// Resident multicomponent backend. Host transfers are explicit at this C ABI.
// No CUDA dependency is added to the existing CPU or single-component paths.
#include <cuda_runtime.h>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <climits>
#include <cstdint>
#include <new>
#ifdef _WIN32
#define API extern "C" __declspec(dllexport)
#else
#define API extern "C"
#endif
namespace {
constexpr int S=64, V=S+4, B=128;
struct Model {
  int n[3],ns,mapped,bc[6]; // 0 periodic,1 reflective,2 dirichlet,3 relaxation,4 MPI
  double origin[3],h[3],xrange[2],nozzle[4],lengths[3],relaxation,length;
  double cfl,diffusion_cfl,chemistry_cfl;
  double thermo[4],species[S][16],transport[2],diff[S],chem[3],reaction[S][3],reference[6][V];
};
struct Grid {
  Model *m; int n[3],ns,nv; size_t ng;
  double *q,*q0,*rhs,*p,*grad,*flux,*io,*halo,*dt,*stats;
  int *error,*steps;
};
struct Context { Grid g{}; Model host{}; size_t cells=0; int uploads=0,downloads=0; };
bool ok(cudaError_t e,const char*where) {
  if(e==cudaSuccess) return true;
  std::fprintf(stderr,"MC CUDA: %s: %s\n",where,cudaGetErrorString(e)); return false;
}
__device__ size_t at(const Grid&g,int i,int j,int k) {return i+size_t(g.n[0])*(j+size_t(g.n[1])*k);}
__device__ void xyz(const Grid&g,size_t t,int*p) {
  p[0]=t%g.n[0];t/=g.n[0];p[1]=t%g.n[1];p[2]=t/g.n[1];
}
__device__ double dot(const double*a,const double*b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
__device__ void bad(Grid g,int code) {atomicCAS(g.error,0,code);}
__device__ double height(const Model&m,double x) {
  if(!m.mapped) return 1;
  double t=m.nozzle[3],a=x<=t?m.nozzle[0]:m.nozzle[2],b=x<=t?m.xrange[0]:m.xrange[1];
  return m.nozzle[1]+(a-m.nozzle[1])*pow((x-t)/(b-t),2);
}
__device__ void heights(const Model&m,int i,double&hl,double&hr) {
  double x=m.origin[0]+(i-2)*m.h[0];hl=height(m,x);hr=height(m,x+m.h[0]);
}
__device__ double volume(const Model&m,int i) {
  double l,r;heights(m,i,l,r);return m.h[0]*m.h[1]*m.h[2]*(l+r)*.5;
}
__device__ void center(const Model&m,const int*p,double*x) {
  double l,r;heights(m,p[0],l,r);
  for(int d=0;d<3;++d)x[d]=m.origin[d]+(p[d]-1.5)*m.h[d];x[1]*=(l+r)*.5;
}
// f is the left cell's padded index, i.e. physical face f[d]-1.
__device__ void area(const Model&m,const int*f,int d,double*a) {
  a[0]=a[1]=a[2]=0;double l,r;heights(m,f[0],l,r);
  if(d==0)a[0]=height(m,m.origin[0]+(f[0]-1)*m.h[0])*m.h[1]*m.h[2];
  if(d==1){a[0]=-(m.origin[1]+(f[1]-1)*m.h[1])*(r-l)*m.h[2];a[1]=m.h[0]*m.h[2];}
  if(d==2)a[2]=m.h[0]*m.h[1]*(l+r)*.5;
}
__device__ double energy(const Model&m,int s,double t) {
  const double*z=m.species[s];const double*a=z+(t<=z[1]?2:9);
  double t2=t*t,t3=t2*t,t4=t3*t,t5=t4*t;
  return m.thermo[0]/z[0]*((a[0]-1)*t+.5*a[1]*t2+a[2]*t3/3+.25*a[3]*t4+.2*a[4]*t5+a[5]);
}
__device__ double cp(const Model&m,int s,double t) {
  const double*z=m.species[s];const double*a=z+(t<=z[1]?2:9);
  return m.thermo[0]/z[0]*(a[0]+a[1]*t+a[2]*t*t+a[3]*t*t*t+a[4]*t*t*t*t);
}
__device__ double mixture_e(const Model&m,const double*y,double t) {
  double e=0;for(int s=0;s<m.ns;++s)e+=y[s]*energy(m,s,t);return e;
}
__device__ double temperature(Grid g,const double*q,double rho,const double*y) {
  const Model&m=*g.m;double target=(q[m.ns+3]-.5*dot(q+m.ns,q+m.ns)/rho)/rho;
  double lo=m.thermo[1],hi=m.thermo[2],el=mixture_e(m,y,lo),eh=mixture_e(m,y,hi);
  double tol=m.thermo[3],scale=fmax(fmax(fabs(target),fabs(el)),fmax(fabs(eh),1.));
  if(!isfinite(target)||target<el-tol*scale||target>eh+tol*scale){bad(g,2);return lo;}
  scale=fmax(fabs(target),1.);
  if(fabs(target-el)<=tol*scale)return lo;
  if(fabs(target-eh)<=tol*scale)return hi;
  for(int it=0;it<100;++it){double mid=.5*(lo+hi),e=mixture_e(m,y,mid);
    if(e<target)lo=mid;else hi=mid;
    if(fabs(e-target)<=tol*scale||hi-lo<=tol*fmax(mid,1.))break;
  }
  return .5*(lo+hi); // Matches the CPU provider's bracket convention.
}
// primitive = u(3),T,Y(ns),rho,p,c
__device__ void primitive(Grid g,const double*q,double*v) {
  const Model&m=*g.m;double rho=0;
  for(int a=0;a<g.nv;++a)if(!isfinite(q[a]))bad(g,1);
  for(int s=0;s<m.ns;++s){rho+=q[s];if(q[s]<-1e-12)bad(g,1);}
  if(!(rho>1e-12)){bad(g,1);rho=1;}
  for(int s=0;s<m.ns;++s)v[4+s]=q[s]/rho;
  for(int d=0;d<3;++d)v[d]=q[m.ns+d]/rho;
  v[3]=temperature(g,q,rho,v+4);double r=0,c=0;
  for(int s=0;s<m.ns;++s){r+=v[4+s]*m.thermo[0]/m.species[s][0];c+=v[4+s]*cp(m,s,v[3]);}
  v[g.nv]=rho;v[g.nv+1]=rho*r*v[3];v[g.nv+2]=sqrt(c/(c-r)*r*v[3]);
  if(!(c>r)||!(v[g.nv+1]>0)||!isfinite(v[g.nv+2]))bad(g,3);
}
__device__ void conservative(Grid g,double rho,const double*u,double p,const double*y,double*q) {
  const Model&m=*g.m;double r=0;for(int s=0;s<m.ns;++s){q[s]=rho*y[s];r+=y[s]*m.thermo[0]/m.species[s][0];}
  double t=p/(rho*r);if(!isfinite(t)||t<m.thermo[1]||t>m.thermo[2])bad(g,2);
  for(int d=0;d<3;++d)q[m.ns+d]=rho*u[d];
  q[m.ns+3]=rho*(mixture_e(m,y,t)+.5*dot(u,u));
}
__device__ double relax(double x,double ref,double eig,double a){return eig>=0?x:(1-a)*x+a*ref;}
__device__ void boundary(Grid g,const double*q,double*out,int face,const int*idx) {
  const Model&m=*g.m;int d=face/2;double normal[3],a[3];int f[3]={idx[0],idx[1],idx[2]};
  f[d]=face%2?m.n[d]+1:1;area(m,f,d,a);double mag=sqrt(dot(a,a));
  for(int j=0;j<3;++j)normal[j]=a[j]/mag;
  for(int v=0;v<g.nv;++v)out[v]=q[v];
  if(m.bc[face]==1){double un=dot(q+m.ns,normal);for(int j=0;j<3;++j)out[m.ns+j]-=2*un*normal[j];return;}
  if(m.bc[face]==2){for(int v=0;v<g.nv;++v)out[v]=m.reference[face][v];return;}
  double vi[V+3],vr[V+3],delta=0,scale=1;primitive(g,q,vi);primitive(g,m.reference[face],vr);
  for(int v=0;v<g.nv;++v){delta=fmax(delta,fabs(q[v]-m.reference[face][v]));scale=fmax(scale,fabs(m.reference[face][v]));}
  if(delta<=100*DBL_EPSILON*scale)return;
  double sign=face%2?1:-1,ui=sign*dot(vi,normal),ur=sign*dot(vr,normal);
  double ci=vi[g.nv+2],cr=vr[g.nv+2];
  double gamma=.5*(ci*ci*vi[g.nv]/vi[g.nv+1]+cr*cr*vr[g.nv]/vr[g.nv+1]);
  double alpha=1-exp(-m.relaxation*.5*m.h[d]/(m.length>0?m.length:m.lengths[d]));
  if(ui+ci<0)alpha=1;
  double jm=relax(ui-2*ci/(gamma-1),ur-2*cr/(gamma-1),ui-ci,alpha);
  double jp=relax(ui+2*ci/(gamma-1),ur+2*cr/(gamma-1),ui+ci,alpha);
  double entropy=relax(vi[g.nv+1]/pow(vi[g.nv],gamma),vr[g.nv+1]/pow(vr[g.nv],gamma),ui,alpha);
  double cb=.25*(gamma-1)*(jp-jm),ub=.5*(jp+jm);
  if(!(cb>0)||!(entropy>0)){bad(g,4);return;}
  double rho=pow(cb*cb/(gamma*entropy),1/(gamma-1)),p=entropy*pow(rho,gamma),u[3],y[S];
  for(int j=0;j<3;++j)u[j]=relax(vi[j]-sign*ui*normal[j],vr[j]-sign*ur*normal[j],ui,alpha)+sign*ub*normal[j];
  for(int s=0;s<m.ns;++s)y[s]=relax(vi[4+s],vr[4+s],ui,alpha);
  conservative(g,rho,u,p,y,out);
}
__global__ void boundary_kernel(Grid g,int d) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);
  if(p[d]>=2&&p[d]<g.n[d]-2)return;
  int face=2*d+(p[d]>=g.n[d]-2);if(g.m->bc[face]==4)return;
  int src[3]={p[0],p[1],p[2]};
  // Clamp transverse indices during early passes; y then z replace valid corners.
  for(int j=d+1;j<3;++j)src[j]=max(2,min(g.n[j]-3,src[j]));
  if(g.m->bc[face]==0)src[d]=2+((p[d]-2)%g.m->n[d]+g.m->n[d])%g.m->n[d];
  else src[d]=face%2?g.n[d]-3:2;
  size_t from=at(g,src[0],src[1],src[2]);double q[V],out[V];
  for(int v=0;v<g.nv;++v)q[v]=g.q[from+g.ng*v];
  if(g.m->bc[face]==0)for(int v=0;v<g.nv;++v)out[v]=q[v];else boundary(g,q,out,face,src);
  for(int v=0;v<g.nv;++v)g.q[t+g.ng*v]=out[v];
}
__global__ void primitive_kernel(Grid g) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;double q[V],v[V+3];
  for(int s=0;s<g.nv;++s)q[s]=g.q[t+g.ng*s];primitive(g,q,v);
  for(int s=0;s<g.nv+3;++s)g.p[t+g.ng*s]=v[s];
}
__global__ void gradient_kernel(Grid g) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);
  if(p[0]==0||p[1]==0||p[2]==0||p[0]==g.n[0]-1||p[1]==g.n[1]-1||p[2]==g.n[2]-1)return;
  size_t stride[3]={1,size_t(g.n[0]),size_t(g.n[0])*g.n[1]};const Model&m=*g.m;
  double hl,hr;heights(m,p[0],hl,hr);double hb=.5*(hl+hr),eta=m.origin[1]+(p[1]-1.5)*m.h[1];
  for(int v=0;v<g.nv;++v){double r[3];for(int d=0;d<3;++d)r[d]=(g.p[t+stride[d]+g.ng*v]-g.p[t-stride[d]+g.ng*v])/(2*m.h[d]);
    g.grad[t+g.ng*(v+g.nv*0)]=r[0]-r[1]*eta*(hr-hl)/(m.h[0]*hb);
    g.grad[t+g.ng*(v+g.nv*1)]=r[1]/hb;g.grad[t+g.ng*(v+g.nv*2)]=r[2];
  }
}
__global__ void face_kernel(Grid g,int d) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);
  for(int a=0;a<3;++a)if(p[a]<(a==d?1:2)||p[a]>g.n[a]-3)return;
  const Model&m=*g.m;int l[3]={p[0],p[1],p[2]},r[3]={p[0],p[1],p[2]};++r[d];
  int wall=-1;if(p[d]==1&&m.bc[2*d]!=0&&m.bc[2*d]!=4)wall=2*d;
  if(p[d]==m.n[d]+1&&m.bc[2*d+1]!=0&&m.bc[2*d+1]!=4)wall=2*d+1;
  size_t il=at(g,l[0],l[1],l[2]),ir=at(g,r[0],r[1],r[2]);int gl[3]={l[0],l[1],l[2]},gr[3]={r[0],r[1],r[2]};
  if(wall>=0){if(wall%2)gr[d]=gl[d];else gl[d]=gr[d];}
  else {
    if(p[d]==1&&m.bc[2*d]==0)gl[d]=m.n[d]+1;
    if(p[d]==m.n[d]+1&&m.bc[2*d+1]==0)gr[d]=2;
  }
  size_t igl=at(g,gl[0],gl[1],gl[2]),igr=at(g,gr[0],gr[1],gr[2]);
  double a[3],normal[3];area(m,p,d,a);double mag=sqrt(dot(a,a));for(int j=0;j<3;++j)normal[j]=a[j]/mag;
  double vl[V],vr[V],v[V],gf[V][3],ql[V],qr[V];
  for(int s=0;s<g.nv;++s){vl[s]=g.p[il+g.ng*s];vr[s]=g.p[ir+g.ng*s];v[s]=.5*(vl[s]+vr[s]);ql[s]=g.q[il+g.ng*s];qr[s]=g.q[ir+g.ng*s];
    for(int j=0;j<3;++j)gf[s][j]=.5*(g.grad[igl+g.ng*(s+g.nv*j)]+g.grad[igr+g.ng*(s+g.nv*j)]);
  }
  if(m.mapped){double delta[3],xl[3],xr[3];center(m,gr,xr);center(m,gl,xl);for(int j=0;j<3;++j)delta[j]=xr[j]-xl[j];
    bool seam=(p[d]==1&&m.bc[2*d]==0)||(p[d]==m.n[d]+1&&m.bc[2*d+1]==0);
    if(wall>=0||seam){double dist=.5*(volume(m,gl[0])+volume(m,gr[0]))/mag;for(int j=0;j<3;++j)delta[j]=normal[j]*dist;}
    double dist=dot(delta,normal);if(!(dist>0)){bad(g,5);return;}
    for(int s=0;s<g.nv;++s){double corr=(vr[s]-vl[s]-dot(gf[s],delta))/dist;for(int j=0;j<3;++j)gf[s][j]+=corr*normal[j];}
  }else if(wall>=0)for(int s=0;s<g.nv;++s)gf[s][d]=(vr[s]-vl[s])/m.h[d];
  double rl=g.p[il+g.ng*g.nv],rr=g.p[ir+g.ng*g.nv],rho=.5*(rl+rr),ul=dot(vl,normal),ur=dot(vr,normal);
  double pl=g.p[il+g.ng*(g.nv+1)],pr=g.p[ir+g.ng*(g.nv+1)];
  double wave=fmax(fabs(ul)+g.p[il+g.ng*(g.nv+2)],fabs(ur)+g.p[ir+g.ng*(g.nv+2)]);
  double js[S],sum=0,cpv=0,stress[3]={0,0,0},visc[V]={0};
  for(int s=0;s<m.ns;++s){js[s]=-rho*m.diff[s]*dot(gf[4+s],normal);sum+=js[s];cpv+=v[4+s]*cp(m,s,v[3]);}
  double residual=0;for(int s=0;s<m.ns;++s){js[s]-=v[4+s]*sum;residual+=js[s];}js[m.ns-1]-=residual;
  double div=gf[0][0]+gf[1][1]+gf[2][2],mu=m.transport[0];
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)stress[i]+=mu*(gf[i][j]+gf[j][i]-(i==j?2*div/3:0))*normal[j];
  for(int s=0;s<m.ns;++s)visc[s]=-js[s];for(int j=0;j<3;++j)visc[m.ns+j]=stress[j];
  visc[m.ns+3]=dot(stress,v)+mu*cpv/m.transport[1]*dot(gf[3],normal);
  for(int s=0;s<m.ns;++s)visc[m.ns+3]-=(energy(m,s,v[3])+m.thermo[0]/m.species[s][0]*v[3])*js[s];
  if(m.mapped&&wall>=0&&m.bc[wall]==1){for(int s=0;s<m.ns;++s)visc[s]=0;double sn=dot(stress,normal);for(int j=0;j<3;++j)visc[m.ns+j]=sn*normal[j];visc[m.ns+3]=0;}
  for(int s=0;s<g.nv;++s){double fl=ql[s]*ul,fr=qr[s]*ur;
    if(s>=m.ns&&s<m.ns+3){fl+=pl*normal[s-m.ns];fr+=pr*normal[s-m.ns];}
    if(s==m.ns+3){fl=(ql[s]+pl)*ul;fr=(qr[s]+pr)*ur;}
    g.flux[t+g.ng*(s+g.nv*d)]=mag*(visc[s]-.5*(fl+fr)+.5*wave*(qr[s]-ql[s]));
  }
}
__device__ bool owned(Grid g,const int*p){return p[0]>=2&&p[0]<g.n[0]-2&&p[1]>=2&&p[1]<g.n[1]-2&&p[2]>=2&&p[2]<g.n[2]-2;}
__global__ void update_kernel(Grid g,double dt,int stage) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);if(!owned(g,p))return;
  size_t stride[3]={1,size_t(g.n[0]),size_t(g.n[0])*g.n[1]};double vol=volume(*g.m,p[0]);
  for(int v=0;v<g.nv;++v){double rhs=0;for(int d=0;d<3;++d)rhs+=g.flux[t+g.ng*(v+g.nv*d)]-g.flux[t-stride[d]+g.ng*(v+g.nv*d)];
    g.rhs[t+g.ng*v]=rhs/vol;if(stage==0)continue;
    double tmp=g.q[t+g.ng*v]+dt*rhs/vol,q0=g.q0[t+g.ng*v];
    g.q[t+g.ng*v]=stage==1?tmp:(stage==2?.75*q0+.25*tmp:q0/3+(2./3)*tmp);
  }
}
__device__ void source(Grid g,const double*q,double*out) {
  const Model&m=*g.m;double v[V+3];primitive(g,q,v);double lograte=log(m.chem[0])+m.chem[1]*log(v[3])-m.chem[2]/v[3];
  for(int s=0;s<g.nv;++s)out[s]=0;
  for(int s=0;s<m.ns;++s)if(m.reaction[s][2]>0){double c=q[s]/m.species[s][0];if(c<=0)return;lograte+=m.reaction[s][2]*log(c);}
  if(lograte<=log(DBL_MIN))return;if(lograte>=log(DBL_MAX)){bad(g,6);return;}
  double rate=exp(lograte),sum=0;int last=-1;
  for(int s=0;s<m.ns;++s){out[s]=m.species[s][0]*(m.reaction[s][1]-m.reaction[s][0])*rate;sum+=out[s];if(m.reaction[s][1]>0)last=s;}
  if(last<0){bad(g,6);return;}out[last]-=sum;
}
__device__ double chemistry_dt(Grid g,const double*q,const double*s,double cap) {
  for(int i=0;i<g.ns;++i)if(s[i]<0)cap=fmin(cap,g.m->chemistry_cfl*q[i]/(-s[i]));return cap;
}
__global__ void chemistry_kernel(Grid g,double duration,int limit) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);if(!owned(g,p))return;
  double q[V],q0[V],tmp[V],s[V],v[V+3],remaining=duration;int steps=0;
  for(int j=0;j<g.nv;++j)q[j]=g.q[t+g.ng*j];
  primitive(g,q,v);
  while(remaining>0){if(steps>=limit){bad(g,7);break;}
    source(g,q,s);double dt=chemistry_dt(g,q,s,remaining);
    if(!(dt>0)||!isfinite(dt)){bad(g,8);break;}
    for(int j=0;j<g.nv;++j){q0[j]=q[j];tmp[j]=q[j]+dt*s[j];}primitive(g,tmp,v);
    source(g,tmp,s);for(int j=0;j<g.nv;++j)tmp[j]=.75*q0[j]+.25*(tmp[j]+dt*s[j]);primitive(g,tmp,v);
    source(g,tmp,s);for(int j=0;j<g.nv;++j)q[j]=q0[j]/3+(2./3)*(tmp[j]+dt*s[j]);primitive(g,q,v);
    ++steps;remaining-=dt;if(remaining<=fmax(16*DBL_EPSILON*duration,DBL_MIN))remaining=0;
    if(*g.error)break;
  }
  for(int j=0;j<g.nv;++j)g.q[t+g.ng*j]=q[j];atomicMax(g.steps,steps);
}
__global__ void dt_kernel(Grid g,int with_chemistry) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);if(!owned(g,p))return;
  const Model&m=*g.m;double q[V],v[V+3];for(int j=0;j<g.nv;++j)q[j]=g.q[t+g.ng*j];primitive(g,q,v);
  double rate=0,metric=0,vol=volume(m,p[0]);
  for(int d=0;d<3;++d)for(int side=0;side<2;++side){int f[3]={p[0],p[1],p[2]};f[d]-=side;double a[3];area(m,f,d,a);double mag=sqrt(dot(a,a));rate+=.5*(fabs(dot(v,a))+v[g.nv+2]*mag)/vol;metric+=.5*mag/vol;}
  double cpv=0,r=0,diff=0;for(int s=0;s<m.ns;++s){cpv+=v[4+s]*cp(m,s,v[3]);r+=v[4+s]*m.thermo[0]/m.species[s][0];diff=fmax(diff,m.diff[s]);}
  diff=fmax(diff,fmax(m.transport[0]/v[g.nv],m.transport[0]*cpv/(m.transport[1]*v[g.nv]*(cpv-r))));
  double geom=m.mapped?metric*metric:1/(m.h[0]*m.h[0])+1/(m.h[1]*m.h[1])+1/(m.h[2]*m.h[2]);
  double dt=fmin((with_chemistry==2?1.:m.cfl)/rate,m.diffusion_cfl/(2*diff*geom));
  if(with_chemistry==1){double s[V];source(g,q,s);dt=chemistry_dt(g,q,s,dt);}
  if(!(dt>0)||!isfinite(dt)){bad(g,8);return;}
  atomicMin(reinterpret_cast<unsigned long long*>(g.dt),__double_as_longlong(dt));
}
__global__ void io_kernel(Grid g,double*data,int direction,int rhs) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x,total=size_t(g.m->n[0])*g.m->n[1]*g.m->n[2];if(t>=total)return;
  size_t x=t;int i=x%g.m->n[0]+2;x/=g.m->n[0];int j=x%g.m->n[1]+2,k=x/g.m->n[1]+2;size_t padded=at(g,i,j,k);
  for(int v=0;v<g.nv;++v){if(direction)g.q[padded+g.ng*v]=data[t+total*v];else data[t+total*v]=(rhs?g.rhs:g.q)[padded+g.ng*v];}
}
__global__ void seed_halo_kernel(Grid g) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;int p[3];xyz(g,t,p);if(owned(g,p))return;
  for(int d=0;d<3;++d)p[d]=max(2,min(g.n[d]-3,p[d]));size_t src=at(g,p[0],p[1],p[2]);
  for(int v=0;v<g.nv;++v)g.q[t+g.ng*v]=g.q[src+g.ng*v];
}
__device__ void minimum(double*address,double value){auto*p=reinterpret_cast<unsigned long long*>(address);auto old=*p;
  while(value<__longlong_as_double(old)){auto expected=old;old=atomicCAS(p,expected,__double_as_longlong(value));if(old==expected)break;}}
__global__ void statistics_kernel(Grid g){size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=g.ng)return;
  int p[3];xyz(g,t,p);if(!owned(g,p))return;double q[V],v[V+3],vol=volume(*g.m,p[0]),minspecies=DBL_MAX;
  for(int s=0;s<g.nv;++s){q[s]=g.q[t+g.ng*s];atomicAdd(g.stats+s,q[s]*vol);if(s<g.ns)minspecies=fmin(minspecies,q[s]);}
  primitive(g,q,v);minimum(g.stats+g.nv,minspecies);minimum(g.stats+g.nv+1,v[g.nv]);
  minimum(g.stats+g.nv+2,v[g.nv+1]);minimum(g.stats+g.nv+3,v[3]);
}
void release(Context*c){if(!c)return;auto&g=c->g;cudaFree(g.m);cudaFree(g.q);cudaFree(g.q0);cudaFree(g.rhs);cudaFree(g.p);cudaFree(g.grad);cudaFree(g.flux);cudaFree(g.io);cudaFree(g.halo);cudaFree(g.dt);cudaFree(g.stats);cudaFree(g.error);cudaFree(g.steps);delete c;}
bool reset(Context*c){return ok(cudaMemset(c->g.error,0,sizeof(int)),"reset error");}
bool finish(Context*c,const char*where){int e=0;if(!ok(cudaGetLastError(),where)||!ok(cudaMemcpy(&e,c->g.error,sizeof(int),cudaMemcpyDeviceToHost),where))return false;
  if(e){std::fprintf(stderr,"MC CUDA %s: numerical error %d (1 state,2 NASA range,3 thermo,4 boundary,5 metric,6 rate,7 substep limit,8 dt)\n",where,e);return false;}return true;}
int blocks(size_t n){return int((n+B-1)/B);}
void boundaries(Context*c){for(int d=0;d<3;++d)boundary_kernel<<<blocks(c->g.ng),B>>>(c->g,d);}
__global__ void halo_kernel(Grid g,int d,int side,int unpack,size_t count) {
  size_t t=blockIdx.x*size_t(blockDim.x)+threadIdx.x;if(t>=count)return;
  int dims[3]={g.n[0],g.n[1],g.n[2]};dims[d]=2;size_t x=t;int p[3];
  for(int j=0;j<3;++j){p[j]=x%dims[j];x/=dims[j];}int v=int(x);
  p[d]+=unpack?(side?g.n[d]-2:0):(side?g.n[d]-4:2);
  size_t offset=at(g,p[0],p[1],p[2])+g.ng*v;
  if(unpack)g.q[offset]=g.halo[t];else g.halo[t]=g.q[offset];
}
}
API int mc_cuda_create(const Model*m,int device,void**handle) {
  *handle=nullptr;if(!m||m->ns<1||m->ns>S)return 1;int count=0;
  if(!ok(cudaGetDeviceCount(&count),"device count")||device<0||device>=count){std::fprintf(stderr,"MC CUDA: invalid device %d\n",device);return 1;}
  if(!ok(cudaSetDevice(device),"select device"))return 1;auto*c=new(std::nothrow) Context;if(!c)return 1;c->host=*m;auto&g=c->g;
  g.ns=m->ns;g.nv=m->ns+4;g.ng=1;c->cells=1;
  for(int d=0;d<3;++d){if(m->n[d]<1||m->n[d]>INT_MAX-4){release(c);return 1;}g.n[d]=m->n[d]+4;g.ng*=g.n[d];c->cells*=m->n[d];}
  if(g.ng>size_t(INT_MAX)*B||g.ng>SIZE_MAX/(sizeof(double)*(10*g.nv+3))){release(c);return 1;}
  size_t free=0,total=0;cudaMemGetInfo(&free,&total);size_t need=(g.ng*(10*g.nv+3)+c->cells*g.nv)*sizeof(double)+sizeof(Model)+32;
  if(need>free){std::fprintf(stderr,"MC CUDA: requires %zu bytes, free %zu bytes\n",need,free);release(c);return 1;}
  auto alloc=[&](double**p,size_t n){return ok(cudaMalloc(reinterpret_cast<void**>(p),n*sizeof(double)),"allocate resident workspace");};
  bool success=ok(cudaMalloc(reinterpret_cast<void**>(&g.m),sizeof(Model)),"allocate model")&&alloc(&g.q,g.ng*g.nv)&&alloc(&g.q0,g.ng*g.nv)&&alloc(&g.rhs,g.ng*g.nv)&&alloc(&g.p,g.ng*(g.nv+3))&&alloc(&g.grad,g.ng*g.nv*3)&&alloc(&g.flux,g.ng*g.nv*3)&&alloc(&g.io,c->cells*g.nv)&&alloc(&g.dt,1)&&ok(cudaMalloc(reinterpret_cast<void**>(&g.error),sizeof(int)),"allocate error")&&ok(cudaMalloc(reinterpret_cast<void**>(&g.steps),sizeof(int)),"allocate counter")&&ok(cudaMemcpy(g.m,m,sizeof(Model),cudaMemcpyHostToDevice),"upload model");
  if(success)success=alloc(&g.stats,g.nv+4);
  if(!success){release(c);return 1;}cudaMemset(g.q,0,g.ng*g.nv*sizeof(double));*handle=c;return 0;
}
API void mc_cuda_destroy(void*handle){release(static_cast<Context*>(handle));}
API int mc_cuda_upload(void*h,const double*q){auto*c=static_cast<Context*>(h);if(!c||!q||!reset(c))return 1;
  ++c->uploads;
  if(!ok(cudaMemcpy(c->g.io,q,c->cells*c->g.nv*sizeof(double),cudaMemcpyHostToDevice),"initial upload"))return 1;
  io_kernel<<<blocks(c->cells),B>>>(c->g,c->g.io,1,0);seed_halo_kernel<<<blocks(c->g.ng),B>>>(c->g);
  boundaries(c);return finish(c,"upload")?0:1;}
API int mc_cuda_download(void*h,double*q,int rhs){auto*c=static_cast<Context*>(h);if(!c||!q)return 1;
  ++c->downloads;
  io_kernel<<<blocks(c->cells),B>>>(c->g,c->g.io,0,rhs);return ok(cudaMemcpy(q,c->g.io,c->cells*c->g.nv*sizeof(double),cudaMemcpyDeviceToHost),"output download")?0:1;}
API int mc_cuda_timestep(void*h,int chem,double*dt){auto*c=static_cast<Context*>(h);if(!c||!reset(c))return 1;double cap=DBL_MAX;
  if(!ok(cudaMemcpy(c->g.dt,&cap,sizeof(double),cudaMemcpyHostToDevice),"reset dt"))return 1;
  dt_kernel<<<blocks(c->g.ng),B>>>(c->g,chem);if(!finish(c,"timestep"))return 1;return ok(cudaMemcpy(dt,c->g.dt,sizeof(double),cudaMemcpyDeviceToHost),"scalar dt")?0:1;}
API int mc_cuda_chemistry(void*h,double duration,int limit,int*steps){auto*c=static_cast<Context*>(h);if(!c||!std::isfinite(duration)||duration<0||limit<1||!reset(c))return 1;
  cudaMemset(c->g.steps,0,sizeof(int));chemistry_kernel<<<blocks(c->g.ng),B>>>(c->g,duration,limit);
  if(!finish(c,"chemistry"))return 1;return ok(cudaMemcpy(steps,c->g.steps,sizeof(int),cudaMemcpyDeviceToHost),"scalar substeps")?0:1;}
API int mc_cuda_begin(void*h){auto*c=static_cast<Context*>(h);if(!c)return 1;return ok(cudaMemcpy(c->g.q0,c->g.q,c->g.ng*c->g.nv*sizeof(double),cudaMemcpyDeviceToDevice),"save resident RK state")?0:1;}
API int mc_cuda_stage(void*h,double dt,int stage){auto*c=static_cast<Context*>(h);if(!c||stage<0||stage>3||!std::isfinite(dt)||dt<0||!reset(c))return 1;
  bool distributed=false;for(int f=0;f<6;++f)distributed|=c->host.bc[f]==4;
  if(!distributed)boundaries(c);
  primitive_kernel<<<blocks(c->g.ng),B>>>(c->g);gradient_kernel<<<blocks(c->g.ng),B>>>(c->g);
  for(int d=0;d<3;++d)face_kernel<<<blocks(c->g.ng),B>>>(c->g,d);
  update_kernel<<<blocks(c->g.ng),B>>>(c->g,dt,stage);
  // Check every RK update, including the final stage, without downloading state.
  if(stage){double cap=DBL_MAX;cudaMemcpy(c->g.dt,&cap,sizeof(double),cudaMemcpyHostToDevice);dt_kernel<<<blocks(c->g.ng),B>>>(c->g,0);}
  return finish(c,"fluid stage")?0:1;
}
API int mc_cuda_local_boundary(void*h){auto*c=static_cast<Context*>(h);if(!c||!reset(c))return 1;boundaries(c);return finish(c,"local boundary")?0:1;}
API int mc_cuda_statistics(void*h,double*stats){auto*c=static_cast<Context*>(h);if(!c||!stats||!reset(c))return 1;
  double init[V+4]={0};for(int i=c->g.nv;i<c->g.nv+4;++i)init[i]=DBL_MAX;
  if(!ok(cudaMemcpy(c->g.stats,init,(c->g.nv+4)*sizeof(double),cudaMemcpyHostToDevice),"reset statistics"))return 1;
  statistics_kernel<<<blocks(c->g.ng),B>>>(c->g);if(!finish(c,"statistics"))return 1;
  return ok(cudaMemcpy(stats,c->g.stats,(c->g.nv+4)*sizeof(double),cudaMemcpyDeviceToHost),"scalar statistics")?0:1;
}
API int mc_cuda_transfer_counts(void*h,int*counts){auto*c=static_cast<Context*>(h);if(!c||!counts)return 1;counts[0]=c->uploads;counts[1]=c->downloads;return 0;}
API int mc_cuda_halo(void*h,int direction,int side,int unpack,double*data){
  auto*c=static_cast<Context*>(h);if(!c||direction<1||direction>2||side<0||side>1||!data)return 1;
  auto&g=c->g;size_t count=g.ng/g.n[direction]*2*g.nv;
  if(!g.halo){size_t maximum=g.ng/(g.n[1]<g.n[2]?g.n[1]:g.n[2])*2*g.nv;
    if(!ok(cudaMalloc(reinterpret_cast<void**>(&g.halo),maximum*sizeof(double)),"allocate persistent halo"))return 1;}
  if(unpack&&!ok(cudaMemcpy(g.halo,data,count*sizeof(double),cudaMemcpyHostToDevice),"halo upload"))return 1;
  halo_kernel<<<blocks(count),B>>>(g,direction,side,unpack,count);
  if(!unpack)return ok(cudaMemcpy(data,g.halo,count*sizeof(double),cudaMemcpyDeviceToHost),"halo download")?0:1;
  return ok(cudaDeviceSynchronize(),"unpack halo")?0:1;
}
