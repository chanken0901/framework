// Conservative orthogonal-grid viscous face flux. Matches CPU fv2, not CENTRAL6.
__device__ inline void fv2_values(const double* field,const GridView& g,const int p[3],double v[4]) {
  const auto c=cell_index(g,p[0],p[1],p[2]);
  for(int a=0;a<4;++a) v[a]=field[c+a*g.cell_count];
}
__device__ inline void fv2_gradient(const double* field,const GridView& g,const int p[3],int axis,double grad[4]) {
  int pm[3]={p[0],p[1],p[2]},pp[3]={p[0],p[1],p[2]};--pm[axis];++pp[axis];
  double vm[4],v0[4],vp[4];fv2_values(field,g,pm,vm);fv2_values(field,g,p,v0);fv2_values(field,g,pp,vp);
  const double dl=g.axis_center[axis][p[axis]]-g.axis_center[axis][pm[axis]];
  const double dr=g.axis_center[axis][pp[axis]]-g.axis_center[axis][p[axis]];
  for(int a=0;a<4;++a) grad[a]=(dr*(v0[a]-vm[a])/dl+dl*(vp[a]-v0[a])/dr)/(dl+dr);
}
__device__ inline void fv2_face(const double* field,const GridView& g,const int p[3],int axis,
    double mu,double heat,double flux[4]) {
  int right[3]={p[0],p[1],p[2]};++right[axis];
  double vl[4],vr[4],vf[4],grad[4][3],gl[4],gr[4];
  fv2_values(field,g,p,vl);fv2_values(field,g,right,vr);
  const double distance=g.axis_center[axis][right[axis]]-g.axis_center[axis][p[axis]];
  const double theta=.5*g.axis_width[axis][p[axis]]/distance;
  for(int a=0;a<4;++a) vf[a]=(1-theta)*vl[a]+theta*vr[a];
  for(int d=0;d<3;++d) {
    if(d==axis) {for(int a=0;a<4;++a) grad[a][d]=(vr[a]-vl[a])/distance;}
    else {
      fv2_gradient(field,g,p,d,gl);fv2_gradient(field,g,right,d,gr);
      for(int a=0;a<4;++a) grad[a][d]=(1-theta)*gl[a]+theta*gr[a];
    }
  }
  const double div=grad[0][0]+grad[1][1]+grad[2][2];
  double tau[3];flux[3]=heat*grad[3][axis];
  for(int a=0;a<3;++a) {
    tau[a]=grad[a][axis]+grad[axis][a]-(a==axis?2.0/3.0*div:0);
    flux[a]=mu*tau[a];flux[3]+=vf[a]*tau[a];
  }
  flux[3]*=mu;
}
__global__ void viscous_fv2_kernel(const double* field,double* rhs,GridView g,
    double mu,double heat,std::size_t count,bool mapped_keep6) {
  const auto linear=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
  if(linear>=count) return;
  const int p[3]={static_cast<int>(linear%g.nx)+g.nghost,
    static_cast<int>((linear/g.nx)%g.ny)+g.nghost,
    static_cast<int>(linear/(static_cast<std::size_t>(g.nx)*g.ny))+g.nghost};
  const auto c=cell_index(g,p[0],p[1],p[2]);double result[4]={0,0,0,0};
  for(int d=0;d<3;++d) {
    int pm[3]={p[0],p[1],p[2]};--pm[d];double plus[4],minus[4];
    fv2_face(field,g,p,d,mu,heat,plus);fv2_face(field,g,pm,d,mu,heat,minus);
    const double width=mapped_keep6?g.axis_keep6_metric[d][p[d]]:g.axis_width[d][p[d]];
    for(int a=0;a<4;++a) result[a]+=(plus[a]-minus[a])/width;
  }
  for(int a=0;a<4;++a) rhs[c+(a+1)*g.cell_count]+=result[a];
}
