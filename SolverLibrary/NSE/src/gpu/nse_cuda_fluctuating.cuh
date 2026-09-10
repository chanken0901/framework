// Landau-Lifshitz pair matching mod_nse_fluctuating, including Philox4x32-10.
// No host staging of stresses, heat fluxes, or stochastic increments.
struct FhView {
  double beta, gamma, mu, kappa, h[3];
  int global_ny, global_nz, y_start, z_start, step, seed;
};

__device__ inline double fh_normal_cuda(unsigned int i, unsigned int j,
    unsigned int k, unsigned int step, unsigned int seed, unsigned int stream) {
  unsigned int c[4]={i,j,k,step}, key[2]={seed,stream};
  for (int r=0;r<10;++r) {
    const unsigned int h0=__umulhi(3528531795u,c[0]), l0=3528531795u*c[0];
    const unsigned int h1=__umulhi(3449720151u,c[2]), l1=3449720151u*c[2];
    const unsigned int n0=h1^c[1]^key[0], n2=h0^c[3]^key[1];
    c[0]=n0; c[1]=l1; c[2]=n2; c[3]=l0;
    key[0]+=2654435769u; key[1]+=3144134277u;
  }
  const double u=(static_cast<double>(c[0])+0.5)/4294967296.0;
  const double v=(static_cast<double>(c[1])+0.5)/4294967296.0;
  return sqrt(-2.0*log(u))*cos(6.2831853071795864769*v);
}

__device__ inline int fh_wrap(int i,int n) { return (i%n+n)%n; }

__device__ inline void fh_primitive(const double* q,GridView g,std::size_t c,
    double gamma,double* u,double& t) {
  const double rho=q[c];
  double kinetic=0;
  for(int v=0;v<3;++v) { u[v]=q[c+(v+1)*g.cell_count]/rho; kinetic+=u[v]*u[v]; }
  t=(gamma-1.0)*(q[c+4*g.cell_count]/rho-0.5*kinetic);
}

__global__ void fh_flux_kernel(const double* q,double* flux,GridView g,FhView f,bool noise,double dt) {
  const std::size_t c=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
  if(c>=g.cell_count) return;
  const int i=c%g.nx_total,j=(c/g.nx_total)%g.ny_total,k=c/(static_cast<std::size_t>(g.nx_total)*g.ny_total);
  if(i<g.nghost-1 || i>g.nghost+g.nx || j<g.nghost-1 || j>g.nghost+g.ny ||
      k<g.nghost-1 || k>g.nghost+g.nz) return;
  const std::size_t stride[3]={1,static_cast<std::size_t>(g.nx_total),
      static_cast<std::size_t>(g.nx_total)*g.ny_total};
  double u[3],t,vel[3][3],grad[3][3],heat[3],s[3][3];
  fh_primitive(q,g,c,f.gamma,u,t);
  for(int d=0;d<3;++d) {
    double up[3],tp; fh_primitive(q,g,c+stride[d],f.gamma,up,tp);
    for(int v=0;v<3;++v) { grad[v][d]=(up[v]-u[v])/f.h[d]; vel[v][d]=0.5*(up[v]+u[v]); }
    heat[d]=f.kappa*(tp-t)/f.h[d];
  }
  if(noise) {
    double z[9];
    for(int v=0;v<9;++v) z[v]=fh_normal_cuda(fh_wrap(i-g.nghost,g.nx),
        fh_wrap(j-g.nghost+f.y_start,f.global_ny),fh_wrap(k-g.nghost+f.z_start,f.global_nz),
        f.step,f.seed,v+1);
    const double scale=2*f.beta/(f.h[0]*f.h[1]*f.h[2]*dt);
    const double a=sqrt(scale*f.mu*t),mean=(z[0]+z[1]+z[2])/3.0;
    for(int v=0;v<3;++v) { s[v][v]=sqrt(2.0)*a*(z[v]-mean); heat[v]=sqrt(scale*f.kappa*t*t)*z[v+6]; }
    s[0][1]=s[1][0]=a*z[3]; s[0][2]=s[2][0]=a*z[4]; s[1][2]=s[2][1]=a*z[5];
  } else {
    const double div=grad[0][0]+grad[1][1]+grad[2][2];
    for(int v=0;v<3;++v) for(int d=0;d<3;++d)
      s[v][d]=f.mu*(grad[v][d]+grad[d][v])-(v==d?2*f.mu*div/3.0:0.0);
  }
  for(int d=0;d<3;++d) {
    double energy=heat[d];
    for(int v=0;v<3;++v) { flux[c+(4*d+v)*g.cell_count]=s[v][d]; energy+=vel[v][d]*s[v][d]; }
    flux[c+(4*d+3)*g.cell_count]=energy;
  }
}

__global__ void fh_divergence_kernel(const double* flux,double* rhs,GridView g,FhView f,bool assign) {
  const std::size_t l=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
  if(l>=static_cast<std::size_t>(g.nx)*g.ny*g.nz) return;
  const auto c=cell_index(g,l%g.nx+g.nghost,(l/g.nx)%g.ny+g.nghost,l/(static_cast<std::size_t>(g.nx)*g.ny)+g.nghost);
  const std::size_t stride[3]={1,static_cast<std::size_t>(g.nx_total),static_cast<std::size_t>(g.nx_total)*g.ny_total};
  if(assign) rhs[c]=0;
  for(int v=0;v<4;++v) {
    double value=0;
    for(int d=0;d<3;++d) value+=(flux[c+(4*d+v)*g.cell_count]-flux[c-stride[d]+(4*d+v)*g.cell_count])/f.h[d];
    if(assign) rhs[c+(v+1)*g.cell_count]=value;
    else rhs[c+(v+1)*g.cell_count]+=value;
  }
}

__global__ void fh_kick_kernel(double* q,const double* noise,GridView g,double dt) {
  const std::size_t l=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
  if(l>=static_cast<std::size_t>(g.nx)*g.ny*g.nz) return;
  const auto c=cell_index(g,l%g.nx+g.nghost,(l/g.nx)%g.ny+g.nghost,l/(static_cast<std::size_t>(g.nx)*g.ny)+g.nghost);
  for(int v=1;v<5;++v) q[c+v*g.cell_count]+=dt*noise[c+v*g.cell_count];
}

bool launch_fh(NseCudaContext* c,bool noise,double dt) {
  FhView f{c->fh_beta,c->gamma,1.0/c->reynolds,
      c->gamma/((c->gamma-1)*c->reynolds*c->prandtl),{c->dx,c->dy,c->dz},
      c->global_ny,c->global_nz,c->global_y_start,c->global_z_start,c->fh_step,c->fh_seed};
  fh_flux_kernel<<<block_count(c->grid.cell_count),256>>>(c->q,c->fh_flux,c->grid,f,noise,dt);
  if(!check_cuda(cudaGetLastError(),"LLNS stress/heat flux")) return false;
  fh_divergence_kernel<<<block_count(c->physical_count),256>>>(c->fh_flux,noise?c->fh_noise:c->rhs,c->grid,f,noise);
  return check_cuda(cudaGetLastError(),"LLNS conservative divergence");
}
