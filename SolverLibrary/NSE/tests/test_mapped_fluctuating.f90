program test_mapped_fluctuating
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config,update_derived_config
  use mod_model_config, only: nse_config
  use mod_grid_fvm, only: build_stretched_grid
  use mod_nse_fluctuating, only: fh_spacing,add_fh_noise,add_fh_transport
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  integer, parameter :: n=6,g=6,nt=2000
  real(dp) :: q(1-g:n+g,1-g:n+g,1-g:n+g,5),rhs(1-g:n+g,1-g:n+g,1-g:n+g,5)
  real(dp) :: w(n,n,n,3),weight(n,n,n),probe(n,n,n),h(3),x,u(3),work,thermal_work
  real(dp) :: dt,expected,variance,average,projection,tvar,tavg,tp,beta,pi
  integer :: i,j,k,v,t,mode
  sim%nx=n;sim%ny=n;sim%nz=n;sim%nghost=g;sim%grid_mapping='sinh'
  sim%grid_stretch=[.4_dp,1._dp,1.5_dp];sim%mapped_keep6=.true.
  call update_derived_config(sim)
  call build_stretched_grid(sim,1,n,1,n,[.true.,.true.,.true.])
  nse%fh_enabled=.true.;nse%fh_boltzmann_number=1.e-5_dp;nse%reynolds=100
  nse%viscous_scheme='central6';dt=.001_dp;beta=nse%fh_boltzmann_number;pi=acos(-1._dp)
  do mode=0,1
    sim%mapped_keep6=mode==1
    do k=1,n
      do j=1,n
        do i=1,n
          h=fh_spacing(i,j,k,sim);weight(i,j,k)=product(h)
          x=2*pi*real(i+2*j+3*k,dp)/n
          w(i,j,k,:)=[sin(x),cos(x),sin(x+.7_dp)];probe(i,j,k)=cos(x)
        end do
      end do
    end do
    do k=1-g,n+g
      do j=1-g,n+g
        do i=1-g,n+g
          u=w(modulo(i-1,n)+1,modulo(j-1,n)+1,modulo(k-1,n)+1,:)
          q(i,j,k,:)=[1._dp,u,1/(nse%gamma-1)+.5_dp*sum(u*u)]
        end do
      end do
    end do
    rhs=0;call add_fh_transport(q,rhs,sim,nse,1,n,1,n)
    do v=1,5
      if(abs(sum(weight*rhs(1:n,1:n,1:n,v)))>1.e-11_dp) error stop 'mapped transport conservation'
    end do
    work=0
    do v=1,3
      work=work+sum(weight*w(:,:,:,v)*rhs(1:n,1:n,1:n,v+1))
    end do
    if(work>=0) error stop 'mapped viscous dissipation'
    expected=-2*beta*work/dt
    do k=1-g,n+g
      do j=1-g,n+g
        do i=1-g,n+g
          q(i,j,k,2:4)=0
          q(i,j,k,5)=(1+.1_dp*probe(modulo(i-1,n)+1,modulo(j-1,n)+1,modulo(k-1,n)+1))/(nse%gamma-1)
        end do
      end do
    end do
    rhs=0;call add_fh_transport(q,rhs,sim,nse,1,n,1,n)
    thermal_work=sum(weight*probe*rhs(1:n,1:n,1:n,5))/.1_dp
    if(thermal_work>=0) error stop 'mapped thermal dissipation'
    q(:,:,:,5)=1/(nse%gamma-1)
    average=0;variance=0;tavg=0;tvar=0
    do t=1,nt
      sim%step=t;rhs=0;call add_fh_noise(q,rhs,dt,sim,nse,1,n,1,n)
      do v=1,5
        if(abs(sum(weight*rhs(1:n,1:n,1:n,v)))>1.e-11_dp) error stop 'mapped noise conservation'
      end do
      projection=0
      do v=1,3
        projection=projection+sum(weight*w(:,:,:,v)*rhs(1:n,1:n,1:n,v+1))
      end do
      tp=sum(weight*probe*rhs(1:n,1:n,1:n,5))
      average=average+projection;variance=variance+projection**2;tavg=tavg+tp;tvar=tvar+tp**2
    end do
    variance=variance/nt-(average/nt)**2;tvar=tvar/nt-(tavg/nt)**2
    print *, 'mapped FDT norm, momentum, heat:',mode,variance/expected,tvar/(-2*beta*thermal_work/dt)
    if(abs(variance/expected-1)>0.12_dp) error stop 'mapped momentum FDT'
    if(abs(tvar/(-2*beta*thermal_work/dt)-1)>0.12_dp) error stop 'mapped thermal FDT'
  end do
  print *, 'Mapped fluctuation dissipation tests passed'
end program
