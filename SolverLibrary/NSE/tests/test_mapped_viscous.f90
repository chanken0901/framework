program test_mapped_viscous
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_model_config, only: nse_config
  use mod_grid_axis, only: build_axis,build_sinh_axis
  use mod_grid_fvm, only: axis_x,axis_y,axis_z
  use mod_viscous_scheme, only: add_viscous_rhs
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:),rhs(:,:,:,:),edges(:)
  real(dp) :: x,s,u,v,w,t,pi,heat,expected(4),err(3),balance(4),ke,weight
  integer :: pass,n,i,j,k
  pi=acos(-1._dp)
  sim%grid_mapping='sinh';sim%mapped_keep6=.true.;sim%nghost=6
  sim%ny=8;sim%nz=8
  nse%viscous_scheme='central6';nse%reynolds=50;nse%prandtl=.72_dp
  heat=nse%gamma/((nse%gamma-1)*nse%prandtl)
  call build_sinh_axis(8,6,0._dp,2*pi,0._dp,.true.,axis_y)
  call build_sinh_axis(8,6,0._dp,2*pi,0._dp,.true.,axis_z)
  do pass=1,3
    n=32*2**(pass-1);sim%nx=n
    allocate(edges(0:n),q(-5:n+6,-5:14,-5:14,5),rhs(-5:n+6,-5:14,-5:14,5))
    do i=0,n
      s=2*pi*i/n;edges(i)=s+.15_dp*sin(s)
    end do
    call build_axis(edges,6,.true.,axis_x)
    do k=-5,14
      do j=-5,14
        do i=-5,n+6
          x=axis_x%center(i);u=.2_dp*sin(x);v=.1_dp*cos(x);w=.1_dp*sin(2*x);t=1+.1_dp*cos(x)
          q(i,j,k,:)=[1._dp,u,v,w,t/(nse%gamma-1)+.5_dp*(u*u+v*v+w*w)]
        end do
      end do
    end do
    rhs=0
    !$OMP PARALLEL DEFAULT(SHARED)
    call add_viscous_rhs(q,rhs,sim,nse,1,8,1,8)
    !$OMP END PARALLEL
    err(pass)=0;balance=0;ke=0
    do i=1,n
      x=axis_x%center(i);weight=axis_x%keep6_metric(i)
      expected=[-(4._dp/3)*.2_dp*sin(x),-.1_dp*cos(x),-.4_dp*sin(2*x), &
        ((4._dp/3)*.04_dp-.01_dp)*cos(2*x)+.04_dp*cos(4*x)-heat*.1_dp*cos(x)]/nse%reynolds
      err(pass)=max(err(pass),maxval(abs(rhs(i,1,1,2:5)-expected)))
      balance=balance+weight*rhs(i,1,1,2:5)
      ke=ke+weight*dot_product(q(i,1,1,2:4),rhs(i,1,1,2:4))
    end do
    if(maxval(abs(balance))>1.e-12_dp) error stop 'Mapped viscous conservation'
    if(ke>=0) error stop 'Mapped viscosity must dissipate kinetic energy'
    if(maxval(abs(rhs(:,:,:,1)))>0) error stop 'Viscosity changed mass'
    deallocate(q,rhs,edges)
  end do
  print *, 'Mapped viscous analytic errors:',err
  print *, 'Mapped viscous convergence ratios:',err(1)/err(2),err(2)/err(3)
  if(min(err(1)/err(2),err(2)/err(3))<32) error stop 'Mapped viscous sixth-order convergence'
end program
