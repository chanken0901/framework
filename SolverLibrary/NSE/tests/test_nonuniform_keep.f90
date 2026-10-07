program test_nonuniform_keep
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_model_config, only: nse_config
  use mod_grid_axis, only: grid_axis,build_sinh_axis
  use mod_convective_keep, only: compute_keep_face_flux
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(grid_axis) :: axis
  real(dp), allocatable :: q(:,:,:,:),face(:,:),rhs(:,:)
  real(dp) :: x,rho,u,p,rate(5),ke_rate,pi,err(2),exact,volume
  integer :: n,i,l,pass,order
  character(len=32) :: mode
  call get_command_argument(1,mode)
  pi=acos(-1._dp);order=2
  if(mode=='mapped') order=6
  if(mode=='reject_keep6') then
    order=6
    nse%convective_scheme='hybrid'
  end if
  ! Nonconstant volume weights with periodic data: conservation and KE balance.
  n=32
  call configure(n)
  do i=1,n
    x=axis%center(i)
    rho=1+.1_dp*sin(2*pi*x);u=.3_dp+.1_dp*cos(2*pi*x);p=1
    q(i,1,1,:)=[rho,rho*u,0._dp,0._dp,p/(nse%gamma-1)+.5_dp*rho*u*u]
  end do
  call evaluate()
  rate=0;ke_rate=0
  do i=1,n
    rho=q(i,1,1,1);u=q(i,1,1,2)/rho;volume=axis%width(i)
    if(order==6) volume=axis%keep6_metric(i)
    rate=rate+volume*rhs(i,:)
    ke_rate=ke_rate+volume*(u*rhs(i,2)-.5_dp*u*u*rhs(i,1))
  end do
  if(maxval(abs(rate))>1.e-12_dp) error stop 'Nonuniform KEEP conserved integral'
  if(abs(ke_rate)>1.e-12_dp) error stop 'Nonuniform KEEP kinetic energy convection'
  q=0;q(:,:,:,1)=1;q(:,:,:,2)=.3_dp;q(:,:,:,5)=1/(nse%gamma-1)+.045_dp
  call evaluate()
  if(maxval(abs(rhs))>1.e-12_dp) error stop 'KEEP free stream preservation'
  deallocate(q,face,rhs)
  ! Interior smooth stretched-grid consistency. Constant velocity and pressure
  ! give mass residual -u*d(rho)/dx. Boundary closure is not tested here.
  do pass=1,2
    n=32*2**(pass-1)
    call configure(n)
    do i=1,n
      x=axis%center(i);rho=1+.1_dp*sin(2*pi*x);u=.3_dp
      q(i,1,1,:)=[rho,rho*u,0._dp,0._dp,1/(nse%gamma-1)+.5_dp*rho*u*u]
    end do
    call evaluate()
    err(pass)=0
    do i=4,n-3
      if(order==6.and.(axis%center(i)<.2_dp.or.axis%center(i)>.8_dp)) cycle
      exact=-.3_dp*.1_dp*2*pi*cos(2*pi*axis%center(i))
      err(pass)=max(err(pass),abs(rhs(i,1)-exact))
    end do
    deallocate(q,face,rhs)
  end do
  if(order==2.and.err(1)/err(2)<3.0_dp) error stop 'KEEP2 smooth-grid convergence'
  if(order==6) print *, 'KEEP6 sinh interior errors and ratio:',err,err(1)/err(2)
  if(order==6.and.err(1)/err(2)<32.0_dp) error stop 'KEEP6 smooth-grid convergence'
  print *, '[OK] stretched KEEP order=',order,': conservation, KE convection, free stream, interior convergence',err
contains
  subroutine configure(n)
    integer, intent(in) :: n
    sim%nx=n;sim%ny=1;sim%nz=1;sim%nghost=3;sim%grid_mapping='sinh'
    call build_sinh_axis(n,3,0._dp,1._dp,1.5_dp,.true.,axis)
    allocate(q(-2:n+3,-2:4,-2:4,5),face(0:n,5),rhs(1:n,5))
    q=0
  end subroutine
  subroutine evaluate()
    do l=-2,n+3
      if(l>=1.and.l<=n) cycle
      q(l,1,1,:)=q(1+modulo(l-1,n),1,1,:)
    end do
    do l=0,n
      call compute_keep_face_flux(q,l,1,1,1,order,sim,nse,1,1,face(l,:))
    end do
    do l=1,n
      rhs(l,:)=-(face(l,:)-face(l-1,:))/axis%width(l)
      if(order==6) rhs(l,:)=-(face(l,:)-face(l-1,:))/axis%keep6_metric(l)
    end do
  end subroutine
end program
