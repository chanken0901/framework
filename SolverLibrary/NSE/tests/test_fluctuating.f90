program test_fluctuating
  use iso_fortran_env, only: int64
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config, update_derived_config
  use mod_model_config, only: nse_config
  use mod_fh_random, only: philox4x32
  use mod_nse_fluctuating, only: fh_sample,add_fh_noise,add_fh_transport,validate_fh
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp) :: q(-2:7,-2:7,-2:7,5),rhs(-2:7,-2:7,-2:7,5),saved(-2:7,-2:7,-2:7,5)
  real(dp) :: local(-2:7,0:6,0:6,5),local_rhs(-2:7,0:6,0:6,5)
  real(dp) :: s(3,3),heat(3),v(4),mean(4),cov(4,4),mu,kappa,a,x,dt,pi,velocity(3)
  real(dp) :: work,variance,projection,w(4,4,4,3),expected,average
  real(dp) :: thermal_probe(4,4,4),thermal_work,thermal_variance,thermal_average,thermal_projection
  integer :: i,j,k,c,d,trial
  integer(int64) :: words(4)
  character(len=32) :: failure
  words=philox4x32([0_int64,0_int64,0_int64,0_int64],[0_int64,0_int64])
  if (any(words/=[int(z'6627e8d5',int64),int(z'e169c58d',int64), &
      int(z'bc57ac4c',int64),int(z'9b00dbd8',int64)])) error stop 'Philox known-answer failed'
  sim%nx=4; sim%ny=4; sim%nz=4; sim%nghost=3
  call update_derived_config(sim)
  nse%fh_enabled=.true.; nse%fh_boltzmann_number=1.e-5_dp
  nse%reynolds=100; nse%viscous_scheme='central6'; dt=0.01_dp
  if (command_argument_count()>0) then
    call get_command_argument(1,failure)
    select case(trim(failure))
    case('adaptive'); sim%use_fixed_dt=.false.
    case('boundary'); nse%boundary_face_type(1)='reflective'
    case('inviscid'); nse%reynolds=0
    case('negative_beta'); nse%fh_boltzmann_number=-1
    case('seed'); nse%fh_seed=-1
    case('weno'); nse%convective_scheme='weno5z_roe'
    case default; stop 0
    end select
    call validate_fh(sim,nse)
    stop 0  ! CTest WILL_FAIL detects a missing rejection.
  end if
  call validate_fh(sim,nse)
  mu=1/nse%reynolds; kappa=mu*nse%gamma/((nse%gamma-1)*nse%prandtl)
  a=2*nse%fh_boltzmann_number/(sim%dx*sim%dy*sim%dz*dt)
  mean=0; cov=0
  do trial=1,30000
    sim%step=trial
    call fh_sample(1,1,1,sim,nse,1.0_dp,dt,s,heat)
    if (abs(s(1,1)+s(2,2)+s(3,3))>1.e-12_dp) error stop 'stress trace'
    if (maxval(abs(s-transpose(s)))>0) error stop 'stress symmetry'
    v=[s(1,1)/sqrt(a*mu),s(2,2)/sqrt(a*mu),s(1,2)/sqrt(a*mu),heat(1)/sqrt(a*kappa)]
    mean=mean+v
    do c=1,4
      do d=1,4
        cov(c,d)=cov(c,d)+v(c)*v(d)
      end do
    end do
  end do
  mean=mean/30000; cov=cov/30000
  if (maxval(abs(mean))>0.03_dp) error stop 'noise mean'
  if (abs(cov(1,1)-4.0_dp/3)>0.05_dp .or. abs(cov(1,2)+2.0_dp/3)>0.05_dp &
      .or. abs(cov(3,3)-1)>0.05_dp .or. abs(cov(4,4)-1)>0.05_dp &
      .or. maxval(abs(cov(1:3,4)))>0.04_dp) error stop 'LL stress/heat covariance'
  q=0; q(:,:,:,1)=1; q(:,:,:,5)=1/(nse%gamma-1)
  rhs=0
  !$OMP PARALLEL
  call add_fh_noise(q,rhs,dt,sim,nse,1,4,1,4)
  !$OMP END PARALLEL
  saved=rhs
  do c=1,5
    if (abs(sum(rhs(1:4,1:4,1:4,c)))>1.e-11_dp) error stop 'periodic conservation'
  end do
  local=q(:,0:6,0:6,:); local_rhs=0
  call add_fh_noise(local,local_rhs,dt,sim,nse,3,3,3,3)
  if (maxval(abs(local_rhs(:,3:3,3:3,:)-rhs(:,3:3,3:3,:)))>1.e-14_dp) &
    error stop 'global index decomposition reproducibility'
  rhs=0
  call add_fh_noise(q,rhs,4*dt,sim,nse,1,4,1,4)
  if (maxval(abs(rhs-0.5_dp*saved))>1.e-13_dp) error stop 'dt scaling'
  nse%fh_enabled=.false.; rhs=0
  call add_fh_noise(q,rhs,dt,sim,nse,1,4,1,4)
  if (maxval(abs(rhs))>0) error stop 'disabled noise changed RHS'
  nse%fh_enabled=.true.
  ! Arbitrary non-axis-aligned velocity probe: test w^T B B^T w = -2 beta T w^T L w / V.
  pi=acos(-1.0_dp)
  do k=-2,7
    do j=-2,7
      do i=-2,7
        x=2*pi*real(i+j+k,dp)/4
        velocity=[sin(x),cos(x),sin(x+0.7_dp)]
        q(i,j,k,2:4)=velocity
        q(i,j,k,5)=1/(nse%gamma-1)+0.5_dp*sum(velocity**2)
      end do
    end do
  end do
  w=q(1:4,1:4,1:4,2:4); rhs=0
  call add_fh_transport(q,rhs,sim,nse,1,4,1,4)
  work=sum(w*rhs(1:4,1:4,1:4,2:4))
  if (work>=0) error stop 'transport must dissipate kinetic energy'
  expected=-a*work
  do k=-2,7
    do j=-2,7
      do i=-2,7
        q(i,j,k,2:4)=0
        q(i,j,k,5)=(1+0.1_dp*sin(2*pi*real(i+j+k,dp)/4))/(nse%gamma-1)
      end do
    end do
  end do
  thermal_probe=((nse%gamma-1)*q(1:4,1:4,1:4,5)-1)/0.1_dp
  rhs=0
  call add_fh_transport(q,rhs,sim,nse,1,4,1,4)
  thermal_work=sum(thermal_probe*rhs(1:4,1:4,1:4,5))/0.1_dp
  if (thermal_work>=0) error stop 'heat conduction must dissipate temperature perturbations'
  q(:,:,:,2:4)=0; q(:,:,:,5)=1/(nse%gamma-1)
  variance=0; average=0; thermal_variance=0; thermal_average=0
  do trial=1,1500
    sim%step=trial; rhs=0
    call add_fh_noise(q,rhs,dt,sim,nse,1,4,1,4)
    projection=sum(w*rhs(1:4,1:4,1:4,2:4))
    variance=variance+projection**2; average=average+projection
    thermal_projection=sum(thermal_probe*rhs(1:4,1:4,1:4,5))
    thermal_variance=thermal_variance+thermal_projection**2
    thermal_average=thermal_average+thermal_projection
  end do
  variance=variance/1500-(average/1500)**2
  if (abs(variance/expected-1)>0.12_dp) error stop 'discrete fluctuation dissipation balance'
  thermal_variance=thermal_variance/1500-(thermal_average/1500)**2
  if (abs(thermal_variance/(-a*thermal_work)-1)>0.12_dp) error stop 'discrete thermal FDT'
  print *, 'Landau-Lifshitz covariance, conservation, decomposition and discrete FDT tests passed'
end program
