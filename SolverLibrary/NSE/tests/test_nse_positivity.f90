program test_nse_positivity
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_grid_fvm, only : build_uniform_grid
  use mod_nse_time_integration, only : advance_nse_ssprk3
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:),q0(:,:,:,:),rhs(:,:,:,:),face(:,:,:,:),initial(:,:,:,:),accepted(:,:,:,:)
  real(dp) :: dt,u,p,total0(5),total1(5)
  integer :: i,j,k,v
  character(len=32) :: mode
  call init_simulation_config(sim)
  call init_nse_config(nse)
  mode=''
  if (command_argument_count()>0) call get_command_argument(1,mode)
  sim%nx=16; sim%ny=8; sim%nz=8; sim%nghost=3
  sim%use_fixed_dt=(mode=='fixed-reject')
  sim%use_mpi=.false.
  nse%convective_scheme='weno5z_roe'; nse%viscous_scheme='none'
  call build_uniform_grid(sim,1,8,1,8)
  allocate(q(-2:19,-2:11,-2:11,5),q0(-2:19,-2:11,-2:11,5), &
    rhs(-2:19,-2:11,-2:11,5),face(0:16,0:8,0:8,5), &
    initial(-2:19,-2:11,-2:11,5),accepted(-2:19,-2:11,-2:11,5))
  q=0; rhs=0
  do k=1,8
    do j=1,8
      do i=1,16
        u=merge(-2.0_dp,2.0_dp,i<=8)
        q(i,j,k,:)=[1.0_dp,u,0.0_dp,0.0_dp,0.001_dp/(nse%gamma-1)+0.5_dp*u*u]
      end do
    end do
  end do
  initial=q
  do v=1,5
    total0(v)=sum(q(1:16,1:8,1:8,v))
  end do
  dt=10.0_dp
  !$OMP PARALLEL DEFAULT(NONE) SHARED(q,q0,rhs,face,dt,sim,nse)
  call advance_nse_ssprk3(q,q0,rhs,face,dt,sim,nse,1,8,1,8)
  !$OMP END PARALLEL
  if (mode=='fixed-reject') stop 0
  if (dt>=10.0_dp .or. dt<=0) error stop 'CPU budget did not reduce dt'
  do k=1,8
    do j=1,8
      do i=1,16
        p=(nse%gamma-1)*(q(i,j,k,5)-0.5_dp*sum(q(i,j,k,2:4)**2)/q(i,j,k,1))
        if (p<nse%small_p .or. q(i,j,k,1)<nse%small_rho) error stop 'CPU negative state'
      end do
    end do
  end do
  do v=1,5
    total1(v)=sum(q(1:16,1:8,1:8,v))
  end do
  if (maxval(abs(total1-total0))>1.e-9_dp) error stop 'CPU retry lost conservation'
  accepted=q
  q=initial
  sim%use_fixed_dt=.true.
  !$OMP PARALLEL DEFAULT(NONE) SHARED(q,q0,rhs,face,dt,sim,nse)
  call advance_nse_ssprk3(q,q0,rhs,face,dt,sim,nse,1,8,1,8)
  !$OMP END PARALLEL
  if (maxval(abs(q-accepted))>1.e-12_dp) error stop 'CPU rollback differs from clean step'
  write(*,'(A,ES16.8)') 'CPU positivity test passed; accepted dt=',dt
end program
