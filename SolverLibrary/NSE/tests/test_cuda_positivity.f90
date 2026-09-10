program test_cuda_positivity
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, nse_gpu_upload, &
    nse_gpu_download, nse_gpu_advance_ssprk3, nse_gpu_finalize
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: gpu
  real(dp), allocatable :: q(:,:,:,:), initial(:,:,:,:), accepted(:,:,:,:)
  real(dp) :: dt, energy0, energy1, p, u, pi
  integer :: i,j,k
  character(len=32) :: mode
  call init_simulation_config(sim)
  call init_nse_config(nse)
  mode = ''
  if (command_argument_count() > 0) call get_command_argument(1,mode)
  sim%nx=8; sim%ny=16; sim%nz=8; sim%nghost=3
  pi=acos(-1.0_dp)
  sim%dx=2*pi/sim%nx; sim%dy=2*pi/sim%ny; sim%dz=2*pi/sim%nz
  nse%convective_scheme='weno5z_roe'
  nse%viscous_scheme='none'
  nse%forcing_scheme='petersen_livescu'
  nse%forcing_spectrum='low_wavenumber'
  nse%forcing_fft_backend='cufft'
  nse%forcing_k_cutoff=2.5_dp
  nse%forcing_target_dissipation=0.1_dp
  nse%forcing_dilatational_ratio=0.0_dp
  nse%forcing_report_interval=0
  allocate(q(-2:11,-2:19,-2:11,5),initial(-2:11,-2:19,-2:11,5), &
    accepted(-2:11,-2:19,-2:11,5))
  q=0.0_dp
  do k=1,sim%nz
    do j=1,sim%ny
      u=0.5_dp*sin(2*pi*(j-1)/sim%ny)
      do i=1,sim%nx
        q(i,j,k,:)=[1.0_dp,u,0.0_dp,0.0_dp,0.01_dp/(nse%gamma-1)+0.5_dp*u*u]
      end do
    end do
  end do
  initial=q
  energy0=sum(q(1:sim%nx,1:sim%ny,1:sim%nz,5))
  call nse_gpu_initialize(gpu,sim,nse)
  call nse_gpu_upload(gpu,q)
  if (mode == 'fixed-reject') then
    call nse_gpu_advance_ssprk3(gpu,0.5_dp)
    ! CTest expects nonzero only when the backend actually rejects this step.
    stop 0
  end if
  call nse_gpu_advance_ssprk3(gpu,0.5_dp,dt)
  if (dt <= 0.0_dp .or. dt >= 0.5_dp) error stop 'forcing budget did not reduce dt'
  call nse_gpu_download(gpu,q)
  if (.not. all(ieee_is_finite(q))) error stop 'nonfinite accepted state'
  do k=1,sim%nz
    do j=1,sim%ny
      do i=1,sim%nx
        p=(nse%gamma-1)*(q(i,j,k,5)-0.5_dp*sum(q(i,j,k,2:4)**2)/q(i,j,k,1))
        if (p < nse%small_p .or. q(i,j,k,1)<nse%small_rho) error stop 'negative accepted state'
      end do
    end do
  end do
  energy1=sum(q(1:sim%nx,1:sim%ny,1:sim%nz,5))
  if (abs(energy1-energy0)>1.e-11_dp) error stop 'retry changed conserved energy'
  accepted=q
  call nse_gpu_upload(gpu,initial)
  call nse_gpu_advance_ssprk3(gpu,dt)
  call nse_gpu_download(gpu,q)
  if (maxval(abs(q-accepted))>1.e-12_dp) error stop 'retry differs from clean accepted step'
  call nse_gpu_finalize(gpu)
  ! Strong double-rarefaction interface exercises the robust shared face flux.
  nse%forcing_scheme='none'
  do k=1,sim%nz
    do j=1,sim%ny
      do i=1,sim%nx
        u=merge(-2.0_dp,2.0_dp,i<=sim%nx/2)
        q(i,j,k,:)=[1.0_dp,u,0.0_dp,0.0_dp,0.001_dp/(nse%gamma-1)+0.5_dp*u*u]
      end do
    end do
  end do
  energy0=sum(q(1:sim%nx,1:sim%ny,1:sim%nz,5))
  call nse_gpu_initialize(gpu,sim,nse)
  call nse_gpu_upload(gpu,q)
  call nse_gpu_advance_ssprk3(gpu,0.5_dp,dt)
  call nse_gpu_download(gpu,q)
  if (.not. all(ieee_is_finite(q))) error stop 'nonfinite rarefaction state'
  do k=1,sim%nz
    do j=1,sim%ny
      do i=1,sim%nx
        p=(nse%gamma-1)*(q(i,j,k,5)-0.5_dp*sum(q(i,j,k,2:4)**2)/q(i,j,k,1))
        if (p<nse%small_p .or. q(i,j,k,1)<nse%small_rho) error stop 'negative rarefaction state'
      end do
    end do
  end do
  energy1=sum(q(1:sim%nx,1:sim%ny,1:sim%nz,5))
  if (abs(energy1-energy0)>1.e-9_dp) error stop 'rarefaction flux lost conservation'
  call nse_gpu_finalize(gpu)
  write(*,'(A,ES16.8)') 'CUDA positivity test passed; accepted dt=',dt
end program
