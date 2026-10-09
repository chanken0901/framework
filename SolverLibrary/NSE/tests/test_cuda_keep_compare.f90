program test_cuda_keep_compare
  use mod_precision, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_convective_scheme, only : compute_convective_flux
  use mod_grid_fvm, only: build_uniform_grid
  use mod_convective_keep, only: keep6_inverse_metric
  use mod_convective_hybrid, only: hybrid_face_weight
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_compute_dt, &
    nse_gpu_advance_ssprk3, nse_gpu_synchronize, nse_gpu_finalize
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: gpu
  real(dp), allocatable :: q_cpu(:,:,:,:), q_gpu(:,:,:,:)
  real(dp), allocatable :: q0(:,:,:,:), rhs(:,:,:,:), fface(:,:,:,:)
  real(dp) :: x, y, z, rho, u, v, w, pressure
  real(dp) :: gpu_dt, cpu_dt, field_error, comparison_tolerance
  integer :: i, j, k, scheme_code, step, steps
  real(dp) :: step_dt
  logical :: shock_test, mapped_test
  integer :: zero_weights,blend_weights,full_weights
  real(dp) :: alpha
  character(len=16) :: scheme_argument, profile_argument

  call init_simulation_config(sim)
  call init_nse_config(nse)
  scheme_code = 6
  call get_command_argument(1, scheme_argument)
  if (len_trim(scheme_argument) > 0) read(scheme_argument,*) scheme_code
  call get_command_argument(2, profile_argument)
  mapped_test = index(trim(profile_argument),'mapped')==1
  shock_test = trim(profile_argument) == 'shock'.or.mapped_test
  if (scheme_code /= 2 .and. scheme_code /= 5 .and. &
      scheme_code /= 6 .and. scheme_code /= 7) then
    error stop 'CUDA comparison scheme code must be 2, 5, 6, or 7'
  end if
  if (scheme_code == 2) then
    nse%convective_scheme = 'keep2'
  else if (scheme_code == 5) then
    nse%convective_scheme = 'weno5z_roe'
  else if (scheme_code == 7) then
    nse%convective_scheme = 'hybrid'
    nse%hybrid_smooth_scheme = 'keep6'
    nse%hybrid_shock_scheme = 'weno5z_roe'
    nse%hybrid_sensor_onset = 1.0e-4_dp
    nse%hybrid_sensor_full = 2.0e-3_dp
  else
    nse%convective_scheme = 'keep6'
  end if
  comparison_tolerance = 2.0e-12_dp
  if (scheme_code == 5 .or. scheme_code == 7) &
    comparison_tolerance = 2.0e-10_dp
  sim%nx = 7
  sim%ny = 6
  sim%nz = 5
  sim%nghost = 3
  sim%dx = 0.17_dp
  sim%dy = 0.19_dp
  sim%dz = 0.23_dp
  nse%cfl = 0.37_dp
  steps = 1
  step_dt = 2.5e-4_dp
  if (shock_test) then
    if (scheme_code /= 7) error stop 'shock regression requires hybrid'
    sim%nx = 128
    sim%ny = 16
    sim%nz = 16
    sim%dx = 0.03125_dp
    sim%dy = sim%dx
    sim%dz = sim%dx
    nse%hybrid_sensor_onset = 0.01_dp
    nse%hybrid_sensor_full = 0.1_dp
    step_dt = 0.002_dp
    steps = 5
  end if
  if(mapped_test) then
    sim%nx=32;sim%ny=12;sim%nz=12
    sim%grid_mapping='sinh';sim%mapped_keep6=.true.;sim%grid_stretch=[.4_dp,1._dp,1.5_dp]
    call build_uniform_grid(sim,1,sim%ny,1,sim%nz,[.true.,.true.,.true.])
    if(trim(profile_argument)=='mapped_blend') nse%hybrid_sensor_full=.9_dp
  end if

  allocate(q_cpu(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost,nse%nv))
  allocate(q_gpu(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost,nse%nv))
  allocate(q0(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost,nse%nv))
  allocate(rhs(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost,nse%nv))
  allocate(fface(0:sim%nx,0:sim%ny,0:sim%nz,nse%nv))
  q_cpu = 0.0_dp
  do k = 1, sim%nz
    z = real(k-1,dp) / real(sim%nz,dp)
    do j = 1, sim%ny
      y = real(j-1,dp) / real(sim%ny,dp)
      do i = 1, sim%nx
        x = real(i-1,dp) / real(sim%nx,dp)
        rho = 1.0_dp + 0.03_dp*sin(2.0_dp*x+3.0_dp*y-z)
        u = 0.12_dp*cos(x+2.0_dp*y)
        v = -0.08_dp*sin(2.0_dp*x-z)
        w = 0.05_dp*cos(y+z)
        pressure = 1.0_dp + 0.04_dp*cos(x-y+2.0_dp*z)
        if (shock_test) then
          rho = 1.0_dp
          u = 0.0_dp
          v = 0.0_dp
          w = 0.0_dp
          pressure = 1.0_dp/nse%gamma
          if (i <= sim%nx/2) pressure = 5.0_dp/nse%gamma
          if(mapped_test) u=merge(.1_dp,-.1_dp,i<=sim%nx/2)
        end if
        q_cpu(i,j,k,1) = rho
        q_cpu(i,j,k,2) = rho*u
        q_cpu(i,j,k,3) = rho*v
        q_cpu(i,j,k,4) = rho*w
        q_cpu(i,j,k,5) = pressure/(nse%gamma-1.0_dp) + &
          0.5_dp*rho*(u*u+v*v+w*w)
      end do
    end do
  end do
  call apply_periodic(q_cpu, sim)
  if(mapped_test) then
    zero_weights=0;blend_weights=0;full_weights=0
    do i=0,sim%nx
      alpha=hybrid_face_weight(q_cpu,i,1,1,1,sim,nse,1,1)
      if(alpha<=0) zero_weights=zero_weights+1
      if(alpha>0.and.alpha<1) blend_weights=blend_weights+1
      if(alpha>=1) full_weights=full_weights+1
    end do
    print *, 'Hybrid zero/blended/full face weights:',zero_weights,blend_weights,full_weights
    if(zero_weights==0) error stop 'Mapped hybrid smooth branch not exercised'
    if(trim(profile_argument)=='mapped_blend') then
      if(blend_weights==0) error stop 'Mapped hybrid blend branch not exercised'
    else
      if(full_weights==0) error stop 'Mapped hybrid WENO branch not exercised'
    end if
  end if
  q_gpu = q_cpu

  call nse_gpu_initialize(gpu, sim, nse)
  call nse_gpu_upload(gpu, q_gpu)
  call nse_gpu_compute_dt(gpu, gpu_dt)
  call reference_dt(q_cpu, cpu_dt, sim, nse)
  if (abs(gpu_dt-cpu_dt) > 2.0e-13_dp) then
    write(*,'(A,3ES24.16)') "CFL comparison failed: ", &
      gpu_dt, cpu_dt, abs(gpu_dt-cpu_dt)
    error stop "CUDA CFL differs from CPU reference"
  end if

  do step = 1, steps
    call reference_ssprk3(q_cpu, q0, rhs, fface, step_dt, sim, nse)
    call nse_gpu_advance_ssprk3(gpu, step_dt)
    call nse_gpu_synchronize(gpu)
    call nse_gpu_download(gpu, q_gpu)
    if (.not. all(ieee_is_finite(q_gpu))) error stop 'nonfinite CUDA state'
    if (.not. all(ieee_is_finite(q_cpu))) error stop 'nonfinite CPU reference'
    field_error = maxval(abs(q_gpu(1:sim%nx,1:sim%ny,1:sim%nz,:) - &
      q_cpu(1:sim%nx,1:sim%ny,1:sim%nz,:)))
    if (field_error > comparison_tolerance) then
      write(*,'(A,ES24.16)') "Convective/SSPRK3 field error: ", field_error
      error stop "CUDA convective step differs from CPU reference"
    end if
  end do

  call nse_gpu_finalize(gpu)
  deallocate(q_cpu, q_gpu, q0, rhs, fface)
  write(*,'(A,I0,A,ES16.8)') &
    "CUDA convective comparison passed; scheme code = ", scheme_code, &
    ", max error = ", field_error

contains

  subroutine apply_periodic(q, sim)
    type(simulation_config), intent(in) :: sim
    real(dp), intent(inout) :: q(1-sim%nghost:,1-sim%nghost:,1-sim%nghost:,:)
    integer :: i, j, k, wrapped_i, wrapped_j, wrapped_k

    do k = 1-sim%nghost, sim%nz+sim%nghost
      wrapped_k = 1 + modulo(k-1, sim%nz)
      do j = 1-sim%nghost, sim%ny+sim%nghost
        wrapped_j = 1 + modulo(j-1, sim%ny)
        do i = 1-sim%nghost, sim%nx+sim%nghost
          if (i >= 1 .and. i <= sim%nx .and. &
              j >= 1 .and. j <= sim%ny .and. &
              k >= 1 .and. k <= sim%nz) cycle
          wrapped_i = 1 + modulo(i-1, sim%nx)
          q(i,j,k,:) = q(wrapped_i,wrapped_j,wrapped_k,:)
        end do
      end do
    end do
  end subroutine apply_periodic

  subroutine reference_dt(q, dt, sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    real(dp), intent(in) :: q(1-sim%nghost:,1-sim%nghost:,1-sim%nghost:,:)
    real(dp), intent(out) :: dt
    real(dp) :: rho, u, v, w, pressure, sound_speed, max_speed
    integer :: i, j, k

    max_speed = 0.0_dp
    do k = 1, sim%nz
      do j = 1, sim%ny
        do i = 1, sim%nx
          rho = max(q(i,j,k,1), nse%small_rho)
          u = q(i,j,k,2)/rho
          v = q(i,j,k,3)/rho
          w = q(i,j,k,4)/rho
          pressure = max((nse%gamma-1.0_dp)*(q(i,j,k,5) - &
            0.5_dp*rho*(u*u+v*v+w*w)), nse%small_p)
          sound_speed = sqrt(nse%gamma*pressure/rho)
          max_speed = max(max_speed, abs(u)+sound_speed, &
            abs(v)+sound_speed, abs(w)+sound_speed)
        end do
      end do
    end do
    dt = nse%cfl*min(sim%dx,sim%dy,sim%dz)/max_speed
    if(sim%mapped_keep6) dt=nse%cfl/(max_speed*(1/sim%dx+1/sim%dy+1/sim%dz))
  end subroutine reference_dt

  subroutine reference_ssprk3(q, q0, rhs, fface, dt, sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    real(dp), intent(inout) :: q(1-sim%nghost:,1-sim%nghost:,1-sim%nghost:,:)
    real(dp), intent(inout) :: q0(1-sim%nghost:,1-sim%nghost:, &
      1-sim%nghost:,:), rhs(1-sim%nghost:,1-sim%nghost:, &
      1-sim%nghost:,:)
    real(dp), intent(inout) :: fface(0:,0:,0:,:)
    real(dp), intent(in) :: dt

    q0 = q
    call reference_rhs(q, rhs, fface, sim, nse)
    q(1:sim%nx,1:sim%ny,1:sim%nz,:) = &
      q0(1:sim%nx,1:sim%ny,1:sim%nz,:) + &
      dt*rhs(1:sim%nx,1:sim%ny,1:sim%nz,:)

    call reference_rhs(q, rhs, fface, sim, nse)
    q(1:sim%nx,1:sim%ny,1:sim%nz,:) = &
      0.75_dp*q0(1:sim%nx,1:sim%ny,1:sim%nz,:) + 0.25_dp*( &
      q(1:sim%nx,1:sim%ny,1:sim%nz,:) + &
      dt*rhs(1:sim%nx,1:sim%ny,1:sim%nz,:))

    call reference_rhs(q, rhs, fface, sim, nse)
    q(1:sim%nx,1:sim%ny,1:sim%nz,:) = &
      (1.0_dp/3.0_dp)*q0(1:sim%nx,1:sim%ny,1:sim%nz,:) + &
      (2.0_dp/3.0_dp)*(q(1:sim%nx,1:sim%ny,1:sim%nz,:) + &
      dt*rhs(1:sim%nx,1:sim%ny,1:sim%nz,:))
  end subroutine reference_ssprk3

  subroutine reference_rhs(q, rhs, fface, sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    real(dp), intent(inout) :: q(1-sim%nghost:,1-sim%nghost:,1-sim%nghost:,:)
    real(dp), intent(out) :: rhs(1-sim%nghost:,1-sim%nghost:, &
      1-sim%nghost:,:), fface(0:,0:,0:,:)
    integer :: i, j, k

    call apply_periodic(q, sim)
    rhs = 0.0_dp
    call compute_convective_flux(q, fface, 1, sim, nse, &
      1, sim%ny, 1, sim%nz)
    do k = 1, sim%nz
      do j = 1, sim%ny
        do i = 1, sim%nx
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (fface(i,j,k,:)-fface(i-1,j,k,:))/sim%dx
          if(sim%mapped_keep6) rhs(i,j,k,:)= &
            -(fface(i,j,k,:)-fface(i-1,j,k,:))*keep6_inverse_metric(i,j,k,1)
        end do
      end do
    end do

    call compute_convective_flux(q, fface, 2, sim, nse, &
      1, sim%ny, 1, sim%nz)
    do k = 1, sim%nz
      do j = 1, sim%ny
        do i = 1, sim%nx
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (fface(i,j,k,:)-fface(i,j-1,k,:))/sim%dy
          if(sim%mapped_keep6) rhs(i,j,k,:)=rhs(i,j,k,:)+ &
            (fface(i,j,k,:)-fface(i,j-1,k,:))*(1._dp/sim%dy-keep6_inverse_metric(i,j,k,2))
        end do
      end do
    end do

    call compute_convective_flux(q, fface, 3, sim, nse, &
      1, sim%ny, 1, sim%nz)
    do k = 1, sim%nz
      do j = 1, sim%ny
        do i = 1, sim%nx
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (fface(i,j,k,:)-fface(i,j,k-1,:))/sim%dz
          if(sim%mapped_keep6) rhs(i,j,k,:)=rhs(i,j,k,:)+ &
            (fface(i,j,k,:)-fface(i,j,k-1,:))*(1._dp/sim%dz-keep6_inverse_metric(i,j,k,3))
        end do
      end do
    end do
  end subroutine reference_rhs

end program test_cuda_keep_compare
