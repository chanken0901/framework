program test_cuda_viscous_compare
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_convective_scheme, only : compute_convective_flux
  use mod_viscous_scheme, only : add_viscous_rhs, viscous_dt_limit
  use mod_nse_fluctuating, only : add_fh_transport, add_fh_noise
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_compute_dt, &
    nse_gpu_advance_ssprk3, nse_gpu_synchronize, nse_gpu_finalize
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: gpu
  real(dp), allocatable :: q_cpu(:,:,:,:), q_gpu(:,:,:,:)
  real(dp), allocatable :: q0(:,:,:,:), rhs(:,:,:,:), fface(:,:,:,:)
  real(dp) :: pi, x, y, z, rho, u, v, w, pressure
  real(dp) :: gpu_dt, cpu_dt, field_error
  integer :: i, j, k, g, iteration, iterations
  character(len=32) :: mode

  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 9
  sim%ny = 8
  sim%nz = 7
  sim%nghost = 3
  pi = acos(-1.0_dp)
  sim%dx = 2.0_dp*pi/real(sim%nx,dp)
  sim%dy = 2.0_dp*pi/real(sim%ny,dp)
  sim%dz = 2.0_dp*pi/real(sim%nz,dp)
  nse%cfl = 0.37_dp
  nse%viscous_scheme = 'central6'
  nse%reynolds = 25.0_dp
  nse%prandtl = 0.72_dp
  iterations=1
  if(command_argument_count()>0) then
    call get_command_argument(1,mode)
    if(trim(mode)=='fluctuating' .or. trim(mode)=='fh_zero') then
      nse%fh_enabled=.true.; nse%fh_boltzmann_number=1.e-4_dp
      if(trim(mode)=='fh_zero') nse%fh_boltzmann_number=0
      iterations=4
    end if
  end if

  g = sim%nghost
  allocate(q_cpu(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,nse%nv))
  allocate(q_gpu(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,nse%nv))
  allocate(q0(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,nse%nv))
  allocate(rhs(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,nse%nv))
  allocate(fface(0:sim%nx,0:sim%ny,0:sim%nz,nse%nv))
  q_cpu = 0.0_dp
  do k = 1, sim%nz
    z = 2.0_dp*pi*real(k-1,dp)/real(sim%nz,dp)
    do j = 1, sim%ny
      y = 2.0_dp*pi*real(j-1,dp)/real(sim%ny,dp)
      do i = 1, sim%nx
        x = 2.0_dp*pi*real(i-1,dp)/real(sim%nx,dp)
        rho = 1.0_dp + 0.03_dp*sin(x+2.0_dp*y-z)
        u = 0.12_dp*cos(x+2.0_dp*y)
        v = -0.08_dp*sin(2.0_dp*x-z)
        w = 0.05_dp*cos(y+z)
        pressure = 1.0_dp + 0.04_dp*cos(x-y+2.0_dp*z)
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
  q_gpu = q_cpu

  call nse_gpu_initialize(gpu, sim, nse)
  call nse_gpu_upload(gpu, q_gpu)
  call nse_gpu_compute_dt(gpu, gpu_dt)
  call reference_dt(q_cpu, cpu_dt, sim, nse)
  if (abs(gpu_dt-cpu_dt) > 2.0e-13_dp) then
    write(*,'(A,3ES24.16)') 'viscous CFL comparison failed: ', &
      gpu_dt, cpu_dt, abs(gpu_dt-cpu_dt)
    error stop 'CUDA viscous CFL differs from CPU reference'
  end if

  do iteration=1,iterations
    call reference_ssprk3(q_cpu, q0, rhs, fface, 1.0e-5_dp, sim, nse)
    call nse_gpu_advance_ssprk3(gpu, 1.0e-5_dp)
    sim%step=sim%step+1
  end do
  call nse_gpu_synchronize(gpu)
  call nse_gpu_download(gpu, q_gpu)
  if (abs(sum(q_gpu(1:sim%nx,1:sim%ny,1:sim%nz,5)) &
      -sum(q0(1:sim%nx,1:sim%ny,1:sim%nz,5)))>3.0e-12_dp) &
    error stop 'periodic CUDA viscous step lost total energy'
  field_error = maxval(abs(q_gpu(1:sim%nx,1:sim%ny,1:sim%nz,:) - &
    q_cpu(1:sim%nx,1:sim%ny,1:sim%nz,:)))
  if (field_error > 3.0e-12_dp) then
    write(*,'(A,ES24.16)') 'CUDA central6 field error: ', field_error
    error stop 'CUDA central6 step differs from CPU reference'
  end if

  call nse_gpu_finalize(gpu)
  if(nse%fh_enabled) write(*,'(A,ES16.8)') 'Resident CUDA LLNS comparison passed; max error=',field_error
  deallocate(q_cpu, q_gpu, q0, rhs, fface)
  write(*,'(A,ES16.8)') &
    'CUDA central6 comparison passed; max error = ', field_error

contains

  subroutine apply_periodic(q, config)
    type(simulation_config), intent(in) :: config
    real(dp), intent(inout) :: q(1-config%nghost:,1-config%nghost:, &
      1-config%nghost:,:)
    integer :: ii, jj, kk, wi, wj, wk

    do kk = 1-config%nghost, config%nz+config%nghost
      wk = 1 + modulo(kk-1, config%nz)
      do jj = 1-config%nghost, config%ny+config%nghost
        wj = 1 + modulo(jj-1, config%ny)
        do ii = 1-config%nghost, config%nx+config%nghost
          if (ii >= 1 .and. ii <= config%nx .and. &
              jj >= 1 .and. jj <= config%ny .and. &
              kk >= 1 .and. kk <= config%nz) cycle
          wi = 1 + modulo(ii-1, config%nx)
          q(ii,jj,kk,:) = q(wi,wj,wk,:)
        end do
      end do
    end do
  end subroutine apply_periodic

  subroutine reference_dt(q, dt, config, fluid)
    type(simulation_config), intent(in) :: config
    type(nse_config), intent(in) :: fluid
    real(dp), intent(in) :: q(1-config%nghost:,1-config%nghost:, &
      1-config%nghost:,:)
    real(dp), intent(out) :: dt
    real(dp) :: local_rho, local_u, local_v, local_w
    real(dp) :: local_p, sound_speed, max_speed, diffusion_dt
    integer :: ii, jj, kk

    max_speed = 0.0_dp
    do kk = 1, config%nz
      do jj = 1, config%ny
        do ii = 1, config%nx
          local_rho = max(q(ii,jj,kk,1), fluid%small_rho)
          local_u = q(ii,jj,kk,2)/local_rho
          local_v = q(ii,jj,kk,3)/local_rho
          local_w = q(ii,jj,kk,4)/local_rho
          local_p = max((fluid%gamma-1.0_dp)*(q(ii,jj,kk,5) - &
            0.5_dp*local_rho*(local_u**2+local_v**2+local_w**2)), &
            fluid%small_p)
          sound_speed = sqrt(fluid%gamma*local_p/local_rho)
          max_speed = max(max_speed, abs(local_u)+sound_speed, &
            abs(local_v)+sound_speed, abs(local_w)+sound_speed)
        end do
      end do
    end do
    dt = fluid%cfl*min(config%dx,config%dy,config%dz)/max_speed
    call viscous_dt_limit(q, diffusion_dt, config, fluid, &
      1, config%ny, 1, config%nz)
    dt = min(dt, diffusion_dt)
  end subroutine reference_dt

  subroutine reference_ssprk3(q, q_initial, local_rhs, flux, dt, &
      config, fluid)
    type(simulation_config), intent(in) :: config
    type(nse_config), intent(in) :: fluid
    real(dp), intent(inout) :: q(1-config%nghost:,1-config%nghost:, &
      1-config%nghost:,:)
    real(dp), intent(inout) :: q_initial(1-config%nghost:, &
      1-config%nghost:,1-config%nghost:,:)
    real(dp), intent(inout) :: local_rhs(1-config%nghost:, &
      1-config%nghost:,1-config%nghost:,:)
    real(dp), intent(inout) :: flux(0:,0:,0:,:)
    real(dp), intent(in) :: dt

    call apply_periodic(q,config)
    q_initial = q
    call reference_rhs(q, local_rhs, flux, config, fluid)
    q(1:config%nx,1:config%ny,1:config%nz,:) = &
      q_initial(1:config%nx,1:config%ny,1:config%nz,:) + &
      dt*local_rhs(1:config%nx,1:config%ny,1:config%nz,:)

    call reference_rhs(q, local_rhs, flux, config, fluid)
    q(1:config%nx,1:config%ny,1:config%nz,:) = &
      0.75_dp*q_initial(1:config%nx,1:config%ny,1:config%nz,:) + &
      0.25_dp*(q(1:config%nx,1:config%ny,1:config%nz,:) + &
      dt*local_rhs(1:config%nx,1:config%ny,1:config%nz,:))

    call reference_rhs(q, local_rhs, flux, config, fluid)
    q(1:config%nx,1:config%ny,1:config%nz,:) = &
      (1.0_dp/3.0_dp)*q_initial(1:config%nx,1:config%ny,1:config%nz,:) + &
      (2.0_dp/3.0_dp)*(q(1:config%nx,1:config%ny,1:config%nz,:) + &
      dt*local_rhs(1:config%nx,1:config%ny,1:config%nz,:))
    if(fluid%fh_enabled) then
      local_rhs=0
      call add_fh_noise(q_initial,local_rhs,dt,config,fluid,1,config%ny,1,config%nz)
      q(1:config%nx,1:config%ny,1:config%nz,:)=q(1:config%nx,1:config%ny,1:config%nz,:) &
        +dt*local_rhs(1:config%nx,1:config%ny,1:config%nz,:)
    end if
  end subroutine reference_ssprk3

  subroutine reference_rhs(q, local_rhs, flux, config, fluid)
    type(simulation_config), intent(in) :: config
    type(nse_config), intent(in) :: fluid
    real(dp), intent(inout) :: q(1-config%nghost:,1-config%nghost:, &
      1-config%nghost:,:)
    real(dp), intent(out) :: local_rhs(1-config%nghost:, &
      1-config%nghost:,1-config%nghost:,:)
    real(dp), intent(out) :: flux(0:,0:,0:,:)
    integer :: ii, jj, kk

    call apply_periodic(q, config)
    local_rhs = 0.0_dp
    call compute_convective_flux(q, flux, 1, config, fluid, &
      1, config%ny, 1, config%nz)
    do kk = 1, config%nz
      do jj = 1, config%ny
        do ii = 1, config%nx
          local_rhs(ii,jj,kk,:) = local_rhs(ii,jj,kk,:) - &
            (flux(ii,jj,kk,:)-flux(ii-1,jj,kk,:))/config%dx
        end do
      end do
    end do

    call compute_convective_flux(q, flux, 2, config, fluid, &
      1, config%ny, 1, config%nz)
    do kk = 1, config%nz
      do jj = 1, config%ny
        do ii = 1, config%nx
          local_rhs(ii,jj,kk,:) = local_rhs(ii,jj,kk,:) - &
            (flux(ii,jj,kk,:)-flux(ii,jj-1,kk,:))/config%dy
        end do
      end do
    end do

    call compute_convective_flux(q, flux, 3, config, fluid, &
      1, config%ny, 1, config%nz)
    do kk = 1, config%nz
      do jj = 1, config%ny
        do ii = 1, config%nx
          local_rhs(ii,jj,kk,:) = local_rhs(ii,jj,kk,:) - &
            (flux(ii,jj,kk,:)-flux(ii,jj,kk-1,:))/config%dz
        end do
      end do
    end do
    if(fluid%fh_enabled) then
      call add_fh_transport(q,local_rhs,config,fluid,1,config%ny,1,config%nz)
    else
      call add_viscous_rhs(q, local_rhs, config, fluid, &
        1, config%ny, 1, config%nz)
    end if
  end subroutine reference_rhs

end program test_cuda_viscous_compare
