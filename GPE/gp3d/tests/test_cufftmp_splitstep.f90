program test_cufftmp_splitstep
  use gp3d_types, only: dp, pi, gp3d_grid_t, gp3d_params_t, gp3d_state_t
  use gp3d_grid, only: gp3d_grid_init
  use gp3d_fft, only: gp3d_fft_plan_t, gp3d_fft_init, gp3d_fft_finalize
  use gp3d_solver, only: gp3d_state_allocate, gp3d_normalize, &
    gp3d_step_split_operator, gp3d_density_norm, gp3d_energy
  use gp3d_gpu, only: gp3d_gpu_context_t, gp3d_gpu_init, gp3d_gpu_upload, &
    gp3d_gpu_download, gp3d_gpu_step, gp3d_gpu_diagnostics, gp3d_gpu_finalize
  use gp3d_mpi, only: gp3d_mpi_t, gp3d_mpi_init, gp3d_mpi_finalize, &
    gp3d_mpi_is_root, gp3d_mpi_max_real
  implicit none

  real(dp), parameter :: field_tolerance = 2.0e-9_dp
  real(dp), parameter :: scalar_tolerance = 2.0e-9_dp
  type(gp3d_mpi_t) :: mpi
  type(gp3d_grid_t) :: grid
  type(gp3d_params_t) :: params
  type(gp3d_state_t) :: cpu_state, gpu_state
  type(gp3d_fft_plan_t) :: fft_plan
  type(gp3d_gpu_context_t) :: gpu
  real(dp) :: local_error, field_error
  real(dp) :: cpu_norm, gpu_norm, cpu_energy, gpu_energy
  integer :: i, j, k, kg

  call gp3d_mpi_init(mpi)
  if (mpi%nprocs > 4) error stop "cuFFTMp test supports at most four ranks"
  call gp3d_grid_init(grid, 4, 4, 4, 2.0_dp * pi, 2.0_dp * pi, &
    2.0_dp * pi, mpi%rank, mpi%nprocs)
  call gp3d_state_allocate(cpu_state, grid)
  do k = 1, grid%local_nz
    kg = grid%k_start + k - 1
    do j = 1, grid%ny
      do i = 1, grid%nx
        cpu_state%psi(i,j,k) = cmplx( &
          1.0_dp + 0.1_dp * cos(grid%x(i) + 2.0_dp * grid%y(j)), &
          0.05_dp * sin(grid%y(j) - grid%z(kg)), kind=dp)
        cpu_state%potential(i,j,k) = 0.02_dp * ( &
          cos(grid%x(i)) + sin(grid%y(j)) + cos(grid%z(kg)))
      end do
    end do
  end do

  params%dt = 2.5e-4_dp
  params%hbar = 1.0_dp
  params%mass = 0.75_dp
  params%g = 1.3_dp
  params%norm = 1.0_dp
  params%imaginary_time = .false.
  call gp3d_normalize(cpu_state, grid, params%norm, mpi)
  gpu_state = cpu_state

  call gp3d_fft_init(fft_plan, grid%nx, grid%ny, grid%nz, &
    mpi%comm, mpi%rank, mpi%nprocs)
  call gp3d_gpu_init(gpu, grid, mpi, need_argle=.false.)
  call gp3d_gpu_upload(gpu, gpu_state, grid)

  call gp3d_step_split_operator(cpu_state, grid, params, fft_plan, mpi)
  call gp3d_gpu_step(gpu, grid, params)
  call gp3d_gpu_download(gpu, gpu_state)
  local_error = maxval(abs(cpu_state%psi - gpu_state%psi))
  call gp3d_mpi_max_real(mpi, local_error, field_error)
  if (field_error > field_tolerance) then
    if (gp3d_mpi_is_root(mpi)) then
      write(*,'(a,es16.8)') "real-time field error: ", field_error
    end if
    error stop "cuFFTMp real-time split step differs from CPU reference"
  end if

  cpu_norm = gp3d_density_norm(cpu_state, grid, mpi)
  cpu_energy = gp3d_energy(cpu_state, grid, params, fft_plan, mpi)
  call gp3d_gpu_diagnostics(gpu, grid, params, gpu_norm, gpu_energy)
  if (abs(cpu_norm - gpu_norm) > scalar_tolerance) then
    error stop "cuFFTMp norm differs from CPU reference"
  end if
  if (abs(cpu_energy - gpu_energy) > scalar_tolerance) then
    if (gp3d_mpi_is_root(mpi)) then
      write(*,'(a,es16.8)') "energy error: ", abs(cpu_energy - gpu_energy)
    end if
    error stop "cuFFTMp energy differs from CPU reference"
  end if

  params%imaginary_time = .true.
  call gp3d_step_split_operator(cpu_state, grid, params, fft_plan, mpi)
  call gp3d_gpu_step(gpu, grid, params)
  call gp3d_gpu_download(gpu, gpu_state)
  local_error = maxval(abs(cpu_state%psi - gpu_state%psi))
  call gp3d_mpi_max_real(mpi, local_error, field_error)
  if (field_error > field_tolerance) then
    if (gp3d_mpi_is_root(mpi)) then
      write(*,'(a,es16.8)') "imaginary-time field error: ", field_error
    end if
    error stop "cuFFTMp imaginary-time split step differs from CPU reference"
  end if

  call gp3d_gpu_finalize(gpu)
  call gp3d_fft_finalize(fft_plan)
  if (gp3d_mpi_is_root(mpi)) then
    write(*,'(a,es16.8)') "cuFFTMp split-step comparison passed; max error = ", &
      field_error
  end if
  call gp3d_mpi_finalize(mpi)
end program test_cufftmp_splitstep
