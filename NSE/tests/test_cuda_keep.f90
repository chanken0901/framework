program test_cuda_keep
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_compute_dt, &
    nse_gpu_advance_ssprk3, nse_gpu_synchronize, nse_gpu_finalize
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: gpu
  real(dp), allocatable :: q(:,:,:,:)
  real(dp) :: state(5), dt, expected_dt, error
  integer :: i, j, k, variable

  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 6
  sim%ny = 5
  sim%nz = 4
  sim%nghost = 3
  sim%dx = 0.2_dp
  sim%dy = 0.25_dp
  sim%dz = 0.3_dp
  sim%cuda_device = 0
  nse%cfl = 0.4_dp

  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost, nse%nv))
  q = 0.0_dp
  state = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    1.0_dp/(nse%gamma-1.0_dp)]
  do k = 1, sim%nz
    do j = 1, sim%ny
      do i = 1, sim%nx
        q(i,j,k,:) = state
      end do
    end do
  end do

  call nse_gpu_initialize(gpu, sim, nse)
  call nse_gpu_upload(gpu, q)
  call nse_gpu_compute_dt(gpu, dt)
  expected_dt = nse%cfl * min(sim%dx, sim%dy, sim%dz) / sqrt(nse%gamma)
  if (abs(dt-expected_dt) > 1.0e-13_dp) then
    write(*,'(A,2ES24.16)') "CUDA CFL mismatch: ", dt, expected_dt
    error stop "CUDA CFL test failed"
  end if

  call nse_gpu_advance_ssprk3(gpu, 1.0e-3_dp)
  call nse_gpu_synchronize(gpu)
  call nse_gpu_download(gpu, q)
  error = 0.0_dp
  do variable = 1, nse%nv
    error = max(error, maxval(abs( &
      q(1:sim%nx,1:sim%ny,1:sim%nz,variable)-state(variable))))
  end do
  if (error > 1.0e-13_dp) then
    write(*,'(A,ES24.16)') "Uniform-state error: ", error
    error stop "CUDA KEEP uniform-state test failed"
  end if
  if (maxval(abs(q(0,1:sim%ny,1:sim%nz,:)- &
      q(sim%nx,1:sim%ny,1:sim%nz,:))) > 1.0e-13_dp) then
    error stop "CUDA periodic halo test failed"
  end if

  call nse_gpu_finalize(gpu)
  deallocate(q)
  write(*,'(A,ES16.8)') "CUDA KEEP test passed; max error = ", error
end program test_cuda_keep
