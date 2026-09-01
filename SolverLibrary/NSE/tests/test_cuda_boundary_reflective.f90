program test_cuda_boundary_reflective
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config, &
    nse_face_x_min, nse_face_x_max
  use mod_nse_boundary, only : apply_nse_boundary, validate_boundary_scheme
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_synchronize, nse_gpu_finalize
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse

  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 8
  sim%ny = 7
  sim%nz = 6
  sim%nghost = 3
  sim%x_min = -0.5_dp
  sim%x_max = 1.5_dp
  sim%y_min = 0.0_dp
  sim%y_max = 1.4_dp
  sim%z_min = -0.3_dp
  sim%z_max = 0.9_dp
  sim%dx = (sim%x_max-sim%x_min)/real(sim%nx,dp)
  sim%dy = (sim%y_max-sim%y_min)/real(sim%ny,dp)
  sim%dz = (sim%z_max-sim%z_min)/real(sim%nz,dp)
  sim%cuda_device = 0
  nse%viscous_scheme = 'central6'
  nse%reynolds = 100.0_dp
  nse%prandtl = 0.72_dp
  nse%boundary_condition = 'mixed'

  nse%boundary_face_type = 'reflective'
  call compare_cpu_and_cuda('all reflective', .true.)

  nse%boundary_face_type = 'periodic'
  nse%boundary_face_type(nse_face_x_min) = 'reflective'
  nse%boundary_face_type(nse_face_x_max) = 'reflective'
  call compare_cpu_and_cuda('x reflective, y/z periodic', .false.)

  write(*,'(A)') 'CUDA reflective boundary comparison passed'

contains

  subroutine compare_cpu_and_cuda(label, check_all_reflective)
    character(len=*), intent(in) :: label
    logical, intent(in) :: check_all_reflective
    type(nse_gpu_context) :: gpu
    real(dp), allocatable :: initial(:,:,:,:), q_cpu(:,:,:,:), q_cuda(:,:,:,:)
    real(dp) :: maximum_error
    integer :: i, j, k

    allocate(initial(1-sim%nghost:sim%nx+sim%nghost, &
      1-sim%nghost:sim%ny+sim%nghost, &
      1-sim%nghost:sim%nz+sim%nghost, nse%nv))
    allocate(q_cpu, mold=initial)
    allocate(q_cuda, mold=initial)
    initial = -huge(1.0_dp)
    do k = 1, sim%nz
      do j = 1, sim%ny
        do i = 1, sim%nx
          call analytic_state(i, j, k, nse%gamma, initial(i,j,k,:))
        end do
      end do
    end do

    call validate_boundary_scheme(sim, nse)
    q_cpu = initial
    call apply_nse_boundary(q_cpu, sim, nse, 1, sim%ny, 1, sim%nz)

    q_cuda = initial
    call nse_gpu_initialize(gpu, sim, nse)
    call nse_gpu_upload(gpu, q_cuda)
    call nse_gpu_synchronize(gpu)
    call nse_gpu_download(gpu, q_cuda)
    call nse_gpu_finalize(gpu)

    maximum_error = maxval(abs(q_cuda-q_cpu))
    if (maximum_error > 1.0e-13_dp) then
      write(*,'(A,A,A,ES24.16)') 'CUDA reflective mismatch (', &
        trim(label), '): ', maximum_error
      error stop 'CUDA and CPU reflective boundaries differ'
    end if
    if (check_all_reflective) call assert_all_reflective(q_cuda)
    deallocate(initial, q_cpu, q_cuda)
  end subroutine compare_cpu_and_cuda

  subroutine assert_all_reflective(q)
    real(dp), intent(in) :: q(1-sim%nghost:,1-sim%nghost:, &
      1-sim%nghost:,:)
    real(dp) :: expected(5), maximum_error
    integer :: i, j, k, source_i, source_j, source_k

    maximum_error = 0.0_dp
    do k = 1-sim%nghost, sim%nz+sim%nghost
      do j = 1-sim%nghost, sim%ny+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          source_i = mirror_index(i, sim%nx)
          source_j = mirror_index(j, sim%ny)
          source_k = mirror_index(k, sim%nz)
          call analytic_state(source_i, source_j, source_k, &
            nse%gamma, expected)
          if (i < 1 .or. i > sim%nx) expected(2) = -expected(2)
          if (j < 1 .or. j > sim%ny) expected(3) = -expected(3)
          if (k < 1 .or. k > sim%nz) expected(4) = -expected(4)
          maximum_error = max(maximum_error, maxval(abs(q(i,j,k,:)-expected)))
        end do
      end do
    end do
    if (maximum_error > 1.0e-13_dp) then
      write(*,'(A,ES24.16)') 'CUDA reflective definition error = ', &
        maximum_error
      error stop 'CUDA reflective face/edge/corner values are incorrect'
    end if
  end subroutine assert_all_reflective

  pure integer function mirror_index(index, extent) result(source)
    integer, intent(in) :: index, extent

    if (index < 1) then
      source = 1-index
    else if (index > extent) then
      source = 2*extent+1-index
    else
      source = index
    end if
  end function mirror_index

  pure subroutine analytic_state(i, j, k, gamma, state)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: state(5)
    real(dp) :: rho, pressure, velocity(3)

    rho = 0.96_dp + 0.004_dp*real(i,dp) + &
      0.002_dp*real(j,dp) + 0.001_dp*real(k,dp)
    velocity = [0.12_dp+0.002_dp*real(i,dp), &
      -0.03_dp+0.001_dp*real(j,dp), &
      0.018_dp-0.0005_dp*real(k,dp)]
    pressure = 1.0_dp/gamma + 0.001_dp*real(i+j+k,dp)
    state(1) = rho
    state(2:4) = rho*velocity
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*sum(velocity*velocity)
  end subroutine analytic_state

end program test_cuda_boundary_reflective
