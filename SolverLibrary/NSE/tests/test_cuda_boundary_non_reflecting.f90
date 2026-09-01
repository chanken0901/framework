program test_cuda_boundary_non_reflecting
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config, &
    nse_boundary_face_count, nse_face_x_min, nse_face_x_max
  use mod_nse_boundary, only : apply_nse_boundary, validate_boundary_scheme
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_synchronize, nse_gpu_finalize
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  integer :: face

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
  nse%boundary_reference_rho = 1.0_dp
  nse%boundary_reference_p = 1.0_dp/nse%gamma
  do face = 1, nse_boundary_face_count
    nse%boundary_reference_velocity(:,face) = &
      [0.18_dp, -0.03_dp, 0.02_dp]
  end do
  nse%boundary_relaxation_strength = 0.15_dp
  nse%boundary_length_scale = -1.0_dp

  nse%boundary_face_type = 'non_reflecting'
  call compare_cpu_and_cuda('all non-reflecting', .false.)

  nse%boundary_face_type = 'periodic'
  nse%boundary_face_type(nse_face_x_min) = 'non_reflecting'
  nse%boundary_face_type(nse_face_x_max) = 'non_reflecting'
  call compare_cpu_and_cuda('x non-reflecting, y/z periodic', .true.)

  write(*,'(A)') 'CUDA non-reflecting boundary comparison passed'

contains

  subroutine compare_cpu_and_cuda(label, require_nonperiodic_x)
    character(len=*), intent(in) :: label
    logical, intent(in) :: require_nonperiodic_x
    type(nse_gpu_context) :: gpu
    real(dp), allocatable :: initial(:,:,:,:), q_cpu(:,:,:,:), q_cuda(:,:,:,:)
    real(dp) :: maximum_error, rho, pressure, velocity(3)
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
          rho = 0.97_dp + 0.004_dp*real(i,dp) + &
            0.002_dp*real(j,dp) + 0.001_dp*real(k,dp)
          velocity = [0.16_dp + 0.002_dp*real(i,dp), &
            -0.025_dp + 0.001_dp*real(j,dp), &
            0.018_dp - 0.0005_dp*real(k,dp)]
          pressure = 1.0_dp/nse%gamma + 0.001_dp*real(i+j+k,dp)
          call make_conserved(rho, velocity, pressure, nse%gamma, &
            initial(i,j,k,:))
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
    if (.not. ieee_is_finite(maximum_error) .or. &
        maximum_error > 5.0e-12_dp) then
      write(*,'(A,A,A,ES24.16)') 'CUDA boundary mismatch (', &
        trim(label), '): ', maximum_error
      error stop 'CUDA and CPU boundary results differ'
    end if
    if (require_nonperiodic_x) then
      if (abs(q_cuda(0,1,1,1)-q_cuda(sim%nx,1,1,1)) < 1.0e-6_dp) then
        error stop 'CUDA x non-reflecting face was replaced by periodic wrapping'
      end if
    end if
    call assert_physical(q_cuda, label)
    deallocate(initial, q_cpu, q_cuda)
  end subroutine compare_cpu_and_cuda

  subroutine assert_physical(q, label)
    real(dp), intent(in) :: q(1-sim%nghost:,1-sim%nghost:, &
      1-sim%nghost:,:)
    character(len=*), intent(in) :: label
    real(dp) :: rho, pressure
    integer :: i, j, k

    do k = 1-sim%nghost, sim%nz+sim%nghost
      do j = 1-sim%nghost, sim%ny+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          rho = q(i,j,k,1)
          pressure = conserved_pressure(q(i,j,k,:), nse%gamma)
          if (.not. ieee_is_finite(rho) .or. rho <= 0.0_dp .or. &
              .not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
            write(*,'(A,A,A,3I6,2ES16.8)') 'invalid CUDA boundary state (', &
              trim(label), '): ', i, j, k, rho, pressure
            error stop 'CUDA boundary generated a nonphysical state'
          end if
        end do
      end do
    end do
  end subroutine assert_physical

  subroutine make_conserved(rho, velocity, pressure, gamma, state)
    real(dp), intent(in) :: rho, velocity(3), pressure, gamma
    real(dp), intent(out) :: state(5)

    state(1) = rho
    state(2:4) = rho*velocity
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*sum(velocity*velocity)
  end subroutine make_conserved

  pure real(dp) function conserved_pressure(state, gamma) result(value)
    real(dp), intent(in) :: state(5), gamma
    real(dp) :: velocity(3)

    velocity = state(2:4)/state(1)
    value = (gamma-1.0_dp) * &
      (state(5)-0.5_dp*state(1)*sum(velocity*velocity))
  end function conserved_pressure

end program test_cuda_boundary_non_reflecting
