program test_boundary_dirichlet
  use module_mpi
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config, &
    nse_face_x_min, nse_face_x_max
  use mod_nse_boundary, only : apply_nse_boundary, validate_boundary_scheme
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:)
  real(dp) :: low_state(5), high_state(5), expected(5), maximum_error
  integer :: i, j, k, source_j, source_k
  integer :: js, je, ks, ke, ierr

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 12
  sim%ny = 12
  sim%nz = 12
  sim%nghost = 3
  sim%x_min = 0.0_dp
  sim%x_max = 1.0_dp
  sim%y_min = 0.0_dp
  sim%y_max = 1.0_dp
  sim%z_min = 0.0_dp
  sim%z_max = 1.0_dp
  nse%boundary_face_type = 'periodic'
  nse%viscous_scheme = 'central6'
  nse%boundary_face_type(nse_face_x_min) = 'dirichlet'
  nse%boundary_face_type(nse_face_x_max) = 'dirichlet'
  nse%boundary_condition = 'mixed'
  nse%boundary_reference_rho(nse_face_x_min) = 2.0_dp
  nse%boundary_reference_velocity(:,nse_face_x_min) = &
    [0.7_dp, 0.1_dp, -0.05_dp]
  nse%boundary_reference_p(nse_face_x_min) = 1.8_dp
  nse%boundary_reference_rho(nse_face_x_max) = 0.9_dp
  nse%boundary_reference_velocity(:,nse_face_x_max) = &
    [-0.2_dp, 0.03_dp, 0.04_dp]
  nse%boundary_reference_p(nse_face_x_max) = 0.65_dp
  call primitive_to_conserved(2.0_dp, [0.7_dp,0.1_dp,-0.05_dp], &
    1.8_dp, nse%gamma, low_state)
  call primitive_to_conserved(0.9_dp, [-0.2_dp,0.03_dp,0.04_dp], &
    0.65_dp, nse%gamma, high_state)

  call validate_boundary_scheme(sim, nse)
  call mp_setup_division(sim%nx, sim%ny, sim%nz)
  js = j_sta
  je = j_end
  ks = k_sta
  ke = k_end
  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    js-sim%nghost:je+sim%nghost, ks-sim%nghost:ke+sim%nghost, nse%nv))
  q = -huge(1.0_dp)
  do k = ks, ke
    do j = js, je
      do i = 1, sim%nx
        call analytic_state(i, j, k, nse%gamma, q(i,j,k,:))
      end do
    end do
  end do

  !$OMP PARALLEL DEFAULT(SHARED)
  call apply_nse_boundary(q, sim, nse, js, je, ks, ke)
  !$OMP END PARALLEL

  maximum_error = 0.0_dp
  do k = ks-sim%nghost, ke+sim%nghost
    do j = js-sim%nghost, je+sim%nghost
      source_j = periodic_index(j, sim%ny)
      source_k = periodic_index(k, sim%nz)
      do i = 1-sim%nghost, sim%nx+sim%nghost
        if (i < 1) then
          expected = low_state
        else if (i > sim%nx) then
          expected = high_state
        else
          call analytic_state(i, source_j, source_k, nse%gamma, expected)
        end if
        maximum_error = max(maximum_error, maxval(abs(q(i,j,k,:)-expected)))
      end do
    end do
  end do
  if (maximum_error > 1.0e-12_dp) then
    write(*,'(A,I0,A,ES24.16)') 'rank ', my_rank, &
      ' Dirichlet boundary error = ', maximum_error
    error stop 'Dirichlet boundary face/edge/corner test failed'
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (my_rank == root) write(*,'(A)') &
    'Dirichlet boundary face/edge/corner test passed'
  deallocate(q)
  call MPI_Finalize(ierr)

contains

  pure integer function periodic_index(index, extent) result(source)
    integer, intent(in) :: index, extent
    source = modulo(index-1, extent) + 1
  end function periodic_index

  pure subroutine primitive_to_conserved(rho, velocity, pressure, gamma, state)
    real(dp), intent(in) :: rho, velocity(3), pressure, gamma
    real(dp), intent(out) :: state(5)
    state(1) = rho
    state(2:4) = rho*velocity
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*sum(velocity*velocity)
  end subroutine primitive_to_conserved

  pure subroutine analytic_state(i, j, k, gamma, state)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: state(5)
    real(dp) :: rho, pressure, velocity(3)
    rho = 0.95_dp + 0.003_dp*real(i,dp) + &
      0.002_dp*real(j,dp) + 0.001_dp*real(k,dp)
    velocity = [0.08_dp+0.001_dp*real(i,dp), &
      -0.04_dp+0.0015_dp*real(j,dp), &
      0.02_dp-0.0007_dp*real(k,dp)]
    pressure = 1.0_dp/gamma + 0.0005_dp*real(i+j+k,dp)
    call primitive_to_conserved(rho, velocity, pressure, gamma, state)
  end subroutine analytic_state

end program test_boundary_dirichlet
