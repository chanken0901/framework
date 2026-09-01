program test_boundary_reflective
  use module_mpi
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_nse_boundary, only : apply_nse_boundary, validate_boundary_scheme
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:)
  real(dp) :: expected(5), maximum_error
  integer :: i, j, k, source_i, source_j, source_k
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
  nse%viscous_scheme = 'central6'
  nse%reynolds = 100.0_dp
  nse%prandtl = 0.72_dp
  nse%boundary_face_type = 'reflective'
  nse%boundary_condition = 'mixed'

  call validate_boundary_scheme(sim, nse)
  call mp_setup_division(sim%nx, sim%ny, sim%nz)
  js = j_sta
  je = j_end
  ks = k_sta
  ke = k_end
  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    js-sim%nghost:je+sim%nghost, &
    ks-sim%nghost:ke+sim%nghost, nse%nv))
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
      do i = 1-sim%nghost, sim%nx+sim%nghost
        source_i = mirror_index(i, sim%nx)
        source_j = mirror_index(j, sim%ny)
        source_k = mirror_index(k, sim%nz)
        call analytic_state(source_i, source_j, source_k, nse%gamma, expected)
        if (i < 1 .or. i > sim%nx) expected(2) = -expected(2)
        if (j < 1 .or. j > sim%ny) expected(3) = -expected(3)
        if (k < 1 .or. k > sim%nz) expected(4) = -expected(4)
        maximum_error = max(maximum_error, maxval(abs(q(i,j,k,:)-expected)))
      end do
    end do
  end do
  if (maximum_error > 1.0e-12_dp) then
    write(*,'(A,I0,A,ES24.16)') 'rank ', my_rank, &
      ' reflective boundary error = ', maximum_error
    error stop 'reflective boundary face/edge/corner test failed'
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (my_rank == root) write(*,'(A)') &
    'reflective boundary face/edge/corner test passed'
  deallocate(q)
  call MPI_Finalize(ierr)

contains

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

    rho = 0.95_dp + 0.003_dp*real(i,dp) + &
      0.002_dp*real(j,dp) + 0.001_dp*real(k,dp)
    velocity = [0.08_dp+0.001_dp*real(i,dp), &
      -0.04_dp+0.0015_dp*real(j,dp), &
      0.02_dp-0.0007_dp*real(k,dp)]
    pressure = 1.0_dp/gamma + 0.0005_dp*real(i+j+k,dp)
    state(1) = rho
    state(2:4) = rho*velocity
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*sum(velocity*velocity)
  end subroutine analytic_state

end program test_boundary_reflective
