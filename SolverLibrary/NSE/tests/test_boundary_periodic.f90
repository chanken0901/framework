program test_boundary_periodic
  use module_mpi
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_nse_boundary, only : apply_nse_boundary
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:)
  real(dp) :: expected, maximum_error
  integer :: i, j, k, l, wi, wj, wk
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
  nse%viscous_scheme = 'central6'

  call mp_setup_division(sim%nx, sim%ny, sim%nz)
  js = j_sta
  je = j_end
  ks = k_sta
  ke = k_end

  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    js-sim%nghost:je+sim%nghost, &
    ks-sim%nghost:ke+sim%nghost, nse%nv))
  q = -huge(1.0_dp)

  do l = 1, nse%nv
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          q(i,j,k,l) = test_value(i, j, k, l)
        end do
      end do
    end do
  end do

  !$OMP PARALLEL DEFAULT(SHARED)
  call apply_nse_boundary(q, sim, nse, js, je, ks, ke)
  !$OMP END PARALLEL

  maximum_error = 0.0_dp
  do l = 1, nse%nv
    do k = ks-sim%nghost, ke+sim%nghost
      wk = 1 + modulo(k-1, sim%nz)
      do j = js-sim%nghost, je+sim%nghost
        wj = 1 + modulo(j-1, sim%ny)
        do i = 1-sim%nghost, sim%nx+sim%nghost
          wi = 1 + modulo(i-1, sim%nx)
          expected = test_value(wi, wj, wk, l)
          maximum_error = max(maximum_error, abs(q(i,j,k,l)-expected))
        end do
      end do
    end do
  end do

  if (maximum_error > 0.0_dp) then
    write(*,'(A,I0,A,ES16.8)') 'rank ', my_rank, &
      ' periodic boundary maximum error = ', maximum_error
    error stop 'periodic boundary face/edge/corner test failed'
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (my_rank == root) write(*,'(A)') &
    'periodic boundary face/edge/corner test passed'
  deallocate(q)
  call MPI_Finalize(ierr)

contains

  pure real(dp) function test_value(ii, jj, kk, ll) result(value)
    integer, intent(in) :: ii, jj, kk, ll

    value = real(ii + 100*jj + 10000*kk + 1000000*ll, dp)
  end function test_value

end program test_boundary_periodic
