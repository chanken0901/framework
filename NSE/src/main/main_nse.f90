program main
  use module_mpi
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config, &
    print_simulation_config, should_output
  use mod_model_config, only : nse_config, init_nse_config, print_nse_config
  use mod_input_reader, only : read_all_inputs
  use mod_grid_fvm, only : build_uniform_grid
  use mod_nse_field, only : allocate_nse_fields, deallocate_nse_fields, &
    Q, Q0, RHS, F
  use mod_nse_initial_conditions, only : initialize_nse_state
  use mod_nse_boundary, only : apply_nse_boundary
  use mod_nse_spatial_operator, only : validate_nse_spatial_configuration
  use mod_nse_forcing, only : initialize_nse_forcing, finalize_nse_forcing
  use mod_nse_time_integration, only : compute_nse_dt, advance_nse_ssprk3, &
    validate_time_integrator
  use mod_slf_output, only : write_nse_conserved_slf, write_meta_json
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  integer :: js, je, ks, ke
  integer :: ierror_local, ierr_local
  character(len=512) :: input_path

  call MPI_Init(ierror_local)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr_local)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr_local)

  input_path = 'input.dat'
  if (command_argument_count() >= 1) call get_command_argument(1, input_path)

  call init_simulation_config(sim)
  call init_nse_config(nse)
  call read_all_inputs(trim(input_path), sim, nse=nse)
  sim%rank = my_rank
  sim%nprocs = nprocs
  sim%use_mpi = .true.

  if (my_rank == root) then
    call print_simulation_config(sim)
    call print_nse_config(nse)
  end if

  call validate_nse_spatial_configuration(sim, nse)
  call validate_time_integrator(nse)
  call mp_setup_division(sim%nx, sim%ny, sim%nz)

  js = j_sta
  je = j_end
  ks = k_sta
  ke = k_end

  call build_uniform_grid(sim, js, je, ks, ke)
  if (my_rank == root) write(*,'(A)') 'Complete build grid'

  call allocate_nse_fields(sim, nse, js, je, ks, ke)
  if (my_rank == root) write(*,'(A)') 'Complete allocate'

  call initialize_nse_state(Q, sim, nse, js, je, ks, ke)
  call apply_nse_boundary(Q, sim, nse, js, je, ks, ke)
  call initialize_nse_forcing(sim, nse, js, je, ks, ke)
  if (my_rank == root) write(*,'(A)') 'Complete initialize'

  call write_meta_json(sim, is=1, ie=sim%nx, js=js, je=je, ks=ks, ke=ke, &
    use_cuda=.false.)
  if (should_output(sim, 0)) then
    call write_nse_conserved_slf(sim, 0, 0.0_dp, Q, rank=my_rank)
  end if

  sim%t = 0.0_dp
  sim%step = 0
  sim%ttotal = 0.0_dp

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(Q,Q0,RHS,F,js,je,ks,ke,my_rank,sim,nse)
  do while (sim%t < sim%t_max .and. sim%step < sim%nsteps)
    !$OMP MASKED
    sim%t1 = MPI_Wtime()
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP MASKED
    if (.not. sim%use_fixed_dt) then
      call compute_nse_dt(Q, sim%dt, sim, nse, js, je, ks, ke)
    end if
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP MASKED
    if (sim%t + sim%dt > sim%t_max) sim%dt = sim%t_max - sim%t
    !$OMP END MASKED
    !$OMP BARRIER

    call advance_nse_ssprk3(Q, Q0, RHS, F, sim%dt, sim, nse, js, je, ks, ke)

    !$OMP MASKED
    sim%t = sim%t + sim%dt
    sim%step = sim%step + 1
    if (should_output(sim, sim%step)) then
      call write_nse_conserved_slf(sim, sim%step, sim%t, Q, rank=my_rank)
    end if
    sim%t2 = MPI_Wtime()
    sim%ttotal = sim%ttotal + (sim%t2 - sim%t1)
    if (my_rank == root) then
      write(*,*) sim%step, sim%t, sim%dt, sim%t2-sim%t1, sim%ttotal
    end if
    call mp_barrier
    !$OMP END MASKED
    !$OMP BARRIER
  end do
  !$OMP END PARALLEL

  call finalize_nse_forcing()
  call deallocate_nse_fields()
  call MPI_Finalize(ierr_local)
  if (ierr_local /= 0) error stop 'MPI_Finalize failed'
  if (my_rank == root) then
    write(*,'(A,I0,A,ES16.8)') &
      'NSE calculation completed successfully: step=', sim%step, &
      ', time=', sim%t
  end if
end program main
