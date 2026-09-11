program main_nse_mpi_cuda
  use module_mpi
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config, &
    print_simulation_config, should_output
  use mod_model_config, only : nse_config, init_nse_config, &
    print_nse_config, resolve_nse_flow_parameters
  use mod_nse_forcing_common, only : forcing_is_enabled
  use mod_input_reader, only : read_all_inputs
  use mod_grid_fvm, only : build_uniform_grid
  use mod_nse_initial_conditions, only : initialize_nse_state
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_select_device, nse_gpu_configure_cufftmp, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_compute_dt, &
    nse_gpu_begin_ssprk3, nse_gpu_advance_ssprk3_stage, nse_gpu_restore_ssprk3, &
    nse_gpu_synchronize, nse_gpu_finalize
  use mod_nse_gpu_mpi, only : nse_gpu_mpi_halo, &
    nse_gpu_mpi_halo_initialize, nse_gpu_mpi_exchange, &
    nse_gpu_mpi_halo_finalize
  use mod_slf_output, only : write_nse_conserved_slf, write_meta_json
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: gpu
  type(nse_gpu_mpi_halo) :: halo
  real(dp), allocatable :: q(:,:,:,:)
  integer :: js, je, ks, ke, stage, selected_device
  integer :: retry, step_status, global_status
  real(dp) :: requested_dt, fh_dt_limit
  integer :: ierror_local, ierr_local, local_cells, minimum_local_cells
  character(len=512) :: input_path

  call MPI_Init(ierror_local)
  if (ierror_local /= 0) error stop 'MPI_Init failed'
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr_local)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr_local)

  input_path = 'input.dat'
  if (command_argument_count() >= 1) call get_command_argument(1, input_path)

  call init_simulation_config(sim)
  call init_nse_config(nse)
  call read_all_inputs(trim(input_path), sim, nse=nse)
  call resolve_nse_flow_parameters(nse, sim%initial_condition)
  sim%backend = 'cuda_mpi'
  sim%use_mpi = .true.
  sim%use_openmp = .false.
  sim%rank = my_rank
  sim%nprocs = nprocs

  call mp_setup_division(sim%nx, sim%ny, sim%nz)
  js = j_sta
  je = j_end
  ks = k_sta
  ke = k_end
  local_cells = min(je-js+1, ke-ks+1)
  call MPI_Allreduce(local_cells, minimum_local_cells, 1, MPI_INTEGER, &
    MPI_MIN, MPI_COMM_WORLD, ierr_local)
  if (ierr_local /= 0) error stop 'MPI local-domain validation failed'
  if (minimum_local_cells < sim%nghost) then
    if (my_rank == root) then
      write(*,'(A,I0,A,I0)') &
        'ERROR: every MPI Y/Z block needs at least nghost cells; minimum=', &
        minimum_local_cells, ', nghost=', sim%nghost
    end if
    call MPI_Abort(MPI_COMM_WORLD, 11, ierr_local)
  end if
  selected_device = resolve_local_cuda_device(sim%cuda_device)
  call nse_gpu_select_device(selected_device)

  if (my_rank == root) then
    call print_simulation_config(sim)
    call print_nse_config(nse)
    write(*,'(A,I0,A,I0,A)') '# MPI+CUDA process grid=', ndiv_ny, 'x', &
      ndiv_nz, ' (YxZ)'
    write(*,'(A)') '# MPI+CUDA halo transport=host_staged'
  end if
  write(*,'(A,I0,A,I0,A,4(I0,1X))') '# rank=', my_rank, &
    ' cuda_device=', selected_device, ' local_range=', js, je, ks, ke

  call build_uniform_grid(sim, js, je, ks, ke)
  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    js-sim%nghost:je+sim%nghost, ks-sim%nghost:ke+sim%nghost, nse%nv))
  q = 0.0_dp
  call initialize_nse_state(q, sim, nse, js, je, ks, ke)

  call nse_gpu_initialize(gpu, sim, nse, local_ny=je-js+1, &
    local_nz=ke-ks+1, distributed_y=ndiv_ny>1, &
    distributed_z=ndiv_nz>1, device=selected_device, &
    global_y_start=js-1, global_z_start=ks-1)
  if (forcing_is_enabled(nse)) then
    call nse_gpu_configure_cufftmp(gpu, sim%ny, sim%nz, js-1, ks-1, &
      MPI_COMM_WORLD)
  end if
  call nse_gpu_upload(gpu, q)
  call nse_gpu_mpi_halo_initialize(halo, gpu, nse)
  call nse_gpu_mpi_exchange(halo, gpu)

  if (sim%write_meta) then
    call write_meta_json(sim, is=1, ie=sim%nx, js=js, je=je, ks=ks, &
      ke=ke, use_cuda=.true.)
  end if
  if (should_output(sim, 0)) then
    call nse_gpu_download(gpu, q)
    call write_nse_conserved_slf(sim, 0, 0.0_dp, q, rank=my_rank)
  end if

  sim%t = 0.0_dp
  sim%step = 0
  sim%ttotal = 0.0_dp
  if (my_rank == root) then
    write(*,'(A)') '# step time dt step_wall_seconds total_wall_seconds'
  end if

  do while (sim%t < sim%t_max .and. sim%step < sim%nsteps)
    sim%t1 = MPI_Wtime()
    if (.not. sim%use_fixed_dt) then
      call nse_gpu_compute_dt(gpu, sim%dt)
      call mp_allminr8(sim%dt)
    end if
    if (sim%t + sim%dt > sim%t_max) sim%dt = sim%t_max - sim%t
    if (sim%dt <= 0.0_dp) exit

    if(nse%fh_enabled) then
      call nse_gpu_compute_dt(gpu,fh_dt_limit)
      call mp_allminr8(fh_dt_limit)
      if(sim%dt>fh_dt_limit) then
        if(my_rank==root) write(*,'(A)') 'ERROR: LLNS dt exceeds global stability bound'
        call MPI_Abort(MPI_COMM_WORLD,15,ierr_local)
      end if
    end if

    call nse_gpu_begin_ssprk3(gpu, sim%dt)
    requested_dt = sim%dt
    do retry = 0, 20
      do stage = 1, 3
        call nse_gpu_mpi_exchange(halo, gpu)
        call nse_gpu_advance_ssprk3_stage(gpu, sim%dt, stage, step_status)
        ! Fatal errors outrank a recoverable positivity rejection.
        if (step_status == 1) step_status = 3
        call MPI_Allreduce(step_status, global_status, 1, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, ierr_local)
        if (ierr_local /= 0) call MPI_Abort(MPI_COMM_WORLD, 12, ierr_local)
        if (global_status /= 0) exit
      end do
      if (global_status == 0) exit
      call nse_gpu_restore_ssprk3(gpu)
      ! LLNS dt is chosen before noise; outcome-dependent retries are forbidden.
      if (global_status /= 2 .or. sim%use_fixed_dt .or. nse%fh_enabled .or. retry == 20) then
        if (my_rank == root) write(*,'(A,I0)') &
          'ERROR: CUDA step rejected; state restored, status=', global_status
        call MPI_Abort(MPI_COMM_WORLD, 13, ierr_local)
        error stop 'CUDA positivity retry failed'
      end if
      sim%dt = 0.5_dp*sim%dt
    end do
    if (my_rank == root .and. retry > 0) write(*,'(A,I0,A,2ES16.8)') &
      '# positivity_retry count=', retry, ' requested/accepted dt=', requested_dt, sim%dt
    if (sim%t + sim%dt <= sim%t) call MPI_Abort(MPI_COMM_WORLD, 14, ierr_local)
    call nse_gpu_synchronize(gpu)

    sim%t = sim%t + sim%dt
    sim%step = sim%step + 1
    if (should_output(sim, sim%step)) then
      call nse_gpu_mpi_exchange(halo, gpu)
      call nse_gpu_download(gpu, q)
      call write_nse_conserved_slf(sim, sim%step, sim%t, q, rank=my_rank)
    end if

    sim%t2 = MPI_Wtime()
    sim%ttotal = sim%ttotal + (sim%t2 - sim%t1)
    if (my_rank == root) then
      write(*,'(I10,1X,4(ES16.8,1X))') sim%step, sim%t, sim%dt, &
        sim%t2-sim%t1, sim%ttotal
    end if
  end do

  call nse_gpu_mpi_halo_finalize(halo)
  call nse_gpu_finalize(gpu)
  deallocate(q)
  call MPI_Finalize(ierr_local)
  if (ierr_local /= 0) error stop 'MPI_Finalize failed'
  if (my_rank == root) then
    write(*,'(A,I0,A,ES16.8)') &
      'NSE MPI+CUDA calculation completed successfully: step=', sim%step, &
      ', time=', sim%t
  end if

contains

  integer function resolve_local_cuda_device(base_device) result(device)
    integer, intent(in) :: base_device
    character(len=32) :: policy
    character(len=1024) :: visible_devices
    integer :: local_rank, local_comm, status, ierr

    policy = 'local_rank'
    call get_environment_variable('NSE_CUDA_DEVICE_POLICY', policy, &
      status=status)
    if (status /= 0 .or. len_trim(policy) == 0) policy = 'local_rank'
    if (trim(adjustl(policy)) == 'fixed') then
      device = max(0, base_device)
      return
    end if
    if (trim(adjustl(policy)) /= 'local_rank') then
      error stop 'NSE_CUDA_DEVICE_POLICY must be local_rank or fixed'
    end if

    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, my_rank, &
      MPI_INFO_NULL, local_comm, ierr)
    if (ierr /= 0) error stop 'MPI_Comm_split_type failed for CUDA binding'
    call MPI_Comm_rank(local_comm, local_rank, ierr)
    if (ierr /= 0) error stop 'MPI local-rank query failed for CUDA binding'
    call MPI_Comm_free(local_comm, ierr)
    if (ierr /= 0) error stop 'MPI local communicator release failed'

    visible_devices = ''
    call get_environment_variable('CUDA_VISIBLE_DEVICES', visible_devices, &
      status=status)
    if (status == 0 .and. len_trim(visible_devices) > 0 .and. &
        index(trim(visible_devices), ',') == 0) then
      ! A scheduler may expose exactly one GPU to each rank.  CUDA then
      ! renumbers that device to zero within the rank.
      device = max(0, base_device)
    else
      device = max(0, base_device) + local_rank
    end if
  end function resolve_local_cuda_device

end program main_nse_mpi_cuda
