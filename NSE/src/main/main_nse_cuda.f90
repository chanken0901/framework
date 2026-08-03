program main_nse_cuda
  use, intrinsic :: iso_fortran_env, only : int64
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config, &
    print_simulation_config, should_output
  use mod_model_config, only : nse_config, init_nse_config, print_nse_config
  use mod_input_reader, only : read_all_inputs
  use mod_grid_fvm, only : build_uniform_grid
  use mod_nse_initial_conditions, only : initialize_nse_state
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_initialize, &
    nse_gpu_upload, nse_gpu_download, nse_gpu_compute_dt, &
    nse_gpu_advance_ssprk3, nse_gpu_synchronize, nse_gpu_finalize
  use mod_slf_output, only : write_nse_conserved_slf, write_meta_json
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: gpu
  real(dp), allocatable :: q(:,:,:,:)
  integer :: js, je, ks, ke
  integer(int64) :: clock_start, clock_end, clock_rate
  character(len=512) :: input_path

  input_path = "input.dat"
  if (command_argument_count() >= 1) call get_command_argument(1, input_path)

  call init_simulation_config(sim)
  call init_nse_config(nse)
  call read_all_inputs(trim(input_path), sim, nse=nse)
  sim%backend = "cuda"
  sim%use_mpi = .false.
  sim%use_openmp = .false.
  sim%rank = 0
  sim%nprocs = 1

  call print_simulation_config(sim)
  call print_nse_config(nse)

  js = 1
  je = sim%ny
  ks = 1
  ke = sim%nz
  call build_uniform_grid(sim, js, je, ks, ke)

  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost, nse%nv))
  q = 0.0_dp
  call initialize_nse_state(q, sim, nse, js, je, ks, ke)

  call nse_gpu_initialize(gpu, sim, nse)
  call nse_gpu_upload(gpu, q)

  if (sim%write_meta) then
    call write_meta_json(sim, is=1, ie=sim%nx, js=js, je=je, ks=ks, &
      ke=ke, use_cuda=.true.)
  end if
  if (should_output(sim, 0)) then
    call nse_gpu_download(gpu, q)
    call write_nse_conserved_slf(sim, 0, 0.0_dp, q, rank=0)
  end if

  sim%t = 0.0_dp
  sim%step = 0
  sim%ttotal = 0.0_dp
  call system_clock(count_rate=clock_rate)
  write(*,'(A)') "# step time dt step_wall_seconds total_wall_seconds"

  do while (sim%t < sim%t_max .and. sim%step < sim%nsteps)
    call system_clock(clock_start)
    if (.not. sim%use_fixed_dt) call nse_gpu_compute_dt(gpu, sim%dt)
    if (sim%t + sim%dt > sim%t_max) sim%dt = sim%t_max - sim%t
    if (sim%dt <= 0.0_dp) exit

    call nse_gpu_advance_ssprk3(gpu, sim%dt)
    call nse_gpu_synchronize(gpu)
    sim%t = sim%t + sim%dt
    sim%step = sim%step + 1

    if (should_output(sim, sim%step)) then
      call nse_gpu_download(gpu, q)
      call write_nse_conserved_slf(sim, sim%step, sim%t, q, rank=0)
    end if

    call system_clock(clock_end)
    sim%t2 = real(clock_end-clock_start, dp) / real(clock_rate, dp)
    sim%ttotal = sim%ttotal + sim%t2
    write(*,'(I10,1X,4(ES16.8,1X))') sim%step, sim%t, sim%dt, &
      sim%t2, sim%ttotal
  end do

  call nse_gpu_finalize(gpu)
  deallocate(q)
  write(*,'(A,I0,A,ES16.8)') &
    'NSE calculation completed successfully: step=', sim%step, &
    ', time=', sim%t
end program main_nse_cuda
