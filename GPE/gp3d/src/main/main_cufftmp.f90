!> MPI + cuFFTMp を用いる複数 GPU 向け GP3D 実行プログラム。
!> 各 MPI rank は一つの GPU と実空間 z スラブを所有する。
program gp3d_cufftmp
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_params_t, gp3d_state_t, &
    gp3d_run_config_t, gp3d_model_config_t
  use gp3d_solver, only: gp3d_step_timing_t, gp3d_report_step_timing
  use gp3d_gpu, only: gp3d_gpu_context_t, gp3d_gpu_init, gp3d_gpu_upload, &
    gp3d_gpu_download, gp3d_gpu_step, gp3d_gpu_relax_taylor_green_argle, &
    gp3d_gpu_diagnostics, gp3d_gpu_finalize
  use gp3d_io, only: gp3d_output_config_t, gp3d_write_meta_json, &
    gp3d_write_gpe_psi_slf
  use gp3d_input, only: gp3d_read_all_inputs, gp3d_configure_problem
  use gp3d_restart, only: gp3d_restart_info_t
  use gp3d_mpi, only: gp3d_mpi_t, gp3d_mpi_init, gp3d_mpi_finalize, &
    gp3d_mpi_is_root, gp3d_mpi_barrier
  implicit none

  type(gp3d_mpi_t) :: mpi
  type(gp3d_run_config_t) :: run_cfg
  type(gp3d_model_config_t) :: model_cfg
  type(gp3d_grid_t) :: grid
  type(gp3d_params_t) :: params
  type(gp3d_state_t) :: state
  type(gp3d_gpu_context_t) :: gpu
  type(gp3d_output_config_t) :: output_cfg, seed_output_cfg
  type(gp3d_step_timing_t) :: step_timing
  type(gp3d_restart_info_t) :: restart_info
  integer :: step, advance_step, start_step
  real(dp) :: norm_value, energy_value, start_time, current_time
  logical :: output_step
  character(len=256) :: input_file

  call gp3d_mpi_init(mpi)

  input_file = "input.nml"
  if (command_argument_count() >= 1) call get_command_argument(1, input_file)

  ! 実行時の rank 数を namelist の指定より優先する。
  call gp3d_read_all_inputs(trim(input_file), run_cfg, model_cfg)
  run_cfg%use_mpi = .true.
  run_cfg%use_cuda = .true.
  run_cfg%cuda_device = -1
  run_cfg%rank = mpi%rank
  run_cfg%nprocs = mpi%nprocs

  call gp3d_configure_problem(run_cfg, model_cfg, grid, params, state, &
    output_cfg, mpi, restart_info)
  start_step = restart_info%step
  start_time = restart_info%time
  if (restart_info%loaded .and. gp3d_mpi_is_root(mpi)) then
    write(*,'(a,a)') "# restart source=", trim(restart_info%source_file)
    write(*,'(a,i0,1x,a,es16.8,1x,a,i0)') "# restart step=", start_step, &
      "time=", start_time, "source_nprocs=", restart_info%source_nprocs
  end if

  ! ARGLE 前の解析的 Taylor-Green seed はホスト側にあるため、そのまま rank 別に保存できる。
  if (.not. restart_info%loaded .and. model_cfg%argle_enabled .and. &
      model_cfg%argle_write_seed) then
    seed_output_cfg = output_cfg
    seed_output_cfg%output_dir = trim(output_cfg%output_dir) // "/argle_seed"
    seed_output_cfg%case_name = trim(output_cfg%case_name) // "_argle_seed"
    if (gp3d_mpi_is_root(mpi)) call gp3d_write_meta_json(seed_output_cfg)
    call gp3d_mpi_barrier(mpi)
    call gp3d_write_gpe_psi_slf(seed_output_cfg, step=0, time=0.0_dp, &
      state=state, rank=mpi%rank)
    call gp3d_mpi_barrier(mpi)
  end if

  call gp3d_gpu_init(gpu, grid, mpi, &
    model_cfg%argle_enabled .and. .not. restart_info%loaded)
  call gp3d_gpu_upload(gpu, state, grid)
  if (.not. restart_info%loaded) then
    call gp3d_gpu_relax_taylor_green_argle(gpu, params, model_cfg)
  end if

  if (gp3d_mpi_is_root(mpi)) then
    if (run_cfg%write_meta) call gp3d_write_meta_json(output_cfg)
  end if
  call gp3d_mpi_barrier(mpi)
  if (run_cfg%write_initial) then
    call gp3d_gpu_download(gpu, state)
    call gp3d_write_gpe_psi_slf(output_cfg, step=start_step, &
      time=start_time, state=state, rank=mpi%rank)
  end if
  call gp3d_mpi_barrier(mpi)

  call gp3d_gpu_diagnostics(gpu, grid, params, norm_value, energy_value)
  if (gp3d_mpi_is_root(mpi)) then
    write(*,'(a)') "# step norm energy"
    write(*,'(i8,1x,2(es24.16,1x))') start_step, norm_value, energy_value
  end if

  ! 出力のないステップでは psi を GPU から戻さない。
  do advance_step = 1, params%nsteps
    step = start_step + advance_step
    current_time = start_time + real(advance_step, dp) * params%dt
    output_step = .false.
    if (params%output_every > 0) output_step = mod(step, params%output_every) == 0
    if (run_cfg%timing_enabled) then
      call gp3d_gpu_step(gpu, grid, params, step_timing)
    else
      call gp3d_gpu_step(gpu, grid, params)
    end if

    if (output_step) then
      call gp3d_gpu_diagnostics(gpu, grid, params, norm_value, energy_value)
      if (gp3d_mpi_is_root(mpi)) then
        write(*,'(i8,1x,2(es24.16,1x))') step, norm_value, energy_value
      end if
      call gp3d_gpu_download(gpu, state)
      call gp3d_write_gpe_psi_slf(output_cfg, step, current_time, state, &
        rank=mpi%rank)
      call gp3d_mpi_barrier(mpi)
    end if
  end do

  if (run_cfg%timing_enabled) call gp3d_report_step_timing(step_timing, mpi)

  ! cuFFTMp/NVSHMEM の解放は必ず MPI_Finalize より先に行う。
  call gp3d_gpu_finalize(gpu)
  call gp3d_mpi_finalize(mpi)
end program gp3d_cufftmp
