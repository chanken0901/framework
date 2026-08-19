!> 単一GPU用GP3D実行プログラム。
!> 初期化後にpsiをGPUへ置き、ARGLEと実時間発展をGPU常駐のまま進め、出力時だけ取得する。
program gp3d_cuda
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
  use gp3d_mpi, only: gp3d_mpi_t, gp3d_mpi_init, gp3d_mpi_finalize
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

  ! CUDA版はMPIスタブを使い、常に単一rankとして初期化する。
  call gp3d_mpi_init(mpi)

  input_file = "input.nml"
  if (command_argument_count() >= 1) call get_command_argument(1, input_file)

  ! namelistを読み、バックエンド設定を単一CUDA実行へ固定する。
  call gp3d_read_all_inputs(trim(input_file), run_cfg, model_cfg)
  run_cfg%use_mpi = .false.
  run_cfg%use_cuda = .true.
  run_cfg%rank = 0
  run_cfg%nprocs = 1

  ! ホスト側で格子・初期条件、または再スタート波動関数を準備する。
  call gp3d_configure_problem(run_cfg, model_cfg, grid, params, state, output_cfg, mpi, restart_info)
  start_step = restart_info%step
  start_time = restart_info%time
  if (restart_info%loaded) then
    write(*,'(a,a)') "# restart source=", trim(restart_info%source_file)
    write(*,'(a,i0,1x,a,es16.8,1x,a,i0)') "# restart step=", start_step, &
      "time=", start_time, "source_nprocs=", restart_info%source_nprocs
  end if

  ! GPUへ転送する前の解析的Taylor-Green種を必要なら保存する。
  if (.not. restart_info%loaded .and. model_cfg%argle_enabled .and. model_cfg%argle_write_seed) then
    seed_output_cfg = output_cfg
    seed_output_cfg%output_dir = trim(output_cfg%output_dir) // "/argle_seed"
    seed_output_cfg%case_name = trim(output_cfg%case_name) // "_argle_seed"
    call gp3d_write_meta_json(seed_output_cfg)
    call gp3d_write_gpe_psi_slf(seed_output_cfg, step=0, time=0.0_dp, state=state)
  end if

  ! ここでGPU常駐領域を確保し、以後は出力時以外psiをホストへ戻さない。
  call gp3d_gpu_init(gpu, grid, model_cfg%argle_enabled .and. .not. restart_info%loaded, &
    run_cfg%cuda_device)
  call gp3d_gpu_upload(gpu, state, grid)
  if (.not. restart_info%loaded) call gp3d_gpu_relax_taylor_green_argle(gpu, params, model_cfg)

  if (run_cfg%write_meta) call gp3d_write_meta_json(output_cfg)
  if (run_cfg%write_initial) then
    call gp3d_gpu_download(gpu, state)
    call gp3d_write_gpe_psi_slf(output_cfg, step=start_step, time=start_time, state=state)
  end if

  call gp3d_gpu_diagnostics(gpu, grid, params, norm_value, energy_value)
  write(*,'(a)') "# step norm energy"
  write(*,'(i8,1x,2(es24.16,1x))') start_step, norm_value, energy_value

  ! cuFFTを含むSplit-stepをGPU上で繰り返し、出力ステップだけダウンロードする。
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
      write(*,'(i8,1x,2(es24.16,1x))') step, norm_value, energy_value
      call gp3d_gpu_download(gpu, state)
      call gp3d_write_gpe_psi_slf(output_cfg, step, current_time, state)
    end if
  end do

  if (run_cfg%timing_enabled) call gp3d_report_step_timing(step_timing, mpi)

  ! cuFFT planと全GPUメモリをまとめて解放する。
  call gp3d_gpu_finalize(gpu)
  call gp3d_mpi_finalize(mpi)
end program gp3d_cuda
