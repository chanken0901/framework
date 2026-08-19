!> CPU逐次版とMPI版で共有するGP3D実行プログラム。
!> 入力、初期化/再開、任意のARGLE、実時間Split-step、SLF出力を順番に制御する。
program gp3d_sequential
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_params_t, gp3d_state_t, &
    gp3d_run_config_t, gp3d_model_config_t
  use gp3d_fft, only: gp3d_fft_plan_t, gp3d_fft_init, gp3d_fft_finalize
  use gp3d_solver, only: gp3d_density_norm, gp3d_energy, gp3d_step_split_operator, &
    gp3d_relax_taylor_green_argle, gp3d_step_timing_t, gp3d_report_step_timing
  use gp3d_io, only: gp3d_output_config_t, gp3d_write_meta_json, gp3d_write_gpe_psi_slf
  use gp3d_input, only: gp3d_read_all_inputs, gp3d_configure_problem
  use gp3d_restart, only: gp3d_restart_info_t
  use gp3d_openmp, only: gp3d_openmp_configure, gp3d_openmp_is_compiled, &
    gp3d_openmp_active, gp3d_openmp_thread_count
  use gp3d_mpi, only: gp3d_mpi_t, gp3d_mpi_init, gp3d_mpi_finalize, &
    gp3d_mpi_is_root, gp3d_mpi_barrier, gp3d_mpi_supports_funneled
  implicit none

  type(gp3d_mpi_t) :: mpi
  type(gp3d_run_config_t) :: run_cfg
  type(gp3d_model_config_t) :: model_cfg
  type(gp3d_grid_t) :: grid
  type(gp3d_params_t) :: params
  type(gp3d_state_t) :: state
  type(gp3d_fft_plan_t) :: fft_plan
  type(gp3d_output_config_t) :: output_cfg, seed_output_cfg
  type(gp3d_step_timing_t) :: step_timing
  type(gp3d_restart_info_t) :: restart_info
  integer :: step, advance_step, start_step
  real(dp) :: norm_value, energy_value, start_time, current_time
  logical :: output_step
  character(len=256) :: input_file

  ! 実MPI版ではrank/sizeを取得し、逐次スタブではrank=0, size=1を設定する。
  call gp3d_mpi_init(mpi)

  input_file = "input.nml"
  if (command_argument_count() >= 1) call get_command_argument(1, input_file)

  ! namelistを読み、実際のMPI状態を実行設定へ反映する。
  call gp3d_read_all_inputs(trim(input_file), run_cfg, model_cfg)
  run_cfg%use_mpi = mpi%enabled
  run_cfg%rank = mpi%rank
  run_cfg%nprocs = mpi%nprocs
  if (run_cfg%use_openmp .and. .not. gp3d_mpi_supports_funneled(mpi)) then
    error stop "MPI implementation does not provide MPI_THREAD_FUNNELED for OpenMP"
  end if
  call gp3d_openmp_configure(run_cfg%use_openmp)
  if (gp3d_mpi_is_root(mpi)) then
    write(*,'(a,l1,a,l1,a,i0)') "# OpenMP compiled=", gp3d_openmp_is_compiled(), &
      " active=", gp3d_openmp_active, " threads_per_rank=", gp3d_openmp_thread_count()
  end if

  ! 格子・ポテンシャル・初期条件、または再スタート波動関数を準備する。
  call gp3d_configure_problem(run_cfg, model_cfg, grid, params, state, output_cfg, mpi, restart_info)
  start_step = restart_info%step
  start_time = restart_info%time
  if (restart_info%loaded .and. gp3d_mpi_is_root(mpi)) then
    write(*,'(a,a)') "# restart source=", trim(restart_info%source_file)
    write(*,'(a,i0,1x,a,es16.8,1x,a,i0)') "# restart step=", start_step, &
      "time=", start_time, "source_nprocs=", restart_info%source_nprocs
  end if
  call gp3d_fft_init(fft_plan, grid%nx, grid%ny, grid%nz, mpi%comm, mpi%rank, mpi%nprocs)

  ! ARGLE前の解析的Taylor-Green種を、比較・再現用に必要なら保存する。
  if (.not. restart_info%loaded .and. model_cfg%argle_enabled .and. model_cfg%argle_write_seed) then
    seed_output_cfg = output_cfg
    seed_output_cfg%output_dir = trim(output_cfg%output_dir) // "/argle_seed"
    seed_output_cfg%case_name = trim(output_cfg%case_name) // "_argle_seed"
    if (gp3d_mpi_is_root(mpi)) call gp3d_write_meta_json(seed_output_cfg)
    call gp3d_mpi_barrier(mpi)
    if (mpi%enabled) then
      call gp3d_write_gpe_psi_slf(seed_output_cfg, step=0, time=0.0_dp, state=state, rank=mpi%rank)
    else
      call gp3d_write_gpe_psi_slf(seed_output_cfg, step=0, time=0.0_dp, state=state)
    end if
    call gp3d_mpi_barrier(mpi)
  end if

  ! 再スタートでない場合だけ、実時間発展前のARGLE緩和を行う。
  if (.not. restart_info%loaded) then
    call gp3d_relax_taylor_green_argle(state, grid, params, model_cfg, fft_plan, mpi)
  end if

  if (gp3d_mpi_is_root(mpi)) then
    if (run_cfg%write_meta) call gp3d_write_meta_json(output_cfg)
  end if
  call gp3d_mpi_barrier(mpi)
  if (run_cfg%write_initial) then
    if (mpi%enabled) then
      call gp3d_write_gpe_psi_slf(output_cfg, step=start_step, time=start_time, state=state, rank=mpi%rank)
    else
      call gp3d_write_gpe_psi_slf(output_cfg, step=start_step, time=start_time, state=state)
    end if
  end if
  call gp3d_mpi_barrier(mpi)

  norm_value = gp3d_density_norm(state, grid, mpi)
  energy_value = gp3d_energy(state, grid, params, fft_plan, mpi)
  if (gp3d_mpi_is_root(mpi)) then
    write(*,'(a)') "# step norm energy"
    write(*,'(i8,1x,2(es24.16,1x))') start_step, norm_value, energy_value
  end if

  ! 実時間または通常の虚時間Split-operatorを進め、指定間隔で診断・SLF出力する。
  do advance_step = 1, params%nsteps
    step = start_step + advance_step
    current_time = start_time + real(advance_step, dp) * params%dt
    output_step = .false.
    if (params%output_every > 0) output_step = mod(step, params%output_every) == 0
    if (run_cfg%timing_enabled) then
      call gp3d_step_split_operator(state, grid, params, fft_plan, mpi, step_timing)
    else
      call gp3d_step_split_operator(state, grid, params, fft_plan, mpi)
    end if

    if (output_step) then
      norm_value = gp3d_density_norm(state, grid, mpi)
      energy_value = gp3d_energy(state, grid, params, fft_plan, mpi)
      if (gp3d_mpi_is_root(mpi)) then
        write(*,'(i8,1x,2(es24.16,1x))') step, norm_value, energy_value
      end if
      if (mpi%enabled) then
        call gp3d_write_gpe_psi_slf(output_cfg, step, current_time, state, rank=mpi%rank)
      else
        call gp3d_write_gpe_psi_slf(output_cfg, step, current_time, state)
      end if
    end if
    if (output_step) call gp3d_mpi_barrier(mpi)
  end do

  if (run_cfg%timing_enabled) call gp3d_report_step_timing(step_timing, mpi)

  ! FFT planとMPI環境は生成と逆順で解放する。
  call gp3d_fft_finalize(fft_plan)
  call gp3d_mpi_finalize(mpi)
end program gp3d_sequential
