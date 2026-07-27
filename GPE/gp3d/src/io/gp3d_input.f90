!> Fortran namelistを読み、格子・物理量・初期波動関数・出力設定を組み立てる。
!> restart指定時はSLFを読み、それ以外はflow名から初期条件を選択して正規化する。
module gp3d_input
  use gp3d_types, only: dp, pi, gp3d_run_config_t, gp3d_model_config_t, gp3d_params_t, &
    gp3d_grid_t, gp3d_state_t
  use gp3d_grid, only: gp3d_grid_init_bounds
  use gp3d_solver, only: gp3d_state_allocate, gp3d_set_harmonic_potential, &
    gp3d_set_gaussian_initial_state, gp3d_normalize
  use gp3d_initial_conditions, only: gp3d_set_uniform_vortex_line, &
    gp3d_set_thomas_fermi_vortex_line, gp3d_imprint_vortex_ring, &
    gp3d_imprint_random_vortex_lines, gp3d_imprint_random_vortex_rings, &
    gp3d_add_random_phase_noise, gp3d_set_quantum_taylor_green
  use gp3d_io, only: gp3d_output_config_t, gp3d_output_config_from_grid
  use gp3d_restart, only: gp3d_restart_info_t, gp3d_restart_requested, gp3d_restart_load
  use gp3d_mpi, only: gp3d_mpi_t
  implicit none
  private

  public :: gp3d_read_common_input
  public :: gp3d_read_gpe_input
  public :: gp3d_read_all_inputs
  public :: gp3d_update_derived_config
  public :: gp3d_configure_problem

contains

  subroutine gp3d_read_common_input(filename, cfg)
    character(len=*), intent(in) :: filename
    type(gp3d_run_config_t), intent(inout) :: cfg

    character(len=32) :: equation, output_format, backend, precision_name
    character(len=256) :: case_name, input_file, output_dir, restart_file
    character(len=64) :: initial_condition
    integer :: nx, ny, nz, nghost, nsteps, output_frequency
    integer :: cuda_device, rank, nprocs
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp) :: dt, t_max, cfl
    logical :: use_fixed_dt, write_initial, write_meta, timing_enabled
    logical :: use_mpi, use_openmp, use_cuda
    integer :: unit, ios
    logical :: exists

    namelist /simulation/ equation, case_name, input_file, initial_condition, restart_file, &
      nx, ny, nz, nghost, &
      x_min, x_max, y_min, y_max, z_min, z_max, &
      dt, t_max, nsteps, cfl, use_fixed_dt, &
      output_frequency, output_dir, output_format, precision_name, &
      write_initial, write_meta, timing_enabled, backend, &
      use_mpi, use_openmp, use_cuda, cuda_device, rank, nprocs

    equation = cfg%equation
    case_name = cfg%case_name
    input_file = filename
    initial_condition = cfg%initial_condition
    restart_file = cfg%restart_file
    nx = cfg%nx
    ny = cfg%ny
    nz = cfg%nz
    nghost = cfg%nghost
    x_min = cfg%x_min
    x_max = cfg%x_max
    y_min = cfg%y_min
    y_max = cfg%y_max
    z_min = cfg%z_min
    z_max = cfg%z_max
    dt = cfg%dt
    t_max = cfg%t_max
    nsteps = cfg%nsteps
    cfl = cfg%cfl
    use_fixed_dt = cfg%use_fixed_dt
    output_frequency = cfg%output_frequency
    output_dir = cfg%output_dir
    output_format = cfg%output_format
    precision_name = cfg%precision_name
    write_initial = cfg%write_initial
    write_meta = cfg%write_meta
    timing_enabled = cfg%timing_enabled
    backend = cfg%backend
    use_mpi = cfg%use_mpi
    use_openmp = cfg%use_openmp
    use_cuda = cfg%use_cuda
    cuda_device = cfg%cuda_device
    rank = cfg%rank
    nprocs = cfg%nprocs

    inquire(file=filename, exist=exists)
    if (.not. exists) then
      write(*,'(A,A)') "WARNING: input file not found. Use defaults: ", trim(filename)
      call gp3d_update_derived_config(cfg)
      return
    end if

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "ERROR: cannot open input file."
    read(unit, nml=simulation, iostat=ios)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') "ERROR: failed to read namelist /simulation/ in ", trim(filename)
      error stop
    end if

    cfg%equation = equation
    cfg%case_name = case_name
    cfg%input_file = input_file
    cfg%initial_condition = initial_condition
    cfg%restart_file = restart_file
    cfg%nx = nx
    cfg%ny = ny
    cfg%nz = nz
    cfg%nghost = nghost
    cfg%x_min = x_min
    cfg%x_max = x_max
    cfg%y_min = y_min
    cfg%y_max = y_max
    cfg%z_min = z_min
    cfg%z_max = z_max
    cfg%dt = dt
    cfg%t_max = t_max
    cfg%nsteps = nsteps
    cfg%cfl = cfl
    cfg%use_fixed_dt = use_fixed_dt
    cfg%output_frequency = output_frequency
    cfg%output_dir = output_dir
    cfg%output_format = output_format
    cfg%precision_name = precision_name
    cfg%write_initial = write_initial
    cfg%write_meta = write_meta
    cfg%timing_enabled = timing_enabled
    cfg%backend = backend
    cfg%use_mpi = use_mpi
    cfg%use_openmp = use_openmp
    cfg%use_cuda = use_cuda
    cfg%cuda_device = cuda_device
    cfg%rank = rank
    cfg%nprocs = nprocs

    call gp3d_update_derived_config(cfg)
  end subroutine gp3d_read_common_input

  subroutine gp3d_read_gpe_input(filename, cfg)
    character(len=*), intent(in) :: filename
    type(gp3d_model_config_t), intent(inout) :: cfg

    real(dp) :: alpha, beta, g, sigma0, wx, wy, wz, hbar, mass, norm, mu, healing_length
    real(dp) :: vortex_x0, vortex_y0, ring_radius, ring_z0, phase_noise
    real(dp) :: ring_radius_min, ring_radius_max, tg_velocity_amplitude
    real(dp) :: argle_dtau, argle_tolerance
    integer :: vortex_charge, tangle_nlines, tangle_nrings, random_seed, tg_winding
    integer :: argle_steps, argle_output_every
    logical :: imaginary_time, argle_enabled, argle_write_seed
    logical :: use_dimensionless_parameters, tg_auto_winding
    integer :: unit, ios
    logical :: exists

    namelist /gpe/ use_dimensionless_parameters, alpha, beta, &
      g, sigma0, wx, wy, wz, hbar, mass, norm, imaginary_time, &
      mu, healing_length, vortex_charge, vortex_x0, vortex_y0, ring_radius, ring_z0, phase_noise, &
      tangle_nlines, tangle_nrings, ring_radius_min, ring_radius_max, random_seed, &
      tg_velocity_amplitude, tg_winding, tg_auto_winding, &
      argle_enabled, argle_write_seed, argle_steps, argle_output_every, &
      argle_dtau, argle_tolerance

    use_dimensionless_parameters = cfg%use_dimensionless_parameters
    alpha = cfg%alpha
    beta = cfg%beta
    g = cfg%g
    sigma0 = cfg%sigma0
    wx = cfg%wx
    wy = cfg%wy
    wz = cfg%wz
    hbar = cfg%hbar
    mass = cfg%mass
    norm = cfg%norm
    imaginary_time = cfg%imaginary_time
    mu = cfg%mu
    healing_length = cfg%healing_length
    vortex_charge = cfg%vortex_charge
    vortex_x0 = cfg%vortex_x0
    vortex_y0 = cfg%vortex_y0
    ring_radius = cfg%ring_radius
    ring_z0 = cfg%ring_z0
    phase_noise = cfg%phase_noise
    tangle_nlines = cfg%tangle_nlines
    tangle_nrings = cfg%tangle_nrings
    ring_radius_min = cfg%ring_radius_min
    ring_radius_max = cfg%ring_radius_max
    random_seed = cfg%random_seed
    tg_velocity_amplitude = cfg%tg_velocity_amplitude
    tg_winding = cfg%tg_winding
    tg_auto_winding = cfg%tg_auto_winding
    argle_enabled = cfg%argle_enabled
    argle_write_seed = cfg%argle_write_seed
    argle_steps = cfg%argle_steps
    argle_output_every = cfg%argle_output_every
    argle_dtau = cfg%argle_dtau
    argle_tolerance = cfg%argle_tolerance

    inquire(file=filename, exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "ERROR: cannot open input file."
    read(unit, nml=gpe, iostat=ios)
    close(unit)
    if (ios /= 0) return

    cfg%use_dimensionless_parameters = use_dimensionless_parameters
    cfg%alpha = alpha
    cfg%beta = beta
    cfg%g = g
    cfg%sigma0 = sigma0
    cfg%wx = wx
    cfg%wy = wy
    cfg%wz = wz
    cfg%hbar = hbar
    cfg%mass = mass
    cfg%norm = norm
    cfg%imaginary_time = imaginary_time
    cfg%mu = mu
    cfg%healing_length = healing_length
    cfg%vortex_charge = vortex_charge
    cfg%vortex_x0 = vortex_x0
    cfg%vortex_y0 = vortex_y0
    cfg%ring_radius = ring_radius
    cfg%ring_z0 = ring_z0
    cfg%phase_noise = phase_noise
    cfg%tangle_nlines = tangle_nlines
    cfg%tangle_nrings = tangle_nrings
    cfg%ring_radius_min = ring_radius_min
    cfg%ring_radius_max = ring_radius_max
    cfg%random_seed = random_seed
    cfg%tg_velocity_amplitude = tg_velocity_amplitude
    cfg%tg_winding = tg_winding
    cfg%tg_auto_winding = tg_auto_winding
    cfg%argle_enabled = argle_enabled
    cfg%argle_write_seed = argle_write_seed
    cfg%argle_steps = argle_steps
    cfg%argle_output_every = argle_output_every
    cfg%argle_dtau = argle_dtau
    cfg%argle_tolerance = argle_tolerance
  end subroutine gp3d_read_gpe_input

  subroutine gp3d_read_all_inputs(filename, run_cfg, model_cfg)
    character(len=*), intent(in) :: filename
    type(gp3d_run_config_t), intent(inout) :: run_cfg
    type(gp3d_model_config_t), intent(inout) :: model_cfg

    call gp3d_read_common_input(filename, run_cfg)
    call gp3d_read_gpe_input(filename, model_cfg)
  end subroutine gp3d_read_all_inputs

  subroutine gp3d_update_derived_config(cfg)
    type(gp3d_run_config_t), intent(inout) :: cfg

    if (cfg%x_max <= cfg%x_min .or. cfg%y_max <= cfg%y_min .or. cfg%z_max <= cfg%z_min) then
      error stop "invalid simulation bounds"
    end if
    if (cfg%nx <= 0 .or. cfg%ny <= 0 .or. cfg%nz <= 0) error stop "grid dimensions must be positive"

    cfg%lx = cfg%x_max - cfg%x_min
    cfg%ly = cfg%y_max - cfg%y_min
    cfg%lz = cfg%z_max - cfg%z_min
    cfg%dx = cfg%lx / real(cfg%nx, dp)
    cfg%dy = cfg%ly / real(cfg%ny, dp)
    cfg%dz = cfg%lz / real(cfg%nz, dp)

    if (cfg%nsteps <= 0 .and. cfg%t_max > 0.0_dp .and. cfg%dt > 0.0_dp) then
      cfg%nsteps = ceiling(cfg%t_max / cfg%dt)
    end if
  end subroutine gp3d_update_derived_config

  subroutine gp3d_configure_problem(run_cfg, model_cfg, grid, params, state, output_cfg, mpi, restart_info)
    ! 読み込み済み設定から、ソルバーが直ちに使える全オブジェクトを一括初期化する。
    type(gp3d_run_config_t), intent(in) :: run_cfg
    type(gp3d_model_config_t), intent(in) :: model_cfg
    type(gp3d_grid_t), intent(out) :: grid
    type(gp3d_params_t), intent(out) :: params
    type(gp3d_state_t), intent(out) :: state
    type(gp3d_output_config_t), intent(out) :: output_cfg
    type(gp3d_mpi_t), intent(in) :: mpi
    type(gp3d_restart_info_t), intent(out), optional :: restart_info

    type(gp3d_restart_info_t) :: loaded_restart

    call gp3d_grid_init_bounds(grid, run_cfg%nx, run_cfg%ny, run_cfg%nz, &
      run_cfg%x_min, run_cfg%x_max, run_cfg%y_min, run_cfg%y_max, run_cfg%z_min, run_cfg%z_max, &
      run_cfg%rank, run_cfg%nprocs)

    params%dt = run_cfg%dt
    params%nsteps = run_cfg%nsteps
    params%output_every = run_cfg%output_frequency
    if (model_cfg%use_dimensionless_parameters) then
      if (model_cfg%alpha <= 0.0_dp .or. model_cfg%beta <= 0.0_dp) then
        error stop "dimensionless alpha and beta must be positive"
      end if
      params%hbar = 1.0_dp
      params%mass = 1.0_dp / (2.0_dp * model_cfg%alpha)
      params%g = model_cfg%beta
    else
      if (model_cfg%hbar <= 0.0_dp .or. model_cfg%mass <= 0.0_dp) then
        error stop "hbar and mass must be positive"
      end if
      params%g = model_cfg%g
      params%mass = model_cfg%mass
      params%hbar = model_cfg%hbar
    end if
    params%norm = model_cfg%norm
    params%imaginary_time = model_cfg%imaginary_time

    call gp3d_state_allocate(state, grid)
    call gp3d_set_harmonic_potential(state, grid, model_cfg%wx, model_cfg%wy, model_cfg%wz)
    loaded_restart = gp3d_restart_info_t()
    if (gp3d_restart_requested(run_cfg)) then
      call gp3d_restart_load(trim(run_cfg%restart_file), state, grid, loaded_restart)
    else
      call set_initial_condition(trim(run_cfg%initial_condition), state, grid, model_cfg)
      if (model_cfg%phase_noise > 0.0_dp) then
        call gp3d_add_random_phase_noise(state, model_cfg%phase_noise, grid, model_cfg%random_seed)
      end if
      if (.not. (is_taylor_green_name(run_cfg%initial_condition) .and. model_cfg%argle_enabled)) then
        call gp3d_normalize(state, grid, params%norm, mpi)
      end if
    end if

    if (.not. loaded_restart%loaded .and. is_taylor_green_name(run_cfg%initial_condition) .and. &
        mpi%rank == mpi%root) then
      call report_taylor_green_parameters(model_cfg, grid)
    end if

    call gp3d_output_config_from_grid(output_cfg, grid, run_cfg%output_dir, run_cfg%case_name)
    output_cfg%equation = run_cfg%equation
    output_cfg%output_format = run_cfg%output_format
    output_cfg%precision_name = run_cfg%precision_name
    output_cfg%nghost = run_cfg%nghost
    output_cfg%use_mpi = run_cfg%use_mpi
    output_cfg%use_openmp = run_cfg%use_openmp
    output_cfg%use_cuda = run_cfg%use_cuda
    output_cfg%cuda_device = run_cfg%cuda_device
    output_cfg%mpi_rank = run_cfg%rank
    output_cfg%mpi_nprocs = run_cfg%nprocs
    if (present(restart_info)) restart_info = loaded_restart
  end subroutine gp3d_configure_problem

  subroutine set_initial_condition(name, state, grid, cfg)
    ! initial_condition文字列を具体的な初期条件生成ルーチンへ振り分ける。
    character(len=*), intent(in) :: name
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_model_config_t), intent(in) :: cfg

    real(dp) :: tg_healing_length, tg_alpha, tg_beta
    integer :: tg_winding

    select case (trim(name))
    case ("gaussian", "Gaussian", "GAUSSIAN")
      call gp3d_set_gaussian_initial_state(state, grid, cfg%sigma0, cfg%sigma0, cfg%sigma0)
    case ("uniform_vortex", "vortex_uniform")
      call gp3d_set_uniform_vortex_line(state, grid, cfg%vortex_charge, cfg%vortex_x0, cfg%vortex_y0, &
        cfg%healing_length, density0=1.0_dp)
    case ("tf_vortex", "thomas_fermi_vortex", "vortex")
      call gp3d_set_thomas_fermi_vortex_line(state, grid, cfg%vortex_charge, cfg%vortex_x0, cfg%vortex_y0, &
        cfg%healing_length, cfg%mu, cfg%g)
    case ("vortex_ring", "ring")
      call gp3d_set_gaussian_initial_state(state, grid, cfg%sigma0, cfg%sigma0, cfg%sigma0)
      call gp3d_imprint_vortex_ring(state, grid, cfg%vortex_charge, cfg%ring_radius, cfg%ring_z0, &
        cfg%healing_length)
    case ("vortex_tangle", "random_vortices", "quantum_turbulence")
      call gp3d_set_gaussian_initial_state(state, grid, cfg%sigma0, cfg%sigma0, cfg%sigma0)
      call gp3d_imprint_random_vortex_lines(state, grid, cfg%tangle_nlines, cfg%healing_length, cfg%random_seed)
    case ("ring_tangle", "vortex_ring_tangle", "random_rings")
      state%psi = cmplx(1.0_dp, 0.0_dp, kind=dp)
      call gp3d_imprint_random_vortex_rings(state, grid, cfg%tangle_nrings, cfg%ring_radius_min, &
        cfg%ring_radius_max, cfg%healing_length, cfg%random_seed)
    case ("quantum_taylor_green", "taylor_green", "tg")
      call resolve_taylor_green_parameters(cfg, tg_alpha, tg_beta, tg_healing_length, tg_winding)
      call gp3d_set_quantum_taylor_green(state, grid, tg_healing_length, tg_winding, density0=1.0_dp)
    case ("restart_slf", "restart", "slf")
      error stop "restart initial condition requires restart_file in /simulation/"
    case default
      write(*,'(A,A,A)') "WARNING: unknown initial_condition='", trim(name), "'. Use gaussian."
      call gp3d_set_gaussian_initial_state(state, grid, cfg%sigma0, cfg%sigma0, cfg%sigma0)
    end select
  end subroutine set_initial_condition

  subroutine resolve_taylor_green_parameters(cfg, alpha, beta, healing_length, winding)
    ! 無次元alpha/betaから渦芯幅と循環巻き数を導出し、設定の物理整合性を確認する。
    type(gp3d_model_config_t), intent(in) :: cfg
    real(dp), intent(out) :: alpha, beta, healing_length
    integer, intent(out) :: winding

    if (cfg%use_dimensionless_parameters) then
      alpha = cfg%alpha
      beta = cfg%beta
      if (alpha <= 0.0_dp .or. beta <= 0.0_dp) then
        error stop "Taylor-Green alpha and beta must be positive"
      end if
      healing_length = sqrt(alpha / beta)
    else
      if (cfg%hbar <= 0.0_dp .or. cfg%mass <= 0.0_dp .or. cfg%g <= 0.0_dp) then
        error stop "Taylor-Green hbar, mass, and g must be positive"
      end if
      alpha = cfg%hbar / (2.0_dp * cfg%mass)
      beta = cfg%g / cfg%hbar
      healing_length = cfg%healing_length
    end if

    if (cfg%tg_velocity_amplitude <= 0.0_dp) then
      error stop "Taylor-Green velocity amplitude must be positive"
    end if
    if (cfg%tg_auto_winding) then
      winding = floor(cfg%tg_velocity_amplitude / (2.0_dp * pi * alpha))
    else
      winding = cfg%tg_winding
    end if
    if (winding < 1) then
      error stop "Taylor-Green parameters give zero circulation winding; reduce alpha or increase velocity"
    end if
  end subroutine resolve_taylor_green_parameters

  subroutine report_taylor_green_parameters(cfg, grid)
    type(gp3d_model_config_t), intent(in) :: cfg
    type(gp3d_grid_t), intent(in) :: grid

    real(dp) :: alpha, beta, healing_length, sound_speed, mach_number
    integer :: winding

    call resolve_taylor_green_parameters(cfg, alpha, beta, healing_length, winding)
    sound_speed = sqrt(2.0_dp * alpha * beta)
    mach_number = cfg%tg_velocity_amplitude / sound_speed
    write(*,'(a)') "# Taylor-Green derived parameters"
    write(*,'(a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
      "# alpha=", alpha, "beta=", beta, "xi=", healing_length
    write(*,'(a,es16.8,1x,a,es16.8,1x,a,i0,1x,a,es16.8)') &
      "# sound_speed=", sound_speed, "Mach=", mach_number, "winding=", winding, &
      "dx/xi=", grid%dx / healing_length
  end subroutine report_taylor_green_parameters

  pure logical function is_taylor_green_name(name) result(is_tg)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("quantum_taylor_green", "taylor_green", "tg")
      is_tg = .true.
    case default
      is_tg = .false.
    end select
  end function is_taylor_green_name

end module gp3d_input
