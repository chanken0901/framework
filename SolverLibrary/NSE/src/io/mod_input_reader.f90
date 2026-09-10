module mod_input_reader
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, update_derived_config
  use mod_model_config, only : gpe_config, nse_config
  implicit none
  private

  public :: read_common_input
  public :: read_gpe_input
  public :: read_nse_input
  public :: read_all_inputs

contains

  subroutine read_common_input(filename, cfg)
    character(len=*), intent(in) :: filename
    type(simulation_config), intent(inout) :: cfg

    character(len=32)  :: equation, output_format, backend, precision_name
    character(len=256) :: case_name, input_file, output_dir
    character(len=64)  :: initial_condition
    integer :: nx, ny, nz, nghost, nsteps, output_frequency
    integer :: rank, nprocs, cuda_device
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp) :: dt, t_max, cfl
    logical :: use_fixed_dt, write_initial, write_meta, use_mpi, use_openmp
    integer :: u, ios
    logical :: exists

    namelist /simulation/ equation, case_name, input_file, initial_condition, &
      nx, ny, nz, nghost, &
      x_min, x_max, y_min, y_max, z_min, z_max, &
      dt, t_max, nsteps, cfl, use_fixed_dt, &
      output_frequency, output_dir, output_format, precision_name, &
      write_initial, write_meta, backend, use_mpi, use_openmp, rank, nprocs, &
      cuda_device

    ! copy defaults from cfg
    equation = cfg%equation
    case_name = cfg%case_name
    input_file = filename
    initial_condition = cfg%initial_condition
    nx = cfg%nx; ny = cfg%ny; nz = cfg%nz
    nghost = cfg%nghost
    x_min = cfg%x_min; x_max = cfg%x_max
    y_min = cfg%y_min; y_max = cfg%y_max
    z_min = cfg%z_min; z_max = cfg%z_max
    dt = cfg%dt; t_max = cfg%t_max; nsteps = cfg%nsteps; cfl = cfg%cfl
    use_fixed_dt = cfg%use_fixed_dt
    output_frequency = cfg%output_frequency
    output_dir = cfg%output_dir
    output_format = cfg%output_format
    precision_name = cfg%precision_name
    write_initial = cfg%write_initial
    write_meta = cfg%write_meta
    backend = cfg%backend
    use_mpi = cfg%use_mpi
    use_openmp = cfg%use_openmp
    rank = cfg%rank
    nprocs = cfg%nprocs
    cuda_device = cfg%cuda_device

    inquire(file=filename, exist=exists)
    if (.not. exists) then
      write(*,'(A,A)') 'WARNING: input file not found. Use defaults: ', trim(filename)
      call update_derived_config(cfg)
      return
    end if

    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot open input file.'

    read(u, nml=simulation, iostat=ios)
    close(u)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: failed to read namelist /simulation/ in ', trim(filename)
      error stop
    end if

    cfg%equation = equation
    cfg%case_name = case_name
    cfg%input_file = input_file
    cfg%initial_condition = initial_condition
    cfg%nx = nx; cfg%ny = ny; cfg%nz = nz
    cfg%nghost = nghost
    cfg%x_min = x_min; cfg%x_max = x_max
    cfg%y_min = y_min; cfg%y_max = y_max
    cfg%z_min = z_min; cfg%z_max = z_max
    cfg%dt = dt; cfg%t_max = t_max; cfg%nsteps = nsteps; cfg%cfl = cfl
    cfg%use_fixed_dt = use_fixed_dt
    cfg%output_frequency = output_frequency
    cfg%output_dir = output_dir
    cfg%output_format = output_format
    cfg%precision_name = precision_name
    cfg%write_initial = write_initial
    cfg%write_meta = write_meta
    cfg%backend = backend
    cfg%use_mpi = use_mpi
    cfg%use_openmp = use_openmp
    cfg%rank = rank
    cfg%nprocs = nprocs
    cfg%cuda_device = cuda_device

    call update_derived_config(cfg)
  end subroutine read_common_input

  subroutine read_gpe_input(filename, cfg)
    character(len=*), intent(in) :: filename
    type(gpe_config), intent(inout) :: cfg

    real(dp) :: g, sigma0, wx, wy, wz, hbar, mass
    integer :: u, ios
    logical :: exists
    namelist /gpe/ g, sigma0, wx, wy, wz, hbar, mass

    g = cfg%g; sigma0 = cfg%sigma0
    wx = cfg%wx; wy = cfg%wy; wz = cfg%wz
    hbar = cfg%hbar; mass = cfg%mass

    inquire(file=filename, exist=exists)
    if (.not. exists) return
    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot open input file.'
    read(u, nml=gpe, iostat=ios)
    close(u)
    if (ios /= 0) return  ! /gpe/ block is optional for non-GPE solver

    cfg%g = g; cfg%sigma0 = sigma0
    cfg%wx = wx; cfg%wy = wy; cfg%wz = wz
    cfg%hbar = hbar; cfg%mass = mass
  end subroutine read_gpe_input

  subroutine read_nse_input(filename, cfg)
    character(len=*), intent(in) :: filename
    type(nse_config), intent(inout) :: cfg

    integer :: nv, nghost
    logical :: fh_enabled
    real(dp) :: fh_boltzmann_number
    integer :: fh_seed
    real(dp) :: gamma, cfl, small_rho, small_p, rho0, mach, reynolds, prandtl
    real(dp) :: hit_turbulent_mach, hit_turbulent_reynolds
    real(dp) :: hit_rms_velocity, hit_peak_wavenumber
    real(dp) :: hit_integral_length, hit_kolmogorov_length
    real(dp) :: hit_johnsen_length_scale_ratio
    real(dp) :: hit_pope_energy_constant, hit_pope_large_scale_constant
    real(dp) :: hit_pope_dissipation_constant
    real(dp) :: hit_pope_large_scale_exponent
    real(dp) :: hit_pope_dissipation_exponent
    real(dp) :: hit_dealias_fraction, hit_isotropy_k_cutoff
    real(dp) :: hit_isotropy_tolerance
    real(dp) :: imported_turbulence_x_start
    real(dp) :: imported_turbulence_x_length
    real(dp) :: imported_turbulence_velocity_offset_x
    real(dp) :: imported_turbulence_velocity_offset_y
    real(dp) :: imported_turbulence_velocity_offset_z
    real(dp) :: imported_turbulence_background_rho
    real(dp) :: imported_turbulence_background_u
    real(dp) :: imported_turbulence_background_v
    real(dp) :: imported_turbulence_background_w
    real(dp) :: imported_turbulence_background_p
    real(dp) :: planar_shock_position, planar_shock_mach
    real(dp) :: planar_shock_upstream_rho, planar_shock_upstream_u
    real(dp) :: planar_shock_upstream_v, planar_shock_upstream_w
    real(dp) :: planar_shock_upstream_p
    real(dp) :: planar_shock_downstream_rho, planar_shock_downstream_u
    real(dp) :: planar_shock_downstream_v, planar_shock_downstream_w
    real(dp) :: planar_shock_downstream_p
    real(dp) :: shock_tube_diaphragm_position
    real(dp) :: shock_tube_driver_rho, shock_tube_driver_u
    real(dp) :: shock_tube_driver_v, shock_tube_driver_w
    real(dp) :: shock_tube_driver_p
    real(dp) :: shock_tube_driven_rho, shock_tube_driven_u
    real(dp) :: shock_tube_driven_v, shock_tube_driven_w
    real(dp) :: shock_tube_driven_p
    real(dp) :: forcing_k_cutoff, forcing_target_dissipation
    real(dp) :: forcing_dilatational_ratio, forcing_denominator_floor
    real(dp) :: forcing_max_coefficient
    real(dp) :: hybrid_sensor_onset, hybrid_sensor_full
    real(dp) :: boundary_x_min_reference_rho, boundary_x_min_reference_u
    real(dp) :: boundary_x_min_reference_v, boundary_x_min_reference_w
    real(dp) :: boundary_x_min_reference_p
    real(dp) :: boundary_x_max_reference_rho, boundary_x_max_reference_u
    real(dp) :: boundary_x_max_reference_v, boundary_x_max_reference_w
    real(dp) :: boundary_x_max_reference_p
    real(dp) :: boundary_y_min_reference_rho, boundary_y_min_reference_u
    real(dp) :: boundary_y_min_reference_v, boundary_y_min_reference_w
    real(dp) :: boundary_y_min_reference_p
    real(dp) :: boundary_y_max_reference_rho, boundary_y_max_reference_u
    real(dp) :: boundary_y_max_reference_v, boundary_y_max_reference_w
    real(dp) :: boundary_y_max_reference_p
    real(dp) :: boundary_z_min_reference_rho, boundary_z_min_reference_u
    real(dp) :: boundary_z_min_reference_v, boundary_z_min_reference_w
    real(dp) :: boundary_z_min_reference_p
    real(dp) :: boundary_z_max_reference_rho, boundary_z_max_reference_u
    real(dp) :: boundary_z_max_reference_v, boundary_z_max_reference_w
    real(dp) :: boundary_z_max_reference_p
    real(dp) :: boundary_relaxation_strength, boundary_length_scale
    integer :: hit_seed, hit_isotropy_max_iterations
    integer :: imported_turbulence_blend_cells
    integer :: forcing_report_interval
    character(len=32) :: convective_scheme, viscous_scheme, hit_spectrum
    character(len=32) :: hybrid_smooth_scheme, hybrid_shock_scheme
    character(len=32) :: hybrid_sensor
    character(len=32) :: hit_isotropy_mode
    character(len=512) :: imported_turbulence_file
    character(len=32) :: imported_turbulence_mode
    character(len=32) :: planar_shock_direction
    character(len=32) :: boundary_condition, time_integrator
    character(len=32) :: boundary_x_min, boundary_x_max
    character(len=32) :: boundary_y_min, boundary_y_max
    character(len=32) :: boundary_z_min, boundary_z_max
    character(len=32) :: forcing_scheme, forcing_spectrum, forcing_fft_backend
    integer :: u, ios
    logical :: exists, any_boundary_face, all_boundary_faces
    namelist /nse/ fh_enabled, fh_boltzmann_number, fh_seed, &
      nv, nghost, gamma, cfl, small_rho, small_p, rho0, mach, &
      reynolds, prandtl, convective_scheme, hybrid_smooth_scheme, &
      hybrid_shock_scheme, hybrid_sensor, hybrid_sensor_onset, &
      hybrid_sensor_full, viscous_scheme, &
      boundary_condition, boundary_x_min, boundary_x_max, &
      boundary_y_min, boundary_y_max, boundary_z_min, boundary_z_max, &
      boundary_x_min_reference_rho, boundary_x_min_reference_u, &
      boundary_x_min_reference_v, boundary_x_min_reference_w, &
      boundary_x_min_reference_p, boundary_x_max_reference_rho, &
      boundary_x_max_reference_u, boundary_x_max_reference_v, &
      boundary_x_max_reference_w, boundary_x_max_reference_p, &
      boundary_y_min_reference_rho, boundary_y_min_reference_u, &
      boundary_y_min_reference_v, boundary_y_min_reference_w, &
      boundary_y_min_reference_p, boundary_y_max_reference_rho, &
      boundary_y_max_reference_u, boundary_y_max_reference_v, &
      boundary_y_max_reference_w, boundary_y_max_reference_p, &
      boundary_z_min_reference_rho, boundary_z_min_reference_u, &
      boundary_z_min_reference_v, boundary_z_min_reference_w, &
      boundary_z_min_reference_p, boundary_z_max_reference_rho, &
      boundary_z_max_reference_u, boundary_z_max_reference_v, &
      boundary_z_max_reference_w, boundary_z_max_reference_p, &
      boundary_relaxation_strength, boundary_length_scale, &
      time_integrator, hit_spectrum, hit_seed, &
      hit_turbulent_mach, hit_turbulent_reynolds, hit_rms_velocity, &
      hit_peak_wavenumber, hit_integral_length, hit_kolmogorov_length, &
      hit_johnsen_length_scale_ratio, hit_pope_energy_constant, &
      hit_pope_large_scale_constant, hit_pope_dissipation_constant, &
      hit_pope_large_scale_exponent, hit_pope_dissipation_exponent, &
      hit_dealias_fraction, hit_isotropy_mode, &
      hit_isotropy_k_cutoff, hit_isotropy_tolerance, &
      hit_isotropy_max_iterations, imported_turbulence_file, &
      imported_turbulence_mode, imported_turbulence_x_start, &
      imported_turbulence_x_length, &
      imported_turbulence_blend_cells, &
      imported_turbulence_velocity_offset_x, &
      imported_turbulence_velocity_offset_y, &
      imported_turbulence_velocity_offset_z, &
      imported_turbulence_background_rho, &
      imported_turbulence_background_u, &
      imported_turbulence_background_v, &
      imported_turbulence_background_w, &
      imported_turbulence_background_p, planar_shock_position, &
      planar_shock_direction, planar_shock_mach, &
      planar_shock_upstream_rho, planar_shock_upstream_u, &
      planar_shock_upstream_v, planar_shock_upstream_w, &
      planar_shock_upstream_p, planar_shock_downstream_rho, &
      planar_shock_downstream_u, planar_shock_downstream_v, &
      planar_shock_downstream_w, planar_shock_downstream_p, &
      shock_tube_diaphragm_position, shock_tube_driver_rho, &
      shock_tube_driver_u, shock_tube_driver_v, shock_tube_driver_w, &
      shock_tube_driver_p, shock_tube_driven_rho, shock_tube_driven_u, &
      shock_tube_driven_v, shock_tube_driven_w, shock_tube_driven_p, &
      forcing_scheme, &
      forcing_spectrum, forcing_fft_backend, forcing_k_cutoff, &
      forcing_target_dissipation, forcing_dilatational_ratio, &
      forcing_denominator_floor, forcing_max_coefficient, &
      forcing_report_interval

    nv = cfg%nv
    nghost = 3
    gamma = cfg%gamma; cfl = cfg%cfl
    small_rho = cfg%small_rho; small_p = cfg%small_p
    rho0 = cfg%rho0; mach = cfg%mach; reynolds = cfg%reynolds; prandtl = cfg%prandtl
    fh_enabled = cfg%fh_enabled
    fh_boltzmann_number = cfg%fh_boltzmann_number
    fh_seed = cfg%fh_seed
    convective_scheme = cfg%convective_scheme
    hybrid_smooth_scheme = cfg%hybrid_smooth_scheme
    hybrid_shock_scheme = cfg%hybrid_shock_scheme
    hybrid_sensor = cfg%hybrid_sensor
    hybrid_sensor_onset = cfg%hybrid_sensor_onset
    hybrid_sensor_full = cfg%hybrid_sensor_full
    viscous_scheme = cfg%viscous_scheme
    boundary_condition = ''
    boundary_x_min = ''
    boundary_x_max = ''
    boundary_y_min = ''
    boundary_y_max = ''
    boundary_z_min = ''
    boundary_z_max = ''
    boundary_x_min_reference_rho = cfg%boundary_reference_rho(1)
    boundary_x_min_reference_u = cfg%boundary_reference_velocity(1,1)
    boundary_x_min_reference_v = cfg%boundary_reference_velocity(2,1)
    boundary_x_min_reference_w = cfg%boundary_reference_velocity(3,1)
    boundary_x_min_reference_p = cfg%boundary_reference_p(1)
    boundary_x_max_reference_rho = cfg%boundary_reference_rho(2)
    boundary_x_max_reference_u = cfg%boundary_reference_velocity(1,2)
    boundary_x_max_reference_v = cfg%boundary_reference_velocity(2,2)
    boundary_x_max_reference_w = cfg%boundary_reference_velocity(3,2)
    boundary_x_max_reference_p = cfg%boundary_reference_p(2)
    boundary_y_min_reference_rho = cfg%boundary_reference_rho(3)
    boundary_y_min_reference_u = cfg%boundary_reference_velocity(1,3)
    boundary_y_min_reference_v = cfg%boundary_reference_velocity(2,3)
    boundary_y_min_reference_w = cfg%boundary_reference_velocity(3,3)
    boundary_y_min_reference_p = cfg%boundary_reference_p(3)
    boundary_y_max_reference_rho = cfg%boundary_reference_rho(4)
    boundary_y_max_reference_u = cfg%boundary_reference_velocity(1,4)
    boundary_y_max_reference_v = cfg%boundary_reference_velocity(2,4)
    boundary_y_max_reference_w = cfg%boundary_reference_velocity(3,4)
    boundary_y_max_reference_p = cfg%boundary_reference_p(4)
    boundary_z_min_reference_rho = cfg%boundary_reference_rho(5)
    boundary_z_min_reference_u = cfg%boundary_reference_velocity(1,5)
    boundary_z_min_reference_v = cfg%boundary_reference_velocity(2,5)
    boundary_z_min_reference_w = cfg%boundary_reference_velocity(3,5)
    boundary_z_min_reference_p = cfg%boundary_reference_p(5)
    boundary_z_max_reference_rho = cfg%boundary_reference_rho(6)
    boundary_z_max_reference_u = cfg%boundary_reference_velocity(1,6)
    boundary_z_max_reference_v = cfg%boundary_reference_velocity(2,6)
    boundary_z_max_reference_w = cfg%boundary_reference_velocity(3,6)
    boundary_z_max_reference_p = cfg%boundary_reference_p(6)
    boundary_relaxation_strength = cfg%boundary_relaxation_strength
    boundary_length_scale = cfg%boundary_length_scale
    time_integrator = cfg%time_integrator
    hit_spectrum = cfg%hit_spectrum
    hit_seed = cfg%hit_seed
    hit_turbulent_mach = cfg%hit_turbulent_mach
    hit_turbulent_reynolds = cfg%hit_turbulent_reynolds
    hit_rms_velocity = cfg%hit_rms_velocity
    hit_peak_wavenumber = cfg%hit_peak_wavenumber
    hit_integral_length = cfg%hit_integral_length
    hit_kolmogorov_length = cfg%hit_kolmogorov_length
    hit_johnsen_length_scale_ratio = &
      cfg%hit_johnsen_length_scale_ratio
    hit_pope_energy_constant = cfg%hit_pope_energy_constant
    hit_pope_large_scale_constant = cfg%hit_pope_large_scale_constant
    hit_pope_dissipation_constant = &
      cfg%hit_pope_dissipation_constant
    hit_pope_large_scale_exponent = &
      cfg%hit_pope_large_scale_exponent
    hit_pope_dissipation_exponent = &
      cfg%hit_pope_dissipation_exponent
    hit_dealias_fraction = cfg%hit_dealias_fraction
    hit_isotropy_mode = cfg%hit_isotropy_mode
    hit_isotropy_k_cutoff = cfg%hit_isotropy_k_cutoff
    hit_isotropy_tolerance = cfg%hit_isotropy_tolerance
    hit_isotropy_max_iterations = cfg%hit_isotropy_max_iterations
    imported_turbulence_file = cfg%imported_turbulence_file
    imported_turbulence_mode = cfg%imported_turbulence_mode
    imported_turbulence_x_start = cfg%imported_turbulence_x_start
    imported_turbulence_x_length = cfg%imported_turbulence_x_length
    imported_turbulence_blend_cells = &
      cfg%imported_turbulence_blend_cells
    imported_turbulence_velocity_offset_x = &
      cfg%imported_turbulence_velocity_offset_x
    imported_turbulence_velocity_offset_y = &
      cfg%imported_turbulence_velocity_offset_y
    imported_turbulence_velocity_offset_z = &
      cfg%imported_turbulence_velocity_offset_z
    imported_turbulence_background_rho = &
      cfg%imported_turbulence_background_rho
    imported_turbulence_background_u = &
      cfg%imported_turbulence_background_u
    imported_turbulence_background_v = &
      cfg%imported_turbulence_background_v
    imported_turbulence_background_w = &
      cfg%imported_turbulence_background_w
    imported_turbulence_background_p = &
      cfg%imported_turbulence_background_p
    planar_shock_position = cfg%planar_shock_position
    planar_shock_direction = cfg%planar_shock_direction
    planar_shock_mach = cfg%planar_shock_mach
    planar_shock_upstream_rho = cfg%planar_shock_upstream_rho
    planar_shock_upstream_u = cfg%planar_shock_upstream_u
    planar_shock_upstream_v = cfg%planar_shock_upstream_v
    planar_shock_upstream_w = cfg%planar_shock_upstream_w
    planar_shock_upstream_p = cfg%planar_shock_upstream_p
    planar_shock_downstream_rho = cfg%planar_shock_downstream_rho
    planar_shock_downstream_u = cfg%planar_shock_downstream_u
    planar_shock_downstream_v = cfg%planar_shock_downstream_v
    planar_shock_downstream_w = cfg%planar_shock_downstream_w
    planar_shock_downstream_p = cfg%planar_shock_downstream_p
    shock_tube_diaphragm_position = cfg%shock_tube_diaphragm_position
    shock_tube_driver_rho = cfg%shock_tube_driver_rho
    shock_tube_driver_u = cfg%shock_tube_driver_u
    shock_tube_driver_v = cfg%shock_tube_driver_v
    shock_tube_driver_w = cfg%shock_tube_driver_w
    shock_tube_driver_p = cfg%shock_tube_driver_p
    shock_tube_driven_rho = cfg%shock_tube_driven_rho
    shock_tube_driven_u = cfg%shock_tube_driven_u
    shock_tube_driven_v = cfg%shock_tube_driven_v
    shock_tube_driven_w = cfg%shock_tube_driven_w
    shock_tube_driven_p = cfg%shock_tube_driven_p
    forcing_scheme = cfg%forcing_scheme
    forcing_spectrum = cfg%forcing_spectrum
    forcing_fft_backend = cfg%forcing_fft_backend
    forcing_k_cutoff = cfg%forcing_k_cutoff
    forcing_target_dissipation = cfg%forcing_target_dissipation
    forcing_dilatational_ratio = cfg%forcing_dilatational_ratio
    forcing_denominator_floor = cfg%forcing_denominator_floor
    forcing_max_coefficient = cfg%forcing_max_coefficient
    forcing_report_interval = cfg%forcing_report_interval

    inquire(file=filename, exist=exists)
    if (.not. exists) return
    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot open input file.'
    read(u, nml=nse, iostat=ios)
    close(u)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: failed to read namelist /nse/ in ', trim(filename)
      error stop
    end if

    cfg%nv = nv
    cfg%gamma = gamma; cfg%cfl = cfl
    cfg%small_rho = small_rho; cfg%small_p = small_p
    cfg%rho0 = rho0; cfg%mach = mach; cfg%reynolds = reynolds; cfg%prandtl = prandtl
    cfg%fh_enabled = fh_enabled
    cfg%fh_boltzmann_number = fh_boltzmann_number
    cfg%fh_seed = fh_seed
    cfg%convective_scheme = convective_scheme
    cfg%hybrid_smooth_scheme = hybrid_smooth_scheme
    cfg%hybrid_shock_scheme = hybrid_shock_scheme
    cfg%hybrid_sensor = hybrid_sensor
    cfg%hybrid_sensor_onset = hybrid_sensor_onset
    cfg%hybrid_sensor_full = hybrid_sensor_full
    cfg%viscous_scheme = viscous_scheme
    any_boundary_face = len_trim(boundary_x_min) > 0 .or. &
      len_trim(boundary_x_max) > 0 .or. len_trim(boundary_y_min) > 0 .or. &
      len_trim(boundary_y_max) > 0 .or. len_trim(boundary_z_min) > 0 .or. &
      len_trim(boundary_z_max) > 0
    all_boundary_faces = len_trim(boundary_x_min) > 0 .and. &
      len_trim(boundary_x_max) > 0 .and. len_trim(boundary_y_min) > 0 .and. &
      len_trim(boundary_y_max) > 0 .and. len_trim(boundary_z_min) > 0 .and. &
      len_trim(boundary_z_max) > 0
    if (any_boundary_face) then
      if (.not. all_boundary_faces) then
        error stop 'all six boundary face types must be specified together'
      end if
      if (len_trim(boundary_condition) > 0) then
        error stop 'boundary_condition cannot be combined with face boundaries'
      end if
      cfg%boundary_face_type(1) = boundary_x_min
      cfg%boundary_face_type(2) = boundary_x_max
      cfg%boundary_face_type(3) = boundary_y_min
      cfg%boundary_face_type(4) = boundary_y_max
      cfg%boundary_face_type(5) = boundary_z_min
      cfg%boundary_face_type(6) = boundary_z_max
    else
      if (len_trim(boundary_condition) == 0) boundary_condition = 'periodic'
      if (trim(adjustl(boundary_condition)) /= 'periodic') then
        error stop 'legacy boundary_condition supports only periodic'
      end if
      cfg%boundary_face_type = 'periodic'
    end if
    if (all(cfg%boundary_face_type == 'periodic')) then
      cfg%boundary_condition = 'periodic'
    else
      cfg%boundary_condition = 'mixed'
    end if
    cfg%boundary_reference_rho = [boundary_x_min_reference_rho, &
      boundary_x_max_reference_rho, boundary_y_min_reference_rho, &
      boundary_y_max_reference_rho, boundary_z_min_reference_rho, &
      boundary_z_max_reference_rho]
    cfg%boundary_reference_velocity(:,1) = [boundary_x_min_reference_u, &
      boundary_x_min_reference_v, boundary_x_min_reference_w]
    cfg%boundary_reference_velocity(:,2) = [boundary_x_max_reference_u, &
      boundary_x_max_reference_v, boundary_x_max_reference_w]
    cfg%boundary_reference_velocity(:,3) = [boundary_y_min_reference_u, &
      boundary_y_min_reference_v, boundary_y_min_reference_w]
    cfg%boundary_reference_velocity(:,4) = [boundary_y_max_reference_u, &
      boundary_y_max_reference_v, boundary_y_max_reference_w]
    cfg%boundary_reference_velocity(:,5) = [boundary_z_min_reference_u, &
      boundary_z_min_reference_v, boundary_z_min_reference_w]
    cfg%boundary_reference_velocity(:,6) = [boundary_z_max_reference_u, &
      boundary_z_max_reference_v, boundary_z_max_reference_w]
    cfg%boundary_reference_p = [boundary_x_min_reference_p, &
      boundary_x_max_reference_p, boundary_y_min_reference_p, &
      boundary_y_max_reference_p, boundary_z_min_reference_p, &
      boundary_z_max_reference_p]
    cfg%boundary_relaxation_strength = boundary_relaxation_strength
    cfg%boundary_length_scale = boundary_length_scale
    cfg%time_integrator = time_integrator
    cfg%hit_spectrum = hit_spectrum
    cfg%hit_seed = hit_seed
    cfg%hit_turbulent_mach = hit_turbulent_mach
    cfg%hit_turbulent_reynolds = hit_turbulent_reynolds
    cfg%hit_rms_velocity = hit_rms_velocity
    cfg%hit_peak_wavenumber = hit_peak_wavenumber
    cfg%hit_integral_length = hit_integral_length
    cfg%hit_kolmogorov_length = hit_kolmogorov_length
    cfg%hit_johnsen_length_scale_ratio = &
      hit_johnsen_length_scale_ratio
    cfg%hit_pope_energy_constant = hit_pope_energy_constant
    cfg%hit_pope_large_scale_constant = hit_pope_large_scale_constant
    cfg%hit_pope_dissipation_constant = &
      hit_pope_dissipation_constant
    cfg%hit_pope_large_scale_exponent = &
      hit_pope_large_scale_exponent
    cfg%hit_pope_dissipation_exponent = &
      hit_pope_dissipation_exponent
    cfg%hit_dealias_fraction = hit_dealias_fraction
    cfg%hit_isotropy_mode = hit_isotropy_mode
    cfg%hit_isotropy_k_cutoff = hit_isotropy_k_cutoff
    cfg%hit_isotropy_tolerance = hit_isotropy_tolerance
    cfg%hit_isotropy_max_iterations = hit_isotropy_max_iterations
    cfg%imported_turbulence_file = imported_turbulence_file
    cfg%imported_turbulence_mode = imported_turbulence_mode
    cfg%imported_turbulence_x_start = imported_turbulence_x_start
    cfg%imported_turbulence_x_length = imported_turbulence_x_length
    cfg%imported_turbulence_blend_cells = &
      imported_turbulence_blend_cells
    cfg%imported_turbulence_velocity_offset_x = &
      imported_turbulence_velocity_offset_x
    cfg%imported_turbulence_velocity_offset_y = &
      imported_turbulence_velocity_offset_y
    cfg%imported_turbulence_velocity_offset_z = &
      imported_turbulence_velocity_offset_z
    cfg%imported_turbulence_background_rho = &
      imported_turbulence_background_rho
    cfg%imported_turbulence_background_u = &
      imported_turbulence_background_u
    cfg%imported_turbulence_background_v = &
      imported_turbulence_background_v
    cfg%imported_turbulence_background_w = &
      imported_turbulence_background_w
    cfg%imported_turbulence_background_p = &
      imported_turbulence_background_p
    cfg%planar_shock_position = planar_shock_position
    cfg%planar_shock_direction = planar_shock_direction
    cfg%planar_shock_mach = planar_shock_mach
    cfg%planar_shock_upstream_rho = planar_shock_upstream_rho
    cfg%planar_shock_upstream_u = planar_shock_upstream_u
    cfg%planar_shock_upstream_v = planar_shock_upstream_v
    cfg%planar_shock_upstream_w = planar_shock_upstream_w
    cfg%planar_shock_upstream_p = planar_shock_upstream_p
    cfg%planar_shock_downstream_rho = planar_shock_downstream_rho
    cfg%planar_shock_downstream_u = planar_shock_downstream_u
    cfg%planar_shock_downstream_v = planar_shock_downstream_v
    cfg%planar_shock_downstream_w = planar_shock_downstream_w
    cfg%planar_shock_downstream_p = planar_shock_downstream_p
    cfg%shock_tube_diaphragm_position = shock_tube_diaphragm_position
    cfg%shock_tube_driver_rho = shock_tube_driver_rho
    cfg%shock_tube_driver_u = shock_tube_driver_u
    cfg%shock_tube_driver_v = shock_tube_driver_v
    cfg%shock_tube_driver_w = shock_tube_driver_w
    cfg%shock_tube_driver_p = shock_tube_driver_p
    cfg%shock_tube_driven_rho = shock_tube_driven_rho
    cfg%shock_tube_driven_u = shock_tube_driven_u
    cfg%shock_tube_driven_v = shock_tube_driven_v
    cfg%shock_tube_driven_w = shock_tube_driven_w
    cfg%shock_tube_driven_p = shock_tube_driven_p
    cfg%forcing_scheme = forcing_scheme
    cfg%forcing_spectrum = forcing_spectrum
    cfg%forcing_fft_backend = forcing_fft_backend
    cfg%forcing_k_cutoff = forcing_k_cutoff
    cfg%forcing_target_dissipation = forcing_target_dissipation
    cfg%forcing_dilatational_ratio = forcing_dilatational_ratio
    cfg%forcing_denominator_floor = forcing_denominator_floor
    cfg%forcing_max_coefficient = forcing_max_coefficient
    cfg%forcing_report_interval = forcing_report_interval
  end subroutine read_nse_input

  subroutine read_all_inputs(filename, sim, gpe, nse)
    character(len=*), intent(in) :: filename
    type(simulation_config), intent(inout) :: sim
    type(gpe_config), intent(inout), optional :: gpe
    type(nse_config), intent(inout), optional :: nse

    call read_common_input(filename, sim)
    if (present(gpe)) call read_gpe_input(filename, gpe)
    if (present(nse)) then
      nse%cfl = sim%cfl
      call read_nse_input(filename, nse)
    end if
  end subroutine read_all_inputs

end module mod_input_reader
