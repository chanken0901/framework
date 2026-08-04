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
    real(dp) :: forcing_k_cutoff, forcing_target_dissipation
    real(dp) :: forcing_dilatational_ratio, forcing_denominator_floor
    real(dp) :: forcing_max_coefficient
    integer :: hit_seed, hit_isotropy_max_iterations
    integer :: forcing_report_interval
    character(len=32) :: convective_scheme, viscous_scheme, hit_spectrum
    character(len=32) :: hit_isotropy_mode
    character(len=32) :: boundary_condition, time_integrator
    character(len=32) :: forcing_scheme, forcing_spectrum, forcing_fft_backend
    integer :: u, ios
    logical :: exists
    namelist /nse/ nv, nghost, gamma, cfl, small_rho, small_p, rho0, mach, &
      reynolds, prandtl, convective_scheme, viscous_scheme, &
      boundary_condition, time_integrator, hit_spectrum, hit_seed, &
      hit_turbulent_mach, hit_turbulent_reynolds, hit_rms_velocity, &
      hit_peak_wavenumber, hit_integral_length, hit_kolmogorov_length, &
      hit_johnsen_length_scale_ratio, hit_pope_energy_constant, &
      hit_pope_large_scale_constant, hit_pope_dissipation_constant, &
      hit_pope_large_scale_exponent, hit_pope_dissipation_exponent, &
      hit_dealias_fraction, hit_isotropy_mode, &
      hit_isotropy_k_cutoff, hit_isotropy_tolerance, &
      hit_isotropy_max_iterations, forcing_scheme, &
      forcing_spectrum, forcing_fft_backend, forcing_k_cutoff, &
      forcing_target_dissipation, forcing_dilatational_ratio, &
      forcing_denominator_floor, forcing_max_coefficient, &
      forcing_report_interval

    nv = cfg%nv
    nghost = 3
    gamma = cfg%gamma; cfl = cfg%cfl
    small_rho = cfg%small_rho; small_p = cfg%small_p
    rho0 = cfg%rho0; mach = cfg%mach; reynolds = cfg%reynolds; prandtl = cfg%prandtl
    convective_scheme = cfg%convective_scheme
    viscous_scheme = cfg%viscous_scheme
    boundary_condition = cfg%boundary_condition
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
    cfg%convective_scheme = convective_scheme
    cfg%viscous_scheme = viscous_scheme
    cfg%boundary_condition = boundary_condition
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
