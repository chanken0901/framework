module mod_mc_euler_config
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_max_species
  implicit none
  private

  integer, parameter :: mc_path_length = 256
  integer, parameter, public :: mc_boundary_face_count = 6
  integer, parameter, public :: mc_face_x_min = 1
  integer, parameter, public :: mc_face_x_max = 2
  integer, parameter, public :: mc_face_y_min = 3
  integer, parameter, public :: mc_face_y_max = 4
  integer, parameter, public :: mc_face_z_min = 5
  integer, parameter, public :: mc_face_z_max = 6

  type, public :: mc_euler_config
    character(len=32) :: geometry = 'cartesian'
    real(dp) :: nozzle_inlet_half_height = 1.0_dp
    real(dp) :: nozzle_throat_half_height = 0.5_dp
    real(dp) :: nozzle_exit_half_height = 1.0_dp
    real(dp) :: nozzle_throat_x = 0.5_dp
    integer :: nx = 64
    integer :: ny = 4
    integer :: nz = 4
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 1.0_dp
    real(dp) :: z_min = 0.0_dp
    real(dp) :: z_max = 1.0_dp
    real(dp) :: gamma = 1.4_dp
    real(dp) :: cfl = 0.35_dp
    real(dp) :: diffusion_cfl = 0.40_dp
    real(dp) :: dt = 0.0_dp
    integer :: nsteps = 20
    character(len=32) :: initial_condition = 'multispecies_sod_x'
    real(dp) :: interface_location = 0.5_dp
    real(dp) :: left_density = 1.0_dp
    real(dp) :: left_velocity(3) = 0.0_dp
    real(dp) :: left_pressure = 1.0_dp
    real(dp) :: left_mass_fractions(mc_max_species) = 0.0_dp
    real(dp) :: right_density = 0.125_dp
    real(dp) :: right_velocity(3) = 0.0_dp
    real(dp) :: right_pressure = 0.1_dp
    real(dp) :: right_mass_fractions(mc_max_species) = 0.0_dp
    real(dp) :: wave_density = 1.0_dp
    real(dp) :: wave_temperature = 300.0_dp
    real(dp) :: wave_velocity(3) = 0.0_dp
    real(dp) :: wave_mean_mass_fractions(mc_max_species) = 0.0_dp
    integer :: wave_positive_species = 1
    integer :: wave_negative_species = 2
    real(dp) :: wave_amplitude = 0.10_dp
    integer :: wave_wavenumber = 1
    character(len=32) :: riemann_solver = 'rusanov1'
    character(len=32) :: boundary_condition = 'periodic'
    character(len=32) :: boundary_face_types(mc_boundary_face_count) = &
      'periodic'
    real(dp) :: boundary_reference_densities(mc_boundary_face_count) = 1.0_dp
    real(dp) :: boundary_reference_velocities(3,mc_boundary_face_count) = &
      0.0_dp
    real(dp) :: boundary_reference_pressures(mc_boundary_face_count) = 1.0_dp
    real(dp) :: boundary_reference_mass_fractions( &
      mc_boundary_face_count,mc_max_species) = 0.0_dp
    real(dp) :: boundary_relaxation_strength = 0.1_dp
    real(dp) :: boundary_length_scale = -1.0_dp
    ! Internal domain metadata; zero means use the local (serial) extent.
    real(dp) :: global_boundary_lengths(3) = 0.0_dp
    character(len=32) :: time_integrator = 'ssprk3'
    logical :: write_final = .true.
    character(len=mc_path_length) :: output_file = &
      'multicomponent_euler_final.csv'
  end type mc_euler_config

  public :: initialize_mc_euler_config
  public :: read_mc_euler_config
  public :: validate_mc_euler_config

contains

  subroutine initialize_mc_euler_config(config, nspecies)
    type(mc_euler_config), intent(out) :: config
    integer, intent(in) :: nspecies

    if (nspecies < 1 .or. nspecies > mc_max_species) then
      error stop 'invalid species count for multicomponent Euler configuration'
    end if
    config = mc_euler_config()
    if (nspecies == 1) then
      config%left_mass_fractions(1) = 1.0_dp
      config%right_mass_fractions(1) = 1.0_dp
      config%wave_mean_mass_fractions(1) = 1.0_dp
    else
      config%left_mass_fractions(1) = 0.8_dp
      config%left_mass_fractions(nspecies) = 0.2_dp
      config%right_mass_fractions(1) = 0.2_dp
      config%right_mass_fractions(nspecies) = 0.8_dp
      config%wave_mean_mass_fractions(1) = 0.5_dp
      config%wave_mean_mass_fractions(nspecies) = 0.5_dp
      config%wave_negative_species = nspecies
    end if
    config%boundary_reference_mass_fractions(:,1) = 1.0_dp
  end subroutine initialize_mc_euler_config

  subroutine read_mc_euler_config(path, nspecies, config)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    type(mc_euler_config), intent(out) :: config
    integer :: unit, ios
    character(len=32) :: geometry
    real(dp) :: nozzle_inlet_half_height, nozzle_throat_half_height
    real(dp) :: nozzle_exit_half_height, nozzle_throat_x
    integer :: nx, ny, nz, nsteps
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp) :: gamma, cfl, diffusion_cfl, dt, interface_location
    real(dp) :: left_density, left_velocity(3), left_pressure
    real(dp) :: left_mass_fractions(mc_max_species)
    real(dp) :: right_density, right_velocity(3), right_pressure
    real(dp) :: right_mass_fractions(mc_max_species)
    real(dp) :: wave_density, wave_temperature, wave_velocity(3)
    real(dp) :: wave_mean_mass_fractions(mc_max_species), wave_amplitude
    integer :: wave_positive_species, wave_negative_species, wave_wavenumber
    character(len=32) :: initial_condition, riemann_solver
    character(len=32) :: boundary_condition, time_integrator
    character(len=32) :: boundary_face_types(mc_boundary_face_count)
    real(dp) :: boundary_reference_densities(mc_boundary_face_count)
    real(dp) :: boundary_reference_velocities(3,mc_boundary_face_count)
    real(dp) :: boundary_reference_pressures(mc_boundary_face_count)
    real(dp) :: boundary_reference_mass_fractions( &
      mc_boundary_face_count,mc_max_species)
    real(dp) :: boundary_relaxation_strength, boundary_length_scale
    logical :: write_final
    character(len=mc_path_length) :: output_file
    character(len=512) :: message
    namelist /multicomponent_euler/ geometry, nozzle_inlet_half_height, &
      nozzle_throat_half_height, nozzle_exit_half_height, nozzle_throat_x, nx, ny, nz, x_min, x_max, y_min, &
      y_max, z_min, z_max, gamma, cfl, diffusion_cfl, dt, nsteps, &
      initial_condition, &
      interface_location, left_density, left_velocity, left_pressure, &
      left_mass_fractions, right_density, right_velocity, right_pressure, &
      right_mass_fractions, wave_density, wave_temperature, wave_velocity, &
      wave_mean_mass_fractions, wave_positive_species, &
      wave_negative_species, wave_amplitude, wave_wavenumber, &
      riemann_solver, boundary_condition, boundary_face_types, &
      boundary_reference_densities, boundary_reference_velocities, &
      boundary_reference_pressures, boundary_reference_mass_fractions, &
      boundary_relaxation_strength, boundary_length_scale, &
      time_integrator, write_final, output_file

    call initialize_mc_euler_config(config, nspecies)
    geometry = config%geometry
    nozzle_inlet_half_height = config%nozzle_inlet_half_height
    nozzle_throat_half_height = config%nozzle_throat_half_height
    nozzle_exit_half_height = config%nozzle_exit_half_height
    nozzle_throat_x = config%nozzle_throat_x
    nx = config%nx
    ny = config%ny
    nz = config%nz
    x_min = config%x_min
    x_max = config%x_max
    y_min = config%y_min
    y_max = config%y_max
    z_min = config%z_min
    z_max = config%z_max
    gamma = config%gamma
    cfl = config%cfl
    diffusion_cfl = config%diffusion_cfl
    dt = config%dt
    nsteps = config%nsteps
    initial_condition = config%initial_condition
    interface_location = config%interface_location
    left_density = config%left_density
    left_velocity = config%left_velocity
    left_pressure = config%left_pressure
    left_mass_fractions = config%left_mass_fractions
    right_density = config%right_density
    right_velocity = config%right_velocity
    right_pressure = config%right_pressure
    right_mass_fractions = config%right_mass_fractions
    wave_density = config%wave_density
    wave_temperature = config%wave_temperature
    wave_velocity = config%wave_velocity
    wave_mean_mass_fractions = config%wave_mean_mass_fractions
    wave_positive_species = config%wave_positive_species
    wave_negative_species = config%wave_negative_species
    wave_amplitude = config%wave_amplitude
    wave_wavenumber = config%wave_wavenumber
    riemann_solver = config%riemann_solver
    boundary_condition = config%boundary_condition
    boundary_face_types = config%boundary_face_types
    boundary_reference_densities = config%boundary_reference_densities
    boundary_reference_velocities = config%boundary_reference_velocities
    boundary_reference_pressures = config%boundary_reference_pressures
    boundary_reference_mass_fractions = &
      config%boundary_reference_mass_fractions
    boundary_relaxation_strength = config%boundary_relaxation_strength
    boundary_length_scale = config%boundary_length_scale
    time_integrator = config%time_integrator
    write_final = config%write_final
    output_file = config%output_file

    open(newunit=unit, file=trim(path), status='old', action='read', &
      iostat=ios, iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open multicomponent Euler input: ', &
        trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open multicomponent Euler input'
    end if
    read(unit, nml=multicomponent_euler, iostat=ios, iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid multicomponent_euler namelist: ', &
        trim(message)
      error stop 'failed to read multicomponent Euler input'
    end if

    config%geometry = geometry
    config%nozzle_inlet_half_height = nozzle_inlet_half_height
    config%nozzle_throat_half_height = nozzle_throat_half_height
    config%nozzle_exit_half_height = nozzle_exit_half_height
    config%nozzle_throat_x = nozzle_throat_x
    if (geometry /= 'cartesian' .and. geometry /= 'planar_nozzle') &
      error stop 'unsupported geometry: use cartesian or planar_nozzle'
    if (geometry == 'planar_nozzle') then
      if (.not. all(ieee_is_finite([nozzle_inlet_half_height, nozzle_throat_half_height, &
          nozzle_exit_half_height, nozzle_throat_x]))) error stop 'non-finite nozzle geometry'
      if (min(nozzle_inlet_half_height,nozzle_throat_half_height,nozzle_exit_half_height) <= 0) &
        error stop 'nozzle half heights must be positive'
      if (nozzle_throat_x <= x_min .or. nozzle_throat_x >= x_max) &
        error stop 'nozzle throat must be inside the x domain'
      if (abs(y_min+1.0_dp) > 1e-12_dp .or. abs(y_max-1.0_dp) > 1e-12_dp) &
        error stop 'planar nozzle computational y domain must be [-1,1]'
    end if
    config%nx = nx
    config%ny = ny
    config%nz = nz
    config%x_min = x_min
    config%x_max = x_max
    config%y_min = y_min
    config%y_max = y_max
    config%z_min = z_min
    config%z_max = z_max
    config%gamma = gamma
    config%cfl = cfl
    config%diffusion_cfl = diffusion_cfl
    config%dt = dt
    config%nsteps = nsteps
    config%initial_condition = adjustl(initial_condition)
    config%interface_location = interface_location
    config%left_density = left_density
    config%left_velocity = left_velocity
    config%left_pressure = left_pressure
    config%left_mass_fractions = left_mass_fractions
    config%right_density = right_density
    config%right_velocity = right_velocity
    config%right_pressure = right_pressure
    config%right_mass_fractions = right_mass_fractions
    config%wave_density = wave_density
    config%wave_temperature = wave_temperature
    config%wave_velocity = wave_velocity
    config%wave_mean_mass_fractions = wave_mean_mass_fractions
    config%wave_positive_species = wave_positive_species
    config%wave_negative_species = wave_negative_species
    config%wave_amplitude = wave_amplitude
    config%wave_wavenumber = wave_wavenumber
    config%riemann_solver = adjustl(riemann_solver)
    config%boundary_condition = adjustl(boundary_condition)
    config%boundary_face_types = adjustl(boundary_face_types)
    config%boundary_reference_densities = boundary_reference_densities
    config%boundary_reference_velocities = boundary_reference_velocities
    config%boundary_reference_pressures = boundary_reference_pressures
    config%boundary_reference_mass_fractions = &
      boundary_reference_mass_fractions
    config%boundary_relaxation_strength = boundary_relaxation_strength
    config%boundary_length_scale = boundary_length_scale
    config%time_integrator = adjustl(time_integrator)
    config%write_final = write_final
    config%output_file = adjustl(output_file)
    call validate_mc_euler_config(config, nspecies)
  end subroutine read_mc_euler_config

  subroutine validate_mc_euler_config(config, nspecies)
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: nspecies
    real(dp), parameter :: fraction_tolerance = 1.0e-12_dp
    integer :: face
    logical :: all_periodic

    if (nspecies < 1 .or. nspecies > mc_max_species) then
      error stop 'invalid multicomponent Euler species count'
    end if
    if (min(config%nx,config%ny,config%nz) < 2) then
      error stop 'multicomponent Euler grid dimensions must be at least two'
    end if
    if (config%x_max <= config%x_min .or. &
        config%y_max <= config%y_min .or. &
        config%z_max <= config%z_min) then
      error stop 'multicomponent Euler domain extents must be positive'
    end if
    if (config%gamma <= 1.0_dp) then
      error stop 'calorically perfect gamma must exceed one'
    end if
    if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      error stop 'multicomponent Euler CFL must be in (0,1]'
    end if
    if (config%diffusion_cfl <= 0.0_dp .or. &
        config%diffusion_cfl > 1.0_dp) then
      error stop 'multicomponent diffusion CFL must be in (0,1]'
    end if
    if (config%dt < 0.0_dp .or. config%nsteps < 0) then
      error stop 'multicomponent Euler time settings are invalid'
    end if
    if (trim(config%initial_condition) /= 'multispecies_sod_x' .and. &
        trim(config%initial_condition) /= 'reactive_shock_tube_x' .and. &
        trim(config%initial_condition) /= 'periodic_species_wave_x') then
      error stop 'unsupported multicomponent initial condition'
    end if
    if (trim(config%initial_condition) == 'multispecies_sod_x' .or. &
        trim(config%initial_condition) == 'reactive_shock_tube_x') then
      if (config%interface_location <= config%x_min .or. &
          config%interface_location >= config%x_max) then
        error stop 'Euler interface location must lie inside the x domain'
      end if
      if (min(config%left_density,config%right_density) <= 0.0_dp .or. &
          min(config%left_pressure,config%right_pressure) <= 0.0_dp) then
        error stop 'Euler initial density and pressure must be positive'
      end if
      if (minval(config%left_mass_fractions(1:nspecies)) < 0.0_dp .or. &
          minval(config%right_mass_fractions(1:nspecies)) < 0.0_dp) then
        error stop 'Euler initial mass fractions must be non-negative'
      end if
      if (abs(sum(config%left_mass_fractions(1:nspecies))-1.0_dp) > &
          fraction_tolerance .or. &
          abs(sum(config%right_mass_fractions(1:nspecies))-1.0_dp) > &
          fraction_tolerance) then
        error stop 'Euler initial mass fractions must sum to one'
      end if
    else
      if (nspecies < 2) then
        error stop 'periodic species wave requires at least two species'
      end if
      if (config%wave_density <= 0.0_dp .or. &
          config%wave_temperature <= 0.0_dp) then
        error stop 'periodic species-wave density and temperature must be positive'
      end if
      if (minval(config%wave_mean_mass_fractions(1:nspecies)) < 0.0_dp .or. &
          abs(sum(config%wave_mean_mass_fractions(1:nspecies))-1.0_dp) > &
          fraction_tolerance) then
        error stop 'periodic species-wave mean fractions are invalid'
      end if
      if (min(config%wave_positive_species,config%wave_negative_species) < 1 .or. &
          max(config%wave_positive_species,config%wave_negative_species) > &
          nspecies .or. &
          config%wave_positive_species == config%wave_negative_species) then
        error stop 'periodic species-wave indices are invalid'
      end if
      if (config%wave_amplitude < 0.0_dp .or. &
          config%wave_amplitude >= min( &
          config%wave_mean_mass_fractions(config%wave_positive_species), &
          config%wave_mean_mass_fractions(config%wave_negative_species))) then
        error stop 'periodic species-wave amplitude violates positivity'
      end if
      if (config%wave_wavenumber < 1) then
        error stop 'periodic species-wave wavenumber must be positive'
      end if
    end if
    if (trim(config%riemann_solver) /= 'rusanov1') then
      error stop 'multicomponent Euler supports riemann_solver=rusanov1'
    end if
    if (.not. ieee_is_finite(config%boundary_relaxation_strength) .or. &
        config%boundary_relaxation_strength < 0.0_dp) then
      error stop 'multicomponent boundary relaxation must be non-negative'
    end if
    if (.not. ieee_is_finite(config%boundary_length_scale) .or. &
        (config%boundary_length_scale <= 0.0_dp .and. &
         abs(config%boundary_length_scale+1.0_dp) > &
         10.0_dp*epsilon(1.0_dp))) then
      error stop 'multicomponent boundary length scale must be AUTO or positive'
    end if
    do face = 1, mc_boundary_face_count
      select case (trim(config%boundary_face_types(face)))
      case ('periodic', 'reflective')
        continue
      case ('non_reflecting', 'dirichlet')
        if (.not. ieee_is_finite( &
            config%boundary_reference_densities(face)) .or. &
            config%boundary_reference_densities(face) <= 0.0_dp .or. &
            .not. ieee_is_finite( &
            config%boundary_reference_pressures(face)) .or. &
            config%boundary_reference_pressures(face) <= 0.0_dp .or. &
            .not. all(ieee_is_finite( &
            config%boundary_reference_velocities(:,face)))) then
          error stop 'invalid multicomponent boundary reference primitive state'
        end if
        if (minval(config%boundary_reference_mass_fractions( &
            face,1:nspecies)) < 0.0_dp .or. &
            abs(sum(config%boundary_reference_mass_fractions( &
            face,1:nspecies))-1.0_dp) > fraction_tolerance) then
          error stop 'invalid multicomponent boundary reference composition'
        end if
      case default
        error stop 'unsupported multicomponent boundary face type'
      end select
    end do
    if ((trim(config%boundary_face_types(mc_face_x_min)) == 'periodic') &
        .neqv. &
        (trim(config%boundary_face_types(mc_face_x_max)) == 'periodic')) then
      error stop 'multicomponent periodic x boundaries must be paired'
    end if
    if ((trim(config%boundary_face_types(mc_face_y_min)) == 'periodic') &
        .neqv. &
        (trim(config%boundary_face_types(mc_face_y_max)) == 'periodic')) then
      error stop 'multicomponent periodic y boundaries must be paired'
    end if
    if ((trim(config%boundary_face_types(mc_face_z_min)) == 'periodic') &
        .neqv. &
        (trim(config%boundary_face_types(mc_face_z_max)) == 'periodic')) then
      error stop 'multicomponent periodic z boundaries must be paired'
    end if
    all_periodic = all(config%boundary_face_types == 'periodic')
    if (all_periodic) then
      if (trim(config%boundary_condition) /= 'periodic') then
        error stop 'multicomponent boundary summary must be periodic'
      end if
    else if (trim(config%boundary_condition) /= 'face_specific') then
      error stop 'multicomponent face boundaries require face_specific summary'
    end if
    if (trim(config%time_integrator) /= 'ssprk3') then
      error stop 'multicomponent Euler supports time_integrator=ssprk3'
    end if
    if (config%write_final .and. len_trim(config%output_file) == 0) then
      error stop 'multicomponent Euler output file must not be empty'
    end if
    if(config%geometry == 'planar_nozzle') then
      if(any(config%boundary_face_types(3:4) == 'periodic')) &
        error stop 'planar nozzle y faces cannot be periodic'
      if(config%boundary_face_types(1) == 'periodic' .and. &
          abs(config%nozzle_inlet_half_height-config%nozzle_exit_half_height)>1e-12_dp) &
        error stop 'periodic nozzle x faces must have matching heights'
    end if
  end subroutine validate_mc_euler_config

end module mod_mc_euler_config
