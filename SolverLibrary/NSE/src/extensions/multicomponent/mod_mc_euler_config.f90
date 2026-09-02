module mod_mc_euler_config
  use mod_precision, only : dp
  use mod_mc_config, only : mc_max_species
  implicit none
  private

  integer, parameter :: mc_path_length = 256

  type, public :: mc_euler_config
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
    character(len=32) :: riemann_solver = 'rusanov1'
    character(len=32) :: boundary_condition = 'periodic'
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
    else
      config%left_mass_fractions(1) = 0.8_dp
      config%left_mass_fractions(nspecies) = 0.2_dp
      config%right_mass_fractions(1) = 0.2_dp
      config%right_mass_fractions(nspecies) = 0.8_dp
    end if
  end subroutine initialize_mc_euler_config

  subroutine read_mc_euler_config(path, nspecies, config)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    type(mc_euler_config), intent(out) :: config
    integer :: unit, ios
    integer :: nx, ny, nz, nsteps
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp) :: gamma, cfl, dt, interface_location
    real(dp) :: left_density, left_velocity(3), left_pressure
    real(dp) :: left_mass_fractions(mc_max_species)
    real(dp) :: right_density, right_velocity(3), right_pressure
    real(dp) :: right_mass_fractions(mc_max_species)
    character(len=32) :: initial_condition, riemann_solver
    character(len=32) :: boundary_condition, time_integrator
    logical :: write_final
    character(len=mc_path_length) :: output_file
    character(len=512) :: message
    namelist /multicomponent_euler/ nx, ny, nz, x_min, x_max, y_min, &
      y_max, z_min, z_max, gamma, cfl, dt, nsteps, initial_condition, &
      interface_location, left_density, left_velocity, left_pressure, &
      left_mass_fractions, right_density, right_velocity, right_pressure, &
      right_mass_fractions, riemann_solver, boundary_condition, &
      time_integrator, write_final, output_file

    call initialize_mc_euler_config(config, nspecies)
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
    riemann_solver = config%riemann_solver
    boundary_condition = config%boundary_condition
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
    config%riemann_solver = adjustl(riemann_solver)
    config%boundary_condition = adjustl(boundary_condition)
    config%time_integrator = adjustl(time_integrator)
    config%write_final = write_final
    config%output_file = adjustl(output_file)
    call validate_mc_euler_config(config, nspecies)
  end subroutine read_mc_euler_config

  subroutine validate_mc_euler_config(config, nspecies)
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: nspecies
    real(dp), parameter :: fraction_tolerance = 1.0e-12_dp

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
    if (config%dt < 0.0_dp .or. config%nsteps < 0) then
      error stop 'multicomponent Euler time settings are invalid'
    end if
    if (trim(config%initial_condition) /= 'multispecies_sod_x') then
      error stop 'multicomponent Euler supports initial_condition=multispecies_sod_x'
    end if
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
    if (trim(config%riemann_solver) /= 'rusanov1') then
      error stop 'multicomponent Euler supports riemann_solver=rusanov1'
    end if
    if (trim(config%boundary_condition) /= 'periodic') then
      error stop 'multicomponent Euler supports periodic boundaries only'
    end if
    if (trim(config%time_integrator) /= 'ssprk3') then
      error stop 'multicomponent Euler supports time_integrator=ssprk3'
    end if
    if (config%write_final .and. len_trim(config%output_file) == 0) then
      error stop 'multicomponent Euler output file must not be empty'
    end if
  end subroutine validate_mc_euler_config

end module mod_mc_euler_config
