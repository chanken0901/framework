module mod_mc_passive_scalar_config
  use mod_precision, only : dp
  implicit none
  private

  integer, parameter :: mc_path_length = 256

  type, public :: mc_passive_scalar_config
    integer :: nx = 32
    integer :: ny = 8
    integer :: nz = 8
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 1.0_dp
    real(dp) :: z_min = 0.0_dp
    real(dp) :: z_max = 1.0_dp
    real(dp) :: velocity(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: cfl = 0.45_dp
    real(dp) :: dt = 0.0_dp
    integer :: nsteps = 20
    character(len=32) :: initial_condition = 'gaussian'
    real(dp) :: tracer_background = 0.05_dp
    real(dp) :: tracer_amplitude = 0.90_dp
    real(dp) :: tracer_center(3) = [0.25_dp, 0.5_dp, 0.5_dp]
    real(dp) :: tracer_width = 0.08_dp
    character(len=32) :: advection_scheme = 'upwind1'
    character(len=32) :: boundary_condition = 'periodic'
    character(len=32) :: time_integrator = 'ssprk3'
    logical :: write_final = .true.
    character(len=mc_path_length) :: output_file = &
      'passive_scalar_final.csv'
  end type mc_passive_scalar_config

  public :: initialize_mc_passive_scalar_config
  public :: read_mc_passive_scalar_config
  public :: validate_mc_passive_scalar_config
  public :: mc_passive_scalar_timestep

contains

  subroutine initialize_mc_passive_scalar_config(config)
    type(mc_passive_scalar_config), intent(out) :: config

    config = mc_passive_scalar_config()
  end subroutine initialize_mc_passive_scalar_config

  subroutine read_mc_passive_scalar_config(path, config)
    character(len=*), intent(in) :: path
    type(mc_passive_scalar_config), intent(out) :: config
    integer :: unit, ios
    integer :: nx, ny, nz, nsteps
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp) :: velocity(3), cfl, dt
    real(dp) :: tracer_background, tracer_amplitude
    real(dp) :: tracer_center(3), tracer_width
    character(len=32) :: initial_condition
    character(len=32) :: advection_scheme, boundary_condition
    character(len=32) :: time_integrator
    logical :: write_final
    character(len=mc_path_length) :: output_file
    character(len=512) :: message
    namelist /passive_scalar/ nx, ny, nz, x_min, x_max, y_min, y_max, &
      z_min, z_max, velocity, cfl, dt, nsteps, initial_condition, &
      tracer_background, tracer_amplitude, tracer_center, tracer_width, &
      advection_scheme, boundary_condition, time_integrator, write_final, &
      output_file

    call initialize_mc_passive_scalar_config(config)
    nx = config%nx
    ny = config%ny
    nz = config%nz
    x_min = config%x_min
    x_max = config%x_max
    y_min = config%y_min
    y_max = config%y_max
    z_min = config%z_min
    z_max = config%z_max
    velocity = config%velocity
    cfl = config%cfl
    dt = config%dt
    nsteps = config%nsteps
    initial_condition = config%initial_condition
    tracer_background = config%tracer_background
    tracer_amplitude = config%tracer_amplitude
    tracer_center = config%tracer_center
    tracer_width = config%tracer_width
    advection_scheme = config%advection_scheme
    boundary_condition = config%boundary_condition
    time_integrator = config%time_integrator
    write_final = config%write_final
    output_file = config%output_file

    open(newunit=unit, file=trim(path), status='old', action='read', &
      iostat=ios, iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open passive-scalar input: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open passive-scalar input'
    end if
    read(unit, nml=passive_scalar, iostat=ios, iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid passive_scalar namelist: ', &
        trim(message)
      error stop 'failed to read passive-scalar input'
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
    config%velocity = velocity
    config%cfl = cfl
    config%dt = dt
    config%nsteps = nsteps
    config%initial_condition = adjustl(initial_condition)
    config%tracer_background = tracer_background
    config%tracer_amplitude = tracer_amplitude
    config%tracer_center = tracer_center
    config%tracer_width = tracer_width
    config%advection_scheme = adjustl(advection_scheme)
    config%boundary_condition = adjustl(boundary_condition)
    config%time_integrator = adjustl(time_integrator)
    config%write_final = write_final
    config%output_file = adjustl(output_file)
    call validate_mc_passive_scalar_config(config)
  end subroutine read_mc_passive_scalar_config

  subroutine validate_mc_passive_scalar_config(config)
    type(mc_passive_scalar_config), intent(in) :: config
    real(dp) :: rate, courant

    if (min(config%nx, config%ny, config%nz) < 2) then
      error stop 'passive-scalar grid dimensions must be at least two'
    end if
    if (config%x_max <= config%x_min .or. &
        config%y_max <= config%y_min .or. &
        config%z_max <= config%z_min) then
      error stop 'passive-scalar domain extents must be positive'
    end if
    if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      error stop 'passive-scalar CFL must be in (0, 1]'
    end if
    if (config%dt < 0.0_dp) then
      error stop 'passive-scalar dt must be non-negative'
    end if
    if (config%nsteps < 0) then
      error stop 'passive-scalar nsteps must be non-negative'
    end if
    if (trim(config%initial_condition) /= 'gaussian') then
      error stop 'stage-1 passive scalar supports initial_condition=gaussian'
    end if
    if (config%tracer_background < 0.0_dp .or. &
        config%tracer_amplitude < 0.0_dp .or. &
        config%tracer_background + config%tracer_amplitude > 1.0_dp) then
      error stop 'passive-scalar tracer values must remain in [0,1]'
    end if
    if (config%tracer_width <= 0.0_dp) then
      error stop 'passive-scalar tracer width must be positive'
    end if
    if (config%tracer_center(1) < config%x_min .or. &
        config%tracer_center(1) > config%x_max .or. &
        config%tracer_center(2) < config%y_min .or. &
        config%tracer_center(2) > config%y_max .or. &
        config%tracer_center(3) < config%z_min .or. &
        config%tracer_center(3) > config%z_max) then
      error stop 'passive-scalar tracer center must lie inside the domain'
    end if
    if (trim(config%advection_scheme) /= 'upwind1') then
      error stop 'stage-1 passive scalar supports advection_scheme=upwind1'
    end if
    if (trim(config%boundary_condition) /= 'periodic') then
      error stop 'stage-1 passive scalar supports periodic boundaries only'
    end if
    if (trim(config%time_integrator) /= 'ssprk3') then
      error stop 'stage-1 passive scalar supports time_integrator=ssprk3'
    end if
    if (config%write_final .and. len_trim(config%output_file) == 0) then
      error stop 'passive-scalar output file must not be empty'
    end if

    rate = advection_rate(config)
    if (config%dt <= 0.0_dp .and. rate <= 0.0_dp) then
      error stop 'zero velocity requires an explicit positive dt'
    end if
    if (config%dt > 0.0_dp) then
      courant = config%dt * rate
      if (courant > 1.0_dp + 100.0_dp*epsilon(1.0_dp)) then
        error stop 'passive-scalar fixed dt violates the upwind CFL limit'
      end if
    end if
  end subroutine validate_mc_passive_scalar_config

  real(dp) function mc_passive_scalar_timestep(config) result(dt)
    type(mc_passive_scalar_config), intent(in) :: config
    real(dp) :: rate

    if (config%dt > 0.0_dp) then
      dt = config%dt
      return
    end if
    rate = advection_rate(config)
    if (rate <= 0.0_dp) then
      error stop 'cannot derive passive-scalar dt from zero velocity'
    end if
    dt = config%cfl / rate
  end function mc_passive_scalar_timestep

  pure real(dp) function advection_rate(config) result(rate)
    type(mc_passive_scalar_config), intent(in) :: config
    real(dp) :: dx, dy, dz

    dx = (config%x_max-config%x_min) / real(config%nx, dp)
    dy = (config%y_max-config%y_min) / real(config%ny, dp)
    dz = (config%z_max-config%z_min) / real(config%nz, dp)
    rate = abs(config%velocity(1))/dx + abs(config%velocity(2))/dy + &
      abs(config%velocity(3))/dz
  end function advection_rate

end module mod_mc_passive_scalar_config
