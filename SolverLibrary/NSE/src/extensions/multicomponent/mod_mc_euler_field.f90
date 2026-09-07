module mod_mc_euler_field
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, mc_pressure, &
    mc_temperature, mc_total_energy_from_primitive, &
    mc_mixture_gas_constant, mc_mixture_gamma
  implicit none
  private

  public :: initialize_mc_euler_state
  public :: set_mc_euler_conservative_state
  public :: validate_mc_euler_state
  public :: compute_mc_euler_totals
  public :: compute_mc_euler_minima

  ! Caller-owned scratch: refreshed from q on every RHS evaluation.
  type, public :: mc_primitive_workspace
    real(dp), allocatable :: density(:,:,:), temperature(:,:,:)
    real(dp), allocatable :: pressure(:,:,:), sound_speed(:,:,:)
    real(dp), allocatable :: velocity(:,:,:,:), mass_fractions(:,:,:,:)
  end type mc_primitive_workspace
  public :: prepare_mc_primitives, evaluate_mc_primitive

contains

  subroutine evaluate_mc_primitive(state,layout,gamma,density,velocity, &
      temperature,mass_fractions,pressure,sound_speed)
    real(dp), intent(in) :: state(:), gamma
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(out) :: density, velocity(3), temperature
    real(dp), intent(out) :: mass_fractions(:), pressure, sound_speed
    real(dp) :: gas_constant, gamma_value

    if (.not. all(ieee_is_finite(state))) &
      error stop 'non-finite multicomponent state'
    if (minval(state(layout%first_species:layout%last_species)) < -1.0e-13_dp) &
      error stop 'negative species partial density in Euler state'
    density = mc_mixture_density(state,layout)
    if (density <= 1.0e-13_dp) &
      error stop 'non-positive mixture density in Euler state'
    velocity = state(layout%momentum)/density
    mass_fractions = state(layout%first_species:layout%last_species)/density
    temperature = mc_temperature(state,layout,gamma)
    gas_constant = mc_mixture_gas_constant(mass_fractions,layout)
    pressure = density*gas_constant*temperature
    if (.not. ieee_is_finite(pressure) .or. pressure <= 1.0e-13_dp) &
      error stop 'non-positive or non-finite pressure in Euler state'
    gamma_value = mc_mixture_gamma(mass_fractions,layout,temperature,gamma)
    sound_speed = sqrt(gamma_value*pressure/density)
    if (.not. ieee_is_finite(sound_speed)) &
      error stop 'non-finite multicomponent sound speed'
  end subroutine evaluate_mc_primitive

  subroutine prepare_mc_primitives(q,layout,config,work)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_primitive_workspace), intent(inout) :: work
    integer :: i,j,k,nx,ny,nz

    nx = config%nx
    ny = config%ny
    nz = config%nz
    if (any(shape(q) /= [nx,ny,nz,layout%nvariables])) &
      error stop 'multicomponent primitive state shape mismatch'
    if (allocated(work%density)) then
      if (any(shape(work%density) /= [nx,ny,nz]) .or. &
          size(work%mass_fractions,4) /= layout%nspecies) then
        deallocate(work%density,work%temperature,work%pressure, &
          work%sound_speed,work%velocity,work%mass_fractions)
      end if
    end if
    if (.not. allocated(work%density)) then
      allocate(work%density(nx,ny,nz),work%temperature(nx,ny,nz), &
        work%pressure(nx,ny,nz),work%sound_speed(nx,ny,nz), &
        work%velocity(nx,ny,nz,3), &
        work%mass_fractions(nx,ny,nz,layout%nspecies))
    end if
    !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call evaluate_mc_primitive(q(i,j,k,:),layout,config%gamma, &
            work%density(i,j,k),work%velocity(i,j,k,:), &
            work%temperature(i,j,k),work%mass_fractions(i,j,k,:), &
            work%pressure(i,j,k),work%sound_speed(i,j,k))
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine prepare_mc_primitives

  subroutine initialize_mc_euler_state(q, layout, config, reference_config)
    real(dp), intent(out) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k
    type(mc_euler_config), intent(in), optional :: reference_config
    real(dp) :: wave_origin, wave_length
    real(dp) :: x, dx, phase, pressure, gas_constant
    real(dp) :: mass_fractions(layout%nspecies)
    real(dp), parameter :: pi = acos(-1.0_dp)

    if (size(q,1) /= config%nx .or. size(q,2) /= config%ny .or. &
        size(q,3) /= config%nz .or. size(q,4) /= layout%nvariables) then
      error stop 'multicomponent Euler state allocation is inconsistent'
    end if
    wave_origin=config%x_min
    wave_length=config%x_max-config%x_min
    if (present(reference_config)) then
      wave_origin=reference_config%x_min
      wave_length=reference_config%x_max-reference_config%x_min
    end if
    dx = (config%x_max-config%x_min) / real(config%nx,dp)
    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          x = config%x_min + (real(i,dp)-0.5_dp)*dx
          select case (trim(config%initial_condition))
          case ('multispecies_sod_x', 'reactive_shock_tube_x')
            if (x < config%interface_location) then
              call set_mc_euler_conservative_state( &
                q(i,j,k,:), layout, config%gamma, config%left_density, &
                config%left_velocity, config%left_pressure, &
                config%left_mass_fractions(1:layout%nspecies))
            else
              call set_mc_euler_conservative_state( &
                q(i,j,k,:), layout, config%gamma, config%right_density, &
                config%right_velocity, config%right_pressure, &
                config%right_mass_fractions(1:layout%nspecies))
            end if
          case ('periodic_species_wave_x')
            phase = 2.0_dp*pi*real(config%wave_wavenumber,dp)* &
              (x-wave_origin)/wave_length
            mass_fractions = &
              config%wave_mean_mass_fractions(1:layout%nspecies)
            mass_fractions(config%wave_positive_species) = &
              mass_fractions(config%wave_positive_species) + &
              config%wave_amplitude*sin(phase)
            mass_fractions(config%wave_negative_species) = &
              mass_fractions(config%wave_negative_species) - &
              config%wave_amplitude*sin(phase)
            gas_constant = mc_mixture_gas_constant(mass_fractions,layout)
            pressure = config%wave_density*gas_constant* &
              config%wave_temperature
            call set_mc_euler_conservative_state( &
              q(i,j,k,:),layout,config%gamma,config%wave_density, &
              config%wave_velocity,pressure,mass_fractions)
          case default
            error stop 'unsupported multicomponent initial condition'
          end select
        end do
      end do
    end do
    call validate_mc_euler_state(q, layout, config)
  end subroutine initialize_mc_euler_state

  subroutine set_mc_euler_conservative_state( &
      state, layout, gamma, density, velocity, pressure, mass_fractions)
    real(dp), intent(out) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma, density, velocity(3), pressure
    real(dp), intent(in) :: mass_fractions(:)

    state = 0.0_dp
    state(layout%first_species:layout%last_species) = &
      density*mass_fractions(1:layout%nspecies)
    state(layout%momentum) = density*velocity
    state(layout%total_energy) = mc_total_energy_from_primitive( &
      layout,gamma,density,velocity,pressure,mass_fractions)
  end subroutine set_mc_euler_conservative_state

  subroutine validate_mc_euler_state(q, layout, config)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k
    real(dp) :: density, pressure
    real(dp), parameter :: small_value = 1.0e-13_dp

    do k = 1, size(q,3)
      do j = 1, size(q,2)
        do i = 1, size(q,1)
          if (minval(q(i,j,k,layout%first_species:layout%last_species)) &
              < -small_value) then
            error stop 'negative species partial density in Euler state'
          end if
          density = mc_mixture_density(q(i,j,k,:), layout)
          if (density <= small_value) then
            error stop 'non-positive mixture density in Euler state'
          end if
          pressure = mc_pressure(q(i,j,k,:), layout, config%gamma)
          if (pressure <= small_value) then
            error stop 'non-positive pressure in Euler state'
          end if
        end do
      end do
    end do
  end subroutine validate_mc_euler_state

  subroutine compute_mc_euler_totals(q, config, totals)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_euler_config), intent(in) :: config
    real(dp), intent(out) :: totals(:)
    integer :: variable
    real(dp) :: cell_volume

    if (size(totals) /= size(q,4)) then
      error stop 'Euler conserved-total diagnostic has the wrong size'
    end if
    cell_volume = (config%x_max-config%x_min) / real(config%nx,dp) * &
      (config%y_max-config%y_min) / real(config%ny,dp) * &
      (config%z_max-config%z_min) / real(config%nz,dp)
    do variable = 1, size(q,4)
      totals(variable) = sum(q(:,:,:,variable))*cell_volume
    end do
  end subroutine compute_mc_euler_totals

  subroutine compute_mc_euler_minima( &
      q, layout, config, species_density, mixture_density, pressure, &
      temperature)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    real(dp), intent(out) :: species_density, mixture_density, pressure
    real(dp), intent(out), optional :: temperature
    integer :: i, j, k

    species_density = minval( &
      q(:,:,:,layout%first_species:layout%last_species))
    mixture_density = huge(1.0_dp)
    pressure = huge(1.0_dp)
    if (present(temperature)) temperature = huge(1.0_dp)
    do k = 1, size(q,3)
      do j = 1, size(q,2)
        do i = 1, size(q,1)
          mixture_density = min(mixture_density, &
            mc_mixture_density(q(i,j,k,:),layout))
          pressure = min(pressure, &
            mc_pressure(q(i,j,k,:),layout,config%gamma))
          if (present(temperature)) then
            temperature = min(temperature, &
              mc_temperature(q(i,j,k,:),layout,config%gamma))
          end if
        end do
      end do
    end do
  end subroutine compute_mc_euler_minima

end module mod_mc_euler_field
