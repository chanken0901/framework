module mod_mc_euler_field
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, mc_pressure, &
    mc_temperature, mc_total_energy_from_primitive
  implicit none
  private

  public :: initialize_mc_euler_state
  public :: set_mc_euler_conservative_state
  public :: validate_mc_euler_state
  public :: compute_mc_euler_totals
  public :: compute_mc_euler_minima

contains

  subroutine initialize_mc_euler_state(q, layout, config)
    real(dp), intent(out) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k
    real(dp) :: x, dx

    if (size(q,1) /= config%nx .or. size(q,2) /= config%ny .or. &
        size(q,3) /= config%nz .or. size(q,4) /= layout%nvariables) then
      error stop 'multicomponent Euler state allocation is inconsistent'
    end if
    dx = (config%x_max-config%x_min) / real(config%nx,dp)
    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          x = config%x_min + (real(i,dp)-0.5_dp)*dx
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
