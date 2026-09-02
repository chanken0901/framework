module mod_mc_thermodynamics_provider
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  implicit none
  private

  character(len=*), parameter, public :: mc_thermodynamics_provider_name = &
    'calorically_perfect'
  logical, parameter, public :: mc_thermodynamics_supports_reactions = .false.

  public :: validate_mc_thermodynamics_provider
  public :: configure_mc_thermodynamics
  public :: mc_mixture_density
  public :: mc_mixture_gas_constant
  public :: mc_mixture_gamma
  public :: mc_pressure
  public :: mc_temperature
  public :: mc_sound_speed
  public :: mc_total_energy_from_primitive

contains

  subroutine validate_mc_thermodynamics_provider(requested_model)
    character(len=*), intent(in) :: requested_model

    if (trim(adjustl(requested_model)) /= mc_thermodynamics_provider_name) then
      write(*,'(A,A,A,A)') 'ERROR: requested thermodynamics provider "', &
        trim(requested_model), '" but this executable contains "', &
        mc_thermodynamics_provider_name
      error stop 'multicomponent thermodynamics provider mismatch'
    end if
  end subroutine validate_mc_thermodynamics_provider

  subroutine configure_mc_thermodynamics(path, nspecies, species_names)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: species_names(:)

    if (len_trim(path) == 0 .or. nspecies < 1 .or. &
        size(species_names) < nspecies) then
      error stop 'invalid calorically-perfect species contract'
    end if
  end subroutine configure_mc_thermodynamics

  pure real(dp) function mc_mixture_density(state, layout) result(density)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout

    density = sum(state(layout%first_species:layout%last_species))
  end function mc_mixture_density

  pure real(dp) function mc_mixture_gas_constant( &
      mass_fractions, layout) result(gas_constant)
    real(dp), intent(in) :: mass_fractions(:)
    type(mc_state_layout), intent(in) :: layout

    gas_constant = 1.0_dp + &
      0.0_dp*sum(mass_fractions(1:layout%nspecies))
  end function mc_mixture_gas_constant

  pure real(dp) function mc_pressure(state, layout, gamma) result(pressure)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    real(dp) :: density, momentum_squared

    density = mc_mixture_density(state, layout)
    momentum_squared = sum(state(layout%momentum)**2)
    pressure = (gamma-1.0_dp) * &
      (state(layout%total_energy)-0.5_dp*momentum_squared/density)
  end function mc_pressure

  pure real(dp) function mc_temperature(state, layout, gamma) result(temperature)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma

    temperature = mc_pressure(state, layout, gamma) / &
      mc_mixture_density(state, layout)
  end function mc_temperature

  pure real(dp) function mc_mixture_gamma( &
      mass_fractions, layout, temperature, gamma) result(gamma_value)
    real(dp), intent(in) :: mass_fractions(:), temperature, gamma
    type(mc_state_layout), intent(in) :: layout

    gamma_value = gamma + 0.0_dp*temperature + &
      0.0_dp*sum(mass_fractions(1:layout%nspecies))
  end function mc_mixture_gamma

  pure real(dp) function mc_sound_speed(state, layout, gamma) result(speed)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma

    speed = sqrt(gamma*mc_pressure(state, layout, gamma) / &
      mc_mixture_density(state, layout))
  end function mc_sound_speed

  pure real(dp) function mc_total_energy_from_primitive( &
      layout, gamma, density, velocity, pressure, mass_fractions) &
      result(total_energy_density)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma, density, velocity(3), pressure
    real(dp), intent(in) :: mass_fractions(:)

    total_energy_density = pressure/(gamma-1.0_dp) + &
      0.5_dp*density*sum(velocity**2) + &
      0.0_dp*sum(mass_fractions(1:layout%nspecies))
  end function mc_total_energy_from_primitive

end module mod_mc_thermodynamics_provider
