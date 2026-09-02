module mod_mc_thermodynamics_provider
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  implicit none
  private

  character(len=*), parameter, public :: mc_thermodynamics_provider_name = &
    'calorically_perfect'
  logical, parameter, public :: mc_thermodynamics_supports_reactions = .false.

  public :: validate_mc_thermodynamics_provider
  public :: mc_mixture_density
  public :: mc_pressure
  public :: mc_sound_speed

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

  pure real(dp) function mc_mixture_density(state, layout) result(density)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout

    density = sum(state(layout%first_species:layout%last_species))
  end function mc_mixture_density

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

  pure real(dp) function mc_sound_speed(state, layout, gamma) result(speed)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma

    speed = sqrt(gamma*mc_pressure(state, layout, gamma) / &
      mc_mixture_density(state, layout))
  end function mc_sound_speed

end module mod_mc_thermodynamics_provider
