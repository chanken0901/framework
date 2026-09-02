module mod_mc_transport_provider
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  implicit none
  private

  character(len=*), parameter, public :: mc_transport_provider_name = 'none'
  logical, parameter, public :: mc_transport_has_species_diffusion = .false.

  public :: validate_mc_transport_provider
  public :: configure_mc_transport
  public :: mc_dynamic_viscosity
  public :: mc_thermal_conductivity
  public :: mc_species_diffusivities

contains

  subroutine validate_mc_transport_provider(requested_model)
    character(len=*), intent(in) :: requested_model

    if (trim(adjustl(requested_model)) /= mc_transport_provider_name) then
      write(*,'(A,A,A,A)') 'ERROR: requested transport provider "', &
        trim(requested_model), '" but this executable contains "', &
        mc_transport_provider_name
      error stop 'multicomponent transport provider mismatch'
    end if
  end subroutine validate_mc_transport_provider

  subroutine configure_mc_transport(path,nspecies,species_names)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: species_names(:)

    if (len_trim(path) == 0 .or. nspecies < 1 .or. &
        size(species_names) < nspecies) then
      error stop 'invalid no-transport species contract'
    end if
  end subroutine configure_mc_transport

  pure real(dp) function mc_dynamic_viscosity( &
      mass_fractions,layout,temperature) result(viscosity)
    real(dp), intent(in) :: mass_fractions(:), temperature
    type(mc_state_layout), intent(in) :: layout

    viscosity = 0.0_dp*(temperature + &
      sum(mass_fractions(1:layout%nspecies)))
  end function mc_dynamic_viscosity

  pure real(dp) function mc_thermal_conductivity( &
      mass_fractions,layout,temperature,mixture_cp) result(conductivity)
    real(dp), intent(in) :: mass_fractions(:), temperature, mixture_cp
    type(mc_state_layout), intent(in) :: layout

    conductivity = 0.0_dp*(temperature+mixture_cp + &
      sum(mass_fractions(1:layout%nspecies)))
  end function mc_thermal_conductivity

  pure subroutine mc_species_diffusivities( &
      mass_fractions,layout,temperature,diffusivities)
    real(dp), intent(in) :: mass_fractions(:), temperature
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(out) :: diffusivities(:)

    if (size(diffusivities) < layout%nspecies) then
      error stop 'species diffusivity output vector is too short'
    end if
    diffusivities(1:layout%nspecies) = 0.0_dp*(temperature + &
      sum(mass_fractions(1:layout%nspecies)))
  end subroutine mc_species_diffusivities

end module mod_mc_transport_provider
