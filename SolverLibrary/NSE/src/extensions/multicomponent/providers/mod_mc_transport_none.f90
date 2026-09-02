module mod_mc_transport_provider
  implicit none
  private

  character(len=*), parameter, public :: mc_transport_provider_name = 'none'
  logical, parameter, public :: mc_transport_has_species_diffusion = .false.

  public :: validate_mc_transport_provider

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

end module mod_mc_transport_provider
