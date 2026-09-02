module mod_mc_chemistry_provider
  implicit none
  private

  character(len=*), parameter, public :: mc_chemistry_provider_name = 'none'
  logical, parameter, public :: mc_chemistry_is_reactive = .false.

  public :: validate_mc_chemistry_provider

contains

  subroutine validate_mc_chemistry_provider(requested_model)
    character(len=*), intent(in) :: requested_model

    if (trim(adjustl(requested_model)) /= mc_chemistry_provider_name) then
      write(*,'(A,A,A,A)') 'ERROR: requested chemistry provider "', &
        trim(requested_model), '" but this executable contains "', &
        mc_chemistry_provider_name
      error stop 'multicomponent chemistry provider mismatch'
    end if
  end subroutine validate_mc_chemistry_provider

end module mod_mc_chemistry_provider
