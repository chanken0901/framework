module mod_mc_state_layout
  implicit none
  private

  type, public :: mc_state_layout
    integer :: nspecies = 0
    integer :: nvariables = 0
    integer :: first_species = 1
    integer :: last_species = 0
    integer :: momentum(3) = 0
    integer :: total_energy = 0
  end type mc_state_layout

  public :: initialize_mc_state_layout
  public :: validate_mc_state_layout
  public :: mc_species_index

contains

  subroutine initialize_mc_state_layout(layout, nspecies)
    type(mc_state_layout), intent(out) :: layout
    integer, intent(in) :: nspecies

    if (nspecies < 1) then
      error stop 'multicomponent state layout requires at least one species'
    end if

    layout%nspecies = nspecies
    layout%nvariables = nspecies + 4
    layout%first_species = 1
    layout%last_species = nspecies
    layout%momentum = [nspecies+1, nspecies+2, nspecies+3]
    layout%total_energy = nspecies + 4
    call validate_mc_state_layout(layout)
  end subroutine initialize_mc_state_layout

  subroutine validate_mc_state_layout(layout)
    type(mc_state_layout), intent(in) :: layout

    if (layout%nspecies < 1) then
      error stop 'invalid multicomponent species count'
    end if
    if (layout%nvariables /= layout%nspecies + 4) then
      error stop 'invalid multicomponent variable count'
    end if
    if (layout%first_species /= 1 .or. &
        layout%last_species /= layout%nspecies) then
      error stop 'invalid multicomponent species range'
    end if
    if (any(layout%momentum /= [layout%nspecies+1, &
        layout%nspecies+2, layout%nspecies+3])) then
      error stop 'invalid multicomponent momentum indices'
    end if
    if (layout%total_energy /= layout%nvariables) then
      error stop 'invalid multicomponent total-energy index'
    end if
  end subroutine validate_mc_state_layout

  integer function mc_species_index(layout, species) result(index)
    type(mc_state_layout), intent(in) :: layout
    integer, intent(in) :: species

    if (species < 1 .or. species > layout%nspecies) then
      error stop 'multicomponent species index is out of range'
    end if
    index = layout%first_species + species - 1
  end function mc_species_index

end module mod_mc_state_layout
