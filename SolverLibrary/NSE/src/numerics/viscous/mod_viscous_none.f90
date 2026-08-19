module mod_viscous_scheme
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  implicit none
  private

  public :: add_viscous_rhs
  public :: validate_viscous_scheme
  public :: viscous_required_ghost_cells
  public :: viscous_scheme_name
  public :: viscous_dt_limit

contains

  subroutine add_viscous_rhs(q, rhs, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)

    ! The "none" backend deliberately leaves the convective RHS unchanged.
    if (size(q,4) /= nse%nv .or. size(rhs,4) /= nse%nv) then
      error stop 'viscous backend received an inconsistent state size'
    end if
    if (je < js .or. ke < ks) error stop 'viscous backend received an empty domain'
  end subroutine add_viscous_rhs

  subroutine viscous_dt_limit(q, dt_limit, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(out) :: dt_limit

    if (size(q,4) /= nse%nv .or. je < js .or. ke < ks) then
      error stop 'viscous time-step backend received an inconsistent domain'
    end if
    dt_limit = huge(1.0_dp)
  end subroutine viscous_dt_limit

  subroutine validate_viscous_scheme(nse)
    type(nse_config), intent(in) :: nse

    if (trim(adjustl(nse%viscous_scheme)) /= viscous_scheme_name()) then
      write(*,'(A,A,A,A)') 'ERROR: executable contains viscous scheme "', &
        viscous_scheme_name(), '", but input requested "', &
        trim(adjustl(nse%viscous_scheme)) // '"'
      error stop
    end if
  end subroutine validate_viscous_scheme

  integer function viscous_required_ghost_cells() result(nghost)
    nghost = 0
  end function viscous_required_ghost_cells

  pure function viscous_scheme_name() result(name)
    character(len=32) :: name
    name = 'none'
  end function viscous_scheme_name

end module mod_viscous_scheme
