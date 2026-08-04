module mod_convective_scheme
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_convective_keep, only : compute_keep_flux, validate_keep_scheme, &
    keep_required_ghost_cells
  use mod_convective_weno5z_roe, only : compute_weno5z_roe_flux, &
    validate_weno5z_roe_scheme, weno5z_roe_required_ghost_cells
  implicit none
  private

  public :: compute_convective_flux
  public :: validate_convective_scheme
  public :: convective_required_ghost_cells

contains

  subroutine compute_convective_flux(q, fface, direction, sim, nse, &
      js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: direction, js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: fface(0:, js-1:, ks-1:, :)

    select case (trim(adjustl(nse%convective_scheme)))
    case ('keep2', 'keep6')
      call compute_keep_flux(q, fface, direction, sim, nse, js, je, ks, ke)
    case ('weno5z_roe')
      call compute_weno5z_roe_flux(q, fface, direction, sim, nse, &
        js, je, ks, ke)
    case default
      write(*,'(A,A,A)') 'ERROR: unsupported convective scheme "', &
        trim(adjustl(nse%convective_scheme)), &
        '"; use keep2, keep6, or weno5z_roe'
      error stop
    end select
  end subroutine compute_convective_flux

  subroutine validate_convective_scheme(nse)
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(nse%convective_scheme)))
    case ('keep2', 'keep6')
      call validate_keep_scheme(nse)
    case ('weno5z_roe')
      call validate_weno5z_roe_scheme(nse)
    case default
      write(*,'(A,A,A)') 'ERROR: unsupported convective scheme "', &
        trim(adjustl(nse%convective_scheme)), &
        '"; use keep2, keep6, or weno5z_roe'
      error stop
    end select
  end subroutine validate_convective_scheme

  integer function convective_required_ghost_cells() result(nghost)
    ! Both current families need three layers. A future KEEP/Roe hybrid can
    ! remain behind this API and return the maximum required by its branches.
    nghost = max(keep_required_ghost_cells(), &
      weno5z_roe_required_ghost_cells())
  end function convective_required_ghost_cells

end module mod_convective_scheme
