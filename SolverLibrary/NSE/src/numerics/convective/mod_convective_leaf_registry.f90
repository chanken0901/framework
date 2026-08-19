module mod_convective_leaf_registry
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_convective_keep, only : compute_keep_face_flux
  use mod_convective_weno5z_roe, only : compute_weno5z_roe_face_flux, &
    validate_weno5z_roe_scheme
  implicit none
  private

  public :: compute_leaf_face_flux
  public :: validate_leaf_scheme
  public :: is_leaf_scheme
  public :: leaf_required_ghost_cells

contains

  pure subroutine compute_leaf_face_flux(scheme, q, i, j, k, direction, &
      sim, nse, js, ks, flux)
    character(len=*), intent(in) :: scheme
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, direction, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: flux(5)

    select case (trim(adjustl(scheme)))
    case ('keep2')
      call compute_keep_face_flux(q, i, j, k, direction, 2, sim, nse, &
        js, ks, flux)
    case ('keep6')
      call compute_keep_face_flux(q, i, j, k, direction, 6, sim, nse, &
        js, ks, flux)
    case ('weno5z_roe')
      call compute_weno5z_roe_face_flux(q, i, j, k, direction, sim, nse, &
        js, ks, flux)
    case default
      flux = 0.0_dp
    end select
  end subroutine compute_leaf_face_flux

  subroutine validate_leaf_scheme(scheme, nse, role)
    character(len=*), intent(in) :: scheme, role
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(scheme)))
    case ('keep2', 'keep6')
      if (nse%nv /= 5) error stop 'KEEP flux requires five variables'
    case ('weno5z_roe')
      call validate_weno5z_roe_scheme(nse)
    case default
      write(*,'(A,A,A,A,A)') 'ERROR: unsupported hybrid ', trim(role), &
        ' scheme "', trim(adjustl(scheme)), &
        '"; use keep2, keep6, or weno5z_roe'
      error stop
    end select
  end subroutine validate_leaf_scheme

  pure logical function is_leaf_scheme(scheme) result(is_leaf)
    character(len=*), intent(in) :: scheme

    select case (trim(adjustl(scheme)))
    case ('keep2', 'keep6', 'weno5z_roe')
      is_leaf = .true.
    case default
      is_leaf = .false.
    end select
  end function is_leaf_scheme

  integer function leaf_required_ghost_cells(scheme) result(nghost)
    character(len=*), intent(in) :: scheme

    if (is_leaf_scheme(scheme)) then
      nghost = 3
    else
      nghost = 0
    end if
  end function leaf_required_ghost_cells

end module mod_convective_leaf_registry
