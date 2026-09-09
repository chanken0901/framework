module mod_convective_weno5z_roe
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_reconstruction_weno5z, only : reconstruct_weno5z_left, &
    reconstruct_weno5z_right
  use mod_riemann_roe, only : rotate_conserved_to_normal, &
    rotate_flux_to_global, roe_eigensystem, roe_numerical_flux, euler_physical_flux
  implicit none
  private

  public :: compute_weno5z_roe_flux
  public :: compute_weno5z_roe_face_flux
  public :: validate_weno5z_roe_scheme
  public :: weno5z_roe_required_ghost_cells
  public :: limit_weno_state

contains

  subroutine compute_weno5z_roe_flux(q, fface, direction, sim, nse, &
      js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: direction, js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: fface(0:, js-1:, ks-1:, :)
    integer :: i, j, k

    if (nse%nv /= 5) error stop 'WENO5-Z/Roe flux requires five variables'
    if (sim%nghost < weno5z_roe_required_ghost_cells()) then
      error stop 'WENO5-Z/Roe flux requires three ghost cells'
    end if

    select case (direction)
    case (1:3)
      !$OMP DO collapse(2) schedule(static)
      do k = ks-1, ke
        do j = js-1, je
          do i = 0, sim%nx
            call compute_weno5z_roe_face_flux(q, i, j, k, direction, &
              sim, nse, &
              js, ks, fface(i,j,k,1:5))
          end do
        end do
      end do
      !$OMP END DO
    case default
      error stop 'convective flux direction must be 1, 2, or 3'
    end select
  end subroutine compute_weno5z_roe_flux

  pure subroutine compute_weno5z_roe_face_flux(q, i, j, k, direction, &
      sim, nse, js, ks, flux)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, direction, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: flux(5)
    real(dp) :: center_left(5), center_right(5)
    real(dp) :: stencil_state(5), characteristic_stencil(5)
    real(dp) :: left_characteristic(5), right_characteristic(5)
    real(dp) :: left_state(5), right_state(5), normal_flux(5)
    real(dp) :: right_matrix(5,5), left_matrix(5,5), eigenvalue(5)
    real(dp) :: fl(5), fr(5), ul, ur, cl, cr, speed
    logical :: limited_left, limited_right
    integer :: characteristic, point

    call normal_state_at_offset(q, i, j, k, direction, 0, sim, js, ks, &
      center_left)
    call normal_state_at_offset(q, i, j, k, direction, 1, sim, js, ks, &
      center_right)
    call roe_eigensystem(center_left, center_right, nse, right_matrix, &
      left_matrix, eigenvalue)

    do characteristic = 1, 5
      do point = 1, 5
        call normal_state_at_offset(q, i, j, k, direction, point-3, &
          sim, js, ks, stencil_state)
        characteristic_stencil(point) = dot_product( &
          left_matrix(characteristic,:), stencil_state)
      end do
      left_characteristic(characteristic) = &
        reconstruct_weno5z_left(characteristic_stencil)

      do point = 1, 5
        call normal_state_at_offset(q, i, j, k, direction, point-2, &
          sim, js, ks, stencil_state)
        characteristic_stencil(point) = dot_product( &
          left_matrix(characteristic,:), stencil_state)
      end do
      right_characteristic(characteristic) = &
        reconstruct_weno5z_right(characteristic_stencil)
    end do

    left_state = matmul(right_matrix, left_characteristic)
    right_state = matmul(right_matrix, right_characteristic)
    call limit_weno_state(center_left, left_state, nse, limited_left)
    call limit_weno_state(center_right, right_state, nse, limited_right)
    ul = left_state(2)/left_state(1)
    ur = right_state(2)/right_state(1)
    cl = sqrt(nse%gamma*state_pressure(left_state,nse)/left_state(1))
    cr = sqrt(nse%gamma*state_pressure(right_state,nse)/right_state(1))
    if (limited_left .or. limited_right .or. ur-ul > 2.0_dp*min(cl,cr)) then
      call euler_physical_flux(left_state, nse, fl)
      call euler_physical_flux(right_state, nse, fr)
      speed = max(abs(ul)+cl, abs(ur)+cr)
      normal_flux = 0.5_dp*(fl+fr-speed*(right_state-left_state))
    else
      call roe_numerical_flux(left_state, right_state, nse, normal_flux)
    end if
    call rotate_flux_to_global(normal_flux, direction, flux)
  end subroutine compute_weno5z_roe_face_flux

  pure real(dp) function state_pressure(state, nse) result(p)
    real(dp), intent(in) :: state(5)
    type(nse_config), intent(in) :: nse
    p = (nse%gamma-1.0_dp)*(state(5)-0.5_dp*sum(state(2:4)**2)/state(1))
  end function

  pure logical function admissible_state(state, nse) result(valid)
    real(dp), intent(in) :: state(5)
    type(nse_config), intent(in) :: nse
    real(dp) :: p
    valid = .false.
    if (.not. all(ieee_is_finite(state))) return
    if (state(1) < nse%small_rho) return
    p = state_pressure(state,nse)
    valid = ieee_is_finite(p) .and. p >= nse%small_p
  end function

  pure subroutine limit_weno_state(center, state, nse, limited)
    real(dp), intent(in) :: center(5)
    real(dp), intent(inout) :: state(5)
    type(nse_config), intent(in) :: nse
    logical, intent(out) :: limited
    real(dp) :: original(5), low, high, theta
    integer :: iteration
    limited = .not. admissible_state(state,nse)
    if (.not. limited) return
    if (.not. all(ieee_is_finite(state))) then
      state = center
      return
    end if
    original = state
    low = 0.0_dp
    high = 1.0_dp
    do iteration = 1, 50
      theta = 0.5_dp*(low+high)
      state = center+theta*(original-center)
      if (admissible_state(state,nse)) then
        low = theta
      else
        high = theta
      end if
    end do
    state = center+(0.99_dp*low)*(original-center)
  end subroutine

  pure subroutine normal_state_at_offset(q, i, j, k, direction, offset, &
      sim, js, ks, state)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, direction, offset, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: state(5)
    real(dp) :: global_state(5)

    select case (direction)
    case (1)
      global_state = q(i+offset,j,k,1:5)
    case (2)
      global_state = q(i,j+offset,k,1:5)
    case (3)
      global_state = q(i,j,k+offset,1:5)
    case default
      global_state = 0.0_dp
    end select
    call rotate_conserved_to_normal(global_state, direction, state)
  end subroutine normal_state_at_offset

  subroutine validate_weno5z_roe_scheme(nse)
    type(nse_config), intent(in) :: nse

    if (nse%nv /= 5) then
      error stop 'WENO5-Z/Roe flux requires five conserved variables'
    end if
    if (nse%gamma <= 1.0_dp) error stop 'WENO5-Z/Roe requires gamma > 1'
    if (nse%small_rho <= 0.0_dp .or. nse%small_p <= 0.0_dp) then
      error stop 'WENO5-Z/Roe requires positive density and pressure floors'
    end if
  end subroutine validate_weno5z_roe_scheme

  integer function weno5z_roe_required_ghost_cells() result(nghost)
    nghost = 3
  end function weno5z_roe_required_ghost_cells

end module mod_convective_weno5z_roe
