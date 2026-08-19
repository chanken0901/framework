module mod_convective_weno5z_roe
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_reconstruction_weno5z, only : reconstruct_weno5z_left, &
    reconstruct_weno5z_right
  use mod_riemann_roe, only : rotate_conserved_to_normal, &
    rotate_flux_to_global, roe_eigensystem, roe_numerical_flux
  implicit none
  private

  public :: compute_weno5z_roe_flux
  public :: compute_weno5z_roe_face_flux
  public :: validate_weno5z_roe_scheme
  public :: weno5z_roe_required_ghost_cells

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
    call roe_numerical_flux(left_state, right_state, nse, normal_flux)
    call rotate_flux_to_global(normal_flux, direction, flux)
  end subroutine compute_weno5z_roe_face_flux

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
