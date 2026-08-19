module mod_nse_spatial_operator
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_grid_fvm, only : vol, area_x, area_y, area_z
  use mod_nse_boundary, only : apply_nse_boundary, validate_boundary_scheme, &
    boundary_required_ghost_cells
  use mod_convective_scheme, only : compute_convective_flux, &
    validate_convective_scheme, convective_required_ghost_cells
  use mod_viscous_scheme, only : add_viscous_rhs, validate_viscous_scheme, &
    viscous_required_ghost_cells
  use mod_nse_forcing, only : add_nse_forcing_rhs, validate_nse_forcing
  implicit none
  private

  public :: compute_nse_rhs
  public :: validate_nse_spatial_configuration
  public :: required_nse_ghost_cells

contains

  subroutine compute_nse_rhs(q, rhs, fface, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(out) :: rhs(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(out) :: fface(0:, js-1:, ks-1:, :)
    integer :: i, j, k
    real(dp) :: flux_divergence(nse%nv)

    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          rhs(i,j,k,:) = 0.0_dp
        end do
      end do
    end do
    !$OMP END DO

    call apply_nse_boundary(q, sim, nse, js, je, ks, ke)

    call compute_convective_flux(q, fface, 1, sim, nse, js, je, ks, ke)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          flux_divergence = (area_x(i,j,k)*fface(i,j,k,:) - &
            area_x(i-1,j,k)*fface(i-1,j,k,:)) / vol(i,j,k)
          rhs(i,j,k,:) = rhs(i,j,k,:) - flux_divergence
        end do
      end do
    end do
    !$OMP END DO

    call compute_convective_flux(q, fface, 2, sim, nse, js, je, ks, ke)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          flux_divergence = (area_y(i,j,k)*fface(i,j,k,:) - &
            area_y(i,j-1,k)*fface(i,j-1,k,:)) / vol(i,j,k)
          rhs(i,j,k,:) = rhs(i,j,k,:) - flux_divergence
        end do
      end do
    end do
    !$OMP END DO

    call compute_convective_flux(q, fface, 3, sim, nse, js, je, ks, ke)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          flux_divergence = (area_z(i,j,k)*fface(i,j,k,:) - &
            area_z(i,j,k-1)*fface(i,j,k-1,:)) / vol(i,j,k)
          rhs(i,j,k,:) = rhs(i,j,k,:) - flux_divergence
        end do
      end do
    end do
    !$OMP END DO

    call add_viscous_rhs(q, rhs, sim, nse, js, je, ks, ke)
    call add_nse_forcing_rhs(q, rhs, sim, nse, js, je, ks, ke)
  end subroutine compute_nse_rhs

  subroutine validate_nse_spatial_configuration(sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer :: required

    call validate_boundary_scheme(sim, nse)
    call validate_convective_scheme(nse)
    call validate_viscous_scheme(nse)
    call validate_nse_forcing(nse)

    required = required_nse_ghost_cells()
    if (sim%nghost < required) then
      write(*,'(A,I0,A,I0)') 'ERROR: configured ghost cells = ', sim%nghost, &
        ', but selected numerical schemes require at least ', required
      error stop
    end if
  end subroutine validate_nse_spatial_configuration

  integer function required_nse_ghost_cells() result(nghost)
    nghost = max(boundary_required_ghost_cells(), &
      convective_required_ghost_cells(), viscous_required_ghost_cells())
  end function required_nse_ghost_cells

end module mod_nse_spatial_operator
