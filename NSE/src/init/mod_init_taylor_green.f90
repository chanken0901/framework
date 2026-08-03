module mod_init_taylor_green
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_grid_fvm, only : x_cell, y_cell, z_cell
  implicit none
  private

  public :: initialize_taylor_green

contains

  subroutine initialize_taylor_green(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    integer :: i, j, k
    real(dp) :: rho, u, v, w, p
    real(dp) :: c1, c2, x2, y2, z2

    if (nse%nv /= 5) error stop 'Taylor-Green initialization requires five conserved variables'

    c1 = 1.0_dp / nse%gamma
    c2 = (nse%rho0 * nse%mach * nse%mach) / 16.0_dp

    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          x2 = 2.0_dp * x_cell(i,j,k)
          y2 = 2.0_dp * y_cell(i,j,k)
          z2 = 2.0_dp * z_cell(i,j,k)

          rho = nse%rho0
          u = nse%mach * sin(x_cell(i,j,k)) * cos(y_cell(i,j,k)) * cos(z_cell(i,j,k))
          v = -nse%mach * cos(x_cell(i,j,k)) * sin(y_cell(i,j,k)) * cos(z_cell(i,j,k))
          w = 0.0_dp
          p = c1 + c2 * (cos(x2) + cos(y2)) * (cos(z2) + 2.0_dp)

          q(i,j,k,1) = rho
          q(i,j,k,2) = rho * u
          q(i,j,k,3) = rho * v
          q(i,j,k,4) = rho * w
          q(i,j,k,5) = p / (nse%gamma - 1.0_dp) + &
            0.5_dp * rho * (u*u + v*v + w*w)
        end do
      end do
    end do
  end subroutine initialize_taylor_green

end module mod_init_taylor_green
