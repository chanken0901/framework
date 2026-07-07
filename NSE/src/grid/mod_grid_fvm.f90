module mod_grid_fvm
  use mod_precision,     only : dp
  use mod_common_config, only : simulation_config
  implicit none

  private

  public :: build_uniform_grid
  public :: x_edge, x_cell
  public :: y_edge, y_cell
  public :: z_edge, z_cell
  public :: vol, area_x, area_y, area_z
  public :: dx_min, dy_min, dz_min

  real(dp), allocatable :: x_edge(:,:,:)
  real(dp), allocatable :: x_cell(:,:,:)
  real(dp), allocatable :: y_edge(:,:,:)
  real(dp), allocatable :: y_cell(:,:,:)
  real(dp), allocatable :: z_edge(:,:,:)
  real(dp), allocatable :: z_cell(:,:,:)

  real(dp), allocatable :: vol(:,:,:)
  real(dp), allocatable :: area_x(:,:,:)
  real(dp), allocatable :: area_y(:,:,:)
  real(dp), allocatable :: area_z(:,:,:)

  real(dp) :: dx_min, dy_min, dz_min

contains

  subroutine build_uniform_grid(sim, js, je, ks, ke)
    type(simulation_config), intent(inout) :: sim
    integer, intent(in) :: js, je, ks, ke

    integer :: i, j, k
    real(dp) :: dx, dy, dz

    allocate(x_edge(-1:sim%nx, js-2:je, ks-2:ke))
    allocate(y_edge(-1:sim%nx, js-2:je, ks-2:ke))
    allocate(z_edge(-1:sim%nx, js-2:je, ks-2:ke))

    allocate(x_cell(-2:sim%nx, js-2:je, ks-2:ke))
    allocate(y_cell(-2:sim%nx, js-2:je, ks-2:ke))
    allocate(z_cell(-2:sim%nx, js-2:je, ks-2:ke))

    allocate(vol   (0:sim%nx, js-1:je, ks-1:ke))
    allocate(area_x(0:sim%nx, js-1:je, ks-1:ke))
    allocate(area_y(0:sim%nx, js-1:je, ks-1:ke))
    allocate(area_z(0:sim%nx, js-1:je, ks-1:ke))

    dx_min = huge(1.0_dp)
    dy_min = huge(1.0_dp)
    dz_min = huge(1.0_dp)

    do k = ks-2, ke
    do j = js-2, je
    do i = -1, sim%nx
      x_edge(i,j,k) = sim%x_min + (sim%x_max - sim%x_min) * real(i,dp) / real(sim%nx,dp)
      y_edge(i,j,k) = sim%y_min + (sim%y_max - sim%y_min) * real(j,dp) / real(sim%ny,dp)
      z_edge(i,j,k) = sim%z_min + (sim%z_max - sim%z_min) * real(k,dp) / real(sim%nz,dp)
    end do
    end do
    end do

    do k = ks-1, ke
    do j = js-1, je
    do i = 0, sim%nx
      x_cell(i,j,k) = 0.5_dp * (x_edge(i-1,j,k) + x_edge(i,j,k))
      y_cell(i,j,k) = 0.5_dp * (y_edge(i,j-1,k) + y_edge(i,j,k))
      z_cell(i,j,k) = 0.5_dp * (z_edge(i,j,k-1) + z_edge(i,j,k))

      dx = x_edge(i,j,k) - x_edge(i-1,j,k)
      dy = y_edge(i,j,k) - y_edge(i,j-1,k)
      dz = z_edge(i,j,k) - z_edge(i,j,k-1)

      dx_min = min(dx_min, dx)
      dy_min = min(dy_min, dy)
      dz_min = min(dz_min, dz)

      area_x(i,j,k) = dy * dz
      area_y(i,j,k) = dz * dx
      area_z(i,j,k) = dx * dy

      vol(i,j,k) = dx * dy * dz
    end do
    end do
    end do

    sim%dx = dx_min
    sim%dy = dy_min
    sim%dz = dz_min

  end subroutine build_uniform_grid

end module mod_grid_fvm