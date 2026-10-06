module mod_grid_fvm
  use mod_precision,     only : dp
  use mod_common_config, only : simulation_config
  use mod_grid_axis, only: grid_axis,build_sinh_axis
  implicit none

  private

  public :: build_uniform_grid
  public :: x_edge, x_cell
  public :: y_edge, y_cell
  public :: z_edge, z_cell
  public :: vol, area_x, area_y, area_z
  public :: dx_min, dy_min, dz_min
  public :: axis_x,axis_y,axis_z,build_stretched_grid
  type(grid_axis) :: axis_x,axis_y,axis_z

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

    if(trim(sim%grid_mapping)=='sinh') then
      call build_stretched_grid(sim,js,je,ks,ke)
      return
    end if
    if(trim(sim%grid_mapping)/='uniform') error stop 'Unsupported grid mapping'

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

  subroutine build_stretched_grid(sim,js,je,ks,ke)
    type(simulation_config), intent(inout) :: sim
    integer, intent(in) :: js,je,ks,ke
    integer :: i,j,k,g
    real(dp) :: dx,dy,dz
    g=sim%nghost
    if(js<1.or.je>sim%ny.or.ks<1.or.ke>sim%nz.or.js>je.or.ks>ke) error stop 'Invalid grid partition'
    if(allocated(vol)) error stop 'Grid already allocated'
    ! Global geometry is generated deterministically before selecting local indices.
    ! Physical-boundary ghost widths are mirrored; periodic mapping is not connected yet.
    call build_sinh_axis(sim%nx,g,sim%x_min,sim%x_max,sim%grid_stretch(1),.false.,axis_x)
    call build_sinh_axis(sim%ny,g,sim%y_min,sim%y_max,sim%grid_stretch(2),.false.,axis_y)
    call build_sinh_axis(sim%nz,g,sim%z_min,sim%z_max,sim%grid_stretch(3),.false.,axis_z)
    allocate(x_edge(-g:sim%nx+g,js-g-1:je+g,ks-g-1:ke+g))
    allocate(y_edge,mold=x_edge);allocate(z_edge,mold=x_edge)
    allocate(x_cell(1-g:sim%nx+g,js-g:je+g,ks-g:ke+g))
    allocate(y_cell,mold=x_cell);allocate(z_cell,mold=x_cell)
    allocate(vol(0:sim%nx,js-1:je,ks-1:ke))
    allocate(area_x,mold=vol);allocate(area_y,mold=vol);allocate(area_z,mold=vol)
    do k=ks-g-1,ke+g
      do j=js-g-1,je+g
        do i=-g,sim%nx+g
          x_edge(i,j,k)=axis_x%edge(i)
          y_edge(i,j,k)=axis_y%edge(j)
          z_edge(i,j,k)=axis_z%edge(k)
        end do
      end do
    end do
    do k=ks-g,ke+g
      do j=js-g,je+g
        do i=1-g,sim%nx+g
          x_cell(i,j,k)=axis_x%center(i)
          y_cell(i,j,k)=axis_y%center(j)
          z_cell(i,j,k)=axis_z%center(k)
        end do
      end do
    end do
    do k=ks-1,ke
      do j=js-1,je
        do i=0,sim%nx
          dx=axis_x%width(i);dy=axis_y%width(j);dz=axis_z%width(k)
          vol(i,j,k)=dx*dy*dz
          area_x(i,j,k)=dy*dz;area_y(i,j,k)=dx*dz;area_z(i,j,k)=dx*dy
        end do
      end do
    end do
    dx_min=axis_x%minimum_width;dy_min=axis_y%minimum_width;dz_min=axis_z%minimum_width
    sim%dx=dx_min;sim%dy=dy_min;sim%dz=dz_min
  end subroutine

end module mod_grid_fvm
