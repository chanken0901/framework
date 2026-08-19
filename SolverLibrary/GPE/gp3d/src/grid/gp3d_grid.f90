!> 周期直交格子、FFT用波数、MPI zスラブ分割を初期化する。
!> 実空間は端点を重複させず、周期長Lをn点で分割する。
module gp3d_grid
  use gp3d_types, only: dp, pi, gp3d_grid_t
  implicit none
  private

  public :: gp3d_grid_init
  public :: gp3d_grid_init_bounds

contains

  subroutine gp3d_grid_init(grid, nx, ny, nz, lx, ly, lz, rank, nprocs)
    type(gp3d_grid_t), intent(out) :: grid
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: lx, ly, lz
    integer, intent(in), optional :: rank, nprocs

    integer :: i

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) error stop "grid dimensions must be positive"
    if (lx <= 0.0_dp .or. ly <= 0.0_dp .or. lz <= 0.0_dp) error stop "domain lengths must be positive"

    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%lx = lx
    grid%ly = ly
    grid%lz = lz
    grid%dx = lx / real(nx, dp)
    grid%dy = ly / real(ny, dp)
    grid%dz = lz / real(nz, dp)
    call set_decomposition(grid, rank, nprocs)

    allocate(grid%x(nx), grid%y(ny), grid%z(nz))
    allocate(grid%kx(nx), grid%ky(ny), grid%kz(nz))

    do i = 1, nx
      grid%x(i) = (real(i - 1, dp) - 0.5_dp * real(nx, dp)) * grid%dx
      grid%kx(i) = wave_number(i, nx, lx)
    end do

    do i = 1, ny
      grid%y(i) = (real(i - 1, dp) - 0.5_dp * real(ny, dp)) * grid%dy
      grid%ky(i) = wave_number(i, ny, ly)
    end do

    do i = 1, nz
      grid%z(i) = (real(i - 1, dp) - 0.5_dp * real(nz, dp)) * grid%dz
      grid%kz(i) = wave_number(i, nz, lz)
    end do
  end subroutine gp3d_grid_init

  subroutine gp3d_grid_init_bounds(grid, nx, ny, nz, x_min, x_max, y_min, y_max, z_min, z_max, rank, nprocs)
    type(gp3d_grid_t), intent(out) :: grid
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    integer, intent(in), optional :: rank, nprocs

    integer :: i
    real(dp) :: lx, ly, lz

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) error stop "grid dimensions must be positive"
    if (x_max <= x_min .or. y_max <= y_min .or. z_max <= z_min) error stop "invalid grid bounds"

    lx = x_max - x_min
    ly = y_max - y_min
    lz = z_max - z_min

    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%lx = lx
    grid%ly = ly
    grid%lz = lz
    grid%dx = lx / real(nx, dp)
    grid%dy = ly / real(ny, dp)
    grid%dz = lz / real(nz, dp)
    call set_decomposition(grid, rank, nprocs)

    allocate(grid%x(nx), grid%y(ny), grid%z(nz))
    allocate(grid%kx(nx), grid%ky(ny), grid%kz(nz))

    do i = 1, nx
      grid%x(i) = x_min + (real(i, dp) - 0.5_dp) * grid%dx
      grid%kx(i) = wave_number(i, nx, lx)
    end do

    do i = 1, ny
      grid%y(i) = y_min + (real(i, dp) - 0.5_dp) * grid%dy
      grid%ky(i) = wave_number(i, ny, ly)
    end do

    do i = 1, nz
      grid%z(i) = z_min + (real(i, dp) - 0.5_dp) * grid%dz
      grid%kz(i) = wave_number(i, nz, lz)
    end do
  end subroutine gp3d_grid_init_bounds

  subroutine set_decomposition(grid, rank, nprocs)
    type(gp3d_grid_t), intent(inout) :: grid
    integer, intent(in), optional :: rank, nprocs
    integer :: rank_value, nprocs_value, base, rest

    rank_value = 0
    nprocs_value = 1
    if (present(rank)) rank_value = rank
    if (present(nprocs)) nprocs_value = nprocs

    if (nprocs_value <= 0) error stop "MPI process count must be positive"
    if (rank_value < 0 .or. rank_value >= nprocs_value) error stop "MPI rank is out of range"
    if (nprocs_value > grid%nz) error stop "MPI process count must not exceed nz"

    base = grid%nz / nprocs_value
    rest = mod(grid%nz, nprocs_value)
    if (rank_value < rest) then
      grid%local_nz = base + 1
      grid%k_start = rank_value * (base + 1) + 1
    else
      grid%local_nz = base
      grid%k_start = rest * (base + 1) + (rank_value - rest) * base + 1
    end if
    grid%k_end = grid%k_start + grid%local_nz - 1
    grid%rank = rank_value
    grid%nprocs = nprocs_value
  end subroutine set_decomposition

  pure real(dp) function wave_number(i, n, length) result(k)
    integer, intent(in) :: i, n
    real(dp), intent(in) :: length
    integer :: mode

    if (i - 1 <= n / 2) then
      mode = i - 1
    else
      mode = i - 1 - n
    end if
    k = 2.0_dp * pi * real(mode, dp) / length
  end function wave_number

end module gp3d_grid
