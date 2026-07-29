!> GPEの初期波動関数へ量子渦、渦輪、渦タングル、Taylor-Green構造を設定する。
!> 振幅の渦芯形状と位相巻き込みを同時に与え、必要に応じて再現可能な位相ノイズを加える。
module gp3d_initial_conditions
  use gp3d_types, only: dp, pi, gp3d_grid_t, gp3d_state_t
  use gp3d_openmp, only: gp3d_openmp_active
  implicit none
  private

  public :: gp3d_set_uniform_vortex_line
  public :: gp3d_set_thomas_fermi_vortex_line
  public :: gp3d_imprint_vortex_line
  public :: gp3d_imprint_vortex_lines
  public :: gp3d_imprint_vortex_ring
  public :: gp3d_imprint_vortex_ring_oriented
  public :: gp3d_imprint_random_vortex_lines
  public :: gp3d_imprint_random_vortex_rings
  public :: gp3d_add_random_phase_noise
  public :: gp3d_set_quantum_taylor_green

contains

  subroutine gp3d_set_quantum_taylor_green(state, grid, healing_length, winding, density0)
    ! 周期Taylor-Green対称性を持つ4本組の素渦を各格子点で合成する。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    real(dp), intent(in) :: healing_length, density0
    integer, intent(in) :: winding

    integer :: i, j, k, kg
    real(dp) :: x_angle, y_angle, z_angle, lambda, mu, root_half
    real(dp) :: x_origin, y_origin, z_origin, cos_z_scale
    complex(dp) :: psi_four

    if (healing_length <= 0.0_dp) error stop "Taylor-Green healing length must be positive"
    if (density0 <= 0.0_dp) error stop "Taylor-Green density must be positive"
    if (winding <= 0) error stop "Taylor-Green winding must be positive"

    x_origin = grid%x(1) - 0.5_dp * grid%dx
    y_origin = grid%y(1) - 0.5_dp * grid%dy
    z_origin = grid%z(1) - 0.5_dp * grid%dz
    root_half = 1.0_dp / sqrt(2.0_dp)

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(kg, x_angle, y_angle, z_angle, cos_z_scale, lambda, mu, psi_four)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          z_angle = 2.0_dp * pi * (grid%z(kg) - z_origin) / grid%lz - pi
          cos_z_scale = sqrt(2.0_dp * abs(cos(z_angle)))
          y_angle = 2.0_dp * pi * (grid%y(j) - y_origin) / grid%ly - pi
          x_angle = 2.0_dp * pi * (grid%x(i) - x_origin) / grid%lx - pi
          lambda = cos(x_angle) * cos_z_scale
          mu = cos(y_angle) * cos_z_scale * sign(1.0_dp, cos(z_angle))

          psi_four = elementary_vortex(lambda - root_half, mu, healing_length) * &
            elementary_vortex(lambda, mu - root_half, healing_length) * &
            elementary_vortex(lambda + root_half, mu, healing_length) * &
            elementary_vortex(lambda, mu + root_half, healing_length)
          state%psi(i,j,k) = sqrt(density0) * psi_four**winding
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine gp3d_set_quantum_taylor_green

  subroutine gp3d_set_uniform_vortex_line(state, grid, charge, x0, y0, healing_length, density0)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: charge
    real(dp), intent(in) :: x0, y0, healing_length, density0

    state%psi = cmplx(sqrt(density0), 0.0_dp, kind=dp)
    call gp3d_imprint_vortex_line(state, grid, charge, x0, y0, healing_length)
  end subroutine gp3d_set_uniform_vortex_line

  subroutine gp3d_set_thomas_fermi_vortex_line(state, grid, charge, x0, y0, healing_length, mu, g)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: charge
    real(dp), intent(in) :: x0, y0, healing_length, mu, g

    integer :: i, j, k
    real(dp) :: density

    if (g <= 0.0_dp) error stop "Thomas-Fermi initial condition requires g > 0"

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) private(density)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          density = max((mu - state%potential(i,j,k)) / g, 0.0_dp)
          state%psi(i,j,k) = cmplx(sqrt(density), 0.0_dp, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

    call gp3d_imprint_vortex_line(state, grid, charge, x0, y0, healing_length)
  end subroutine gp3d_set_thomas_fermi_vortex_line

  subroutine gp3d_imprint_vortex_line(state, grid, charge, x0, y0, healing_length)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: charge
    real(dp), intent(in) :: x0, y0, healing_length

    integer :: charges(1)
    real(dp) :: x_positions(1), y_positions(1)

    charges(1) = charge
    x_positions(1) = x0
    y_positions(1) = y0
    call gp3d_imprint_vortex_lines(state, grid, charges, x_positions, y_positions, healing_length)
  end subroutine gp3d_imprint_vortex_line

  subroutine gp3d_imprint_vortex_lines(state, grid, charges, x0, y0, healing_length)
    ! 既存の背景波動関数へ、z方向直線渦の芯振幅と位相因子を乗算する。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: charges(:)
    real(dp), intent(in) :: x0(:), y0(:), healing_length

    integer :: i, j, k, n
    real(dp) :: dx, dy, radius, phase, core
    complex(dp) :: factor

    if (size(charges) /= size(x0) .or. size(charges) /= size(y0)) then
      error stop "vortex arrays must have the same length"
    end if
    if (healing_length <= 0.0_dp) error stop "healing length must be positive"

    do n = 1, size(charges)
      if (charges(n) == 0) cycle

      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
      !$omp& private(dx, dy, radius, phase, core, factor)
      do k = 1, grid%local_nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            dx = grid%x(i) - x0(n)
            dy = grid%y(j) - y0(n)
            radius = sqrt(dx * dx + dy * dy)
            phase = real(charges(n), dp) * atan2(dy, dx)
            core = vortex_core_profile(radius, healing_length, abs(charges(n)))
            factor = core * cmplx(cos(phase), sin(phase), kind=dp)
            state%psi(i,j,k) = state%psi(i,j,k) * factor
          end do
        end do
      end do
      !$omp end parallel do
    end do
  end subroutine gp3d_imprint_vortex_lines

  subroutine gp3d_imprint_vortex_ring(state, grid, charge, ring_radius, z0, healing_length)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: charge
    real(dp), intent(in) :: ring_radius, z0, healing_length

    integer :: i, j, k, kg
    real(dp) :: rho, dr, dz, distance, phase, core
    complex(dp) :: factor

    if (charge == 0) return
    if (ring_radius <= 0.0_dp) error stop "ring radius must be positive"
    if (healing_length <= 0.0_dp) error stop "healing length must be positive"

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(kg, rho, dr, dz, distance, phase, core, factor)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          rho = sqrt(grid%x(i)**2 + grid%y(j)**2)
          dr = rho - ring_radius
          dz = grid%z(kg) - z0
          distance = sqrt(dr * dr + dz * dz)
          phase = real(charge, dp) * atan2(dz, dr)
          core = vortex_core_profile(distance, healing_length, abs(charge))
          factor = core * cmplx(cos(phase), sin(phase), kind=dp)
          state%psi(i,j,k) = state%psi(i,j,k) * factor
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine gp3d_imprint_vortex_ring

  subroutine gp3d_imprint_vortex_ring_oriented(state, grid, charge, center, normal, ring_radius, healing_length)
    ! 任意の中心・法線を持つ局所座標系へ写像し、3次元の渦輪を刻印する。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: charge
    real(dp), intent(in) :: center(3), normal(3), ring_radius, healing_length

    integer :: i, j, k, kg
    real(dp) :: nvec(3), rvec(3), parallel, in_plane2, rho, dr, distance, phase, core, norm_normal
    complex(dp) :: factor

    if (charge == 0) return
    if (ring_radius <= 0.0_dp) error stop "ring radius must be positive"
    if (healing_length <= 0.0_dp) error stop "healing length must be positive"

    norm_normal = sqrt(sum(normal * normal))
    if (norm_normal <= 0.0_dp) error stop "ring normal must be nonzero"
    nvec = normal / norm_normal

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(kg, rvec, parallel, in_plane2, rho, dr, distance, phase, core, factor)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          rvec = [grid%x(i) - center(1), grid%y(j) - center(2), grid%z(kg) - center(3)]
          parallel = sum(rvec * nvec)
          in_plane2 = max(sum(rvec * rvec) - parallel * parallel, 0.0_dp)
          rho = sqrt(in_plane2)
          dr = rho - ring_radius
          distance = sqrt(dr * dr + parallel * parallel)
          phase = real(charge, dp) * atan2(parallel, dr)
          core = vortex_core_profile(distance, healing_length, abs(charge))
          factor = core * cmplx(cos(phase), sin(phase), kind=dp)
          state%psi(i,j,k) = state%psi(i,j,k) * factor
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine gp3d_imprint_vortex_ring_oriented

  subroutine gp3d_imprint_random_vortex_lines(state, grid, nlines, healing_length, seed)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: nlines
    real(dp), intent(in) :: healing_length
    integer, intent(in), optional :: seed

    integer, allocatable :: charges(:)
    real(dp), allocatable :: x0(:), y0(:)
    integer :: n
    real(dp) :: rx, ry, margin

    if (nlines <= 0) return
    if (healing_length <= 0.0_dp) error stop "healing length must be positive"

    if (present(seed)) call set_random_seed(seed)

    allocate(charges(nlines), x0(nlines), y0(nlines))
    margin = 2.0_dp * healing_length

    do n = 1, nlines
      call random_number(rx)
      call random_number(ry)
      x0(n) = (minval(grid%x) + margin) + rx * max(grid%lx - 2.0_dp * margin, grid%dx)
      y0(n) = (minval(grid%y) + margin) + ry * max(grid%ly - 2.0_dp * margin, grid%dy)
      if (mod(n, 2) == 0) then
        charges(n) = -1
      else
        charges(n) = 1
      end if
    end do

    call gp3d_imprint_vortex_lines(state, grid, charges, x0, y0, healing_length)

    deallocate(charges, x0, y0)
  end subroutine gp3d_imprint_random_vortex_lines

  subroutine gp3d_imprint_random_vortex_rings(state, grid, nrings, radius_min, radius_max, healing_length, seed)
    ! 位置、半径、法線を乱数で選び、符号を交互にした渦輪群を作る。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: nrings
    real(dp), intent(in) :: radius_min, radius_max, healing_length
    integer, intent(in), optional :: seed

    integer :: n, charge
    real(dp) :: rx, ry, rz, rr, ru, rphi
    real(dp) :: center(3), normal(3), radius, margin, zeta, azimuth

    if (nrings <= 0) return
    if (radius_min <= 0.0_dp .or. radius_max < radius_min) error stop "invalid random ring radius range"
    if (healing_length <= 0.0_dp) error stop "healing length must be positive"

    if (present(seed)) call set_random_seed(seed)

    margin = radius_max + 4.0_dp * healing_length
    do n = 1, nrings
      call random_number(rx)
      call random_number(ry)
      call random_number(rz)
      call random_number(rr)
      call random_number(ru)
      call random_number(rphi)

      center(1) = (minval(grid%x) + margin) + rx * max(grid%lx - 2.0_dp * margin, grid%dx)
      center(2) = (minval(grid%y) + margin) + ry * max(grid%ly - 2.0_dp * margin, grid%dy)
      center(3) = (minval(grid%z) + margin) + rz * max(grid%lz - 2.0_dp * margin, grid%dz)
      radius = radius_min + rr * (radius_max - radius_min)

      zeta = 2.0_dp * ru - 1.0_dp
      azimuth = 2.0_dp * pi * rphi
      normal = [sqrt(max(1.0_dp - zeta * zeta, 0.0_dp)) * cos(azimuth), &
                sqrt(max(1.0_dp - zeta * zeta, 0.0_dp)) * sin(azimuth), zeta]

      if (mod(n, 2) == 0) then
        charge = -1
      else
        charge = 1
      end if
      call gp3d_imprint_vortex_ring_oriented(state, grid, charge, center, normal, radius, healing_length)
    end do
  end subroutine gp3d_imprint_random_vortex_rings

  subroutine gp3d_add_random_phase_noise(state, amplitude, grid, seed)
    ! MPI分割数によらず同じ全体格子点に同じ位相を与える決定論的ノイズを使う。
    type(gp3d_state_t), intent(inout) :: state
    real(dp), intent(in) :: amplitude
    type(gp3d_grid_t), intent(in), optional :: grid
    integer, intent(in), optional :: seed

    integer :: i, j, k, kg, seed_value
    real(dp) :: random_value, phase, key

    if (amplitude < 0.0_dp) error stop "phase noise amplitude must be non-negative"

    seed_value = 0
    if (present(seed)) seed_value = seed
    do k = 1, size(state%psi, 3)
      kg = k
      if (present(grid)) kg = grid%k_start + k - 1
      do j = 1, size(state%psi, 2)
        do i = 1, size(state%psi, 1)
          if (present(grid)) then
            key = 12.9898_dp * real(i, dp) + 78.233_dp * real(j, dp) + &
              37.719_dp * real(kg, dp) + 0.12345_dp * real(seed_value, dp)
            random_value = modulo(sin(key) * 43758.5453_dp, 1.0_dp)
          else
            call random_number(random_value)
          end if
          phase = amplitude * (2.0_dp * random_value - 1.0_dp)
          state%psi(i,j,k) = state%psi(i,j,k) * cmplx(cos(phase), sin(phase), kind=dp)
        end do
      end do
    end do
  end subroutine gp3d_add_random_phase_noise

  pure real(dp) function vortex_core_profile(radius, healing_length, charge_abs) result(core)
    real(dp), intent(in) :: radius, healing_length
    integer, intent(in) :: charge_abs
    real(dp) :: scaled_radius

    if (charge_abs <= 0) then
      core = 1.0_dp
      return
    end if

    scaled_radius = radius / healing_length
    core = (scaled_radius / sqrt(1.0_dp + scaled_radius * scaled_radius))**charge_abs
  end function vortex_core_profile

  pure complex(dp) function elementary_vortex(lambda, mu, healing_length) result(psi)
    real(dp), intent(in) :: lambda, mu, healing_length
    real(dp) :: radius, core

    radius = sqrt(lambda * lambda + mu * mu)
    if (radius <= tiny(1.0_dp)) then
      psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
      return
    end if

    core = tanh(radius / (sqrt(2.0_dp) * healing_length))
    psi = core * cmplx(lambda / radius, mu / radius, kind=dp)
  end function elementary_vortex

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: values(:)

    call random_seed(size=n)
    allocate(values(n))
    do i = 1, n
      values(i) = seed + 104729 * (i - 1)
    end do
    call random_seed(put=values)
    deallocate(values)
  end subroutine set_random_seed

end module gp3d_initial_conditions
