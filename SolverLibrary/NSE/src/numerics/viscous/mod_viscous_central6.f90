module mod_viscous_scheme
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  implicit none
  private

  real(dp), parameter :: d1_positive(3) = [ &
    3.0_dp/4.0_dp, -3.0_dp/20.0_dp, 1.0_dp/60.0_dp ]
  real(dp), parameter :: d2_spectral_radius = 272.0_dp / 45.0_dp
  real(dp), parameter :: diffusion_stability_radius = 2.0_dp
  real(dp), allocatable, save :: primitive(:,:,:,:)

  public :: add_viscous_rhs
  public :: validate_viscous_scheme
  public :: viscous_required_ghost_cells
  public :: viscous_scheme_name
  public :: viscous_dt_limit

contains

  subroutine add_viscous_rhs(q, rhs, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    integer :: i, j, k
    real(dp) :: rho, u, v, w, pressure
    real(dp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
    real(dp) :: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz
    real(dp) :: uxy, uxz, vxy, vyz, wxz, wyz
    real(dp) :: momentum_x, momentum_y, momentum_z
    real(dp) :: div_velocity, tau_xx, tau_yy, tau_zz
    real(dp) :: tau_xy, tau_xz, tau_yz, dissipation
    real(dp) :: lap_temperature, inverse_reynolds, heat_coefficient

    if (.not. viscosity_is_enabled(nse)) return

    !$OMP MASKED
    call ensure_primitive_workspace(sim, js, je, ks, ke)
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          rho = max(q(i,j,k,1), nse%small_rho)
          u = q(i,j,k,2) / rho
          v = q(i,j,k,3) / rho
          w = q(i,j,k,4) / rho
          pressure = max((nse%gamma-1.0_dp) * &
            (q(i,j,k,5)-0.5_dp*rho*(u*u+v*v+w*w)), nse%small_p)
          primitive(i,j,k,1) = u
          primitive(i,j,k,2) = v
          primitive(i,j,k,3) = w
          primitive(i,j,k,4) = pressure / rho
        end do
      end do
    end do
    !$OMP END DO

    inverse_reynolds = 1.0_dp / nse%reynolds
    heat_coefficient = nse%gamma / &
      ((nse%gamma-1.0_dp) * nse%prandtl)

    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          u = primitive(i,j,k,1)
          v = primitive(i,j,k,2)
          w = primitive(i,j,k,3)

          ux = derivative_x(primitive(:,:,:,1), i, j, k, sim, js, ks)
          uy = derivative_y(primitive(:,:,:,1), i, j, k, sim, js, ks)
          uz = derivative_z(primitive(:,:,:,1), i, j, k, sim, js, ks)
          vx = derivative_x(primitive(:,:,:,2), i, j, k, sim, js, ks)
          vy = derivative_y(primitive(:,:,:,2), i, j, k, sim, js, ks)
          vz = derivative_z(primitive(:,:,:,2), i, j, k, sim, js, ks)
          wx = derivative_x(primitive(:,:,:,3), i, j, k, sim, js, ks)
          wy = derivative_y(primitive(:,:,:,3), i, j, k, sim, js, ks)
          wz = derivative_z(primitive(:,:,:,3), i, j, k, sim, js, ks)

          uxx = second_x(primitive(:,:,:,1), i, j, k, sim, js, ks)
          uyy = second_y(primitive(:,:,:,1), i, j, k, sim, js, ks)
          uzz = second_z(primitive(:,:,:,1), i, j, k, sim, js, ks)
          vxx = second_x(primitive(:,:,:,2), i, j, k, sim, js, ks)
          vyy = second_y(primitive(:,:,:,2), i, j, k, sim, js, ks)
          vzz = second_z(primitive(:,:,:,2), i, j, k, sim, js, ks)
          wxx = second_x(primitive(:,:,:,3), i, j, k, sim, js, ks)
          wyy = second_y(primitive(:,:,:,3), i, j, k, sim, js, ks)
          wzz = second_z(primitive(:,:,:,3), i, j, k, sim, js, ks)

          uxy = mixed_xy(primitive(:,:,:,1), i, j, k, sim, js, ks)
          uxz = mixed_xz(primitive(:,:,:,1), i, j, k, sim, js, ks)
          vxy = mixed_xy(primitive(:,:,:,2), i, j, k, sim, js, ks)
          vyz = mixed_yz(primitive(:,:,:,2), i, j, k, sim, js, ks)
          wxz = mixed_xz(primitive(:,:,:,3), i, j, k, sim, js, ks)
          wyz = mixed_yz(primitive(:,:,:,3), i, j, k, sim, js, ks)

          momentum_x = (4.0_dp/3.0_dp)*uxx + uyy + uzz + &
            (vxy+wxz)/3.0_dp
          momentum_y = vxx + (4.0_dp/3.0_dp)*vyy + vzz + &
            (uxy+wyz)/3.0_dp
          momentum_z = wxx + wyy + (4.0_dp/3.0_dp)*wzz + &
            (uxz+vyz)/3.0_dp

          div_velocity = ux + vy + wz
          tau_xx = 2.0_dp*ux - (2.0_dp/3.0_dp)*div_velocity
          tau_yy = 2.0_dp*vy - (2.0_dp/3.0_dp)*div_velocity
          tau_zz = 2.0_dp*wz - (2.0_dp/3.0_dp)*div_velocity
          tau_xy = uy + vx
          tau_xz = uz + wx
          tau_yz = vz + wy
          dissipation = tau_xx*ux + tau_yy*vy + tau_zz*wz + &
            tau_xy*(uy+vx) + tau_xz*(uz+wx) + tau_yz*(vz+wy)

          lap_temperature = second_x(primitive(:,:,:,4), i, j, k, &
            sim, js, ks) + second_y(primitive(:,:,:,4), i, j, k, &
            sim, js, ks) + second_z(primitive(:,:,:,4), i, j, k, &
            sim, js, ks)

          rhs(i,j,k,2) = rhs(i,j,k,2) + inverse_reynolds*momentum_x
          rhs(i,j,k,3) = rhs(i,j,k,3) + inverse_reynolds*momentum_y
          rhs(i,j,k,4) = rhs(i,j,k,4) + inverse_reynolds*momentum_z
          rhs(i,j,k,5) = rhs(i,j,k,5) + inverse_reynolds * &
            (u*momentum_x + v*momentum_y + w*momentum_z + &
             dissipation + heat_coefficient*lap_temperature)
        end do
      end do
    end do
    !$OMP END DO
  end subroutine add_viscous_rhs

  subroutine viscous_dt_limit(q, dt_limit, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: dt_limit
    integer :: i, j, k
    real(dp) :: minimum_density, maximum_diffusivity, inverse_spacing_square

    dt_limit = huge(1.0_dp)
    if (.not. viscosity_is_enabled(nse)) return

    minimum_density = huge(1.0_dp)
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          minimum_density = min(minimum_density, &
            max(q(i,j,k,1), nse%small_rho))
        end do
      end do
    end do

    maximum_diffusivity = max(4.0_dp/3.0_dp, &
      nse%gamma/nse%prandtl) / (nse%reynolds*minimum_density)
    inverse_spacing_square = 1.0_dp/(sim%dx*sim%dx) + &
      1.0_dp/(sim%dy*sim%dy) + 1.0_dp/(sim%dz*sim%dz)
    dt_limit = diffusion_stability_radius / &
      (d2_spectral_radius*maximum_diffusivity*inverse_spacing_square)
  end subroutine viscous_dt_limit

  subroutine validate_viscous_scheme(nse)
    type(nse_config), intent(in) :: nse
    character(len=32) :: requested

    requested = trim(adjustl(nse%viscous_scheme))
    if (requested /= 'none' .and. requested /= viscous_scheme_name()) then
      write(*,'(A,A,A)') 'ERROR: unsupported viscous scheme "', &
        trim(requested), '"; available schemes: none, central6'
      error stop
    end if
    if (requested == viscous_scheme_name()) then
      if (nse%reynolds <= 0.0_dp) then
        error stop 'central6 viscosity requires reynolds > 0'
      end if
      if (nse%prandtl <= 0.0_dp) then
        error stop 'central6 viscosity requires prandtl > 0'
      end if
      if (nse%gamma <= 1.0_dp) then
        error stop 'central6 viscosity requires gamma > 1'
      end if
    end if
  end subroutine validate_viscous_scheme

  integer function viscous_required_ghost_cells() result(nghost)
    nghost = 3
  end function viscous_required_ghost_cells

  pure function viscous_scheme_name() result(name)
    character(len=32) :: name
    name = 'central6'
  end function viscous_scheme_name

  pure logical function viscosity_is_enabled(nse) result(enabled)
    type(nse_config), intent(in) :: nse
    enabled = trim(adjustl(nse%viscous_scheme)) == viscous_scheme_name()
  end function viscosity_is_enabled

  subroutine ensure_primitive_workspace(sim, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: js, je, ks, ke
    logical :: correct_shape

    correct_shape = allocated(primitive)
    if (correct_shape) then
      correct_shape = lbound(primitive,1) == 1-sim%nghost .and. &
        ubound(primitive,1) == sim%nx+sim%nghost .and. &
        lbound(primitive,2) == js-sim%nghost .and. &
        ubound(primitive,2) == je+sim%nghost .and. &
        lbound(primitive,3) == ks-sim%nghost .and. &
        ubound(primitive,3) == ke+sim%nghost
    end if
    if (correct_shape) return

    if (allocated(primitive)) deallocate(primitive)
    allocate(primitive(1-sim%nghost:sim%nx+sim%nghost, &
      js-sim%nghost:je+sim%nghost, ks-sim%nghost:ke+sim%nghost, 4))
  end subroutine ensure_primitive_workspace

  pure real(dp) function derivative_x(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)

    value = (-f(i-3,j,k) + 9.0_dp*f(i-2,j,k) - 45.0_dp*f(i-1,j,k) + &
      45.0_dp*f(i+1,j,k) - 9.0_dp*f(i+2,j,k) + f(i+3,j,k)) / &
      (60.0_dp*sim%dx)
  end function derivative_x

  pure real(dp) function derivative_y(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)

    value = (-f(i,j-3,k) + 9.0_dp*f(i,j-2,k) - 45.0_dp*f(i,j-1,k) + &
      45.0_dp*f(i,j+1,k) - 9.0_dp*f(i,j+2,k) + f(i,j+3,k)) / &
      (60.0_dp*sim%dy)
  end function derivative_y

  pure real(dp) function derivative_z(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)

    value = (-f(i,j,k-3) + 9.0_dp*f(i,j,k-2) - 45.0_dp*f(i,j,k-1) + &
      45.0_dp*f(i,j,k+1) - 9.0_dp*f(i,j,k+2) + f(i,j,k+3)) / &
      (60.0_dp*sim%dz)
  end function derivative_z

  pure real(dp) function second_x(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)

    value = (2.0_dp*f(i-3,j,k) - 27.0_dp*f(i-2,j,k) + &
      270.0_dp*f(i-1,j,k) - 490.0_dp*f(i,j,k) + &
      270.0_dp*f(i+1,j,k) - 27.0_dp*f(i+2,j,k) + &
      2.0_dp*f(i+3,j,k)) / (180.0_dp*sim%dx*sim%dx)
  end function second_x

  pure real(dp) function second_y(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)

    value = (2.0_dp*f(i,j-3,k) - 27.0_dp*f(i,j-2,k) + &
      270.0_dp*f(i,j-1,k) - 490.0_dp*f(i,j,k) + &
      270.0_dp*f(i,j+1,k) - 27.0_dp*f(i,j+2,k) + &
      2.0_dp*f(i,j+3,k)) / (180.0_dp*sim%dy*sim%dy)
  end function second_y

  pure real(dp) function second_z(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)

    value = (2.0_dp*f(i,j,k-3) - 27.0_dp*f(i,j,k-2) + &
      270.0_dp*f(i,j,k-1) - 490.0_dp*f(i,j,k) + &
      270.0_dp*f(i,j,k+1) - 27.0_dp*f(i,j,k+2) + &
      2.0_dp*f(i,j,k+3)) / (180.0_dp*sim%dz*sim%dz)
  end function second_z

  pure real(dp) function mixed_xy(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)
    integer :: a, b

    value = 0.0_dp
    do b = 1, 3
      do a = 1, 3
        value = value + d1_positive(a)*d1_positive(b) * &
          (f(i+a,j+b,k)-f(i+a,j-b,k)-f(i-a,j+b,k)+f(i-a,j-b,k))
      end do
    end do
    value = value / (sim%dx*sim%dy)
  end function mixed_xy

  pure real(dp) function mixed_xz(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)
    integer :: a, b

    value = 0.0_dp
    do b = 1, 3
      do a = 1, 3
        value = value + d1_positive(a)*d1_positive(b) * &
          (f(i+a,j,k+b)-f(i+a,j,k-b)-f(i-a,j,k+b)+f(i-a,j,k-b))
      end do
    end do
    value = value / (sim%dx*sim%dz)
  end function mixed_xz

  pure real(dp) function mixed_yz(f, i, j, k, sim, js, ks) result(value)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: f(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:)
    integer :: a, b

    value = 0.0_dp
    do b = 1, 3
      do a = 1, 3
        value = value + d1_positive(a)*d1_positive(b) * &
          (f(i,j+a,k+b)-f(i,j+a,k-b)-f(i,j-a,k+b)+f(i,j-a,k-b))
      end do
    end do
    value = value / (sim%dy*sim%dz)
  end function mixed_yz

end module mod_viscous_scheme
