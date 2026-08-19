module mod_convective_hybrid
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_convective_leaf_registry, only : compute_leaf_face_flux, &
    validate_leaf_scheme, leaf_required_ghost_cells
  implicit none
  private

  public :: compute_hybrid_flux
  public :: validate_hybrid_scheme
  public :: hybrid_required_ghost_cells
  public :: hybrid_face_weight

contains

  subroutine compute_hybrid_flux(q, fface, direction, sim, nse, &
      js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: direction, js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: fface(0:, js-1:, ks-1:, :)
    real(dp) :: smooth_flux(5), shock_flux(5), alpha
    integer :: i, j, k

    if (direction < 1 .or. direction > 3) then
      error stop 'convective flux direction must be 1, 2, or 3'
    end if
    if (sim%nghost < hybrid_required_ghost_cells(nse)) then
      error stop 'hybrid convective flux requires three ghost cells'
    end if

    !$OMP DO collapse(2) schedule(static)
    do k = ks-1, ke
      do j = js-1, je
        do i = 0, sim%nx
          alpha = hybrid_face_weight(q, i, j, k, direction, sim, nse, &
            js, ks)
          if (alpha <= 0.0_dp) then
            call compute_leaf_face_flux(nse%hybrid_smooth_scheme, q, &
              i, j, k, direction, sim, nse, js, ks, &
              fface(i,j,k,1:5))
          else if (alpha >= 1.0_dp) then
            call compute_leaf_face_flux(nse%hybrid_shock_scheme, q, &
              i, j, k, direction, sim, nse, js, ks, &
              fface(i,j,k,1:5))
          else
            call compute_leaf_face_flux(nse%hybrid_smooth_scheme, q, &
              i, j, k, direction, sim, nse, js, ks, smooth_flux)
            call compute_leaf_face_flux(nse%hybrid_shock_scheme, q, &
              i, j, k, direction, sim, nse, js, ks, shock_flux)
            fface(i,j,k,1:5) = (1.0_dp-alpha)*smooth_flux + &
              alpha*shock_flux
          end if
        end do
      end do
    end do
    !$OMP END DO
  end subroutine compute_hybrid_flux

  pure real(dp) function hybrid_face_weight(q, i, j, k, direction, &
      sim, nse, js, ks) result(alpha)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, direction, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    integer :: ip, jp, kp
    real(dp) :: raw_sensor, scaled

    ip = i
    jp = j
    kp = k
    select case (direction)
    case (1)
      ip = i+1
    case (2)
      jp = j+1
    case (3)
      kp = k+1
    case default
      alpha = 0.0_dp
      return
    end select

    select case (trim(adjustl(nse%hybrid_sensor)))
    case ('ducros_pressure')
      raw_sensor = max(ducros_pressure_cell_sensor(q, i, j, k, &
        direction, sim, nse, js, ks), &
        ducros_pressure_cell_sensor(q, ip, jp, kp, direction, sim, nse, &
        js, ks))
    case default
      raw_sensor = 0.0_dp
    end select

    scaled = (raw_sensor-nse%hybrid_sensor_onset) / &
      (nse%hybrid_sensor_full-nse%hybrid_sensor_onset)
    scaled = max(0.0_dp, min(1.0_dp, scaled))
    alpha = scaled*scaled*(3.0_dp-2.0_dp*scaled)
  end function hybrid_face_weight

  pure real(dp) function ducros_pressure_cell_sensor(q, i, j, k, &
      direction, sim, nse, js, ks) result(sensor)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, direction, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), parameter :: sensor_epsilon = 1.0e-30_dp
    real(dp) :: velocity_xm(3), velocity_xp(3)
    real(dp) :: velocity_ym(3), velocity_yp(3)
    real(dp) :: velocity_zm(3), velocity_zp(3)
    real(dp) :: pressure_minus, pressure_center, pressure_plus
    real(dp) :: divergence, scaled_divergence, compression, vorticity(3)
    real(dp) :: scaled_vorticity(3), rate_scale
    real(dp) :: vorticity_squared, pressure_curvature, ducros_factor

    call velocity_pressure_at(q, i-1, j, k, sim, nse, js, ks, &
      velocity_xm, pressure_minus)
    call velocity_pressure_at(q, i+1, j, k, sim, nse, js, ks, &
      velocity_xp, pressure_plus)
    call velocity_pressure_at(q, i, j-1, k, sim, nse, js, ks, &
      velocity_ym, pressure_center)
    call velocity_pressure_at(q, i, j+1, k, sim, nse, js, ks, &
      velocity_yp, pressure_center)
    call velocity_pressure_at(q, i, j, k-1, sim, nse, js, ks, &
      velocity_zm, pressure_center)
    call velocity_pressure_at(q, i, j, k+1, sim, nse, js, ks, &
      velocity_zp, pressure_center)

    divergence = (velocity_xp(1)-velocity_xm(1))/(2.0_dp*sim%dx) + &
      (velocity_yp(2)-velocity_ym(2))/(2.0_dp*sim%dy) + &
      (velocity_zp(3)-velocity_zm(3))/(2.0_dp*sim%dz)
    vorticity(1) = (velocity_yp(3)-velocity_ym(3))/(2.0_dp*sim%dy) - &
      (velocity_zp(2)-velocity_zm(2))/(2.0_dp*sim%dz)
    vorticity(2) = (velocity_zp(1)-velocity_zm(1))/(2.0_dp*sim%dz) - &
      (velocity_xp(3)-velocity_xm(3))/(2.0_dp*sim%dx)
    vorticity(3) = (velocity_xp(2)-velocity_xm(2))/(2.0_dp*sim%dx) - &
      (velocity_yp(1)-velocity_ym(1))/(2.0_dp*sim%dy)
    rate_scale = max(abs(divergence), maxval(abs(vorticity)))
    if (rate_scale <= sqrt(sensor_epsilon)) then
      ducros_factor = 0.0_dp
    else
      scaled_divergence = divergence/rate_scale
      scaled_vorticity = vorticity/rate_scale
      vorticity_squared = dot_product(scaled_vorticity, scaled_vorticity)
      compression = min(scaled_divergence, 0.0_dp)
      ducros_factor = compression*compression / &
        (scaled_divergence*scaled_divergence+vorticity_squared + &
        sensor_epsilon)
    end if

    call directional_pressures(q, i, j, k, direction, sim, nse, &
      js, ks, pressure_minus, pressure_center, pressure_plus)
    pressure_curvature = abs(pressure_plus-2.0_dp*pressure_center + &
      pressure_minus) / (pressure_plus+2.0_dp*pressure_center + &
      pressure_minus+sensor_epsilon)
    sensor = max(0.0_dp, min(1.0_dp, pressure_curvature*ducros_factor))
  end function ducros_pressure_cell_sensor

  pure subroutine directional_pressures(q, i, j, k, direction, sim, nse, &
      js, ks, pressure_minus, pressure_center, pressure_plus)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, direction, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: pressure_minus, pressure_center, pressure_plus
    real(dp) :: unused_velocity(3)
    integer :: di, dj, dk

    di = 0
    dj = 0
    dk = 0
    if (direction == 1) di = 1
    if (direction == 2) dj = 1
    if (direction == 3) dk = 1
    call velocity_pressure_at(q, i-di, j-dj, k-dk, sim, nse, js, ks, &
      unused_velocity, pressure_minus)
    call velocity_pressure_at(q, i, j, k, sim, nse, js, ks, &
      unused_velocity, pressure_center)
    call velocity_pressure_at(q, i+di, j+dj, k+dk, sim, nse, js, ks, &
      unused_velocity, pressure_plus)
  end subroutine directional_pressures

  pure subroutine velocity_pressure_at(q, i, j, k, sim, nse, js, ks, &
      velocity, pressure)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: velocity(3), pressure
    real(dp) :: density

    density = max(q(i,j,k,1), nse%small_rho)
    velocity = q(i,j,k,2:4)/density
    pressure = max((nse%gamma-1.0_dp)*(q(i,j,k,5) - &
      0.5_dp*density*dot_product(velocity, velocity)), nse%small_p)
  end subroutine velocity_pressure_at

  subroutine validate_hybrid_scheme(nse)
    type(nse_config), intent(in) :: nse

    call validate_leaf_scheme(nse%hybrid_smooth_scheme, nse, 'smooth')
    call validate_leaf_scheme(nse%hybrid_shock_scheme, nse, 'shock')
    if (trim(adjustl(nse%hybrid_sensor)) /= 'ducros_pressure') then
      write(*,'(A,A,A)') 'ERROR: unsupported hybrid sensor "', &
        trim(adjustl(nse%hybrid_sensor)), '"; use ducros_pressure'
      error stop
    end if
    if (nse%hybrid_sensor_onset < 0.0_dp .or. &
        nse%hybrid_sensor_full <= nse%hybrid_sensor_onset) then
      error stop 'hybrid sensor requires 0 <= onset < full'
    end if
  end subroutine validate_hybrid_scheme

  integer function hybrid_required_ghost_cells(nse) result(nghost)
    type(nse_config), intent(in) :: nse

    nghost = max(3, leaf_required_ghost_cells(nse%hybrid_smooth_scheme), &
      leaf_required_ghost_cells(nse%hybrid_shock_scheme))
  end function hybrid_required_ghost_cells

end module mod_convective_hybrid
