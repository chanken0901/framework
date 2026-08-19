module mod_riemann_roe
  use mod_precision, only : dp
  use mod_model_config, only : nse_config
  implicit none
  private

  public :: rotate_conserved_to_normal
  public :: rotate_flux_to_global
  public :: roe_eigensystem
  public :: roe_numerical_flux
  public :: euler_physical_flux

contains

  pure subroutine rotate_conserved_to_normal(global_state, direction, &
      normal_state)
    real(dp), intent(in) :: global_state(5)
    integer, intent(in) :: direction
    real(dp), intent(out) :: normal_state(5)

    normal_state(1) = global_state(1)
    normal_state(5) = global_state(5)
    select case (direction)
    case (1)
      normal_state(2:4) = global_state(2:4)
    case (2)
      normal_state(2) = global_state(3)
      normal_state(3) = global_state(2)
      normal_state(4) = global_state(4)
    case (3)
      normal_state(2) = global_state(4)
      normal_state(3) = global_state(2)
      normal_state(4) = global_state(3)
    case default
      normal_state = 0.0_dp
    end select
  end subroutine rotate_conserved_to_normal

  pure subroutine rotate_flux_to_global(normal_flux, direction, global_flux)
    real(dp), intent(in) :: normal_flux(5)
    integer, intent(in) :: direction
    real(dp), intent(out) :: global_flux(5)

    global_flux(1) = normal_flux(1)
    global_flux(5) = normal_flux(5)
    select case (direction)
    case (1)
      global_flux(2:4) = normal_flux(2:4)
    case (2)
      global_flux(2) = normal_flux(3)
      global_flux(3) = normal_flux(2)
      global_flux(4) = normal_flux(4)
    case (3)
      global_flux(2) = normal_flux(3)
      global_flux(3) = normal_flux(4)
      global_flux(4) = normal_flux(2)
    case default
      global_flux = 0.0_dp
    end select
  end subroutine rotate_flux_to_global

  pure subroutine roe_numerical_flux(left_state, right_state, nse, flux)
    real(dp), intent(in) :: left_state(5), right_state(5)
    type(nse_config), intent(in) :: nse
    real(dp), intent(out) :: flux(5)
    real(dp) :: left_flux(5), right_flux(5), right_matrix(5,5)
    real(dp) :: left_matrix(5,5), eigenvalue(5), absolute_eigenvalue(5)
    real(dp) :: wave_strength(5), state_jump(5), dissipation(5)
    real(dp) :: rho_l, u_l, v_l, w_l, pressure_l, enthalpy_l
    real(dp) :: rho_r, u_r, v_r, w_r, pressure_r, enthalpy_r
    real(dp) :: sound_l, sound_r, entropy_delta(5)
    integer :: wave

    call primitive_state(left_state, nse, rho_l, u_l, v_l, w_l, &
      pressure_l, enthalpy_l)
    call primitive_state(right_state, nse, rho_r, u_r, v_r, w_r, &
      pressure_r, enthalpy_r)
    call euler_physical_flux(left_state, nse, left_flux)
    call euler_physical_flux(right_state, nse, right_flux)
    call roe_eigensystem(left_state, right_state, nse, right_matrix, &
      left_matrix, eigenvalue)

    sound_l = sqrt(nse%gamma*pressure_l/rho_l)
    sound_r = sqrt(nse%gamma*pressure_r/rho_r)
    entropy_delta(1) = max(0.0_dp, (u_r-sound_r)-(u_l-sound_l))
    entropy_delta(2:4) = max(0.0_dp, u_r-u_l)
    entropy_delta(5) = max(0.0_dp, (u_r+sound_r)-(u_l+sound_l))
    do wave = 1, 5
      absolute_eigenvalue(wave) = harten_hyman_fix( &
        eigenvalue(wave), entropy_delta(wave))
    end do

    state_jump = right_state - left_state
    wave_strength = matmul(left_matrix, state_jump)
    dissipation = matmul(right_matrix, absolute_eigenvalue*wave_strength)
    flux = 0.5_dp*(left_flux+right_flux-dissipation)
  end subroutine roe_numerical_flux

  pure subroutine euler_physical_flux(state, nse, flux)
    real(dp), intent(in) :: state(5)
    type(nse_config), intent(in) :: nse
    real(dp), intent(out) :: flux(5)
    real(dp) :: rho, u, v, w, pressure, enthalpy

    call primitive_state(state, nse, rho, u, v, w, pressure, enthalpy)
    flux(1) = rho*u
    flux(2) = rho*u*u + pressure
    flux(3) = rho*u*v
    flux(4) = rho*u*w
    flux(5) = u*(state(5)+pressure)
  end subroutine euler_physical_flux

  pure subroutine roe_eigensystem(left_state, right_state, nse, &
      right_matrix, left_matrix, eigenvalue)
    real(dp), intent(in) :: left_state(5), right_state(5)
    type(nse_config), intent(in) :: nse
    real(dp), intent(out) :: right_matrix(5,5), left_matrix(5,5)
    real(dp), intent(out) :: eigenvalue(5)
    real(dp) :: rho_l, u_l, v_l, w_l, pressure_l, enthalpy_l
    real(dp) :: rho_r, u_r, v_r, w_r, pressure_r, enthalpy_r
    real(dp) :: sqrt_l, sqrt_r, denominator, roe_density
    real(dp) :: u, v, w, enthalpy, sound_speed, sound_speed_squared
    real(dp) :: velocity_squared, gamma_minus_one, b1, b2

    call primitive_state(left_state, nse, rho_l, u_l, v_l, w_l, &
      pressure_l, enthalpy_l)
    call primitive_state(right_state, nse, rho_r, u_r, v_r, w_r, &
      pressure_r, enthalpy_r)

    sqrt_l = sqrt(rho_l)
    sqrt_r = sqrt(rho_r)
    denominator = sqrt_l + sqrt_r
    roe_density = max(sqrt_l*sqrt_r, nse%small_rho)
    u = (sqrt_l*u_l + sqrt_r*u_r) / denominator
    v = (sqrt_l*v_l + sqrt_r*v_r) / denominator
    w = (sqrt_l*w_l + sqrt_r*w_r) / denominator
    enthalpy = (sqrt_l*enthalpy_l + sqrt_r*enthalpy_r) / denominator
    velocity_squared = u*u + v*v + w*w
    sound_speed_squared = max((nse%gamma-1.0_dp) * &
      (enthalpy-0.5_dp*velocity_squared), nse%small_p/roe_density)
    sound_speed = sqrt(sound_speed_squared)

    eigenvalue = [u-sound_speed, u, u, u, u+sound_speed]
    right_matrix = 0.0_dp
    right_matrix(:,1) = [1.0_dp, u-sound_speed, v, w, &
      enthalpy-u*sound_speed]
    right_matrix(:,2) = [1.0_dp, u, v, w, 0.5_dp*velocity_squared]
    right_matrix(:,3) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, v]
    right_matrix(:,4) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, w]
    right_matrix(:,5) = [1.0_dp, u+sound_speed, v, w, &
      enthalpy+u*sound_speed]

    gamma_minus_one = nse%gamma - 1.0_dp
    b1 = gamma_minus_one*velocity_squared / (2.0_dp*sound_speed_squared)
    b2 = gamma_minus_one / sound_speed_squared
    left_matrix = 0.0_dp
    left_matrix(1,:) = 0.5_dp * [b1+u/sound_speed, &
      -b2*u-1.0_dp/sound_speed, -b2*v, -b2*w, b2]
    left_matrix(2,:) = [1.0_dp-b1, b2*u, b2*v, b2*w, -b2]
    left_matrix(3,:) = [-v, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    left_matrix(4,:) = [-w, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    left_matrix(5,:) = 0.5_dp * [b1-u/sound_speed, &
      -b2*u+1.0_dp/sound_speed, -b2*v, -b2*w, b2]
  end subroutine roe_eigensystem

  pure subroutine primitive_state(state, nse, rho, u, v, w, pressure, &
      enthalpy)
    real(dp), intent(in) :: state(5)
    type(nse_config), intent(in) :: nse
    real(dp), intent(out) :: rho, u, v, w, pressure, enthalpy

    rho = max(state(1), nse%small_rho)
    u = state(2) / rho
    v = state(3) / rho
    w = state(4) / rho
    pressure = max((nse%gamma-1.0_dp) * &
      (state(5)-0.5_dp*rho*(u*u+v*v+w*w)), nse%small_p)
    enthalpy = (state(5)+pressure) / rho
  end subroutine primitive_state

  pure real(dp) function harten_hyman_fix(eigenvalue, delta) result(value)
    real(dp), intent(in) :: eigenvalue, delta

    if (delta > 0.0_dp .and. abs(eigenvalue) < delta) then
      value = 0.5_dp*(eigenvalue*eigenvalue/delta+delta)
    else
      value = abs(eigenvalue)
    end if
  end function harten_hyman_fix

end module mod_riemann_roe
