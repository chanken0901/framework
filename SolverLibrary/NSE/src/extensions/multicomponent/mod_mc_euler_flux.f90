module mod_mc_euler_flux
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_euler_field, only : validate_mc_euler_state
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_pressure, mc_sound_speed
  implicit none
  private

  public :: compute_mc_euler_physical_flux
  public :: compute_mc_euler_rusanov_flux
  public :: compute_mc_euler_rhs
  public :: compute_mc_euler_timestep
  public :: advance_mc_euler_ssprk3

contains

  subroutine compute_mc_euler_physical_flux( &
      state, layout, gamma, direction, flux)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    integer, intent(in) :: direction
    real(dp), intent(out) :: flux(:)
    real(dp) :: density, velocity(3), normal_velocity, pressure

    density = mc_mixture_density(state, layout)
    velocity = state(layout%momentum)/density
    normal_velocity = velocity(direction)
    pressure = mc_pressure(state,layout,gamma)
    flux = 0.0_dp
    flux(layout%first_species:layout%last_species) = &
      state(layout%first_species:layout%last_species)*normal_velocity
    flux(layout%momentum) = state(layout%momentum)*normal_velocity
    flux(layout%momentum(direction)) = &
      flux(layout%momentum(direction)) + pressure
    flux(layout%total_energy) = &
      (state(layout%total_energy)+pressure)*normal_velocity
  end subroutine compute_mc_euler_physical_flux

  subroutine compute_mc_euler_rusanov_flux( &
      left, right, layout, gamma, direction, flux)
    real(dp), intent(in) :: left(:), right(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    integer, intent(in) :: direction
    real(dp), intent(out) :: flux(:)
    real(dp) :: left_flux(size(left)), right_flux(size(right))
    real(dp) :: left_speed, right_speed, wave_speed

    call compute_mc_euler_physical_flux( &
      left,layout,gamma,direction,left_flux)
    call compute_mc_euler_physical_flux( &
      right,layout,gamma,direction,right_flux)
    left_speed = abs(left(layout%momentum(direction)) / &
      mc_mixture_density(left,layout)) + mc_sound_speed(left,layout,gamma)
    right_speed = abs(right(layout%momentum(direction)) / &
      mc_mixture_density(right,layout)) + mc_sound_speed(right,layout,gamma)
    wave_speed = max(left_speed,right_speed)
    flux = 0.5_dp*(left_flux+right_flux) - &
      0.5_dp*wave_speed*(right-left)
  end subroutine compute_mc_euler_rusanov_flux

  subroutine compute_mc_euler_rhs(q, rhs, layout, config)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k, im, ip, jm, jp, km, kp
    real(dp) :: dx, dy, dz
    real(dp) :: positive_flux(layout%nvariables)
    real(dp) :: negative_flux(layout%nvariables)

    if (any(shape(rhs) /= shape(q))) then
      error stop 'multicomponent Euler RHS allocation does not match state'
    end if
    call validate_mc_euler_state(q,layout,config)
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    rhs = 0.0_dp
    do k = 1, config%nz
      km = merge(config%nz,k-1,k == 1)
      kp = merge(1,k+1,k == config%nz)
      do j = 1, config%ny
        jm = merge(config%ny,j-1,j == 1)
        jp = merge(1,j+1,j == config%ny)
        do i = 1, config%nx
          im = merge(config%nx,i-1,i == 1)
          ip = merge(1,i+1,i == config%nx)
          call compute_mc_euler_rusanov_flux( &
            q(i,j,k,:),q(ip,j,k,:),layout,config%gamma,1,positive_flux)
          call compute_mc_euler_rusanov_flux( &
            q(im,j,k,:),q(i,j,k,:),layout,config%gamma,1,negative_flux)
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (positive_flux-negative_flux)/dx
          call compute_mc_euler_rusanov_flux( &
            q(i,j,k,:),q(i,jp,k,:),layout,config%gamma,2,positive_flux)
          call compute_mc_euler_rusanov_flux( &
            q(i,jm,k,:),q(i,j,k,:),layout,config%gamma,2,negative_flux)
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (positive_flux-negative_flux)/dy
          call compute_mc_euler_rusanov_flux( &
            q(i,j,k,:),q(i,j,kp,:),layout,config%gamma,3,positive_flux)
          call compute_mc_euler_rusanov_flux( &
            q(i,j,km,:),q(i,j,k,:),layout,config%gamma,3,negative_flux)
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (positive_flux-negative_flux)/dz
        end do
      end do
    end do
  end subroutine compute_mc_euler_rhs

  real(dp) function compute_mc_euler_timestep(q,layout,config) result(dt)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k
    real(dp) :: density, velocity(3), sound_speed, rate, maximum_rate
    real(dp) :: dx, dy, dz

    call validate_mc_euler_state(q,layout,config)
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    maximum_rate = 0.0_dp
    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          density = mc_mixture_density(q(i,j,k,:),layout)
          velocity = q(i,j,k,layout%momentum)/density
          sound_speed = mc_sound_speed(q(i,j,k,:),layout,config%gamma)
          rate = (abs(velocity(1))+sound_speed)/dx + &
            (abs(velocity(2))+sound_speed)/dy + &
            (abs(velocity(3))+sound_speed)/dz
          maximum_rate = max(maximum_rate,rate)
        end do
      end do
    end do
    if (config%dt > 0.0_dp) then
      if (config%dt*maximum_rate > 1.0_dp+100.0_dp*epsilon(1.0_dp)) then
        error stop 'fixed multicomponent Euler dt violates the CFL limit'
      end if
      dt = config%dt
    else
      dt = config%cfl/maximum_rate
    end if
  end function compute_mc_euler_timestep

  subroutine advance_mc_euler_ssprk3(q,q0,rhs,dt,layout,config)
    real(dp), intent(inout) :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config

    if (any(shape(q0) /= shape(q)) .or. any(shape(rhs) /= shape(q))) then
      error stop 'multicomponent Euler SSPRK3 work arrays do not match state'
    end if
    q0 = q
    call compute_mc_euler_rhs(q,rhs,layout,config)
    q = q0 + dt*rhs
    call compute_mc_euler_rhs(q,rhs,layout,config)
    q = 0.75_dp*q0 + 0.25_dp*(q+dt*rhs)
    call compute_mc_euler_rhs(q,rhs,layout,config)
    q = (1.0_dp/3.0_dp)*q0 + (2.0_dp/3.0_dp)*(q+dt*rhs)
    call validate_mc_euler_state(q,layout,config)
  end subroutine advance_mc_euler_ssprk3

end module mod_mc_euler_flux
