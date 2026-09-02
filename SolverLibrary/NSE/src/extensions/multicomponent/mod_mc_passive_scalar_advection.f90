module mod_mc_passive_scalar_advection
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_passive_scalar_config, only : mc_passive_scalar_config
  implicit none
  private

  public :: compute_mc_passive_scalar_rhs
  public :: advance_mc_passive_scalar_ssprk3

contains

  subroutine compute_mc_passive_scalar_rhs(q, rhs, layout, config)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_passive_scalar_config), intent(in) :: config
    integer :: i, j, k, variable
    integer :: im, ip, jm, jp, km, kp
    real(dp) :: dx, dy, dz
    real(dp) :: positive_velocity(3), negative_velocity(3)

    if (any(shape(rhs) /= shape(q))) then
      error stop 'passive-scalar RHS allocation does not match state'
    end if
    dx = (config%x_max-config%x_min) / real(config%nx, dp)
    dy = (config%y_max-config%y_min) / real(config%ny, dp)
    dz = (config%z_max-config%z_min) / real(config%nz, dp)
    positive_velocity = max(config%velocity, 0.0_dp)
    negative_velocity = min(config%velocity, 0.0_dp)
    rhs = 0.0_dp

    do variable = layout%first_species, layout%last_species
      do k = 1, config%nz
        km = merge(config%nz, k-1, k == 1)
        kp = merge(1, k+1, k == config%nz)
        do j = 1, config%ny
          jm = merge(config%ny, j-1, j == 1)
          jp = merge(1, j+1, j == config%ny)
          do i = 1, config%nx
            im = merge(config%nx, i-1, i == 1)
            ip = merge(1, i+1, i == config%nx)
            rhs(i,j,k,variable) = -( &
              positive_velocity(1)*(q(i,j,k,variable)-q(im,j,k,variable))/dx + &
              negative_velocity(1)*(q(ip,j,k,variable)-q(i,j,k,variable))/dx + &
              positive_velocity(2)*(q(i,j,k,variable)-q(i,jm,k,variable))/dy + &
              negative_velocity(2)*(q(i,jp,k,variable)-q(i,j,k,variable))/dy + &
              positive_velocity(3)*(q(i,j,k,variable)-q(i,j,km,variable))/dz + &
              negative_velocity(3)*(q(i,j,kp,variable)-q(i,j,k,variable))/dz)
          end do
        end do
      end do
    end do
  end subroutine compute_mc_passive_scalar_rhs

  subroutine advance_mc_passive_scalar_ssprk3(q, q0, rhs, dt, layout, config)
    real(dp), intent(inout) :: q(:,:,:,:)
    real(dp), intent(inout) :: q0(:,:,:,:)
    real(dp), intent(inout) :: rhs(:,:,:,:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout
    type(mc_passive_scalar_config), intent(in) :: config
    integer :: first, last

    if (any(shape(q0) /= shape(q)) .or. any(shape(rhs) /= shape(q))) then
      error stop 'passive-scalar SSPRK3 work arrays do not match state'
    end if
    first = layout%first_species
    last = layout%last_species
    q0 = q

    call compute_mc_passive_scalar_rhs(q, rhs, layout, config)
    q(:,:,:,first:last) = q0(:,:,:,first:last) + &
      dt*rhs(:,:,:,first:last)

    call compute_mc_passive_scalar_rhs(q, rhs, layout, config)
    q(:,:,:,first:last) = 0.75_dp*q0(:,:,:,first:last) + &
      0.25_dp*(q(:,:,:,first:last) + dt*rhs(:,:,:,first:last))

    call compute_mc_passive_scalar_rhs(q, rhs, layout, config)
    q(:,:,:,first:last) = (1.0_dp/3.0_dp)*q0(:,:,:,first:last) + &
      (2.0_dp/3.0_dp)*(q(:,:,:,first:last) + &
      dt*rhs(:,:,:,first:last))
  end subroutine advance_mc_passive_scalar_ssprk3

end module mod_mc_passive_scalar_advection
