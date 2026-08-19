module mod_nse_time_integration
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_nse_spatial_operator, only : compute_nse_rhs
  use mod_viscous_scheme, only : viscous_dt_limit
  use module_mpi, only : mp_barrier, mp_allminr8
  implicit none
  private

  public :: compute_nse_dt
  public :: advance_nse_ssprk3
  public :: validate_time_integrator

contains

  subroutine compute_nse_dt(q, dt, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(out) :: dt
    integer :: i, j, k
    real(dp) :: rho, u, v, w, p, sound_speed, max_speed, diffusion_dt

    max_speed = 0.0_dp
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          rho = max(q(i,j,k,1), nse%small_rho)
          u = q(i,j,k,2) / rho
          v = q(i,j,k,3) / rho
          w = q(i,j,k,4) / rho
          p = max((nse%gamma-1.0_dp) * &
            (q(i,j,k,5)-0.5_dp*rho*(u*u+v*v+w*w)), nse%small_p)
          sound_speed = sqrt(nse%gamma * p / rho)
          max_speed = max(max_speed, abs(u)+sound_speed, &
            abs(v)+sound_speed, abs(w)+sound_speed)
        end do
      end do
    end do

    dt = nse%cfl * min(sim%dx/max_speed, sim%dy/max_speed, sim%dz/max_speed)
    call viscous_dt_limit(q, diffusion_dt, sim, nse, js, je, ks, ke)
    dt = min(dt, diffusion_dt)
    call mp_barrier
    call mp_allminr8(dt)
  end subroutine compute_nse_dt

  subroutine advance_nse_ssprk3(q, q0, rhs, fface, dt, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(inout) :: q0(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(inout) :: fface(0:, js-1:, ks-1:, :)
    real(dp), intent(in) :: dt
    integer :: i, j, k

    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          q0(i,j,k,:) = q(i,j,k,:)
        end do
      end do
    end do
    !$OMP END DO

    call compute_nse_rhs(q, rhs, fface, sim, nse, js, je, ks, ke)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          q(i,j,k,:) = q0(i,j,k,:) + dt * rhs(i,j,k,:)
        end do
      end do
    end do
    !$OMP END DO

    call compute_nse_rhs(q, rhs, fface, sim, nse, js, je, ks, ke)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          q(i,j,k,:) = 0.75_dp*q0(i,j,k,:) + &
            0.25_dp*(q(i,j,k,:) + dt*rhs(i,j,k,:))
        end do
      end do
    end do
    !$OMP END DO

    call compute_nse_rhs(q, rhs, fface, sim, nse, js, je, ks, ke)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          q(i,j,k,:) = (1.0_dp/3.0_dp)*q0(i,j,k,:) + &
            (2.0_dp/3.0_dp)*(q(i,j,k,:) + dt*rhs(i,j,k,:))
        end do
      end do
    end do
    !$OMP END DO
  end subroutine advance_nse_ssprk3

  subroutine validate_time_integrator(nse)
    type(nse_config), intent(in) :: nse

    if (trim(adjustl(nse%time_integrator)) /= 'ssprk3') then
      write(*,'(A,A)') 'ERROR: unsupported NSE time integrator: ', &
        trim(adjustl(nse%time_integrator))
      error stop
    end if
  end subroutine validate_time_integrator

end module mod_nse_time_integration
