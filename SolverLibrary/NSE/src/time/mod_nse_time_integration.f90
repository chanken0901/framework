module mod_nse_time_integration
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_nse_spatial_operator, only : compute_nse_rhs
  use mod_viscous_scheme, only : viscous_dt_limit
  use module_mpi, only : mp_barrier, mp_allminr8
  implicit none
  private
  ! Shared only across the OpenMP team advancing one solver state.
  real(dp), save :: budget_status = 0.0_dp

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
    real(dp), intent(inout) :: dt
    integer :: i, j, k, retry, stage

    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
        do i = 1-sim%nghost, sim%nx+sim%nghost
          q0(i,j,k,:) = q(i,j,k,:)
        end do
      end do
    end do
    !$OMP END DO

    do retry = 0,20
      do stage = 1,3
        call check_step_budget(q,rhs,0.0_dp,sim,nse,js,je,ks,ke,.true.)
        if (budget_status /= 0.0_dp) exit
        call compute_nse_rhs(q,rhs,fface,sim,nse,js,je,ks,ke)
        call check_step_budget(q,rhs,dt,sim,nse,js,je,ks,ke,.false.)
        if (budget_status /= 0.0_dp) exit
        !$OMP DO collapse(2) schedule(static)
        do k = ks,ke
          do j = js,je
            do i = 1,sim%nx
              select case(stage)
              case(1)
                q(i,j,k,:) = q0(i,j,k,:) + dt*rhs(i,j,k,:)
              case(2)
                q(i,j,k,:) = 0.75_dp*q0(i,j,k,:) + 0.25_dp*(q(i,j,k,:)+dt*rhs(i,j,k,:))
              case(3)
                q(i,j,k,:) = q0(i,j,k,:)/3.0_dp + (2.0_dp/3.0_dp)*(q(i,j,k,:)+dt*rhs(i,j,k,:))
              end select
            end do
          end do
        end do
        !$OMP END DO
      end do
      if (budget_status == 0.0_dp) then
        call check_step_budget(q,rhs,0.0_dp,sim,nse,js,je,ks,ke,.true.)
        if (budget_status == 0.0_dp) exit
      end if
      !$OMP DO collapse(2) schedule(static)
      do k = ks-sim%nghost,ke+sim%nghost
        do j = js-sim%nghost,je+sim%nghost
          do i = 1-sim%nghost,sim%nx+sim%nghost
            q(i,j,k,:) = q0(i,j,k,:)
          end do
        end do
      end do
      !$OMP END DO
      !$OMP MASKED
      if (budget_status > 1.0_dp .or. sim%use_fixed_dt .or. retry == 20) &
        error stop 'NSE density/internal-energy budget failed; original state restored; review dt/forcing/resolution'
      dt = 0.5_dp*dt
      !$OMP END MASKED
      !$OMP BARRIER
    end do
    !$OMP MASKED
    if (sim%t+dt <= sim%t) error stop 'time step cannot advance time'
    if (retry > 0 .and. sim%rank == 0) write(*,'(A,I0,A,ES16.8)') &
      '# positivity_retry count=',retry,' accepted_dt=',dt
    !$OMP END MASKED
    !$OMP BARRIER
  end subroutine advance_nse_ssprk3

  subroutine check_step_budget(q,rhs,dt,sim,nse,js,je,ks,ke,state_only)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js,je,ks,ke
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(in) :: rhs(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(in) :: dt
    logical, intent(in) :: state_only
    real(dp) :: s(5),trial(5),p
    integer :: i,j,k
    !$OMP MASKED
    budget_status=0.0_dp
    !$OMP END MASKED
    !$OMP BARRIER
    !$OMP DO collapse(2) schedule(static) reduction(max:budget_status)
    do k=ks,ke
      do j=js,je
        do i=1,sim%nx
          s=q(i,j,k,1:5)
          if (.not. admissible(s,nse%small_rho,nse%small_p,nse%gamma)) then
            budget_status=3.0_dp
            cycle
          end if
          if (state_only) cycle
          if (.not. all(ieee_is_finite(rhs(i,j,k,1:5)))) then
            budget_status=3.0_dp
            cycle
          end if
          p=(nse%gamma-1)*(s(5)-0.5_dp*sum(s(2:4)**2)/s(1))
          trial=s+dt*rhs(i,j,k,1:5)
          if (.not. admissible(trial,max(nse%small_rho,0.1_dp*s(1)), &
              max(nse%small_p,0.1_dp*p),nse%gamma)) budget_status=max(budget_status,1.0_dp)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP MASKED
    budget_status=-budget_status
    call mp_allminr8(budget_status)
    budget_status=-budget_status
    !$OMP END MASKED
    !$OMP BARRIER
  end subroutine

  pure logical function admissible(s,rho_floor,p_floor,gamma) result(valid)
    real(dp), intent(in) :: s(5),rho_floor,p_floor,gamma
    real(dp) :: p
    valid=.false.
    if (.not. all(ieee_is_finite(s))) return
    if (s(1)<rho_floor) return
    p=(gamma-1)*(s(5)-0.5_dp*sum(s(2:4)**2)/s(1))
    valid=ieee_is_finite(p) .and. p>=p_floor
  end function

  subroutine validate_time_integrator(nse)
    type(nse_config), intent(in) :: nse

    if (trim(adjustl(nse%time_integrator)) /= 'ssprk3') then
      write(*,'(A,A)') 'ERROR: unsupported NSE time integrator: ', &
        trim(adjustl(nse%time_integrator))
      error stop
    end if
  end subroutine validate_time_integrator

end module mod_nse_time_integration
