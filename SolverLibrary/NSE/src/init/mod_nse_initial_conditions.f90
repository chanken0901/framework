module mod_nse_initial_conditions
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_init_taylor_green, only : initialize_taylor_green
  use mod_init_hit_spectral, only : initialize_hit_spectral
  use mod_init_imported_turbulence, only : initialize_imported_turbulence
  use mod_init_imported_turbulence, only : initialize_restart
  use mod_init_shock_turbulence, only : initialize_shock_turbulence
  use mod_init_shock_tube_turbulence, only : &
    initialize_shock_tube_turbulence
  implicit none
  private

  public :: initialize_nse_state

contains

  subroutine initialize_nse_state(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(inout) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    character(len=:), allocatable :: initial_condition
    real(dp) :: state(5), rho, pressure, velocity(3)
    integer :: variable

    if (len_trim(sim%restart_file) > 0) then
      call initialize_restart(q, sim, nse, js, je, ks, ke)
      return
    end if
    sim%t = 0.0_dp
    sim%step = 0
    initial_condition = trim(adjustl(sim%initial_condition))
    select case (initial_condition)
    case ('uniform_flow')
      if (nse%nv/=5) error stop 'Uniform flow requires five conserved variables'
      rho = nse%uniform_state(1)
      velocity = nse%uniform_state(2:4)
      pressure = nse%uniform_state(5)
      if (.not.all(ieee_is_finite(nse%uniform_state))) error stop 'Nonfinite uniform initial state'
      if (rho<=nse%small_rho.or.pressure<=nse%small_p.or.nse%gamma<=1.0_dp) &
        error stop 'Invalid uniform initial density, pressure or gamma'
      state(1) = rho
      state(2:4) = rho*velocity
      state(5) = pressure/(nse%gamma-1.0_dp)+0.5_dp*rho*sum(velocity**2)
      if (.not.all(ieee_is_finite(state))) error stop 'Uniform conserved state overflow'
      do variable=1,5
        q(1:sim%nx,js:je,ks:ke,variable) = state(variable)
      end do
    case ('default', 'taylor_green', 'taylor_green_vortex')
      call initialize_taylor_green(q, sim, nse, js, je, ks, ke)
    case ('hit_spectral', 'homogeneous_isotropic_turbulence')
      call initialize_hit_spectral(q, sim, nse, js, je, ks, ke)
    case ('imported_turbulence')
      call initialize_imported_turbulence(q, sim, nse, js, je, ks, ke)
    case ('shock_turbulence_interaction')
      call initialize_shock_turbulence(q, sim, nse, js, je, ks, ke)
    case ('shock_tube_turbulence_interaction')
      call initialize_shock_tube_turbulence(q, sim, nse, js, je, ks, ke)
    case default
      write(*,'(A,A)') 'ERROR: unsupported NSE initial condition: ', initial_condition
      error stop
    end select
  end subroutine initialize_nse_state

end module mod_nse_initial_conditions
