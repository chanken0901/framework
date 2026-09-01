module mod_nse_initial_conditions
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_init_taylor_green, only : initialize_taylor_green
  use mod_init_hit_spectral, only : initialize_hit_spectral
  use mod_init_imported_turbulence, only : initialize_imported_turbulence
  use mod_init_shock_turbulence, only : initialize_shock_turbulence
  use mod_init_shock_tube_turbulence, only : &
    initialize_shock_tube_turbulence
  implicit none
  private

  public :: initialize_nse_state

contains

  subroutine initialize_nse_state(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    character(len=:), allocatable :: initial_condition

    initial_condition = trim(adjustl(sim%initial_condition))
    select case (initial_condition)
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
