module mod_nse_forcing
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_nse_forcing_common, only : forcing_is_enabled
  implicit none
  private

  public :: initialize_nse_forcing
  public :: add_nse_forcing_rhs
  public :: finalize_nse_forcing
  public :: validate_nse_forcing

contains

  subroutine validate_nse_forcing(nse)
    type(nse_config), intent(in) :: nse

    if (forcing_is_enabled(nse)) then
      error stop 'Forcing is enabled, but this build has no forcing FFT backend'
    end if
  end subroutine validate_nse_forcing

  subroutine initialize_nse_forcing(sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke

    call validate_nse_forcing(nse)
  end subroutine initialize_nse_forcing

  subroutine add_nse_forcing_rhs(q, rhs, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
  end subroutine add_nse_forcing_rhs

  subroutine finalize_nse_forcing()
  end subroutine finalize_nse_forcing

end module mod_nse_forcing
