module mod_init_hit_spectral
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  implicit none
  private

  public :: initialize_hit_spectral

contains

  subroutine initialize_hit_spectral(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)

    write(*,'(A)') 'ERROR: hit_spectral requires the 2decomp_fftw init backend.'
    write(*,'(A)') 'Configure with NSE_INIT_FFT_BACKEND=2decomp_fftw.'
    write(*,'(A,A)') 'Requested spectrum: ', trim(nse%hit_spectrum)
    q(1,js,ks,1) = q(1,js,ks,1)
    error stop
  end subroutine initialize_hit_spectral

end module mod_init_hit_spectral
