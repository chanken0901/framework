program test_hit_parameter_resolution
  use mod_precision, only : dp
  use mod_model_config, only : nse_config, init_nse_config, &
    resolve_nse_flow_parameters
  implicit none

  type(nse_config) :: cfg
  real(dp), parameter :: tolerance = 1.0e-12_dp
  real(dp) :: expected_re_l, expected_rms, expected_nu

  call init_nse_config(cfg)
  cfg%hit_spectrum = 'johnsen'
  cfg%hit_turbulent_mach = 0.3_dp
  cfg%hit_turbulent_reynolds = 40.0_dp
  cfg%hit_integral_length = 2.0_dp
  cfg%hit_johnsen_length_scale_ratio = 2.5_dp

  call resolve_nse_flow_parameters(cfg, 'hit_spectral')

  expected_re_l = 240.0_dp
  expected_rms = 0.3_dp / sqrt(3.0_dp)
  expected_nu = sqrt(1.5_dp) * expected_rms * 2.0_dp / expected_re_l
  call assert_close(cfg%hit_integral_reynolds, expected_re_l, 'Re_L')
  call assert_close(cfg%hit_rms_velocity, expected_rms, 'component RMS')
  call assert_close(cfg%hit_kinematic_viscosity, expected_nu, 'nu')
  call assert_close(cfg%reynolds, 1.0_dp/expected_nu, 'solver Reynolds')
  call assert_close(cfg%hit_taylor_microscale, &
    2.0_dp*sqrt(10.0_dp/expected_re_l), 'Taylor microscale')
  call assert_close(cfg%hit_kolmogorov_length, &
    2.0_dp*expected_re_l**(-0.75_dp), 'Kolmogorov length')
  call assert_close(cfg%hit_peak_wavenumber, 2.5_dp, &
    'Johnsen peak wavenumber')

  call init_nse_config(cfg)
  cfg%reynolds = 123.0_dp
  call resolve_nse_flow_parameters(cfg, 'taylor_green')
  call assert_close(cfg%reynolds, 123.0_dp, 'non-HIT Reynolds')

  write(*,'(A)') 'HIT parameter resolution test passed'

contains

  subroutine assert_close(actual, expected, label)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label

    if (abs(actual-expected) > tolerance*max(1.0_dp,abs(expected))) then
      write(*,'(A,A,2(1X,ES24.16))') 'Mismatch for ', trim(label), &
        actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_hit_parameter_resolution
