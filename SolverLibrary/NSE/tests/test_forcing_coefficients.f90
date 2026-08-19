program test_forcing_coefficients
  use mod_precision, only : dp
  use mod_model_config, only : nse_config, init_nse_config
  use mod_nse_forcing_common, only : compute_forcing_coefficients
  implicit none

  type(nse_config) :: nse
  real(dp) :: coefficient_s, coefficient_d, target_s, target_d

  call init_nse_config(nse)
  nse%forcing_target_dissipation = 0.3_dp
  nse%forcing_dilatational_ratio = 0.5_dp
  nse%forcing_denominator_floor = 1.0e-14_dp
  call compute_forcing_coefficients(nse, 0.4_dp, 0.2_dp, 0.02_dp, &
    coefficient_s, coefficient_d, target_s, target_d)

  call assert_close(target_s, 0.2_dp, 'solenoidal target')
  call assert_close(target_d, 0.1_dp, 'dilatational target')
  call assert_close(coefficient_s, 0.5_dp, 'solenoidal coefficient')
  call assert_close(coefficient_d, 0.4_dp, 'dilatational coefficient')

  nse%forcing_max_coefficient = 0.25_dp
  call compute_forcing_coefficients(nse, 0.4_dp, 0.2_dp, 0.02_dp, &
    coefficient_s, coefficient_d, target_s, target_d)
  call assert_close(coefficient_s, 0.25_dp, 'clipped solenoidal coefficient')
  call assert_close(coefficient_d, 0.25_dp, 'clipped dilatational coefficient')
  write(*,'(A)') 'Petersen-Livescu forcing coefficient test passed'

contains

  subroutine assert_close(actual, expected, label)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label

    if (abs(actual-expected) > 1.0e-12_dp) then
      write(*,'(A,A,2(1X,ES16.8))') 'ERROR: ', trim(label), actual, expected
      error stop
    end if
  end subroutine assert_close

end program test_forcing_coefficients
