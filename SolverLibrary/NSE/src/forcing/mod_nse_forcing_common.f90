module mod_nse_forcing_common
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_model_config, only : nse_config
  implicit none
  private

  public :: forcing_is_enabled
  public :: validate_forcing_parameters
  public :: compute_forcing_coefficients
  public :: compute_watanabe_coefficient

contains

  logical function forcing_is_enabled(nse) result(enabled)
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(nse%forcing_scheme)))
    case ('none')
      enabled = .false.
    case ('petersen_livescu', 'watanabe2021')
      enabled = .true.
    case default
      write(*,'(A,A)') 'ERROR: unsupported forcing_scheme: ', &
        trim(nse%forcing_scheme)
      error stop
    end select
  end function forcing_is_enabled

  subroutine validate_forcing_parameters(nse, compiled_backend)
    type(nse_config), intent(in) :: nse
    character(len=*), intent(in) :: compiled_backend
    character(len=32) :: requested_backend

    if (.not. forcing_is_enabled(nse)) return

    if (trim(nse%forcing_scheme) == 'watanabe2021') then
      if (nse%nv /= 5) error stop 'watanabe2021 requires single-component NSE'
      if (trim(compiled_backend) /= '2decomp_fftw' .and. trim(compiled_backend) /= 'cufft' .and. &
          trim(compiled_backend) /= 'cufftmp') error stop 'watanabe2021 requires a forcing FFT backend'
      if (trim(nse%forcing_spectrum) /= 'full_spectrum') &
        error stop 'watanabe2021 requires full_spectrum'
      if (.not. all(ieee_is_finite([nse%forcing_target_velocity_variance, &
          nse%forcing_feedback_exponent, nse%forcing_target_dissipation, nse%reynolds]))) &
        error stop 'non-finite watanabe2021 parameter'
      if (nse%forcing_target_velocity_variance <= 0.0_dp .or. &
          nse%forcing_feedback_exponent < 0.0_dp) error stop 'invalid watanabe2021 feedback parameter'
      if (nse%forcing_dilatational_ratio /= 0.0_dp .or. nse%forcing_max_coefficient /= 0.0_dp) &
        error stop 'watanabe2021 does not allow dilatational forcing or coefficient clipping'
      if (nse%reynolds <= 0.0_dp .or. trim(nse%viscous_scheme) /= 'central6' .or. nse%fh_enabled) &
        error stop 'watanabe2021 currently requires CENTRAL6 viscosity without fluctuating hydrodynamics'
    end if

    if (trim(adjustl(nse%boundary_condition)) /= 'periodic') then
      error stop 'Petersen-Livescu forcing requires periodic boundaries'
    end if
    select case (trim(adjustl(nse%forcing_spectrum)))
    case ('full_spectrum')
    case ('low_wavenumber')
      if (nse%forcing_k_cutoff <= 0.0_dp) then
        error stop 'low_wavenumber forcing requires forcing_k_cutoff > 0'
      end if
    case default
      error stop 'forcing_spectrum must be full_spectrum or low_wavenumber'
    end select
    if (nse%forcing_target_dissipation <= 0.0_dp) then
      error stop 'forcing_target_dissipation must be positive'
    end if
    if (nse%forcing_dilatational_ratio < 0.0_dp) then
      error stop 'forcing_dilatational_ratio must be non-negative'
    end if
    if (nse%forcing_denominator_floor <= 0.0_dp) then
      error stop 'forcing_denominator_floor must be positive'
    end if
    if (nse%forcing_max_coefficient < 0.0_dp) then
      error stop 'forcing_max_coefficient must be non-negative'
    end if
    if (nse%forcing_report_interval < 0) then
      error stop 'forcing_report_interval must be non-negative'
    end if

    requested_backend = trim(adjustl(nse%forcing_fft_backend))
    if (requested_backend /= 'auto' .and. &
        requested_backend /= trim(compiled_backend)) then
      write(*,'(A,A,A,A)') 'ERROR: forcing_fft_backend=', &
        trim(requested_backend), ' but compiled backend is ', &
        trim(compiled_backend)
      error stop
    end if
  end subroutine validate_forcing_parameters

  subroutine compute_forcing_coefficients(nse, solenoidal_denominator, &
      dilatational_denominator, pressure_dilatation, coefficient_s, &
      coefficient_d, target_s, target_d)
    type(nse_config), intent(in) :: nse
    real(dp), intent(in) :: solenoidal_denominator
    real(dp), intent(in) :: dilatational_denominator
    real(dp), intent(in) :: pressure_dilatation
    real(dp), intent(out) :: coefficient_s, coefficient_d
    real(dp), intent(out) :: target_s, target_d
    real(dp) :: ratio, dilatational_numerator

    if (.not. all(ieee_is_finite([solenoidal_denominator, &
        dilatational_denominator, pressure_dilatation]))) &
      error stop 'non-finite forcing denominator or pressure dilatation'
    ratio = nse%forcing_dilatational_ratio
    target_s = nse%forcing_target_dissipation / (1.0_dp + ratio)
    target_d = nse%forcing_target_dissipation - target_s
    dilatational_numerator = target_d - pressure_dilatation

    if (solenoidal_denominator <= nse%forcing_denominator_floor) then
      error stop 'Petersen-Livescu solenoidal forcing denominator is too small'
    end if
    coefficient_s = target_s / solenoidal_denominator

    if (dilatational_denominator <= nse%forcing_denominator_floor) then
      if (target_d <= nse%forcing_denominator_floor) then
        ! A purely solenoidal target cannot correct pressure dilatation when
        ! the selected band contains no dilatational energy.
        coefficient_d = 0.0_dp
      else
        error stop 'Petersen-Livescu dilatational forcing denominator is too small'
      end if
    else
      coefficient_d = dilatational_numerator / dilatational_denominator
    end if

    if (.not. all(ieee_is_finite([coefficient_s,coefficient_d]))) &
      error stop 'non-finite forcing coefficient'
    if (nse%forcing_max_coefficient > 0.0_dp) then
      coefficient_s = max(-nse%forcing_max_coefficient, &
        min(nse%forcing_max_coefficient, coefficient_s))
      coefficient_d = max(-nse%forcing_max_coefficient, &
        min(nse%forcing_max_coefficient, coefficient_d))
    end if
  end subroutine compute_forcing_coefficients

  subroutine compute_watanabe_coefficient(nse, denominator_s, pressure_dilatation, &
      variance, coefficient_s)
    type(nse_config), intent(in) :: nse
    real(dp), intent(in) :: denominator_s, pressure_dilatation, variance
    real(dp), intent(out) :: coefficient_s
    real(dp) :: alpha

    if (.not. all(ieee_is_finite([denominator_s, pressure_dilatation, variance]))) &
      error stop 'non-finite watanabe2021 statistics'
    if (denominator_s <= nse%forcing_denominator_floor .or. variance <= 0.0_dp) &
      error stop 'watanabe2021 requires nonzero solenoidal energy and velocity variance'
    ! Eq. (10),(12), in solver units. <u.grad(p)> = -<p.div(u)> on a periodic grid.
    alpha = (variance / nse%forcing_target_velocity_variance)**(-nse%forcing_feedback_exponent)
    coefficient_s = alpha * (nse%forcing_target_dissipation-pressure_dilatation) / denominator_s
    if (.not. ieee_is_finite(coefficient_s)) error stop 'non-finite watanabe2021 coefficient'
  end subroutine compute_watanabe_coefficient

end module mod_nse_forcing_common
