module mod_nse_forcing_common
  use mod_precision, only : dp
  use mod_model_config, only : nse_config
  implicit none
  private

  public :: forcing_is_enabled
  public :: validate_forcing_parameters
  public :: compute_forcing_coefficients

contains

  logical function forcing_is_enabled(nse) result(enabled)
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(nse%forcing_scheme)))
    case ('none')
      enabled = .false.
    case ('petersen_livescu')
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

    if (nse%forcing_max_coefficient > 0.0_dp) then
      coefficient_s = max(-nse%forcing_max_coefficient, &
        min(nse%forcing_max_coefficient, coefficient_s))
      coefficient_d = max(-nse%forcing_max_coefficient, &
        min(nse%forcing_max_coefficient, coefficient_d))
    end if
  end subroutine compute_forcing_coefficients

end module mod_nse_forcing_common
