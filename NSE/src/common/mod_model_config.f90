module mod_model_config
  use mod_precision, only : dp
  implicit none
  private

  public :: gpe_config, nse_config
  public :: init_gpe_config, init_nse_config
  public :: print_gpe_config, print_nse_config

  type :: nse_config
    integer :: nv = 5
    real(dp) :: gamma = 1.4_dp
    real(dp) :: cfl = 0.5_dp
    real(dp) :: small_rho = 1.0e-12_dp
    real(dp) :: small_p   = 1.0e-12_dp
    real(dp) :: rho0 = 1.0_dp
    real(dp) :: mach = 0.5_dp
    real(dp) :: reynolds = 0.0_dp
    real(dp) :: prandtl = 0.72_dp
    character(len=32) :: convective_scheme = 'keep6'
    character(len=32) :: viscous_scheme = 'none'
    character(len=32) :: boundary_condition = 'periodic'
    character(len=32) :: time_integrator = 'ssprk3'
    character(len=32) :: hit_spectrum = 'johnsen'
    integer :: hit_seed = 13579
    real(dp) :: hit_rms_velocity = -1.0_dp
    real(dp) :: hit_peak_wavenumber = 4.0_dp
    real(dp) :: hit_integral_length = 1.0_dp
    real(dp) :: hit_kolmogorov_length = 0.02_dp
    real(dp) :: hit_dealias_fraction = 2.0_dp / 3.0_dp
    character(len=32) :: forcing_scheme = 'none'
    character(len=32) :: forcing_spectrum = 'low_wavenumber'
    character(len=32) :: forcing_fft_backend = 'auto'
    real(dp) :: forcing_k_cutoff = 2.5_dp
    real(dp) :: forcing_target_dissipation = 0.0_dp
    real(dp) :: forcing_dilatational_ratio = 0.0_dp
    real(dp) :: forcing_denominator_floor = 1.0e-14_dp
    real(dp) :: forcing_max_coefficient = 0.0_dp
    integer :: forcing_report_interval = 100
  end type nse_config

  type :: gpe_config
    real(dp) :: g = 50.0_dp
    real(dp) :: sigma0 = 1.0_dp
    real(dp) :: wx = 1.0_dp
    real(dp) :: wy = 1.0_dp
    real(dp) :: wz = 1.0_dp
    real(dp) :: hbar = 1.0_dp
    real(dp) :: mass = 1.0_dp
  end type gpe_config

contains

  subroutine init_gpe_config(cfg)
    type(gpe_config), intent(out) :: cfg
    cfg = gpe_config()
  end subroutine init_gpe_config

  subroutine init_nse_config(cfg)
    type(nse_config), intent(out) :: cfg
    cfg = nse_config()
  end subroutine init_nse_config

  subroutine print_gpe_config(cfg, unit)
    type(gpe_config), intent(in) :: cfg
    integer, intent(in), optional :: unit
    integer :: u
    u = 6
    if (present(unit)) u = unit
    write(u,'(A)') '--- gpe_config ---'
    write(u,'(A,ES16.8)') 'g      = ', cfg%g
    write(u,'(A,ES16.8)') 'sigma0 = ', cfg%sigma0
    write(u,'(A,3ES16.8)') 'wx, wy, wz = ', cfg%wx, cfg%wy, cfg%wz
  end subroutine print_gpe_config

  subroutine print_nse_config(cfg, unit)
    type(nse_config), intent(in) :: cfg
    integer, intent(in), optional :: unit
    integer :: u
    u = 6
    if (present(unit)) u = unit
    write(u,'(A)') '--- nse_config ---'
    write(u,'(A,I10)') 'nv     = ', cfg%nv
    write(u,'(A,ES16.8)') 'gamma  = ', cfg%gamma
    write(u,'(A,ES16.8)') 'cfl    = ', cfg%cfl
    write(u,'(A,ES16.8)') 'rho0   = ', cfg%rho0
    write(u,'(A,ES16.8)') 'mach   = ', cfg%mach
    write(u,'(A,ES16.8)') 'reynolds   = ', cfg%reynolds
    write(u,'(A,ES16.8)') 'prandtl   = ', cfg%prandtl
    write(u,'(A,A)') 'convective_scheme = ', trim(cfg%convective_scheme)
    write(u,'(A,A)') 'viscous_scheme    = ', trim(cfg%viscous_scheme)
    write(u,'(A,A)') 'boundary_condition = ', trim(cfg%boundary_condition)
    write(u,'(A,A)') 'time_integrator   = ', trim(cfg%time_integrator)
    write(u,'(A,A)') 'hit_spectrum      = ', trim(cfg%hit_spectrum)
    write(u,'(A,I10)') 'hit_seed          = ', cfg%hit_seed
    write(u,'(A,ES16.8)') 'hit_rms_velocity  = ', cfg%hit_rms_velocity
    write(u,'(A,ES16.8)') 'hit_peak_wavenumber = ', cfg%hit_peak_wavenumber
    write(u,'(A,ES16.8)') 'hit_integral_length = ', cfg%hit_integral_length
    write(u,'(A,ES16.8)') 'hit_kolmogorov_length = ', &
      cfg%hit_kolmogorov_length
    write(u,'(A,ES16.8)') 'hit_dealias_fraction = ', &
      cfg%hit_dealias_fraction
    write(u,'(A,A)') 'forcing_scheme    = ', trim(cfg%forcing_scheme)
    write(u,'(A,A)') 'forcing_spectrum  = ', trim(cfg%forcing_spectrum)
    write(u,'(A,A)') 'forcing_fft_backend = ', trim(cfg%forcing_fft_backend)
    write(u,'(A,ES16.8)') 'forcing_k_cutoff = ', cfg%forcing_k_cutoff
    write(u,'(A,ES16.8)') 'forcing_target_dissipation = ', &
      cfg%forcing_target_dissipation
    write(u,'(A,ES16.8)') 'forcing_dilatational_ratio = ', &
      cfg%forcing_dilatational_ratio
    write(u,'(A,ES16.8)') 'forcing_denominator_floor = ', &
      cfg%forcing_denominator_floor
    write(u,'(A,ES16.8)') 'forcing_max_coefficient = ', &
      cfg%forcing_max_coefficient
    write(u,'(A,I10)') 'forcing_report_interval = ', &
      cfg%forcing_report_interval
  end subroutine print_nse_config

end module mod_model_config
