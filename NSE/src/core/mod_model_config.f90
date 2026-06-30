module mod_model_config
  use mod_precision, only : dp
  implicit none
  private

  public :: gpe_config, nse_config
  public :: init_gpe_config, init_nse_config
  public :: print_gpe_config, print_nse_config

  type :: gpe_config
    real(dp) :: g = 50.0_dp
    real(dp) :: sigma0 = 1.0_dp
    real(dp) :: wx = 1.0_dp
    real(dp) :: wy = 1.0_dp
    real(dp) :: wz = 1.0_dp
    real(dp) :: hbar = 1.0_dp
    real(dp) :: mass = 1.0_dp
  end type gpe_config

  type :: nse_config
    integer :: nv = 5
    integer :: nghost = 3
    real(dp) :: gamma = 1.4_dp
    real(dp) :: small_rho = 1.0e-12_dp
    real(dp) :: small_p   = 1.0e-12_dp
    real(dp) :: rho0 = 1.0_dp
    real(dp) :: mach = 0.5_dp
    real(dp) :: reynolds = 0.0_dp
    real(dp) :: prandtl = 0.72_dp
  end type nse_config

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
    write(u,'(A,I10)') 'nghost = ', cfg%nghost
    write(u,'(A,ES16.8)') 'gamma  = ', cfg%gamma
    write(u,'(A,ES16.8)') 'mach   = ', cfg%mach
    write(u,'(A,ES16.8)') 'reynolds   = ', cfg%reynolds
    write(u,'(A,ES16.8)') 'prandtl   = ', cfg%prandtl
  end subroutine print_nse_config

end module mod_model_config
