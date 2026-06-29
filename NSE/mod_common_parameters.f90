module mod_common_parameters
  use mod_precision, only: dp
  implicit none
  private

  public :: case_id, case_label, physics_model
  public :: nx, ny, nz
  public :: x_min, x_max, x_center
  public :: y_min, y_max, y_center
  public :: z_min, z_max, z_center
  public :: t_max, dt
  public :: output_frequency
  public :: output_location, restart_location
  public :: pi
  public :: update_common_derived_parameters

  !==========================
  ! case information
  !==========================
  character(len=64)  :: case_id       = ""
  character(len=64)  :: case_label    = ""
  character(len=64)  :: physics_model = ""

  !==========================
  ! grid
  !==========================
  integer :: nx = 64
  integer :: ny = 64
  integer :: nz = 64

  real(dp) :: x_min = 0.0_dp
  real(dp) :: x_max = 2.0_dp * acos(-1.0_dp)
  real(dp) :: x_center = 0.0_dp

  real(dp) :: y_min = 0.0_dp
  real(dp) :: y_max = 2.0_dp * acos(-1.0_dp)
  real(dp) :: y_center = 0.0_dp

  real(dp) :: z_min = 0.0_dp
  real(dp) :: z_max = 2.0_dp * acos(-1.0_dp)
  real(dp) :: z_center = 0.0_dp

  !==========================
  ! time / output
  !==========================
  real(dp) :: dt    = 0.0_dp
  real(dp) :: t_max = 100.0_dp

  integer :: output_frequency = 10

  character(len=256) :: output_location  = "output"
  character(len=256) :: restart_location = "restart"

  !==========================
  ! constants
  !==========================
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  subroutine update_common_derived_parameters()
    x_center = 0.5_dp * (x_min + x_max)
    y_center = 0.5_dp * (y_min + y_max)
    z_center = 0.5_dp * (z_min + z_max)
  end subroutine update_common_derived_parameters

end module mod_common_parameters
