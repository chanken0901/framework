module param
use mod_precision, only: dp
  implicit none

  real(dp),parameter :: rho_0 = 1.0_dp
  real(dp),parameter :: Mach  = 0.5_dp

  integer :: nx,ny,nz
  integer :: nghost         ! WENO5 needs 3 ghosts

  real(dp) :: x_min, x_max, x_center
  real(dp) :: y_min, y_max, y_center
  real(dp) :: z_min, z_max, z_center

  real(dp) :: t_max = 100.0_dp
  real(dp) :: cfl   = 0.50_dp   ! a bit safer for sharper capture

  real(dp) :: gamma = 1.4_dp
  real(dp) :: small_rho = 1.0d-12
  real(dp) :: small_p   = 1.0d-12
  integer  :: output_frequency =  10 ! Output every 20 steps for Paraview

  real(dp), parameter :: pi = acos(-1.0_dp)

end module