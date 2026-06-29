module mod_constants
  use mod_precision, only : dp
  implicit none
  private
  real(dp), parameter :: pi = acos(-1.0_dp)
end module mod_precision