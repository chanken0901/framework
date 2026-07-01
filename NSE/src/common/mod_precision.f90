module mod_precision
  use, intrinsic :: iso_fortran_env, only : real64, int32, int64
  implicit none
  public

  integer, parameter :: dp = real64
  integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64
end module mod_precision
