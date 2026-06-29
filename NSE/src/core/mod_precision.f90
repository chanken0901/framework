module mod_precision
  implicit none
  private
  public :: dp
  integer, parameter :: dp = selected_real_kind(15,307)
end module mod_precision