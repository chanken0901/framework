module mod_openmp_runtime
  implicit none
  private

  public :: nse_max_threads

contains

  integer function nse_max_threads() result(count)
    count = 1
  end function nse_max_threads

end module mod_openmp_runtime
