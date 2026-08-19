module mod_openmp_runtime
  use omp_lib, only : omp_get_max_threads
  implicit none
  private

  public :: nse_max_threads

contains

  integer function nse_max_threads() result(count)
    count = omp_get_max_threads()
  end function nse_max_threads

end module mod_openmp_runtime
