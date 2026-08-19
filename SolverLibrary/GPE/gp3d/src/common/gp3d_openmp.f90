!> Runtime OpenMP control shared by the CPU solver modules.
!> The same OpenMP-capable executable can run with one thread or many threads.
module gp3d_openmp
  !$ use omp_lib, only: omp_get_max_threads, omp_set_dynamic, omp_set_num_threads
  implicit none
  private

  logical, save, public :: gp3d_openmp_active = .false.
  logical, save :: openmp_compiled = .false.
  integer, save :: active_threads = 1

  public :: gp3d_openmp_configure
  public :: gp3d_openmp_is_compiled
  public :: gp3d_openmp_thread_count

contains

  subroutine gp3d_openmp_configure(requested)
    logical, intent(in) :: requested

    openmp_compiled = .false.
    !$ openmp_compiled = .true.
    if (requested .and. .not. openmp_compiled) then
      error stop "OpenMP was requested at run time, but this executable was built without OpenMP"
    end if

    gp3d_openmp_active = requested .and. openmp_compiled
    active_threads = 1
    !$ call omp_set_dynamic(.false.)
    if (gp3d_openmp_active) then
      !$ active_threads = omp_get_max_threads()
    else
      !$ call omp_set_num_threads(1)
    end if
  end subroutine gp3d_openmp_configure

  logical function gp3d_openmp_is_compiled() result(compiled)
    compiled = openmp_compiled
  end function gp3d_openmp_is_compiled

  integer function gp3d_openmp_thread_count() result(thread_count)
    thread_count = active_threads
  end function gp3d_openmp_thread_count

end module gp3d_openmp
