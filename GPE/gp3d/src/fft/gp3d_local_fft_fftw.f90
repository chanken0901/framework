!> MPI分散FFT内の各1次元変換をFFTW3で処理する局所FFTバックエンド。
!> gp3d_local_fft_dftと同じ公開APIを持ち、CMakeでどちらか一方が選択される。
module gp3d_local_fft
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_loc
  use gp3d_types, only: dp
  implicit none
  private

  public :: gp3d_local_fft_plan_t
  public :: gp3d_local_fft_init
  public :: gp3d_local_fft_execute
  public :: gp3d_local_fft_finalize

  integer(c_int), parameter :: FFTW_FORWARD = -1_c_int
  integer(c_int), parameter :: FFTW_BACKWARD = 1_c_int
  integer(c_int), parameter :: FFTW_ESTIMATE = 64_c_int

  type :: gp3d_local_fft_plan_t
    integer :: n = 0
    type(c_ptr) :: forward_plan = c_null_ptr
    type(c_ptr) :: inverse_plan = c_null_ptr
  end type gp3d_local_fft_plan_t

  interface
    function fftw_plan_dft_1d(n, input, output, sign, flags) bind(C, name="fftw_plan_dft_1d")
      import :: c_int, c_ptr
      integer(c_int), value :: n, sign, flags
      type(c_ptr), value :: input, output
      type(c_ptr) :: fftw_plan_dft_1d
    end function fftw_plan_dft_1d

    subroutine fftw_execute_dft(plan, input, output) bind(C, name="fftw_execute_dft")
      import :: c_ptr
      type(c_ptr), value :: plan, input, output
    end subroutine fftw_execute_dft

    subroutine fftw_destroy_plan(plan) bind(C, name="fftw_destroy_plan")
      import :: c_ptr
      type(c_ptr), value :: plan
    end subroutine fftw_destroy_plan
  end interface

contains

  subroutine gp3d_local_fft_init(plan, n)
    type(gp3d_local_fft_plan_t), intent(out) :: plan
    integer, intent(in) :: n
    complex(dp), allocatable, target :: work_in(:), work_out(:)

    if (n <= 0) error stop "local FFT length must be positive"
    plan%n = n
    allocate(work_in(n), work_out(n))
    work_in = (0.0_dp, 0.0_dp)
    work_out = (0.0_dp, 0.0_dp)
    plan%forward_plan = fftw_plan_dft_1d(int(n, c_int), c_loc(work_in), c_loc(work_out), &
      FFTW_FORWARD, FFTW_ESTIMATE)
    plan%inverse_plan = fftw_plan_dft_1d(int(n, c_int), c_loc(work_in), c_loc(work_out), &
      FFTW_BACKWARD, FFTW_ESTIMATE)
    if (.not. c_associated(plan%forward_plan)) error stop "failed to create FFTW forward plan"
    if (.not. c_associated(plan%inverse_plan)) error stop "failed to create FFTW inverse plan"
  end subroutine gp3d_local_fft_init

  subroutine gp3d_local_fft_execute(plan, input, output, sign)
    type(gp3d_local_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in), target :: input(plan%n)
    complex(dp), intent(out), target :: output(plan%n)
    integer, intent(in) :: sign

    select case (sign)
    case (-1)
      call fftw_execute_dft(plan%forward_plan, c_loc(input), c_loc(output))
    case (1)
      call fftw_execute_dft(plan%inverse_plan, c_loc(input), c_loc(output))
    case default
      error stop "local FFT sign must be -1 or 1"
    end select
  end subroutine gp3d_local_fft_execute

  subroutine gp3d_local_fft_finalize(plan)
    type(gp3d_local_fft_plan_t), intent(inout) :: plan

    if (c_associated(plan%forward_plan)) call fftw_destroy_plan(plan%forward_plan)
    if (c_associated(plan%inverse_plan)) call fftw_destroy_plan(plan%inverse_plan)
    plan%forward_plan = c_null_ptr
    plan%inverse_plan = c_null_ptr
    plan%n = 0
  end subroutine gp3d_local_fft_finalize

end module gp3d_local_fft
