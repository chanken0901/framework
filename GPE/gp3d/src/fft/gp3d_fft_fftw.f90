!> FFTW3をISO_C_BINDING経由で呼び出す、逐次3次元FFTバックエンド。
!> gp3d_fft参照実装と同じ公開APIを持ち、CMakeプロファイルで排他的に選択される。
module gp3d_fft
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_loc
  use gp3d_types, only: dp
  implicit none
  private

  public :: gp3d_fft_plan_t
  public :: gp3d_fft_init
  public :: gp3d_fft_forward
  public :: gp3d_fft_inverse
  public :: gp3d_fft_finalize

  integer(c_int), parameter :: FFTW_FORWARD = -1_c_int
  integer(c_int), parameter :: FFTW_BACKWARD = 1_c_int
  integer(c_int), parameter :: FFTW_ESTIMATE = 64_c_int

  type :: gp3d_fft_plan_t
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    type(c_ptr) :: forward_plan = c_null_ptr
    type(c_ptr) :: inverse_plan = c_null_ptr
  end type gp3d_fft_plan_t

  interface
    function fftw_plan_dft_3d(n0, n1, n2, input, output, sign, flags) bind(C, name="fftw_plan_dft_3d")
      import :: c_int, c_ptr
      integer(c_int), value :: n0, n1, n2
      type(c_ptr), value :: input
      type(c_ptr), value :: output
      integer(c_int), value :: sign
      integer(c_int), value :: flags
      type(c_ptr) :: fftw_plan_dft_3d
    end function fftw_plan_dft_3d

    subroutine fftw_execute_dft(plan, input, output) bind(C, name="fftw_execute_dft")
      import :: c_ptr
      type(c_ptr), value :: plan
      type(c_ptr), value :: input
      type(c_ptr), value :: output
    end subroutine fftw_execute_dft

    subroutine fftw_destroy_plan(plan) bind(C, name="fftw_destroy_plan")
      import :: c_ptr
      type(c_ptr), value :: plan
    end subroutine fftw_destroy_plan
  end interface

contains

  subroutine gp3d_fft_init(plan, nx, ny, nz, comm, rank, nprocs)
    type(gp3d_fft_plan_t), intent(out) :: plan
    integer, intent(in) :: nx, ny, nz
    integer, intent(in), optional :: comm, rank, nprocs

    complex(dp), allocatable, target :: work_in(:,:,:)
    complex(dp), allocatable, target :: work_out(:,:,:)

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) error stop "FFT dimensions must be positive"

    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    if (present(comm)) continue
    if (present(rank)) continue
    if (present(nprocs)) continue

    allocate(work_in(nx, ny, nz), work_out(nx, ny, nz))
    work_in = (0.0_dp, 0.0_dp)
    work_out = (0.0_dp, 0.0_dp)

    ! Fortran arrays are column-major.  Reversing the dimensions gives FFTW
    ! the same contiguous memory layout as psi(i,j,k).
    plan%forward_plan = fftw_plan_dft_3d( &
      int(nz, c_int), int(ny, c_int), int(nx, c_int), &
      c_loc(work_in), c_loc(work_out), FFTW_FORWARD, FFTW_ESTIMATE)
    plan%inverse_plan = fftw_plan_dft_3d( &
      int(nz, c_int), int(ny, c_int), int(nx, c_int), &
      c_loc(work_in), c_loc(work_out), FFTW_BACKWARD, FFTW_ESTIMATE)

    if (.not. c_associated(plan%forward_plan)) error stop "failed to create FFTW forward plan"
    if (.not. c_associated(plan%inverse_plan)) error stop "failed to create FFTW inverse plan"
  end subroutine gp3d_fft_init

  subroutine gp3d_fft_forward(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in), target :: input(plan%nx, plan%ny, plan%nz)
    complex(dp), intent(out), target :: output(plan%nx, plan%ny, plan%nz)

    call fftw_execute_dft(plan%forward_plan, c_loc(input), c_loc(output))
  end subroutine gp3d_fft_forward

  subroutine gp3d_fft_inverse(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in), target :: input(plan%nx, plan%ny, plan%nz)
    complex(dp), intent(out), target :: output(plan%nx, plan%ny, plan%nz)

    call fftw_execute_dft(plan%inverse_plan, c_loc(input), c_loc(output))
    output = output / real(plan%nx * plan%ny * plan%nz, dp)
  end subroutine gp3d_fft_inverse

  subroutine gp3d_fft_finalize(plan)
    type(gp3d_fft_plan_t), intent(inout) :: plan

    if (c_associated(plan%forward_plan)) call fftw_destroy_plan(plan%forward_plan)
    if (c_associated(plan%inverse_plan)) call fftw_destroy_plan(plan%inverse_plan)

    plan%forward_plan = c_null_ptr
    plan%inverse_plan = c_null_ptr
    plan%nx = 0
    plan%ny = 0
    plan%nz = 0
  end subroutine gp3d_fft_finalize

end module gp3d_fft
