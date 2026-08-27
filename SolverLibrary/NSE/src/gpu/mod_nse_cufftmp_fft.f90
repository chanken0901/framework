module mod_nse_cufftmp_fft
  use, intrinsic :: iso_c_binding, only : c_char, c_double_complex, c_int, &
    c_null_char
  use, intrinsic :: iso_fortran_env, only : error_unit
  use mod_precision, only : dp
  use module_mpi, only : jjsta, jjend, kksta, kkend, MPI_COMM_WORLD
  implicit none
  private

  integer, parameter, public :: mytype = dp
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1

  type, public :: decomp_info
    integer :: xst(3) = 0
    integer :: xen(3) = -1
  end type decomp_info

  integer, public :: xstart(3) = 0
  integer, public :: xend(3) = -1
  integer, public :: zstart(3) = 0
  integer, public :: zend(3) = -1

  type(decomp_info), target :: physical_layout
  integer :: global_shape(3) = 0
  logical :: initialized = .false.

  public :: decomp_2d_init, decomp_2d_finalize
  public :: decomp_2d_fft_init, decomp_2d_fft_finalize
  public :: decomp_2d_fft_3d, decomp_2d_fft_get_ph
  public :: alloc_x, alloc_z

  interface
    function c_cufftmp_fft_initialize(nx, ny, nz, ylo, yhi, zlo, zhi, &
        communicator) bind(C, name='nse_cufftmp_fft_initialize') result(status)
      import :: c_int
      integer(c_int), value :: nx, ny, nz, ylo, yhi, zlo, zhi
      integer(c_int), value :: communicator
      integer(c_int) :: status
    end function c_cufftmp_fft_initialize

    function c_cufftmp_fft_execute(input, output, direction) &
        bind(C, name='nse_cufftmp_fft_execute') result(status)
      import :: c_double_complex, c_int
      complex(c_double_complex), intent(in) :: input(*)
      complex(c_double_complex), intent(out) :: output(*)
      integer(c_int), value :: direction
      integer(c_int) :: status
    end function c_cufftmp_fft_execute

    subroutine c_cufftmp_fft_finalize() &
        bind(C, name='nse_cufftmp_fft_finalize')
    end subroutine c_cufftmp_fft_finalize

    subroutine c_cuda_get_last_error(buffer, buffer_size) &
        bind(C, name='nse_cuda_get_last_error')
      import :: c_char, c_int
      character(c_char), intent(out) :: buffer(*)
      integer(c_int), value :: buffer_size
    end subroutine c_cuda_get_last_error
  end interface

contains

  subroutine decomp_2d_init(nx, ny, nz, p_row, p_col, complex_pool)
    integer, intent(in) :: nx, ny, nz, p_row, p_col
    logical, intent(in), optional :: complex_pool

    if (initialized) error stop 'cuFFTMp decomposition is already initialized'
    if (p_row <= 0 .or. p_col <= 0) then
      error stop 'cuFFTMp requires a positive MPI process grid'
    end if
    ! The argument is retained for source compatibility with 2DECOMP&FFT.
    if (present(complex_pool)) then
      if (.not. complex_pool) continue
    end if

    global_shape = [nx, ny, nz]
    xstart = [1, jjsta, kksta]
    xend = [nx, jjend, kkend]
    ! cuFFTMp is configured with equal input and output pencils.  This lets the
    ! existing distributed HIT mathematics use its spectral zstart/zend names
    ! without an extra redistribution or a root gather.
    zstart = xstart
    zend = xend
    physical_layout%xst = xstart
    physical_layout%xen = xend
  end subroutine decomp_2d_init

  subroutine decomp_2d_fft_init(layout)
    integer, intent(in) :: layout
    integer(c_int) :: status

    if (layout /= PHYSICAL_IN_X) then
      error stop 'cuFFTMp NSE backend requires PHYSICAL_IN_X'
    end if
    status = c_cufftmp_fft_initialize(int(global_shape(1),c_int), &
      int(global_shape(2),c_int), int(global_shape(3),c_int), &
      int(xstart(2),c_int), int(xend(2),c_int), &
      int(xstart(3),c_int), int(xend(3),c_int), &
      int(MPI_COMM_WORLD,c_int))
    call require_success(status, 'initialize cuFFTMp HIT plan')
    initialized = .true.
  end subroutine decomp_2d_fft_init

  function decomp_2d_fft_get_ph() result(layout)
    type(decomp_info), pointer :: layout
    layout => physical_layout
  end function decomp_2d_fft_get_ph

  subroutine alloc_x(field, layout, complex_data)
    complex(mytype), allocatable, intent(out) :: field(:,:,:)
    type(decomp_info), intent(in) :: layout
    logical, intent(in), optional :: complex_data

    if (present(complex_data)) then
      if (.not. complex_data) continue
    end if
    allocate(field(layout%xst(1):layout%xen(1), &
      layout%xst(2):layout%xen(2), layout%xst(3):layout%xen(3)))
  end subroutine alloc_x

  subroutine alloc_z(field, layout, complex_data)
    complex(mytype), allocatable, intent(out) :: field(:,:,:)
    type(decomp_info), intent(in) :: layout
    logical, intent(in), optional :: complex_data

    call alloc_x(field, layout, complex_data)
  end subroutine alloc_z

  subroutine decomp_2d_fft_3d(input, output, direction)
    complex(mytype), contiguous, intent(in) :: input(:,:,:)
    complex(mytype), contiguous, intent(out) :: output(:,:,:)
    integer, intent(in) :: direction
    integer(c_int) :: status

    if (.not. initialized) error stop 'cuFFTMp plan is not initialized'
    if (size(input) /= size(output)) then
      error stop 'cuFFTMp input/output local sizes differ'
    end if
    if (direction /= DECOMP_2D_FFT_FORWARD .and. &
        direction /= DECOMP_2D_FFT_BACKWARD) then
      error stop 'invalid cuFFTMp transform direction'
    end if
    status = c_cufftmp_fft_execute(input, output, int(direction,c_int))
    call require_success(status, 'execute distributed cuFFTMp transform')
  end subroutine decomp_2d_fft_3d

  subroutine decomp_2d_fft_finalize()
    if (.not. initialized) return
    call c_cufftmp_fft_finalize()
    initialized = .false.
  end subroutine decomp_2d_fft_finalize

  subroutine decomp_2d_finalize()
    global_shape = 0
    xstart = 0
    xend = -1
    zstart = 0
    zend = -1
  end subroutine decomp_2d_finalize

  subroutine require_success(status, operation)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: operation
    character(c_char) :: raw(1024)
    character(len=1024) :: message
    integer :: i

    if (status == 0_c_int) return
    raw = c_null_char
    call c_cuda_get_last_error(raw, int(size(raw),c_int))
    message = ''
    do i = 1, size(raw)
      if (raw(i) == c_null_char) exit
      message(i:i) = raw(i)
    end do
    write(error_unit,'(A)') trim(operation)//': '//trim(message)
    error stop 'cuFFTMp operation failed'
  end subroutine require_success

end module mod_nse_cufftmp_fft
