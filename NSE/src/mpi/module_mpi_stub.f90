module module_mpi
  use mod_precision, only : dp
  implicit none
  private

  integer, public :: nprocs = 1
  integer, public :: my_rank = 0
  integer, parameter, public :: root = 0
  integer, parameter, public :: MPI_INTEGER = 1
  integer, parameter, public :: MPI_DOUBLE_PRECISION = 2
  integer, parameter, public :: MPI_COMM_WORLD = 0
  integer, parameter, public :: MPI_STATUS_SIZE = 5
  integer, public :: ndiv_nx = 1
  integer, public :: ndiv_ny = 1
  integer, public :: ndiv_nz = 1
  integer, public :: i_sta = 1, i_end = 1
  integer, public :: j_sta = 1, j_end = 1
  integer, public :: k_sta = 1, k_end = 1
  integer, allocatable, public :: itable(:,:)
  integer, allocatable, public :: jjsta(:), jjend(:)
  integer, allocatable, public :: kksta(:), kkend(:)

  public :: MPI_Init
  public :: MPI_Comm_size
  public :: MPI_Comm_rank
  public :: MPI_Finalize
  public :: MPI_Wtime
  public :: MPI_Gather
  public :: mp_setup_division
  public :: mp_send_recv_pre_r8_Vec
  public :: mp_sendrecv_r8
  public :: mp_barrier
  public :: mp_allminr8
  public :: mp_stop

contains

  subroutine MPI_Init(ierr)
    integer, intent(out) :: ierr

    ierr = 0
  end subroutine MPI_Init

  subroutine MPI_Comm_size(communicator, size_out, ierr)
    integer, intent(in) :: communicator
    integer, intent(out) :: size_out, ierr

    size_out = 1
    ierr = 0
  end subroutine MPI_Comm_size

  subroutine MPI_Comm_rank(communicator, rank_out, ierr)
    integer, intent(in) :: communicator
    integer, intent(out) :: rank_out, ierr

    rank_out = 0
    ierr = 0
  end subroutine MPI_Comm_rank

  subroutine MPI_Finalize(ierr)
    integer, intent(out) :: ierr

    ierr = 0
  end subroutine MPI_Finalize

  real(dp) function MPI_Wtime() result(seconds)
    integer :: count, rate

    call system_clock(count, rate)
    seconds = real(count, dp) / real(rate, dp)
  end function MPI_Wtime

  subroutine mp_setup_division(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz

    ndiv_nx = 1
    ndiv_ny = 1
    ndiv_nz = 1
    i_sta = 1
    i_end = nx
    j_sta = 1
    j_end = ny
    k_sta = 1
    k_end = nz

    if (allocated(itable)) deallocate(itable)
    if (allocated(jjsta)) deallocate(jjsta, jjend)
    if (allocated(kksta)) deallocate(kksta, kkend)
    allocate(itable(0:0,0:0), jjsta(0:0), jjend(0:0))
    allocate(kksta(0:0), kkend(0:0))
    itable(0,0) = 0
    jjsta(0) = 1
    jjend(0) = ny
    kksta(0) = 1
    kkend(0) = nz
  end subroutine mp_setup_division

  subroutine MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, &
      recvtype, root_rank, communicator, ierr)
    integer, intent(in) :: sendbuf(:)
    integer, intent(in) :: sendcount, sendtype, recvcount, recvtype
    integer, intent(in) :: root_rank, communicator
    integer, intent(out) :: recvbuf(:,:)
    integer, intent(out) :: ierr

    if (sendcount /= recvcount) error stop "serial MPI_Gather count mismatch"
    if (size(recvbuf,1) < sendcount .or. size(recvbuf,2) < 1) then
      error stop "serial MPI_Gather receive buffer is too small"
    end if
    recvbuf(1:sendcount,1) = sendbuf(1:sendcount)
    ierr = 0
  end subroutine MPI_Gather

  subroutine mp_sendrecv_r8(sendbuf, sendcount, sendtype, destination, &
      sendtag, recvbuf, recvcount, recvtype, source, recvtag, communicator, &
      status, ierr)
    real(dp), intent(in) :: sendbuf(*)
    integer, intent(in) :: sendcount, sendtype, destination, sendtag
    real(dp), intent(out) :: recvbuf(*)
    integer, intent(in) :: recvcount, recvtype, source, recvtag, communicator
    integer, intent(out) :: status(*)
    integer, intent(out) :: ierr
    integer :: item

    if (sendcount /= recvcount) error stop 'serial send/receive count mismatch'
    do item = 1, sendcount
      recvbuf(item) = sendbuf(item)
    end do
    status(1:MPI_STATUS_SIZE) = 0
    ierr = 0
  end subroutine mp_sendrecv_r8

  subroutine mp_send_recv_pre_r8_Vec(field, mode, nx1, nx2, ny1, ny2, &
      nz1, nz2)
    integer, intent(in) :: mode, nx1, nx2, ny1, ny2, nz1, nz2
    real(dp), intent(inout) :: field(nx1:nx2,ny1:ny2,nz1:nz2,5)

    ! Periodic ghost cells are filled locally by mod_nse_boundary.
  end subroutine mp_send_recv_pre_r8_Vec

  subroutine mp_barrier()
  end subroutine mp_barrier

  subroutine mp_allminr8(value)
    real(dp), intent(inout) :: value
  end subroutine mp_allminr8

  subroutine mp_stop(code)
    integer, intent(in) :: code

    if (code == 0) stop
    error stop "serial MPI compatibility layer requested an error stop"
  end subroutine mp_stop

end module module_mpi
