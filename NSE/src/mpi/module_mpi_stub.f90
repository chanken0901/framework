module module_mpi
  use mod_precision, only : dp
  implicit none
  private

  integer, public :: nprocs = 1
  integer, public :: my_rank = 0
  integer, parameter, public :: root = 0
  integer, parameter, public :: MPI_INTEGER = 1
  integer, parameter, public :: MPI_COMM_WORLD = 0

  public :: MPI_Gather
  public :: mp_barrier
  public :: mp_allminr8
  public :: mp_stop

contains

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
