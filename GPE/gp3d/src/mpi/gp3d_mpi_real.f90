!> Microsoft MPIを含む実MPI実行で使用する通信コンテキストと集団通信ラッパー。
!> mainやsolverがMPI APIへ直接依存しないよう、初期化・同期・sum/maxを共通化する。
module gp3d_mpi
  implicit none
  include 'mpif.h'
  private

  public :: gp3d_mpi_t
  public :: gp3d_mpi_init
  public :: gp3d_mpi_finalize
  public :: gp3d_mpi_is_root
  public :: gp3d_mpi_barrier
  public :: gp3d_mpi_z_range
  public :: gp3d_mpi_sum_real
  public :: gp3d_mpi_max_real

  type :: gp3d_mpi_t
    logical :: enabled = .true.
    integer :: comm = MPI_COMM_WORLD
    integer :: rank = 0
    integer :: nprocs = 1
    integer :: root = 0
  end type gp3d_mpi_t

contains

  subroutine gp3d_mpi_init(ctx)
    type(gp3d_mpi_t), intent(out) :: ctx
    integer :: ierr
    logical :: initialized

    call MPI_Initialized(initialized, ierr)
    if (.not. initialized) call MPI_Init(ierr)

    ctx%enabled = .true.
    ctx%comm = MPI_COMM_WORLD
    ctx%root = 0
    call MPI_Comm_rank(ctx%comm, ctx%rank, ierr)
    call MPI_Comm_size(ctx%comm, ctx%nprocs, ierr)
  end subroutine gp3d_mpi_init

  subroutine gp3d_mpi_finalize(ctx)
    type(gp3d_mpi_t), intent(inout) :: ctx
    integer :: ierr
    logical :: finalized

    call MPI_Finalized(finalized, ierr)
    if (.not. finalized) call MPI_Finalize(ierr)
    ctx%enabled = .false.
  end subroutine gp3d_mpi_finalize

  pure logical function gp3d_mpi_is_root(ctx) result(is_root)
    type(gp3d_mpi_t), intent(in) :: ctx

    is_root = (ctx%rank == ctx%root)
  end function gp3d_mpi_is_root

  subroutine gp3d_mpi_barrier(ctx)
    type(gp3d_mpi_t), intent(in) :: ctx
    integer :: ierr

    call MPI_Barrier(ctx%comm, ierr)
  end subroutine gp3d_mpi_barrier

  subroutine gp3d_mpi_sum_real(ctx, local_value, global_value)
    type(gp3d_mpi_t), intent(in) :: ctx
    real(kind(1.0d0)), intent(in) :: local_value
    real(kind(1.0d0)), intent(out) :: global_value
    integer :: ierr

    call MPI_Allreduce(local_value, global_value, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ctx%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Allreduce failed"
  end subroutine gp3d_mpi_sum_real

  subroutine gp3d_mpi_max_real(ctx, local_value, global_value)
    type(gp3d_mpi_t), intent(in) :: ctx
    real(kind(1.0d0)), intent(in) :: local_value
    real(kind(1.0d0)), intent(out) :: global_value
    integer :: ierr

    call MPI_Allreduce(local_value, global_value, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ctx%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Allreduce failed"
  end subroutine gp3d_mpi_max_real

  pure subroutine gp3d_mpi_z_range(nz, rank, nprocs, k_start, k_end)
    integer, intent(in) :: nz, rank, nprocs
    integer, intent(out) :: k_start, k_end
    integer :: base, rest

    base = nz / nprocs
    rest = mod(nz, nprocs)
    if (rank < rest) then
      k_start = rank * (base + 1) + 1
      k_end = k_start + base
    else
      k_start = rest * (base + 1) + (rank - rest) * base + 1
      k_end = k_start + base - 1
    end if
  end subroutine gp3d_mpi_z_range

end module gp3d_mpi
