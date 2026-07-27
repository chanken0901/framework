!> MPIをリンクしない逐次実行用の互換スタブ。
!> 実MPI版と同じ公開APIを持ち、rank=0、nprocs=1として通信を無処理化する。
module gp3d_mpi
  implicit none
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
    logical :: enabled = .false.
    integer :: comm = 0
    integer :: rank = 0
    integer :: nprocs = 1
    integer :: root = 0
  end type gp3d_mpi_t

contains

  subroutine gp3d_mpi_init(ctx)
    type(gp3d_mpi_t), intent(out) :: ctx

    ctx%enabled = .false.
    ctx%comm = 0
    ctx%rank = 0
    ctx%nprocs = 1
    ctx%root = 0
  end subroutine gp3d_mpi_init

  subroutine gp3d_mpi_finalize(ctx)
    type(gp3d_mpi_t), intent(inout) :: ctx

    ctx%enabled = .false.
  end subroutine gp3d_mpi_finalize

  pure logical function gp3d_mpi_is_root(ctx) result(is_root)
    type(gp3d_mpi_t), intent(in) :: ctx

    is_root = (ctx%rank == ctx%root)
  end function gp3d_mpi_is_root

  subroutine gp3d_mpi_barrier(ctx)
    type(gp3d_mpi_t), intent(in) :: ctx

    if (ctx%enabled) continue
  end subroutine gp3d_mpi_barrier

  subroutine gp3d_mpi_sum_real(ctx, local_value, global_value)
    type(gp3d_mpi_t), intent(in) :: ctx
    real(kind(1.0d0)), intent(in) :: local_value
    real(kind(1.0d0)), intent(out) :: global_value

    if (ctx%enabled) continue
    global_value = local_value
  end subroutine gp3d_mpi_sum_real

  subroutine gp3d_mpi_max_real(ctx, local_value, global_value)
    type(gp3d_mpi_t), intent(in) :: ctx
    real(kind(1.0d0)), intent(in) :: local_value
    real(kind(1.0d0)), intent(out) :: global_value

    if (ctx%enabled) continue
    global_value = local_value
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
