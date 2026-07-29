!> MPI zスラブ分割上で3次元FFTを行う分散バックエンド。
!> x-y局所変換後にAlltoallvでyスラブへ転置してz変換し、zスラブへ戻す。
!> 全領域をrank 0へ集約せず、局所1次元FFTはgp3d_local_fftへ委譲する。
module gp3d_fft
  use gp3d_types, only: dp
  use gp3d_local_fft, only: gp3d_local_fft_plan_t, gp3d_local_fft_init, &
    gp3d_local_fft_execute, gp3d_local_fft_finalize
  use gp3d_openmp, only: gp3d_openmp_active
  implicit none
  include 'mpif.h'
  private

  public :: gp3d_fft_plan_t
  public :: gp3d_fft_init
  public :: gp3d_fft_forward
  public :: gp3d_fft_inverse
  public :: gp3d_fft_finalize

  type :: gp3d_fft_plan_t
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: local_nz = 0
    integer :: local_ny = 0
    integer :: k_start = 1
    integer :: j_start = 1
    integer :: comm = MPI_COMM_WORLD
    integer :: rank = 0
    integer :: nprocs = 1
    integer, allocatable :: z_starts(:), z_counts(:)
    integer, allocatable :: y_starts(:), y_counts(:)
    integer, allocatable :: zy_sendcounts(:), zy_senddispls(:)
    integer, allocatable :: zy_recvcounts(:), zy_recvdispls(:)
    type(gp3d_local_fft_plan_t) :: x_plan
    type(gp3d_local_fft_plan_t) :: y_plan
    type(gp3d_local_fft_plan_t) :: z_plan
  end type gp3d_fft_plan_t

contains

  subroutine gp3d_fft_init(plan, nx, ny, nz, comm, rank, nprocs)
    type(gp3d_fft_plan_t), intent(out) :: plan
    integer, intent(in) :: nx, ny, nz, comm, rank, nprocs
    integer :: r

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) error stop "FFT dimensions must be positive"
    if (nprocs <= 0 .or. rank < 0 .or. rank >= nprocs) error stop "invalid MPI FFT topology"
    if (nprocs > ny .or. nprocs > nz) then
      error stop "MPI process count must not exceed ny or nz for slab FFT"
    end if

    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    plan%comm = comm
    plan%rank = rank
    plan%nprocs = nprocs
    allocate(plan%z_starts(nprocs), plan%z_counts(nprocs))
    allocate(plan%y_starts(nprocs), plan%y_counts(nprocs))
    allocate(plan%zy_sendcounts(nprocs), plan%zy_senddispls(nprocs))
    allocate(plan%zy_recvcounts(nprocs), plan%zy_recvdispls(nprocs))

    do r = 0, nprocs - 1
      call block_range(nz, r, nprocs, plan%z_starts(r + 1), plan%z_counts(r + 1))
      call block_range(ny, r, nprocs, plan%y_starts(r + 1), plan%y_counts(r + 1))
    end do
    plan%local_nz = plan%z_counts(rank + 1)
    plan%local_ny = plan%y_counts(rank + 1)
    plan%k_start = plan%z_starts(rank + 1)
    plan%j_start = plan%y_starts(rank + 1)

    do r = 1, nprocs
      plan%zy_sendcounts(r) = nx * plan%y_counts(r) * plan%local_nz
      plan%zy_recvcounts(r) = nx * plan%local_ny * plan%z_counts(r)
    end do
    call make_displacements(plan%zy_sendcounts, plan%zy_senddispls)
    call make_displacements(plan%zy_recvcounts, plan%zy_recvdispls)

    call gp3d_local_fft_init(plan%x_plan, nx)
    call gp3d_local_fft_init(plan%y_plan, ny)
    call gp3d_local_fft_init(plan%z_plan, nz)
  end subroutine gp3d_fft_init

  subroutine gp3d_fft_forward(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%local_nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%local_nz)

    call distributed_transform(plan, input, output, -1, .false.)
  end subroutine gp3d_fft_forward

  subroutine gp3d_fft_inverse(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%local_nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%local_nz)

    call distributed_transform(plan, input, output, 1, .true.)
  end subroutine gp3d_fft_inverse

  subroutine gp3d_fft_finalize(plan)
    type(gp3d_fft_plan_t), intent(inout) :: plan

    call gp3d_local_fft_finalize(plan%x_plan)
    call gp3d_local_fft_finalize(plan%y_plan)
    call gp3d_local_fft_finalize(plan%z_plan)
    if (allocated(plan%z_starts)) deallocate(plan%z_starts, plan%z_counts)
    if (allocated(plan%y_starts)) deallocate(plan%y_starts, plan%y_counts)
    if (allocated(plan%zy_sendcounts)) deallocate(plan%zy_sendcounts, plan%zy_senddispls)
    if (allocated(plan%zy_recvcounts)) deallocate(plan%zy_recvcounts, plan%zy_recvdispls)
    plan%nx = 0
    plan%ny = 0
    plan%nz = 0
    plan%local_nz = 0
    plan%local_ny = 0
  end subroutine gp3d_fft_finalize

  subroutine distributed_transform(plan, input, output, sign, normalize)
    ! 1) x-y変換、2) z->y分散転置、3) z変換、4) y->z逆転置の順に処理する。
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%local_nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%local_nz)
    integer, intent(in) :: sign
    logical, intent(in) :: normalize
    complex(dp), allocatable :: z_slab(:,:,:), y_slab(:,:,:)
    integer :: i, j, k

    allocate(z_slab(plan%nx, plan%ny, plan%local_nz))
    allocate(y_slab(plan%nx, plan%local_ny, plan%nz))
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active)
    do k = 1, plan%local_nz
      do j = 1, plan%ny
        do i = 1, plan%nx
          z_slab(i,j,k) = input(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

    call transform_xy(plan, z_slab, sign)
    call transpose_z_to_y(plan, z_slab, y_slab)
    call transform_z(plan, y_slab, sign)
    call transpose_y_to_z(plan, y_slab, output)

    if (normalize) then
      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active)
      do k = 1, plan%local_nz
        do j = 1, plan%ny
          do i = 1, plan%nx
            output(i,j,k) = output(i,j,k) / real(plan%nx * plan%ny * plan%nz, dp)
          end do
        end do
      end do
      !$omp end parallel do
    end if
    deallocate(z_slab, y_slab)
  end subroutine distributed_transform

  subroutine transform_xy(plan, slab, sign)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(inout) :: slab(plan%nx, plan%ny, plan%local_nz)
    integer, intent(in) :: sign
    complex(dp), allocatable :: line_in(:), line_out(:)
    integer :: i, j, k

    !$omp parallel if(gp3d_openmp_active) private(i, j, k, line_in, line_out)
    allocate(line_in(max(plan%nx, plan%ny)), line_out(max(plan%nx, plan%ny)))
    !$omp do collapse(2) schedule(static)
    do k = 1, plan%local_nz
      do j = 1, plan%ny
        line_in(1:plan%nx) = slab(:,j,k)
        call gp3d_local_fft_execute(plan%x_plan, line_in(1:plan%nx), line_out(1:plan%nx), sign)
        slab(:,j,k) = line_out(1:plan%nx)
      end do
    end do
    !$omp end do
    !$omp do collapse(2) schedule(static)
    do k = 1, plan%local_nz
      do i = 1, plan%nx
        line_in(1:plan%ny) = slab(i,:,k)
        call gp3d_local_fft_execute(plan%y_plan, line_in(1:plan%ny), line_out(1:plan%ny), sign)
        slab(i,:,k) = line_out(1:plan%ny)
      end do
    end do
    !$omp end do
    deallocate(line_in, line_out)
    !$omp end parallel
  end subroutine transform_xy

  subroutine transform_z(plan, slab, sign)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(inout) :: slab(plan%nx, plan%local_ny, plan%nz)
    integer, intent(in) :: sign
    complex(dp), allocatable :: line_in(:), line_out(:)
    integer :: i, j

    !$omp parallel if(gp3d_openmp_active) private(i, j, line_in, line_out)
    allocate(line_in(plan%nz), line_out(plan%nz))
    !$omp do collapse(2) schedule(static)
    do j = 1, plan%local_ny
      do i = 1, plan%nx
        line_in(1:plan%nz) = slab(i,j,:)
        call gp3d_local_fft_execute(plan%z_plan, line_in(1:plan%nz), line_out(1:plan%nz), sign)
        slab(i,j,:) = line_out(1:plan%nz)
      end do
    end do
    !$omp end do
    deallocate(line_in, line_out)
    !$omp end parallel
  end subroutine transform_z

  subroutine transpose_z_to_y(plan, z_slab, y_slab)
    ! 各rankのzスラブを、z方向の線が局所的に揃うyスラブ配置へAlltoallvする。
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: z_slab(plan%nx, plan%ny, plan%local_nz)
    complex(dp), intent(out) :: y_slab(plan%nx, plan%local_ny, plan%nz)
    complex(dp), allocatable :: sendbuf(:), recvbuf(:)
    integer :: dest, source, i, j, k, offset, ierr

    allocate(sendbuf(sum(plan%zy_sendcounts)), recvbuf(sum(plan%zy_recvcounts)))
    !$omp parallel do schedule(static) if(gp3d_openmp_active) private(offset, i, j, k)
    do dest = 1, plan%nprocs
      offset = plan%zy_senddispls(dest) + 1
      do k = 1, plan%local_nz
        do j = plan%y_starts(dest), plan%y_starts(dest) + plan%y_counts(dest) - 1
          do i = 1, plan%nx
            sendbuf(offset) = z_slab(i,j,k)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call MPI_Alltoallv(sendbuf, plan%zy_sendcounts, plan%zy_senddispls, MPI_DOUBLE_COMPLEX, &
      recvbuf, plan%zy_recvcounts, plan%zy_recvdispls, MPI_DOUBLE_COMPLEX, plan%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Alltoallv z-to-y transpose failed"

    !$omp parallel do schedule(static) if(gp3d_openmp_active) private(offset, i, j, k)
    do source = 1, plan%nprocs
      offset = plan%zy_recvdispls(source) + 1
      do k = plan%z_starts(source), plan%z_starts(source) + plan%z_counts(source) - 1
        do j = 1, plan%local_ny
          do i = 1, plan%nx
            y_slab(i,j,k) = recvbuf(offset)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(sendbuf, recvbuf)
  end subroutine transpose_z_to_y

  subroutine transpose_y_to_z(plan, y_slab, z_slab)
    ! スペクトルデータをソルバーが所有する元のzスラブ配置へ戻す。
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: y_slab(plan%nx, plan%local_ny, plan%nz)
    complex(dp), intent(out) :: z_slab(plan%nx, plan%ny, plan%local_nz)
    complex(dp), allocatable :: sendbuf(:), recvbuf(:)
    integer :: dest, source, i, j, k, offset, ierr

    allocate(sendbuf(sum(plan%zy_recvcounts)), recvbuf(sum(plan%zy_sendcounts)))
    !$omp parallel do schedule(static) if(gp3d_openmp_active) private(offset, i, j, k)
    do dest = 1, plan%nprocs
      offset = plan%zy_recvdispls(dest) + 1
      do k = plan%z_starts(dest), plan%z_starts(dest) + plan%z_counts(dest) - 1
        do j = 1, plan%local_ny
          do i = 1, plan%nx
            sendbuf(offset) = y_slab(i,j,k)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call MPI_Alltoallv(sendbuf, plan%zy_recvcounts, plan%zy_recvdispls, MPI_DOUBLE_COMPLEX, &
      recvbuf, plan%zy_sendcounts, plan%zy_senddispls, MPI_DOUBLE_COMPLEX, plan%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Alltoallv y-to-z transpose failed"

    !$omp parallel do schedule(static) if(gp3d_openmp_active) private(offset, i, j, k)
    do source = 1, plan%nprocs
      offset = plan%zy_senddispls(source) + 1
      do k = 1, plan%local_nz
        do j = plan%y_starts(source), plan%y_starts(source) + plan%y_counts(source) - 1
          do i = 1, plan%nx
            z_slab(i,j,k) = recvbuf(offset)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(sendbuf, recvbuf)
  end subroutine transpose_y_to_z

  pure subroutine block_range(n, rank, nprocs, start_index, count)
    integer, intent(in) :: n, rank, nprocs
    integer, intent(out) :: start_index, count
    integer :: base, rest

    base = n / nprocs
    rest = mod(n, nprocs)
    if (rank < rest) then
      count = base + 1
      start_index = rank * (base + 1) + 1
    else
      count = base
      start_index = rest * (base + 1) + (rank - rest) * base + 1
    end if
  end subroutine block_range

  pure subroutine make_displacements(counts, displacements)
    integer, intent(in) :: counts(:)
    integer, intent(out) :: displacements(size(counts))
    integer :: i

    displacements(1) = 0
    do i = 2, size(counts)
      displacements(i) = displacements(i - 1) + counts(i - 1)
    end do
  end subroutine make_displacements

end module gp3d_fft
