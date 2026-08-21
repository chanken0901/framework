!> 2DECOMP&FFTを用いたMPIペンシル分割FFTバックエンド。
!>
!> GP3D本体とSLF入出力が保持するzスラブ配列を公開APIに残し、FFTの入口で
!> Xペンシルへ、出口でZペンシルからzスラブへ再分配する。これにより従来の
!> スラブバックエンドとファイル互換性を保ったまま、FFT内部の2次元プロセス
!> 格子を選択できる。
module gp3d_fft
  use gp3d_types, only: dp
  use gp3d_openmp, only: gp3d_openmp_active
  use decomp_2d, only: decomp_info, decomp_2d_init, decomp_2d_finalize, &
    decomp_log, get_decomp_dims, xstart, xend, zstart, zend, alloc_x, alloc_z
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_finalize, &
    decomp_2d_fft_get_ph, decomp_2d_fft_3d
  use decomp_2d_constants, only: mytype, complex_type, D2D_LOG_QUIET, &
    PHYSICAL_IN_X, DECOMP_2D_FFT_FORWARD, DECOMP_2D_FFT_BACKWARD
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
    integer :: k_start = 1
    integer :: comm = MPI_COMM_WORLD
    integer :: rank = 0
    integer :: nprocs = 1
    integer :: p_row = 0
    integer :: p_col = 0
    logical :: initialized = .false.
    integer, allocatable :: slab_starts(:), slab_counts(:)
    integer, allocatable :: x_starts(:,:), x_ends(:,:)
    integer, allocatable :: z_starts(:,:), z_ends(:,:)
    type(decomp_info), pointer :: ph => null()
  end type gp3d_fft_plan_t

contains

  subroutine gp3d_fft_init(plan, nx, ny, nz, comm, rank, nprocs)
    type(gp3d_fft_plan_t), intent(out) :: plan
    integer, intent(in) :: nx, ny, nz, comm, rank, nprocs
    integer :: ierr, r
    integer :: process_grid(2)

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) error stop "FFT dimensions must be positive"
    if (nprocs <= 0 .or. rank < 0 .or. rank >= nprocs) error stop "invalid MPI FFT topology"
    ! GP3Dの保存配列は互換性維持のためzスラブのままなので、この制約は残る。
    if (nprocs > nz) error stop "MPI process count must not exceed nz for GP3D slab storage"
    if (kind(0.0_dp) /= mytype) then
      error stop "GP3D and 2DECOMP&FFT floating-point precision do not match"
    end if

    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    plan%comm = comm
    plan%rank = rank
    plan%nprocs = nprocs
    allocate(plan%slab_starts(nprocs), plan%slab_counts(nprocs))
    allocate(plan%x_starts(3,nprocs), plan%x_ends(3,nprocs))
    allocate(plan%z_starts(3,nprocs), plan%z_ends(3,nprocs))
    do r = 0, nprocs - 1
      call block_range(nz, r, nprocs, plan%slab_starts(r + 1), plan%slab_counts(r + 1))
    end do
    plan%local_nz = plan%slab_counts(rank + 1)
    plan%k_start = plan%slab_starts(rank + 1)

    ! p_row=p_col=0 lets 2DECOMP&FFT choose the closest valid factorisation.
    plan%p_row = 0
    plan%p_col = 0
    decomp_log = D2D_LOG_QUIET
    call decomp_2d_init(nx, ny, nz, plan%p_row, plan%p_col, &
      comm=comm, complex_pool=.true.)
    call decomp_2d_fft_init(PHYSICAL_IN_X)
    plan%ph => decomp_2d_fft_get_ph()

    call MPI_Allgather(xstart, 3, MPI_INTEGER, plan%x_starts, 3, MPI_INTEGER, comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Allgather for X-pencil starts failed"
    call MPI_Allgather(xend, 3, MPI_INTEGER, plan%x_ends, 3, MPI_INTEGER, comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Allgather for X-pencil ends failed"
    call MPI_Allgather(zstart, 3, MPI_INTEGER, plan%z_starts, 3, MPI_INTEGER, comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Allgather for Z-pencil starts failed"
    call MPI_Allgather(zend, 3, MPI_INTEGER, plan%z_ends, 3, MPI_INTEGER, comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Allgather for Z-pencil ends failed"

    process_grid = get_decomp_dims()
    if (any(process_grid /= [plan%p_row, plan%p_col])) then
      error stop "2DECOMP&FFT returned an inconsistent process grid"
    end if
    plan%initialized = .true.
    if (rank == 0) then
      write(*,'(a,i0,a,i0)') "# FFT decomposition=pencil backend=2decomp_fftw process_grid=", &
        plan%p_row, "x", plan%p_col
    end if
  end subroutine gp3d_fft_init

  subroutine gp3d_fft_forward(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%local_nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%local_nz)
    complex(mytype), allocatable :: x_pencil(:,:,:), z_pencil(:,:,:)

    if (.not. plan%initialized) error stop "pencil FFT plan is not initialized"
    call alloc_x(x_pencil, plan%ph, .true.)
    call alloc_z(z_pencil, plan%ph, .true.)
    call slab_to_pencil(plan, input, x_pencil, plan%x_starts, plan%x_ends)
    call decomp_2d_fft_3d(x_pencil, z_pencil, DECOMP_2D_FFT_FORWARD)
    call pencil_to_slab(plan, z_pencil, output, plan%z_starts, plan%z_ends)
    deallocate(x_pencil, z_pencil)
  end subroutine gp3d_fft_forward

  subroutine gp3d_fft_inverse(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%local_nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%local_nz)
    complex(mytype), allocatable :: x_pencil(:,:,:), z_pencil(:,:,:)
    integer :: i, j, k
    real(dp) :: scale

    if (.not. plan%initialized) error stop "pencil FFT plan is not initialized"
    call alloc_x(x_pencil, plan%ph, .true.)
    call alloc_z(z_pencil, plan%ph, .true.)
    call slab_to_pencil(plan, input, z_pencil, plan%z_starts, plan%z_ends)
    call decomp_2d_fft_3d(z_pencil, x_pencil, DECOMP_2D_FFT_BACKWARD)
    call pencil_to_slab(plan, x_pencil, output, plan%x_starts, plan%x_ends)

    scale = 1.0_dp / (real(plan%nx, dp) * real(plan%ny, dp) * real(plan%nz, dp))
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active)
    do k = 1, plan%local_nz
      do j = 1, plan%ny
        do i = 1, plan%nx
          output(i,j,k) = output(i,j,k) * scale
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(x_pencil, z_pencil)
  end subroutine gp3d_fft_inverse

  subroutine gp3d_fft_finalize(plan)
    type(gp3d_fft_plan_t), intent(inout) :: plan

    if (plan%initialized) then
      nullify(plan%ph)
      call decomp_2d_fft_finalize
      call decomp_2d_finalize
    end if
    if (allocated(plan%slab_starts)) deallocate(plan%slab_starts, plan%slab_counts)
    if (allocated(plan%x_starts)) deallocate(plan%x_starts, plan%x_ends)
    if (allocated(plan%z_starts)) deallocate(plan%z_starts, plan%z_ends)
    plan%nx = 0
    plan%ny = 0
    plan%nz = 0
    plan%local_nz = 0
    plan%initialized = .false.
  end subroutine gp3d_fft_finalize

  subroutine slab_to_pencil(plan, slab, pencil, starts, ends)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: slab(plan%nx, plan%ny, plan%local_nz)
    complex(mytype), intent(out) :: pencil(:,:,:)
    integer, intent(in) :: starts(3,plan%nprocs), ends(3,plan%nprocs)
    complex(dp), allocatable :: sendbuf(:), recvbuf(:)
    integer, allocatable :: sendcounts(:), senddispls(:), recvcounts(:), recvdispls(:)
    integer :: dest, source, i, j, k, offset, ierr
    integer :: k_first, k_last, local_start(3), local_end(3)

    allocate(sendcounts(plan%nprocs), senddispls(plan%nprocs))
    allocate(recvcounts(plan%nprocs), recvdispls(plan%nprocs))
    local_start = starts(:,plan%rank + 1)
    local_end = ends(:,plan%rank + 1)

    do dest = 1, plan%nprocs
      k_first = max(plan%k_start, starts(3,dest))
      k_last = min(plan%k_start + plan%local_nz - 1, ends(3,dest))
      sendcounts(dest) = range_size(starts(1,dest), ends(1,dest)) * &
        range_size(starts(2,dest), ends(2,dest)) * range_size(k_first, k_last)
    end do
    do source = 1, plan%nprocs
      k_first = max(plan%slab_starts(source), local_start(3))
      k_last = min(plan%slab_starts(source) + plan%slab_counts(source) - 1, local_end(3))
      recvcounts(source) = range_size(local_start(1), local_end(1)) * &
        range_size(local_start(2), local_end(2)) * range_size(k_first, k_last)
    end do
    call make_displacements(sendcounts, senddispls)
    call make_displacements(recvcounts, recvdispls)
    allocate(sendbuf(sum(sendcounts)), recvbuf(sum(recvcounts)))

    !$omp parallel do schedule(static) if(gp3d_openmp_active) &
    !$omp& private(offset, k_first, k_last, i, j, k)
    do dest = 1, plan%nprocs
      offset = senddispls(dest) + 1
      k_first = max(plan%k_start, starts(3,dest))
      k_last = min(plan%k_start + plan%local_nz - 1, ends(3,dest))
      do k = k_first, k_last
        do j = starts(2,dest), ends(2,dest)
          do i = starts(1,dest), ends(1,dest)
            sendbuf(offset) = slab(i,j,k - plan%k_start + 1)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call MPI_Alltoallv(sendbuf, sendcounts, senddispls, complex_type, &
      recvbuf, recvcounts, recvdispls, complex_type, plan%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Alltoallv slab-to-pencil redistribution failed"

    !$omp parallel do schedule(static) if(gp3d_openmp_active) &
    !$omp& private(offset, k_first, k_last, i, j, k)
    do source = 1, plan%nprocs
      offset = recvdispls(source) + 1
      k_first = max(plan%slab_starts(source), local_start(3))
      k_last = min(plan%slab_starts(source) + plan%slab_counts(source) - 1, local_end(3))
      do k = k_first, k_last
        do j = local_start(2), local_end(2)
          do i = local_start(1), local_end(1)
            pencil(i - local_start(1) + 1, j - local_start(2) + 1, &
              k - local_start(3) + 1) = recvbuf(offset)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(sendbuf, recvbuf, sendcounts, senddispls, recvcounts, recvdispls)
  end subroutine slab_to_pencil

  subroutine pencil_to_slab(plan, pencil, slab, starts, ends)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(mytype), intent(in) :: pencil(:,:,:)
    complex(dp), intent(out) :: slab(plan%nx, plan%ny, plan%local_nz)
    integer, intent(in) :: starts(3,plan%nprocs), ends(3,plan%nprocs)
    complex(dp), allocatable :: sendbuf(:), recvbuf(:)
    integer, allocatable :: sendcounts(:), senddispls(:), recvcounts(:), recvdispls(:)
    integer :: dest, source, i, j, k, offset, ierr
    integer :: k_first, k_last, local_start(3), local_end(3)

    allocate(sendcounts(plan%nprocs), senddispls(plan%nprocs))
    allocate(recvcounts(plan%nprocs), recvdispls(plan%nprocs))
    local_start = starts(:,plan%rank + 1)
    local_end = ends(:,plan%rank + 1)

    do dest = 1, plan%nprocs
      k_first = max(plan%slab_starts(dest), local_start(3))
      k_last = min(plan%slab_starts(dest) + plan%slab_counts(dest) - 1, local_end(3))
      sendcounts(dest) = range_size(local_start(1), local_end(1)) * &
        range_size(local_start(2), local_end(2)) * range_size(k_first, k_last)
    end do
    do source = 1, plan%nprocs
      k_first = max(plan%k_start, starts(3,source))
      k_last = min(plan%k_start + plan%local_nz - 1, ends(3,source))
      recvcounts(source) = range_size(starts(1,source), ends(1,source)) * &
        range_size(starts(2,source), ends(2,source)) * range_size(k_first, k_last)
    end do
    call make_displacements(sendcounts, senddispls)
    call make_displacements(recvcounts, recvdispls)
    allocate(sendbuf(sum(sendcounts)), recvbuf(sum(recvcounts)))

    !$omp parallel do schedule(static) if(gp3d_openmp_active) &
    !$omp& private(offset, k_first, k_last, i, j, k)
    do dest = 1, plan%nprocs
      offset = senddispls(dest) + 1
      k_first = max(plan%slab_starts(dest), local_start(3))
      k_last = min(plan%slab_starts(dest) + plan%slab_counts(dest) - 1, local_end(3))
      do k = k_first, k_last
        do j = local_start(2), local_end(2)
          do i = local_start(1), local_end(1)
            sendbuf(offset) = pencil(i - local_start(1) + 1, &
              j - local_start(2) + 1, k - local_start(3) + 1)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call MPI_Alltoallv(sendbuf, sendcounts, senddispls, complex_type, &
      recvbuf, recvcounts, recvdispls, complex_type, plan%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Alltoallv pencil-to-slab redistribution failed"

    !$omp parallel do schedule(static) if(gp3d_openmp_active) &
    !$omp& private(offset, k_first, k_last, i, j, k)
    do source = 1, plan%nprocs
      offset = recvdispls(source) + 1
      k_first = max(plan%k_start, starts(3,source))
      k_last = min(plan%k_start + plan%local_nz - 1, ends(3,source))
      do k = k_first, k_last
        do j = starts(2,source), ends(2,source)
          do i = starts(1,source), ends(1,source)
            slab(i,j,k - plan%k_start + 1) = recvbuf(offset)
            offset = offset + 1
          end do
        end do
      end do
    end do
    !$omp end parallel do
    deallocate(sendbuf, recvbuf, sendcounts, senddispls, recvcounts, recvdispls)
  end subroutine pencil_to_slab

  pure integer function range_size(first, last) result(count)
    integer, intent(in) :: first, last

    count = max(0, last - first + 1)
  end function range_size

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
