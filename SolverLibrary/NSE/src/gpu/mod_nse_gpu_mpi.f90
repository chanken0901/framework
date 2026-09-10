module mod_nse_gpu_mpi
  use mod_precision, only : dp
  use mod_model_config, only : nse_config, nse_face_y_min, nse_face_z_min
  use mod_nse_gpu, only : nse_gpu_context, nse_gpu_apply_local_boundary, &
    nse_gpu_halo_count, nse_gpu_pack_halo, nse_gpu_unpack_halo, &
    nse_gpu_device_mpi_available, nse_gpu_exchange_device_halo, &
    nse_gpu_halo_y, nse_gpu_halo_z, nse_gpu_halo_low, nse_gpu_halo_high
  use module_mpi, only : ndiv_ny, ndiv_nz, j_myrank, k_myrank, itable, &
    mp_sendrecv_r8, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE, &
    MPI_PROC_NULL, mp_allmaxi, mp_stop, my_rank
  implicit none
  private

  integer, parameter :: tag_y_low = 7101
  integer, parameter :: tag_y_high = 7102
  integer, parameter :: tag_z_low = 7201
  integer, parameter :: tag_z_high = 7202

  type, public :: nse_gpu_mpi_halo
    private
    integer :: y_count = 0
    integer :: z_count = 0
    logical :: initialized = .false.
    logical :: device_transport = .false.
    logical :: periodic_y = .true.
    logical :: periodic_z = .true.
    real(dp), allocatable :: send_low(:)
    real(dp), allocatable :: send_high(:)
    real(dp), allocatable :: recv_low(:)
    real(dp), allocatable :: recv_high(:)
  end type nse_gpu_mpi_halo

  public :: nse_gpu_mpi_halo_initialize
  public :: nse_gpu_mpi_exchange
  public :: nse_gpu_mpi_halo_finalize

contains

  subroutine nse_gpu_mpi_halo_initialize(halo, context, nse)
    type(nse_gpu_mpi_halo), intent(inout) :: halo
    type(nse_gpu_context), intent(in) :: context
    type(nse_config), intent(in) :: nse
    integer :: max_count, env_status, mode, mode_max, mode_neg, unavailable
    character(len=32) :: transport

    if (halo%initialized) then
      error stop 'NSE MPI+CUDA halo workspace is already initialized'
    end if
    halo%y_count = 0
    halo%z_count = 0
    halo%periodic_y = &
      trim(adjustl(nse%boundary_face_type(nse_face_y_min))) == 'periodic'
    halo%periodic_z = &
      trim(adjustl(nse%boundary_face_type(nse_face_z_min))) == 'periodic'
    if (ndiv_ny > 1) halo%y_count = nse_gpu_halo_count(context, nse_gpu_halo_y)
    if (ndiv_nz > 1) halo%z_count = nse_gpu_halo_count(context, nse_gpu_halo_z)
    max_count = max(halo%y_count, halo%z_count)
    if (max_count <= 0) then
      error stop 'MPI+CUDA requires at least one distributed direction'
    end if
    transport='auto'
    call get_environment_variable('NSE_CUDA_MPI_TRANSPORT',transport,status=env_status)
    if(env_status==1) transport='auto'
    mode=-1
    if(env_status>=0 .or. env_status==1) then
      select case(trim(transport))
      case('auto'); mode=0
      case('staged'); mode=1
      case('device'); mode=2
      end select
    end if
    mode_max=mode; mode_neg=-mode
    call mp_allmaxi(mode_max)
    call mp_allmaxi(mode_neg)
    if(mode<0 .or. mode_max/=-mode_neg) then
      write(*,'(a)') 'NSE_CUDA_MPI_TRANSPORT must be auto/staged/device and identical on all ranks'
      call mp_stop(16)
    end if
    unavailable=0
    if(.not. nse_gpu_device_mpi_available()) unavailable=1
    call mp_allmaxi(unavailable)
    if(mode==2 .and. unavailable/=0) then
      write(*,'(a)') 'Device MPI requested but unavailable: enable CUDA-aware build and CUDA-aware Open MPI'
      call mp_stop(16)
    end if
    halo%device_transport=mode/=1 .and. unavailable==0
    if(halo%device_transport) then
      if(my_rank==0) write(*,'(a)') 'NSE CUDA MPI halo transport: device (CUDA-aware MPI)'
    else
      if(my_rank==0) write(*,'(a)') 'NSE CUDA MPI halo transport: staged (host buffers)'
      allocate(halo%send_low(max_count), halo%send_high(max_count))
      allocate(halo%recv_low(max_count), halo%recv_high(max_count))
    end if
    halo%initialized=.true.
  end subroutine nse_gpu_mpi_halo_initialize

  subroutine nse_gpu_mpi_exchange(halo, context)
    type(nse_gpu_mpi_halo), intent(inout) :: halo
    type(nse_gpu_context), intent(in) :: context

    if (.not. halo%initialized) then
      error stop 'NSE MPI+CUDA halo workspace is not initialized'
    end if

    ! Physical periodic/non-reflecting/reflective faces and every
    ! non-decomposed direction are completed on the GPU. Y is exchanged before
    ! Z so the Z slabs carry valid X-Y edges and Y-Z corners.
    call nse_gpu_apply_local_boundary(context)
    if (ndiv_ny > 1) call exchange_y(halo, context)
    if (ndiv_nz > 1) call exchange_z(halo, context)
  end subroutine nse_gpu_mpi_exchange

  subroutine exchange_y(halo, context)
    type(nse_gpu_mpi_halo), intent(inout) :: halo
    type(nse_gpu_context), intent(in) :: context
    integer :: rank_low, rank_high, ierr
    integer :: status(MPI_STATUS_SIZE)

    if (halo%periodic_y) then
      rank_low = itable(modulo(j_myrank-1, ndiv_ny), k_myrank)
      rank_high = itable(modulo(j_myrank+1, ndiv_ny), k_myrank)
    else
      rank_low = itable(j_myrank-1, k_myrank)
      rank_high = itable(j_myrank+1, k_myrank)
    end if
    if(halo%device_transport) then
      call nse_gpu_exchange_device_halo(context,nse_gpu_halo_y, &
        merge(-1,rank_low,rank_low==MPI_PROC_NULL),merge(-1,rank_high,rank_high==MPI_PROC_NULL), &
        MPI_COMM_WORLD,tag_y_low)
      return
    end if
    call nse_gpu_pack_halo(context, nse_gpu_halo_y, nse_gpu_halo_low, &
      halo%send_low(1:halo%y_count))
    call nse_gpu_pack_halo(context, nse_gpu_halo_y, nse_gpu_halo_high, &
      halo%send_high(1:halo%y_count))

    call mp_sendrecv_r8(halo%send_low, halo%y_count, MPI_DOUBLE_PRECISION, &
      rank_low, tag_y_low, halo%recv_high, halo%y_count, &
      MPI_DOUBLE_PRECISION, rank_high, tag_y_low, MPI_COMM_WORLD, status, ierr)
    if (ierr /= 0) error stop 'MPI+CUDA Y-low/Y-high halo exchange failed'
    call mp_sendrecv_r8(halo%send_high, halo%y_count, MPI_DOUBLE_PRECISION, &
      rank_high, tag_y_high, halo%recv_low, halo%y_count, &
      MPI_DOUBLE_PRECISION, rank_low, tag_y_high, MPI_COMM_WORLD, status, ierr)
    if (ierr /= 0) error stop 'MPI+CUDA Y-high/Y-low halo exchange failed'

    if (rank_high /= MPI_PROC_NULL) then
      call nse_gpu_unpack_halo(context, nse_gpu_halo_y, nse_gpu_halo_high, &
        halo%recv_high(1:halo%y_count))
    end if
    if (rank_low /= MPI_PROC_NULL) then
      call nse_gpu_unpack_halo(context, nse_gpu_halo_y, nse_gpu_halo_low, &
        halo%recv_low(1:halo%y_count))
    end if
  end subroutine exchange_y

  subroutine exchange_z(halo, context)
    type(nse_gpu_mpi_halo), intent(inout) :: halo
    type(nse_gpu_context), intent(in) :: context
    integer :: rank_low, rank_high, ierr
    integer :: status(MPI_STATUS_SIZE)

    if (halo%periodic_z) then
      rank_low = itable(j_myrank, modulo(k_myrank-1, ndiv_nz))
      rank_high = itable(j_myrank, modulo(k_myrank+1, ndiv_nz))
    else
      rank_low = itable(j_myrank, k_myrank-1)
      rank_high = itable(j_myrank, k_myrank+1)
    end if
    if(halo%device_transport) then
      call nse_gpu_exchange_device_halo(context,nse_gpu_halo_z, &
        merge(-1,rank_low,rank_low==MPI_PROC_NULL),merge(-1,rank_high,rank_high==MPI_PROC_NULL), &
        MPI_COMM_WORLD,tag_z_low)
      return
    end if
    call nse_gpu_pack_halo(context, nse_gpu_halo_z, nse_gpu_halo_low, &
      halo%send_low(1:halo%z_count))
    call nse_gpu_pack_halo(context, nse_gpu_halo_z, nse_gpu_halo_high, &
      halo%send_high(1:halo%z_count))

    call mp_sendrecv_r8(halo%send_low, halo%z_count, MPI_DOUBLE_PRECISION, &
      rank_low, tag_z_low, halo%recv_high, halo%z_count, &
      MPI_DOUBLE_PRECISION, rank_high, tag_z_low, MPI_COMM_WORLD, status, ierr)
    if (ierr /= 0) error stop 'MPI+CUDA Z-low/Z-high halo exchange failed'
    call mp_sendrecv_r8(halo%send_high, halo%z_count, MPI_DOUBLE_PRECISION, &
      rank_high, tag_z_high, halo%recv_low, halo%z_count, &
      MPI_DOUBLE_PRECISION, rank_low, tag_z_high, MPI_COMM_WORLD, status, ierr)
    if (ierr /= 0) error stop 'MPI+CUDA Z-high/Z-low halo exchange failed'

    if (rank_high /= MPI_PROC_NULL) then
      call nse_gpu_unpack_halo(context, nse_gpu_halo_z, nse_gpu_halo_high, &
        halo%recv_high(1:halo%z_count))
    end if
    if (rank_low /= MPI_PROC_NULL) then
      call nse_gpu_unpack_halo(context, nse_gpu_halo_z, nse_gpu_halo_low, &
        halo%recv_low(1:halo%z_count))
    end if
  end subroutine exchange_z

  subroutine nse_gpu_mpi_halo_finalize(halo)
    type(nse_gpu_mpi_halo), intent(inout) :: halo

    if (allocated(halo%send_low)) deallocate(halo%send_low)
    if (allocated(halo%send_high)) deallocate(halo%send_high)
    if (allocated(halo%recv_low)) deallocate(halo%recv_low)
    if (allocated(halo%recv_high)) deallocate(halo%recv_high)
    halo%y_count = 0
    halo%z_count = 0
    halo%periodic_y = .true.
    halo%periodic_z = .true.
    halo%initialized = .false.
    halo%device_transport = .false.
  end subroutine nse_gpu_mpi_halo_finalize

end module mod_nse_gpu_mpi
