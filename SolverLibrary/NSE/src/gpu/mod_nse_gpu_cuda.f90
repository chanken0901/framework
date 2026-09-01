module mod_nse_gpu
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated, &
    c_int, c_double, c_char, c_null_char, c_size_t
  use, intrinsic :: iso_fortran_env, only : error_unit
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config, nse_boundary_face_count
  use mod_nse_boundary, only : validate_boundary_scheme
  use mod_nse_forcing_common, only : forcing_is_enabled, &
    validate_forcing_parameters
  implicit none
  private

  integer(c_int), parameter :: cuda_convective_keep2 = 1_c_int
  integer(c_int), parameter :: cuda_convective_keep6 = 2_c_int
  integer(c_int), parameter :: cuda_convective_weno5z_roe = 3_c_int
  integer(c_int), parameter :: cuda_convective_hybrid = 4_c_int
  integer(c_int), parameter :: cuda_sensor_ducros_pressure = 1_c_int
  integer(c_int), parameter :: cuda_boundary_periodic = 0_c_int
  integer(c_int), parameter :: cuda_boundary_non_reflecting = 1_c_int
  integer(c_int), parameter :: cuda_boundary_reflective = 2_c_int
  integer(c_int), parameter :: cuda_boundary_dirichlet = 3_c_int
  integer, parameter, public :: nse_gpu_halo_y = 1
  integer, parameter, public :: nse_gpu_halo_z = 2
  integer, parameter, public :: nse_gpu_halo_low = -1
  integer, parameter, public :: nse_gpu_halo_high = 1

  type, public :: nse_gpu_context
    private
    type(c_ptr) :: handle = c_null_ptr
  end type nse_gpu_context

  public :: nse_gpu_initialize
  public :: nse_gpu_select_device
  public :: nse_gpu_configure_cufftmp
  public :: nse_gpu_upload
  public :: nse_gpu_download
  public :: nse_gpu_compute_dt
  public :: nse_gpu_advance_ssprk3
  public :: nse_gpu_begin_ssprk3
  public :: nse_gpu_advance_ssprk3_stage
  public :: nse_gpu_apply_local_boundary
  public :: nse_gpu_apply_local_periodic
  public :: nse_gpu_halo_count
  public :: nse_gpu_pack_halo
  public :: nse_gpu_unpack_halo
  public :: nse_gpu_synchronize
  public :: nse_gpu_finalize
  public :: validate_nse_gpu_configuration

  interface
    function c_nse_cuda_set_device(device) &
        bind(C, name="nse_cuda_set_device") result(status)
      import :: c_int
      integer(c_int), value :: device
      integer(c_int) :: status
    end function c_nse_cuda_set_device

    function c_nse_cuda_create(handle, nx, ny, nz, nghost, nvar, device, &
        distributed_y, distributed_z, &
        global_ny, global_nz, global_y_start, global_z_start, &
        boundary_type, boundary_reference, boundary_relaxation_strength, &
        boundary_length_scale, &
        convective_scheme, hybrid_smooth_scheme, hybrid_shock_scheme, &
        hybrid_sensor, viscous_enabled, forcing_enabled, forcing_spectrum, &
        forcing_report_interval, gamma, cfl, small_rho, small_p, reynolds, &
        prandtl, dx, dy, dz, forcing_k_cutoff, forcing_target_dissipation, &
        forcing_dilatational_ratio, forcing_denominator_floor, &
        forcing_max_coefficient, hybrid_sensor_onset, hybrid_sensor_full) &
        bind(C, name="nse_cuda_create") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), intent(out) :: handle
      integer(c_int), value :: nx, ny, nz, nghost, nvar, device
      integer(c_int), value :: distributed_y, distributed_z
      integer(c_int), value :: global_ny, global_nz
      integer(c_int), value :: global_y_start, global_z_start
      integer(c_int), intent(in) :: boundary_type(6)
      real(c_double), intent(in) :: boundary_reference(5,6)
      real(c_double), value :: boundary_relaxation_strength
      real(c_double), value :: boundary_length_scale
      integer(c_int), value :: convective_scheme
      integer(c_int), value :: hybrid_smooth_scheme, hybrid_shock_scheme
      integer(c_int), value :: hybrid_sensor, viscous_enabled
      integer(c_int), value :: forcing_enabled, forcing_spectrum
      integer(c_int), value :: forcing_report_interval
      real(c_double), value :: gamma, cfl, small_rho, small_p
      real(c_double), value :: reynolds, prandtl, dx, dy, dz
      real(c_double), value :: forcing_k_cutoff, forcing_target_dissipation
      real(c_double), value :: forcing_dilatational_ratio
      real(c_double), value :: forcing_denominator_floor
      real(c_double), value :: forcing_max_coefficient
      real(c_double), value :: hybrid_sensor_onset, hybrid_sensor_full
      integer(c_int) :: status
    end function c_nse_cuda_create

    function c_nse_cuda_configure_cufftmp(handle, global_ny, global_nz, &
        global_y_start, global_z_start, communicator) &
        bind(C, name="nse_cuda_configure_cufftmp") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int), value :: global_ny, global_nz
      integer(c_int), value :: global_y_start, global_z_start
      integer(c_int), value :: communicator
      integer(c_int) :: status
    end function c_nse_cuda_configure_cufftmp

    function c_nse_cuda_upload(handle, q) &
        bind(C, name="nse_cuda_upload") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: q(*)
      integer(c_int) :: status
    end function c_nse_cuda_upload

    function c_nse_cuda_download(handle, q) &
        bind(C, name="nse_cuda_download") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: q(*)
      integer(c_int) :: status
    end function c_nse_cuda_download

    function c_nse_cuda_apply_local_boundary(handle) &
        bind(C, name="nse_cuda_apply_local_boundary") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: status
    end function c_nse_cuda_apply_local_boundary

    function c_nse_cuda_halo_count(handle, direction, count) &
        bind(C, name="nse_cuda_halo_count") result(status)
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: handle
      integer(c_int), value :: direction
      integer(c_size_t), intent(out) :: count
      integer(c_int) :: status
    end function c_nse_cuda_halo_count

    function c_nse_cuda_pack_halo(handle, direction, side, buffer) &
        bind(C, name="nse_cuda_pack_halo") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      integer(c_int), value :: direction, side
      real(c_double), intent(out) :: buffer(*)
      integer(c_int) :: status
    end function c_nse_cuda_pack_halo

    function c_nse_cuda_unpack_halo(handle, direction, side, buffer) &
        bind(C, name="nse_cuda_unpack_halo") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      integer(c_int), value :: direction, side
      real(c_double), intent(in) :: buffer(*)
      integer(c_int) :: status
    end function c_nse_cuda_unpack_halo

    function c_nse_cuda_compute_dt(handle, dt) &
        bind(C, name="nse_cuda_compute_dt") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: dt
      integer(c_int) :: status
    end function c_nse_cuda_compute_dt

    function c_nse_cuda_advance(handle, dt) &
        bind(C, name="nse_cuda_advance_ssprk3") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: dt
      integer(c_int) :: status
    end function c_nse_cuda_advance

    function c_nse_cuda_begin_ssprk3(handle, dt) &
        bind(C, name="nse_cuda_begin_ssprk3") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: dt
      integer(c_int) :: status
    end function c_nse_cuda_begin_ssprk3

    function c_nse_cuda_advance_stage(handle, dt, stage) &
        bind(C, name="nse_cuda_advance_ssprk3_stage") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: dt
      integer(c_int), value :: stage
      integer(c_int) :: status
    end function c_nse_cuda_advance_stage

    function c_nse_cuda_synchronize(handle) &
        bind(C, name="nse_cuda_synchronize") result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: status
    end function c_nse_cuda_synchronize

    subroutine c_nse_cuda_destroy(handle) bind(C, name="nse_cuda_destroy")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine c_nse_cuda_destroy

    subroutine c_nse_cuda_get_last_error(buffer, buffer_size) &
        bind(C, name="nse_cuda_get_last_error")
      import :: c_char, c_int
      character(c_char), intent(out) :: buffer(*)
      integer(c_int), value :: buffer_size
    end subroutine c_nse_cuda_get_last_error
  end interface

contains

  subroutine validate_nse_gpu_configuration(sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer(c_int) :: selected_scheme

    if (nse%nv /= 5) error stop "CUDA backend requires five variables"
    if (sim%nghost < 3) then
      error stop "CUDA convective schemes require three ghost cells"
    end if
    selected_scheme = requested_cuda_convective_scheme(nse)
    if (selected_scheme == 0_c_int) then
      error stop &
        "CUDA convective scheme must be keep2, keep6, weno5z_roe, or hybrid"
    end if
    if (selected_scheme == cuda_convective_hybrid) then
      if (requested_cuda_leaf_scheme(nse%hybrid_smooth_scheme) == 0_c_int) &
        error stop "invalid CUDA hybrid_smooth_scheme"
      if (requested_cuda_leaf_scheme(nse%hybrid_shock_scheme) == 0_c_int) &
        error stop "invalid CUDA hybrid_shock_scheme"
      if (requested_cuda_sensor(nse%hybrid_sensor) == 0_c_int) &
        error stop "CUDA hybrid_sensor must be ducros_pressure"
      if (nse%hybrid_sensor_onset < 0.0_dp .or. &
          nse%hybrid_sensor_full <= nse%hybrid_sensor_onset) then
        error stop "CUDA hybrid sensor requires 0 <= onset < full"
      end if
    end if
    if (trim(adjustl(nse%viscous_scheme)) /= "none" .and. &
        trim(adjustl(nse%viscous_scheme)) /= "central6") then
      error stop "CUDA backend supports viscous_scheme=none or central6"
    end if
    if (trim(adjustl(nse%viscous_scheme)) == "central6") then
      if (sim%nghost < 3) then
        error stop "CUDA central6 viscosity requires at least three ghost cells"
      end if
      if (nse%reynolds <= 0.0_dp) then
        error stop "CUDA central6 viscosity requires reynolds > 0"
      end if
      if (nse%prandtl <= 0.0_dp) then
        error stop "CUDA central6 viscosity requires prandtl > 0"
      end if
    end if
    call validate_boundary_scheme(sim, nse)
    if (trim(adjustl(nse%time_integrator)) /= "ssprk3") then
      error stop "CUDA backend currently supports time_integrator=ssprk3"
    end if
    if (forcing_is_enabled(nse)) then
      if (sim%use_mpi) then
        call validate_forcing_parameters(nse, 'cufftmp')
      else
        call validate_forcing_parameters(nse, 'cufft')
      end if
    end if
  end subroutine validate_nse_gpu_configuration

  subroutine nse_gpu_select_device(device)
    integer, intent(in) :: device
    integer(c_int) :: status

    status = c_nse_cuda_set_device(int(device,c_int))
    call require_success(status, "select CUDA device")
  end subroutine nse_gpu_select_device

  subroutine nse_gpu_initialize(context, sim, nse, local_ny, local_nz, &
      distributed_y, distributed_z, device, global_y_start, global_z_start)
    type(nse_gpu_context), intent(inout) :: context
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in), optional :: local_ny, local_nz, device
    integer, intent(in), optional :: global_y_start, global_z_start
    logical, intent(in), optional :: distributed_y, distributed_z
    integer(c_int) :: status, viscous_enabled, convective_scheme
    integer(c_int) :: forcing_enabled, forcing_spectrum
    integer(c_int) :: hybrid_smooth_scheme, hybrid_shock_scheme, hybrid_sensor
    integer(c_int) :: cuda_distributed_y, cuda_distributed_z
    integer(c_int) :: boundary_type(nse_boundary_face_count)
    real(c_double) :: boundary_reference(5,nse_boundary_face_count)
    integer :: cuda_ny, cuda_nz, cuda_device
    integer :: cuda_global_y_start, cuda_global_z_start, face

    call validate_nse_gpu_configuration(sim, nse)
    if (kind(1.0_dp) /= c_double) then
      error stop "CUDA backend requires double-precision dp"
    end if
    if (c_associated(context%handle)) then
      error stop "NSE CUDA context is already initialized"
    end if

    cuda_ny = sim%ny
    cuda_nz = sim%nz
    cuda_device = sim%cuda_device
    cuda_global_y_start = 0
    cuda_global_z_start = 0
    if (present(local_ny)) cuda_ny = local_ny
    if (present(local_nz)) cuda_nz = local_nz
    if (present(device)) cuda_device = device
    if (present(global_y_start)) cuda_global_y_start = global_y_start
    if (present(global_z_start)) cuda_global_z_start = global_z_start
    if (cuda_ny <= 0 .or. cuda_nz <= 0) then
      error stop "local CUDA domain dimensions must be positive"
    end if
    if (cuda_global_y_start < 0 .or. cuda_global_z_start < 0 .or. &
        cuda_global_y_start+cuda_ny > sim%ny .or. &
        cuda_global_z_start+cuda_nz > sim%nz) then
      error stop "local CUDA domain lies outside the global Y-Z domain"
    end if
    cuda_distributed_y = 0_c_int
    cuda_distributed_z = 0_c_int
    if (present(distributed_y)) then
      if (distributed_y) cuda_distributed_y = 1_c_int
    end if
    if (present(distributed_z)) then
      if (distributed_z) cuda_distributed_z = 1_c_int
    end if

    boundary_reference = 0.0_c_double
    do face = 1, nse_boundary_face_count
      select case (trim(adjustl(nse%boundary_face_type(face))))
      case ('periodic')
        boundary_type(face) = cuda_boundary_periodic
      case ('non_reflecting')
        boundary_type(face) = cuda_boundary_non_reflecting
      case ('reflective')
        boundary_type(face) = cuda_boundary_reflective
      case ('dirichlet')
        boundary_type(face) = cuda_boundary_dirichlet
      case default
        error stop "unsupported CUDA boundary face type"
      end select
      boundary_reference(1,face) = nse%boundary_reference_rho(face)
      boundary_reference(2:4,face) = &
        nse%boundary_reference_velocity(:,face)
      boundary_reference(5,face) = nse%boundary_reference_p(face)
    end do

    viscous_enabled = 0_c_int
    if (trim(adjustl(nse%viscous_scheme)) == "central6") then
      viscous_enabled = 1_c_int
    end if
    convective_scheme = requested_cuda_convective_scheme(nse)
    hybrid_smooth_scheme = requested_cuda_leaf_scheme( &
      nse%hybrid_smooth_scheme)
    hybrid_shock_scheme = requested_cuda_leaf_scheme( &
      nse%hybrid_shock_scheme)
    hybrid_sensor = requested_cuda_sensor(nse%hybrid_sensor)
    forcing_enabled = 0_c_int
    forcing_spectrum = 0_c_int
    if (forcing_is_enabled(nse)) then
      forcing_enabled = 1_c_int
      select case (trim(adjustl(nse%forcing_spectrum)))
      case ('full_spectrum')
        forcing_spectrum = 1_c_int
      case ('low_wavenumber')
        forcing_spectrum = 2_c_int
      end select
    end if

    status = c_nse_cuda_create(context%handle, int(sim%nx, c_int), &
      int(cuda_ny, c_int), int(cuda_nz, c_int), int(sim%nghost, c_int), &
      int(nse%nv, c_int), int(cuda_device, c_int), &
      cuda_distributed_y, cuda_distributed_z, &
      int(sim%ny,c_int), int(sim%nz,c_int), &
      int(cuda_global_y_start,c_int), int(cuda_global_z_start,c_int), &
      boundary_type, boundary_reference, nse%boundary_relaxation_strength, &
      nse%boundary_length_scale, &
      convective_scheme, hybrid_smooth_scheme, hybrid_shock_scheme, &
      hybrid_sensor, viscous_enabled, forcing_enabled, forcing_spectrum, &
      int(nse%forcing_report_interval, c_int), &
      nse%gamma, nse%cfl, nse%small_rho, nse%small_p, nse%reynolds, &
      nse%prandtl, sim%dx, sim%dy, sim%dz, nse%forcing_k_cutoff, &
      nse%forcing_target_dissipation, nse%forcing_dilatational_ratio, &
      nse%forcing_denominator_floor, nse%forcing_max_coefficient, &
      nse%hybrid_sensor_onset, nse%hybrid_sensor_full)
    call require_success(status, "initialize NSE CUDA context")
  end subroutine nse_gpu_initialize

  subroutine nse_gpu_configure_cufftmp(context, global_ny, global_nz, &
      global_y_start, global_z_start, communicator)
    type(nse_gpu_context), intent(in) :: context
    integer, intent(in) :: global_ny, global_nz
    integer, intent(in) :: global_y_start, global_z_start
    integer, intent(in) :: communicator
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_configure_cufftmp(context%handle, &
      int(global_ny,c_int), int(global_nz,c_int), &
      int(global_y_start,c_int), int(global_z_start,c_int), &
      int(communicator,c_int))
    call require_success(status, "configure distributed cuFFTMp forcing")
  end subroutine nse_gpu_configure_cufftmp

  pure integer(c_int) function requested_cuda_convective_scheme(nse) &
      result(scheme)
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(nse%convective_scheme)))
    case ("keep2")
      scheme = cuda_convective_keep2
    case ("keep6")
      scheme = cuda_convective_keep6
    case ("weno5z_roe")
      scheme = cuda_convective_weno5z_roe
    case ("hybrid")
      scheme = cuda_convective_hybrid
    case default
      scheme = 0_c_int
    end select
  end function requested_cuda_convective_scheme

  pure integer(c_int) function requested_cuda_leaf_scheme(name) &
      result(scheme)
    character(len=*), intent(in) :: name

    select case (trim(adjustl(name)))
    case ("keep2")
      scheme = cuda_convective_keep2
    case ("keep6")
      scheme = cuda_convective_keep6
    case ("weno5z_roe")
      scheme = cuda_convective_weno5z_roe
    case default
      scheme = 0_c_int
    end select
  end function requested_cuda_leaf_scheme

  pure integer(c_int) function requested_cuda_sensor(name) result(sensor)
    character(len=*), intent(in) :: name

    select case (trim(adjustl(name)))
    case ("ducros_pressure")
      sensor = cuda_sensor_ducros_pressure
    case default
      sensor = 0_c_int
    end select
  end function requested_cuda_sensor

  subroutine nse_gpu_upload(context, q)
    type(nse_gpu_context), intent(in) :: context
    real(dp), contiguous, intent(in) :: q(:,:,:,:)
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_upload(context%handle, q)
    call require_success(status, "upload NSE state")
  end subroutine nse_gpu_upload

  subroutine nse_gpu_download(context, q)
    type(nse_gpu_context), intent(in) :: context
    real(dp), contiguous, intent(out) :: q(:,:,:,:)
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_download(context%handle, q)
    call require_success(status, "download NSE state")
  end subroutine nse_gpu_download

  subroutine nse_gpu_apply_local_boundary(context)
    type(nse_gpu_context), intent(in) :: context
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_apply_local_boundary(context%handle)
    call require_success(status, "apply local CUDA boundaries")
  end subroutine nse_gpu_apply_local_boundary

  subroutine nse_gpu_apply_local_periodic(context)
    type(nse_gpu_context), intent(in) :: context

    ! Backward-compatible API name. The runtime kernel now applies the
    ! configured periodic/non-reflecting boundary on each physical face.
    call nse_gpu_apply_local_boundary(context)
  end subroutine nse_gpu_apply_local_periodic

  integer function nse_gpu_halo_count(context, direction) result(count)
    type(nse_gpu_context), intent(in) :: context
    integer, intent(in) :: direction
    integer(c_size_t) :: c_count
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_halo_count(context%handle, &
      int(direction, c_int), c_count)
    call require_success(status, "query CUDA MPI halo size")
    if (c_count > int(huge(count), c_size_t)) then
      error stop "CUDA MPI halo exceeds the MPI default-integer count limit"
    end if
    count = int(c_count)
  end function nse_gpu_halo_count

  subroutine nse_gpu_pack_halo(context, direction, side, buffer)
    type(nse_gpu_context), intent(in) :: context
    integer, intent(in) :: direction, side
    real(dp), contiguous, intent(out) :: buffer(:)
    integer(c_int) :: status

    call require_context(context)
    if (size(buffer) < nse_gpu_halo_count(context, direction)) then
      error stop "host send buffer is too small for the CUDA MPI halo"
    end if
    status = c_nse_cuda_pack_halo(context%handle, int(direction, c_int), &
      int(side, c_int), buffer)
    call require_success(status, "pack CUDA MPI halo")
  end subroutine nse_gpu_pack_halo

  subroutine nse_gpu_unpack_halo(context, direction, side, buffer)
    type(nse_gpu_context), intent(in) :: context
    integer, intent(in) :: direction, side
    real(dp), contiguous, intent(in) :: buffer(:)
    integer(c_int) :: status

    call require_context(context)
    if (size(buffer) < nse_gpu_halo_count(context, direction)) then
      error stop "host receive buffer is too small for the CUDA MPI halo"
    end if
    status = c_nse_cuda_unpack_halo(context%handle, &
      int(direction, c_int), int(side, c_int), buffer)
    call require_success(status, "unpack CUDA MPI halo")
  end subroutine nse_gpu_unpack_halo

  subroutine nse_gpu_compute_dt(context, dt)
    type(nse_gpu_context), intent(in) :: context
    real(dp), intent(out) :: dt
    real(c_double) :: c_dt
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_compute_dt(context%handle, c_dt)
    call require_success(status, "compute CUDA CFL time step")
    dt = c_dt
  end subroutine nse_gpu_compute_dt

  subroutine nse_gpu_advance_ssprk3(context, dt)
    type(nse_gpu_context), intent(in) :: context
    real(dp), intent(in) :: dt
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_advance(context%handle, dt)
    call require_success(status, "advance CUDA SSPRK3 step")
  end subroutine nse_gpu_advance_ssprk3

  subroutine nse_gpu_begin_ssprk3(context, dt)
    type(nse_gpu_context), intent(in) :: context
    real(dp), intent(in) :: dt
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_begin_ssprk3(context%handle, dt)
    call require_success(status, "begin distributed CUDA SSPRK3 step")
  end subroutine nse_gpu_begin_ssprk3

  subroutine nse_gpu_advance_ssprk3_stage(context, dt, stage)
    type(nse_gpu_context), intent(in) :: context
    real(dp), intent(in) :: dt
    integer, intent(in) :: stage
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_advance_stage(context%handle, dt, int(stage, c_int))
    call require_success(status, "advance distributed CUDA SSPRK3 stage")
  end subroutine nse_gpu_advance_ssprk3_stage

  subroutine nse_gpu_synchronize(context)
    type(nse_gpu_context), intent(in) :: context
    integer(c_int) :: status

    call require_context(context)
    status = c_nse_cuda_synchronize(context%handle)
    call require_success(status, "synchronize NSE CUDA backend")
  end subroutine nse_gpu_synchronize

  subroutine nse_gpu_finalize(context)
    type(nse_gpu_context), intent(inout) :: context

    if (c_associated(context%handle)) call c_nse_cuda_destroy(context%handle)
    context%handle = c_null_ptr
  end subroutine nse_gpu_finalize

  subroutine require_context(context)
    type(nse_gpu_context), intent(in) :: context

    if (.not. c_associated(context%handle)) then
      error stop "NSE CUDA context is not initialized"
    end if
  end subroutine require_context

  subroutine require_success(status, operation)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: operation
    character(c_char) :: buffer(512)
    character(len=512) :: message
    integer :: i

    if (status == 0_c_int) return
    buffer = c_null_char
    call c_nse_cuda_get_last_error(buffer, int(size(buffer), c_int))
    message = ""
    do i = 1, size(buffer)
      if (buffer(i) == c_null_char) exit
      message(i:i) = buffer(i)
    end do
    write(error_unit,'(3a)') "CUDA error during ", trim(operation), &
      ": " // trim(message)
    error stop "NSE CUDA backend failure"
  end subroutine require_success

end module mod_nse_gpu
