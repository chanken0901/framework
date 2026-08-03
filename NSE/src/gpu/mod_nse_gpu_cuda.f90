module mod_nse_gpu
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated, &
    c_int, c_double, c_char, c_null_char
  use, intrinsic :: iso_fortran_env, only : error_unit
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_nse_forcing_common, only : forcing_is_enabled, &
    validate_forcing_parameters
  implicit none
  private

  type, public :: nse_gpu_context
    private
    type(c_ptr) :: handle = c_null_ptr
  end type nse_gpu_context

  public :: nse_gpu_initialize
  public :: nse_gpu_upload
  public :: nse_gpu_download
  public :: nse_gpu_compute_dt
  public :: nse_gpu_advance_ssprk3
  public :: nse_gpu_synchronize
  public :: nse_gpu_finalize
  public :: validate_nse_gpu_configuration

  interface
    function c_nse_cuda_create(handle, nx, ny, nz, nghost, nvar, device, &
        keep_order, viscous_enabled, forcing_enabled, forcing_spectrum, &
        forcing_report_interval, gamma, cfl, small_rho, small_p, reynolds, &
        prandtl, dx, dy, dz, forcing_k_cutoff, forcing_target_dissipation, &
        forcing_dilatational_ratio, forcing_denominator_floor, &
        forcing_max_coefficient) &
        bind(C, name="nse_cuda_create") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), intent(out) :: handle
      integer(c_int), value :: nx, ny, nz, nghost, nvar, device
      integer(c_int), value :: keep_order, viscous_enabled
      integer(c_int), value :: forcing_enabled, forcing_spectrum
      integer(c_int), value :: forcing_report_interval
      real(c_double), value :: gamma, cfl, small_rho, small_p
      real(c_double), value :: reynolds, prandtl, dx, dy, dz
      real(c_double), value :: forcing_k_cutoff, forcing_target_dissipation
      real(c_double), value :: forcing_dilatational_ratio
      real(c_double), value :: forcing_denominator_floor
      real(c_double), value :: forcing_max_coefficient
      integer(c_int) :: status
    end function c_nse_cuda_create

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
    integer :: selected_order

    if (nse%nv /= 5) error stop "CUDA KEEP backend requires five variables"
    if (sim%nghost < 3) then
      error stop "selectable-order CUDA KEEP backend requires three ghost cells"
    end if
    selected_order = requested_keep_order(nse)
    if (selected_order /= 2 .and. selected_order /= 6) then
      error stop "CUDA KEEP scheme must be keep2 or keep6"
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
    if (trim(adjustl(nse%boundary_condition)) /= "periodic") then
      error stop "CUDA backend currently supports periodic boundaries"
    end if
    if (trim(adjustl(nse%time_integrator)) /= "ssprk3") then
      error stop "CUDA backend currently supports time_integrator=ssprk3"
    end if
    if (forcing_is_enabled(nse)) then
      call validate_forcing_parameters(nse, 'cufft')
    end if
  end subroutine validate_nse_gpu_configuration

  subroutine nse_gpu_initialize(context, sim, nse)
    type(nse_gpu_context), intent(inout) :: context
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer(c_int) :: status, viscous_enabled, keep_order
    integer(c_int) :: forcing_enabled, forcing_spectrum

    call validate_nse_gpu_configuration(sim, nse)
    if (kind(1.0_dp) /= c_double) then
      error stop "CUDA backend requires double-precision dp"
    end if
    if (c_associated(context%handle)) then
      error stop "NSE CUDA context is already initialized"
    end if

    viscous_enabled = 0_c_int
    if (trim(adjustl(nse%viscous_scheme)) == "central6") then
      viscous_enabled = 1_c_int
    end if
    keep_order = int(requested_keep_order(nse), c_int)
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
      int(sim%ny, c_int), int(sim%nz, c_int), int(sim%nghost, c_int), &
      int(nse%nv, c_int), int(sim%cuda_device, c_int), &
      keep_order, viscous_enabled, forcing_enabled, forcing_spectrum, &
      int(nse%forcing_report_interval, c_int), &
      nse%gamma, nse%cfl, nse%small_rho, nse%small_p, nse%reynolds, &
      nse%prandtl, sim%dx, sim%dy, sim%dz, nse%forcing_k_cutoff, &
      nse%forcing_target_dissipation, nse%forcing_dilatational_ratio, &
      nse%forcing_denominator_floor, nse%forcing_max_coefficient)
    call require_success(status, "initialize NSE CUDA context")
  end subroutine nse_gpu_initialize

  pure integer function requested_keep_order(nse) result(order)
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(nse%convective_scheme)))
    case ("keep2")
      order = 2
    case ("keep6")
      order = 6
    case default
      order = 0
    end select
  end function requested_keep_order

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
