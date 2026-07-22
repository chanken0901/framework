!> CUDA CブリッジをFortranから呼ぶためのISO_C_BINDINGラッパー。
!> GPUコンテキストの寿命、ホスト-GPU転送、Split-step、ARGLE、診断量を型安全に公開する。
module gp3d_gpu
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, &
    c_double, c_double_complex, c_char, c_null_char
  use, intrinsic :: iso_fortran_env, only: error_unit
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_params_t, gp3d_state_t, &
    gp3d_model_config_t
  use gp3d_solver, only: gp3d_step_timing_t
  implicit none
  private

  type, public :: gp3d_gpu_context_t
    private
    type(c_ptr) :: handle = c_null_ptr
  end type gp3d_gpu_context_t

  public :: gp3d_gpu_init
  public :: gp3d_gpu_upload
  public :: gp3d_gpu_download
  public :: gp3d_gpu_step
  public :: gp3d_gpu_relax_taylor_green_argle
  public :: gp3d_gpu_diagnostics
  public :: gp3d_gpu_finalize

  interface
    function c_gpu_create(handle, nx, ny, nz, need_argle, device) &
        bind(C, name="gp3d_cuda_create") result(status)
      import :: c_ptr, c_int
      type(c_ptr), intent(out) :: handle
      integer(c_int), value :: nx, ny, nz, need_argle, device
      integer(c_int) :: status
    end function c_gpu_create

    function c_gpu_upload(handle, psi, potential, kx, ky, kz) &
        bind(C, name="gp3d_cuda_upload") result(status)
      import :: c_ptr, c_int, c_double, c_double_complex
      type(c_ptr), value :: handle
      complex(c_double_complex), intent(in) :: psi(*)
      real(c_double), intent(in) :: potential(*), kx(*), ky(*), kz(*)
      integer(c_int) :: status
    end function c_gpu_upload

    function c_gpu_download(handle, psi) &
        bind(C, name="gp3d_cuda_download") result(status)
      import :: c_ptr, c_int, c_double_complex
      type(c_ptr), value :: handle
      complex(c_double_complex), intent(out) :: psi(*)
      integer(c_int) :: status
    end function c_gpu_download

    function c_gpu_step(handle, dt, hbar, mass, g, target_norm, &
        volume_element, imaginary_time, measure_timing, nonlinear_seconds, &
        fft_seconds, kinetic_seconds, other_seconds, total_seconds) &
        bind(C, name="gp3d_cuda_step") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: dt, hbar, mass, g, target_norm, volume_element
      integer(c_int), value :: imaginary_time, measure_timing
      real(c_double), intent(out) :: nonlinear_seconds, fft_seconds
      real(c_double), intent(out) :: kinetic_seconds, other_seconds, total_seconds
      integer(c_int) :: status
    end function c_gpu_step

    function c_gpu_argle_step(handle, alpha, reaction_constant, &
        nonlinear_coefficient, potential_scale, dtau, velocity_amplitude, &
        compute_metrics, max_delta_rate, mean_density, min_density) &
        bind(C, name="gp3d_cuda_argle_step") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: alpha, reaction_constant, nonlinear_coefficient
      real(c_double), value :: potential_scale, dtau, velocity_amplitude
      integer(c_int), value :: compute_metrics
      real(c_double), intent(out) :: max_delta_rate, mean_density, min_density
      integer(c_int) :: status
    end function c_gpu_argle_step

    function c_gpu_diagnostics(handle, hbar, mass, g, volume_element, &
        norm_value, energy_value) bind(C, name="gp3d_cuda_diagnostics") &
        result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: hbar, mass, g, volume_element
      real(c_double), intent(out) :: norm_value, energy_value
      integer(c_int) :: status
    end function c_gpu_diagnostics

    subroutine c_gpu_destroy(handle) bind(C, name="gp3d_cuda_destroy")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine c_gpu_destroy

    subroutine c_gpu_get_last_error(buffer, buffer_size) &
        bind(C, name="gp3d_cuda_get_last_error")
      import :: c_char, c_int
      character(c_char), intent(out) :: buffer(*)
      integer(c_int), value :: buffer_size
    end subroutine c_gpu_get_last_error
  end interface

contains

  subroutine gp3d_gpu_init(context, grid, need_argle, device)
    ! CUDAデバイス、cuFFT plan、常駐配列を一度だけ確保する。
    type(gp3d_gpu_context_t), intent(inout) :: context
    type(gp3d_grid_t), intent(in) :: grid
    logical, intent(in) :: need_argle
    integer, intent(in), optional :: device

    integer :: selected_device
    integer(c_int) :: status

    if (kind(1.0_dp) /= c_double) error stop "CUDA backend requires double-precision dp"
    if (grid%local_nz /= grid%nz) error stop "single-GPU backend requires the full z domain"
    if (c_associated(context%handle)) error stop "CUDA context is already initialized"
    selected_device = 0
    if (present(device)) selected_device = device
    status = c_gpu_create(context%handle, int(grid%nx, c_int), &
      int(grid%ny, c_int), int(grid%nz, c_int), &
      merge(1_c_int, 0_c_int, need_argle), int(selected_device, c_int))
    call require_success(status, "initialize CUDA/cuFFT context")
  end subroutine gp3d_gpu_init

  subroutine gp3d_gpu_upload(context, state, grid)
    ! 初期psi、ポテンシャル、波数配列をホストからGPUへ転送する。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_state_t), intent(in) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer(c_int) :: status

    call require_context(context)
    status = c_gpu_upload(context%handle, state%psi, state%potential, &
      grid%kx, grid%ky, grid%kz)
    call require_success(status, "upload GPE state to GPU")
  end subroutine gp3d_gpu_upload

  subroutine gp3d_gpu_download(context, state)
    ! 出力または終了時だけ、現在のpsiをGPUからホストへ戻す。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_state_t), intent(inout) :: state
    integer(c_int) :: status

    call require_context(context)
    status = c_gpu_download(context%handle, state%psi)
    call require_success(status, "download GPE state from GPU")
  end subroutine gp3d_gpu_download

  subroutine gp3d_gpu_step(context, grid, params, timing)
    ! 1ステップ全体をGPU上で実行し、ループ中のpsi転送を発生させない。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    type(gp3d_step_timing_t), intent(inout), optional :: timing

    real(c_double) :: nonlinear_seconds, fft_seconds, kinetic_seconds
    real(c_double) :: other_seconds, total_seconds
    real(dp) :: volume_element
    integer(c_int) :: status, measure_timing

    call require_context(context)
    volume_element = grid%dx * grid%dy * grid%dz
    measure_timing = merge(1_c_int, 0_c_int, present(timing))
    status = c_gpu_step(context%handle, params%dt, params%hbar, params%mass, &
      params%g, params%norm, volume_element, &
      merge(1_c_int, 0_c_int, params%imaginary_time), measure_timing, &
      nonlinear_seconds, fft_seconds, kinetic_seconds, other_seconds, total_seconds)
    call require_success(status, "advance one CUDA split-step")
    if (present(timing)) then
      timing%steps = timing%steps + 1
      timing%nonlinear_seconds = timing%nonlinear_seconds + nonlinear_seconds
      timing%fft_seconds = timing%fft_seconds + fft_seconds
      timing%kinetic_seconds = timing%kinetic_seconds + kinetic_seconds
      timing%other_seconds = timing%other_seconds + other_seconds
      timing%total_seconds = timing%total_seconds + total_seconds
    end if
  end subroutine gp3d_gpu_step

  subroutine gp3d_gpu_relax_taylor_green_argle(context, params, model_cfg)
    ! 初期条件緩和もGPU常駐データとcuFFTを使って実行する。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_params_t), intent(in) :: params
    type(gp3d_model_config_t), intent(in) :: model_cfg

    integer :: step
    integer(c_int) :: status, compute_metrics
    real(dp) :: alpha, reaction_constant, nonlinear_coefficient, potential_scale
    real(c_double) :: max_delta_rate, mean_density, min_density
    real(dp) :: pseudo_time
    logical :: converged, report_step, need_metrics

    if (.not. model_cfg%argle_enabled) return
    if (model_cfg%argle_steps <= 0) error stop "ARGLE steps must be positive"
    if (model_cfg%argle_dtau <= 0.0_dp) error stop "ARGLE time step must be positive"
    if (model_cfg%argle_tolerance < 0.0_dp) error stop "ARGLE tolerance must be non-negative"

    alpha = params%hbar / (2.0_dp * params%mass)
    if (model_cfg%use_dimensionless_parameters) then
      reaction_constant = model_cfg%beta
      nonlinear_coefficient = model_cfg%beta
      potential_scale = 1.0_dp
    else
      reaction_constant = model_cfg%mu / params%hbar
      nonlinear_coefficient = params%g / params%hbar
      potential_scale = 1.0_dp / params%hbar
    end if

    write(*,'(a)') "# CUDA ARGLE: unconstrained minimization of the driven energy"
    write(*,'(a)') "# step pseudo_time max_delta_rate mean_density min_density"
    do step = 1, model_cfg%argle_steps
      report_step = step == 1 .or. step == model_cfg%argle_steps
      if (model_cfg%argle_output_every > 0) then
        report_step = report_step .or. mod(step, model_cfg%argle_output_every) == 0
      end if
      need_metrics = report_step .or. model_cfg%argle_tolerance > 0.0_dp
      compute_metrics = merge(1_c_int, 0_c_int, need_metrics)
      status = c_gpu_argle_step(context%handle, alpha, reaction_constant, &
        nonlinear_coefficient, potential_scale, model_cfg%argle_dtau, &
        model_cfg%tg_velocity_amplitude, compute_metrics, max_delta_rate, &
        mean_density, min_density)
      call require_success(status, "advance one CUDA ARGLE step")
      converged = model_cfg%argle_tolerance > 0.0_dp .and. &
        max_delta_rate < model_cfg%argle_tolerance
      report_step = report_step .or. converged
      if (report_step) then
        pseudo_time = real(step, dp) * model_cfg%argle_dtau
        write(*,'(i8,1x,4(es16.8,1x))') step, pseudo_time, &
          max_delta_rate, mean_density, min_density
      end if
      if (converged) exit
    end do
  end subroutine gp3d_gpu_relax_taylor_green_argle

  subroutine gp3d_gpu_diagnostics(context, grid, params, norm_value, energy_value)
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    real(dp), intent(out) :: norm_value, energy_value

    real(c_double) :: c_norm, c_energy
    integer(c_int) :: status

    call require_context(context)
    status = c_gpu_diagnostics(context%handle, params%hbar, params%mass, &
      params%g, grid%dx * grid%dy * grid%dz, c_norm, c_energy)
    call require_success(status, "compute CUDA GPE diagnostics")
    norm_value = c_norm
    energy_value = c_energy
  end subroutine gp3d_gpu_diagnostics

  subroutine gp3d_gpu_finalize(context)
    type(gp3d_gpu_context_t), intent(inout) :: context

    if (c_associated(context%handle)) call c_gpu_destroy(context%handle)
    context%handle = c_null_ptr
  end subroutine gp3d_gpu_finalize

  subroutine require_context(context)
    type(gp3d_gpu_context_t), intent(in) :: context

    if (.not. c_associated(context%handle)) error stop "CUDA context is not initialized"
  end subroutine require_context

  subroutine require_success(status, operation)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: operation

    character(c_char) :: buffer(512)
    character(len=512) :: message
    integer :: i

    if (status == 0_c_int) return
    buffer = c_null_char
    call c_gpu_get_last_error(buffer, int(size(buffer), c_int))
    message = ""
    do i = 1, size(buffer)
      if (buffer(i) == c_null_char) exit
      message(i:i) = buffer(i)
    end do
    write(error_unit,'(3a)') "CUDA backend error during ", trim(operation), ": " // trim(message)
    error stop "CUDA backend failure"
  end subroutine require_success

end module gp3d_gpu
