!> MPI + cuFFTMp バックエンドを Fortran から利用するための C バインディング。
!> 波動関数は各 rank の GPU に常駐し、出力時だけ局所 z スラブをホストへ戻す。
module gp3d_gpu
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, &
    c_double, c_double_complex, c_char, c_null_char
  use, intrinsic :: iso_fortran_env, only: error_unit
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_params_t, gp3d_state_t, &
    gp3d_model_config_t
  use gp3d_solver, only: gp3d_step_timing_t
  use gp3d_mpi, only: gp3d_mpi_t
  implicit none
  private

  type, public :: gp3d_gpu_context_t
    private
    type(c_ptr) :: handle = c_null_ptr
    integer(c_int) :: comm = 0_c_int
    integer :: rank = 0
  end type gp3d_gpu_context_t

  public :: gp3d_gpu_init
  public :: gp3d_gpu_upload
  public :: gp3d_gpu_download
  public :: gp3d_gpu_step
  public :: gp3d_gpu_relax_taylor_green_argle
  public :: gp3d_gpu_diagnostics
  public :: gp3d_gpu_finalize

  interface
    function c_gpu_create(handle, nx, ny, nz, local_nz, k_start, &
        need_argle, comm) bind(C, name="gp3d_cufftmp_create") result(status)
      import :: c_ptr, c_int
      type(c_ptr), intent(out) :: handle
      integer(c_int), value :: nx, ny, nz, local_nz, k_start
      integer(c_int), value :: need_argle, comm
      integer(c_int) :: status
    end function c_gpu_create

    function c_gpu_upload(handle, psi, potential, kx, ky, kz) &
        bind(C, name="gp3d_cufftmp_upload") result(status)
      import :: c_ptr, c_int, c_double, c_double_complex
      type(c_ptr), value :: handle
      complex(c_double_complex), intent(in) :: psi(*)
      real(c_double), intent(in) :: potential(*), kx(*), ky(*), kz(*)
      integer(c_int) :: status
    end function c_gpu_upload

    function c_gpu_download(handle, psi) &
        bind(C, name="gp3d_cufftmp_download") result(status)
      import :: c_ptr, c_int, c_double_complex
      type(c_ptr), value :: handle
      complex(c_double_complex), intent(out) :: psi(*)
      integer(c_int) :: status
    end function c_gpu_download

    function c_gpu_step(handle, dt, hbar, mass, g, target_norm, &
        volume_element, imaginary_time, measure_timing, nonlinear_seconds, &
        fft_seconds, kinetic_seconds, other_seconds, total_seconds) &
        bind(C, name="gp3d_cufftmp_step") result(status)
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
        bind(C, name="gp3d_cufftmp_argle_step") result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: alpha, reaction_constant, nonlinear_coefficient
      real(c_double), value :: potential_scale, dtau, velocity_amplitude
      integer(c_int), value :: compute_metrics
      real(c_double), intent(out) :: max_delta_rate, mean_density, min_density
      integer(c_int) :: status
    end function c_gpu_argle_step

    function c_gpu_diagnostics(handle, hbar, mass, g, volume_element, &
        norm_value, energy_value) bind(C, name="gp3d_cufftmp_diagnostics") &
        result(status)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: hbar, mass, g, volume_element
      real(c_double), intent(out) :: norm_value, energy_value
      integer(c_int) :: status
    end function c_gpu_diagnostics

    subroutine c_gpu_destroy(handle) bind(C, name="gp3d_cufftmp_destroy")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine c_gpu_destroy

    subroutine c_gpu_get_last_error(buffer, buffer_size) &
        bind(C, name="gp3d_cufftmp_get_last_error")
      import :: c_char, c_int
      character(c_char), intent(out) :: buffer(*)
      integer(c_int), value :: buffer_size
    end subroutine c_gpu_get_last_error

    subroutine c_gpu_abort(comm, error_code) bind(C, name="gp3d_cufftmp_abort")
      import :: c_int
      integer(c_int), value :: comm, error_code
    end subroutine c_gpu_abort
  end interface

contains

  subroutine gp3d_gpu_init(context, grid, mpi, need_argle)
    ! cuFFTMp plan と全 rank で対称な分散 GPU バッファを一度だけ確保する。
    type(gp3d_gpu_context_t), intent(inout) :: context
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_mpi_t), intent(in) :: mpi
    logical, intent(in) :: need_argle

    integer(c_int) :: status

    if (kind(1.0_dp) /= c_double) error stop "cuFFTMp backend requires double-precision dp"
    if (.not. mpi%enabled) error stop "cuFFTMp backend requires MPI"
    if (grid%rank /= mpi%rank .or. grid%nprocs /= mpi%nprocs) then
      error stop "cuFFTMp grid and MPI decomposition disagree"
    end if
    if (c_associated(context%handle)) error stop "cuFFTMp context is already initialized"

    context%comm = int(mpi%comm, c_int)
    context%rank = mpi%rank
    status = c_gpu_create(context%handle, int(grid%nx, c_int), &
      int(grid%ny, c_int), int(grid%nz, c_int), int(grid%local_nz, c_int), &
      int(grid%k_start, c_int), merge(1_c_int, 0_c_int, need_argle), &
      context%comm)
    call require_success(context, status, "initialize MPI/cuFFTMp context")
  end subroutine gp3d_gpu_init

  subroutine gp3d_gpu_upload(context, state, grid)
    ! 現 rank が所有する実空間 z スラブと全方向の波数表を GPU へ転送する。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_state_t), intent(in) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer(c_int) :: status

    call require_context(context)
    if (size(state%psi, 3) /= grid%local_nz) then
      error stop "cuFFTMp upload shape does not match the local z slab"
    end if
    status = c_gpu_upload(context%handle, state%psi, state%potential, &
      grid%kx, grid%ky, grid%kz)
    call require_success(context, status, "upload distributed GPE state")
  end subroutine gp3d_gpu_upload

  subroutine gp3d_gpu_download(context, state)
    ! SLF 出力や終了処理のため、現在の局所 z スラブだけをホストへ戻す。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_state_t), intent(inout) :: state
    integer(c_int) :: status

    call require_context(context)
    status = c_gpu_download(context%handle, state%psi)
    call require_success(context, status, "download distributed GPE state")
  end subroutine gp3d_gpu_download

  subroutine gp3d_gpu_step(context, grid, params, timing)
    ! 局所非線形項、分散 FFT、Fourier 空間運動項を GPU 常駐のまま進める。
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
    call require_success(context, status, "advance one distributed CUDA split step")
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
    ! ARGLE の複数回の微分 FFT も cuFFTMp で集団実行する。
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

    if (context%rank == 0) then
      write(*,'(a)') "# cuFFTMp ARGLE: unconstrained distributed minimization"
      write(*,'(a)') "# step pseudo_time max_delta_rate mean_density min_density"
    end if
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
      call require_success(context, status, "advance one distributed CUDA ARGLE step")
      converged = model_cfg%argle_tolerance > 0.0_dp .and. &
        max_delta_rate < model_cfg%argle_tolerance
      report_step = report_step .or. converged
      if (report_step .and. context%rank == 0) then
        pseudo_time = real(step, dp) * model_cfg%argle_dtau
        write(*,'(i8,1x,4(es16.8,1x))') step, pseudo_time, &
          max_delta_rate, mean_density, min_density
      end if
      if (converged) exit
    end do
  end subroutine gp3d_gpu_relax_taylor_green_argle

  subroutine gp3d_gpu_diagnostics(context, grid, params, norm_value, energy_value)
    ! GPU 上の局所和を計算し、C ブリッジ内部の MPI_Allreduce で全体値にする。
    type(gp3d_gpu_context_t), intent(in) :: context
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    real(dp), intent(out) :: norm_value, energy_value

    real(c_double) :: c_norm, c_energy
    integer(c_int) :: status

    call require_context(context)
    status = c_gpu_diagnostics(context%handle, params%hbar, params%mass, &
      params%g, grid%dx * grid%dy * grid%dz, c_norm, c_energy)
    call require_success(context, status, "compute distributed CUDA diagnostics")
    norm_value = c_norm
    energy_value = c_energy
  end subroutine gp3d_gpu_diagnostics

  subroutine gp3d_gpu_finalize(context)
    ! 全 rank が MPI_Finalize より前に分散バッファと plan を解放する。
    type(gp3d_gpu_context_t), intent(inout) :: context

    if (c_associated(context%handle)) call c_gpu_destroy(context%handle)
    context%handle = c_null_ptr
    context%comm = 0_c_int
    context%rank = 0
  end subroutine gp3d_gpu_finalize

  subroutine require_context(context)
    type(gp3d_gpu_context_t), intent(in) :: context

    if (.not. c_associated(context%handle)) error stop "cuFFTMp context is not initialized"
  end subroutine require_context

  subroutine require_success(context, status, operation)
    type(gp3d_gpu_context_t), intent(in) :: context
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: operation

    character(c_char) :: buffer(1024)
    character(len=1024) :: message
    integer :: i

    if (status == 0_c_int) return
    buffer = c_null_char
    call c_gpu_get_last_error(buffer, int(size(buffer), c_int))
    message = ""
    do i = 1, size(buffer)
      if (buffer(i) == c_null_char) exit
      message(i:i) = buffer(i)
    end do
    write(error_unit,'(a,i0,3a)') "cuFFTMp rank ", context%rank, &
      " failed during ", trim(operation), ": " // trim(message)
    call c_gpu_abort(context%comm, status)
    error stop "cuFFTMp backend failure"
  end subroutine require_success

end module gp3d_gpu
