!> GPEの状態初期化、診断量、Split-operator時間積分、ARGLE緩和を実装する。
!> FFTとMPIは抽象化された同一APIを呼ぶため、CPUバックエンドを差し替えても本体は変わらない。
module gp3d_solver
  use, intrinsic :: iso_fortran_env, only: int64
  use gp3d_types, only: dp, pi, gp3d_grid_t, gp3d_params_t, gp3d_state_t, &
    gp3d_model_config_t
  use gp3d_fft, only: gp3d_fft_plan_t, gp3d_fft_forward, gp3d_fft_inverse
  use gp3d_mpi, only: gp3d_mpi_t, gp3d_mpi_sum_real, gp3d_mpi_max_real
  use gp3d_openmp, only: gp3d_openmp_active
  implicit none
  private

  type, public :: gp3d_step_timing_t
    integer :: steps = 0
    real(dp) :: nonlinear_seconds = 0.0_dp
    real(dp) :: fft_seconds = 0.0_dp
    real(dp) :: kinetic_seconds = 0.0_dp
    real(dp) :: other_seconds = 0.0_dp
    real(dp) :: total_seconds = 0.0_dp
  end type gp3d_step_timing_t

  public :: gp3d_state_allocate
  public :: gp3d_set_harmonic_potential
  public :: gp3d_set_gaussian_initial_state
  public :: gp3d_normalize
  public :: gp3d_density_norm
  public :: gp3d_energy
  public :: gp3d_step_split_operator
  public :: gp3d_report_step_timing
  public :: gp3d_relax_taylor_green_argle

contains

  subroutine gp3d_state_allocate(state, grid)
    type(gp3d_state_t), intent(out) :: state
    type(gp3d_grid_t), intent(in) :: grid

    allocate(state%psi(grid%nx, grid%ny, grid%local_nz))
    allocate(state%potential(grid%nx, grid%ny, grid%local_nz))
    state%psi = (0.0_dp, 0.0_dp)
    state%potential = 0.0_dp
  end subroutine gp3d_state_allocate

  subroutine gp3d_set_harmonic_potential(state, grid, omega_x, omega_y, omega_z)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    real(dp), intent(in) :: omega_x, omega_y, omega_z

    integer :: i, j, k, kg

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) private(kg)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          state%potential(i,j,k) = 0.5_dp * ( &
            (omega_x * grid%x(i))**2 + &
            (omega_y * grid%y(j))**2 + &
            (omega_z * grid%z(kg))**2)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine gp3d_set_harmonic_potential

  subroutine gp3d_set_gaussian_initial_state(state, grid, sigma_x, sigma_y, sigma_z)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    real(dp), intent(in) :: sigma_x, sigma_y, sigma_z

    integer :: i, j, k, kg
    real(dp) :: exponent_value

    if (sigma_x <= 0.0_dp .or. sigma_y <= 0.0_dp .or. sigma_z <= 0.0_dp) then
      error stop "Gaussian widths must be positive"
    end if

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) private(kg, exponent_value)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          exponent_value = -0.5_dp * ( &
            (grid%x(i) / sigma_x)**2 + &
            (grid%y(j) / sigma_y)**2 + &
            (grid%z(kg) / sigma_z)**2)
          state%psi(i,j,k) = cmplx(exp(exponent_value), 0.0_dp, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine gp3d_set_gaussian_initial_state

  real(dp) function gp3d_density_norm(state, grid, mpi) result(norm_value)
    type(gp3d_state_t), intent(in) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_mpi_t), intent(in), optional :: mpi
    real(dp) :: local_norm
    integer :: i, j, k

    local_norm = 0.0_dp
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) reduction(+:local_norm)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          local_norm = local_norm + abs(state%psi(i,j,k))**2
        end do
      end do
    end do
    !$omp end parallel do
    local_norm = local_norm * grid%dx * grid%dy * grid%dz
    if (present(mpi)) then
      call gp3d_mpi_sum_real(mpi, local_norm, norm_value)
    else
      norm_value = local_norm
    end if
  end function gp3d_density_norm

  subroutine gp3d_normalize(state, grid, target_norm, mpi)
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    real(dp), intent(in) :: target_norm
    type(gp3d_mpi_t), intent(in), optional :: mpi

    integer :: i, j, k
    real(dp) :: current_norm, scale

    if (present(mpi)) then
      current_norm = gp3d_density_norm(state, grid, mpi)
    else
      current_norm = gp3d_density_norm(state, grid)
    end if
    if (current_norm <= 0.0_dp) error stop "cannot normalize a zero wave function"
    scale = sqrt(target_norm / current_norm)
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          state%psi(i,j,k) = state%psi(i,j,k) * scale
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine gp3d_normalize

  subroutine gp3d_step_split_operator(state, grid, params, fft_plan, mpi, timing)
    ! Strang分割: 局所半ステップ -> FFT -> 運動項1ステップ -> 逆FFT -> 局所半ステップ。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    type(gp3d_fft_plan_t), intent(in) :: fft_plan
    type(gp3d_mpi_t), intent(in), optional :: mpi
    type(gp3d_step_timing_t), intent(inout), optional :: timing

    complex(dp), allocatable :: psi_k(:,:,:)
    real(dp) :: start_time = 0.0_dp, step_start = 0.0_dp
    real(dp) :: nonlinear_elapsed = 0.0_dp, fft_elapsed = 0.0_dp
    real(dp) :: kinetic_elapsed = 0.0_dp, total_elapsed = 0.0_dp

    if (present(timing)) then
      step_start = wall_time_seconds()
      nonlinear_elapsed = 0.0_dp
      fft_elapsed = 0.0_dp
      kinetic_elapsed = 0.0_dp
      start_time = wall_time_seconds()
    end if
    call apply_local_half_step(state, grid, params)
    if (present(timing)) nonlinear_elapsed = nonlinear_elapsed + wall_time_seconds() - start_time

    allocate(psi_k(grid%nx, grid%ny, grid%local_nz))
    if (present(timing)) start_time = wall_time_seconds()
    call gp3d_fft_forward(fft_plan, state%psi, psi_k)
    if (present(timing)) fft_elapsed = fft_elapsed + wall_time_seconds() - start_time
    if (present(timing)) start_time = wall_time_seconds()
    call apply_kinetic_step(psi_k, grid, params)
    if (present(timing)) kinetic_elapsed = kinetic_elapsed + wall_time_seconds() - start_time
    if (present(timing)) start_time = wall_time_seconds()
    call gp3d_fft_inverse(fft_plan, psi_k, state%psi)
    if (present(timing)) fft_elapsed = fft_elapsed + wall_time_seconds() - start_time
    deallocate(psi_k)

    if (present(timing)) start_time = wall_time_seconds()
    call apply_local_half_step(state, grid, params)
    if (present(timing)) nonlinear_elapsed = nonlinear_elapsed + wall_time_seconds() - start_time

    if (params%imaginary_time) then
      if (present(mpi)) then
        call gp3d_normalize(state, grid, params%norm, mpi)
      else
        call gp3d_normalize(state, grid, params%norm)
      end if
    end if

    if (present(timing)) then
      total_elapsed = wall_time_seconds() - step_start
      timing%steps = timing%steps + 1
      timing%nonlinear_seconds = timing%nonlinear_seconds + nonlinear_elapsed
      timing%fft_seconds = timing%fft_seconds + fft_elapsed
      timing%kinetic_seconds = timing%kinetic_seconds + kinetic_elapsed
      timing%other_seconds = timing%other_seconds + &
        max(0.0_dp, total_elapsed - nonlinear_elapsed - fft_elapsed - kinetic_elapsed)
      timing%total_seconds = timing%total_seconds + total_elapsed
    end if
  end subroutine gp3d_step_split_operator

  subroutine gp3d_report_step_timing(timing, mpi)
    type(gp3d_step_timing_t), intent(in) :: timing
    type(gp3d_mpi_t), intent(in) :: mpi

    character(len=24), parameter :: labels(5) = [character(len=24) :: &
      "nonlinear_local", "fft_forward_inverse", "kinetic_spectral", &
      "allocation_and_other", "step_total"]
    real(dp) :: local_values(5), sum_values(5), max_values(5)
    real(dp) :: mean_values(5), percent
    integer :: i

    if (timing%steps <= 0) then
      if (mpi%rank == mpi%root) write(*,'(a)') "# split-step timing: no steps measured"
      return
    end if

    local_values = [timing%nonlinear_seconds, timing%fft_seconds, &
      timing%kinetic_seconds, timing%other_seconds, timing%total_seconds]
    do i = 1, size(local_values)
      call gp3d_mpi_sum_real(mpi, local_values(i), sum_values(i))
      call gp3d_mpi_max_real(mpi, local_values(i), max_values(i))
    end do
    mean_values = sum_values / real(max(1, mpi%nprocs), dp)

    if (mpi%rank /= mpi%root) return
    write(*,'(a,i0)') "# split-step timing steps=", timing%steps
    write(*,'(a)') "# region rank_mean_s rank_max_s max_s_per_step mean_percent"
    do i = 1, size(local_values)
      percent = 0.0_dp
      if (mean_values(5) > 0.0_dp) percent = 100.0_dp * mean_values(i) / mean_values(5)
      write(*,'(a24,1x,3(es16.8,1x),f10.3)') adjustl(labels(i)), &
        mean_values(i), max_values(i), max_values(i) / real(timing%steps, dp), percent
    end do
    write(*,'(a)') "# FFT timing contains one forward and one inverse transform per step."
    write(*,'(a)') "# Nonlinear timing contains both local half steps."
  end subroutine gp3d_report_step_timing

  subroutine gp3d_relax_taylor_green_argle(state, grid, params, model_cfg, fft_plan, mpi)
    ! Taylor-Green速度場を拘束として加え、無拘束・半陰的ARGLEで音波成分を緩和する。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    type(gp3d_model_config_t), intent(in) :: model_cfg
    type(gp3d_fft_plan_t), intent(in) :: fft_plan
    type(gp3d_mpi_t), intent(in) :: mpi

    complex(dp), allocatable :: psi_old(:,:,:), psi_k(:,:,:), derivative_k(:,:,:)
    complex(dp), allocatable :: grad_x(:,:,:), grad_y(:,:,:), rhs(:,:,:), rhs_k(:,:,:)
    complex(dp), allocatable :: next_k(:,:,:)
    real(dp), allocatable :: velocity_x(:,:,:), velocity_y(:,:,:), velocity2(:,:,:)
    integer :: i, j, k, kg, step
    real(dp) :: alpha, beta, k2, denominator, explicit_factor, reaction
    real(dp) :: local_error, point_error, global_error, local_min_density, global_min_density
    real(dp) :: neg_global_min_density, norm_value, mean_density, pseudo_time
    logical :: converged, report_step

    if (.not. model_cfg%argle_enabled) return
    if (model_cfg%argle_steps <= 0) error stop "ARGLE steps must be positive"
    if (model_cfg%argle_dtau <= 0.0_dp) error stop "ARGLE time step must be positive"
    if (model_cfg%argle_tolerance < 0.0_dp) error stop "ARGLE tolerance must be non-negative"
    if (params%mass <= 0.0_dp .or. params%hbar <= 0.0_dp) then
      error stop "ARGLE requires positive mass and hbar"
    end if

    alpha = params%hbar / (2.0_dp * params%mass)
    beta = params%g / params%hbar
    allocate(psi_old(grid%nx, grid%ny, grid%local_nz))
    allocate(psi_k(grid%nx, grid%ny, grid%local_nz))
    allocate(derivative_k(grid%nx, grid%ny, grid%local_nz))
    allocate(grad_x(grid%nx, grid%ny, grid%local_nz))
    allocate(grad_y(grid%nx, grid%ny, grid%local_nz))
    allocate(rhs(grid%nx, grid%ny, grid%local_nz))
    allocate(rhs_k(grid%nx, grid%ny, grid%local_nz))
    allocate(next_k(grid%nx, grid%ny, grid%local_nz))
    allocate(velocity_x(grid%nx, grid%ny, grid%local_nz))
    allocate(velocity_y(grid%nx, grid%ny, grid%local_nz))
    allocate(velocity2(grid%nx, grid%ny, grid%local_nz))

    call fill_taylor_green_velocity(grid, model_cfg%tg_velocity_amplitude, &
      velocity_x, velocity_y, velocity2)

    if (mpi%rank == mpi%root) then
      write(*,'(a)') "# ARGLE: unconstrained minimization of the driven energy"
      write(*,'(a)') "# step pseudo_time max_delta_rate mean_density min_density"
    end if
    do step = 1, model_cfg%argle_steps
      psi_old = state%psi
      call gp3d_fft_forward(fft_plan, psi_old, psi_k)

      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active)
      do k = 1, grid%local_nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            derivative_k(i,j,k) = cmplx(0.0_dp, grid%kx(i), kind=dp) * psi_k(i,j,k)
          end do
        end do
      end do
      !$omp end parallel do
      call gp3d_fft_inverse(fft_plan, derivative_k, grad_x)

      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active)
      do k = 1, grid%local_nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            derivative_k(i,j,k) = cmplx(0.0_dp, grid%ky(j), kind=dp) * psi_k(i,j,k)
          end do
        end do
      end do
      !$omp end parallel do
      call gp3d_fft_inverse(fft_plan, derivative_k, grad_y)

      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) private(reaction)
      do k = 1, grid%local_nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            if (model_cfg%use_dimensionless_parameters) then
              reaction = beta * (1.0_dp - abs(psi_old(i,j,k))**2) - &
                state%potential(i,j,k) - velocity2(i,j,k) / (4.0_dp * alpha)
            else
              reaction = (model_cfg%mu - state%potential(i,j,k) - &
                params%g * abs(psi_old(i,j,k))**2) / params%hbar - &
                velocity2(i,j,k) / (4.0_dp * alpha)
            end if
            rhs(i,j,k) = reaction * psi_old(i,j,k) - &
              cmplx(0.0_dp, 1.0_dp, kind=dp) * &
              (velocity_x(i,j,k) * grad_x(i,j,k) + velocity_y(i,j,k) * grad_y(i,j,k))
          end do
        end do
      end do
      !$omp end parallel do
      call gp3d_fft_forward(fft_plan, rhs, rhs_k)

      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
      !$omp& private(kg, k2, denominator, explicit_factor)
      do k = 1, grid%local_nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            kg = grid%k_start + k - 1
            k2 = grid%kx(i)**2 + grid%ky(j)**2 + grid%kz(kg)**2
            denominator = 1.0_dp + 0.5_dp * model_cfg%argle_dtau * alpha * k2
            explicit_factor = 1.0_dp - 0.5_dp * model_cfg%argle_dtau * alpha * k2
            next_k(i,j,k) = (explicit_factor * psi_k(i,j,k) + &
              model_cfg%argle_dtau * rhs_k(i,j,k)) / denominator
          end do
        end do
      end do
      !$omp end parallel do
      call gp3d_fft_inverse(fft_plan, next_k, state%psi)

      local_error = 0.0_dp
      !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
      !$omp& private(point_error) reduction(max:local_error)
      do k = 1, grid%local_nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            point_error = abs(state%psi(i,j,k) - psi_old(i,j,k))
            local_error = max(local_error, point_error)
          end do
        end do
      end do
      !$omp end parallel do
      local_error = local_error / model_cfg%argle_dtau
      call gp3d_mpi_max_real(mpi, local_error, global_error)
      converged = model_cfg%argle_tolerance > 0.0_dp .and. &
        global_error < model_cfg%argle_tolerance
      report_step = step == 1 .or. step == model_cfg%argle_steps .or. converged
      if (model_cfg%argle_output_every > 0) then
        report_step = report_step .or. mod(step, model_cfg%argle_output_every) == 0
      end if
      if (report_step) then
        norm_value = gp3d_density_norm(state, grid, mpi)
        mean_density = norm_value / (grid%lx * grid%ly * grid%lz)
        local_min_density = huge(1.0_dp)
        !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
        !$omp& reduction(min:local_min_density)
        do k = 1, grid%local_nz
          do j = 1, grid%ny
            do i = 1, grid%nx
              local_min_density = min(local_min_density, abs(state%psi(i,j,k))**2)
            end do
          end do
        end do
        !$omp end parallel do
        call gp3d_mpi_max_real(mpi, -local_min_density, neg_global_min_density)
        global_min_density = -neg_global_min_density
        pseudo_time = real(step, dp) * model_cfg%argle_dtau
      end if
      if (report_step .and. mpi%rank == mpi%root) then
        write(*,'(i8,1x,4(es16.8,1x))') step, pseudo_time, global_error, &
          mean_density, global_min_density
      end if
      if (converged) exit
    end do

    deallocate(psi_old, psi_k, derivative_k, grad_x, grad_y, rhs, rhs_k, next_k)
    deallocate(velocity_x, velocity_y, velocity2)
  end subroutine gp3d_relax_taylor_green_argle

  real(dp) function gp3d_energy(state, grid, params, fft_plan, mpi) result(energy)
    type(gp3d_state_t), intent(in) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    type(gp3d_fft_plan_t), intent(in) :: fft_plan
    type(gp3d_mpi_t), intent(in), optional :: mpi

    complex(dp), allocatable :: psi_k(:,:,:)
    integer :: i, j, k, kg
    real(dp) :: volume_element, spectral_scale, k2, kinetic, interaction, local_energy

    allocate(psi_k(grid%nx, grid%ny, grid%local_nz))
    call gp3d_fft_forward(fft_plan, state%psi, psi_k)

    kinetic = 0.0_dp
    interaction = 0.0_dp
    spectral_scale = grid%dx * grid%dy * grid%dz / real(grid%nx * grid%ny * grid%nz, dp)
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(kg, k2) reduction(+:kinetic, interaction)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          k2 = grid%kx(i)**2 + grid%ky(j)**2 + grid%kz(kg)**2
          kinetic = kinetic + 0.5_dp * params%hbar**2 / params%mass * k2 * abs(psi_k(i,j,k))**2
          interaction = interaction + state%potential(i,j,k) * abs(state%psi(i,j,k))**2 + &
            0.5_dp * params%g * abs(state%psi(i,j,k))**4
        end do
      end do
    end do
    !$omp end parallel do
    kinetic = kinetic * spectral_scale

    volume_element = grid%dx * grid%dy * grid%dz
    local_energy = kinetic + interaction * volume_element
    if (present(mpi)) then
      call gp3d_mpi_sum_real(mpi, local_energy, energy)
    else
      energy = local_energy
    end if

    deallocate(psi_k)
  end function gp3d_energy

  subroutine apply_local_half_step(state, grid, params)
    ! V + g|psi|^2 は格子点ごとに独立なので、実空間で指数演算子を直接乗算する。
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params

    integer :: i, j, k
    real(dp) :: local_energy, tau
    complex(dp) :: factor

    tau = 0.5_dp * params%dt / params%hbar
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(local_energy, factor)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          local_energy = state%potential(i,j,k) + params%g * abs(state%psi(i,j,k))**2
          if (params%imaginary_time) then
            factor = cmplx(exp(-tau * local_energy), 0.0_dp, kind=dp)
          else
            factor = exp(cmplx(0.0_dp, -tau * local_energy, kind=dp))
          end if
          state%psi(i,j,k) = state%psi(i,j,k) * factor
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine apply_local_half_step

  subroutine apply_kinetic_step(psi_k, grid, params)
    ! 運動項は波数空間で対角化され、各Fourier係数へ位相因子を乗算する。
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_params_t), intent(in) :: params
    complex(dp), intent(inout) :: psi_k(grid%nx, grid%ny, grid%local_nz)

    integer :: i, j, k, kg
    real(dp) :: kinetic_energy, tau
    complex(dp) :: factor

    tau = params%dt / params%hbar
    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(kg, kinetic_energy, factor)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          kinetic_energy = 0.5_dp * params%hbar**2 / params%mass * &
            (grid%kx(i)**2 + grid%ky(j)**2 + grid%kz(kg)**2)
          if (params%imaginary_time) then
            factor = cmplx(exp(-tau * kinetic_energy), 0.0_dp, kind=dp)
          else
            factor = exp(cmplx(0.0_dp, -tau * kinetic_energy, kind=dp))
          end if
          psi_k(i,j,k) = psi_k(i,j,k) * factor
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine apply_kinetic_step

  subroutine fill_taylor_green_velocity(grid, amplitude, velocity_x, velocity_y, velocity2)
    type(gp3d_grid_t), intent(in) :: grid
    real(dp), intent(in) :: amplitude
    real(dp), intent(out) :: velocity_x(grid%nx, grid%ny, grid%local_nz)
    real(dp), intent(out) :: velocity_y(grid%nx, grid%ny, grid%local_nz)
    real(dp), intent(out) :: velocity2(grid%nx, grid%ny, grid%local_nz)

    integer :: i, j, k, kg
    real(dp) :: x_angle, y_angle, z_angle
    real(dp) :: x_origin, y_origin, z_origin

    x_origin = grid%x(1) - 0.5_dp * grid%dx
    y_origin = grid%y(1) - 0.5_dp * grid%dy
    z_origin = grid%z(1) - 0.5_dp * grid%dz

    !$omp parallel do collapse(3) schedule(static) if(gp3d_openmp_active) &
    !$omp& private(kg, x_angle, y_angle, z_angle)
    do k = 1, grid%local_nz
      do j = 1, grid%ny
        do i = 1, grid%nx
          kg = grid%k_start + k - 1
          z_angle = 2.0_dp * pi * (grid%z(kg) - z_origin) / grid%lz - pi
          y_angle = 2.0_dp * pi * (grid%y(j) - y_origin) / grid%ly - pi
          x_angle = 2.0_dp * pi * (grid%x(i) - x_origin) / grid%lx - pi
          velocity_x(i,j,k) = amplitude * sin(x_angle) * cos(y_angle) * cos(z_angle)
          velocity_y(i,j,k) = -amplitude * cos(x_angle) * sin(y_angle) * cos(z_angle)
          velocity2(i,j,k) = velocity_x(i,j,k)**2 + velocity_y(i,j,k)**2
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine fill_taylor_green_velocity

  real(dp) function wall_time_seconds() result(seconds)
    integer(int64) :: count, count_rate

    call system_clock(count=count, count_rate=count_rate)
    if (count_rate > 0_int64) then
      seconds = real(count, dp) / real(count_rate, dp)
    else
      seconds = 0.0_dp
    end if
  end function wall_time_seconds

end module gp3d_solver
