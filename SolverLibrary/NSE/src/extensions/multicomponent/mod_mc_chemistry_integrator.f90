module mod_mc_chemistry_integrator
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_pressure, mc_temperature
  use mod_mc_chemistry_provider, only : compute_mc_chemistry_source, &
    compute_mc_chemistry_timestep
  implicit none
  private

  real(dp), parameter :: state_tolerance = 1.0e-12_dp

  public :: validate_mc_chemistry_state
  public :: advance_mc_chemistry_ssprk3
  public :: advance_mc_chemistry_interval

contains

  subroutine validate_mc_chemistry_state(state,layout,gamma)
    real(dp), intent(in) :: state(:), gamma
    type(mc_state_layout), intent(in) :: layout
    real(dp) :: density, pressure, temperature

    if (size(state) /= layout%nvariables .or. &
        .not. all(ieee_is_finite(state))) then
      error stop 'chemistry state does not match the state layout'
    end if
    if (minval(state(layout%first_species:layout%last_species)) < &
        -state_tolerance) then
      error stop 'negative species partial density in chemistry state'
    end if
    density = mc_mixture_density(state,layout)
    if (.not. ieee_is_finite(density) .or. &
        density <= state_tolerance) then
      error stop 'non-positive density in chemistry state'
    end if
    temperature = mc_temperature(state,layout,gamma)
    pressure = mc_pressure(state,layout,gamma)
    if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_dp .or. &
        .not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
      error stop 'non-positive thermodynamic chemistry state'
    end if
  end subroutine validate_mc_chemistry_state

  subroutine advance_mc_chemistry_ssprk3(state,dt,layout,gamma)
    real(dp), intent(inout) :: state(:)
    real(dp), intent(in) :: dt, gamma
    type(mc_state_layout), intent(in) :: layout
    real(dp) :: initial(layout%nvariables)
    real(dp) :: stage(layout%nvariables)
    real(dp) :: source(layout%nvariables)

    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      error stop 'chemistry SSPRK3 dt must be finite and positive'
    end if
    initial = state
    call compute_mc_chemistry_source(initial,layout,gamma,source)
    stage = initial+dt*source
    call validate_mc_chemistry_state(stage,layout,gamma)

    call compute_mc_chemistry_source(stage,layout,gamma,source)
    stage = 0.75_dp*initial+0.25_dp*(stage+dt*source)
    call validate_mc_chemistry_state(stage,layout,gamma)

    call compute_mc_chemistry_source(stage,layout,gamma,source)
    state = initial/3.0_dp+2.0_dp*(stage+dt*source)/3.0_dp
    call validate_mc_chemistry_state(state,layout,gamma)
  end subroutine advance_mc_chemistry_ssprk3

  subroutine advance_mc_chemistry_interval( &
      state,duration,layout,gamma,chemistry_cfl,maximum_substeps, &
      substeps_used)
    real(dp), intent(inout) :: state(:)
    real(dp), intent(in) :: duration, gamma, chemistry_cfl
    type(mc_state_layout), intent(in) :: layout
    integer, intent(in) :: maximum_substeps
    integer, intent(out), optional :: substeps_used
    real(dp) :: remaining, stable_dt, substep_dt, time_tolerance
    integer :: substeps

    if (.not. ieee_is_finite(duration) .or. duration < 0.0_dp) then
      error stop 'chemistry interval must be finite and nonnegative'
    end if
    if (.not. ieee_is_finite(chemistry_cfl) .or. &
        chemistry_cfl <= 0.0_dp .or. chemistry_cfl > 1.0_dp) then
      error stop 'chemistry CFL must be in (0,1]'
    end if
    if (maximum_substeps < 1) then
      error stop 'maximum chemistry substeps must be positive'
    end if
    call validate_mc_chemistry_state(state,layout,gamma)

    remaining = duration
    substeps = 0
    time_tolerance = max( &
      16.0_dp*epsilon(1.0_dp)*duration,tiny(1.0_dp))
    do while (remaining > 0.0_dp)
      if (substeps >= maximum_substeps) then
        error stop 'maximum chemistry substeps exceeded'
      end if
      stable_dt = compute_mc_chemistry_timestep( &
        state,layout,gamma,chemistry_cfl,remaining)
      substep_dt = min(remaining,stable_dt)
      call advance_mc_chemistry_ssprk3( &
        state,substep_dt,layout,gamma)
      substeps = substeps+1
      remaining = remaining-substep_dt
      if (remaining <= time_tolerance) remaining = 0.0_dp
    end do
    if (present(substeps_used)) substeps_used = substeps
  end subroutine advance_mc_chemistry_interval

end module mod_mc_chemistry_integrator
