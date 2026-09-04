module mod_mc_reactive_solver
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    validate_mc_euler_state, compute_mc_euler_totals, &
    compute_mc_euler_minima
  use mod_mc_euler_solver, only : write_mc_euler_csv
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_timestep, &
    advance_mc_navier_stokes_ssprk3
  use mod_mc_chemistry_provider, only : compute_mc_chemistry_timestep
  use mod_mc_chemistry_integrator, only : advance_mc_chemistry_interval
  use mod_mc_reactive_config, only : mc_reactive_config
  implicit none
  private

  public :: compute_mc_reactive_timestep
  public :: apply_mc_field_chemistry
  public :: advance_mc_reactive_strang
  public :: run_mc_reactive

contains

  real(dp) function compute_mc_reactive_timestep( &
      q,layout,numerics,reactive) result(dt)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: numerics
    type(mc_reactive_config), intent(in) :: reactive
    type(mc_euler_config) :: automatic_numerics
    real(dp) :: fluid_dt, chemistry_dt
    integer :: i, j, k

    if (numerics%dt > 0.0_dp) then
      dt = compute_mc_navier_stokes_timestep(q,layout,numerics)
      return
    end if

    automatic_numerics = numerics
    automatic_numerics%dt = 0.0_dp
    fluid_dt = compute_mc_navier_stokes_timestep( &
      q,layout,automatic_numerics)
    chemistry_dt = fluid_dt
    do k = 1, size(q,3)
      do j = 1, size(q,2)
        do i = 1, size(q,1)
          chemistry_dt = min(chemistry_dt, &
            compute_mc_chemistry_timestep( &
            q(i,j,k,:),layout,numerics%gamma,reactive%chemistry_cfl, &
            fluid_dt))
        end do
      end do
    end do
    dt = min(fluid_dt,chemistry_dt)
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      error stop 'reactive Navier-Stokes timestep is invalid'
    end if
  end function compute_mc_reactive_timestep

  subroutine apply_mc_field_chemistry( &
      q,duration,layout,numerics,reactive,maximum_substeps_used)
    real(dp), intent(inout) :: q(:,:,:,:)
    real(dp), intent(in) :: duration
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: numerics
    type(mc_reactive_config), intent(in) :: reactive
    integer, intent(out), optional :: maximum_substeps_used
    integer :: i, j, k, substeps, maximum_used

    maximum_used = 0
    do k = 1, size(q,3)
      do j = 1, size(q,2)
        do i = 1, size(q,1)
          call advance_mc_chemistry_interval( &
            q(i,j,k,:),duration,layout,numerics%gamma, &
            reactive%chemistry_cfl, &
            reactive%maximum_chemistry_substeps,substeps)
          maximum_used = max(maximum_used,substeps)
        end do
      end do
    end do
    call validate_mc_euler_state(q,layout,numerics)
    if (present(maximum_substeps_used)) then
      maximum_substeps_used = maximum_used
    end if
  end subroutine apply_mc_field_chemistry

  subroutine advance_mc_reactive_strang( &
      q,q0,rhs,dt,layout,numerics,reactive,maximum_substeps_used)
    real(dp), intent(inout) :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: numerics
    type(mc_reactive_config), intent(in) :: reactive
    integer, intent(out), optional :: maximum_substeps_used
    type(mc_euler_config) :: fixed_step_numerics
    real(dp) :: checked_dt
    integer :: first_half_substeps, second_half_substeps

    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      error stop 'reactive Strang timestep must be finite and positive'
    end if
    call apply_mc_field_chemistry( &
      q,0.5_dp*dt,layout,numerics,reactive,first_half_substeps)

    fixed_step_numerics = numerics
    fixed_step_numerics%dt = dt
    checked_dt = compute_mc_navier_stokes_timestep( &
      q,layout,fixed_step_numerics)
    call advance_mc_navier_stokes_ssprk3( &
      q,q0,rhs,checked_dt,layout,numerics)

    call apply_mc_field_chemistry( &
      q,0.5_dp*dt,layout,numerics,reactive,second_half_substeps)
    if (present(maximum_substeps_used)) then
      maximum_substeps_used = max( &
        first_half_substeps,second_half_substeps)
    end if
  end subroutine advance_mc_reactive_strang

  subroutine run_mc_reactive(model,layout,numerics,reactive)
    type(mc_config), intent(in) :: model
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: numerics
    type(mc_reactive_config), intent(in) :: reactive
    real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), allocatable :: initial_totals(:), final_totals(:)
    real(dp) :: dt, time, initial_mass, final_mass, scale
    real(dp) :: conservation_error, minimum_species, minimum_density
    real(dp) :: minimum_pressure, minimum_temperature
    integer :: step, direction, step_substeps, maximum_substeps_used

    allocate(q(numerics%nx,numerics%ny,numerics%nz,layout%nvariables))
    allocate(q0,mold=q)
    allocate(rhs,mold=q)
    allocate(initial_totals(layout%nvariables),final_totals(layout%nvariables))
    call initialize_mc_euler_state(q,layout,numerics)
    call compute_mc_euler_totals(q,numerics,initial_totals)

    time = 0.0_dp
    maximum_substeps_used = 0
    write(*,'(A)') '--- stage-6 reactive multicomponent Navier-Stokes ---'
    write(*,'(A,3(I0,1X))') 'grid = ', &
      numerics%nx,numerics%ny,numerics%nz
    write(*,'(A,I0)') 'nsteps = ', numerics%nsteps
    do step = 1, numerics%nsteps
      dt = compute_mc_reactive_timestep(q,layout,numerics,reactive)
      call advance_mc_reactive_strang( &
        q,q0,rhs,dt,layout,numerics,reactive,step_substeps)
      maximum_substeps_used = max(maximum_substeps_used,step_substeps)
      time = time+dt
    end do

    call compute_mc_euler_totals(q,numerics,final_totals)
    initial_mass = sum(initial_totals( &
      layout%first_species:layout%last_species))
    final_mass = sum(final_totals( &
      layout%first_species:layout%last_species))
    conservation_error = abs(final_mass-initial_mass)/ &
      max(abs(initial_mass),1.0_dp)
    do direction = 1, 3
      scale = max(abs(initial_totals(layout%momentum(direction))),1.0_dp)
      conservation_error = max(conservation_error, &
        abs(final_totals(layout%momentum(direction))- &
        initial_totals(layout%momentum(direction)))/scale)
    end do
    scale = max(abs(initial_totals(layout%total_energy)),1.0_dp)
    conservation_error = max(conservation_error, &
      abs(final_totals(layout%total_energy)- &
      initial_totals(layout%total_energy))/scale)
    call compute_mc_euler_minima( &
      q,layout,numerics,minimum_species,minimum_density,minimum_pressure, &
      minimum_temperature)
    if (conservation_error > 1.0e-10_dp) then
      error stop 'reactive Navier-Stokes conservation check failed'
    end if
    if (minimum_species < -1.0e-12_dp .or. &
        minimum_density <= 0.0_dp .or. minimum_pressure <= 0.0_dp .or. &
        minimum_temperature <= 0.0_dp) then
      error stop 'reactive Navier-Stokes positivity check failed'
    end if
    if (numerics%write_final) then
      call write_mc_euler_csv( &
        trim(numerics%output_file),q,model,layout,numerics)
    end if

    write(*,'(A,ES12.4)') 'final time = ', time
    write(*,'(A,I0)') 'maximum chemistry substeps per half-step = ', &
      maximum_substeps_used
    write(*,'(A,ES12.4)') 'maximum relative conservation error = ', &
      conservation_error
    write(*,'(A,ES12.4)') 'minimum species partial density = ', &
      minimum_species
    write(*,'(A,ES12.4)') 'minimum mixture density = ', minimum_density
    write(*,'(A,ES12.4)') 'minimum pressure = ', minimum_pressure
    write(*,'(A,ES12.4)') 'minimum temperature = ', minimum_temperature
    write(*,'(A)') &
      'Reactive multicomponent Navier-Stokes calculation completed successfully.'

    deallocate(q,q0,rhs,initial_totals,final_totals)
  end subroutine run_mc_reactive

end module mod_mc_reactive_solver
