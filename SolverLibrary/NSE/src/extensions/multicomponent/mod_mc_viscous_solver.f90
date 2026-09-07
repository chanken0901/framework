module mod_mc_viscous_solver
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    compute_mc_euler_totals, compute_mc_euler_minima
  use mod_mc_euler_solver, only : write_mc_euler_csv
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_timestep, &
    advance_mc_navier_stokes_ssprk3, mc_navier_stokes_workspace
  implicit none
  private

  public :: run_mc_viscous

contains

  subroutine run_mc_viscous(config,layout,numerics)
    type(mc_config), intent(in) :: config
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: numerics
    real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), allocatable :: initial_totals(:), final_totals(:)
    real(dp) :: dt, time, conservation_error, scale
    real(dp) :: minimum_species, minimum_density, minimum_pressure
    real(dp) :: minimum_temperature
    type(mc_navier_stokes_workspace) :: workspace
    integer :: step, variable

    allocate(q(numerics%nx,numerics%ny,numerics%nz,layout%nvariables))
    allocate(q0,mold=q)
    allocate(rhs,mold=q)
    allocate(initial_totals(layout%nvariables),final_totals(layout%nvariables))
    call initialize_mc_euler_state(q,layout,numerics)
    call compute_mc_euler_totals(q,numerics,initial_totals)

    time = 0.0_dp
    write(*,'(A)') '--- stage-4 multicomponent Navier-Stokes ---'
    write(*,'(A,3(I0,1X))') 'grid = ', &
      numerics%nx,numerics%ny,numerics%nz
    write(*,'(A,I0)') 'nsteps = ', numerics%nsteps
    do step = 1, numerics%nsteps
      dt = compute_mc_navier_stokes_timestep(q,layout,numerics)
      call advance_mc_navier_stokes_ssprk3( &
        q,q0,rhs,dt,layout,numerics,workspace)
      time = time+dt
    end do

    call compute_mc_euler_totals(q,numerics,final_totals)
    conservation_error = 0.0_dp
    do variable = 1, layout%nvariables
      scale = max(abs(initial_totals(variable)),1.0_dp)
      conservation_error = max(conservation_error, &
        abs(final_totals(variable)-initial_totals(variable))/scale)
    end do
    call compute_mc_euler_minima( &
      q,layout,numerics,minimum_species,minimum_density,minimum_pressure, &
      minimum_temperature)
    if (conservation_error > 1.0e-10_dp) then
      error stop 'multicomponent Navier-Stokes conservation check failed'
    end if
    if (minimum_species < -1.0e-12_dp .or. &
        minimum_density <= 0.0_dp .or. minimum_pressure <= 0.0_dp .or. &
        minimum_temperature <= 0.0_dp) then
      error stop 'multicomponent Navier-Stokes positivity check failed'
    end if
    if (numerics%write_final) then
      call write_mc_euler_csv( &
        trim(numerics%output_file),q,config,layout,numerics)
    end if

    write(*,'(A,ES12.4)') 'final time = ', time
    write(*,'(A,ES12.4)') 'maximum relative conservation error = ', &
      conservation_error
    write(*,'(A,ES12.4)') 'minimum species partial density = ', &
      minimum_species
    write(*,'(A,ES12.4)') 'minimum mixture density = ', minimum_density
    write(*,'(A,ES12.4)') 'minimum pressure = ', minimum_pressure
    write(*,'(A,ES12.4)') 'minimum temperature = ', minimum_temperature
    write(*,'(A)') &
      'Multicomponent Navier-Stokes calculation completed successfully.'

    deallocate(q,q0,rhs,initial_totals,final_totals)
  end subroutine run_mc_viscous

end module mod_mc_viscous_solver
