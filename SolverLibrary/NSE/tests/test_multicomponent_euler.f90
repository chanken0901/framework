program test_multicomponent_euler
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, &
    initialize_mc_euler_config, validate_mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    set_mc_euler_conservative_state, validate_mc_euler_state, &
    compute_mc_euler_totals, compute_mc_euler_minima
  use mod_mc_euler_flux, only : compute_mc_euler_physical_flux, &
    compute_mc_euler_rusanov_flux, compute_mc_euler_rhs, &
    compute_mc_euler_timestep, advance_mc_euler_ssprk3
  implicit none

  type(mc_state_layout) :: layout
  type(mc_state_layout) :: single_layout
  type(mc_euler_config) :: config
  real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
  real(dp), allocatable :: initial_totals(:), final_totals(:)
  real(dp), allocatable :: state(:), physical_flux(:), numerical_flux(:)
  real(dp), allocatable :: single_state(:), single_flux(:)
  real(dp) :: mass_fractions(2), velocity(3), dt
  real(dp) :: minimum_species, minimum_density, minimum_pressure
  integer :: i, j, k, step

  call initialize_mc_state_layout(layout,2)
  call initialize_mc_euler_config(config,layout%nspecies)
  config%nx = 24
  config%ny = 6
  config%nz = 4
  config%cfl = 0.35_dp
  config%nsteps = 8
  config%write_final = .false.
  call validate_mc_euler_config(config,layout%nspecies)

  allocate(state(layout%nvariables),physical_flux(layout%nvariables), &
    numerical_flux(layout%nvariables))
  mass_fractions = [0.4_dp,0.6_dp]
  velocity = [0.3_dp,-0.1_dp,0.2_dp]
  call set_mc_euler_conservative_state( &
    state,layout,config%gamma,1.2_dp,velocity,0.9_dp,mass_fractions)
  call compute_mc_euler_physical_flux( &
    state,layout,config%gamma,1,physical_flux)
  call compute_mc_euler_rusanov_flux( &
    state,state,layout,config%gamma,1,numerical_flux)
  call assert_true(maxval(abs(physical_flux-numerical_flux)) < &
    10.0_dp*epsilon(1.0_dp), &
    'Rusanov flux must equal physical flux for identical states')
  call assert_true(abs(sum(physical_flux( &
    layout%first_species:layout%last_species))-1.2_dp*velocity(1)) < &
    10.0_dp*epsilon(1.0_dp), &
    'species flux sum must equal total mass flux')

  allocate(q(config%nx,config%ny,config%nz,layout%nvariables))
  allocate(q0,mold=q)
  allocate(rhs,mold=q)
  allocate(initial_totals(layout%nvariables),final_totals(layout%nvariables))
  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        q(i,j,k,:) = state
      end do
    end do
  end do
  call compute_mc_euler_rhs(q,rhs,layout,config)
  call assert_true(maxval(abs(rhs)) < 10.0_dp*epsilon(1.0_dp), &
    'uniform periodic Euler state must have zero RHS')

  call initialize_mc_euler_state(q,layout,config)
  call compute_mc_euler_totals(q,config,initial_totals)
  do step = 1, config%nsteps
    dt = compute_mc_euler_timestep(q,layout,config)
    call advance_mc_euler_ssprk3(q,q0,rhs,dt,layout,config)
  end do
  call validate_mc_euler_state(q,layout,config)
  call compute_mc_euler_totals(q,config,final_totals)
  call assert_true(maxval(abs(final_totals-initial_totals)) < 5.0e-13_dp, &
    'periodic Euler update must conserve every state variable')
  call compute_mc_euler_minima( &
    q,layout,config,minimum_species,minimum_density,minimum_pressure)
  call assert_true(minimum_species >= -1.0e-13_dp, &
    'Euler update must preserve non-negative partial densities')
  call assert_true(minimum_density > 0.0_dp .and. minimum_pressure > 0.0_dp, &
    'Euler update must preserve positive density and pressure')

  call initialize_mc_state_layout(single_layout,1)
  allocate(single_state(single_layout%nvariables))
  allocate(single_flux(single_layout%nvariables))
  call set_mc_euler_conservative_state( &
    single_state,single_layout,config%gamma,1.1_dp, &
    [0.25_dp,-0.05_dp,0.1_dp],0.8_dp,[1.0_dp])
  call compute_mc_euler_physical_flux( &
    single_state,single_layout,config%gamma,1,single_flux)
  call assert_true(single_layout%nvariables == 5, &
    'one-species Euler limit must have five variables')
  call assert_true(abs(single_flux(single_layout%first_species) - &
    1.1_dp*0.25_dp) < 10.0_dp*epsilon(1.0_dp), &
    'one-species Euler mass flux is inconsistent')
  call assert_true(abs(single_flux(single_layout%momentum(1)) - &
    (1.1_dp*0.25_dp**2+0.8_dp)) < 10.0_dp*epsilon(1.0_dp), &
    'one-species Euler momentum flux is inconsistent')

  write(*,'(A)') 'Multicomponent Euler tests passed.'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'multicomponent Euler test failure'
    end if
  end subroutine assert_true

end program test_multicomponent_euler
