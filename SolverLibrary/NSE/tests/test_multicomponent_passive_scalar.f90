program test_multicomponent_passive_scalar
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, initialize_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_passive_scalar_config, only : mc_passive_scalar_config, &
    initialize_mc_passive_scalar_config, validate_mc_passive_scalar_config, &
    mc_passive_scalar_timestep
  use mod_mc_passive_scalar_field, only : &
    initialize_mc_passive_scalar_state, compute_mc_species_masses, &
    compute_mc_species_sum_error, compute_mc_tracer_bounds
  use mod_mc_passive_scalar_advection, only : &
    compute_mc_passive_scalar_rhs, advance_mc_passive_scalar_ssprk3
  implicit none

  type(mc_config) :: config
  type(mc_state_layout) :: layout
  type(mc_passive_scalar_config) :: passive
  real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
  real(dp) :: initial_mass(2), final_mass(2), dt
  real(dp) :: tracer_minimum, tracer_maximum
  integer :: step

  call initialize_mc_config(config)
  config%nspecies = 2
  config%species_names = ''
  config%species_names(1) = 'tracer'
  config%species_names(2) = 'carrier'
  config%simulation_mode = 'passive_scalar'
  call initialize_mc_state_layout(layout, config%nspecies)

  call initialize_mc_passive_scalar_config(passive)
  passive%nx = 24
  passive%ny = 12
  passive%nz = 8
  passive%velocity = [0.7_dp, -0.2_dp, 0.1_dp]
  passive%cfl = 0.45_dp
  passive%nsteps = 25
  passive%write_final = .false.
  call validate_mc_passive_scalar_config(passive)
  dt = mc_passive_scalar_timestep(passive)

  allocate(q(passive%nx,passive%ny,passive%nz,layout%nvariables))
  allocate(q0, mold=q)
  allocate(rhs, mold=q)
  call initialize_mc_passive_scalar_state(q, layout, passive)
  call compute_mc_species_masses(q, layout, passive, initial_mass)
  do step = 1, passive%nsteps
    call advance_mc_passive_scalar_ssprk3(q, q0, rhs, dt, layout, passive)
  end do
  call compute_mc_species_masses(q, layout, passive, final_mass)

  call assert_true(maxval(abs(final_mass-initial_mass)) < 2.0e-13_dp, &
    'periodic advection must conserve every species mass')
  call assert_true(compute_mc_species_sum_error(q, layout) < 2.0e-13_dp, &
    'species partial densities must continue to sum to one')
  call compute_mc_tracer_bounds(q, layout, tracer_minimum, tracer_maximum)
  call assert_true(tracer_minimum >= -1.0e-14_dp .and. &
    tracer_maximum <= 1.0_dp+1.0e-14_dp, &
    'upwind SSPRK3 transport must preserve tracer bounds')
  call assert_true(maxval(abs(q(:,:,:,layout%momentum(1)) - &
    passive%velocity(1))) < epsilon(1.0_dp), &
    'passive transport must not change carrier momentum')
  call assert_true(maxval(abs(q(:,:,:,layout%total_energy) - &
    (1.0_dp+0.5_dp*sum(passive%velocity**2)))) < epsilon(1.0_dp), &
    'passive transport must not change carrier energy')

  q(:,:,:,layout%first_species) = 0.25_dp
  q(:,:,:,layout%last_species) = 0.75_dp
  call compute_mc_passive_scalar_rhs(q, rhs, layout, passive)
  call assert_true(maxval(abs(rhs)) < epsilon(1.0_dp), &
    'uniform species fields must have exactly zero periodic RHS')

  write(*,'(A)') 'Multicomponent passive-scalar tests passed.'

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'multicomponent passive-scalar test failure'
    end if
  end subroutine assert_true

end program test_multicomponent_passive_scalar
