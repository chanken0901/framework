program test_multicomponent_reactive
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, read_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, read_mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    validate_mc_euler_state, compute_mc_euler_totals
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_timestep, &
    advance_mc_navier_stokes_ssprk3
  use mod_mc_reactive_config, only : mc_reactive_config, &
    read_mc_reactive_config
  use mod_mc_reactive_solver, only : compute_mc_reactive_timestep, &
    advance_mc_reactive_strang
  use mod_mc_chemistry_integrator, only : advance_mc_chemistry_interval
  use mod_mc_thermodynamics_provider, only : configure_mc_thermodynamics
  use mod_mc_transport_provider, only : configure_mc_transport
  use mod_mc_chemistry_provider, only : configure_mc_chemistry
  implicit none

  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: numerics, test_numerics
  type(mc_reactive_config) :: reactive
  character(len=512) :: input_path
  real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
  real(dp), allocatable :: reference(:,:,:,:), initial_totals(:)
  real(dp), allocatable :: final_totals(:), state_reference(:)
  real(dp) :: dt, initial_mass, final_mass, scale, error
  real(dp) :: initial_fuel, initial_product, final_fuel, final_product
  integer :: i, j, k, direction, step, substeps, maximum_substeps

  if (command_argument_count() < 1) then
    error stop 'stage-6 test requires its input file path'
  end if
  call get_command_argument(1,input_path)
  call read_mc_config(trim(input_path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  call configure_mc_thermodynamics( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call configure_mc_transport( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call configure_mc_chemistry( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call read_mc_euler_config(trim(input_path),model%nspecies,numerics)
  call read_mc_reactive_config(trim(input_path),reactive)
  numerics%write_final = .false.

  allocate(q(numerics%nx,numerics%ny,numerics%nz,layout%nvariables))
  allocate(q0,mold=q)
  allocate(rhs,mold=q)
  allocate(reference,mold=q)
  allocate(initial_totals(layout%nvariables),final_totals(layout%nvariables))
  allocate(state_reference(layout%nvariables))

  call initialize_mc_euler_state(q,layout,numerics)
  call compute_mc_euler_totals(q,numerics,initial_totals)
  initial_mass = sum(initial_totals( &
    layout%first_species:layout%last_species))
  initial_fuel = initial_totals(layout%first_species)
  initial_product = initial_totals(layout%first_species+2)
  maximum_substeps = 0
  do step = 1, 4
    dt = compute_mc_reactive_timestep(q,layout,numerics,reactive)
    call assert_true(dt > 0.0_dp, &
      'reactive timestep must be positive')
    call advance_mc_reactive_strang( &
      q,q0,rhs,dt,layout,numerics,reactive,substeps)
    maximum_substeps = max(maximum_substeps,substeps)
  end do
  call validate_mc_euler_state(q,layout,numerics)
  call compute_mc_euler_totals(q,numerics,final_totals)
  final_mass = sum(final_totals( &
    layout%first_species:layout%last_species))
  final_fuel = final_totals(layout%first_species)
  final_product = final_totals(layout%first_species+2)
  call assert_true(final_fuel < initial_fuel, &
    'reactive flow must consume fuel')
  call assert_true(final_product > initial_product, &
    'reactive flow must create product')
  call assert_true(maximum_substeps >= 1, &
    'each reactive half-step must advance chemistry')
  call assert_close(final_mass,initial_mass,1.0e-11_dp, &
    'reactive flow must conserve total mass')
  do direction = 1, 3
    call assert_close( &
      final_totals(layout%momentum(direction)), &
      initial_totals(layout%momentum(direction)),1.0e-11_dp, &
      'reactive flow must conserve momentum')
  end do
  call assert_close( &
    final_totals(layout%total_energy), &
    initial_totals(layout%total_energy),1.0e-11_dp, &
    'reactive flow must conserve total energy')

  test_numerics = numerics
  test_numerics%wave_mean_mass_fractions(1:3) = &
    [0.5_dp,0.5_dp,0.0_dp]
  test_numerics%wave_amplitude = 0.0_dp
  call initialize_mc_euler_state(q,layout,test_numerics)
  state_reference = q(1,1,1,:)
  dt = compute_mc_reactive_timestep(q,layout,test_numerics,reactive)
  call advance_mc_chemistry_interval( &
    state_reference,0.5_dp*dt,layout,test_numerics%gamma, &
    reactive%chemistry_cfl,reactive%maximum_chemistry_substeps)
  call advance_mc_chemistry_interval( &
    state_reference,0.5_dp*dt,layout,test_numerics%gamma, &
    reactive%chemistry_cfl,reactive%maximum_chemistry_substeps)
  call advance_mc_reactive_strang( &
    q,q0,rhs,dt,layout,test_numerics,reactive)
  error = 0.0_dp
  scale = max(maxval(abs(state_reference)),1.0_dp)
  do k = 1, test_numerics%nz
    do j = 1, test_numerics%ny
      do i = 1, test_numerics%nx
        error = max(error,maxval(abs(q(i,j,k,:)-state_reference))/scale)
      end do
    end do
  end do
  call assert_true(error < 2.0e-12_dp, &
    'uniform reactive flow must reduce to homogeneous chemistry')

  test_numerics%wave_mean_mass_fractions(1:3) = &
    [0.5_dp,0.0_dp,0.5_dp]
  test_numerics%wave_positive_species = 1
  test_numerics%wave_negative_species = 3
  test_numerics%wave_amplitude = 0.1_dp
  call initialize_mc_euler_state(q,layout,test_numerics)
  reference = q
  dt = compute_mc_navier_stokes_timestep( &
    reference,layout,test_numerics)
  call advance_mc_navier_stokes_ssprk3( &
    reference,q0,rhs,dt,layout,test_numerics)
  call advance_mc_reactive_strang( &
    q,q0,rhs,dt,layout,test_numerics,reactive)
  scale = max(maxval(abs(reference)),1.0_dp)
  error = maxval(abs(q-reference))/scale
  call assert_true(error < 2.0e-12_dp, &
    'zero-rate reactive flow must reduce to stage-4 transport')

  write(*,'(A)') &
    'Reactive multicomponent Strang-splitting tests passed.'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'reactive multicomponent test failure'
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,tolerance,message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    real(dp) :: comparison_scale

    comparison_scale = max(abs(expected),1.0_dp)
    call assert_true(abs(actual-expected) <= tolerance*comparison_scale, &
      message)
  end subroutine assert_close

end program test_multicomponent_reactive
