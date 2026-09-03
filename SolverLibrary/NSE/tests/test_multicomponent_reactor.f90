program test_multicomponent_reactor
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, read_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_reactor_config, only : mc_reactor_config, &
    read_mc_reactor_config
  use mod_mc_reactor_solver, only : initialize_mc_reactor_state, &
    advance_mc_reactor_ssprk3, validate_mc_reactor_state
  use mod_mc_thermodynamics_provider, only : &
    configure_mc_thermodynamics, mc_mixture_density, mc_temperature
  use mod_mc_chemistry_provider, only : configure_mc_chemistry, &
    compute_mc_chemistry_source, compute_mc_chemistry_timestep
  implicit none

  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_reactor_config) :: config
  character(len=512) :: input_path
  real(dp), allocatable :: state(:), source(:)
  real(dp) :: initial_density, initial_energy, initial_temperature
  real(dp) :: dt, progress_rate, density_scale, energy_scale
  integer :: step

  if (command_argument_count() < 1) then
    error stop 'stage-5 test requires its input file path'
  end if
  call get_command_argument(1,input_path)
  call read_mc_config(trim(input_path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  call configure_mc_thermodynamics( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call configure_mc_chemistry( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call read_mc_reactor_config( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies),config)

  allocate(state(layout%nvariables),source(layout%nvariables))
  call initialize_mc_reactor_state(state,layout,config)
  call compute_mc_chemistry_source( &
    state,layout,1.4_dp,source,progress_rate)
  call assert_true(progress_rate > 0.0_dp, &
    'reactive initial state must have a positive progress rate')
  call assert_true(source(1) < 0.0_dp .and. source(2) < 0.0_dp, &
    'the one-step source must consume both reactants')
  call assert_true(source(3) > 0.0_dp, &
    'the one-step source must produce the product')
  call assert_true(abs(sum(source(1:model%nspecies))) < 1.0e-14_dp, &
    'the chemistry source must conserve total species mass')
  call assert_true(maxval(abs(source(layout%momentum))) < tiny(1.0_dp), &
    'homogeneous chemistry must not create momentum')
  call assert_true(abs(source(layout%total_energy)) < tiny(1.0_dp), &
    'adiabatic chemistry must conserve total energy')

  initial_density = mc_mixture_density(state,layout)
  initial_energy = state(layout%total_energy)
  initial_temperature = mc_temperature(state,layout,1.4_dp)
  density_scale = max(abs(initial_density),1.0_dp)
  energy_scale = max(abs(initial_energy),1.0_dp)
  do step = 1, 100
    dt = compute_mc_chemistry_timestep( &
      state,layout,1.4_dp,config%chemistry_cfl,config%maximum_dt)
    call assert_true(dt > 0.0_dp .and. dt <= config%maximum_dt, &
      'adaptive chemistry timestep must be positive and bounded')
    call advance_mc_reactor_ssprk3(state,dt,layout)
  end do
  call validate_mc_reactor_state(state,layout)
  call assert_true(state(1) < 0.5_dp*initial_density, &
    'finite-rate integration must consume fuel')
  call assert_true(state(3) > 0.0_dp, &
    'finite-rate integration must create product')
  call assert_true(mc_temperature(state,layout,1.4_dp) > &
    initial_temperature, &
    'exothermic formation energy must raise the reactor temperature')
  call assert_true(abs(mc_mixture_density(state,layout)-initial_density) < &
    1.0e-11_dp*density_scale, &
    'reactor integration must conserve total mass')
  call assert_true(abs(state(layout%total_energy)-initial_energy) < &
    1.0e-13_dp*energy_scale, &
    'reactor integration must conserve total energy')

  write(*,'(A)') 'Homogeneous finite-rate chemistry tests passed.'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'homogeneous finite-rate chemistry test failure'
    end if
  end subroutine assert_true

end program test_multicomponent_reactor
