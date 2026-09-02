program test_multicomponent_thermally_perfect
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, read_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, read_mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    set_mc_euler_conservative_state, validate_mc_euler_state, &
    compute_mc_euler_totals
  use mod_mc_euler_flux, only : compute_mc_euler_physical_flux, &
    compute_mc_euler_rusanov_flux, compute_mc_euler_rhs, &
    compute_mc_euler_timestep, advance_mc_euler_ssprk3
  use mod_mc_thermodynamics_provider, only : &
    configure_mc_thermodynamics, mc_mixture_gas_constant, &
    mc_mixture_gamma, mc_pressure, mc_temperature
  implicit none

  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: config
  character(len=512) :: input_path
  real(dp), allocatable :: state(:), physical_flux(:), numerical_flux(:)
  real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
  real(dp), allocatable :: initial_totals(:), final_totals(:)
  real(dp) :: fractions(2), velocity(3), gas_constant
  real(dp) :: gamma_300, gamma_2500, pressure, temperature, dt, scale
  integer :: i, j, k, step, variable

  if (command_argument_count() < 1) then
    error stop 'thermally-perfect test requires its input file path'
  end if
  call get_command_argument(1,input_path)
  call read_mc_config(trim(input_path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  call configure_mc_thermodynamics( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call read_mc_euler_config(trim(input_path),model%nspecies,config)
  config%write_final = .false.

  fractions = [0.767_dp,0.233_dp]
  velocity = [35.0_dp,-4.0_dp,2.0_dp]
  gas_constant = mc_mixture_gas_constant(fractions,layout)
  allocate(state(layout%nvariables),physical_flux(layout%nvariables), &
    numerical_flux(layout%nvariables))
  call set_mc_euler_conservative_state( &
    state,layout,config%gamma,1.17197031944841_dp,velocity, &
    101325.0_dp,fractions)
  temperature = mc_temperature(state,layout,config%gamma)
  pressure = mc_pressure(state,layout,config%gamma)
  call assert_close(temperature,300.0_dp,1.0e-7_dp, &
    'primitive-to-conservative temperature recovery failed')
  call assert_close(pressure,101325.0_dp,1.0e-4_dp, &
    'thermally-perfect pressure recovery failed')
  call assert_close( &
    101325.0_dp/(1.17197031944841_dp*gas_constant),300.0_dp, &
    2.0e-8_dp,'mixture gas constant is inconsistent')

  gamma_300 = mc_mixture_gamma( &
    fractions,layout,300.0_dp,config%gamma)
  gamma_2500 = mc_mixture_gamma( &
    fractions,layout,2500.0_dp,config%gamma)
  call assert_true(gamma_300 > 1.0_dp .and. gamma_2500 > 1.0_dp, &
    'mixture gamma must exceed one')
  call assert_true(abs(gamma_300-gamma_2500) > 1.0e-3_dp, &
    'NASA-7 mixture gamma must depend on temperature')

  call compute_mc_euler_physical_flux( &
    state,layout,config%gamma,1,physical_flux)
  call compute_mc_euler_rusanov_flux( &
    state,state,layout,config%gamma,1,numerical_flux)
  call assert_true(maxval(abs(physical_flux-numerical_flux)) < &
    1.0e-10_dp*max(maxval(abs(physical_flux)),1.0_dp), &
    'identical-state Rusanov flux must equal the physical flux')

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
  call assert_true(maxval(abs(rhs)) < &
    1.0e-10_dp*max(maxval(abs(state)),1.0_dp), &
    'uniform thermally-perfect state must have zero periodic RHS')

  call initialize_mc_euler_state(q,layout,config)
  call compute_mc_euler_totals(q,config,initial_totals)
  do step = 1, config%nsteps
    dt = compute_mc_euler_timestep(q,layout,config)
    call advance_mc_euler_ssprk3(q,q0,rhs,dt,layout,config)
  end do
  call validate_mc_euler_state(q,layout,config)
  call compute_mc_euler_totals(q,config,final_totals)
  do variable = 1, layout%nvariables
    scale = max(abs(initial_totals(variable)),1.0_dp)
    call assert_true(abs(final_totals(variable)-initial_totals(variable)) < &
      2.0e-10_dp*scale, &
      'periodic stage-3 update must conserve every state variable')
  end do

  write(*,'(A)') 'Thermally-perfect multicomponent tests passed.'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'thermally-perfect multicomponent test failure'
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,tolerance,message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    if (abs(actual-expected) > tolerance) then
      write(*,'(A,ES24.16)') 'actual = ', actual
      write(*,'(A,ES24.16)') 'expected = ', expected
      write(*,'(A,ES24.16)') 'tolerance = ', tolerance
    end if
    call assert_true(abs(actual-expected) <= tolerance,message)
  end subroutine assert_close

end program test_multicomponent_thermally_perfect
