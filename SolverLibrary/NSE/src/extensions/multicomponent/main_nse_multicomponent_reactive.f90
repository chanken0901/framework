program main_nse_multicomponent_reactive
  use mod_mc_config, only : mc_config, read_mc_config, print_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_provider_registry, only : validate_mc_providers, &
    print_mc_providers
  use mod_mc_euler_config, only : mc_euler_config, read_mc_euler_config
  use mod_mc_reactive_config, only : mc_reactive_config, &
    read_mc_reactive_config
  use mod_mc_reactive_solver, only : run_mc_reactive
  use mod_mc_thermodynamics_provider, only : &
    configure_mc_thermodynamics, mc_thermodynamics_provider_name
  use mod_mc_transport_provider, only : configure_mc_transport, &
    mc_transport_provider_name
  use mod_mc_chemistry_provider, only : configure_mc_chemistry, &
    mc_chemistry_provider_name
  implicit none

  type(mc_config) :: config
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: numerics
  type(mc_reactive_config) :: reactive
  character(len=512) :: input_path
  logical :: input_exists

  input_path = 'input.dat'
  if (command_argument_count() >= 1) call get_command_argument(1,input_path)
  inquire(file=trim(input_path),exist=input_exists)
  if (.not. input_exists) then
    write(*,'(A,A)') 'ERROR: reactive input does not exist: ', &
      trim(input_path)
    error stop 'missing reactive Navier-Stokes input'
  end if

  call read_mc_config(trim(input_path),config)
  if (trim(config%simulation_mode) /= 'reactive_navier_stokes') then
    error stop 'reactive executable requires mode=reactive_navier_stokes'
  end if
  call read_mc_euler_config(trim(input_path),config%nspecies,numerics)
  if (trim(numerics%initial_condition) /= 'periodic_species_wave_x' .and. &
      trim(numerics%initial_condition) /= 'reactive_shock_tube_x') then
    error stop 'reactive executable received unsupported initial condition'
  end if
  call read_mc_reactive_config(trim(input_path),reactive)
  call initialize_mc_state_layout(layout,config%nspecies)
  call validate_mc_providers(config)
  call configure_mc_thermodynamics( &
    trim(input_path),config%nspecies, &
    config%species_names(1:config%nspecies))
  call configure_mc_transport( &
    trim(input_path),config%nspecies, &
    config%species_names(1:config%nspecies))
  call configure_mc_chemistry( &
    trim(input_path),config%nspecies, &
    config%species_names(1:config%nspecies))
  if (mc_thermodynamics_provider_name /= 'thermally_perfect' .or. &
      mc_transport_provider_name /= 'mixture_averaged' .or. &
      mc_chemistry_provider_name /= 'one_step_arrhenius') then
    error stop 'reactive executable contains incompatible providers'
  end if
  call print_mc_config(config)
  call print_mc_providers()
  write(*,'(A,I0)') 'conservative variables = ', layout%nvariables
  call run_mc_reactive(config,layout,numerics,reactive)
end program main_nse_multicomponent_reactive
