program main_nse_multicomponent_inviscid
  use mod_mc_config, only : mc_config, read_mc_config, print_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_provider_registry, only : validate_mc_providers, &
    print_mc_providers
  use mod_mc_euler_config, only : mc_euler_config, read_mc_euler_config
  use mod_mc_euler_solver, only : run_mc_euler
  use mod_mc_thermodynamics_provider, only : &
    configure_mc_thermodynamics, mc_thermodynamics_provider_name
  implicit none

  type(mc_config) :: config
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: euler
  character(len=512) :: input_path
  logical :: input_exists

  input_path = 'input.dat'
  if (command_argument_count() >= 1) call get_command_argument(1,input_path)
  inquire(file=trim(input_path),exist=input_exists)
  if (.not. input_exists) then
    write(*,'(A,A)') 'ERROR: multicomponent Euler input does not exist: ', &
      trim(input_path)
    error stop 'missing multicomponent Euler input'
  end if

  call read_mc_config(trim(input_path),config)
  if (trim(config%simulation_mode) /= 'inviscid_euler' .and. &
      trim(config%simulation_mode) /= 'thermally_perfect_euler') then
    error stop 'Euler executable requires an Euler simulation mode'
  end if
  call read_mc_euler_config(trim(input_path),config%nspecies,euler)
  call initialize_mc_state_layout(layout,config%nspecies)
  call validate_mc_providers(config)
  call configure_mc_thermodynamics( &
    trim(input_path),config%nspecies, &
    config%species_names(1:config%nspecies))
  if (trim(config%simulation_mode) == 'thermally_perfect_euler' .and. &
      mc_thermodynamics_provider_name /= 'thermally_perfect') then
    error stop 'stage-3 mode requires thermally-perfect thermodynamics'
  end if
  call print_mc_config(config)
  call print_mc_providers()
  write(*,'(A,I0)') 'conservative variables = ', layout%nvariables
  call run_mc_euler(config,layout,euler)
end program main_nse_multicomponent_inviscid
