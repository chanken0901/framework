program main_nse_multicomponent
  use mod_mc_config, only : mc_config, initialize_mc_config, read_mc_config, &
    print_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_provider_registry, only : validate_mc_providers, &
    print_mc_providers
  implicit none

  type(mc_config) :: config
  type(mc_state_layout) :: layout
  character(len=512) :: input_path
  logical :: input_exists

  input_path = 'input.dat'
  if (command_argument_count() >= 1) call get_command_argument(1, input_path)
  inquire(file=trim(input_path), exist=input_exists)
  if (input_exists) then
    call read_mc_config(trim(input_path), config)
  else if (command_argument_count() >= 1) then
    write(*,'(A,A)') 'ERROR: multicomponent input does not exist: ', &
      trim(input_path)
    error stop 'missing multicomponent input'
  else
    call initialize_mc_config(config)
  end if

  if (trim(config%simulation_mode) /= 'foundation') then
    error stop 'stage-0 executable requires simulation_mode=foundation'
  end if

  call initialize_mc_state_layout(layout, config%nspecies)
  call validate_mc_providers(config)
  call print_mc_config(config)
  call print_mc_providers()
  write(*,'(A,I0)') 'conservative variables = ', layout%nvariables
  write(*,'(A)') 'Multicomponent NSE stage-0 foundation initialized successfully.'
  write(*,'(A)') 'No flow advancement is performed by the stage-0 executable.'
end program main_nse_multicomponent
