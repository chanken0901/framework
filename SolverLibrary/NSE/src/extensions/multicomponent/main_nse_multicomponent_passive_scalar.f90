program main_nse_multicomponent_passive_scalar
  use mod_mc_config, only : mc_config, read_mc_config, print_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_provider_registry, only : validate_mc_providers, &
    print_mc_providers
  use mod_mc_passive_scalar_config, only : mc_passive_scalar_config, &
    read_mc_passive_scalar_config
  use mod_mc_passive_scalar_solver, only : run_mc_passive_scalar
  implicit none

  type(mc_config) :: config
  type(mc_state_layout) :: layout
  type(mc_passive_scalar_config) :: passive
  character(len=512) :: input_path
  logical :: input_exists

  input_path = 'input.dat'
  if (command_argument_count() >= 1) call get_command_argument(1, input_path)
  inquire(file=trim(input_path), exist=input_exists)
  if (.not. input_exists) then
    write(*,'(A,A)') 'ERROR: passive-scalar input does not exist: ', &
      trim(input_path)
    error stop 'missing passive-scalar input'
  end if

  call read_mc_config(trim(input_path), config)
  if (trim(config%simulation_mode) /= 'passive_scalar') then
    write(*,'(A,A)') 'ERROR: stage-1 executable requires simulation_mode=', &
      'passive_scalar'
    error stop 'incorrect multicomponent simulation mode'
  end if
  if (config%nspecies < 2) then
    error stop 'passive-scalar mode requires at least two species'
  end if
  call read_mc_passive_scalar_config(trim(input_path), passive)
  call initialize_mc_state_layout(layout, config%nspecies)
  call validate_mc_providers(config)
  call print_mc_config(config)
  call print_mc_providers()
  write(*,'(A,I0)') 'conservative variables = ', layout%nvariables
  call run_mc_passive_scalar(config, layout, passive)
end program main_nse_multicomponent_passive_scalar
