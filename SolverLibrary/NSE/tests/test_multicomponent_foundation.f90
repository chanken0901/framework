program test_multicomponent_foundation
  use mod_mc_config, only : mc_config, initialize_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout, mc_species_index
  use mod_mc_provider_registry, only : validate_mc_providers
  use mod_mc_thermodynamics_provider, only : &
    mc_thermodynamics_provider_name
  use mod_mc_transport_provider, only : mc_transport_provider_name
  use mod_mc_chemistry_provider, only : mc_chemistry_provider_name
  implicit none

  type(mc_config) :: config
  type(mc_state_layout) :: layout

  call initialize_mc_config(config)
  call assert_true(trim(config%simulation_mode) == 'foundation', &
    'stage-0 default simulation mode must remain foundation')
  call initialize_mc_state_layout(layout, config%nspecies)
  config%thermodynamics_model = mc_thermodynamics_provider_name
  config%transport_model = mc_transport_provider_name
  config%chemistry_model = mc_chemistry_provider_name
  call assert_true(layout%nvariables == 5, &
    'one species must reduce to five conservative variables')
  call assert_true(mc_species_index(layout, 1) == 1, &
    'first species index must be one')
  call assert_true(all(layout%momentum == [2, 3, 4]), &
    'one-species momentum indices are invalid')
  call assert_true(layout%total_energy == 5, &
    'one-species energy index must be five')
  call validate_mc_providers(config)

  call initialize_mc_state_layout(layout, 3)
  call assert_true(layout%nvariables == 7, &
    'three species must have seven conservative variables')
  call assert_true(mc_species_index(layout, 3) == 3, &
    'third species index is invalid')
  call assert_true(all(layout%momentum == [4, 5, 6]), &
    'three-species momentum indices are invalid')
  call assert_true(layout%total_energy == 7, &
    'three-species energy index is invalid')

  write(*,'(A)') 'Multicomponent foundation contract tests passed.'

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'multicomponent foundation contract failure'
    end if
  end subroutine assert_true

end program test_multicomponent_foundation
