module mod_mc_provider_registry
  use mod_mc_config, only : mc_config
  use mod_mc_thermodynamics_provider, only : &
    mc_thermodynamics_provider_name, &
    mc_thermodynamics_supports_reactions, &
    validate_mc_thermodynamics_provider
  use mod_mc_transport_provider, only : mc_transport_provider_name, &
    mc_transport_has_species_diffusion, validate_mc_transport_provider
  use mod_mc_chemistry_provider, only : mc_chemistry_provider_name, &
    mc_chemistry_is_reactive, validate_mc_chemistry_provider
  implicit none
  private

  public :: validate_mc_providers
  public :: print_mc_providers

contains

  subroutine validate_mc_providers(config)
    type(mc_config), intent(in) :: config

    call validate_mc_thermodynamics_provider(config%thermodynamics_model)
    call validate_mc_transport_provider(config%transport_model)
    call validate_mc_chemistry_provider(config%chemistry_model)
    if (mc_chemistry_is_reactive .and. &
        .not. mc_thermodynamics_supports_reactions) then
      error stop 'reactive chemistry requires a reaction-capable thermodynamics provider'
    end if
  end subroutine validate_mc_providers

  subroutine print_mc_providers(unit)
    integer, intent(in), optional :: unit
    integer :: output_unit

    output_unit = 6
    if (present(unit)) output_unit = unit
    write(output_unit,'(A,A)') 'compiled thermodynamics provider = ', &
      mc_thermodynamics_provider_name
    write(output_unit,'(A,A)') 'compiled transport provider = ', &
      mc_transport_provider_name
    write(output_unit,'(A,L1)') 'species diffusion available = ', &
      mc_transport_has_species_diffusion
    write(output_unit,'(A,A)') 'compiled chemistry provider = ', &
      mc_chemistry_provider_name
    write(output_unit,'(A,L1)') 'reactive chemistry available = ', &
      mc_chemistry_is_reactive
  end subroutine print_mc_providers

end module mod_mc_provider_registry
