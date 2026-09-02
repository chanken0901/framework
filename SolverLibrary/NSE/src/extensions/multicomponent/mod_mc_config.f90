module mod_mc_config
  implicit none
  private

  integer, parameter, public :: mc_max_species = 64
  integer, parameter :: mc_name_length = 32

  type, public :: mc_config
    integer :: nspecies = 1
    character(len=mc_name_length) :: species_names(mc_max_species) = ''
    character(len=mc_name_length) :: simulation_mode = 'foundation'
    character(len=mc_name_length) :: thermodynamics_model = &
      'calorically_perfect'
    character(len=mc_name_length) :: transport_model = 'none'
    character(len=mc_name_length) :: chemistry_model = 'none'
  end type mc_config

  public :: initialize_mc_config
  public :: read_mc_config
  public :: validate_mc_config
  public :: print_mc_config

contains

  subroutine initialize_mc_config(config)
    type(mc_config), intent(out) :: config

    config = mc_config()
    config%species_names(1) = 'mixture'
  end subroutine initialize_mc_config

  subroutine read_mc_config(path, config)
    character(len=*), intent(in) :: path
    type(mc_config), intent(out) :: config
    integer :: unit, ios
    integer :: nspecies
    character(len=mc_name_length) :: species_names(mc_max_species)
    character(len=mc_name_length) :: simulation_mode
    character(len=mc_name_length) :: thermodynamics_model
    character(len=mc_name_length) :: transport_model
    character(len=mc_name_length) :: chemistry_model
    character(len=512) :: message
    namelist /multicomponent/ nspecies, species_names, simulation_mode, &
      thermodynamics_model, transport_model, chemistry_model

    call initialize_mc_config(config)
    nspecies = config%nspecies
    species_names = config%species_names
    simulation_mode = config%simulation_mode
    thermodynamics_model = config%thermodynamics_model
    transport_model = config%transport_model
    chemistry_model = config%chemistry_model

    open(newunit=unit, file=trim(path), status='old', action='read', &
      iostat=ios, iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open multicomponent input: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open multicomponent input'
    end if
    read(unit, nml=multicomponent, iostat=ios, iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid multicomponent namelist: ', &
        trim(message)
      error stop 'failed to read multicomponent input'
    end if

    config%nspecies = nspecies
    config%species_names = species_names
    config%simulation_mode = adjustl(simulation_mode)
    config%thermodynamics_model = adjustl(thermodynamics_model)
    config%transport_model = adjustl(transport_model)
    config%chemistry_model = adjustl(chemistry_model)
    call validate_mc_config(config)
  end subroutine read_mc_config

  subroutine validate_mc_config(config)
    type(mc_config), intent(in) :: config
    integer :: species, other

    if (config%nspecies < 1 .or. config%nspecies > mc_max_species) then
      error stop 'multicomponent nspecies is outside the supported range'
    end if
    do species = 1, config%nspecies
      if (len_trim(config%species_names(species)) == 0) then
        error stop 'every multicomponent species must have a name'
      end if
      do other = species + 1, config%nspecies
        if (trim(config%species_names(species)) == &
            trim(config%species_names(other))) then
          error stop 'multicomponent species names must be unique'
        end if
      end do
    end do
    if (len_trim(config%thermodynamics_model) == 0) then
      error stop 'thermodynamics model must not be empty'
    end if
    if (trim(config%simulation_mode) /= 'foundation' .and. &
        trim(config%simulation_mode) /= 'passive_scalar' .and. &
        trim(config%simulation_mode) /= 'inviscid_euler' .and. &
        trim(config%simulation_mode) /= 'thermally_perfect_euler' .and. &
        trim(config%simulation_mode) /= 'viscous_navier_stokes') then
      error stop 'unsupported multicomponent simulation mode'
    end if
    if (len_trim(config%transport_model) == 0) then
      error stop 'transport model must not be empty'
    end if
    if (len_trim(config%chemistry_model) == 0) then
      error stop 'chemistry model must not be empty'
    end if
  end subroutine validate_mc_config

  subroutine print_mc_config(config, unit)
    type(mc_config), intent(in) :: config
    integer, intent(in), optional :: unit
    integer :: output_unit, species

    output_unit = 6
    if (present(unit)) output_unit = unit
    write(output_unit,'(A)') '--- multicomponent configuration ---'
    write(output_unit,'(A,I0)') 'nspecies = ', config%nspecies
    write(output_unit,'(A,A)') 'simulation mode = ', &
      trim(config%simulation_mode)
    do species = 1, config%nspecies
      write(output_unit,'(A,I0,A,A)') 'species(', species, ') = ', &
        trim(config%species_names(species))
    end do
    write(output_unit,'(A,A)') 'thermodynamics = ', &
      trim(config%thermodynamics_model)
    write(output_unit,'(A,A)') 'transport = ', trim(config%transport_model)
    write(output_unit,'(A,A)') 'chemistry = ', trim(config%chemistry_model)
  end subroutine print_mc_config

end module mod_mc_config
