module mod_mc_reactor_config
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_max_species
  implicit none
  private

  integer, parameter :: reactor_name_length = 32
  integer, parameter :: reactor_path_length = 256

  type, public :: mc_reactor_config
    real(dp) :: initial_density = 1.0_dp
    real(dp) :: initial_temperature = 1000.0_dp
    real(dp) :: initial_mass_fractions(mc_max_species) = 0.0_dp
    real(dp) :: dt = 0.0_dp
    real(dp) :: maximum_dt = 1.0e-3_dp
    real(dp) :: chemistry_cfl = 0.1_dp
    integer :: nsteps = 100
    integer :: output_every = 1
    logical :: write_history = .true.
    character(len=reactor_path_length) :: output_file = &
      'homogeneous_reactor.csv'
  end type mc_reactor_config

  public :: read_mc_reactor_config
  public :: validate_mc_reactor_config

contains

  subroutine read_mc_reactor_config( &
      path,nspecies,species_names,config)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: species_names(:)
    type(mc_reactor_config), intent(out) :: config
    integer :: unit, ios, species
    character(len=512) :: message
    character(len=reactor_name_length) :: &
      reactor_species_names(mc_max_species)
    real(dp) :: initial_density, initial_temperature
    real(dp) :: initial_mass_fractions(mc_max_species)
    real(dp) :: dt, maximum_dt, chemistry_cfl
    integer :: nsteps, output_every
    logical :: write_history
    character(len=reactor_path_length) :: output_file
    namelist /homogeneous_reactor/ reactor_species_names, &
      initial_density, initial_temperature, initial_mass_fractions, &
      dt, maximum_dt, chemistry_cfl, nsteps, output_every, &
      write_history, output_file

    if (nspecies < 2 .or. nspecies > mc_max_species .or. &
        size(species_names) < nspecies) then
      error stop 'invalid homogeneous-reactor species contract'
    end if
    config = mc_reactor_config()
    reactor_species_names = ''
    initial_density = config%initial_density
    initial_temperature = config%initial_temperature
    initial_mass_fractions = config%initial_mass_fractions
    dt = config%dt
    maximum_dt = config%maximum_dt
    chemistry_cfl = config%chemistry_cfl
    nsteps = config%nsteps
    output_every = config%output_every
    write_history = config%write_history
    output_file = config%output_file

    open(newunit=unit,file=trim(path),status='old',action='read', &
      iostat=ios,iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open homogeneous-reactor input: ', &
        trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open homogeneous-reactor input'
    end if
    read(unit,nml=homogeneous_reactor,iostat=ios,iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid homogeneous_reactor namelist: ', &
        trim(message)
      error stop 'failed to read homogeneous-reactor input'
    end if

    do species = 1, nspecies
      if (trim(reactor_species_names(species)) /= &
          trim(species_names(species))) then
        error stop 'reactor species order does not match physics species'
      end if
    end do
    config%initial_density = initial_density
    config%initial_temperature = initial_temperature
    config%initial_mass_fractions = initial_mass_fractions
    config%dt = dt
    config%maximum_dt = maximum_dt
    config%chemistry_cfl = chemistry_cfl
    config%nsteps = nsteps
    config%output_every = output_every
    config%write_history = write_history
    config%output_file = adjustl(output_file)
    call validate_mc_reactor_config(config,nspecies)
  end subroutine read_mc_reactor_config

  subroutine validate_mc_reactor_config(config,nspecies)
    type(mc_reactor_config), intent(in) :: config
    integer, intent(in) :: nspecies

    if (.not. ieee_is_finite(config%initial_density) .or. &
        config%initial_density <= 0.0_dp) then
      error stop 'homogeneous-reactor density must be finite and positive'
    end if
    if (.not. ieee_is_finite(config%initial_temperature) .or. &
        config%initial_temperature <= 0.0_dp) then
      error stop 'homogeneous-reactor temperature must be finite and positive'
    end if
    if (.not. all(ieee_is_finite( &
        config%initial_mass_fractions(1:nspecies))) .or. &
        minval(config%initial_mass_fractions(1:nspecies)) < 0.0_dp .or. &
        abs(sum(config%initial_mass_fractions(1:nspecies))-1.0_dp) > &
        1.0e-12_dp) then
      error stop 'homogeneous-reactor mass fractions must be nonnegative and sum to one'
    end if
    if (.not. ieee_is_finite(config%dt) .or. config%dt < 0.0_dp .or. &
        .not. ieee_is_finite(config%maximum_dt) .or. &
        config%maximum_dt <= 0.0_dp) then
      error stop 'homogeneous-reactor timestep settings are invalid'
    end if
    if (config%dt > config%maximum_dt) then
      error stop 'homogeneous-reactor dt must not exceed maximum_dt'
    end if
    if (.not. ieee_is_finite(config%chemistry_cfl) .or. &
        config%chemistry_cfl <= 0.0_dp .or. &
        config%chemistry_cfl > 1.0_dp) then
      error stop 'homogeneous-reactor chemistry_cfl must be in (0,1]'
    end if
    if (config%nsteps < 0 .or. config%output_every < 1) then
      error stop 'homogeneous-reactor step counts are invalid'
    end if
    if (config%write_history .and. len_trim(config%output_file) == 0) then
      error stop 'homogeneous-reactor output file must not be empty'
    end if
  end subroutine validate_mc_reactor_config

end module mod_mc_reactor_config
