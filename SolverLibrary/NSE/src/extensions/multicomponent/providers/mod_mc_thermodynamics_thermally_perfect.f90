module mod_mc_thermodynamics_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_max_species
  use mod_mc_state_layout, only : mc_state_layout
  implicit none
  private

  integer, parameter :: nasa_coefficient_count = 7
  integer, parameter :: thermo_name_length = 32
  real(dp), parameter :: default_universal_gas_constant = &
    8.31446261815324e3_dp

  character(len=*), parameter, public :: mc_thermodynamics_provider_name = &
    'thermally_perfect'
  logical, parameter, public :: mc_thermodynamics_supports_reactions = .false.

  integer, save :: configured_species = 0
  real(dp), save :: universal_gas_constant = default_universal_gas_constant
  real(dp), save :: temperature_min = 200.0_dp
  real(dp), save :: temperature_max = 6000.0_dp
  real(dp), save :: temperature_tolerance = 1.0e-10_dp
  real(dp), save :: molecular_weights(mc_max_species) = 0.0_dp
  real(dp), save :: temperature_midpoints(mc_max_species) = 0.0_dp
  real(dp), save :: nasa_low_coefficients( &
    nasa_coefficient_count,mc_max_species) = 0.0_dp
  real(dp), save :: nasa_high_coefficients( &
    nasa_coefficient_count,mc_max_species) = 0.0_dp

  public :: validate_mc_thermodynamics_provider
  public :: configure_mc_thermodynamics
  public :: mc_mixture_density
  public :: mc_mixture_gas_constant
  public :: mc_mixture_cp
  public :: mc_mixture_gamma
  public :: mc_species_enthalpies
  public :: mc_pressure
  public :: mc_temperature
  public :: mc_sound_speed
  public :: mc_total_energy_from_primitive

contains

  subroutine validate_mc_thermodynamics_provider(requested_model)
    character(len=*), intent(in) :: requested_model

    if (trim(adjustl(requested_model)) /= mc_thermodynamics_provider_name) then
      write(*,'(A,A,A,A)') 'ERROR: requested thermodynamics provider "', &
        trim(requested_model), '" but this executable contains "', &
        mc_thermodynamics_provider_name
      error stop 'multicomponent thermodynamics provider mismatch'
    end if
  end subroutine validate_mc_thermodynamics_provider

  subroutine configure_mc_thermodynamics(path,nspecies,species_names)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: species_names(:)
    integer :: unit, ios, species
    character(len=512) :: message
    character(len=thermo_name_length) :: thermo_species_names(mc_max_species)
    namelist /thermally_perfect/ thermo_species_names, &
      universal_gas_constant, temperature_min, temperature_max, &
      temperature_tolerance, molecular_weights, temperature_midpoints, &
      nasa_low_coefficients, nasa_high_coefficients

    if (nspecies < 1 .or. nspecies > mc_max_species .or. &
        size(species_names) < nspecies) then
      error stop 'invalid species contract for thermally-perfect provider'
    end if

    configured_species = 0
    thermo_species_names = ''
    universal_gas_constant = default_universal_gas_constant
    temperature_min = 200.0_dp
    temperature_max = 6000.0_dp
    temperature_tolerance = 1.0e-10_dp
    molecular_weights = 0.0_dp
    temperature_midpoints = 0.0_dp
    nasa_low_coefficients = 0.0_dp
    nasa_high_coefficients = 0.0_dp

    open(newunit=unit,file=trim(path),status='old',action='read', &
      iostat=ios,iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open thermodynamics input: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open thermally-perfect input'
    end if
    read(unit,nml=thermally_perfect,iostat=ios,iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid thermally_perfect namelist: ', &
        trim(message)
      error stop 'failed to read thermally-perfect input'
    end if

    do species = 1, nspecies
      if (trim(thermo_species_names(species)) /= &
          trim(species_names(species))) then
        write(*,'(A,I0,A,A,A,A)') 'ERROR: thermodynamics species(', species, &
          ')="', trim(thermo_species_names(species)), &
          '" does not match physics species "', trim(species_names(species))
        error stop 'thermodynamics species order mismatch'
      end if
    end do
    configured_species = nspecies
    call validate_thermally_perfect_data()
  end subroutine configure_mc_thermodynamics

  subroutine validate_thermally_perfect_data()
    integer :: species, sample
    real(dp) :: temperature, cp_value, cv_value, gas_constant
    real(dp) :: cp_low, cp_high, energy_low, energy_high, scale

    if (.not. ieee_is_finite(universal_gas_constant) .or. &
        universal_gas_constant <= 0.0_dp) then
      error stop 'universal gas constant must be finite and positive'
    end if
    if (.not. ieee_is_finite(temperature_min) .or. &
        .not. ieee_is_finite(temperature_max) .or. &
        temperature_min <= 0.0_dp .or. temperature_max <= temperature_min) then
      error stop 'thermally-perfect temperature range is invalid'
    end if
    if (.not. ieee_is_finite(temperature_tolerance) .or. &
        temperature_tolerance <= 0.0_dp .or. &
        temperature_tolerance >= 1.0e-3_dp) then
      error stop 'thermally-perfect temperature tolerance is invalid'
    end if

    do species = 1, configured_species
      if (.not. ieee_is_finite(molecular_weights(species)) .or. &
          molecular_weights(species) <= 0.0_dp) then
        error stop 'every species molecular weight must be finite and positive'
      end if
      if (.not. ieee_is_finite(temperature_midpoints(species)) .or. &
          temperature_midpoints(species) <= temperature_min .or. &
          temperature_midpoints(species) >= temperature_max) then
        error stop 'every NASA midpoint must lie inside the temperature range'
      end if
      if (.not. all(ieee_is_finite( &
          nasa_low_coefficients(:,species))) .or. &
          .not. all(ieee_is_finite( &
          nasa_high_coefficients(:,species)))) then
        error stop 'NASA-7 coefficients must be finite'
      end if

      gas_constant = species_gas_constant(species)
      do sample = 0, 8
        temperature = temperature_min + &
          (temperature_midpoints(species)-temperature_min)* &
          real(sample,dp)/8.0_dp
        cp_value = species_cp_with_coefficients( &
          species,temperature,nasa_low_coefficients(:,species))
        cv_value = cp_value-gas_constant
        if (.not. ieee_is_finite(cp_value) .or. cv_value <= 0.0_dp) then
          error stop 'low-temperature NASA-7 coefficients produce invalid cv'
        end if
        temperature = temperature_midpoints(species) + &
          (temperature_max-temperature_midpoints(species))* &
          real(sample,dp)/8.0_dp
        cp_value = species_cp_with_coefficients( &
          species,temperature,nasa_high_coefficients(:,species))
        cv_value = cp_value-gas_constant
        if (.not. ieee_is_finite(cp_value) .or. cv_value <= 0.0_dp) then
          error stop 'high-temperature NASA-7 coefficients produce invalid cv'
        end if
      end do

      cp_low = species_cp_with_coefficients( &
        species,temperature_midpoints(species), &
        nasa_low_coefficients(:,species))
      cp_high = species_cp_with_coefficients( &
        species,temperature_midpoints(species), &
        nasa_high_coefficients(:,species))
      scale = max(abs(cp_low),abs(cp_high),1.0_dp)
      if (abs(cp_high-cp_low) > 1.0e-3_dp*scale) then
        error stop 'NASA-7 cp is discontinuous at a species midpoint'
      end if
      energy_low = species_internal_energy_with_coefficients( &
        species,temperature_midpoints(species), &
        nasa_low_coefficients(:,species))
      energy_high = species_internal_energy_with_coefficients( &
        species,temperature_midpoints(species), &
        nasa_high_coefficients(:,species))
      scale = max(abs(energy_low),abs(energy_high),1.0_dp)
      if (abs(energy_high-energy_low) > 1.0e-3_dp*scale) then
        error stop 'NASA-7 internal energy is discontinuous at a midpoint'
      end if
    end do
  end subroutine validate_thermally_perfect_data

  subroutine ensure_configured(layout)
    type(mc_state_layout), intent(in) :: layout

    if (configured_species == 0) then
      error stop 'thermally-perfect provider has not been configured'
    end if
    if (layout%nspecies /= configured_species) then
      error stop 'state layout and thermodynamics species counts differ'
    end if
  end subroutine ensure_configured

  real(dp) function species_gas_constant(species) result(gas_constant)
    integer, intent(in) :: species

    gas_constant = universal_gas_constant/molecular_weights(species)
  end function species_gas_constant

  real(dp) function species_cp_with_coefficients( &
      species,temperature,coefficients) result(cp_value)
    integer, intent(in) :: species
    real(dp), intent(in) :: temperature
    real(dp), intent(in) :: coefficients(nasa_coefficient_count)
    real(dp) :: cp_over_gas_constant

    cp_over_gas_constant = coefficients(1) + temperature*( &
      coefficients(2) + temperature*(coefficients(3) + temperature*( &
      coefficients(4) + temperature*coefficients(5))))
    cp_value = species_gas_constant(species)*cp_over_gas_constant
  end function species_cp_with_coefficients

  real(dp) function species_cp(species,temperature) result(cp_value)
    integer, intent(in) :: species
    real(dp), intent(in) :: temperature

    if (temperature <= temperature_midpoints(species)) then
      cp_value = species_cp_with_coefficients( &
        species,temperature,nasa_low_coefficients(:,species))
    else
      cp_value = species_cp_with_coefficients( &
        species,temperature,nasa_high_coefficients(:,species))
    end if
  end function species_cp

  real(dp) function species_internal_energy_with_coefficients( &
      species,temperature,coefficients) result(internal_energy)
    integer, intent(in) :: species
    real(dp), intent(in) :: temperature
    real(dp), intent(in) :: coefficients(nasa_coefficient_count)
    real(dp) :: temperature2, temperature3, temperature4, temperature5

    temperature2 = temperature*temperature
    temperature3 = temperature2*temperature
    temperature4 = temperature3*temperature
    temperature5 = temperature4*temperature
    internal_energy = species_gas_constant(species)*( &
      (coefficients(1)-1.0_dp)*temperature + &
      0.5_dp*coefficients(2)*temperature2 + &
      coefficients(3)*temperature3/3.0_dp + &
      0.25_dp*coefficients(4)*temperature4 + &
      0.2_dp*coefficients(5)*temperature5 + coefficients(6))
  end function species_internal_energy_with_coefficients

  real(dp) function species_internal_energy( &
      species,temperature) result(internal_energy)
    integer, intent(in) :: species
    real(dp), intent(in) :: temperature

    if (temperature <= temperature_midpoints(species)) then
      internal_energy = species_internal_energy_with_coefficients( &
        species,temperature,nasa_low_coefficients(:,species))
    else
      internal_energy = species_internal_energy_with_coefficients( &
        species,temperature,nasa_high_coefficients(:,species))
    end if
  end function species_internal_energy

  real(dp) function species_enthalpy( &
      species,temperature) result(enthalpy)
    integer, intent(in) :: species
    real(dp), intent(in) :: temperature

    enthalpy = species_internal_energy(species,temperature) + &
      species_gas_constant(species)*temperature
  end function species_enthalpy

  real(dp) function mixture_internal_energy( &
      mass_fractions,temperature) result(internal_energy)
    real(dp), intent(in) :: mass_fractions(:), temperature
    integer :: species

    internal_energy = 0.0_dp
    do species = 1, configured_species
      internal_energy = internal_energy + mass_fractions(species)* &
        species_internal_energy(species,temperature)
    end do
  end function mixture_internal_energy

  real(dp) function recover_temperature( &
      mass_fractions,target_energy) result(temperature)
    real(dp), intent(in) :: mass_fractions(:), target_energy
    integer :: iteration
    real(dp) :: lower_temperature, upper_temperature, midpoint
    real(dp) :: lower_energy, upper_energy, midpoint_energy
    real(dp) :: bound_scale, convergence_scale

    lower_temperature = temperature_min
    upper_temperature = temperature_max
    lower_energy = mixture_internal_energy( &
      mass_fractions,lower_temperature)
    upper_energy = mixture_internal_energy( &
      mass_fractions,upper_temperature)
    bound_scale = max( &
      abs(target_energy),abs(lower_energy),abs(upper_energy),1.0_dp)
    convergence_scale = max(abs(target_energy),1.0_dp)
    if (target_energy < lower_energy-temperature_tolerance*bound_scale .or. &
        target_energy > upper_energy+temperature_tolerance*bound_scale) then
      error stop 'mixture internal energy is outside the NASA temperature range'
    end if
    if (abs(target_energy-lower_energy) <= &
        temperature_tolerance*convergence_scale) then
      temperature = lower_temperature
      return
    end if
    if (abs(target_energy-upper_energy) <= &
        temperature_tolerance*convergence_scale) then
      temperature = upper_temperature
      return
    end if

    do iteration = 1, 100
      midpoint = 0.5_dp*(lower_temperature+upper_temperature)
      midpoint_energy = mixture_internal_energy(mass_fractions,midpoint)
      if (midpoint_energy < target_energy) then
        lower_temperature = midpoint
      else
        upper_temperature = midpoint
      end if
      if (abs(midpoint_energy-target_energy) <= &
          temperature_tolerance*convergence_scale .or. &
          upper_temperature-lower_temperature <= &
          temperature_tolerance*max(midpoint,1.0_dp)) exit
    end do
    temperature = 0.5_dp*(lower_temperature+upper_temperature)
  end function recover_temperature

  real(dp) function mc_mixture_density(state,layout) result(density)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout

    density = sum(state(layout%first_species:layout%last_species))
  end function mc_mixture_density

  real(dp) function mc_mixture_gas_constant( &
      mass_fractions,layout) result(gas_constant)
    real(dp), intent(in) :: mass_fractions(:)
    type(mc_state_layout), intent(in) :: layout
    integer :: species

    call ensure_configured(layout)
    if (size(mass_fractions) < layout%nspecies) then
      error stop 'mixture mass-fraction vector is too short'
    end if
    gas_constant = 0.0_dp
    do species = 1, layout%nspecies
      gas_constant = gas_constant + &
        mass_fractions(species)*species_gas_constant(species)
    end do
  end function mc_mixture_gas_constant

  real(dp) function mc_mixture_cp( &
      mass_fractions,layout,temperature) result(cp_value)
    real(dp), intent(in) :: mass_fractions(:), temperature
    type(mc_state_layout), intent(in) :: layout
    integer :: species

    call ensure_configured(layout)
    if (size(mass_fractions) < layout%nspecies) then
      error stop 'mixture mass-fraction vector is too short'
    end if
    if (temperature < temperature_min .or. temperature > temperature_max) then
      error stop 'mixture cp temperature is outside the NASA range'
    end if
    cp_value = 0.0_dp
    do species = 1, layout%nspecies
      cp_value = cp_value + &
        mass_fractions(species)*species_cp(species,temperature)
    end do
  end function mc_mixture_cp

  real(dp) function mc_mixture_gamma( &
      mass_fractions,layout,temperature,gamma) result(gamma_value)
    real(dp), intent(in) :: mass_fractions(:), temperature, gamma
    type(mc_state_layout), intent(in) :: layout
    real(dp) :: cp_value, gas_constant

    if (gamma <= 0.0_dp) error stop 'fallback gamma must be positive'
    cp_value = mc_mixture_cp(mass_fractions,layout,temperature)
    gas_constant = mc_mixture_gas_constant(mass_fractions,layout)
    gamma_value = cp_value/(cp_value-gas_constant)
  end function mc_mixture_gamma

  subroutine mc_species_enthalpies( &
      layout,temperature,enthalpies)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: enthalpies(:)
    integer :: species

    call ensure_configured(layout)
    if (temperature < temperature_min .or. temperature > temperature_max .or. &
        size(enthalpies) < layout%nspecies) then
      error stop 'species-enthalpy request is outside the provider contract'
    end if
    do species = 1, layout%nspecies
      enthalpies(species) = species_enthalpy(species,temperature)
    end do
  end subroutine mc_species_enthalpies

  real(dp) function mc_temperature(state,layout,gamma) result(temperature)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    real(dp) :: density, kinetic_energy, target_energy
    real(dp) :: mass_fractions(layout%nspecies)

    if (gamma <= 0.0_dp) error stop 'fallback gamma must be positive'
    call ensure_configured(layout)
    density = mc_mixture_density(state,layout)
    if (density <= 0.0_dp) error stop 'mixture density must be positive'
    mass_fractions = &
      state(layout%first_species:layout%last_species)/density
    kinetic_energy = 0.5_dp*sum(state(layout%momentum)**2)/density
    target_energy = (state(layout%total_energy)-kinetic_energy)/density
    temperature = recover_temperature(mass_fractions,target_energy)
  end function mc_temperature

  real(dp) function mc_pressure(state,layout,gamma) result(pressure)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    real(dp) :: density, temperature
    real(dp) :: mass_fractions(layout%nspecies)

    density = mc_mixture_density(state,layout)
    temperature = mc_temperature(state,layout,gamma)
    mass_fractions = &
      state(layout%first_species:layout%last_species)/density
    pressure = density*mc_mixture_gas_constant( &
      mass_fractions,layout)*temperature
  end function mc_pressure

  real(dp) function mc_sound_speed(state,layout,gamma) result(speed)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    real(dp) :: density, pressure, temperature, gamma_value
    real(dp) :: mass_fractions(layout%nspecies)

    density = mc_mixture_density(state,layout)
    temperature = mc_temperature(state,layout,gamma)
    mass_fractions = &
      state(layout%first_species:layout%last_species)/density
    pressure = density*mc_mixture_gas_constant( &
      mass_fractions,layout)*temperature
    gamma_value = mc_mixture_gamma( &
      mass_fractions,layout,temperature,gamma)
    speed = sqrt(gamma_value*pressure/density)
  end function mc_sound_speed

  real(dp) function mc_total_energy_from_primitive( &
      layout,gamma,density,velocity,pressure,mass_fractions) &
      result(total_energy_density)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma, density, velocity(3), pressure
    real(dp), intent(in) :: mass_fractions(:)
    real(dp) :: gas_constant, temperature, internal_energy

    if (gamma <= 0.0_dp) error stop 'fallback gamma must be positive'
    call ensure_configured(layout)
    if (density <= 0.0_dp .or. pressure <= 0.0_dp) then
      error stop 'primitive density and pressure must be positive'
    end if
    if (size(mass_fractions) < layout%nspecies .or. &
        minval(mass_fractions(1:layout%nspecies)) < 0.0_dp .or. &
        abs(sum(mass_fractions(1:layout%nspecies))-1.0_dp) > 1.0e-12_dp) then
      error stop 'primitive mass fractions are invalid'
    end if
    gas_constant = mc_mixture_gas_constant( &
      mass_fractions(1:layout%nspecies),layout)
    temperature = pressure/(density*gas_constant)
    if (temperature < temperature_min .or. temperature > temperature_max) then
      error stop 'primitive state is outside the NASA temperature range'
    end if
    internal_energy = mixture_internal_energy( &
      mass_fractions(1:layout%nspecies),temperature)
    total_energy_density = density*( &
      internal_energy+0.5_dp*sum(velocity**2))
  end function mc_total_energy_from_primitive

end module mod_mc_thermodynamics_provider
