module mod_mc_chemistry_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_max_species
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_thermodynamics_provider, only : mc_temperature, &
    mc_get_species_molecular_weights
  implicit none
  private

  integer, parameter :: chemistry_name_length = 32
  character(len=*), parameter, public :: mc_chemistry_provider_name = &
    'one_step_arrhenius'
  logical, parameter, public :: mc_chemistry_is_reactive = .true.

  integer, save :: configured_species = 0
  real(dp), save :: reactant_stoich(mc_max_species) = 0.0_dp
  real(dp), save :: product_stoich(mc_max_species) = 0.0_dp
  real(dp), save :: reaction_orders(mc_max_species) = 0.0_dp
  real(dp), save :: molecular_weights(mc_max_species) = 0.0_dp
  real(dp), save :: pre_exponential_factor = 0.0_dp
  real(dp), save :: temperature_exponent = 0.0_dp
  real(dp), save :: activation_temperature = 0.0_dp

  public :: validate_mc_chemistry_provider
  public :: configure_mc_chemistry
  public :: compute_mc_chemistry_source
  public :: compute_mc_chemistry_timestep
  public :: mc_export_chemistry

contains

  subroutine mc_export_chemistry(n,controls,data)
    integer,intent(in)::n
    real(dp),intent(out)::controls(3),data(3,n)
    if(n/=configured_species .or. n<1) error stop 'chemistry export before configuration'
    controls=[pre_exponential_factor,temperature_exponent,activation_temperature]
    data(1,:)=reactant_stoich(1:n)
    data(2,:)=product_stoich(1:n)
    data(3,:)=reaction_orders(1:n)
  end subroutine

  subroutine validate_mc_chemistry_provider(requested_model)
    character(len=*), intent(in) :: requested_model

    if (trim(adjustl(requested_model)) /= mc_chemistry_provider_name) then
      write(*,'(A,A,A,A)') 'ERROR: requested chemistry provider "', &
        trim(requested_model), '" but this executable contains "', &
        mc_chemistry_provider_name
      error stop 'multicomponent chemistry provider mismatch'
    end if
  end subroutine validate_mc_chemistry_provider

  subroutine configure_mc_chemistry(path,nspecies,species_names)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: species_names(:)
    integer :: unit, ios, species
    real(dp) :: mass_scale, mass_residual
    character(len=512) :: message
    character(len=chemistry_name_length) :: &
      chemistry_species_names(mc_max_species)
    namelist /one_step_arrhenius/ chemistry_species_names, &
      reactant_stoich, product_stoich, reaction_orders, &
      pre_exponential_factor, temperature_exponent, activation_temperature

    if (nspecies < 2 .or. nspecies > mc_max_species .or. &
        size(species_names) < nspecies) then
      error stop 'invalid species contract for one-step chemistry provider'
    end if

    configured_species = 0
    chemistry_species_names = ''
    reactant_stoich = 0.0_dp
    product_stoich = 0.0_dp
    reaction_orders = -1.0_dp
    molecular_weights = 0.0_dp
    pre_exponential_factor = 0.0_dp
    temperature_exponent = 0.0_dp
    activation_temperature = 0.0_dp

    open(newunit=unit,file=trim(path),status='old',action='read', &
      iostat=ios,iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open chemistry input: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open one-step chemistry input'
    end if
    read(unit,nml=one_step_arrhenius,iostat=ios,iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid one_step_arrhenius namelist: ', &
        trim(message)
      error stop 'failed to read one-step chemistry input'
    end if

    do species = 1, nspecies
      if (trim(chemistry_species_names(species)) /= &
          trim(species_names(species))) then
        error stop 'chemistry species order does not match physics species'
      end if
      if (reaction_orders(species) < 0.0_dp) then
        reaction_orders(species) = reactant_stoich(species)
      end if
    end do
    configured_species = nspecies
    call mc_get_species_molecular_weights( &
      configured_species,molecular_weights(1:configured_species))

    if (.not. ieee_is_finite(pre_exponential_factor) .or. &
        pre_exponential_factor <= 0.0_dp) then
      error stop 'chemistry pre-exponential factor must be finite and positive'
    end if
    if (.not. ieee_is_finite(temperature_exponent)) then
      error stop 'chemistry temperature exponent must be finite'
    end if
    if (.not. ieee_is_finite(activation_temperature) .or. &
        activation_temperature < 0.0_dp) then
      error stop 'chemistry activation temperature must be finite and nonnegative'
    end if
    if (.not. all(ieee_is_finite( &
        reactant_stoich(1:configured_species))) .or. &
        .not. all(ieee_is_finite( &
        product_stoich(1:configured_species))) .or. &
        .not. all(ieee_is_finite( &
        reaction_orders(1:configured_species)))) then
      error stop 'chemistry stoichiometry and orders must be finite'
    end if
    if (minval(reactant_stoich(1:configured_species)) < 0.0_dp .or. &
        minval(product_stoich(1:configured_species)) < 0.0_dp .or. &
        minval(reaction_orders(1:configured_species)) < 0.0_dp) then
      error stop 'chemistry stoichiometry and orders must be nonnegative'
    end if
    if (maxval(reactant_stoich(1:configured_species)) <= 0.0_dp .or. &
        maxval(product_stoich(1:configured_species)) <= 0.0_dp) then
      error stop 'one-step chemistry requires reactants and products'
    end if
    if (any(reactant_stoich(1:configured_species) > 0.0_dp .and. &
        product_stoich(1:configured_species) > 0.0_dp)) then
      error stop 'one-step chemistry does not support species on both sides'
    end if
    mass_residual = sum(molecular_weights(1:configured_species)*( &
      product_stoich(1:configured_species)- &
      reactant_stoich(1:configured_species)))
    mass_scale = max(sum(molecular_weights(1:configured_species)*( &
      product_stoich(1:configured_species)+ &
      reactant_stoich(1:configured_species))),1.0_dp)
    if (abs(mass_residual) > 1.0e-12_dp*mass_scale) then
      error stop 'one-step reaction stoichiometry does not conserve mass'
    end if
  end subroutine configure_mc_chemistry

  subroutine ensure_configured(layout)
    type(mc_state_layout), intent(in) :: layout

    if (configured_species == 0 .or. &
        layout%nspecies /= configured_species) then
      error stop 'one-step chemistry provider has not been configured'
    end if
  end subroutine ensure_configured

  subroutine compute_mc_chemistry_source( &
      state,layout,gamma,source,progress_rate)
    real(dp), intent(in) :: state(:), gamma
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(out) :: source(:)
    real(dp), intent(out), optional :: progress_rate
    integer :: species, state_index, correction_index
    real(dp) :: temperature, concentration, log_rate, rate, residual
    real(dp) :: log_tiny, log_huge

    call ensure_configured(layout)
    if (size(state) < layout%nvariables .or. &
        size(source) < layout%nvariables) then
      error stop 'chemistry source vectors do not match the state layout'
    end if
    source = 0.0_dp
    temperature = mc_temperature(state,layout,gamma)
    log_rate = log(pre_exponential_factor) + &
      temperature_exponent*log(temperature) - &
      activation_temperature/temperature
    do species = 1, configured_species
      if (reaction_orders(species) <= 0.0_dp) cycle
      state_index = layout%first_species+species-1
      concentration = state(state_index)/molecular_weights(species)
      if (concentration <= 0.0_dp) then
        if (present(progress_rate)) progress_rate = 0.0_dp
        return
      end if
      log_rate = log_rate + reaction_orders(species)*log(concentration)
    end do

    log_tiny = log(tiny(1.0_dp))
    log_huge = log(huge(1.0_dp))
    if (log_rate <= log_tiny) then
      rate = 0.0_dp
    else if (log_rate >= log_huge) then
      error stop 'one-step Arrhenius reaction rate overflow'
    else
      rate = exp(log_rate)
    end if
    if (.not. ieee_is_finite(rate) .or. rate < 0.0_dp) then
      error stop 'one-step Arrhenius reaction rate is invalid'
    end if

    correction_index = 0
    do species = 1, configured_species
      state_index = layout%first_species+species-1
      source(state_index) = molecular_weights(species)*( &
        product_stoich(species)-reactant_stoich(species))*rate
      if (product_stoich(species) > 0.0_dp) correction_index = state_index
    end do
    residual = sum(source(layout%first_species:layout%last_species))
    source(correction_index) = source(correction_index)-residual
    if (present(progress_rate)) progress_rate = rate
  end subroutine compute_mc_chemistry_source

  real(dp) function compute_mc_chemistry_timestep( &
      state,layout,gamma,safety,maximum_dt) result(dt)
    real(dp), intent(in) :: state(:), gamma, safety, maximum_dt
    type(mc_state_layout), intent(in) :: layout
    real(dp) :: source(layout%nvariables), depletion_time
    integer :: species, state_index
    logical :: has_consuming_species

    if (.not. ieee_is_finite(safety) .or. safety <= 0.0_dp .or. &
        safety > 1.0_dp .or. .not. ieee_is_finite(maximum_dt) .or. &
        maximum_dt <= 0.0_dp) then
      error stop 'chemistry timestep controls are invalid'
    end if
    call compute_mc_chemistry_source(state,layout,gamma,source)
    depletion_time = huge(1.0_dp)
    has_consuming_species = .false.
    do species = 1, configured_species
      state_index = layout%first_species+species-1
      if (source(state_index) < 0.0_dp) then
        has_consuming_species = .true.
        depletion_time = min( &
          depletion_time,state(state_index)/(-source(state_index)))
      end if
    end do
    if (.not. has_consuming_species) then
      dt = maximum_dt
    else
      dt = min(maximum_dt,safety*depletion_time)
    end if
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
      error stop 'computed chemistry timestep is invalid'
    end if
  end function compute_mc_chemistry_timestep

end module mod_mc_chemistry_provider
