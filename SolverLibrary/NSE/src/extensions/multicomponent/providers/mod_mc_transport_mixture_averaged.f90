module mod_mc_transport_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_max_species
  use mod_mc_state_layout, only : mc_state_layout
  implicit none
  private

  integer, parameter :: transport_name_length = 32

  character(len=*), parameter, public :: mc_transport_provider_name = &
    'mixture_averaged'
  logical, parameter, public :: mc_transport_has_species_diffusion = .true.

  integer, save :: configured_species = 0
  real(dp), save :: reference_dynamic_viscosity = 1.8e-5_dp
  real(dp), save :: prandtl_number = 0.72_dp
  real(dp), save :: species_diffusivities(mc_max_species) = 0.0_dp

  public :: validate_mc_transport_provider
  public :: configure_mc_transport
  public :: mc_dynamic_viscosity
  public :: mc_thermal_conductivity
  public :: mc_species_diffusivities
  public :: mc_export_transport

contains

  subroutine mc_export_transport(n,controls,data)
    integer,intent(in)::n
    real(dp),intent(out)::controls(2),data(n)
    if(n/=configured_species .or. n<1) error stop 'transport export before configuration'
    controls=[reference_dynamic_viscosity,prandtl_number]
    data=species_diffusivities(1:n)
  end subroutine

  subroutine validate_mc_transport_provider(requested_model)
    character(len=*), intent(in) :: requested_model

    if (trim(adjustl(requested_model)) /= mc_transport_provider_name) then
      write(*,'(A,A,A,A)') 'ERROR: requested transport provider "', &
        trim(requested_model), '" but this executable contains "', &
        mc_transport_provider_name
      error stop 'multicomponent transport provider mismatch'
    end if
  end subroutine validate_mc_transport_provider

  subroutine configure_mc_transport(path,nspecies,species_names)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: species_names(:)
    integer :: unit, ios, species
    character(len=512) :: message
    character(len=transport_name_length) :: &
      transport_species_names(mc_max_species)
    namelist /mixture_averaged_transport/ transport_species_names, &
      reference_dynamic_viscosity, prandtl_number, species_diffusivities

    if (nspecies < 1 .or. nspecies > mc_max_species .or. &
        size(species_names) < nspecies) then
      error stop 'invalid species contract for mixture-averaged transport'
    end if
    configured_species = 0
    transport_species_names = ''
    reference_dynamic_viscosity = 1.8e-5_dp
    prandtl_number = 0.72_dp
    species_diffusivities = 0.0_dp

    open(newunit=unit,file=trim(path),status='old',action='read', &
      iostat=ios,iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open transport input: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open mixture-averaged transport input'
    end if
    read(unit,nml=mixture_averaged_transport,iostat=ios,iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: invalid mixture_averaged_transport namelist: ', &
        trim(message)
      error stop 'failed to read mixture-averaged transport input'
    end if

    do species = 1, nspecies
      if (trim(transport_species_names(species)) /= &
          trim(species_names(species))) then
        error stop 'transport species order does not match physics species'
      end if
    end do
    if (.not. ieee_is_finite(reference_dynamic_viscosity) .or. &
        reference_dynamic_viscosity <= 0.0_dp) then
      error stop 'reference dynamic viscosity must be finite and positive'
    end if
    if (.not. ieee_is_finite(prandtl_number) .or. &
        prandtl_number <= 0.0_dp) then
      error stop 'Prandtl number must be finite and positive'
    end if
    if (.not. all(ieee_is_finite(species_diffusivities(1:nspecies))) .or. &
        minval(species_diffusivities(1:nspecies)) <= 0.0_dp) then
      error stop 'species diffusivities must be finite and positive'
    end if
    configured_species = nspecies
  end subroutine configure_mc_transport

  subroutine ensure_transport_configured(layout,mass_fractions)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: mass_fractions(:)

    if (configured_species == 0) then
      error stop 'mixture-averaged transport has not been configured'
    end if
    if (layout%nspecies /= configured_species .or. &
        size(mass_fractions) < layout%nspecies) then
      error stop 'state layout and transport species contract differ'
    end if
  end subroutine ensure_transport_configured

  real(dp) function mc_dynamic_viscosity( &
      mass_fractions,layout,temperature) result(viscosity)
    real(dp), intent(in) :: mass_fractions(:), temperature
    type(mc_state_layout), intent(in) :: layout

    call ensure_transport_configured(layout,mass_fractions)
    if (temperature <= 0.0_dp) error stop 'transport temperature must be positive'
    viscosity = reference_dynamic_viscosity
  end function mc_dynamic_viscosity

  real(dp) function mc_thermal_conductivity( &
      mass_fractions,layout,temperature,mixture_cp) result(conductivity)
    real(dp), intent(in) :: mass_fractions(:), temperature, mixture_cp
    type(mc_state_layout), intent(in) :: layout

    call ensure_transport_configured(layout,mass_fractions)
    if (temperature <= 0.0_dp .or. mixture_cp <= 0.0_dp) then
      error stop 'thermal-conductivity state is invalid'
    end if
    conductivity = reference_dynamic_viscosity*mixture_cp/prandtl_number
  end function mc_thermal_conductivity

  subroutine mc_species_diffusivities( &
      mass_fractions,layout,temperature,diffusivities)
    real(dp), intent(in) :: mass_fractions(:), temperature
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(out) :: diffusivities(:)

    call ensure_transport_configured(layout,mass_fractions)
    if (temperature <= 0.0_dp .or. size(diffusivities) < layout%nspecies) then
      error stop 'species-diffusivity request is invalid'
    end if
    diffusivities(1:layout%nspecies) = &
      species_diffusivities(1:layout%nspecies)
  end subroutine mc_species_diffusivities

end module mod_mc_transport_provider
