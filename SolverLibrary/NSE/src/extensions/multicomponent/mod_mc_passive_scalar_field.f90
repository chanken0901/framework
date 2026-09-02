module mod_mc_passive_scalar_field
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_passive_scalar_config, only : mc_passive_scalar_config
  implicit none
  private

  public :: initialize_mc_passive_scalar_state
  public :: compute_mc_species_masses
  public :: compute_mc_species_sum_error
  public :: compute_mc_tracer_bounds

contains

  subroutine initialize_mc_passive_scalar_state(q, layout, config)
    type(mc_state_layout), intent(in) :: layout
    type(mc_passive_scalar_config), intent(in) :: config
    real(dp), intent(out) :: q(:,:,:,:)
    integer :: i, j, k, carrier
    real(dp) :: x, y, z, radius_squared, tracer
    real(dp) :: dx, dy, dz, kinetic_energy

    if (layout%nspecies < 2) then
      error stop 'passive-scalar mode requires tracer and carrier species'
    end if
    if (size(q,1) /= config%nx .or. size(q,2) /= config%ny .or. &
        size(q,3) /= config%nz .or. size(q,4) /= layout%nvariables) then
      error stop 'passive-scalar state allocation does not match configuration'
    end if

    dx = (config%x_max-config%x_min) / real(config%nx, dp)
    dy = (config%y_max-config%y_min) / real(config%ny, dp)
    dz = (config%z_max-config%z_min) / real(config%nz, dp)
    kinetic_energy = 0.5_dp * sum(config%velocity**2)
    carrier = layout%last_species
    q = 0.0_dp
    do k = 1, config%nz
      z = config%z_min + (real(k,dp)-0.5_dp)*dz
      do j = 1, config%ny
        y = config%y_min + (real(j,dp)-0.5_dp)*dy
        do i = 1, config%nx
          x = config%x_min + (real(i,dp)-0.5_dp)*dx
          radius_squared = periodic_distance_squared(x, y, z, config)
          tracer = config%tracer_background + config%tracer_amplitude * &
            exp(-0.5_dp*radius_squared/config%tracer_width**2)
          q(i,j,k,layout%first_species) = tracer
          q(i,j,k,carrier) = 1.0_dp - tracer
          q(i,j,k,layout%momentum) = config%velocity
          q(i,j,k,layout%total_energy) = 1.0_dp + kinetic_energy
        end do
      end do
    end do
  end subroutine initialize_mc_passive_scalar_state

  subroutine compute_mc_species_masses(q, layout, config, masses)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_passive_scalar_config), intent(in) :: config
    real(dp), intent(out) :: masses(:)
    integer :: species, variable
    real(dp) :: cell_volume

    if (size(masses) /= layout%nspecies) then
      error stop 'species-mass diagnostic has the wrong size'
    end if
    cell_volume = (config%x_max-config%x_min) / real(config%nx,dp) * &
      (config%y_max-config%y_min) / real(config%ny,dp) * &
      (config%z_max-config%z_min) / real(config%nz,dp)
    do species = 1, layout%nspecies
      variable = layout%first_species + species - 1
      masses(species) = sum(q(:,:,:,variable)) * cell_volume
    end do
  end subroutine compute_mc_species_masses

  real(dp) function compute_mc_species_sum_error(q, layout) result(error)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    integer :: i, j, k
    real(dp) :: density

    error = 0.0_dp
    do k = 1, size(q,3)
      do j = 1, size(q,2)
        do i = 1, size(q,1)
          density = sum(q(i,j,k,layout%first_species:layout%last_species))
          error = max(error, abs(density-1.0_dp))
        end do
      end do
    end do
  end function compute_mc_species_sum_error

  subroutine compute_mc_tracer_bounds(q, layout, minimum, maximum)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(out) :: minimum, maximum

    minimum = minval(q(:,:,:,layout%first_species))
    maximum = maxval(q(:,:,:,layout%first_species))
  end subroutine compute_mc_tracer_bounds

  pure real(dp) function periodic_distance_squared(x, y, z, config) &
      result(distance_squared)
    real(dp), intent(in) :: x, y, z
    type(mc_passive_scalar_config), intent(in) :: config
    real(dp) :: distance(3), length(3)

    length = [config%x_max-config%x_min, config%y_max-config%y_min, &
      config%z_max-config%z_min]
    distance = abs([x,y,z] - config%tracer_center)
    distance = min(distance, length-distance)
    distance_squared = sum(distance**2)
  end function periodic_distance_squared

end module mod_mc_passive_scalar_field
