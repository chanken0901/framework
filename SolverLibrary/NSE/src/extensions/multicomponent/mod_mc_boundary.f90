module mod_mc_boundary
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, &
    mc_face_x_min, mc_face_x_max, &
    mc_face_y_min, mc_face_y_max, mc_face_z_min, mc_face_z_max
  use mod_mc_euler_field, only : set_mc_euler_conservative_state
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_mixture_gamma, mc_pressure, &
    mc_temperature, mc_sound_speed
  implicit none
  private

  public :: mc_boundary_state
  public :: mc_all_boundaries_periodic
  public :: mc_boundary_is_periodic
  public :: mc_boundary_face_geometry

contains

  pure logical function mc_boundary_is_periodic(config,face) result(periodic)
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: face

    periodic = trim(config%boundary_face_types(face)) == 'periodic'
  end function mc_boundary_is_periodic

  pure logical function mc_all_boundaries_periodic(config) result(periodic)
    type(mc_euler_config), intent(in) :: config

    periodic = all(config%boundary_face_types == 'periodic')
  end function mc_all_boundaries_periodic

  subroutine mc_boundary_face_geometry(face,normal_axis,outward_sign)
    integer, intent(in) :: face
    integer, intent(out) :: normal_axis
    real(dp), intent(out) :: outward_sign

    select case (face)
    case (mc_face_x_min)
      normal_axis = 1
      outward_sign = -1.0_dp
    case (mc_face_x_max)
      normal_axis = 1
      outward_sign = 1.0_dp
    case (mc_face_y_min)
      normal_axis = 2
      outward_sign = -1.0_dp
    case (mc_face_y_max)
      normal_axis = 2
      outward_sign = 1.0_dp
    case (mc_face_z_min)
      normal_axis = 3
      outward_sign = -1.0_dp
    case (mc_face_z_max)
      normal_axis = 3
      outward_sign = 1.0_dp
    case default
      error stop 'invalid multicomponent boundary face index'
    end select
  end subroutine mc_boundary_face_geometry

  subroutine mc_boundary_state(interior,ghost,layout,config,face)
    real(dp), intent(in) :: interior(:)
    real(dp), intent(out) :: ghost(:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: face
    integer :: normal_axis
    real(dp) :: outward_sign

    call mc_boundary_face_geometry(face,normal_axis,outward_sign)
    select case (trim(config%boundary_face_types(face)))
    case ('reflective')
      ghost = interior
      ghost(layout%momentum(normal_axis)) = &
        -ghost(layout%momentum(normal_axis))
    case ('dirichlet')
      call reference_conservative_state(ghost,layout,config,face)
    case ('non_reflecting')
      call characteristic_relaxation_state( &
        interior,ghost,layout,config,face,normal_axis,outward_sign)
    case ('periodic')
      error stop 'periodic multicomponent boundary needs the opposite cell'
    case default
      error stop 'unsupported multicomponent boundary face type'
    end select
  end subroutine mc_boundary_state

  subroutine reference_conservative_state(state,layout,config,face)
    real(dp), intent(out) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: face

    call set_mc_euler_conservative_state( &
      state,layout,config%gamma,config%boundary_reference_densities(face), &
      config%boundary_reference_velocities(:,face), &
      config%boundary_reference_pressures(face), &
      config%boundary_reference_mass_fractions(face,1:layout%nspecies))
  end subroutine reference_conservative_state

  subroutine characteristic_relaxation_state( &
      interior,ghost,layout,config,face,normal_axis,outward_sign)
    real(dp), intent(in) :: interior(:)
    real(dp), intent(out) :: ghost(:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: face, normal_axis
    real(dp), intent(in) :: outward_sign
    real(dp) :: density_i, pressure_i, temperature_i, velocity_i(3)
    real(dp) :: density_r, pressure_r, temperature_r, velocity_r(3)
    real(dp) :: normal_i, normal_r, normal_b, velocity_b(3)
    real(dp) :: sound_i, sound_r, sound_b, gamma_i, gamma_r, gamma_frozen
    real(dp) :: jminus_i, jplus_i, jminus_r, jplus_r
    real(dp) :: jminus_b, jplus_b, entropy_i, entropy_r, entropy_b
    real(dp) :: density_b, pressure_b, alpha, spacing, length_scale
    real(dp) :: mass_fractions_i(layout%nspecies)
    real(dp) :: mass_fractions_r(layout%nspecies)
    real(dp) :: mass_fractions_b(layout%nspecies)
    real(dp) :: reference_state(layout%nvariables)
    integer :: component

    density_i = mc_mixture_density(interior,layout)
    pressure_i = mc_pressure(interior,layout,config%gamma)
    temperature_i = mc_temperature(interior,layout,config%gamma)
    velocity_i = interior(layout%momentum)/density_i
    mass_fractions_i = &
      interior(layout%first_species:layout%last_species)/density_i
    call reference_conservative_state(reference_state,layout,config,face)
    if (maxval(abs(interior-reference_state)) <= &
        100.0_dp*epsilon(1.0_dp)* &
        max(maxval(abs(reference_state)),1.0_dp)) then
      ghost = interior
      return
    end if
    density_r = config%boundary_reference_densities(face)
    pressure_r = config%boundary_reference_pressures(face)
    velocity_r = config%boundary_reference_velocities(:,face)
    mass_fractions_r = &
      config%boundary_reference_mass_fractions(face,1:layout%nspecies)
    temperature_r = mc_temperature(reference_state,layout,config%gamma)
    sound_i = mc_sound_speed(interior,layout,config%gamma)
    sound_r = mc_sound_speed(reference_state,layout,config%gamma)
    gamma_i = mc_mixture_gamma( &
      mass_fractions_i,layout,temperature_i,config%gamma)
    gamma_r = mc_mixture_gamma( &
      mass_fractions_r,layout,temperature_r,config%gamma)
    gamma_frozen = 0.5_dp*(gamma_i+gamma_r)
    normal_i = outward_sign*velocity_i(normal_axis)
    normal_r = outward_sign*velocity_r(normal_axis)

    jminus_i = normal_i-2.0_dp*sound_i/(gamma_frozen-1.0_dp)
    jplus_i = normal_i+2.0_dp*sound_i/(gamma_frozen-1.0_dp)
    jminus_r = normal_r-2.0_dp*sound_r/(gamma_frozen-1.0_dp)
    jplus_r = normal_r+2.0_dp*sound_r/(gamma_frozen-1.0_dp)
    entropy_i = pressure_i/density_i**gamma_frozen
    entropy_r = pressure_r/density_r**gamma_frozen

    call boundary_spacing_and_length(config,normal_axis,spacing,length_scale)
    if (config%boundary_length_scale > 0.0_dp) then
      length_scale = config%boundary_length_scale
    end if
    alpha = 1.0_dp-exp(-config%boundary_relaxation_strength* &
      0.5_dp*spacing/length_scale)
    alpha = max(0.0_dp,min(1.0_dp,alpha))
    if (normal_i+sound_i < 0.0_dp) alpha = 1.0_dp

    jminus_b = outgoing_or_relaxed( &
      jminus_i,jminus_r,normal_i-sound_i,alpha)
    jplus_b = outgoing_or_relaxed( &
      jplus_i,jplus_r,normal_i+sound_i,alpha)
    entropy_b = outgoing_or_relaxed( &
      entropy_i,entropy_r,normal_i,alpha)
    velocity_b = velocity_i
    do component = 1, 3
      if (component == normal_axis) cycle
      velocity_b(component) = outgoing_or_relaxed( &
        velocity_i(component),velocity_r(component),normal_i,alpha)
    end do
    if (normal_i >= 0.0_dp) then
      mass_fractions_b = mass_fractions_i
    else
      mass_fractions_b = (1.0_dp-alpha)*mass_fractions_i+ &
        alpha*mass_fractions_r
    end if

    normal_b = 0.5_dp*(jplus_b+jminus_b)
    sound_b = 0.25_dp*(gamma_frozen-1.0_dp)*(jplus_b-jminus_b)
    if (.not. ieee_is_finite(sound_b) .or. sound_b <= 0.0_dp .or. &
        .not. ieee_is_finite(entropy_b) .or. entropy_b <= 0.0_dp) then
      error stop 'invalid multicomponent characteristic boundary state'
    end if
    density_b = (sound_b*sound_b/(gamma_frozen*entropy_b))** &
      (1.0_dp/(gamma_frozen-1.0_dp))
    pressure_b = entropy_b*density_b**gamma_frozen
    velocity_b(normal_axis) = outward_sign*normal_b
    call set_mc_euler_conservative_state( &
      ghost,layout,config%gamma,density_b,velocity_b,pressure_b, &
      mass_fractions_b)
  end subroutine characteristic_relaxation_state

  subroutine boundary_spacing_and_length( &
      config,normal_axis,spacing,length_scale)
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: normal_axis
    real(dp), intent(out) :: spacing, length_scale

    select case (normal_axis)
    case (1)
      length_scale = config%x_max-config%x_min
      spacing = length_scale/real(config%nx,dp)
    case (2)
      length_scale = config%y_max-config%y_min
      spacing = length_scale/real(config%ny,dp)
    case (3)
      length_scale = config%z_max-config%z_min
      spacing = length_scale/real(config%nz,dp)
    case default
      error stop 'invalid multicomponent boundary normal axis'
    end select
    if (config%global_boundary_lengths(normal_axis) > 0.0_dp) &
      length_scale=config%global_boundary_lengths(normal_axis)
  end subroutine boundary_spacing_and_length

  pure real(dp) function outgoing_or_relaxed( &
      interior,reference,eigenvalue,alpha) result(value)
    real(dp), intent(in) :: interior, reference, eigenvalue, alpha

    if (eigenvalue >= 0.0_dp) then
      value = interior
    else
      value = (1.0_dp-alpha)*interior+alpha*reference
    end if
  end function outgoing_or_relaxed

end module mod_mc_boundary
