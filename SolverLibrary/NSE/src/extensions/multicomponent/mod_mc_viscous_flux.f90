module mod_mc_viscous_flux
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, &
    mc_face_x_min, mc_face_x_max, mc_face_y_min, mc_face_y_max, &
    mc_face_z_min, mc_face_z_max
  use mod_mc_euler_field, only : validate_mc_euler_state, &
    mc_primitive_workspace, prepare_mc_primitives
  use mod_mc_euler_flux, only : compute_mc_euler_rhs_prepared, &
    compute_mc_euler_timestep
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_mixture_gas_constant, mc_mixture_cp, mc_temperature, &
    mc_species_enthalpies
  use mod_mc_transport_provider, only : mc_dynamic_viscosity, &
    mc_thermal_conductivity, mc_species_diffusivities
  use mod_mc_boundary, only : mc_boundary_state, mc_boundary_is_periodic
  implicit none
  private

  public :: compute_mc_transport_rhs
  public :: compute_mc_navier_stokes_rhs
  public :: compute_mc_navier_stokes_timestep
  public :: advance_mc_navier_stokes_ssprk3

  type, public :: mc_navier_stokes_workspace
    type(mc_primitive_workspace) :: primitive
    real(dp), allocatable :: velocity_gradient(:,:,:,:,:)
    real(dp), allocatable :: temperature_gradient(:,:,:,:)
    real(dp), allocatable :: fraction_gradient(:,:,:,:,:)
  end type mc_navier_stokes_workspace

contains

  subroutine prepare_workspace(q,layout,config,work)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_navier_stokes_workspace), intent(inout) :: work
    integer :: nx,ny,nz
    nx=config%nx
    ny=config%ny
    nz=config%nz
    call prepare_mc_primitives(q,layout,config,work%primitive)
    if (allocated(work%fraction_gradient)) then
      if (any(shape(work%fraction_gradient) /= [nx,ny,nz,layout%nspecies,3])) then
        deallocate(work%velocity_gradient,work%temperature_gradient, &
          work%fraction_gradient)
      end if
    end if
    if (.not. allocated(work%fraction_gradient)) then
      allocate(work%velocity_gradient(nx,ny,nz,3,3), &
        work%temperature_gradient(nx,ny,nz,3), &
        work%fraction_gradient(nx,ny,nz,layout%nspecies,3))
    end if
  end subroutine prepare_workspace

  subroutine compute_mc_transport_rhs(q,rhs,layout,config,workspace)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_navier_stokes_workspace), intent(inout), optional, target :: workspace
    type(mc_navier_stokes_workspace), target :: local_workspace
    type(mc_navier_stokes_workspace), pointer :: work
    work => local_workspace
    if (present(workspace)) work => workspace
    call prepare_workspace(q,layout,config,work)
    rhs=0.0_dp
    call add_transport_rhs(q,rhs,layout,config,work)
  end subroutine compute_mc_transport_rhs

  subroutine add_transport_rhs(q,rhs,layout,config,work)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(inout) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_navier_stokes_workspace), intent(inout) :: work
    real(dp) :: dx, dy, dz
    real(dp) :: face_flux(layout%nvariables)
    integer :: i, j, k, ip, jp, kp

    if (any(shape(rhs) /= shape(q))) then
      error stop 'transport RHS allocation does not match state'
    end if
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    associate(density => work%primitive%density, &
      temperature => work%primitive%temperature, &
      velocity => work%primitive%velocity, &
      mass_fractions => work%primitive%mass_fractions, &
      velocity_gradient => work%velocity_gradient, &
      temperature_gradient => work%temperature_gradient, &
      fraction_gradient => work%fraction_gradient)
    call compute_mc_gradients( &
      q,layout,config,velocity,temperature,mass_fractions,dx,dy,dz, &
      velocity_gradient,temperature_gradient,fraction_gradient)

    do k = 1, config%nz
      kp = merge(1,k+1,k == config%nz)
      do j = 1, config%ny
        jp = merge(1,j+1,j == config%ny)
        do i = 1, config%nx
          ip = merge(1,i+1,i == config%nx)
          if (i < config%nx .or. &
              mc_boundary_is_periodic(config,mc_face_x_max)) then
            call compute_transport_face_flux( &
              0.5_dp*(density(i,j,k)+density(ip,j,k)), &
              0.5_dp*(velocity(i,j,k,:)+velocity(ip,j,k,:)), &
              0.5_dp*(temperature(i,j,k)+temperature(ip,j,k)), &
              0.5_dp*(mass_fractions(i,j,k,:)+mass_fractions(ip,j,k,:)), &
              0.5_dp*(velocity_gradient(i,j,k,:,:)+ &
                velocity_gradient(ip,j,k,:,:)), &
              0.5_dp*(temperature_gradient(i,j,k,:)+ &
                temperature_gradient(ip,j,k,:)), &
              0.5_dp*(fraction_gradient(i,j,k,:,:)+ &
                fraction_gradient(ip,j,k,:,:)),layout,1,face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)+face_flux/dx
            rhs(ip,j,k,:) = rhs(ip,j,k,:)-face_flux/dx
          else
            call compute_boundary_transport_face_flux( &
              q(i,j,k,:),density(i,j,k),velocity(i,j,k,:), &
              temperature(i,j,k),mass_fractions(i,j,k,:), &
              velocity_gradient(i,j,k,:,:),temperature_gradient(i,j,k,:), &
              fraction_gradient(i,j,k,:,:),layout,config,mc_face_x_max,1, &
              face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)+face_flux/dx
          end if
          if (i == 1 .and. &
              .not. mc_boundary_is_periodic(config,mc_face_x_min)) then
            call compute_boundary_transport_face_flux( &
              q(i,j,k,:),density(i,j,k),velocity(i,j,k,:), &
              temperature(i,j,k),mass_fractions(i,j,k,:), &
              velocity_gradient(i,j,k,:,:),temperature_gradient(i,j,k,:), &
              fraction_gradient(i,j,k,:,:),layout,config,mc_face_x_min,1, &
              face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)-face_flux/dx
          end if

          if (j < config%ny .or. &
              mc_boundary_is_periodic(config,mc_face_y_max)) then
            call compute_transport_face_flux( &
              0.5_dp*(density(i,j,k)+density(i,jp,k)), &
              0.5_dp*(velocity(i,j,k,:)+velocity(i,jp,k,:)), &
              0.5_dp*(temperature(i,j,k)+temperature(i,jp,k)), &
              0.5_dp*(mass_fractions(i,j,k,:)+mass_fractions(i,jp,k,:)), &
              0.5_dp*(velocity_gradient(i,j,k,:,:)+ &
                velocity_gradient(i,jp,k,:,:)), &
              0.5_dp*(temperature_gradient(i,j,k,:)+ &
                temperature_gradient(i,jp,k,:)), &
              0.5_dp*(fraction_gradient(i,j,k,:,:)+ &
                fraction_gradient(i,jp,k,:,:)),layout,2,face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)+face_flux/dy
            rhs(i,jp,k,:) = rhs(i,jp,k,:)-face_flux/dy
          else
            call compute_boundary_transport_face_flux( &
              q(i,j,k,:),density(i,j,k),velocity(i,j,k,:), &
              temperature(i,j,k),mass_fractions(i,j,k,:), &
              velocity_gradient(i,j,k,:,:),temperature_gradient(i,j,k,:), &
              fraction_gradient(i,j,k,:,:),layout,config,mc_face_y_max,2, &
              face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)+face_flux/dy
          end if
          if (j == 1 .and. &
              .not. mc_boundary_is_periodic(config,mc_face_y_min)) then
            call compute_boundary_transport_face_flux( &
              q(i,j,k,:),density(i,j,k),velocity(i,j,k,:), &
              temperature(i,j,k),mass_fractions(i,j,k,:), &
              velocity_gradient(i,j,k,:,:),temperature_gradient(i,j,k,:), &
              fraction_gradient(i,j,k,:,:),layout,config,mc_face_y_min,2, &
              face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)-face_flux/dy
          end if

          if (k < config%nz .or. &
              mc_boundary_is_periodic(config,mc_face_z_max)) then
            call compute_transport_face_flux( &
              0.5_dp*(density(i,j,k)+density(i,j,kp)), &
              0.5_dp*(velocity(i,j,k,:)+velocity(i,j,kp,:)), &
              0.5_dp*(temperature(i,j,k)+temperature(i,j,kp)), &
              0.5_dp*(mass_fractions(i,j,k,:)+mass_fractions(i,j,kp,:)), &
              0.5_dp*(velocity_gradient(i,j,k,:,:)+ &
                velocity_gradient(i,j,kp,:,:)), &
              0.5_dp*(temperature_gradient(i,j,k,:)+ &
                temperature_gradient(i,j,kp,:)), &
              0.5_dp*(fraction_gradient(i,j,k,:,:)+ &
                fraction_gradient(i,j,kp,:,:)),layout,3,face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)+face_flux/dz
            rhs(i,j,kp,:) = rhs(i,j,kp,:)-face_flux/dz
          else
            call compute_boundary_transport_face_flux( &
              q(i,j,k,:),density(i,j,k),velocity(i,j,k,:), &
              temperature(i,j,k),mass_fractions(i,j,k,:), &
              velocity_gradient(i,j,k,:,:),temperature_gradient(i,j,k,:), &
              fraction_gradient(i,j,k,:,:),layout,config,mc_face_z_max,3, &
              face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)+face_flux/dz
          end if
          if (k == 1 .and. &
              .not. mc_boundary_is_periodic(config,mc_face_z_min)) then
            call compute_boundary_transport_face_flux( &
              q(i,j,k,:),density(i,j,k),velocity(i,j,k,:), &
              temperature(i,j,k),mass_fractions(i,j,k,:), &
              velocity_gradient(i,j,k,:,:),temperature_gradient(i,j,k,:), &
              fraction_gradient(i,j,k,:,:),layout,config,mc_face_z_min,3, &
              face_flux)
            rhs(i,j,k,:) = rhs(i,j,k,:)-face_flux/dz
          end if
        end do
      end do
    end do
    end associate
  end subroutine add_transport_rhs

  subroutine compute_mc_gradients( &
      q,layout,config,velocity,temperature,mass_fractions,dx,dy,dz, &
      velocity_gradient,temperature_gradient,fraction_gradient)
    real(dp), intent(in) :: q(:,:,:,:), velocity(:,:,:,:)
    real(dp), intent(in) :: temperature(:,:,:), mass_fractions(:,:,:,:)
    real(dp), intent(in) :: dx, dy, dz
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    real(dp), intent(out) :: velocity_gradient(:,:,:,:,:)
    real(dp), intent(out) :: temperature_gradient(:,:,:,:)
    real(dp), intent(out) :: fraction_gradient(:,:,:,:,:)
    real(dp) :: velocity_minus(3), velocity_plus(3)
    real(dp) :: temperature_minus, temperature_plus
    real(dp) :: fractions_minus(layout%nspecies)
    real(dp) :: fractions_plus(layout%nspecies)
    integer :: i, j, k

    do k = 1, size(temperature,3)
      do j = 1, size(temperature,2)
        do i = 1, size(temperature,1)
          call sample_mc_neighbor_primitive( &
            q,layout,config,velocity,temperature,mass_fractions, &
            i,j,k,1,-1,velocity_minus,temperature_minus, &
            fractions_minus)
          call sample_mc_neighbor_primitive( &
            q,layout,config,velocity,temperature,mass_fractions, &
            i,j,k,1,1,velocity_plus,temperature_plus,fractions_plus)
          velocity_gradient(i,j,k,:,1) = &
            (velocity_plus-velocity_minus)/(2.0_dp*dx)
          temperature_gradient(i,j,k,1) = &
            (temperature_plus-temperature_minus)/(2.0_dp*dx)
          fraction_gradient(i,j,k,:,1) = &
            (fractions_plus-fractions_minus)/(2.0_dp*dx)

          call sample_mc_neighbor_primitive( &
            q,layout,config,velocity,temperature,mass_fractions, &
            i,j,k,2,-1,velocity_minus,temperature_minus, &
            fractions_minus)
          call sample_mc_neighbor_primitive( &
            q,layout,config,velocity,temperature,mass_fractions, &
            i,j,k,2,1,velocity_plus,temperature_plus,fractions_plus)
          velocity_gradient(i,j,k,:,2) = &
            (velocity_plus-velocity_minus)/(2.0_dp*dy)
          temperature_gradient(i,j,k,2) = &
            (temperature_plus-temperature_minus)/(2.0_dp*dy)
          fraction_gradient(i,j,k,:,2) = &
            (fractions_plus-fractions_minus)/(2.0_dp*dy)

          call sample_mc_neighbor_primitive( &
            q,layout,config,velocity,temperature,mass_fractions, &
            i,j,k,3,-1,velocity_minus,temperature_minus, &
            fractions_minus)
          call sample_mc_neighbor_primitive( &
            q,layout,config,velocity,temperature,mass_fractions, &
            i,j,k,3,1,velocity_plus,temperature_plus,fractions_plus)
          velocity_gradient(i,j,k,:,3) = &
            (velocity_plus-velocity_minus)/(2.0_dp*dz)
          temperature_gradient(i,j,k,3) = &
            (temperature_plus-temperature_minus)/(2.0_dp*dz)
          fraction_gradient(i,j,k,:,3) = &
            (fractions_plus-fractions_minus)/(2.0_dp*dz)
        end do
      end do
    end do
  end subroutine compute_mc_gradients

  subroutine sample_mc_neighbor_primitive( &
      q,layout,config,velocity,temperature,mass_fractions, &
      i,j,k,direction,offset,neighbor_velocity,neighbor_temperature, &
      neighbor_fractions)
    real(dp), intent(in) :: q(:,:,:,:), velocity(:,:,:,:)
    real(dp), intent(in) :: temperature(:,:,:), mass_fractions(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: i, j, k, direction, offset
    real(dp), intent(out) :: neighbor_velocity(3), neighbor_temperature
    real(dp), intent(out) :: neighbor_fractions(layout%nspecies)
    real(dp) :: ghost(layout%nvariables), ghost_density
    integer :: ni, nj, nk, face

    ni = i
    nj = j
    nk = k
    select case (direction)
    case (1)
      ni = i+offset
      face = merge(mc_face_x_min,mc_face_x_max,offset < 0)
      if (ni < 1) ni = config%nx
      if (ni > config%nx) ni = 1
    case (2)
      nj = j+offset
      face = merge(mc_face_y_min,mc_face_y_max,offset < 0)
      if (nj < 1) nj = config%ny
      if (nj > config%ny) nj = 1
    case (3)
      nk = k+offset
      face = merge(mc_face_z_min,mc_face_z_max,offset < 0)
      if (nk < 1) nk = config%nz
      if (nk > config%nz) nk = 1
    case default
      error stop 'invalid multicomponent gradient direction'
    end select

    if ((direction == 1 .and. i+offset >= 1 .and. &
         i+offset <= config%nx) .or. &
        (direction == 2 .and. j+offset >= 1 .and. &
         j+offset <= config%ny) .or. &
        (direction == 3 .and. k+offset >= 1 .and. &
         k+offset <= config%nz) .or. mc_boundary_is_periodic(config,face)) then
      neighbor_velocity = velocity(ni,nj,nk,:)
      neighbor_temperature = temperature(ni,nj,nk)
      neighbor_fractions = mass_fractions(ni,nj,nk,:)
    else
      call mc_boundary_state(q(i,j,k,:),ghost,layout,config,face)
      ghost_density = mc_mixture_density(ghost,layout)
      neighbor_velocity = ghost(layout%momentum)/ghost_density
      neighbor_temperature = mc_temperature(ghost,layout,config%gamma)
      neighbor_fractions = &
        ghost(layout%first_species:layout%last_species)/ghost_density
    end if
  end subroutine sample_mc_neighbor_primitive

  subroutine compute_boundary_transport_face_flux( &
      interior_state,interior_density,interior_velocity, &
      interior_temperature,interior_fractions,interior_velocity_gradient, &
      interior_temperature_gradient,interior_fraction_gradient, &
      layout,config,face,direction,flux)
    real(dp), intent(in) :: interior_state(:), interior_density
    real(dp), intent(in) :: interior_velocity(3), interior_temperature
    real(dp), intent(in) :: interior_fractions(:)
    real(dp), intent(in) :: interior_velocity_gradient(3,3)
    real(dp), intent(in) :: interior_temperature_gradient(3)
    real(dp), intent(in) :: interior_fraction_gradient(:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer, intent(in) :: face, direction
    real(dp), intent(out) :: flux(:)
    real(dp) :: ghost(layout%nvariables), ghost_density
    real(dp) :: ghost_velocity(3), ghost_temperature
    real(dp) :: ghost_fractions(layout%nspecies)
    real(dp) :: face_velocity_gradient(3,3)
    real(dp) :: face_temperature_gradient(3)
    real(dp) :: face_fraction_gradient(layout%nspecies,3)
    real(dp) :: spacing, orientation

    call mc_boundary_state(interior_state,ghost,layout,config,face)
    ghost_density = mc_mixture_density(ghost,layout)
    ghost_velocity = ghost(layout%momentum)/ghost_density
    ghost_temperature = mc_temperature(ghost,layout,config%gamma)
    ghost_fractions = &
      ghost(layout%first_species:layout%last_species)/ghost_density
    select case (direction)
    case (1)
      spacing = (config%x_max-config%x_min)/real(config%nx,dp)
    case (2)
      spacing = (config%y_max-config%y_min)/real(config%ny,dp)
    case (3)
      spacing = (config%z_max-config%z_min)/real(config%nz,dp)
    case default
      error stop 'invalid multicomponent transport boundary direction'
    end select
    if (face == mc_face_x_min .or. face == mc_face_y_min .or. &
        face == mc_face_z_min) then
      orientation = 1.0_dp
    else
      orientation = -1.0_dp
    end if
    face_velocity_gradient = interior_velocity_gradient
    face_temperature_gradient = interior_temperature_gradient
    face_fraction_gradient = interior_fraction_gradient
    face_velocity_gradient(:,direction) = orientation* &
      (interior_velocity-ghost_velocity)/spacing
    face_temperature_gradient(direction) = orientation* &
      (interior_temperature-ghost_temperature)/spacing
    face_fraction_gradient(:,direction) = orientation* &
      (interior_fractions(1:layout%nspecies)-ghost_fractions)/spacing

    call compute_transport_face_flux( &
      0.5_dp*(interior_density+ghost_density), &
      0.5_dp*(interior_velocity+ghost_velocity), &
      0.5_dp*(interior_temperature+ghost_temperature), &
      0.5_dp*(interior_fractions(1:layout%nspecies)+ghost_fractions), &
      face_velocity_gradient,face_temperature_gradient, &
      face_fraction_gradient,layout,direction,flux)
  end subroutine compute_boundary_transport_face_flux

  subroutine compute_transport_face_flux( &
      density,velocity,temperature,mass_fractions,velocity_gradient, &
      temperature_gradient,fraction_gradient,layout,direction,flux)
    real(dp), intent(in) :: density, velocity(3), temperature
    real(dp), intent(in) :: mass_fractions(:)
    real(dp), intent(in) :: velocity_gradient(3,3)
    real(dp), intent(in) :: temperature_gradient(3)
    real(dp), intent(in) :: fraction_gradient(:,:)
    type(mc_state_layout), intent(in) :: layout
    integer, intent(in) :: direction
    real(dp), intent(out) :: flux(:)
    real(dp) :: viscosity, conductivity, cp_value, divergence
    real(dp) :: diffusivities(layout%nspecies)
    real(dp) :: preliminary_flux(layout%nspecies)
    real(dp) :: diffusion_flux(layout%nspecies)
    real(dp) :: enthalpies(layout%nspecies), stress_normal(3)
    real(dp) :: correction
    integer :: species, component

    viscosity = mc_dynamic_viscosity( &
      mass_fractions,layout,temperature)
    cp_value = mc_mixture_cp(mass_fractions,layout,temperature)
    conductivity = mc_thermal_conductivity( &
      mass_fractions,layout,temperature,cp_value)
    call mc_species_diffusivities( &
      mass_fractions,layout,temperature,diffusivities)
    call mc_species_enthalpies(layout,temperature,enthalpies)

    do species = 1, layout%nspecies
      preliminary_flux(species) = -density*diffusivities(species)* &
        fraction_gradient(species,direction)
    end do
    correction = sum(preliminary_flux)
    diffusion_flux = preliminary_flux-mass_fractions*correction
    diffusion_flux(layout%nspecies) = &
      diffusion_flux(layout%nspecies)-sum(diffusion_flux)

    divergence = velocity_gradient(1,1) + velocity_gradient(2,2) + &
      velocity_gradient(3,3)
    do component = 1, 3
      stress_normal(component) = viscosity*( &
        velocity_gradient(component,direction) + &
        velocity_gradient(direction,component))
      if (component == direction) then
        stress_normal(component) = stress_normal(component) - &
          (2.0_dp/3.0_dp)*viscosity*divergence
      end if
    end do

    flux = 0.0_dp
    flux(layout%first_species:layout%last_species) = -diffusion_flux
    flux(layout%momentum) = stress_normal
    flux(layout%total_energy) = dot_product(stress_normal,velocity) + &
      conductivity*temperature_gradient(direction) - &
      dot_product(enthalpies,diffusion_flux)
  end subroutine compute_transport_face_flux

  subroutine compute_mc_navier_stokes_rhs(q,rhs,layout,config,workspace)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_navier_stokes_workspace), intent(inout), optional, target :: workspace
    type(mc_navier_stokes_workspace), target :: local_workspace
    type(mc_navier_stokes_workspace), pointer :: work

    work => local_workspace
    if (present(workspace)) work => workspace
    call prepare_workspace(q,layout,config,work)
    call compute_mc_euler_rhs_prepared(q,rhs,layout,config,work%primitive)
    call add_transport_rhs(q,rhs,layout,config,work)
  end subroutine compute_mc_navier_stokes_rhs

  real(dp) function compute_mc_navier_stokes_timestep( &
      q,layout,config) result(dt)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    real(dp) :: convective_dt, transport_dt, maximum_rate, rate
    real(dp) :: density, temperature, cp_value, gas_constant
    real(dp) :: viscosity, conductivity, maximum_diffusivity
    real(dp) :: dx, dy, dz
    real(dp) :: mass_fractions(layout%nspecies)
    real(dp) :: diffusivities(layout%nspecies)
    integer :: i, j, k

    call validate_mc_euler_state(q,layout,config)
    convective_dt = compute_mc_euler_timestep(q,layout,config)
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    maximum_rate = 0.0_dp
    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          density = mc_mixture_density(q(i,j,k,:),layout)
          temperature = mc_temperature(q(i,j,k,:),layout,config%gamma)
          mass_fractions = &
            q(i,j,k,layout%first_species:layout%last_species)/density
          cp_value = mc_mixture_cp(mass_fractions,layout,temperature)
          gas_constant = mc_mixture_gas_constant(mass_fractions,layout)
          viscosity = mc_dynamic_viscosity( &
            mass_fractions,layout,temperature)
          conductivity = mc_thermal_conductivity( &
            mass_fractions,layout,temperature,cp_value)
          call mc_species_diffusivities( &
            mass_fractions,layout,temperature,diffusivities)
          maximum_diffusivity = max( &
            maxval(diffusivities),viscosity/density, &
            conductivity/(density*(cp_value-gas_constant)))
          rate = 2.0_dp*maximum_diffusivity*( &
            1.0_dp/dx**2+1.0_dp/dy**2+1.0_dp/dz**2)
          maximum_rate = max(maximum_rate,rate)
        end do
      end do
    end do
    if (maximum_rate <= 0.0_dp) then
      error stop 'transport stability rate must be positive'
    end if
    transport_dt = config%diffusion_cfl/maximum_rate
    if (config%dt > 0.0_dp) then
      if (config%dt > transport_dt*(1.0_dp+100.0_dp*epsilon(1.0_dp))) then
        error stop 'fixed multicomponent dt violates the diffusion limit'
      end if
      dt = config%dt
    else
      dt = min(convective_dt,transport_dt)
    end if
  end function compute_mc_navier_stokes_timestep

  subroutine advance_mc_navier_stokes_ssprk3( &
      q,q0,rhs,dt,layout,config,workspace)
    real(dp), intent(inout) :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_navier_stokes_workspace), intent(inout), optional, target :: workspace
    type(mc_navier_stokes_workspace), target :: local_workspace
    type(mc_navier_stokes_workspace), pointer :: work

    if (any(shape(q0) /= shape(q)) .or. any(shape(rhs) /= shape(q))) then
      error stop 'multicomponent Navier-Stokes work arrays do not match state'
    end if
    work => local_workspace
    if (present(workspace)) work => workspace
    q0 = q
    call compute_mc_navier_stokes_rhs(q,rhs,layout,config,work)
    q = q0+dt*rhs
    call compute_mc_navier_stokes_rhs(q,rhs,layout,config,work)
    q = 0.75_dp*q0+0.25_dp*(q+dt*rhs)
    call compute_mc_navier_stokes_rhs(q,rhs,layout,config,work)
    q = (1.0_dp/3.0_dp)*q0+(2.0_dp/3.0_dp)*(q+dt*rhs)
    call validate_mc_euler_state(q,layout,config)
  end subroutine advance_mc_navier_stokes_ssprk3

end module mod_mc_viscous_flux
