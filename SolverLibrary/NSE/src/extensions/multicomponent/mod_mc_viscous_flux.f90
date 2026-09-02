module mod_mc_viscous_flux
  use mod_precision, only : dp
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_euler_field, only : validate_mc_euler_state
  use mod_mc_euler_flux, only : compute_mc_euler_rhs, &
    compute_mc_euler_timestep
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_mixture_gas_constant, mc_mixture_cp, mc_temperature, &
    mc_species_enthalpies
  use mod_mc_transport_provider, only : mc_dynamic_viscosity, &
    mc_thermal_conductivity, mc_species_diffusivities
  implicit none
  private

  public :: compute_mc_transport_rhs
  public :: compute_mc_navier_stokes_rhs
  public :: compute_mc_navier_stokes_timestep
  public :: advance_mc_navier_stokes_ssprk3

contains

  subroutine compute_mc_transport_rhs(q,rhs,layout,config)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    real(dp), allocatable :: density(:,:,:), temperature(:,:,:)
    real(dp), allocatable :: velocity(:,:,:,:), mass_fractions(:,:,:,:)
    real(dp), allocatable :: velocity_gradient(:,:,:,:,:)
    real(dp), allocatable :: temperature_gradient(:,:,:,:)
    real(dp), allocatable :: fraction_gradient(:,:,:,:,:)
    real(dp) :: dx, dy, dz
    real(dp) :: face_flux(layout%nvariables)
    integer :: i, j, k, ip, jp, kp

    if (any(shape(rhs) /= shape(q))) then
      error stop 'transport RHS allocation does not match state'
    end if
    call validate_mc_euler_state(q,layout,config)
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    allocate(density(config%nx,config%ny,config%nz))
    allocate(temperature(config%nx,config%ny,config%nz))
    allocate(velocity(config%nx,config%ny,config%nz,3))
    allocate(mass_fractions( &
      config%nx,config%ny,config%nz,layout%nspecies))
    allocate(velocity_gradient( &
      config%nx,config%ny,config%nz,3,3))
    allocate(temperature_gradient(config%nx,config%ny,config%nz,3))
    allocate(fraction_gradient( &
      config%nx,config%ny,config%nz,layout%nspecies,3))

    call compute_primitive_fields( &
      q,layout,config,density,velocity,temperature,mass_fractions)
    call compute_periodic_gradients( &
      velocity,temperature,mass_fractions,dx,dy,dz, &
      velocity_gradient,temperature_gradient,fraction_gradient)

    rhs = 0.0_dp
    do k = 1, config%nz
      kp = merge(1,k+1,k == config%nz)
      do j = 1, config%ny
        jp = merge(1,j+1,j == config%ny)
        do i = 1, config%nx
          ip = merge(1,i+1,i == config%nx)
          call compute_transport_face_flux( &
            0.5_dp*(density(i,j,k)+density(ip,j,k)), &
            0.5_dp*(velocity(i,j,k,:)+velocity(ip,j,k,:)), &
            0.5_dp*(temperature(i,j,k)+temperature(ip,j,k)), &
            0.5_dp*(mass_fractions(i,j,k,:)+ &
              mass_fractions(ip,j,k,:)), &
            0.5_dp*(velocity_gradient(i,j,k,:,:)+ &
              velocity_gradient(ip,j,k,:,:)), &
            0.5_dp*(temperature_gradient(i,j,k,:)+ &
              temperature_gradient(ip,j,k,:)), &
            0.5_dp*(fraction_gradient(i,j,k,:,:)+ &
              fraction_gradient(ip,j,k,:,:)),layout,1,face_flux)
          rhs(i,j,k,:) = rhs(i,j,k,:) + face_flux/dx
          rhs(ip,j,k,:) = rhs(ip,j,k,:) - face_flux/dx

          call compute_transport_face_flux( &
            0.5_dp*(density(i,j,k)+density(i,jp,k)), &
            0.5_dp*(velocity(i,j,k,:)+velocity(i,jp,k,:)), &
            0.5_dp*(temperature(i,j,k)+temperature(i,jp,k)), &
            0.5_dp*(mass_fractions(i,j,k,:)+ &
              mass_fractions(i,jp,k,:)), &
            0.5_dp*(velocity_gradient(i,j,k,:,:)+ &
              velocity_gradient(i,jp,k,:,:)), &
            0.5_dp*(temperature_gradient(i,j,k,:)+ &
              temperature_gradient(i,jp,k,:)), &
            0.5_dp*(fraction_gradient(i,j,k,:,:)+ &
              fraction_gradient(i,jp,k,:,:)),layout,2,face_flux)
          rhs(i,j,k,:) = rhs(i,j,k,:) + face_flux/dy
          rhs(i,jp,k,:) = rhs(i,jp,k,:) - face_flux/dy

          call compute_transport_face_flux( &
            0.5_dp*(density(i,j,k)+density(i,j,kp)), &
            0.5_dp*(velocity(i,j,k,:)+velocity(i,j,kp,:)), &
            0.5_dp*(temperature(i,j,k)+temperature(i,j,kp)), &
            0.5_dp*(mass_fractions(i,j,k,:)+ &
              mass_fractions(i,j,kp,:)), &
            0.5_dp*(velocity_gradient(i,j,k,:,:)+ &
              velocity_gradient(i,j,kp,:,:)), &
            0.5_dp*(temperature_gradient(i,j,k,:)+ &
              temperature_gradient(i,j,kp,:)), &
            0.5_dp*(fraction_gradient(i,j,k,:,:)+ &
              fraction_gradient(i,j,kp,:,:)),layout,3,face_flux)
          rhs(i,j,k,:) = rhs(i,j,k,:) + face_flux/dz
          rhs(i,j,kp,:) = rhs(i,j,kp,:) - face_flux/dz
        end do
      end do
    end do
    deallocate(density,temperature,velocity,mass_fractions, &
      velocity_gradient,temperature_gradient,fraction_gradient)
  end subroutine compute_mc_transport_rhs

  subroutine compute_primitive_fields( &
      q,layout,config,density,velocity,temperature,mass_fractions)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    real(dp), intent(out) :: density(:,:,:), velocity(:,:,:,:)
    real(dp), intent(out) :: temperature(:,:,:)
    real(dp), intent(out) :: mass_fractions(:,:,:,:)
    integer :: i, j, k

    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          density(i,j,k) = mc_mixture_density(q(i,j,k,:),layout)
          velocity(i,j,k,:) = q(i,j,k,layout%momentum)/density(i,j,k)
          temperature(i,j,k) = &
            mc_temperature(q(i,j,k,:),layout,config%gamma)
          mass_fractions(i,j,k,:) = &
            q(i,j,k,layout%first_species:layout%last_species)/density(i,j,k)
        end do
      end do
    end do
  end subroutine compute_primitive_fields

  subroutine compute_periodic_gradients( &
      velocity,temperature,mass_fractions,dx,dy,dz, &
      velocity_gradient,temperature_gradient,fraction_gradient)
    real(dp), intent(in) :: velocity(:,:,:,:), temperature(:,:,:)
    real(dp), intent(in) :: mass_fractions(:,:,:,:), dx, dy, dz
    real(dp), intent(out) :: velocity_gradient(:,:,:,:,:)
    real(dp), intent(out) :: temperature_gradient(:,:,:,:)
    real(dp), intent(out) :: fraction_gradient(:,:,:,:,:)
    integer :: i, j, k, im, ip, jm, jp, km, kp

    do k = 1, size(temperature,3)
      km = merge(size(temperature,3),k-1,k == 1)
      kp = merge(1,k+1,k == size(temperature,3))
      do j = 1, size(temperature,2)
        jm = merge(size(temperature,2),j-1,j == 1)
        jp = merge(1,j+1,j == size(temperature,2))
        do i = 1, size(temperature,1)
          im = merge(size(temperature,1),i-1,i == 1)
          ip = merge(1,i+1,i == size(temperature,1))
          velocity_gradient(i,j,k,:,1) = &
            (velocity(ip,j,k,:)-velocity(im,j,k,:))/(2.0_dp*dx)
          velocity_gradient(i,j,k,:,2) = &
            (velocity(i,jp,k,:)-velocity(i,jm,k,:))/(2.0_dp*dy)
          velocity_gradient(i,j,k,:,3) = &
            (velocity(i,j,kp,:)-velocity(i,j,km,:))/(2.0_dp*dz)
          temperature_gradient(i,j,k,1) = &
            (temperature(ip,j,k)-temperature(im,j,k))/(2.0_dp*dx)
          temperature_gradient(i,j,k,2) = &
            (temperature(i,jp,k)-temperature(i,jm,k))/(2.0_dp*dy)
          temperature_gradient(i,j,k,3) = &
            (temperature(i,j,kp)-temperature(i,j,km))/(2.0_dp*dz)
          fraction_gradient(i,j,k,:,1) = &
            (mass_fractions(ip,j,k,:)-mass_fractions(im,j,k,:))/ &
            (2.0_dp*dx)
          fraction_gradient(i,j,k,:,2) = &
            (mass_fractions(i,jp,k,:)-mass_fractions(i,jm,k,:))/ &
            (2.0_dp*dy)
          fraction_gradient(i,j,k,:,3) = &
            (mass_fractions(i,j,kp,:)-mass_fractions(i,j,km,:))/ &
            (2.0_dp*dz)
        end do
      end do
    end do
  end subroutine compute_periodic_gradients

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

  subroutine compute_mc_navier_stokes_rhs(q,rhs,layout,config)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    real(dp), allocatable :: transport_rhs(:,:,:,:)

    allocate(transport_rhs,mold=q)
    call compute_mc_euler_rhs(q,rhs,layout,config)
    call compute_mc_transport_rhs(q,transport_rhs,layout,config)
    rhs = rhs+transport_rhs
    deallocate(transport_rhs)
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
      q,q0,rhs,dt,layout,config)
    real(dp), intent(inout) :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config

    if (any(shape(q0) /= shape(q)) .or. any(shape(rhs) /= shape(q))) then
      error stop 'multicomponent Navier-Stokes work arrays do not match state'
    end if
    q0 = q
    call compute_mc_navier_stokes_rhs(q,rhs,layout,config)
    q = q0+dt*rhs
    call compute_mc_navier_stokes_rhs(q,rhs,layout,config)
    q = 0.75_dp*q0+0.25_dp*(q+dt*rhs)
    call compute_mc_navier_stokes_rhs(q,rhs,layout,config)
    q = (1.0_dp/3.0_dp)*q0+(2.0_dp/3.0_dp)*(q+dt*rhs)
    call validate_mc_euler_state(q,layout,config)
  end subroutine advance_mc_navier_stokes_ssprk3

end module mod_mc_viscous_flux
