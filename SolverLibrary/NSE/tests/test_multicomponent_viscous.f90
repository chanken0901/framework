program test_multicomponent_viscous
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, read_mc_config
  use mod_mc_state_layout, only : mc_state_layout, &
    initialize_mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, read_mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    set_mc_euler_conservative_state, validate_mc_euler_state, &
    compute_mc_euler_totals
  use mod_mc_viscous_flux, only : compute_mc_transport_rhs, &
    compute_mc_navier_stokes_timestep, advance_mc_navier_stokes_ssprk3
  use mod_mc_thermodynamics_provider, only : &
    configure_mc_thermodynamics, mc_mixture_gas_constant
  use mod_mc_transport_provider, only : configure_mc_transport
  implicit none

  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: config
  character(len=512) :: input_path
  real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
  real(dp), allocatable :: initial_totals(:), final_totals(:)
  real(dp), allocatable :: state(:)
  real(dp) :: fractions(2), velocity(3), pressure, gas_constant
  real(dp) :: x, dx, phase, damping, total, scale, dt
  real(dp), parameter :: pi = acos(-1.0_dp)
  integer :: i, j, k, variable, step

  if (command_argument_count() < 1) then
    error stop 'stage-4 test requires its input file path'
  end if
  call get_command_argument(1,input_path)
  call read_mc_config(trim(input_path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  call configure_mc_thermodynamics( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call configure_mc_transport( &
    trim(input_path),model%nspecies, &
    model%species_names(1:model%nspecies))
  call read_mc_euler_config(trim(input_path),model%nspecies,config)
  config%write_final = .false.

  allocate(q(config%nx,config%ny,config%nz,layout%nvariables))
  allocate(q0,mold=q)
  allocate(rhs,mold=q)
  allocate(state(layout%nvariables))
  allocate(initial_totals(layout%nvariables),final_totals(layout%nvariables))
  call initialize_mc_euler_state(q,layout,config)
  call compute_mc_transport_rhs(q,rhs,layout,config)

  call assert_true(maxval(abs(sum( &
    rhs(:,:,:,layout%first_species:layout%last_species),dim=4))) < 1.0e-12_dp, &
    'corrected species diffusion fluxes must sum to zero in every cell')
  do variable = 1, layout%nvariables
    total = sum(rhs(:,:,:,variable))
    call assert_true(abs(total) < &
      1.0e-11_dp*max(sum(abs(rhs(:,:,:,variable))),1.0_dp), &
      'periodic transport RHS must conserve every state variable')
  end do
  damping = 0.0_dp
  dx = (config%x_max-config%x_min)/real(config%nx,dp)
  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        x = config%x_min+(real(i,dp)-0.5_dp)*dx
        phase = 2.0_dp*pi*real(config%wave_wavenumber,dp)* &
          (x-config%x_min)/(config%x_max-config%x_min)
        damping = damping + sin(phase)* &
          rhs(i,j,k,layout%first_species)
      end do
    end do
  end do
  call assert_true(damping < 0.0_dp, &
    'species diffusion must damp the periodic composition wave')

  fractions = [0.5_dp,0.5_dp]
  velocity = 0.0_dp
  gas_constant = mc_mixture_gas_constant(fractions,layout)
  pressure = config%wave_density*gas_constant*config%wave_temperature
  call set_mc_euler_conservative_state( &
    state,layout,config%gamma,config%wave_density,velocity,pressure,fractions)
  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        q(i,j,k,:) = state
      end do
    end do
  end do
  call compute_mc_transport_rhs(q,rhs,layout,config)
  call assert_true(maxval(abs(rhs)) < 1.0e-12_dp, &
    'uniform state must have zero viscous and diffusive RHS')

  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        x = config%x_min+(real(i,dp)-0.5_dp)*dx
        phase = 2.0_dp*pi*(x-config%x_min)/ &
          (config%x_max-config%x_min)
        velocity = [10.0_dp*sin(phase),0.0_dp,0.0_dp]
        call set_mc_euler_conservative_state( &
          q(i,j,k,:),layout,config%gamma,config%wave_density,velocity, &
          pressure,fractions)
      end do
    end do
  end do
  call compute_mc_transport_rhs(q,rhs,layout,config)
  damping = 0.0_dp
  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        damping = damping + &
          q(i,j,k,layout%momentum(1))/config%wave_density * &
          rhs(i,j,k,layout%momentum(1))
      end do
    end do
  end do
  call assert_true(damping < 0.0_dp, &
    'viscous stress must damp the periodic velocity wave')

  velocity = 0.0_dp
  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        x = config%x_min+(real(i,dp)-0.5_dp)*dx
        phase = 2.0_dp*pi*(x-config%x_min)/ &
          (config%x_max-config%x_min)
        pressure = config%wave_density*gas_constant* &
          (config%wave_temperature+20.0_dp*sin(phase))
        call set_mc_euler_conservative_state( &
          q(i,j,k,:),layout,config%gamma,config%wave_density,velocity, &
          pressure,fractions)
      end do
    end do
  end do
  call compute_mc_transport_rhs(q,rhs,layout,config)
  damping = 0.0_dp
  do k = 1, config%nz
    do j = 1, config%ny
      do i = 1, config%nx
        x = config%x_min+(real(i,dp)-0.5_dp)*dx
        phase = 2.0_dp*pi*(x-config%x_min)/ &
          (config%x_max-config%x_min)
        damping = damping + &
          sin(phase)*rhs(i,j,k,layout%total_energy)
      end do
    end do
  end do
  call assert_true(damping < 0.0_dp, &
    'Fourier heat conduction must damp the periodic temperature wave')

  call initialize_mc_euler_state(q,layout,config)
  call compute_mc_euler_totals(q,config,initial_totals)
  do step = 1, config%nsteps
    dt = compute_mc_navier_stokes_timestep(q,layout,config)
    call assert_true(dt > 0.0_dp, &
      'stage-4 stable timestep must be positive')
    call advance_mc_navier_stokes_ssprk3(q,q0,rhs,dt,layout,config)
  end do
  call validate_mc_euler_state(q,layout,config)
  call compute_mc_euler_totals(q,config,final_totals)
  do variable = 1, layout%nvariables
    scale = max(abs(initial_totals(variable)),1.0_dp)
    call assert_true(abs(final_totals(variable)-initial_totals(variable)) < &
      2.0e-10_dp*scale, &
      'periodic stage-4 update must conserve every state variable')
  end do

  write(*,'(A)') 'Multicomponent viscous transport tests passed.'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'multicomponent viscous transport test failure'
    end if
  end subroutine assert_true

end program test_multicomponent_viscous
