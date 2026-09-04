module mod_mc_reactor_solver
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_reactor_config, only : mc_reactor_config
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_mixture_gas_constant, mc_pressure, mc_temperature, &
    mc_total_energy_from_primitive
  use mod_mc_chemistry_provider, only : compute_mc_chemistry_timestep
  use mod_mc_chemistry_integrator, only : validate_mc_chemistry_state, &
    advance_mc_chemistry_ssprk3
  implicit none
  private

  public :: initialize_mc_reactor_state
  public :: validate_mc_reactor_state
  public :: advance_mc_reactor_ssprk3
  public :: run_mc_reactor

contains

  subroutine initialize_mc_reactor_state(state,layout,config)
    real(dp), intent(out) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_reactor_config), intent(in) :: config
    real(dp) :: gas_constant, pressure
    real(dp) :: velocity(3)

    if (size(state) /= layout%nvariables) then
      error stop 'homogeneous-reactor state does not match layout'
    end if
    velocity = 0.0_dp
    gas_constant = mc_mixture_gas_constant( &
      config%initial_mass_fractions(1:layout%nspecies),layout)
    pressure = config%initial_density*gas_constant* &
      config%initial_temperature
    state = 0.0_dp
    state(layout%first_species:layout%last_species) = &
      config%initial_density* &
      config%initial_mass_fractions(1:layout%nspecies)
    state(layout%total_energy) = mc_total_energy_from_primitive( &
      layout,1.4_dp,config%initial_density,velocity,pressure, &
      config%initial_mass_fractions(1:layout%nspecies))
    call validate_mc_reactor_state(state,layout)
  end subroutine initialize_mc_reactor_state

  subroutine validate_mc_reactor_state(state,layout)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout

    call validate_mc_chemistry_state(state,layout,1.4_dp)
  end subroutine validate_mc_reactor_state

  subroutine advance_mc_reactor_ssprk3(state,dt,layout)
    real(dp), intent(inout) :: state(:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout

    call advance_mc_chemistry_ssprk3(state,dt,layout,1.4_dp)
  end subroutine advance_mc_reactor_ssprk3

  subroutine write_history_header(unit,model)
    integer, intent(in) :: unit
    type(mc_config), intent(in) :: model
    integer :: species

    write(unit,'(A)',advance='no') 'step,time,dt,density,temperature,pressure'
    do species = 1, model%nspecies
      write(unit,'(A)',advance='no') ',Y_'// &
        trim(model%species_names(species))
    end do
    write(unit,'(A)') ''
  end subroutine write_history_header

  subroutine write_history_row(unit,step,time,dt,state,layout)
    integer, intent(in) :: unit, step
    real(dp), intent(in) :: time, dt, state(:)
    type(mc_state_layout), intent(in) :: layout
    integer :: species
    real(dp) :: density, temperature, pressure

    density = mc_mixture_density(state,layout)
    temperature = mc_temperature(state,layout,1.4_dp)
    pressure = mc_pressure(state,layout,1.4_dp)
    write(unit,'(I0,5(",",ES24.16))',advance='no') &
      step,time,dt,density,temperature,pressure
    do species = 1, layout%nspecies
      write(unit,'(",",ES24.16)',advance='no') state(species)/density
    end do
    write(unit,'(A)') ''
  end subroutine write_history_row

  subroutine run_mc_reactor(model,layout,config)
    type(mc_config), intent(in) :: model
    type(mc_state_layout), intent(in) :: layout
    type(mc_reactor_config), intent(in) :: config
    real(dp) :: state(layout%nvariables)
    real(dp) :: dt, stable_dt, time, initial_density, initial_energy
    real(dp) :: density_scale, energy_scale
    integer :: unit, ios, step
    character(len=512) :: message

    call initialize_mc_reactor_state(state,layout,config)
    initial_density = mc_mixture_density(state,layout)
    initial_energy = state(layout%total_energy)
    density_scale = max(abs(initial_density),1.0_dp)
    energy_scale = max(abs(initial_energy),1.0_dp)
    time = 0.0_dp

    write(*,'(A)') '--- stage-5 homogeneous finite-rate reactor ---'
    write(*,'(A,ES14.6)') 'initial temperature = ', &
      mc_temperature(state,layout,1.4_dp)
    if (config%write_history) then
      open(newunit=unit,file=trim(config%output_file),status='replace', &
        action='write',iostat=ios,iomsg=message)
      if (ios /= 0) then
        write(*,'(A,A)') 'ERROR: cannot open reactor history: ', &
          trim(config%output_file)
        write(*,'(A,A)') 'ERROR: ', trim(message)
        error stop 'failed to open reactor history'
      end if
      call write_history_header(unit,model)
      call write_history_row(unit,0,time,0.0_dp,state,layout)
    end if

    do step = 1, config%nsteps
      stable_dt = compute_mc_chemistry_timestep( &
        state,layout,1.4_dp,config%chemistry_cfl,config%maximum_dt)
      if (config%dt > 0.0_dp) then
        if (config%dt > stable_dt*(1.0_dp+1.0e-12_dp)) then
          error stop 'fixed homogeneous-reactor dt violates chemistry CFL'
        end if
        dt = config%dt
      else
        dt = stable_dt
      end if
      call advance_mc_reactor_ssprk3(state,dt,layout)
      time = time + dt
      if (abs(mc_mixture_density(state,layout)-initial_density) > &
          1.0e-11_dp*density_scale) then
        error stop 'homogeneous reaction failed total-mass conservation'
      end if
      if (abs(state(layout%total_energy)-initial_energy) > &
          1.0e-13_dp*energy_scale) then
        error stop 'homogeneous reaction failed total-energy conservation'
      end if
      if (config%write_history .and. &
          (mod(step,config%output_every) == 0 .or. &
          step == config%nsteps)) then
        call write_history_row(unit,step,time,dt,state,layout)
      end if
    end do
    if (config%write_history) close(unit)

    write(*,'(A,I0,A,ES14.6,A,ES14.6)') &
      'Homogeneous reactor calculation completed successfully: steps=', &
      config%nsteps, ', time=', time, ', temperature=', &
      mc_temperature(state,layout,1.4_dp)
  end subroutine run_mc_reactor

end module mod_mc_reactor_solver
