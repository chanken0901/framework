module mod_mc_euler_solver
  use mod_mc_euler_field, only : mc_primitive_workspace
  use mod_precision, only : dp
  use mod_mc_geometry, only: mc_cell_center
  use mod_mc_config, only : mc_config
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    compute_mc_euler_totals, compute_mc_euler_minima
  use mod_mc_euler_flux, only : compute_mc_euler_timestep, &
    advance_mc_euler_ssprk3
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, mc_pressure, &
    mc_temperature, mc_thermodynamics_provider_name
  implicit none
  private

  public :: run_mc_euler
  public :: write_mc_euler_csv

contains

  subroutine run_mc_euler(config,layout,euler)
    type(mc_config), intent(in) :: config
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: euler
    type(mc_primitive_workspace) :: workspace
    real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), allocatable :: initial_totals(:), final_totals(:)
    real(dp) :: dt, time, conservation_error, scale
    real(dp) :: minimum_species, minimum_density, minimum_pressure
    real(dp) :: minimum_temperature
    integer :: step, variable

    allocate(q(euler%nx,euler%ny,euler%nz,layout%nvariables))
    allocate(q0,mold=q)
    allocate(rhs,mold=q)
    allocate(initial_totals(layout%nvariables),final_totals(layout%nvariables))
    call initialize_mc_euler_state(q,layout,euler)
    call compute_mc_euler_totals(q,euler,initial_totals)

    time = 0.0_dp
    if (trim(config%simulation_mode) == 'thermally_perfect_euler') then
      write(*,'(A)') '--- stage-3 thermally-perfect multicomponent Euler ---'
    else
      write(*,'(A)') '--- stage-2 multicomponent inviscid Euler ---'
    end if
    write(*,'(A,3(I0,1X))') 'grid = ', euler%nx, euler%ny, euler%nz
    if (mc_thermodynamics_provider_name == 'calorically_perfect') then
      write(*,'(A,F8.4)') 'gamma = ', euler%gamma
    end if
    write(*,'(A,I0)') 'nsteps = ', euler%nsteps
    do step = 1, euler%nsteps
      dt = compute_mc_euler_timestep(q,layout,euler)
      call advance_mc_euler_ssprk3(q,q0,rhs,dt,layout,euler,workspace)
      time = time + dt
    end do

    call compute_mc_euler_totals(q,euler,final_totals)
    conservation_error = 0.0_dp
    do variable = 1, layout%nvariables
      scale = max(abs(initial_totals(variable)),1.0_dp)
      conservation_error = max(conservation_error, &
        abs(final_totals(variable)-initial_totals(variable))/scale)
    end do
    call compute_mc_euler_minima( &
      q,layout,euler,minimum_species,minimum_density,minimum_pressure, &
      minimum_temperature)
    if (conservation_error > 1.0e-10_dp) then
      error stop 'multicomponent Euler conservation check failed'
    end if
    if (minimum_species < -1.0e-12_dp .or. &
        minimum_density <= 0.0_dp .or. minimum_pressure <= 0.0_dp) then
      error stop 'multicomponent Euler positivity check failed'
    end if
    if (euler%write_final) then
      call write_mc_euler_csv(trim(euler%output_file),q,config,layout,euler)
    end if

    write(*,'(A,ES12.4)') 'final time = ', time
    write(*,'(A,ES12.4)') 'maximum relative conservation error = ', &
      conservation_error
    write(*,'(A,ES12.4)') 'minimum species partial density = ', &
      minimum_species
    write(*,'(A,ES12.4)') 'minimum mixture density = ', minimum_density
    write(*,'(A,ES12.4)') 'minimum pressure = ', minimum_pressure
    if (mc_thermodynamics_provider_name == 'thermally_perfect') then
      write(*,'(A,ES12.4)') 'minimum temperature = ', minimum_temperature
    end if
    write(*,'(A)') 'Multicomponent Euler calculation completed successfully.'

    deallocate(q,q0,rhs,initial_totals,final_totals)
  end subroutine run_mc_euler

  subroutine write_mc_euler_csv(path,q,config,layout,euler)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_config), intent(in) :: config
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: euler
    integer :: unit, ios, i, j, k, species, variable
    real(dp) :: x, y, z, dx, dy, dz, density, velocity(3), pressure
    real(dp) :: temperature, position(3)
    character(len=512) :: message

    open(newunit=unit,file=trim(path),status='replace',action='write', &
      iostat=ios,iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open multicomponent Euler output: ', &
        trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open multicomponent Euler output'
    end if
    write(unit,'(A)',advance='no') 'x,y,z'
    do species = 1, layout%nspecies
      write(unit,'(A,A)',advance='no') ',rho_', &
        trim(config%species_names(species))
    end do
    if (mc_thermodynamics_provider_name == 'thermally_perfect') then
      write(unit,'(A)') ',rho,u,v,w,p,T,rhoE'
    else
      write(unit,'(A)') ',rho,u,v,w,p,rhoE'
    end if

    dx = (euler%x_max-euler%x_min)/real(euler%nx,dp)
    dy = (euler%y_max-euler%y_min)/real(euler%ny,dp)
    dz = (euler%z_max-euler%z_min)/real(euler%nz,dp)
    do k = 1, euler%nz
      z = euler%z_min + (real(k,dp)-0.5_dp)*dz
      do j = 1, euler%ny
        y = euler%y_min + (real(j,dp)-0.5_dp)*dy
        do i = 1, euler%nx
          position=mc_cell_center(euler,i,j,k)
          x=position(1)
          y=position(2)
          z=position(3)
          density = mc_mixture_density(q(i,j,k,:),layout)
          velocity = q(i,j,k,layout%momentum)/density
          pressure = mc_pressure(q(i,j,k,:),layout,euler%gamma)
          if (mc_thermodynamics_provider_name == 'thermally_perfect') then
            temperature = mc_temperature(q(i,j,k,:),layout,euler%gamma)
          end if
          write(unit,'(ES24.16)',advance='no') x
          write(unit,'(A,ES24.16)',advance='no') ',', y
          write(unit,'(A,ES24.16)',advance='no') ',', z
          do species = 1, layout%nspecies
            variable = layout%first_species + species - 1
            write(unit,'(A,ES24.16)',advance='no') ',', q(i,j,k,variable)
          end do
          write(unit,'(A,ES24.16)',advance='no') ',', density
          write(unit,'(A,ES24.16)',advance='no') ',', velocity(1)
          write(unit,'(A,ES24.16)',advance='no') ',', velocity(2)
          write(unit,'(A,ES24.16)',advance='no') ',', velocity(3)
          write(unit,'(A,ES24.16)',advance='no') ',', pressure
          if (mc_thermodynamics_provider_name == 'thermally_perfect') then
            write(unit,'(A,ES24.16)',advance='no') ',', temperature
          end if
          write(unit,'(A,ES24.16)') ',', q(i,j,k,layout%total_energy)
        end do
      end do
    end do
    close(unit)
    write(*,'(A,A)') 'multicomponent Euler output = ', trim(path)
  end subroutine write_mc_euler_csv

end module mod_mc_euler_solver
