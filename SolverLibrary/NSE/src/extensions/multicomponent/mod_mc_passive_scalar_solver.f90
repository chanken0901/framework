module mod_mc_passive_scalar_solver
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_passive_scalar_config, only : mc_passive_scalar_config, &
    mc_passive_scalar_timestep
  use mod_mc_passive_scalar_field, only : &
    initialize_mc_passive_scalar_state, compute_mc_species_masses, &
    compute_mc_species_sum_error, compute_mc_tracer_bounds
  use mod_mc_passive_scalar_advection, only : &
    advance_mc_passive_scalar_ssprk3
  implicit none
  private

  public :: run_mc_passive_scalar

contains

  subroutine run_mc_passive_scalar(config, layout, passive)
    type(mc_config), intent(in) :: config
    type(mc_state_layout), intent(in) :: layout
    type(mc_passive_scalar_config), intent(in) :: passive
    real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), allocatable :: initial_mass(:), final_mass(:)
    real(dp) :: dt, time, mass_drift, sum_error
    real(dp) :: tracer_minimum, tracer_maximum, scale
    integer :: step, species

    allocate(q(passive%nx,passive%ny,passive%nz,layout%nvariables))
    allocate(q0, mold=q)
    allocate(rhs, mold=q)
    allocate(initial_mass(layout%nspecies), final_mass(layout%nspecies))
    call initialize_mc_passive_scalar_state(q, layout, passive)
    call compute_mc_species_masses(q, layout, passive, initial_mass)
    dt = mc_passive_scalar_timestep(passive)

    write(*,'(A)') '--- stage-1 passive-scalar advection ---'
    write(*,'(A,3(I0,1X))') 'grid = ', passive%nx, passive%ny, passive%nz
    write(*,'(A,3(ES12.4,1X))') 'velocity = ', passive%velocity
    write(*,'(A,ES12.4)') 'dt = ', dt
    write(*,'(A,I0)') 'nsteps = ', passive%nsteps
    do step = 1, passive%nsteps
      call advance_mc_passive_scalar_ssprk3(q, q0, rhs, dt, layout, passive)
    end do
    time = real(passive%nsteps,dp) * dt

    call compute_mc_species_masses(q, layout, passive, final_mass)
    mass_drift = 0.0_dp
    do species = 1, layout%nspecies
      scale = max(abs(initial_mass(species)), 1.0_dp)
      mass_drift = max(mass_drift, &
        abs(final_mass(species)-initial_mass(species))/scale)
    end do
    sum_error = compute_mc_species_sum_error(q, layout)
    call compute_mc_tracer_bounds(q, layout, tracer_minimum, tracer_maximum)

    if (mass_drift > 1.0e-11_dp) then
      error stop 'passive-scalar species mass conservation failed'
    end if
    if (sum_error > 1.0e-11_dp) then
      error stop 'passive-scalar species sum constraint failed'
    end if
    if (tracer_minimum < -1.0e-12_dp .or. &
        tracer_maximum > 1.0_dp+1.0e-12_dp) then
      error stop 'passive-scalar boundedness check failed'
    end if
    if (passive%write_final) then
      call write_mc_passive_scalar_csv( &
        trim(passive%output_file), q, config, layout, passive)
    end if

    write(*,'(A,ES12.4)') 'final time = ', time
    write(*,'(A,ES12.4)') 'maximum relative species-mass drift = ', mass_drift
    write(*,'(A,ES12.4)') 'maximum species-sum error = ', sum_error
    write(*,'(A,2(ES12.4,1X))') 'tracer bounds = ', &
      tracer_minimum, tracer_maximum
    write(*,'(A)') 'Passive-scalar calculation completed successfully.'

    deallocate(q, q0, rhs, initial_mass, final_mass)
  end subroutine run_mc_passive_scalar

  subroutine write_mc_passive_scalar_csv(path, q, config, layout, passive)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_config), intent(in) :: config
    type(mc_state_layout), intent(in) :: layout
    type(mc_passive_scalar_config), intent(in) :: passive
    integer :: unit, ios, i, j, k, species, variable
    real(dp) :: x, y, z, dx, dy, dz
    character(len=512) :: message

    open(newunit=unit, file=trim(path), status='replace', action='write', &
      iostat=ios, iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open passive-scalar output: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open passive-scalar output'
    end if
    write(unit,'(A)',advance='no') 'x,y,z'
    do species = 1, layout%nspecies
      write(unit,'(A,A)',advance='no') ',rho_', &
        trim(config%species_names(species))
    end do
    write(unit,'()')

    dx = (passive%x_max-passive%x_min) / real(passive%nx,dp)
    dy = (passive%y_max-passive%y_min) / real(passive%ny,dp)
    dz = (passive%z_max-passive%z_min) / real(passive%nz,dp)
    do k = 1, passive%nz
      z = passive%z_min + (real(k,dp)-0.5_dp)*dz
      do j = 1, passive%ny
        y = passive%y_min + (real(j,dp)-0.5_dp)*dy
        do i = 1, passive%nx
          x = passive%x_min + (real(i,dp)-0.5_dp)*dx
          write(unit,'(ES24.16)',advance='no') x
          write(unit,'(A,ES24.16)',advance='no') ',', y
          write(unit,'(A,ES24.16)',advance='no') ',', z
          do species = 1, layout%nspecies
            variable = layout%first_species + species - 1
            write(unit,'(A,ES24.16)',advance='no') ',', q(i,j,k,variable)
          end do
          write(unit,'()')
        end do
      end do
    end do
    close(unit)
    write(*,'(A,A)') 'passive-scalar output = ', trim(path)
  end subroutine write_mc_passive_scalar_csv

end module mod_mc_passive_scalar_solver
