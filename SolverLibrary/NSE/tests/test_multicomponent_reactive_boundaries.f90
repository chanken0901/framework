program test_multicomponent_reactive_boundaries
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, read_mc_config
  use mod_mc_state_layout, only : mc_state_layout, initialize_mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, read_mc_euler_config, &
    validate_mc_euler_config, mc_face_x_min, mc_face_x_max
  use mod_mc_euler_field, only : initialize_mc_euler_state, &
    set_mc_euler_conservative_state, validate_mc_euler_state
  use mod_mc_boundary, only : mc_boundary_state
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_rhs
  use mod_mc_reactive_config, only : mc_reactive_config, &
    read_mc_reactive_config
  use mod_mc_reactive_solver, only : compute_mc_reactive_timestep, &
    advance_mc_reactive_strang
  use mod_mc_thermodynamics_provider, only : configure_mc_thermodynamics
  use mod_mc_transport_provider, only : configure_mc_transport
  use mod_mc_chemistry_provider, only : configure_mc_chemistry
  implicit none

  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: numerics, reflected
  type(mc_reactive_config) :: reactive
  character(len=512) :: input_path
  real(dp), allocatable :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
  real(dp) :: ghost(7), expected(7), dt, density
  integer :: step

  if (command_argument_count() < 1) then
    error stop 'stage-7 test requires its input file path'
  end if
  call get_command_argument(1,input_path)
  call read_mc_config(trim(input_path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  call configure_mc_thermodynamics( &
    trim(input_path),model%nspecies,model%species_names(1:model%nspecies))
  call configure_mc_transport( &
    trim(input_path),model%nspecies,model%species_names(1:model%nspecies))
  call configure_mc_chemistry( &
    trim(input_path),model%nspecies,model%species_names(1:model%nspecies))
  call read_mc_euler_config(trim(input_path),model%nspecies,numerics)
  call read_mc_reactive_config(trim(input_path),reactive)
  numerics%write_final = .false.
  reactive%write_snapshots = .false.
  reactive%write_history = .false.

  allocate(q(numerics%nx,numerics%ny,numerics%nz,layout%nvariables))
  allocate(q0,mold=q)
  allocate(rhs,mold=q)
  call initialize_mc_euler_state(q,layout,numerics)
  call assert_true(all(ieee_is_finite(q)), &
    'shock-tube initial state must be finite')

  call mc_boundary_state( &
    q(1,1,1,:),ghost,layout,numerics,mc_face_x_min)
  call set_mc_euler_conservative_state( &
    expected,layout,numerics%gamma, &
    numerics%boundary_reference_densities(mc_face_x_min), &
    numerics%boundary_reference_velocities(:,mc_face_x_min), &
    numerics%boundary_reference_pressures(mc_face_x_min), &
    numerics%boundary_reference_mass_fractions( &
      mc_face_x_min,1:layout%nspecies))
  call assert_close(maxval(abs(ghost-expected)),0.0_dp,1.0e-12_dp, &
    'Dirichlet boundary must reproduce its reference state')

  call mc_boundary_state( &
    q(numerics%nx,1,1,:),ghost,layout,numerics,mc_face_x_max)
  call assert_close( &
    maxval(abs(ghost-q(numerics%nx,1,1,:)))/ &
    max(maxval(abs(q(numerics%nx,1,1,:))),1.0_dp), &
    0.0_dp,1.0e-12_dp, &
    'uniform non-reflecting boundary must preserve the far field')

  reflected = numerics
  reflected%boundary_face_types(mc_face_x_min) = 'reflective'
  reflected%boundary_face_types(mc_face_x_max) = 'reflective'
  reflected%boundary_condition = 'face_specific'
  call validate_mc_euler_config(reflected,model%nspecies)
  density = sum(q(1,1,1,layout%first_species:layout%last_species))
  q(1,1,1,layout%momentum(1)) = 0.2_dp*density
  call mc_boundary_state( &
    q(1,1,1,:),ghost,layout,reflected,mc_face_x_min)
  call assert_close(ghost(layout%momentum(1)), &
    -q(1,1,1,layout%momentum(1)),1.0e-12_dp, &
    'reflective boundary must reverse normal momentum')
  call assert_close(ghost(layout%total_energy), &
    q(1,1,1,layout%total_energy),1.0e-12_dp, &
    'reflective boundary must preserve total energy')

  call initialize_mc_euler_state(q,layout,numerics)
  call compute_mc_navier_stokes_rhs(q,rhs,layout,numerics)
  call assert_true(all(ieee_is_finite(rhs)), &
    'face-specific Navier-Stokes RHS must be finite')
  do step = 1, 3
    dt = compute_mc_reactive_timestep(q,layout,numerics,reactive)
    call advance_mc_reactive_strang( &
      q,q0,rhs,dt,layout,numerics,reactive)
  end do
  call validate_mc_euler_state(q,layout,numerics)

  write(*,'(A)') &
    'Reactive multicomponent boundary and shock-tube tests passed.'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(message)
      error stop 'reactive multicomponent boundary test failure'
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected_value,tolerance,message)
    real(dp), intent(in) :: actual, expected_value, tolerance
    character(len=*), intent(in) :: message
    real(dp) :: scale

    scale = max(abs(expected_value),1.0_dp)
    call assert_true(abs(actual-expected_value) <= tolerance*scale,message)
  end subroutine assert_close

end program test_multicomponent_reactive_boundaries
