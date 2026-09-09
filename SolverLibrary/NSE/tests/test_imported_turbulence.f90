program test_imported_turbulence
  use, intrinsic :: iso_fortran_env, only : int32
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config, &
    update_derived_config
  use mod_model_config, only : nse_config, init_nse_config, &
    nse_face_x_min, nse_face_x_max
  use mod_init_imported_turbulence, only : &
    initialize_imported_turbulence, imported_turbulence_weight
  use mod_init_shock_turbulence, only : initialize_shock_turbulence
  use mod_nse_initial_conditions, only : initialize_nse_state
  implicit none

  character(len=*), parameter :: source_file = &
    'test_imported_turbulence_source.slf'
  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:)
  real(dp) :: expected(5), background_test(5)
  integer :: unit, ios, first_i, last_i

  call write_source_slf(source_file)
  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 8
  sim%ny = 2
  sim%nz = 2
  sim%nghost = 3
  sim%x_min = 0.0_dp
  sim%x_max = 8.0_dp
  sim%y_min = 0.0_dp
  sim%y_max = 2.0_dp
  sim%z_min = 0.0_dp
  sim%z_max = 2.0_dp
  sim%rank = 0
  call update_derived_config(sim)

  nse%imported_turbulence_file = source_file
  nse%imported_turbulence_mode = 'embed'
  nse%imported_turbulence_x_start = 2.0_dp
  nse%imported_turbulence_blend_cells = 0
  nse%imported_turbulence_velocity_offset_x = 1.0_dp
  nse%imported_turbulence_background_rho = 1.0_dp
  nse%imported_turbulence_background_u = 10.0_dp
  nse%imported_turbulence_background_p = 2.0_dp

  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost, 5))
  q = -999.0_dp
  call initialize_imported_turbulence(q, sim, nse, 1, sim%ny, 1, sim%nz)

  call primitive_to_conserved(1.0_dp, 10.0_dp, 0.0_dp, 0.0_dp, &
    2.0_dp, nse%gamma, expected)
  call assert_vector_close(q(2,1,1,:), expected, 'embed background left')
  call assert_vector_close(q(7,2,2,:), expected, 'embed background right')
  call primitive_to_conserved(2.0_dp, 2.0_dp, 0.1_dp, 0.2_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(3,1,1,:), expected, 'embed first source cell')
  call primitive_to_conserved(2.0_dp, 5.0_dp, 0.2_dp, 0.4_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(6,2,2,:), expected, 'embed last source cell')

  nse%imported_turbulence_mode = 'tile'
  nse%imported_turbulence_x_start = 0.0_dp
  nse%imported_turbulence_background_u = 0.0_dp
  q = -999.0_dp
  call initialize_imported_turbulence(q, sim, nse, 1, sim%ny, 1, sim%nz)
  call primitive_to_conserved(2.0_dp, 2.0_dp, 0.1_dp, 0.2_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(1,1,1,:), expected, 'tile first source cell')
  call assert_vector_close(q(5,1,1,:), expected, 'tile repeated source cell')

  deallocate(q)
  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    2-sim%nghost:2+sim%nghost, 2-sim%nghost:2+sim%nghost, 5))
  q = -999.0_dp
  call initialize_imported_turbulence(q, sim, nse, 2, 2, 2, 2)
  call primitive_to_conserved(2.0_dp, 2.0_dp, 0.2_dp, 0.4_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(1,2,2,:), expected, &
    'decomposed y-z source selection')

  ! A six-cell window of a four-cell periodic source: 1,2,3,4,1,2.
  nse%imported_turbulence_mode = 'periodic_embed'
  nse%imported_turbulence_x_start = 1.0_dp
  nse%imported_turbulence_x_length = 6.0_dp
  call initialize_imported_turbulence(q,sim,nse,2,2,2,2,first_i,last_i)
  call assert_true(first_i == 2 .and. last_i == 7,'periodic window bounds')
  call primitive_to_conserved(1.0_dp,0.0_dp,0.0_dp,0.0_dp,2.0_dp,nse%gamma,background_test)
  call assert_vector_close(q(1,2,2,:),background_test,'periodic window left background')
  call assert_vector_close(q(8,2,2,:),background_test,'periodic window right background')
  call primitive_to_conserved(2.0_dp,2.0_dp,0.2_dp,0.4_dp,3.0_dp,nse%gamma,expected)
  call assert_vector_close(q(2,2,2,:),expected,'periodic window first cell')
  call assert_vector_close(q(6,2,2,:),expected,'periodic window repeated cell')
  call primitive_to_conserved(2.0_dp,3.0_dp,0.2_dp,0.4_dp,3.0_dp,nse%gamma,expected)
  call assert_vector_close(q(7,2,2,:),expected,'periodic window partial last block')
  nse%imported_turbulence_blend_cells=1
  call initialize_imported_turbulence(q,sim,nse,2,2,2,2)
  call assert_vector_close(q(7,2,2,:),0.5_dp*(background_test+expected),'window outer edge blend')
  call primitive_to_conserved(2.0_dp,2.0_dp,0.2_dp,0.4_dp,3.0_dp,nse%gamma,expected)
  call assert_vector_close(q(6,2,2,:),expected,'no blend at internal periodic seam')
  nse%imported_turbulence_blend_cells=0
  nse%imported_turbulence_x_length=2.0_dp
  call initialize_imported_turbulence(q,sim,nse,2,2,2,2,first_i,last_i)
  call assert_true(first_i==2 .and. last_i==3,'window shorter than source')
  call assert_vector_close(q(4,2,2,:),background_test,'short window background')
  nse%imported_turbulence_x_length=-1.0_dp

  nse%imported_turbulence_mode = 'embed'
  nse%imported_turbulence_x_start = 4.0_dp
  nse%imported_turbulence_velocity_offset_x = 0.0_dp
  nse%imported_turbulence_background_rho = 1.0_dp
  nse%imported_turbulence_background_u = 0.0_dp
  nse%imported_turbulence_background_v = 0.0_dp
  nse%imported_turbulence_background_w = 0.0_dp
  nse%imported_turbulence_background_p = 1.0_dp/nse%gamma
  nse%planar_shock_position = 2.0_dp
  nse%planar_shock_direction = 'positive_x'
  nse%planar_shock_upstream_rho = 1.0_dp
  nse%planar_shock_upstream_u = 0.0_dp
  nse%planar_shock_upstream_v = 0.0_dp
  nse%planar_shock_upstream_w = 0.0_dp
  nse%planar_shock_upstream_p = 1.0_dp/nse%gamma
  nse%planar_shock_downstream_rho = 2.0_dp
  nse%planar_shock_downstream_u = 0.75_dp
  nse%planar_shock_downstream_v = 0.0_dp
  nse%planar_shock_downstream_w = 0.0_dp
  nse%planar_shock_downstream_p = 1.5_dp
  nse%boundary_face_type = 'periodic'
  nse%boundary_face_type(nse_face_x_min) = 'dirichlet'
  nse%boundary_face_type(nse_face_x_max) = 'non_reflecting'
  nse%boundary_condition = 'mixed'
  nse%boundary_reference_rho(nse_face_x_min) = 2.0_dp
  nse%boundary_reference_velocity(:,nse_face_x_min) = &
    [0.75_dp,0.0_dp,0.0_dp]
  nse%boundary_reference_p(nse_face_x_min) = 1.5_dp
  q = -999.0_dp
  call initialize_shock_turbulence(q, sim, nse, 2, 2, 2, 2)
  call primitive_to_conserved(2.0_dp, 0.75_dp, 0.0_dp, 0.0_dp, &
    1.5_dp, nse%gamma, expected)
  call assert_vector_close(q(1,2,2,:), expected, &
    'shock downstream region')
  call assert_vector_close(q(2,2,2,:), expected, &
    'shock downstream last cell')
  call primitive_to_conserved(1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    1.0_dp/nse%gamma, nse%gamma, expected)
  call assert_vector_close(q(3,2,2,:), expected, &
    'shock upstream background gap')
  call primitive_to_conserved(2.0_dp, 1.0_dp, 0.2_dp, 0.4_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(5,2,2,:), expected, &
    'turbulence remains ahead of shock')

  nse%imported_turbulence_x_start = 0.0_dp
  nse%planar_shock_position = 6.0_dp
  nse%planar_shock_direction = 'negative_x'
  nse%planar_shock_downstream_u = -0.75_dp
  nse%boundary_face_type(nse_face_x_min) = 'non_reflecting'
  nse%boundary_face_type(nse_face_x_max) = 'dirichlet'
  nse%boundary_reference_rho(nse_face_x_max) = 2.0_dp
  nse%boundary_reference_velocity(:,nse_face_x_max) = &
    [-0.75_dp,0.0_dp,0.0_dp]
  nse%boundary_reference_p(nse_face_x_max) = 1.5_dp
  q = -999.0_dp
  call initialize_shock_turbulence(q, sim, nse, 2, 2, 2, 2)
  call primitive_to_conserved(2.0_dp, -0.75_dp, 0.0_dp, 0.0_dp, &
    1.5_dp, nse%gamma, expected)
  call assert_vector_close(q(7,2,2,:), expected, &
    'negative-x shock downstream first cell')
  call assert_vector_close(q(8,2,2,:), expected, &
    'negative-x shock downstream region')
  call primitive_to_conserved(1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    1.0_dp/nse%gamma, nse%gamma, expected)
  call assert_vector_close(q(5,2,2,:), expected, &
    'negative-x shock upstream background gap')
  call primitive_to_conserved(2.0_dp, 4.0_dp, 0.2_dp, 0.4_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(4,2,2,:), expected, &
    'turbulence remains ahead of negative-x shock')

  nse%imported_turbulence_mode='periodic_embed'
  nse%imported_turbulence_x_length=5.0_dp
  call initialize_shock_turbulence(q,sim,nse,2,2,2,2)
  call primitive_to_conserved(2.0_dp,1.0_dp,0.2_dp,0.4_dp,3.0_dp,nse%gamma,expected)
  call assert_vector_close(q(5,2,2,:),expected,'negative-x shock periodic window tail')
  nse%imported_turbulence_mode='embed'
  nse%imported_turbulence_x_length=-1.0_dp

  nse%imported_turbulence_x_start = 4.0_dp
  nse%imported_turbulence_background_rho = 1.0_dp
  nse%imported_turbulence_background_u = 0.0_dp
  nse%imported_turbulence_background_v = 0.0_dp
  nse%imported_turbulence_background_w = 0.0_dp
  nse%imported_turbulence_background_p = 1.0_dp/nse%gamma
  nse%shock_tube_diaphragm_position = 2.0_dp
  nse%shock_tube_driver_rho = 1.0_dp
  nse%shock_tube_driver_u = 0.0_dp
  nse%shock_tube_driver_v = 0.0_dp
  nse%shock_tube_driver_w = 0.0_dp
  nse%shock_tube_driver_p = 5.0_dp/nse%gamma
  nse%shock_tube_driven_rho = 1.0_dp
  nse%shock_tube_driven_u = 0.0_dp
  nse%shock_tube_driven_v = 0.0_dp
  nse%shock_tube_driven_w = 0.0_dp
  nse%shock_tube_driven_p = 1.0_dp/nse%gamma
  nse%boundary_face_type(nse_face_x_min) = 'reflective'
  nse%boundary_face_type(nse_face_x_max) = 'non_reflecting'
  nse%boundary_reference_rho(nse_face_x_max) = 1.0_dp
  nse%boundary_reference_velocity(:,nse_face_x_max) = 0.0_dp
  nse%boundary_reference_p(nse_face_x_max) = 1.0_dp/nse%gamma
  sim%initial_condition = 'shock_tube_turbulence_interaction'
  q = -999.0_dp
  call initialize_nse_state(q, sim, nse, 2, 2, 2, 2)
  call primitive_to_conserved(1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    5.0_dp/nse%gamma, nse%gamma, expected)
  call assert_vector_close(q(1,2,2,:), expected, &
    'shock-tube driver first cell')
  call assert_vector_close(q(2,2,2,:), expected, &
    'shock-tube driver last cell')
  call primitive_to_conserved(1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    1.0_dp/nse%gamma, nse%gamma, expected)
  call assert_vector_close(q(3,2,2,:), expected, &
    'shock-tube driven background gap')
  call primitive_to_conserved(2.0_dp, 1.0_dp, 0.2_dp, 0.4_dp, &
    3.0_dp, nse%gamma, expected)
  call assert_vector_close(q(5,2,2,:), expected, &
    'shock-tube turbulence remains in driven region')

  nse%imported_turbulence_mode='periodic_embed'
  nse%imported_turbulence_x_start=3.0_dp
  nse%imported_turbulence_x_length=5.0_dp
  call initialize_nse_state(q,sim,nse,2,2,2,2)
  call assert_vector_close(q(8,2,2,:),expected,'shock-tube periodic window repeated tail')

  call assert_true(imported_turbulence_weight(1, 8, 2) > 0.0_dp, &
    'blend edge weight is positive')
  call assert_true(imported_turbulence_weight(1, 8, 2) < 1.0_dp, &
    'blend edge weight is below one')
  call assert_close(imported_turbulence_weight(1, 8, 2), &
    imported_turbulence_weight(8, 8, 2), 'blend symmetry')
  call assert_close(imported_turbulence_weight(3, 8, 2), 1.0_dp, &
    'blend interior weight')

  deallocate(q)
  open(newunit=unit, file=source_file, status='old', iostat=ios)
  if (ios == 0) close(unit, status='delete')
  write(*,'(A)') 'Imported turbulence initialization tests passed'

contains

  subroutine write_source_slf(path)
    character(len=*), intent(in) :: path
    character(len=8) :: magic
    character(len=32) :: names(5)
    integer(int32) :: version, dtype_code, ndim, shape4(4), metadata(8)
    integer(int32) :: nvar
    real(dp) :: field(4,2,2,5), bounds(6), time
    real(dp) :: state(5)
    integer :: i, j, k, output_unit

    do k = 1, 2
      do j = 1, 2
        do i = 1, 4
          call primitive_to_conserved(2.0_dp, real(i,dp), &
            0.1_dp*real(j,dp), 0.2_dp*real(k,dp), 3.0_dp, 1.4_dp, state)
          field(i,j,k,:) = state
        end do
      end do
    end do

    magic = 'SLF1'//char(0)//char(0)//char(0)//char(0)
    version = 1_int32
    dtype_code = 2_int32
    ndim = 4_int32
    shape4 = [4_int32, 2_int32, 2_int32, 5_int32]
    metadata = [0_int32, 0_int32, 4_int32, 2_int32, 2_int32, &
      0_int32, 1_int32, 0_int32]
    time = 0.0_dp
    bounds = [0.0_dp, 4.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp]
    nvar = 5_int32
    names = [character(len=32) :: &
      'rho', 'rho_u', 'rho_v', 'rho_w', 'rho_E']

    open(newunit=output_unit, file=path, access='stream', &
      form='unformatted', status='replace', action='write', &
      convert='little_endian')
    write(output_unit) magic
    write(output_unit) version
    write(output_unit) dtype_code
    write(output_unit) ndim
    write(output_unit) shape4
    write(output_unit) metadata
    write(output_unit) time
    write(output_unit) bounds
    write(output_unit) nvar
    do i = 1, 5
      write(output_unit) names(i)
    end do
    write(output_unit) field
    close(output_unit)
  end subroutine write_source_slf

  pure subroutine primitive_to_conserved(rho, u, v, w, p, gamma, state)
    real(dp), intent(in) :: rho, u, v, w, p, gamma
    real(dp), intent(out) :: state(5)

    state(1) = rho
    state(2) = rho*u
    state(3) = rho*v
    state(4) = rho*w
    state(5) = p/(gamma-1.0_dp) + 0.5_dp*rho*(u*u+v*v+w*w)
  end subroutine primitive_to_conserved

  subroutine assert_vector_close(actual, reference, label)
    real(dp), intent(in) :: actual(:), reference(:)
    character(len=*), intent(in) :: label

    if (maxval(abs(actual-reference)) > 1.0e-12_dp) then
      write(*,'(A,A)') 'FAILED: ', trim(label)
      write(*,'(A,5ES16.8)') 'actual:   ', actual
      write(*,'(A,5ES16.8)') 'reference:', reference
      error stop 'imported turbulence vector assertion failed'
    end if
  end subroutine assert_vector_close

  subroutine assert_close(actual, reference, label)
    real(dp), intent(in) :: actual, reference
    character(len=*), intent(in) :: label

    if (abs(actual-reference) > 1.0e-12_dp) then
      write(*,'(A,A)') 'FAILED: ', trim(label)
      error stop 'imported turbulence scalar assertion failed'
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*,'(A,A)') 'FAILED: ', trim(label)
      error stop 'imported turbulence logical assertion failed'
    end if
  end subroutine assert_true

end program test_imported_turbulence
