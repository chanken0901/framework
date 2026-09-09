module mod_init_shock_tube_turbulence
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config, nse_face_x_min, nse_face_x_max
  use mod_init_imported_turbulence, only : initialize_imported_turbulence
  implicit none
  private

  integer, parameter :: nconserved = 5

  public :: initialize_shock_tube_turbulence

contains

  subroutine initialize_shock_tube_turbulence(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, &
      js-sim%nghost:, ks-sim%nghost:, :)

    real(dp) :: driver(nconserved), driven(nconserved)
    real(dp) :: aligned_position, scale, velocity_scale
    integer :: first_i, last_i, face_index
    integer :: i, j, k

    if (trim(adjustl(nse%imported_turbulence_mode)) /= 'embed' .and. &
        trim(adjustl(nse%imported_turbulence_mode)) /= 'periodic_embed') then
      error stop 'shock-tube turbulence interaction requires embed or periodic_embed mode'
    end if
    if (nse%nv /= nconserved) then
      error stop 'shock-tube turbulence interaction requires five variables'
    end if
    if (nse%gamma <= 1.0_dp) then
      error stop 'shock-tube turbulence interaction requires gamma > 1'
    end if

    call primitive_to_conserved(nse%shock_tube_driver_rho, &
      nse%shock_tube_driver_u, nse%shock_tube_driver_v, &
      nse%shock_tube_driver_w, nse%shock_tube_driver_p, nse%gamma, &
      driver, 'shock-tube driver')
    call primitive_to_conserved(nse%shock_tube_driven_rho, &
      nse%shock_tube_driven_u, nse%shock_tube_driven_v, &
      nse%shock_tube_driven_w, nse%shock_tube_driven_p, nse%gamma, &
      driven, 'shock-tube driven region')
    if (nse%shock_tube_driver_p <= nse%shock_tube_driven_p) then
      error stop 'shock-tube driver pressure must exceed driven pressure'
    end if
    velocity_scale = max(1.0_dp, abs(nse%shock_tube_driver_u), &
      abs(nse%shock_tube_driver_v), abs(nse%shock_tube_driver_w))
    if (abs(nse%shock_tube_driver_u) > 1.0e-12_dp*velocity_scale) then
      error stop 'reflective shock-tube driver requires zero x velocity'
    end if
    call validate_driven_background(nse)

    face_index = nint((nse%shock_tube_diaphragm_position-sim%x_min) / &
      sim%dx)
    aligned_position = sim%x_min + real(face_index,dp)*sim%dx
    scale = max(1.0_dp, abs(nse%shock_tube_diaphragm_position), &
      abs(aligned_position))
    if (abs(nse%shock_tube_diaphragm_position-aligned_position) > &
        1.0e-10_dp*scale) then
      error stop 'shock-tube diaphragm must lie on an x-cell boundary'
    end if
    if (face_index < 1 .or. face_index >= sim%nx) then
      error stop 'shock-tube diaphragm must lie strictly inside x domain'
    end if

    call initialize_imported_turbulence(q, sim, nse, js, je, ks, ke, &
      first_i, last_i)
    if (face_index >= first_i) then
      error stop 'shock-tube driver overlaps imported turbulence'
    end if

    do k = ks, ke
      do j = js, je
        do i = 1, face_index
          q(i,j,k,1:nconserved) = driver
        end do
      end do
    end do

    call validate_boundaries(nse)
    if (sim%rank == 0) then
      write(*,'(A,ES16.8,A,I0,A,I0,A,I0)') &
        'Shock-tube turbulence initialized: diaphragm=', aligned_position, &
        ', face_index=', face_index, ', turbulence_i=', first_i, ':', last_i
    end if
  end subroutine initialize_shock_tube_turbulence

  subroutine validate_driven_background(nse)
    type(nse_config), intent(in) :: nse
    real(dp) :: background(5), driven(5)

    background = [nse%imported_turbulence_background_rho, &
      nse%imported_turbulence_background_u, &
      nse%imported_turbulence_background_v, &
      nse%imported_turbulence_background_w, &
      nse%imported_turbulence_background_p]
    driven = [nse%shock_tube_driven_rho, nse%shock_tube_driven_u, &
      nse%shock_tube_driven_v, nse%shock_tube_driven_w, &
      nse%shock_tube_driven_p]
    if (.not. states_match(background, driven)) then
      error stop 'imported turbulence background must equal shock-tube driven state'
    end if
  end subroutine validate_driven_background

  subroutine validate_boundaries(nse)
    type(nse_config), intent(in) :: nse
    real(dp) :: boundary_state(5), driven(5)

    if (trim(adjustl(nse%boundary_face_type(nse_face_x_min))) /= &
        'reflective') then
      error stop 'shock-tube x_min boundary must be reflective'
    end if
    if (trim(adjustl(nse%boundary_face_type(nse_face_x_max))) /= &
        'non_reflecting') then
      error stop 'shock-tube x_max boundary must be non_reflecting'
    end if
    boundary_state = [nse%boundary_reference_rho(nse_face_x_max), &
      nse%boundary_reference_velocity(:,nse_face_x_max), &
      nse%boundary_reference_p(nse_face_x_max)]
    driven = [nse%shock_tube_driven_rho, nse%shock_tube_driven_u, &
      nse%shock_tube_driven_v, nse%shock_tube_driven_w, &
      nse%shock_tube_driven_p]
    if (.not. states_match(boundary_state, driven)) then
      error stop 'shock-tube x_max reference must equal driven state'
    end if
  end subroutine validate_boundaries

  subroutine primitive_to_conserved(rho, u, v, w, pressure, gamma, &
      state, label)
    real(dp), intent(in) :: rho, u, v, w, pressure, gamma
    real(dp), intent(out) :: state(nconserved)
    character(len=*), intent(in) :: label

    if (.not. ieee_is_finite(rho) .or. rho <= 0.0_dp .or. &
        .not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp .or. &
        .not. all(ieee_is_finite([u,v,w]))) then
      write(*,'(A,A)') 'ERROR: invalid primitive state: ', trim(label)
      error stop 'invalid shock-tube primitive state'
    end if
    state(1) = rho
    state(2) = rho*u
    state(3) = rho*v
    state(4) = rho*w
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*(u*u+v*v+w*w)
  end subroutine primitive_to_conserved

  pure logical function states_match(left, right) result(matches)
    real(dp), intent(in) :: left(5), right(5)
    real(dp) :: scale(5)

    scale = max(1.0_dp, abs(left), abs(right))
    matches = all(abs(left-right) <= 1.0e-11_dp*scale)
  end function states_match

end module mod_init_shock_tube_turbulence
