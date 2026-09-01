module mod_init_shock_turbulence
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config, nse_face_x_min, nse_face_x_max
  use mod_init_imported_turbulence, only : initialize_imported_turbulence
  implicit none
  private

  integer, parameter :: nconserved = 5

  public :: initialize_shock_turbulence

contains

  subroutine initialize_shock_turbulence(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, &
      js-sim%nghost:, ks-sim%nghost:, :)

    character(len=32) :: direction
    real(dp) :: upstream(nconserved), downstream(nconserved)
    real(dp) :: aligned_position, scale
    integer :: first_i, last_i, face_index, driver_face
    integer :: i, j, k

    if (trim(adjustl(nse%imported_turbulence_mode)) /= 'embed') then
      error stop 'shock-turbulence interaction requires embed mode'
    end if
    if (nse%nv /= nconserved) then
      error stop 'shock-turbulence interaction requires five variables'
    end if
    if (nse%gamma <= 1.0_dp) then
      error stop 'shock-turbulence interaction requires gamma > 1'
    end if

    call primitive_to_conserved(nse%planar_shock_upstream_rho, &
      nse%planar_shock_upstream_u, nse%planar_shock_upstream_v, &
      nse%planar_shock_upstream_w, nse%planar_shock_upstream_p, &
      nse%gamma, upstream, 'planar shock upstream')
    call primitive_to_conserved(nse%planar_shock_downstream_rho, &
      nse%planar_shock_downstream_u, nse%planar_shock_downstream_v, &
      nse%planar_shock_downstream_w, nse%planar_shock_downstream_p, &
      nse%gamma, downstream, 'planar shock downstream')
    call validate_upstream_background(nse)

    face_index = nint((nse%planar_shock_position-sim%x_min)/sim%dx)
    aligned_position = sim%x_min + real(face_index,dp)*sim%dx
    scale = max(1.0_dp, abs(nse%planar_shock_position), &
      abs(aligned_position))
    if (abs(nse%planar_shock_position-aligned_position) > &
        1.0e-10_dp*scale) then
      error stop 'planar shock position must lie on an x-cell boundary'
    end if
    if (face_index < 1 .or. face_index >= sim%nx) then
      error stop 'planar shock position must lie strictly inside x domain'
    end if

    call initialize_imported_turbulence(q, sim, nse, js, je, ks, ke, &
      first_i, last_i)

    direction = lowercase(trim(adjustl(nse%planar_shock_direction)))
    select case (direction)
    case ('positive_x')
      driver_face = nse_face_x_min
      if (face_index >= first_i) then
        error stop 'positive-x shock initially overlaps imported turbulence'
      end if
      do k = ks, ke
        do j = js, je
          do i = 1, face_index
            q(i,j,k,1:nconserved) = downstream
          end do
        end do
      end do
    case ('negative_x')
      driver_face = nse_face_x_max
      if (face_index < last_i) then
        error stop 'negative-x shock initially overlaps imported turbulence'
      end if
      do k = ks, ke
        do j = js, je
          do i = face_index+1, sim%nx
            q(i,j,k,1:nconserved) = downstream
          end do
        end do
      end do
    case default
      error stop 'planar shock direction must be positive_x or negative_x'
    end select

    call validate_driver_boundary(nse, driver_face)
    if (sim%rank == 0) then
      write(*,'(A,A,A,ES16.8,A,I0)') &
        'Planar shock initialized: direction=', trim(direction), &
        ', position=', aligned_position, ', face_index=', face_index
    end if
  end subroutine initialize_shock_turbulence

  subroutine validate_upstream_background(nse)
    type(nse_config), intent(in) :: nse
    real(dp) :: background(5), upstream(5)

    background = [nse%imported_turbulence_background_rho, &
      nse%imported_turbulence_background_u, &
      nse%imported_turbulence_background_v, &
      nse%imported_turbulence_background_w, &
      nse%imported_turbulence_background_p]
    upstream = [nse%planar_shock_upstream_rho, &
      nse%planar_shock_upstream_u, nse%planar_shock_upstream_v, &
      nse%planar_shock_upstream_w, nse%planar_shock_upstream_p]
    if (.not. states_match(background, upstream)) then
      error stop 'imported turbulence background must equal shock upstream state'
    end if
  end subroutine validate_upstream_background

  subroutine validate_driver_boundary(nse, face)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face
    real(dp) :: boundary_state(5), downstream(5)

    if (trim(adjustl(nse%boundary_face_type(face))) /= 'dirichlet') then
      error stop 'shock driver face must use dirichlet boundary'
    end if
    boundary_state = [nse%boundary_reference_rho(face), &
      nse%boundary_reference_velocity(:,face), &
      nse%boundary_reference_p(face)]
    downstream = [nse%planar_shock_downstream_rho, &
      nse%planar_shock_downstream_u, nse%planar_shock_downstream_v, &
      nse%planar_shock_downstream_w, nse%planar_shock_downstream_p]
    if (.not. states_match(boundary_state, downstream)) then
      error stop 'shock driver dirichlet state must equal downstream state'
    end if
  end subroutine validate_driver_boundary

  subroutine primitive_to_conserved(rho, u, v, w, pressure, gamma, state, label)
    real(dp), intent(in) :: rho, u, v, w, pressure, gamma
    real(dp), intent(out) :: state(nconserved)
    character(len=*), intent(in) :: label

    if (.not. ieee_is_finite(rho) .or. rho <= 0.0_dp .or. &
        .not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp .or. &
        .not. all(ieee_is_finite([u,v,w]))) then
      write(*,'(A,A)') 'ERROR: invalid primitive state: ', trim(label)
      error stop 'invalid planar shock primitive state'
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

  pure function lowercase(value) result(lower)
    character(len=*), intent(in) :: value
    character(len=len(value)) :: lower
    integer :: code, index

    lower = value
    do index = 1, len(value)
      code = iachar(value(index:index))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(index:index) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase

end module mod_init_shock_turbulence
