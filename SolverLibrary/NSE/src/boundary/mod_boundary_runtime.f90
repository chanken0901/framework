module mod_nse_boundary
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_grid_fvm, only : axis_x,axis_y,axis_z
  use mod_grid_axis, only : grid_axis
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config, nse_boundary_face_count, &
    nse_face_x_min, nse_face_x_max, nse_face_y_min, nse_face_y_max, &
    nse_face_z_min, nse_face_z_max
  use module_mpi, only : ndiv_ny, ndiv_nz, itable, jjsta, kksta, j_sta, &
    mp_send_recv_pre_r8_Vec, mp_sendrecv_r8, MPI_COMM_WORLD, &
    MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE
  implicit none
  private

  character(len=5), parameter :: face_name(nse_boundary_face_count) = [ &
    character(len=5) :: 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']

  public :: apply_nse_boundary
  public :: validate_boundary_scheme
  public :: boundary_required_ghost_cells
  public :: boundary_scheme_name
  public :: boundary_ghost_distance

contains

  real(dp) function boundary_ghost_distance(axis,high,layer) result(distance)
    type(grid_axis), intent(in) :: axis
    logical, intent(in) :: high
    integer, intent(in) :: layer
    if(.not.allocated(axis%center)) error stop 'Boundary grid coordinates are missing'
    if(layer<1.or.layer>axis%ng) error stop 'Invalid boundary ghost layer'
    if(high) then
      distance=axis%center(axis%n+layer)-axis%center(axis%n)
    else
      distance=axis%center(1)-axis%center(1-layer)
    end if
    if(.not.ieee_is_finite(distance).or.distance<=0) error stop 'Invalid boundary ghost distance'
  end function

  subroutine apply_nse_boundary(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)

    call apply_boundary_x(q, sim, nse, js, je, ks, ke)

    !$OMP MASKED
    call apply_boundary_y(q, sim, nse, js, je, ks, ke)
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP MASKED
    call apply_boundary_z(q, sim, nse, js, je, ks, ke)
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP MASKED
    if(sim%nghost==6) then
      call exchange_wide_halo(q,sim,js,je,ks,ke)
      call exchange_wide_halo(q,sim,js,je,ks,ke)
    else
    call mp_send_recv_pre_r8_Vec(q, sim%nghost, &
      1-sim%nghost, sim%nx+sim%nghost, &
      js-sim%nghost, je+sim%nghost, ks-sim%nghost, ke+sim%nghost)
    ! The second staged exchange propagates completed physical-face values
    ! through transverse MPI halos for sixth-order mixed derivatives.
    if (trim(adjustl(nse%viscous_scheme)) == 'central6' .or. &
        trim(adjustl(nse%viscous_scheme)) == 'fv2') then
      call mp_send_recv_pre_r8_Vec(q, sim%nghost, &
        1-sim%nghost, sim%nx+sim%nghost, &
        js-sim%nghost, je+sim%nghost, ks-sim%nghost, ke+sim%nghost)
    end if
    end if
    !$OMP END MASKED
    !$OMP BARRIER
  end subroutine apply_nse_boundary

  subroutine exchange_wide_halo(q,sim,js,je,ks,ke)
    ! Six-layer exchange for mapped D6(D6). Preserve the legacy three-layer MPI path.
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: js,je,ks,ke
    real(dp), intent(inout) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), allocatable :: sendbuf(:),recvbuf(:)
    integer :: row,col,r,c,phase,partner,count,ierr,lo,hi,jlo,jhi,klo,khi
    integer :: status(MPI_STATUS_SIZE)
    row=0;col=0
    do r=0,ndiv_ny-1
      if(jjsta(r)==js) row=r
    end do
    do c=0,ndiv_nz-1
      if(kksta(c)==ks) col=c
    end do
    lo=1-6;hi=sim%nx+6;jlo=js-6;jhi=je+6;klo=ks-6;khi=ke+6
    count=(hi-lo+1)*6*(khi-klo+1)*5
    allocate(sendbuf(count),recvbuf(count))
    ! Two disjoint edge colours: exchange all independent neighbours together.
    do phase=0,1
      r=row
      if(modulo(row,2)/=phase) r=row-1
      if(r<0.or.r>=ndiv_ny-1) cycle
      if(row==r) then
        partner=itable(r+1,col);sendbuf=reshape(q(:,je-5:je,:,:),[count])
      else
        partner=itable(r,col);sendbuf=reshape(q(:,js:js+5,:,:),[count])
      end if
      call mp_sendrecv_r8(sendbuf(1),count,MPI_DOUBLE_PRECISION,partner,71, &
        recvbuf(1),count,MPI_DOUBLE_PRECISION,partner,71,MPI_COMM_WORLD,status,ierr)
      if(ierr/=0) error stop 'Six-layer Y halo exchange failed'
      if(row==r) then
        q(:,je+1:je+6,:,:)=reshape(recvbuf,[hi-lo+1,6,khi-klo+1,5])
      else
        q(:,js-6:js-1,:,:)=reshape(recvbuf,[hi-lo+1,6,khi-klo+1,5])
      end if
    end do
    deallocate(sendbuf,recvbuf)
    count=(hi-lo+1)*(jhi-jlo+1)*6*5
    allocate(sendbuf(count),recvbuf(count))
    do phase=0,1
      c=col
      if(modulo(col,2)/=phase) c=col-1
      if(c<0.or.c>=ndiv_nz-1) cycle
      if(col==c) then
        partner=itable(row,c+1);sendbuf=reshape(q(:,:,ke-5:ke,:),[count])
      else
        partner=itable(row,c);sendbuf=reshape(q(:,:,ks:ks+5,:),[count])
      end if
      call mp_sendrecv_r8(sendbuf(1),count,MPI_DOUBLE_PRECISION,partner,72, &
        recvbuf(1),count,MPI_DOUBLE_PRECISION,partner,72,MPI_COMM_WORLD,status,ierr)
      if(ierr/=0) error stop 'Six-layer Z halo exchange failed'
      if(col==c) then
        q(:,:,ke+1:ke+6,:)=reshape(recvbuf,[hi-lo+1,jhi-jlo+1,6,5])
      else
        q(:,:,ks-6:ks-1,:)=reshape(recvbuf,[hi-lo+1,jhi-jlo+1,6,5])
      end if
    end do
  end subroutine

  subroutine apply_boundary_x(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: ghost_state(nse%nv)
    integer :: g, j, k

    !$OMP DO collapse(2) schedule(static) private(g,ghost_state)
    do k = ks, ke
      do j = js, je
        do g = 1, sim%nghost
          if (is_periodic(nse, nse_face_x_min)) then
            q(1-g,j,k,:) = q(sim%nx+1-g,j,k,:)
          else if (is_reflective(nse, nse_face_x_min)) then
            call reflective_state(q(g,j,k,:), ghost_state, 1)
            q(1-g,j,k,:) = ghost_state
          else if (is_dirichlet(nse, nse_face_x_min)) then
            call reference_state(nse, nse_face_x_min, ghost_state)
            q(1-g,j,k,:) = ghost_state
          else
            call non_reflecting_state(q(1,j,k,:), ghost_state, sim, nse, &
              nse_face_x_min, 1, -1.0_dp, g)
            q(1-g,j,k,:) = ghost_state
          end if
          if (is_periodic(nse, nse_face_x_max)) then
            q(sim%nx+g,j,k,:) = q(g,j,k,:)
          else if (is_reflective(nse, nse_face_x_max)) then
            call reflective_state(q(sim%nx+1-g,j,k,:), ghost_state, 1)
            q(sim%nx+g,j,k,:) = ghost_state
          else if (is_dirichlet(nse, nse_face_x_max)) then
            call reference_state(nse, nse_face_x_max, ghost_state)
            q(sim%nx+g,j,k,:) = ghost_state
          else
            call non_reflecting_state(q(sim%nx,j,k,:), ghost_state, sim, &
              nse, nse_face_x_max, 1, 1.0_dp, g)
            q(sim%nx+g,j,k,:) = ghost_state
          end if
        end do
      end do
    end do
    !$OMP END DO
  end subroutine apply_boundary_x

  subroutine apply_boundary_y(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: ghost_state(nse%nv)
    integer :: g, i, k, ilo, ihi

    if (is_periodic(nse, nse_face_y_min)) then
      call apply_periodic_y(q, sim, nse, js, je, ks, ke)
      return
    end if

    ilo = 1 - sim%nghost
    ihi = sim%nx + sim%nghost
    if (js == 1) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            if (is_reflective(nse, nse_face_y_min)) then
              call reflective_state(q(i,g,k,:), ghost_state, 2)
            else if (is_dirichlet(nse, nse_face_y_min)) then
              call reference_state(nse, nse_face_y_min, ghost_state)
            else
              call non_reflecting_state(q(i,1,k,:), ghost_state, sim, nse, &
                nse_face_y_min, 2, -1.0_dp, g)
            end if
            q(i,1-g,k,:) = ghost_state
          end do
        end do
      end do
    end if
    if (je == sim%ny) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            if (is_reflective(nse, nse_face_y_max)) then
              call reflective_state(q(i,sim%ny+1-g,k,:), ghost_state, 2)
            else if (is_dirichlet(nse, nse_face_y_max)) then
              call reference_state(nse, nse_face_y_max, ghost_state)
            else
              call non_reflecting_state(q(i,sim%ny,k,:), ghost_state, sim, &
                nse, nse_face_y_max, 2, 1.0_dp, g)
            end if
            q(i,sim%ny+g,k,:) = ghost_state
          end do
        end do
      end do
    end if
  end subroutine apply_boundary_y

  subroutine apply_boundary_z(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: ghost_state(nse%nv)
    integer :: g, i, j, ilo, ihi, jlo, jhi

    if (is_periodic(nse, nse_face_z_min)) then
      call apply_periodic_z(q, sim, nse, js, je, ks, ke)
      return
    end if

    ilo = 1 - sim%nghost
    ihi = sim%nx + sim%nghost
    ! Transverse y halos at internal rank boundaries are not valid until the
    ! staged MPI exchange below.  Include only physical y ghosts, which were
    ! completed by apply_boundary_y, and let that exchange propagate z-face
    ! values to internal y halos.
    jlo = js
    jhi = je
    if (js == 1) jlo = js - sim%nghost
    if (je == sim%ny) jhi = je + sim%nghost
    if (ks == 1) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            if (is_reflective(nse, nse_face_z_min)) then
              call reflective_state(q(i,j,g,:), ghost_state, 3)
            else if (is_dirichlet(nse, nse_face_z_min)) then
              call reference_state(nse, nse_face_z_min, ghost_state)
            else
              call non_reflecting_state(q(i,j,1,:), ghost_state, sim, nse, &
                nse_face_z_min, 3, -1.0_dp, g)
            end if
            q(i,j,1-g,:) = ghost_state
          end do
        end do
      end do
    end if
    if (ke == sim%nz) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            if (is_reflective(nse, nse_face_z_max)) then
              call reflective_state(q(i,j,sim%nz+1-g,:), ghost_state, 3)
            else if (is_dirichlet(nse, nse_face_z_max)) then
              call reference_state(nse, nse_face_z_max, ghost_state)
            else
              call non_reflecting_state(q(i,j,sim%nz,:), ghost_state, sim, &
                nse, nse_face_z_max, 3, 1.0_dp, g)
            end if
            q(i,j,sim%nz+g,:) = ghost_state
          end do
        end do
      end do
    end if
  end subroutine apply_boundary_z

  pure subroutine reflective_state(q_inside, q_ghost, normal_axis)
    real(dp), intent(in) :: q_inside(:)
    real(dp), intent(out) :: q_ghost(:)
    integer, intent(in) :: normal_axis

    q_ghost = q_inside
    q_ghost(normal_axis+1) = -q_ghost(normal_axis+1)
  end subroutine reflective_state

  subroutine reference_state(nse, face, q_ghost)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face
    real(dp), intent(out) :: q_ghost(:)

    call primitive_to_conserved(nse%boundary_reference_rho(face), &
      nse%boundary_reference_velocity(:,face), &
      nse%boundary_reference_p(face), q_ghost, nse, face)
  end subroutine reference_state

  subroutine non_reflecting_state(q_inside, q_ghost, sim, nse, face, &
      normal_axis, outward_sign, ghost_layer)
    real(dp), intent(in) :: q_inside(:)
    real(dp), intent(out) :: q_ghost(:)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face, normal_axis, ghost_layer
    real(dp), intent(in) :: outward_sign
    real(dp) :: rho_i, pressure_i, velocity_i(3), sound_i
    real(dp) :: rho_r, pressure_r, velocity_r(3), sound_r
    real(dp) :: normal_i, normal_r, normal_b, velocity_b(3)
    real(dp) :: jminus_i, jplus_i, jminus_r, jplus_r
    real(dp) :: jminus_b, jplus_b, entropy_i, entropy_r, entropy_b
    real(dp) :: sound_b, rho_b, pressure_b, alpha
    real(dp) :: spacing, length_scale, distance
    integer :: component

    call conserved_to_primitive(q_inside, rho_i, velocity_i, pressure_i, &
      nse, face)
    rho_r = nse%boundary_reference_rho(face)
    velocity_r = nse%boundary_reference_velocity(:,face)
    pressure_r = nse%boundary_reference_p(face)
    sound_i = sqrt(nse%gamma * pressure_i / rho_i)
    sound_r = sqrt(nse%gamma * pressure_r / rho_r)
    normal_i = outward_sign * velocity_i(normal_axis)
    normal_r = outward_sign * velocity_r(normal_axis)
    jminus_i = normal_i - 2.0_dp*sound_i/(nse%gamma-1.0_dp)
    jplus_i = normal_i + 2.0_dp*sound_i/(nse%gamma-1.0_dp)
    jminus_r = normal_r - 2.0_dp*sound_r/(nse%gamma-1.0_dp)
    jplus_r = normal_r + 2.0_dp*sound_r/(nse%gamma-1.0_dp)
    entropy_i = pressure_i / rho_i**nse%gamma
    entropy_r = pressure_r / rho_r**nse%gamma

    select case (normal_axis)
    case (1)
      spacing = (sim%x_max-sim%x_min) / real(sim%nx,dp)
      length_scale = sim%x_max-sim%x_min
    case (2)
      spacing = (sim%y_max-sim%y_min) / real(sim%ny,dp)
      length_scale = sim%y_max-sim%y_min
    case (3)
      spacing = (sim%z_max-sim%z_min) / real(sim%nz,dp)
      length_scale = sim%z_max-sim%z_min
    case default
      error stop 'invalid non-reflecting boundary normal axis'
    end select
    distance=real(ghost_layer,dp)*spacing
    if(sim%grid_mapping=='sinh') then
      select case(normal_axis)
      case(1);distance=boundary_ghost_distance(axis_x,outward_sign>0,ghost_layer)
      case(2);distance=boundary_ghost_distance(axis_y,outward_sign>0,ghost_layer)
      case(3);distance=boundary_ghost_distance(axis_z,outward_sign>0,ghost_layer)
      end select
    end if
    if (nse%boundary_length_scale > 0.0_dp) then
      length_scale = nse%boundary_length_scale
    end if
    alpha = 1.0_dp - exp(-nse%boundary_relaxation_strength * &
      distance / length_scale)
    alpha = max(0.0_dp, min(1.0_dp, alpha))
    if (normal_i + sound_i < 0.0_dp) alpha = 1.0_dp

    jminus_b = outgoing_or_relaxed(jminus_i, jminus_r, &
      normal_i-sound_i, alpha)
    jplus_b = outgoing_or_relaxed(jplus_i, jplus_r, &
      normal_i+sound_i, alpha)
    entropy_b = outgoing_or_relaxed(entropy_i, entropy_r, normal_i, alpha)
    velocity_b = velocity_i
    do component = 1, 3
      if (component == normal_axis) cycle
      velocity_b(component) = outgoing_or_relaxed(velocity_i(component), &
        velocity_r(component), normal_i, alpha)
    end do

    normal_b = 0.5_dp * (jplus_b+jminus_b)
    sound_b = 0.25_dp * (nse%gamma-1.0_dp) * (jplus_b-jminus_b)
    if (.not. ieee_is_finite(sound_b) .or. sound_b <= 0.0_dp .or. &
        .not. ieee_is_finite(entropy_b) .or. entropy_b <= 0.0_dp) then
      call invalid_boundary_state(face, 'invalid characteristic sound speed/entropy')
    end if
    rho_b = (sound_b*sound_b/(nse%gamma*entropy_b))** &
      (1.0_dp/(nse%gamma-1.0_dp))
    pressure_b = entropy_b * rho_b**nse%gamma
    velocity_b(normal_axis) = outward_sign * normal_b
    call primitive_to_conserved(rho_b, velocity_b, pressure_b, q_ghost, &
      nse, face)
  end subroutine non_reflecting_state

  pure real(dp) function outgoing_or_relaxed(interior, reference, &
      eigenvalue, alpha) result(value)
    real(dp), intent(in) :: interior, reference, eigenvalue, alpha

    if (eigenvalue >= 0.0_dp) then
      value = interior
    else
      value = (1.0_dp-alpha)*interior + alpha*reference
    end if
  end function outgoing_or_relaxed

  subroutine conserved_to_primitive(qstate, rho, velocity, pressure, nse, face)
    real(dp), intent(in) :: qstate(:)
    real(dp), intent(out) :: rho, velocity(3), pressure
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face
    real(dp) :: kinetic

    rho = qstate(1)
    if (.not. ieee_is_finite(rho) .or. rho <= nse%small_rho) then
      call invalid_boundary_state(face, 'non-positive interior density')
    end if
    velocity = qstate(2:4) / rho
    kinetic = 0.5_dp*rho*sum(velocity*velocity)
    pressure = (nse%gamma-1.0_dp) * (qstate(5)-kinetic)
    if (.not. all(ieee_is_finite(velocity)) .or. &
        .not. ieee_is_finite(pressure) .or. pressure <= nse%small_p) then
      call invalid_boundary_state(face, 'invalid interior velocity/pressure')
    end if
  end subroutine conserved_to_primitive

  subroutine primitive_to_conserved(rho, velocity, pressure, qstate, nse, face)
    real(dp), intent(in) :: rho, velocity(3), pressure
    real(dp), intent(out) :: qstate(:)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face

    if (.not. ieee_is_finite(rho) .or. rho <= nse%small_rho .or. &
        .not. all(ieee_is_finite(velocity)) .or. &
        .not. ieee_is_finite(pressure) .or. pressure <= nse%small_p) then
      call invalid_boundary_state(face, 'invalid reconstructed primitive state')
    end if
    qstate(1) = rho
    qstate(2:4) = rho*velocity
    qstate(5) = pressure/(nse%gamma-1.0_dp) + &
      0.5_dp*rho*sum(velocity*velocity)
  end subroutine primitive_to_conserved

  subroutine invalid_boundary_state(face, reason)
    integer, intent(in) :: face
    character(len=*), intent(in) :: reason

    write(*,'(A,A,A,A)') 'ERROR: boundary ', &
      trim(face_name(face)), ': ', trim(reason)
    error stop 'invalid boundary state'
  end subroutine invalid_boundary_state

  logical function is_periodic(nse, face) result(periodic)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face

    periodic = trim(adjustl(nse%boundary_face_type(face))) == 'periodic'
  end function is_periodic

  logical function is_reflective(nse, face) result(reflective)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face

    reflective = trim(adjustl(nse%boundary_face_type(face))) == 'reflective'
  end function is_reflective

  logical function is_dirichlet(nse, face) result(dirichlet)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: face

    dirichlet = trim(adjustl(nse%boundary_face_type(face))) == 'dirichlet'
  end function is_dirichlet

  subroutine apply_periodic_y(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), allocatable :: sendbuf(:,:,:,:), recv_high(:,:,:,:)
    real(dp), allocatable :: recv_low(:,:,:,:)
    integer :: g, i, k, column, partner, count, ierr
    integer :: ilo, ihi
    integer :: status(MPI_STATUS_SIZE)

    ilo = 1 - sim%nghost
    ihi = sim%nx + sim%nghost
    if (ndiv_ny == 1) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,sim%ny+g,k,:) = q(i,g,k,:)
            q(i,1-g,k,:) = q(i,sim%ny+1-g,k,:)
          end do
        end do
      end do
      return
    end if
    if (js /= 1 .and. je /= sim%ny) return

    allocate(sendbuf(ilo:ihi,1:sim%nghost,ks:ke,1:nse%nv))
    allocate(recv_high(ilo:ihi,1:sim%nghost,ks:ke,1:nse%nv))
    allocate(recv_low(ilo:ihi,1:sim%nghost,ks:ke,1:nse%nv))
    count = (ihi-ilo+1) * sim%nghost * (ke-ks+1) * nse%nv
    do column = 0, ndiv_nz-1
      if (ks /= kksta(column)) cycle
      if (je == sim%ny) then
        partner = itable(0,column)
        do k = ks, ke
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,g,k,:) = q(i,sim%ny+g-sim%nghost,k,:)
            end do
          end do
        end do
        call mp_sendrecv_r8(sendbuf(ilo,1,ks,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, recv_high(ilo,1,ks,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic y-high exchange failed'
      end if
      if (js == 1) then
        partner = itable(ndiv_ny-1,column)
        do k = ks, ke
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,g,k,:) = q(i,g,k,:)
            end do
          end do
        end do
        call mp_sendrecv_r8(sendbuf(ilo,1,ks,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, recv_low(ilo,1,ks,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic y-low exchange failed'
      end if
    end do
    if (je == sim%ny) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,sim%ny+g,k,:) = recv_high(i,g,k,:)
          end do
        end do
      end do
    end if
    if (js == 1) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,1-g,k,:) = recv_low(i,sim%nghost-g+1,k,:)
          end do
        end do
      end do
    end if
    deallocate(sendbuf, recv_high, recv_low)
  end subroutine apply_periodic_y

  subroutine apply_periodic_z(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), allocatable :: sendbuf(:,:,:,:), recv_high(:,:,:,:)
    real(dp), allocatable :: recv_low(:,:,:,:)
    integer :: g, i, j, row, partner, count, ierr
    integer :: ilo, ihi, jlo, jhi
    integer :: status(MPI_STATUS_SIZE)

    ilo = 1 - sim%nghost
    ihi = sim%nx + sim%nghost
    jlo = js - sim%nghost
    jhi = je + sim%nghost
    if (ndiv_nz == 1) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,j,sim%nz+g,:) = q(i,j,g,:)
            q(i,j,1-g,:) = q(i,j,sim%nz+1-g,:)
          end do
        end do
      end do
      return
    end if
    if (ks /= 1 .and. ke /= sim%nz) return

    allocate(sendbuf(ilo:ihi,jlo:jhi,1:sim%nghost,1:nse%nv))
    allocate(recv_high(ilo:ihi,jlo:jhi,1:sim%nghost,1:nse%nv))
    allocate(recv_low(ilo:ihi,jlo:jhi,1:sim%nghost,1:nse%nv))
    count = (ihi-ilo+1) * (jhi-jlo+1) * sim%nghost * nse%nv
    do row = 0, ndiv_ny-1
      if (j_sta /= jjsta(row)) cycle
      if (ke == sim%nz) then
        partner = itable(row,0)
        do j = jlo, jhi
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,j,g,:) = q(i,j,sim%nz+g-sim%nghost,:)
            end do
          end do
        end do
        call mp_sendrecv_r8(sendbuf(ilo,jlo,1,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, recv_high(ilo,jlo,1,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic z-high exchange failed'
      end if
      if (ks == 1) then
        partner = itable(row,ndiv_nz-1)
        do j = jlo, jhi
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,j,g,:) = q(i,j,g,:)
            end do
          end do
        end do
        call mp_sendrecv_r8(sendbuf(ilo,jlo,1,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, recv_low(ilo,jlo,1,1), count, &
          MPI_DOUBLE_PRECISION, partner, 1, MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic z-low exchange failed'
      end if
    end do
    if (ke == sim%nz) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,j,sim%nz+g,:) = recv_high(i,j,g,:)
          end do
        end do
      end do
    end if
    if (ks == 1) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,j,1-g,:) = recv_low(i,j,sim%nghost-g+1,:)
          end do
        end do
      end do
    end if
    deallocate(sendbuf, recv_high, recv_low)
  end subroutine apply_periodic_z

  subroutine validate_boundary_scheme(sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer :: face
    logical :: all_periodic

    if (sim%nghost /= boundary_required_ghost_cells().and. &
        .not.(sim%mapped_keep6.and.sim%nghost==6.and.nse%viscous_scheme=='central6')) then
      error stop 'Use three ghosts, or six for mapped CENTRAL6'
    end if
    if (nse%nv /= 5) then
      error stop 'MPI halo exchange currently requires five conserved variables'
    end if
    if (nse%gamma <= 1.0_dp) error stop 'boundary treatment requires gamma > 1'
    if (.not. ieee_is_finite(nse%boundary_relaxation_strength) .or. &
        nse%boundary_relaxation_strength < 0.0_dp) then
      error stop 'boundary relaxation strength must be finite and non-negative'
    end if
    if (.not. ieee_is_finite(nse%boundary_length_scale) .or. &
        abs(nse%boundary_length_scale) <= tiny(1.0_dp) .or. &
        (nse%boundary_length_scale < 0.0_dp .and. &
         abs(nse%boundary_length_scale+1.0_dp) > 10.0_dp*epsilon(1.0_dp))) then
      error stop 'boundary length scale must be AUTO or positive'
    end if

    do face = 1, nse_boundary_face_count
      select case (trim(adjustl(nse%boundary_face_type(face))))
      case ('periodic')
        continue
      case ('reflective')
        continue
      case ('non_reflecting', 'dirichlet')
        if (.not. ieee_is_finite(nse%boundary_reference_rho(face)) .or. &
            nse%boundary_reference_rho(face) <= nse%small_rho) then
          call invalid_boundary_state(face, 'reference density must be positive')
        end if
        if (.not. all(ieee_is_finite( &
            nse%boundary_reference_velocity(:,face)))) then
          call invalid_boundary_state(face, 'reference velocity must be finite')
        end if
        if (.not. ieee_is_finite(nse%boundary_reference_p(face)) .or. &
            nse%boundary_reference_p(face) <= nse%small_p) then
          call invalid_boundary_state(face, 'reference pressure must be positive')
        end if
      case default
        write(*,'(A,A,A)') 'ERROR: unsupported boundary type on ', &
          trim(face_name(face)), ': ' // trim(nse%boundary_face_type(face))
        error stop 'unsupported boundary type'
      end select
    end do
    call validate_periodic_pair(nse, nse_face_x_min, nse_face_x_max, 'x')
    call validate_periodic_pair(nse, nse_face_y_min, nse_face_y_max, 'y')
    call validate_periodic_pair(nse, nse_face_z_min, nse_face_z_max, 'z')

    all_periodic = all(nse%boundary_face_type == 'periodic')
    if (all_periodic .neqv. &
        (trim(adjustl(nse%boundary_condition)) == 'periodic')) then
      error stop 'boundary summary is inconsistent with six face types'
    end if
  end subroutine validate_boundary_scheme

  subroutine validate_periodic_pair(nse, lower_face, upper_face, direction)
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: lower_face, upper_face
    character(len=*), intent(in) :: direction

    if (is_periodic(nse, lower_face) .neqv. is_periodic(nse, upper_face)) then
      write(*,'(A,A,A)') 'ERROR: periodic ', trim(direction), &
        ' boundaries must be specified on both physical faces'
      error stop 'unpaired periodic boundary'
    end if
  end subroutine validate_periodic_pair

  integer function boundary_required_ghost_cells() result(nghost)
    nghost = 3
  end function boundary_required_ghost_cells

  pure function boundary_scheme_name() result(name)
    character(len=32) :: name
    name = 'runtime'
  end function boundary_scheme_name

end module mod_nse_boundary
