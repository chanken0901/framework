module mod_mc_euler_flux
  use mod_precision, only : dp
  use mod_mc_mapped_flux, only: mc_mapped_rhs, mc_mapped_rate
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config, &
    mc_face_x_min, mc_face_x_max, mc_face_y_min, mc_face_y_max, &
    mc_face_z_min, mc_face_z_max
  use mod_mc_euler_field, only : validate_mc_euler_state, &
    mc_primitive_workspace, prepare_mc_primitives, evaluate_mc_primitive
  use mod_mc_boundary, only : mc_boundary_state, mc_boundary_is_periodic
  use mod_mc_thermodynamics_provider, only : mc_mixture_density, &
    mc_pressure, mc_sound_speed
  implicit none
  private

  public :: compute_mc_euler_physical_flux
  public :: compute_mc_euler_rusanov_flux
  public :: compute_mc_euler_rhs, compute_mc_euler_rhs_prepared
  public :: compute_mc_euler_timestep
  public :: advance_mc_euler_ssprk3

contains

  subroutine compute_mc_euler_physical_flux( &
      state, layout, gamma, direction, flux)
    real(dp), intent(in) :: state(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    integer, intent(in) :: direction
    real(dp), intent(out) :: flux(:)
    real(dp) :: density, velocity(3), normal_velocity, pressure

    density = mc_mixture_density(state, layout)
    velocity = state(layout%momentum)/density
    normal_velocity = velocity(direction)
    pressure = mc_pressure(state,layout,gamma)
    flux = 0.0_dp
    flux(layout%first_species:layout%last_species) = &
      state(layout%first_species:layout%last_species)*normal_velocity
    flux(layout%momentum) = state(layout%momentum)*normal_velocity
    flux(layout%momentum(direction)) = &
      flux(layout%momentum(direction)) + pressure
    flux(layout%total_energy) = &
      (state(layout%total_energy)+pressure)*normal_velocity
  end subroutine compute_mc_euler_physical_flux

  subroutine compute_mc_euler_rusanov_flux( &
      left, right, layout, gamma, direction, flux)
    real(dp), intent(in) :: left(:), right(:)
    type(mc_state_layout), intent(in) :: layout
    real(dp), intent(in) :: gamma
    integer, intent(in) :: direction
    real(dp), intent(out) :: flux(:)
    real(dp) :: left_flux(size(left)), right_flux(size(right))
    real(dp) :: left_speed, right_speed, wave_speed

    call compute_mc_euler_physical_flux( &
      left,layout,gamma,direction,left_flux)
    call compute_mc_euler_physical_flux( &
      right,layout,gamma,direction,right_flux)
    left_speed = abs(left(layout%momentum(direction)) / &
      mc_mixture_density(left,layout)) + mc_sound_speed(left,layout,gamma)
    right_speed = abs(right(layout%momentum(direction)) / &
      mc_mixture_density(right,layout)) + mc_sound_speed(right,layout,gamma)
    wave_speed = max(left_speed,right_speed)
    flux = 0.5_dp*(left_flux+right_flux) - &
      0.5_dp*wave_speed*(right-left)
  end subroutine compute_mc_euler_rusanov_flux

  subroutine compute_mc_euler_rhs(q,rhs,layout,config,workspace)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_primitive_workspace), intent(inout), optional, target :: workspace
    type(mc_primitive_workspace), target :: local_workspace
    type(mc_primitive_workspace), pointer :: work

    work => local_workspace
    if (present(workspace)) work => workspace
    call prepare_mc_primitives(q,layout,config,work)
    call compute_mc_euler_rhs_prepared(q,rhs,layout,config,work)
  end subroutine compute_mc_euler_rhs

  ! Low-level kernel: work MUST have been prepared from this exact q.
  ! Each line retains its previous face, including the periodic seam.
  ! Thus every face is evaluated once without full-volume flux buffers.
  subroutine compute_mc_euler_rhs_prepared(q,rhs,layout,config,work)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_primitive_workspace), intent(in) :: work
    integer :: direction, a, b, ia, ib, n, extent(3), left(3), right(3)
    real(dp) :: spacing(3), previous(layout%nvariables)
    real(dp) :: current(layout%nvariables), first(layout%nvariables)
    real(dp) :: ghost(layout%nvariables), rho, vel(3), temp, p, sound
    real(dp) :: fractions(layout%nspecies)
    logical :: periodic

    if (any(shape(rhs) /= shape(q))) &
      error stop 'multicomponent Euler RHS allocation does not match state'
    if(config%geometry /= 'cartesian') then
      call mc_mapped_rhs(q,rhs,layout,config,.false.)
      return
    end if
    extent = [config%nx,config%ny,config%nz]
    spacing = [config%x_max-config%x_min,config%y_max-config%y_min, &
      config%z_max-config%z_min]/real(extent,dp)
    rhs = 0.0_dp
    do direction=1,3
      a = mod(direction,3)+1
      b = mod(direction+1,3)+1
      periodic = mc_boundary_is_periodic(config,2*direction-1)
      !$omp parallel do collapse(2) default(shared) schedule(static) &
      !$omp private(ib,ia,n,left,right,previous,current,first,ghost,rho,vel,temp,p,sound,fractions)
      do ib=1,extent(b)
        do ia=1,extent(a)
          right = 1
          right(a)=ia
          right(b)=ib
          left=right
          if (periodic) then
            left(direction)=extent(direction)
            call interior_flux(left,right,direction,previous)
          else
            call mc_boundary_state(q(right(1),right(2),right(3),:), &
              ghost,layout,config,2*direction-1)
            call evaluate_mc_primitive(ghost,layout,config%gamma,rho,vel, &
              temp,fractions,p,sound)
            call cached_rusanov(ghost,q(right(1),right(2),right(3),:), &
              vel(direction),work%velocity(right(1),right(2),right(3),direction), &
              p,work%pressure(right(1),right(2),right(3)), &
              sound,work%sound_speed(right(1),right(2),right(3)), &
              layout,direction,previous)
          end if
          first=previous
          do n=1,extent(direction)
            left=right
            left(direction)=n
            right=left
            if (n < extent(direction)) then
              right(direction)=n+1
              call interior_flux(left,right,direction,current)
            else if (periodic) then
              current=first
            else
              call mc_boundary_state(q(left(1),left(2),left(3),:), &
                ghost,layout,config,2*direction)
              call evaluate_mc_primitive(ghost,layout,config%gamma,rho,vel, &
                temp,fractions,p,sound)
              call cached_rusanov(q(left(1),left(2),left(3),:),ghost, &
                work%velocity(left(1),left(2),left(3),direction),vel(direction), &
                work%pressure(left(1),left(2),left(3)),p, &
                work%sound_speed(left(1),left(2),left(3)),sound, &
                layout,direction,current)
            end if
            rhs(left(1),left(2),left(3),:) = rhs(left(1),left(2),left(3),:) - &
              (current-previous)/spacing(direction)
            previous=current
          end do
        end do
      end do
      !$omp end parallel do
    end do
  contains
    subroutine interior_flux(l,r,d,flux)
      integer, intent(in) :: l(3),r(3),d
      real(dp), intent(out) :: flux(:)
      call cached_rusanov(q(l(1),l(2),l(3),:),q(r(1),r(2),r(3),:), &
        work%velocity(l(1),l(2),l(3),d),work%velocity(r(1),r(2),r(3),d), &
        work%pressure(l(1),l(2),l(3)),work%pressure(r(1),r(2),r(3)), &
        work%sound_speed(l(1),l(2),l(3)),work%sound_speed(r(1),r(2),r(3)), &
        layout,d,flux)
    end subroutine interior_flux
  end subroutine compute_mc_euler_rhs_prepared

  subroutine cached_rusanov(left,right,ul,ur,pl,pr,cl,cr,layout,d,flux)
    real(dp), intent(in) :: left(:),right(:),ul,ur,pl,pr,cl,cr
    type(mc_state_layout), intent(in) :: layout
    integer, intent(in) :: d
    real(dp), intent(out) :: flux(:)
    real(dp) :: fl(size(left)),fr(size(right))
    fl=left*ul
    fr=right*ur
    fl(layout%momentum(d))=fl(layout%momentum(d))+pl
    fr(layout%momentum(d))=fr(layout%momentum(d))+pr
    fl(layout%total_energy)=(left(layout%total_energy)+pl)*ul
    fr(layout%total_energy)=(right(layout%total_energy)+pr)*ur
    flux=0.5_dp*(fl+fr)-0.5_dp*max(abs(ul)+cl,abs(ur)+cr)*(right-left)
  end subroutine cached_rusanov

  real(dp) function compute_mc_euler_timestep(q,layout,config) result(dt)
    real(dp), intent(in) :: q(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k
    real(dp) :: density, velocity(3), sound_speed, rate, maximum_rate
    real(dp) :: dx, dy, dz

    call validate_mc_euler_state(q,layout,config)
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    maximum_rate = 0.0_dp
    !$omp parallel do collapse(3) default(shared) schedule(static) &
    !$omp private(i,j,k,density,velocity,sound_speed,rate) reduction(max:maximum_rate)
    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          density = mc_mixture_density(q(i,j,k,:),layout)
          velocity = q(i,j,k,layout%momentum)/density
          sound_speed = mc_sound_speed(q(i,j,k,:),layout,config%gamma)
          rate = (abs(velocity(1))+sound_speed)/dx + &
            (abs(velocity(2))+sound_speed)/dy + &
            (abs(velocity(3))+sound_speed)/dz
          if(config%geometry /= 'cartesian') rate=mc_mapped_rate(config,i,j,k,velocity,sound_speed)
          maximum_rate = max(maximum_rate,rate)
        end do
      end do
    end do
    !$omp end parallel do
    if (config%dt > 0.0_dp) then
      if (config%dt*maximum_rate > 1.0_dp+100.0_dp*epsilon(1.0_dp)) then
        error stop 'fixed multicomponent Euler dt violates the CFL limit'
      end if
      dt = config%dt
    else
      dt = config%cfl/maximum_rate
    end if
  end function compute_mc_euler_timestep

  subroutine advance_mc_euler_ssprk3(q,q0,rhs,dt,layout,config,workspace)
    real(dp), intent(inout) :: q(:,:,:,:), q0(:,:,:,:), rhs(:,:,:,:)
    real(dp), intent(in) :: dt
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    type(mc_primitive_workspace), intent(inout), optional, target :: workspace
    type(mc_primitive_workspace), target :: local_workspace
    type(mc_primitive_workspace), pointer :: work

    if (any(shape(q0) /= shape(q)) .or. any(shape(rhs) /= shape(q))) then
      error stop 'multicomponent Euler SSPRK3 work arrays do not match state'
    end if
    work => local_workspace
    if (present(workspace)) work => workspace
    q0 = q
    call compute_mc_euler_rhs(q,rhs,layout,config,work)
    q = q0 + dt*rhs
    call compute_mc_euler_rhs(q,rhs,layout,config,work)
    q = 0.75_dp*q0 + 0.25_dp*(q+dt*rhs)
    call compute_mc_euler_rhs(q,rhs,layout,config,work)
    q = (1.0_dp/3.0_dp)*q0 + (2.0_dp/3.0_dp)*(q+dt*rhs)
    call validate_mc_euler_state(q,layout,config)
  end subroutine advance_mc_euler_ssprk3

end module mod_mc_euler_flux
