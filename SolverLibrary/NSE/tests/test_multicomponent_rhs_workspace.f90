program test_multicomponent_rhs_workspace
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config, read_mc_config
  use mod_mc_state_layout, only : mc_state_layout, initialize_mc_state_layout
  use mod_mc_euler_config
  use mod_mc_euler_field, only : set_mc_euler_conservative_state, validate_mc_euler_state
  use mod_mc_boundary, only : mc_boundary_state, mc_boundary_is_periodic
  use mod_mc_euler_flux, only : compute_mc_euler_rhs, compute_mc_euler_rusanov_flux
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_rhs, &
    compute_mc_transport_rhs, mc_navier_stokes_workspace
  use mod_mc_thermodynamics_provider, only : configure_mc_thermodynamics, mc_mixture_gas_constant
  use mod_mc_transport_provider, only : configure_mc_transport
  implicit none
  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: config
  type(mc_navier_stokes_workspace) :: workspace
  real(dp), allocatable :: q(:,:,:,:), actual(:,:,:,:), expected(:,:,:,:), transport(:,:,:,:)
  real(dp) :: fractions(3), velocity(3), rho, temperature, pressure, phase
  real(dp) :: reference_seconds,cached_seconds
  integer :: grid,boundary,revision,i,j,k,repeat
  integer(kind=8) :: start,finish,clock_rate
  character(len=512) :: path
  character(len=16), parameter :: kinds(4) = &
    [character(len=16) :: 'periodic','reflective','dirichlet','non_reflecting']

  call get_command_argument(1,path)
  call read_mc_config(trim(path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  if (layout%nspecies /= 3) error stop 'workspace fixture requires three species'
  call configure_mc_thermodynamics(trim(path),model%nspecies,model%species_names)
  call configure_mc_transport(trim(path),model%nspecies,model%species_names)
  call read_mc_euler_config(trim(path),model%nspecies,config)
  config%boundary_reference_densities=1.0_dp
  config%boundary_reference_pressures=300000.0_dp
  config%boundary_reference_velocities=0.0_dp
  config%boundary_reference_mass_fractions=0.0_dp
  config%boundary_reference_mass_fractions(:,1)=0.3_dp
  config%boundary_reference_mass_fractions(:,2)=0.4_dp
  config%boundary_reference_mass_fractions(:,3)=0.3_dp

  ! Exercise reuse after state changes and resizing, including single-cell axes.
  do grid=1,2
    if (grid == 1) then
      config%nx=5
      config%ny=4
      config%nz=3
    else
      config%nx=1
      config%ny=2
      config%nz=1
    end if
    allocate(q(config%nx,config%ny,config%nz,layout%nvariables))
    allocate(actual,mold=q)
    allocate(expected,mold=q)
    allocate(transport,mold=q)
    do boundary=1,5
      config%boundary_condition='face_specific'
      if (boundary <= 4) then
        config%boundary_face_types=kinds(min(boundary,4))
      else
        config%boundary_face_types = [character(len=32) :: &
          'reflective','non_reflecting','dirichlet','reflective','non_reflecting','dirichlet']
      end if
      do revision=1,2
        do k=1,config%nz
          do j=1,config%ny
            do i=1,config%nx
              phase=real(i+2*j+3*k+revision,dp)
              rho=1.0_dp+0.05_dp*sin(phase)
              fractions=[0.3_dp+0.02_dp*sin(phase),0.4_dp,0.3_dp-0.02_dp*sin(phase)]
              temperature=900.0_dp+200.0_dp*sin(phase)
              velocity=[2.0_dp*sin(phase),3.0_dp*cos(phase),sin(2.0_dp*phase)]
              pressure=rho*mc_mixture_gas_constant(fractions,layout)*temperature
              call set_mc_euler_conservative_state(q(i,j,k,:),layout,config%gamma, &
                rho,velocity,pressure,fractions)
            end do
          end do
        end do
        call reference_euler_rhs(q,expected,layout,config)
        call compute_mc_euler_rhs(q,actual,layout,config,workspace%primitive)
        call assert_fields(actual,expected,'cached Euler / original six-face kernel')
        call compute_mc_transport_rhs(q,transport,layout,config)
        expected=expected+transport
        call compute_mc_navier_stokes_rhs(q,actual,layout,config,workspace)
        call assert_fields(actual,expected,'shared cache / direct accumulation')
        call compute_mc_navier_stokes_rhs(q,actual,layout,config)
        call assert_fields(actual,expected,'legacy API fallback')
      end do
    end do
    if (grid == 1) then
      config%boundary_face_types='periodic'
      call system_clock(start,clock_rate)
      do repeat=1,30
        call reference_euler_rhs(q,expected,layout,config)
      end do
      call system_clock(finish)
      reference_seconds=real(finish-start,dp)/real(clock_rate,dp)
      call system_clock(start)
      do repeat=1,30
        call compute_mc_euler_rhs(q,actual,layout,config,workspace%primitive)
      end do
      call system_clock(finish)
      cached_seconds=real(finish-start,dp)/real(clock_rate,dp)
      write(*,'(A,2ES14.5)') 'Euler reference / cached seconds: ',reference_seconds,cached_seconds
    end if
    deallocate(q,actual,expected,transport)
  end do
  write(*,'(A)') 'Multicomponent RHS workspace regression tests passed.'
contains
  subroutine assert_fields(actual,expected,message)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp), intent(in) :: actual(:,:,:,:),expected(:,:,:,:)
    character(len=*), intent(in) :: message
    integer :: v
    if (.not. all(ieee_is_finite(actual))) error stop 'non-finite workspace result'
    do v=1,size(actual,4)
      if (maxval(abs(actual(:,:,:,v)-expected(:,:,:,v))) > &
          2.0e-11_dp*max(1.0_dp,maxval(abs(expected(:,:,:,v))))) then
        write(*,*) trim(message),v
        error stop 'workspace regression mismatch'
      end if
    end do
  end subroutine assert_fields

  ! Frozen pre-optimization reference, intentionally retains duplicate face work.
  subroutine reference_euler_rhs(q, rhs, layout, config)
    real(dp), intent(in) :: q(:,:,:,:)
    real(dp), intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout), intent(in) :: layout
    type(mc_euler_config), intent(in) :: config
    integer :: i, j, k, im, ip, jm, jp, km, kp
    real(dp) :: dx, dy, dz
    real(dp) :: positive_flux(layout%nvariables)
    real(dp) :: negative_flux(layout%nvariables)
    real(dp) :: ghost_state(layout%nvariables)

    if (any(shape(rhs) /= shape(q))) then
      error stop 'multicomponent Euler RHS allocation does not match state'
    end if
    call validate_mc_euler_state(q,layout,config)
    dx = (config%x_max-config%x_min)/real(config%nx,dp)
    dy = (config%y_max-config%y_min)/real(config%ny,dp)
    dz = (config%z_max-config%z_min)/real(config%nz,dp)
    rhs = 0.0_dp
    do k = 1, config%nz
      km = merge(config%nz,k-1,k == 1)
      kp = merge(1,k+1,k == config%nz)
      do j = 1, config%ny
        jm = merge(config%ny,j-1,j == 1)
        jp = merge(1,j+1,j == config%ny)
        do i = 1, config%nx
          im = merge(config%nx,i-1,i == 1)
          ip = merge(1,i+1,i == config%nx)
          if (i == config%nx .and. &
              .not. mc_boundary_is_periodic(config,mc_face_x_max)) then
            call mc_boundary_state( &
              q(i,j,k,:),ghost_state,layout,config,mc_face_x_max)
            call compute_mc_euler_rusanov_flux( &
              q(i,j,k,:),ghost_state,layout,config%gamma,1,positive_flux)
          else
            call compute_mc_euler_rusanov_flux( &
              q(i,j,k,:),q(ip,j,k,:),layout,config%gamma,1,positive_flux)
          end if
          if (i == 1 .and. &
              .not. mc_boundary_is_periodic(config,mc_face_x_min)) then
            call mc_boundary_state( &
              q(i,j,k,:),ghost_state,layout,config,mc_face_x_min)
            call compute_mc_euler_rusanov_flux( &
              ghost_state,q(i,j,k,:),layout,config%gamma,1,negative_flux)
          else
            call compute_mc_euler_rusanov_flux( &
              q(im,j,k,:),q(i,j,k,:),layout,config%gamma,1,negative_flux)
          end if
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (positive_flux-negative_flux)/dx
          if (j == config%ny .and. &
              .not. mc_boundary_is_periodic(config,mc_face_y_max)) then
            call mc_boundary_state( &
              q(i,j,k,:),ghost_state,layout,config,mc_face_y_max)
            call compute_mc_euler_rusanov_flux( &
              q(i,j,k,:),ghost_state,layout,config%gamma,2,positive_flux)
          else
            call compute_mc_euler_rusanov_flux( &
              q(i,j,k,:),q(i,jp,k,:),layout,config%gamma,2,positive_flux)
          end if
          if (j == 1 .and. &
              .not. mc_boundary_is_periodic(config,mc_face_y_min)) then
            call mc_boundary_state( &
              q(i,j,k,:),ghost_state,layout,config,mc_face_y_min)
            call compute_mc_euler_rusanov_flux( &
              ghost_state,q(i,j,k,:),layout,config%gamma,2,negative_flux)
          else
            call compute_mc_euler_rusanov_flux( &
              q(i,jm,k,:),q(i,j,k,:),layout,config%gamma,2,negative_flux)
          end if
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (positive_flux-negative_flux)/dy
          if (k == config%nz .and. &
              .not. mc_boundary_is_periodic(config,mc_face_z_max)) then
            call mc_boundary_state( &
              q(i,j,k,:),ghost_state,layout,config,mc_face_z_max)
            call compute_mc_euler_rusanov_flux( &
              q(i,j,k,:),ghost_state,layout,config%gamma,3,positive_flux)
          else
            call compute_mc_euler_rusanov_flux( &
              q(i,j,k,:),q(i,j,kp,:),layout,config%gamma,3,positive_flux)
          end if
          if (k == 1 .and. &
              .not. mc_boundary_is_periodic(config,mc_face_z_min)) then
            call mc_boundary_state( &
              q(i,j,k,:),ghost_state,layout,config,mc_face_z_min)
            call compute_mc_euler_rusanov_flux( &
              ghost_state,q(i,j,k,:),layout,config%gamma,3,negative_flux)
          else
            call compute_mc_euler_rusanov_flux( &
              q(i,j,km,:),q(i,j,k,:),layout,config%gamma,3,negative_flux)
          end if
          rhs(i,j,k,:) = rhs(i,j,k,:) - &
            (positive_flux-negative_flux)/dz
        end do
      end do
    end do
  end subroutine reference_euler_rhs
end program test_multicomponent_rhs_workspace
