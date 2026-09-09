program test_weno5z_roe
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_reconstruction_weno5z, only : reconstruct_weno5z_left, &
    reconstruct_weno5z_right
  use mod_riemann_roe, only : rotate_conserved_to_normal, &
    rotate_flux_to_global, roe_eigensystem, roe_numerical_flux, &
    euler_physical_flux
  use mod_convective_scheme, only : compute_convective_flux, &
    validate_convective_scheme, convective_required_ghost_cells
  use mod_convective_weno5z_roe, only : limit_weno_state
  implicit none

  call test_scalar_reconstruction_order()
  call test_roe_eigensystem_and_uniform_flux()
  call test_directional_flux_rotation()
  call test_entropy_wave_order_and_conservation()
  call test_sod_discontinuity_is_finite()
  call test_positive_reconstruction()

  write(*,'(A)') 'WENO5-Z/Roe tests passed'

contains

  subroutine test_positive_reconstruction()
    type(nse_config) :: nse
    real(dp) :: center(5), trial(5), original(5), p
    logical :: limited
    integer :: mode
    call init_nse_config(nse)
    center=[1.0_dp,0.3_dp,0.1_dp,0.0_dp,2.55_dp]
    trial=center
    call limit_weno_state(center,trial,nse,limited)
    if (limited .or. any(trial/=center)) error stop 'limiter changed admissible state'
    do mode=1,2
      trial=center
      if (mode==1) trial(1)=-1.0_dp
      if (mode==2) trial(5)=0.0_dp
      original=trial
      call limit_weno_state(center,trial,nse,limited)
      if (.not. limited) error stop 'invalid reconstruction was not limited'
      if (trial(1)<nse%small_rho) error stop 'limiter left negative density'
      p=(nse%gamma-1)*(trial(5)-0.5_dp*sum(trial(2:4)**2)/trial(1))
      if (p<nse%small_p) error stop 'limiter left negative pressure'
      if (any(trial<min(center,original)) .or. any(trial>max(center,original))) &
        error stop 'limiter is not a convex reconstruction'
    end do
  end subroutine

  subroutine test_scalar_reconstruction_order()
    real(dp), parameter :: spacing(3) = [0.2_dp, 0.1_dp, 0.05_dp]
    real(dp), parameter :: left_offset(5) = &
      [-2.5_dp, -1.5_dp, -0.5_dp, 0.5_dp, 1.5_dp]
    real(dp), parameter :: right_offset(5) = &
      [-1.5_dp, -0.5_dp, 0.5_dp, 1.5_dp, 2.5_dp]
    real(dp) :: left_value(5), right_value(5), error(3), rate(2)
    integer :: level, point

    do level = 1, size(spacing)
      do point = 1, 5
        left_value(point) = exponential_cell_average( &
          left_offset(point)*spacing(level), spacing(level))
        right_value(point) = exponential_cell_average( &
          right_offset(point)*spacing(level), spacing(level))
      end do
      error(level) = max( &
        abs(reconstruct_weno5z_left(left_value)-1.0_dp), &
        abs(reconstruct_weno5z_right(right_value)-1.0_dp))
    end do
    rate(1) = log(error(1)/error(2))/log(2.0_dp)
    rate(2) = log(error(2)/error(3))/log(2.0_dp)
    if (minval(rate) < 4.5_dp) then
      write(*,'(A,3ES16.8)') 'WENO5-Z scalar errors: ', error
      write(*,'(A,2F10.5)') 'WENO5-Z scalar rates: ', rate
      error stop 'WENO5-Z scalar reconstruction did not attain fifth order'
    end if
  end subroutine test_scalar_reconstruction_order

  pure real(dp) function exponential_cell_average(center, spacing) &
      result(value)
    real(dp), intent(in) :: center, spacing

    value = (exp(center+0.5_dp*spacing) - &
      exp(center-0.5_dp*spacing)) / spacing
  end function exponential_cell_average

  subroutine test_roe_eigensystem_and_uniform_flux()
    type(nse_config) :: nse
    real(dp) :: state(5), numerical_flux(5), physical_flux(5)
    real(dp) :: right_matrix(5,5), left_matrix(5,5), eigenvalue(5)
    real(dp) :: identity(5,5), expected_identity(5,5)
    integer :: component

    call init_nse_config(nse)
    call conserved_state(1.2_dp, 0.4_dp, -0.2_dp, 0.1_dp, 1.0_dp, &
      nse%gamma, state)
    call roe_eigensystem(state, state, nse, right_matrix, left_matrix, &
      eigenvalue)
    identity = matmul(left_matrix, right_matrix)
    expected_identity = 0.0_dp
    do component = 1, 5
      expected_identity(component,component) = 1.0_dp
    end do
    if (maxval(abs(identity-expected_identity)) > 2.0e-12_dp) then
      error stop 'Roe left and right eigenvectors are not inverses'
    end if

    call roe_numerical_flux(state, state, nse, numerical_flux)
    call euler_physical_flux(state, nse, physical_flux)
    if (maxval(abs(numerical_flux-physical_flux)) > 2.0e-12_dp) then
      error stop 'Roe flux is inconsistent for equal left/right states'
    end if
  end subroutine test_roe_eigensystem_and_uniform_flux

  subroutine test_directional_flux_rotation()
    type(nse_config) :: nse
    real(dp) :: global_state(5), normal_state(5), normal_flux(5)
    real(dp) :: global_flux(5), expected(5)
    real(dp) :: rho, u, v, w, pressure, energy, normal_velocity
    integer :: direction

    call init_nse_config(nse)
    rho = 1.3_dp
    u = 0.4_dp
    v = -0.25_dp
    w = 0.15_dp
    pressure = 0.9_dp
    call conserved_state(rho, u, v, w, pressure, nse%gamma, global_state)
    energy = global_state(5)

    do direction = 1, 3
      call rotate_conserved_to_normal(global_state, direction, normal_state)
      call euler_physical_flux(normal_state, nse, normal_flux)
      call rotate_flux_to_global(normal_flux, direction, global_flux)
      select case (direction)
      case (1)
        normal_velocity = u
      case (2)
        normal_velocity = v
      case (3)
        normal_velocity = w
      end select
      expected = [rho*normal_velocity, rho*normal_velocity*u, &
        rho*normal_velocity*v, rho*normal_velocity*w, &
        normal_velocity*(energy+pressure)]
      expected(1+direction) = expected(1+direction) + pressure
      if (maxval(abs(global_flux-expected)) > 2.0e-12_dp) then
        error stop 'normal-to-global Roe flux rotation is incorrect'
      end if
    end do
  end subroutine test_directional_flux_rotation

  subroutine test_entropy_wave_order_and_conservation()
    integer, parameter :: resolution(3) = [24, 48, 96]
    real(dp) :: error(3), conservation(3), rate(2)
    integer :: level

    if (convective_required_ghost_cells() /= 3) then
      error stop 'WENO5-Z/Roe must require three ghost cells'
    end if
    do level = 1, size(resolution)
      call entropy_wave_error(resolution(level), error(level), &
        conservation(level))
    end do
    rate(1) = log(error(1)/error(2))/log(2.0_dp)
    rate(2) = log(error(2)/error(3))/log(2.0_dp)
    if (minval(rate) < 4.3_dp) then
      write(*,'(A,3ES16.8)') 'WENO5-Z/Roe entropy-wave errors: ', error
      write(*,'(A,2F10.5)') 'WENO5-Z/Roe entropy-wave rates: ', rate
      error stop 'WENO5-Z/Roe flux divergence did not attain fifth order'
    end if
    if (maxval(conservation) > 3.0e-12_dp) then
      write(*,'(A,3ES16.8)') 'WENO5-Z/Roe conservation errors: ', &
        conservation
      error stop 'WENO5-Z/Roe face flux difference is not conservative'
    end if
  end subroutine test_entropy_wave_order_and_conservation

  subroutine entropy_wave_error(nx, l2_error, conservation_error)
    integer, intent(in) :: nx
    real(dp), intent(out) :: l2_error, conservation_error
    type(simulation_config) :: sim
    type(nse_config) :: nse
    real(dp), allocatable :: q(:,:,:,:), fface(:,:,:,:)
    real(dp) :: pi, x, density, density_derivative, average_factor
    real(dp) :: u, v, w, pressure, speed_squared, exact(5)
    real(dp) :: divergence(5), residual_sum(5), error_sum
    integer :: i, j, k, nghost

    call init_simulation_config(sim)
    call init_nse_config(nse)
    nse%convective_scheme = 'weno5z_roe'
    call validate_convective_scheme(nse)
    nghost = 3
    pi = acos(-1.0_dp)
    sim%nx = nx
    sim%ny = 1
    sim%nz = 1
    sim%nghost = nghost
    sim%dx = 2.0_dp*pi/real(nx,dp)
    sim%dy = 1.0_dp
    sim%dz = 1.0_dp
    u = 0.35_dp
    v = -0.12_dp
    w = 0.08_dp
    pressure = 1.0_dp
    speed_squared = u*u + v*v + w*w
    average_factor = sin(0.5_dp*sim%dx)/(0.5_dp*sim%dx)

    allocate(q(1-nghost:nx+nghost,1-nghost:1+nghost, &
      1-nghost:1+nghost,5))
    allocate(fface(0:nx,0:1,0:1,5))
    do k = 1-nghost, 1+nghost
      do j = 1-nghost, 1+nghost
        do i = 1-nghost, nx+nghost
          x = (real(i,dp)-0.5_dp)*sim%dx
          density = 1.0_dp + 0.1_dp*average_factor*sin(x)
          call conserved_state(density, u, v, w, pressure, nse%gamma, &
            q(i,j,k,:))
        end do
      end do
    end do

    !$OMP PARALLEL DEFAULT(shared)
    call compute_convective_flux(q, fface, 1, sim, nse, 1, 1, 1, 1)
    !$OMP END PARALLEL

    error_sum = 0.0_dp
    residual_sum = 0.0_dp
    do i = 1, nx
      x = (real(i,dp)-0.5_dp)*sim%dx
      density_derivative = 0.1_dp*average_factor*cos(x)
      exact = density_derivative * [u, u*u, u*v, u*w, &
        0.5_dp*u*speed_squared]
      divergence = (fface(i,1,1,:)-fface(i-1,1,1,:))/sim%dx
      error_sum = error_sum + sum((divergence-exact)**2)
      residual_sum = residual_sum + divergence*sim%dx
    end do
    l2_error = sqrt(error_sum/real(5*nx,dp))
    conservation_error = maxval(abs(residual_sum))
    deallocate(q, fface)
  end subroutine entropy_wave_error

  subroutine test_sod_discontinuity_is_finite()
    type(simulation_config) :: sim
    type(nse_config) :: nse
    real(dp), allocatable :: q(:,:,:,:), fface(:,:,:,:)
    real(dp) :: density, pressure
    integer :: i, j, k, nghost

    call init_simulation_config(sim)
    call init_nse_config(nse)
    nse%convective_scheme = 'weno5z_roe'
    nghost = 3
    sim%nx = 32
    sim%ny = 1
    sim%nz = 1
    sim%nghost = nghost
    allocate(q(1-nghost:sim%nx+nghost,1-nghost:1+nghost, &
      1-nghost:1+nghost,5))
    allocate(fface(0:sim%nx,0:1,0:1,5))
    do k = 1-nghost, 1+nghost
      do j = 1-nghost, 1+nghost
        do i = 1-nghost, sim%nx+nghost
          if (i <= sim%nx/2) then
            density = 1.0_dp
            pressure = 1.0_dp
          else
            density = 0.125_dp
            pressure = 0.1_dp
          end if
          call conserved_state(density, 0.0_dp, 0.0_dp, 0.0_dp, &
            pressure, nse%gamma, q(i,j,k,:))
        end do
      end do
    end do

    !$OMP PARALLEL DEFAULT(shared)
    call compute_convective_flux(q, fface, 1, sim, nse, 1, 1, 1, 1)
    !$OMP END PARALLEL
    if (.not. all(ieee_is_finite(fface))) then
      error stop 'WENO5-Z/Roe produced a non-finite Sod interface flux'
    end if
    deallocate(q, fface)
  end subroutine test_sod_discontinuity_is_finite

  pure subroutine conserved_state(rho, u, v, w, pressure, gamma, state)
    real(dp), intent(in) :: rho, u, v, w, pressure, gamma
    real(dp), intent(out) :: state(5)

    state(1) = rho
    state(2) = rho*u
    state(3) = rho*v
    state(4) = rho*w
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*(u*u+v*v+w*w)
  end subroutine conserved_state

end program test_weno5z_roe
