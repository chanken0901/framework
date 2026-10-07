program test_keep6_accuracy
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_grid_axis, only : build_axis
  use mod_grid_fvm, only : axis_x
  use mod_convective_keep, only : keep6_inverse_metric
  use mod_convective_scheme, only : compute_convective_flux, &
    convective_required_ghost_cells
  implicit none

  integer, parameter :: resolution(3) = [24, 48, 96]
  real(dp) :: error2(3), error6(3), conservation2(3), conservation6(3)
  real(dp) :: rate2(2), rate6(2)
  integer :: level
  character(len=32) :: mode
  logical :: mapped
  call get_command_argument(1,mode)
  mapped=mode=='mapped'

  if (convective_required_ghost_cells() /= 3) then
    error stop 'sixth-order KEEP must require three ghost cells'
  end if

  do level = 1, size(resolution)
    if(.not.mapped) then
      call measure_error(resolution(level), 2, error2(level), conservation2(level))
    else
      error2(level)=1._dp;conservation2(level)=0 ! KEEP2 is covered by the separate uniform test.
    end if
    call measure_error(resolution(level), 6, error6(level), &
      conservation6(level))
  end do
  rate2(1) = log(error2(1)/error2(2))/log(2.0_dp)
  rate2(2) = log(error2(2)/error2(3))/log(2.0_dp)
  rate6(1) = log(error6(1)/error6(2))/log(2.0_dp)
  rate6(2) = log(error6(2)/error6(3))/log(2.0_dp)

  if (.not.mapped .and. minval(rate2) < 1.8_dp) then
    write(*,'(A,3ES16.8)') 'KEEP2 errors: ', error2
    write(*,'(A,2F10.5)') 'KEEP2 rates: ', rate2
    error stop 'KEEP convective derivative did not attain second order'
  end if
  if (minval(rate6) < 5.5_dp) then
    write(*,'(A,3ES16.8)') 'KEEP6 errors: ', error6
    write(*,'(A,2F10.5)') 'KEEP6 rates: ', rate6
    error stop 'KEEP convective derivative did not attain sixth order'
  end if
  if (max(maxval(conservation2),maxval(conservation6)) > 2.0e-12_dp) then
    write(*,'(A,3ES16.8)') 'KEEP2 conservation errors: ', conservation2
    write(*,'(A,3ES16.8)') 'KEEP6 conservation errors: ', conservation6
    error stop 'KEEP flux difference is not conservative'
  end if

  if(mapped) then
    write(*,'(A,2F10.5,A,ES12.4)') 'Mapped KEEP6 rates = ',rate6, &
      ', metric-weighted conservation error = ',maxval(conservation6)
  else
  write(*,'(A,2F10.5,A,2F10.5,A,ES12.4)') &
    'KEEP selectable-order accuracy test passed; KEEP2 rates = ', rate2, &
    ', KEEP6 rates = ', rate6, ', conservation error = ', &
    max(maxval(conservation2),maxval(conservation6))
  end if

contains

  subroutine measure_error(nx, keep_order, l2_error, conservation_error)
    integer, intent(in) :: nx, keep_order
    real(dp), intent(out) :: l2_error, conservation_error
    type(simulation_config) :: sim
    type(nse_config) :: nse
    real(dp), allocatable :: q(:,:,:,:), fface(:,:,:,:)
    real(dp) :: pi, x, rho, u, v, w, pressure
    real(dp) :: drho, du, dv, dw, dpressure
    real(dp) :: energy, denergy, speed_squared
    real(dp) :: exact(5), divergence(5), residual_sum(5)
    real(dp) :: error_sum,metric,s
    real(dp), allocatable :: edges(:)
    integer :: i, j, k, nghost

    call init_simulation_config(sim)
    call init_nse_config(nse)
    nghost = 3
    pi = acos(-1.0_dp)
    sim%nx = nx
    sim%ny = 1
    sim%nz = 1
    sim%nghost = nghost
    sim%dx = 2.0_dp*pi/real(nx,dp)
    sim%dy = 1.0_dp
    sim%dz = 1.0_dp
    if(mapped) then
      sim%grid_mapping='sinh'
      allocate(edges(0:nx))
      do i=0,nx
        s=2*pi*real(i,dp)/nx
        edges(i)=s+.2_dp*sin(s)
      end do
      call build_axis(edges,nghost,.true.,axis_x)
    end if
    if (keep_order == 2) then
      nse%convective_scheme = 'keep2'
    else
      nse%convective_scheme = 'keep6'
    end if

    allocate(q(1-nghost:nx+nghost,1-nghost:1+nghost, &
      1-nghost:1+nghost,5))
    allocate(fface(0:nx,0:1,0:1,5))

    do k = 1-nghost, 1+nghost
      do j = 1-nghost, 1+nghost
        do i = 1-nghost, nx+nghost
          x = (real(i,dp)-0.5_dp)*sim%dx
          if(mapped) x=axis_x%center(i)
          call analytic_state(x, nse%gamma, q(i,j,k,1), q(i,j,k,2), &
            q(i,j,k,3), q(i,j,k,4), q(i,j,k,5))
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
      if(mapped) x=axis_x%center(i)
      rho = 1.0_dp + 0.1_dp*sin(x)
      u = 0.3_dp + 0.05_dp*cos(2.0_dp*x)
      v = 0.1_dp*sin(3.0_dp*x)
      w = -0.05_dp*cos(x)
      pressure = 1.0_dp + 0.08_dp*cos(x)
      drho = 0.1_dp*cos(x)
      du = -0.1_dp*sin(2.0_dp*x)
      dv = 0.3_dp*cos(3.0_dp*x)
      dw = 0.05_dp*sin(x)
      dpressure = -0.08_dp*sin(x)
      speed_squared = u*u + v*v + w*w
      energy = pressure/(nse%gamma-1.0_dp) + &
        0.5_dp*rho*speed_squared
      denergy = dpressure/(nse%gamma-1.0_dp) + &
        0.5_dp*drho*speed_squared + rho*(u*du+v*dv+w*dw)

      exact(1) = drho*u + rho*du
      exact(2) = drho*u*u + 2.0_dp*rho*u*du + dpressure
      exact(3) = drho*u*v + rho*du*v + rho*u*dv
      exact(4) = drho*u*w + rho*du*w + rho*u*dw
      exact(5) = (denergy+dpressure)*u + (energy+pressure)*du
      metric=sim%dx
      if(mapped) metric=1._dp/keep6_inverse_metric(i,1,1,1)
      divergence = (fface(i,1,1,:)-fface(i-1,1,1,:))/metric
      error_sum = error_sum + sum((divergence-exact)**2)
      residual_sum = residual_sum + divergence*metric
    end do
    l2_error = sqrt(error_sum/real(5*nx,dp))
    conservation_error = maxval(abs(residual_sum))

    deallocate(q, fface)
  end subroutine measure_error

  pure subroutine analytic_state(x, gamma, density, momentum_x, &
      momentum_y, momentum_z, total_energy)
    real(dp), intent(in) :: x, gamma
    real(dp), intent(out) :: density, momentum_x, momentum_y
    real(dp), intent(out) :: momentum_z, total_energy
    real(dp) :: rho, u, v, w, pressure

    rho = 1.0_dp + 0.1_dp*sin(x)
    u = 0.3_dp + 0.05_dp*cos(2.0_dp*x)
    v = 0.1_dp*sin(3.0_dp*x)
    w = -0.05_dp*cos(x)
    pressure = 1.0_dp + 0.08_dp*cos(x)
    density = rho
    momentum_x = rho*u
    momentum_y = rho*v
    momentum_z = rho*w
    total_energy = pressure/(gamma-1.0_dp) + &
      0.5_dp*rho*(u*u+v*v+w*w)
  end subroutine analytic_state

end program test_keep6_accuracy
