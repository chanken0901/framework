program test_viscous_central6
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_viscous_scheme, only : add_viscous_rhs, &
    validate_viscous_scheme, viscous_dt_limit
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:), rhs(:,:,:,:)
  real(dp) :: pi, x, rho, u, pressure
  real(dp) :: expected_momentum, expected_energy
  real(dp) :: momentum_error, energy_error, dt_limit, expected_dt
  integer :: i, j, k, g

  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 32
  sim%ny = 8
  sim%nz = 8
  sim%nghost = 3
  pi = acos(-1.0_dp)
  sim%dx = 2.0_dp*pi/real(sim%nx,dp)
  sim%dy = 2.0_dp*pi/real(sim%ny,dp)
  sim%dz = 2.0_dp*pi/real(sim%nz,dp)
  nse%viscous_scheme = 'central6'
  nse%reynolds = 50.0_dp
  nse%prandtl = 0.72_dp
  call validate_viscous_scheme(nse)

  g = sim%nghost
  allocate(q(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,nse%nv))
  allocate(rhs(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,nse%nv))
  q = 0.0_dp
  rhs = 0.0_dp
  rho = 1.0_dp
  pressure = 1.0_dp
  do k = 1, sim%nz
    do j = 1, sim%ny
      do i = 1, sim%nx
        x = real(i-1,dp)*sim%dx
        u = sin(x)
        q(i,j,k,1) = rho
        q(i,j,k,2) = rho*u
        q(i,j,k,3) = 0.0_dp
        q(i,j,k,4) = 0.0_dp
        q(i,j,k,5) = pressure/(nse%gamma-1.0_dp) + 0.5_dp*rho*u*u
      end do
    end do
  end do
  call apply_periodic(q, sim)

  !$OMP PARALLEL DEFAULT(SHARED)
  call add_viscous_rhs(q, rhs, sim, nse, 1, sim%ny, 1, sim%nz)
  !$OMP END PARALLEL

  momentum_error = 0.0_dp
  energy_error = 0.0_dp
  do k = 1, sim%nz
    do j = 1, sim%ny
      do i = 1, sim%nx
        x = real(i-1,dp)*sim%dx
        expected_momentum = -(4.0_dp/(3.0_dp*nse%reynolds))*sin(x)
        expected_energy = (4.0_dp/(3.0_dp*nse%reynolds))*cos(2.0_dp*x)
        momentum_error = max(momentum_error, &
          abs(rhs(i,j,k,2)-expected_momentum))
        energy_error = max(energy_error, abs(rhs(i,j,k,5)-expected_energy))
      end do
    end do
  end do

  if (momentum_error > 2.0e-7_dp .or. energy_error > 4.0e-7_dp) then
    write(*,'(A,2ES16.8)') 'central6 analytic errors: ', &
      momentum_error, energy_error
    error stop 'sixth-order viscous analytic test failed'
  end if
  if (maxval(abs(rhs(1:sim%nx,1:sim%ny,1:sim%nz,1))) > 1.0e-14_dp .or. &
      maxval(abs(rhs(1:sim%nx,1:sim%ny,1:sim%nz,3:4))) > 1.0e-14_dp) then
    error stop 'viscous operator changed mass or transverse momentum'
  end if

  call viscous_dt_limit(q, dt_limit, sim, nse, 1, sim%ny, 1, sim%nz)
  expected_dt = 2.0_dp / ((272.0_dp/45.0_dp) * &
    max(4.0_dp/3.0_dp,nse%gamma/nse%prandtl) / nse%reynolds * &
    (1.0_dp/sim%dx**2 + 1.0_dp/sim%dy**2 + 1.0_dp/sim%dz**2))
  if (abs(dt_limit-expected_dt) > 1.0e-14_dp) then
    error stop 'viscous time-step limit test failed'
  end if

  deallocate(q, rhs)
  write(*,'(A,2ES16.8)') 'central6 viscous test passed; errors = ', &
    momentum_error, energy_error

contains

  subroutine apply_periodic(field, config)
    type(simulation_config), intent(in) :: config
    real(dp), intent(inout) :: field(1-config%nghost:, &
      1-config%nghost:, 1-config%nghost:, :)
    integer :: ii, jj, kk, wi, wj, wk

    do kk = 1-config%nghost, config%nz+config%nghost
      wk = 1 + modulo(kk-1, config%nz)
      do jj = 1-config%nghost, config%ny+config%nghost
        wj = 1 + modulo(jj-1, config%ny)
        do ii = 1-config%nghost, config%nx+config%nghost
          if (ii >= 1 .and. ii <= config%nx .and. &
              jj >= 1 .and. jj <= config%ny .and. &
              kk >= 1 .and. kk <= config%nz) cycle
          wi = 1 + modulo(ii-1, config%nx)
          field(ii,jj,kk,:) = field(wi,wj,wk,:)
        end do
      end do
    end do
  end subroutine apply_periodic

end program test_viscous_central6
