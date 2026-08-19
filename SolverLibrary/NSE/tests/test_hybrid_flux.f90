program test_hybrid_flux
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config
  use mod_convective_scheme, only : compute_convective_flux, &
    validate_convective_scheme
  use mod_convective_hybrid, only : hybrid_face_weight
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:), hybrid_flux(:,:,:,:)
  real(dp), allocatable :: reference_flux(:,:,:,:)
  real(dp) :: flux_error, alpha, conservation(5)
  integer :: i, j, k, middle

  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 12
  sim%ny = 4
  sim%nz = 3
  sim%nghost = 3
  sim%dx = 1.0_dp/real(sim%nx,dp)
  sim%dy = 1.0_dp/real(sim%ny,dp)
  sim%dz = 1.0_dp/real(sim%nz,dp)
  nse%convective_scheme = 'hybrid'
  nse%hybrid_smooth_scheme = 'keep6'
  nse%hybrid_shock_scheme = 'weno5z_roe'
  nse%hybrid_sensor = 'ducros_pressure'
  call validate_convective_scheme(nse)

  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    1-sim%nghost:sim%ny+sim%nghost, &
    1-sim%nghost:sim%nz+sim%nghost,5))
  allocate(hybrid_flux(0:sim%nx,0:sim%ny,0:sim%nz,5))
  allocate(reference_flux(0:sim%nx,0:sim%ny,0:sim%nz,5))

  do k = 1, sim%nz
    do j = 1, sim%ny
      do i = 1, sim%nx
        call set_state(q(i,j,k,:), &
          1.0_dp+0.02_dp*sin(real(i+j+k,dp)), &
          0.1_dp, -0.04_dp, 0.02_dp, 1.0_dp, nse%gamma)
      end do
    end do
  end do
  call apply_periodic(q, sim)

  !$OMP PARALLEL DEFAULT(shared)
  call compute_convective_flux(q, hybrid_flux, 1, sim, nse, &
    1, sim%ny, 1, sim%nz)
  !$OMP END PARALLEL
  nse%convective_scheme = 'keep6'
  !$OMP PARALLEL DEFAULT(shared)
  call compute_convective_flux(q, reference_flux, 1, sim, nse, &
    1, sim%ny, 1, sim%nz)
  !$OMP END PARALLEL
  flux_error = maxval(abs(hybrid_flux-reference_flux))
  if (flux_error > 5.0e-14_dp) then
    write(*,'(A,ES24.16)') 'Smooth-limit flux error: ', flux_error
    error stop 'hybrid flux did not reduce to its smooth scheme'
  end if

  middle = sim%nx/2
  do k = 1, sim%nz
    do j = 1, sim%ny
      do i = 1, sim%nx
        if (i <= middle) then
          call set_state(q(i,j,k,:), 1.0_dp, 0.3_dp, 0.0_dp, 0.0_dp, &
            1.5_dp, nse%gamma)
        else
          call set_state(q(i,j,k,:), 1.2_dp, -0.3_dp, 0.0_dp, 0.0_dp, &
            0.7_dp, nse%gamma)
        end if
      end do
    end do
  end do
  call apply_periodic(q, sim)
  nse%convective_scheme = 'hybrid'
  nse%hybrid_sensor_onset = 0.0_dp
  nse%hybrid_sensor_full = 1.0e-8_dp
  alpha = hybrid_face_weight(q, middle, 1, 1, 1, sim, nse, 1, 1)
  if (abs(alpha-1.0_dp) > 1.0e-14_dp) then
    write(*,'(A,ES24.16)') 'Shock-face WENO weight: ', alpha
    error stop 'hybrid sensor did not select the shock scheme'
  end if

  !$OMP PARALLEL DEFAULT(shared)
  call compute_convective_flux(q, hybrid_flux, 1, sim, nse, &
    1, sim%ny, 1, sim%nz)
  !$OMP END PARALLEL
  nse%convective_scheme = 'weno5z_roe'
  !$OMP PARALLEL DEFAULT(shared)
  call compute_convective_flux(q, reference_flux, 1, sim, nse, &
    1, sim%ny, 1, sim%nz)
  !$OMP END PARALLEL
  flux_error = maxval(abs(hybrid_flux(middle,1:sim%ny,1:sim%nz,:) - &
    reference_flux(middle,1:sim%ny,1:sim%nz,:)))
  if (flux_error > 5.0e-13_dp) then
    write(*,'(A,ES24.16)') 'Shock-limit flux error: ', flux_error
    error stop 'hybrid flux did not reduce to its shock scheme'
  end if

  conservation = 0.0_dp
  do k = 1, sim%nz
    do j = 1, sim%ny
      do i = 1, sim%nx
        conservation = conservation + hybrid_flux(i,j,k,:) - &
          hybrid_flux(i-1,j,k,:)
      end do
    end do
  end do
  if (maxval(abs(conservation)) > 5.0e-13_dp) then
    write(*,'(A,5ES16.8)') 'Hybrid conservation error: ', conservation
    error stop 'hybrid face flux is not conservative'
  end if

  deallocate(q, hybrid_flux, reference_flux)
  write(*,'(A)') 'Hybrid KEEP/WENO flux tests passed'

contains

  pure subroutine set_state(state, density, u, v, w, pressure, gamma)
    real(dp), intent(out) :: state(5)
    real(dp), intent(in) :: density, u, v, w, pressure, gamma

    state(1) = density
    state(2) = density*u
    state(3) = density*v
    state(4) = density*w
    state(5) = pressure/(gamma-1.0_dp) + &
      0.5_dp*density*(u*u+v*v+w*w)
  end subroutine set_state

  subroutine apply_periodic(state, sim)
    type(simulation_config), intent(in) :: sim
    real(dp), intent(inout) :: state(1-sim%nghost:,1-sim%nghost:, &
      1-sim%nghost:,:)
    integer :: i, j, k, wrapped_i, wrapped_j, wrapped_k

    do k = 1-sim%nghost, sim%nz+sim%nghost
      wrapped_k = 1+modulo(k-1,sim%nz)
      do j = 1-sim%nghost, sim%ny+sim%nghost
        wrapped_j = 1+modulo(j-1,sim%ny)
        do i = 1-sim%nghost, sim%nx+sim%nghost
          if (i >= 1 .and. i <= sim%nx .and. j >= 1 .and. &
              j <= sim%ny .and. k >= 1 .and. k <= sim%nz) cycle
          wrapped_i = 1+modulo(i-1,sim%nx)
          state(i,j,k,:) = state(wrapped_i,wrapped_j,wrapped_k,:)
        end do
      end do
    end do
  end subroutine apply_periodic

end program test_hybrid_flux
