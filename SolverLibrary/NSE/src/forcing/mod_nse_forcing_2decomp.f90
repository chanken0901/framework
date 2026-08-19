module mod_nse_forcing
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_nse_forcing_common, only : forcing_is_enabled, &
    validate_forcing_parameters, compute_forcing_coefficients
  use module_mpi, only : my_rank, root, nprocs, ndiv_ny, ndiv_nz, &
    jjsta, jjend, kksta, kkend, mp_allsumr8, MPI_COMM_WORLD, &
    MPI_INTEGER, MPI_DOUBLE_PRECISION
  use decomp_2d, only : decomp_info, decomp_2d_init, decomp_2d_finalize, &
    alloc_x, alloc_z, xstart, xend, zstart, zend
  use decomp_2d_fft, only : decomp_2d_fft_init, decomp_2d_fft_finalize, &
    decomp_2d_fft_3d, decomp_2d_fft_get_ph
  use decomp_2d_constants, only : mytype, PHYSICAL_IN_X, &
    DECOMP_2D_FFT_FORWARD, DECOMP_2D_FFT_BACKWARD
  implicit none
  private

  logical :: initialized = .false.
  integer :: evaluation_count = 0
  complex(mytype), allocatable :: work_x(:,:,:)
  complex(mytype), allocatable :: ws_x(:,:,:), ws_y(:,:,:), ws_z(:,:,:)
  complex(mytype), allocatable :: wd_x(:,:,:), wd_y(:,:,:), wd_z(:,:,:)
  real(dp), allocatable :: local_component(:,:,:)

  public :: initialize_nse_forcing
  public :: add_nse_forcing_rhs
  public :: finalize_nse_forcing
  public :: validate_nse_forcing

contains

  subroutine validate_nse_forcing(nse)
    type(nse_config), intent(in) :: nse

    call validate_forcing_parameters(nse, '2decomp_fftw')
  end subroutine validate_nse_forcing

  subroutine initialize_nse_forcing(sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    type(decomp_info), pointer :: ph

    call validate_nse_forcing(nse)
    if (.not. forcing_is_enabled(nse)) return
    if (initialized) error stop 'NSE forcing backend is already initialized'
    if (storage_size(0.0_mytype) /= storage_size(0.0_dp)) then
      error stop '2DECOMP&FFT precision must match mod_precision dp'
    end if

    call decomp_2d_init(sim%nx, sim%ny, sim%nz, ndiv_ny, ndiv_nz, &
      complex_pool=.true.)
    call decomp_2d_fft_init(PHYSICAL_IN_X)
    ph => decomp_2d_fft_get_ph()
    if (xstart(1) /= 1 .or. xend(1) /= sim%nx) then
      error stop 'Forcing requires a 2DECOMP physical x-pencil'
    end if
    call alloc_x(work_x, ph, .true.)
    call alloc_z(ws_x, ph, .true.)
    call alloc_z(ws_y, ph, .true.)
    call alloc_z(ws_z, ph, .true.)
    call alloc_z(wd_x, ph, .true.)
    call alloc_z(wd_y, ph, .true.)
    call alloc_z(wd_z, ph, .true.)
    allocate(local_component(1:sim%nx, js:je, ks:ke))
    nullify(ph)

    initialized = .true.
    evaluation_count = 0
    if (my_rank == root) then
      write(*,'(A)') '# Petersen-Livescu forcing initialized'
      write(*,'(A)') '# forcing FFT backend: 2DECOMP&FFT + FFTW3'
      write(*,'(A,A)') '# forcing spectrum: ', trim(nse%forcing_spectrum)
      write(*,'(A,I0,A,I0)') '# forcing process grid: ', ndiv_ny, ' x ', ndiv_nz
    end if
  end subroutine initialize_nse_forcing

  subroutine add_nse_forcing_rhs(q, rhs, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)

    if (.not. forcing_is_enabled(nse)) return
    if (.not. initialized) error stop 'NSE forcing backend is not initialized'
    !$OMP BARRIER
    !$OMP MASTER
    call evaluate_and_add_forcing(q, rhs, sim, nse, js, je, ks, ke)
    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine add_nse_forcing_rhs

  subroutine finalize_nse_forcing()
    if (.not. initialized) return
    deallocate(work_x, ws_x, ws_y, ws_z, wd_x, wd_y, wd_z)
    deallocate(local_component)
    call decomp_2d_fft_finalize
    call decomp_2d_finalize
    initialized = .false.
  end subroutine finalize_nse_forcing

  subroutine evaluate_and_add_forcing(q, rhs, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: denominator_s, denominator_d, pressure_dilatation
    real(dp) :: coefficient_s, coefficient_d, target_s, target_d

    call forward_weighted_velocity(q, 2, ws_x, sim, nse, js, je, ks, ke)
    call forward_weighted_velocity(q, 3, ws_y, sim, nse, js, je, ks, ke)
    call forward_weighted_velocity(q, 4, ws_z, sim, nse, js, je, ks, ke)
    call project_helmholtz(sim, nse, denominator_s, denominator_d)
    pressure_dilatation = compute_pressure_dilatation(q, sim, nse, &
      js, je, ks, ke)
    call compute_forcing_coefficients(nse, denominator_s, denominator_d, &
      pressure_dilatation, coefficient_s, coefficient_d, target_s, target_d)

    call inverse_and_add(ws_x, rhs, 2, coefficient_s, q, sim, nse, &
      js, je, ks, ke)
    call inverse_and_add(ws_y, rhs, 3, coefficient_s, q, sim, nse, &
      js, je, ks, ke)
    call inverse_and_add(ws_z, rhs, 4, coefficient_s, q, sim, nse, &
      js, je, ks, ke)
    call inverse_and_add(wd_x, rhs, 2, coefficient_d, q, sim, nse, &
      js, je, ks, ke)
    call inverse_and_add(wd_y, rhs, 3, coefficient_d, q, sim, nse, &
      js, je, ks, ke)
    call inverse_and_add(wd_z, rhs, 4, coefficient_d, q, sim, nse, &
      js, je, ks, ke)

    evaluation_count = evaluation_count + 1
    if (my_rank == root .and. nse%forcing_report_interval > 0) then
      if (evaluation_count == 1 .or. &
          modulo(evaluation_count, nse%forcing_report_interval) == 0) then
        write(*,'(A,I0,7(1X,ES13.5))') '# forcing', evaluation_count, &
          coefficient_s, coefficient_d, denominator_s, denominator_d, &
          pressure_dilatation, target_s, target_d
      end if
    end if
  end subroutine evaluate_and_add_forcing

  subroutine forward_weighted_velocity(q, variable, field_hat, sim, nse, &
      js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: variable, js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    complex(mytype), intent(out) :: field_hat(zstart(1):, zstart(2):, &
      zstart(3):)

    call redistribute_nse_to_decomp_x(q, variable, work_x, sim, nse, &
      js, je, ks, ke)
    call decomp_2d_fft_3d(work_x, field_hat, DECOMP_2D_FFT_FORWARD)
  end subroutine forward_weighted_velocity

  subroutine project_helmholtz(sim, nse, denominator_s, denominator_d)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    real(dp), intent(out) :: denominator_s, denominator_d
    complex(mytype) :: dot_product, wx, wy, wz
    real(dp) :: pi, kx, ky, kz, k_squared, k_magnitude, point_count
    logical :: retain_mode, low_wavenumber
    integer :: i, j, k, mx, my, mz

    pi = acos(-1.0_dp)
    denominator_s = 0.0_dp
    denominator_d = 0.0_dp
    low_wavenumber = &
      trim(adjustl(nse%forcing_spectrum)) == 'low_wavenumber'

    do k = zstart(3), zend(3)
      mz = signed_mode(k, sim%nz)
      kz = 2.0_dp*pi*real(mz,dp) / sim%lz
      do j = zstart(2), zend(2)
        my = signed_mode(j, sim%ny)
        ky = 2.0_dp*pi*real(my,dp) / sim%ly
        do i = zstart(1), zend(1)
          mx = signed_mode(i, sim%nx)
          kx = 2.0_dp*pi*real(mx,dp) / sim%lx
          k_squared = kx*kx + ky*ky + kz*kz
          k_magnitude = sqrt(k_squared)
          retain_mode = k_squared > 0.0_dp
          if (low_wavenumber) then
            retain_mode = retain_mode .and. &
              k_magnitude < nse%forcing_k_cutoff
          end if

          wx = ws_x(i,j,k)
          wy = ws_y(i,j,k)
          wz = ws_z(i,j,k)
          if (retain_mode) then
            dot_product = (kx*wx + ky*wy + kz*wz) / k_squared
            wd_x(i,j,k) = kx * dot_product
            wd_y(i,j,k) = ky * dot_product
            wd_z(i,j,k) = kz * dot_product
            ws_x(i,j,k) = wx - wd_x(i,j,k)
            ws_y(i,j,k) = wy - wd_y(i,j,k)
            ws_z(i,j,k) = wz - wd_z(i,j,k)
          else
            ws_x(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            ws_y(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            ws_z(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            wd_x(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            wd_y(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            wd_z(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
          end if
          denominator_s = denominator_s + &
            real(ws_x(i,j,k)*conjg(ws_x(i,j,k)) + &
                 ws_y(i,j,k)*conjg(ws_y(i,j,k)) + &
                 ws_z(i,j,k)*conjg(ws_z(i,j,k)), dp)
          denominator_d = denominator_d + &
            real(wd_x(i,j,k)*conjg(wd_x(i,j,k)) + &
                 wd_y(i,j,k)*conjg(wd_y(i,j,k)) + &
                 wd_z(i,j,k)*conjg(wd_z(i,j,k)), dp)
        end do
      end do
    end do
    call mp_allsumr8(denominator_s)
    call mp_allsumr8(denominator_d)
    point_count = real(sim%nx,dp)*real(sim%ny,dp)*real(sim%nz,dp)
    denominator_s = denominator_s / (point_count*point_count)
    denominator_d = denominator_d / (point_count*point_count)
  end subroutine project_helmholtz

  real(dp) function compute_pressure_dilatation(q, sim, nse, js, je, &
      ks, ke) result(value)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: rho, u, v, w, pressure, divergence, point_count
    integer :: i, j, k

    value = 0.0_dp
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          rho = max(q(i,j,k,1), nse%small_rho)
          u = q(i,j,k,2) / rho
          v = q(i,j,k,3) / rho
          w = q(i,j,k,4) / rho
          pressure = max((nse%gamma-1.0_dp) * &
            (q(i,j,k,5)-0.5_dp*rho*(u*u+v*v+w*w)), nse%small_p)
          divergence = derivative_velocity(q, 2, i, j, k, 1, sim%dx, &
            sim, nse, js, ks) + &
            derivative_velocity(q, 3, i, j, k, 2, sim%dy, &
            sim, nse, js, ks) + &
            derivative_velocity(q, 4, i, j, k, 3, sim%dz, &
            sim, nse, js, ks)
          value = value + pressure * divergence
        end do
      end do
    end do
    call mp_allsumr8(value)
    point_count = real(sim%nx,dp)*real(sim%ny,dp)*real(sim%nz,dp)
    value = value / point_count
  end function compute_pressure_dilatation

  real(dp) function derivative_velocity(q, variable, i, j, k, direction, &
      spacing, sim, nse, js, ks) result(derivative)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: variable, i, j, k, direction, js, ks
    real(dp), intent(in) :: spacing
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    integer :: di, dj, dk

    di = merge(1, 0, direction == 1)
    dj = merge(1, 0, direction == 2)
    dk = merge(1, 0, direction == 3)
    derivative = (-velocity(q,variable,i-3*di,j-3*dj,k-3*dk, &
      sim,nse,js,ks) + 9.0_dp*velocity(q,variable,i-2*di,j-2*dj,k-2*dk, &
      sim,nse,js,ks) - 45.0_dp*velocity(q,variable,i-di,j-dj,k-dk, &
      sim,nse,js,ks) + 45.0_dp*velocity(q,variable,i+di,j+dj,k+dk, &
      sim,nse,js,ks) - 9.0_dp*velocity(q,variable,i+2*di,j+2*dj,k+2*dk, &
      sim,nse,js,ks) + velocity(q,variable,i+3*di,j+3*dj,k+3*dk, &
      sim,nse,js,ks)) / (60.0_dp*spacing)
  end function derivative_velocity

  real(dp) function velocity(q, variable, i, j, k, sim, nse, js, ks) &
      result(component)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: variable, i, j, k, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)

    component = q(i,j,k,variable) / max(q(i,j,k,1), nse%small_rho)
  end function velocity

  subroutine inverse_and_add(field_hat, rhs, variable, coefficient, q, &
      sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: variable, js, je, ks, ke
    real(dp), intent(in) :: coefficient
    complex(mytype), intent(inout) :: field_hat(zstart(1):, zstart(2):, &
      zstart(3):)
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(inout) :: rhs(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: normalization
    integer :: i, j, k

    call decomp_2d_fft_3d(field_hat, work_x, DECOMP_2D_FFT_BACKWARD)
    normalization = 1.0_dp / &
      (real(sim%nx,dp)*real(sim%ny,dp)*real(sim%nz,dp))
    call redistribute_decomp_x_to_local(work_x, local_component, sim, &
      js, je, ks, ke, normalization)
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          rhs(i,j,k,variable) = rhs(i,j,k,variable) + coefficient * &
            sqrt(max(q(i,j,k,1), nse%small_rho)) * local_component(i,j,k)
        end do
      end do
    end do
  end subroutine inverse_and_add

  subroutine redistribute_nse_to_decomp_x(q, variable, output, sim, nse, &
      js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: variable, js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    complex(mytype), intent(out) :: output(xstart(1):, xstart(2):, &
      xstart(3):)
    integer, allocatable :: send_counts(:), recv_counts(:)
    integer, allocatable :: send_displs(:), recv_displs(:)
    integer, allocatable :: decomp_bounds(:,:)
    real(dp), allocatable :: send_buffer(:), recv_buffer(:)
    integer :: local_bounds(4), dest, source, source_j, source_k
    integer :: ylo, yhi, zlo, zhi, i, j, k, position, ierr
    integer :: send_total, recv_total

    allocate(send_counts(0:nprocs-1), recv_counts(0:nprocs-1))
    allocate(send_displs(0:nprocs-1), recv_displs(0:nprocs-1))
    allocate(decomp_bounds(4,0:nprocs-1))
    local_bounds = [xstart(2), xend(2), xstart(3), xend(3)]
    call MPI_ALLGATHER(local_bounds, 4, MPI_INTEGER, decomp_bounds, 4, &
      MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) error stop 'MPI_Allgather failed in forcing redistribution'

    do dest = 0, nprocs - 1
      ylo = max(js, decomp_bounds(1,dest))
      yhi = min(je, decomp_bounds(2,dest))
      zlo = max(ks, decomp_bounds(3,dest))
      zhi = min(ke, decomp_bounds(4,dest))
      send_counts(dest) = overlap_size(sim%nx, ylo, yhi, zlo, zhi)
    end do
    do source = 0, nprocs - 1
      source_j = modulo(source, ndiv_ny)
      source_k = source / ndiv_ny
      ylo = max(xstart(2), jjsta(source_j))
      yhi = min(xend(2), jjend(source_j))
      zlo = max(xstart(3), kksta(source_k))
      zhi = min(xend(3), kkend(source_k))
      recv_counts(source) = overlap_size(sim%nx, ylo, yhi, zlo, zhi)
    end do
    call make_displacements(send_counts, send_displs)
    call make_displacements(recv_counts, recv_displs)
    send_total = sum(send_counts)
    recv_total = sum(recv_counts)
    if (send_total /= sim%nx*(je-js+1)*(ke-ks+1) .or. &
        recv_total /= size(output)) then
      error stop 'Invalid NSE-to-2DECOMP forcing overlap'
    end if
    allocate(send_buffer(max(1,send_total)), recv_buffer(max(1,recv_total)))

    do dest = 0, nprocs - 1
      ylo = max(js, decomp_bounds(1,dest))
      yhi = min(je, decomp_bounds(2,dest))
      zlo = max(ks, decomp_bounds(3,dest))
      zhi = min(ke, decomp_bounds(4,dest))
      position = send_displs(dest)
      do k = zlo, zhi
        do j = ylo, yhi
          do i = 1, sim%nx
            position = position + 1
            send_buffer(position) = q(i,j,k,variable) / &
              sqrt(max(q(i,j,k,1), nse%small_rho))
          end do
        end do
      end do
    end do
    call MPI_ALLTOALLV(send_buffer, send_counts, send_displs, &
      MPI_DOUBLE_PRECISION, recv_buffer, recv_counts, recv_displs, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) error stop 'MPI_Alltoallv failed in forcing forward path'

    do source = 0, nprocs - 1
      source_j = modulo(source, ndiv_ny)
      source_k = source / ndiv_ny
      ylo = max(xstart(2), jjsta(source_j))
      yhi = min(xend(2), jjend(source_j))
      zlo = max(xstart(3), kksta(source_k))
      zhi = min(xend(3), kkend(source_k))
      position = recv_displs(source)
      do k = zlo, zhi
        do j = ylo, yhi
          do i = 1, sim%nx
            position = position + 1
            output(i,j,k) = cmplx(recv_buffer(position), 0.0_dp, kind=mytype)
          end do
        end do
      end do
    end do
    deallocate(send_counts, recv_counts, send_displs, recv_displs)
    deallocate(decomp_bounds, send_buffer, recv_buffer)
  end subroutine redistribute_nse_to_decomp_x

  subroutine redistribute_decomp_x_to_local(input, output, sim, js, je, &
      ks, ke, normalization)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(in) :: normalization
    complex(mytype), intent(in) :: input(xstart(1):, xstart(2):, xstart(3):)
    real(dp), intent(out) :: output(1:, js:, ks:)
    integer, allocatable :: send_counts(:), recv_counts(:)
    integer, allocatable :: send_displs(:), recv_displs(:)
    integer, allocatable :: source_bounds(:,:)
    real(dp), allocatable :: send_buffer(:), recv_buffer(:)
    integer :: local_bounds(4), dest, source, dest_j, dest_k
    integer :: ylo, yhi, zlo, zhi, i, j, k, position, ierr
    integer :: send_total, recv_total

    allocate(send_counts(0:nprocs-1), recv_counts(0:nprocs-1))
    allocate(send_displs(0:nprocs-1), recv_displs(0:nprocs-1))
    allocate(source_bounds(4,0:nprocs-1))
    local_bounds = [xstart(2), xend(2), xstart(3), xend(3)]
    call MPI_ALLGATHER(local_bounds, 4, MPI_INTEGER, source_bounds, 4, &
      MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) error stop 'MPI_Allgather failed in forcing inverse path'

    do dest = 0, nprocs - 1
      dest_j = modulo(dest, ndiv_ny)
      dest_k = dest / ndiv_ny
      ylo = max(xstart(2), jjsta(dest_j))
      yhi = min(xend(2), jjend(dest_j))
      zlo = max(xstart(3), kksta(dest_k))
      zhi = min(xend(3), kkend(dest_k))
      send_counts(dest) = overlap_size(sim%nx, ylo, yhi, zlo, zhi)
    end do
    do source = 0, nprocs - 1
      ylo = max(js, source_bounds(1,source))
      yhi = min(je, source_bounds(2,source))
      zlo = max(ks, source_bounds(3,source))
      zhi = min(ke, source_bounds(4,source))
      recv_counts(source) = overlap_size(sim%nx, ylo, yhi, zlo, zhi)
    end do
    call make_displacements(send_counts, send_displs)
    call make_displacements(recv_counts, recv_displs)
    send_total = sum(send_counts)
    recv_total = sum(recv_counts)
    if (send_total /= size(input) .or. &
        recv_total /= sim%nx*(je-js+1)*(ke-ks+1)) then
      error stop 'Invalid 2DECOMP-to-NSE forcing overlap'
    end if
    allocate(send_buffer(max(1,send_total)), recv_buffer(max(1,recv_total)))

    do dest = 0, nprocs - 1
      dest_j = modulo(dest, ndiv_ny)
      dest_k = dest / ndiv_ny
      ylo = max(xstart(2), jjsta(dest_j))
      yhi = min(xend(2), jjend(dest_j))
      zlo = max(xstart(3), kksta(dest_k))
      zhi = min(xend(3), kkend(dest_k))
      position = send_displs(dest)
      do k = zlo, zhi
        do j = ylo, yhi
          do i = 1, sim%nx
            position = position + 1
            send_buffer(position) = real(input(i,j,k),dp) * normalization
          end do
        end do
      end do
    end do
    call MPI_ALLTOALLV(send_buffer, send_counts, send_displs, &
      MPI_DOUBLE_PRECISION, recv_buffer, recv_counts, recv_displs, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) error stop 'MPI_Alltoallv failed in forcing inverse path'

    do source = 0, nprocs - 1
      ylo = max(js, source_bounds(1,source))
      yhi = min(je, source_bounds(2,source))
      zlo = max(ks, source_bounds(3,source))
      zhi = min(ke, source_bounds(4,source))
      position = recv_displs(source)
      do k = zlo, zhi
        do j = ylo, yhi
          do i = 1, sim%nx
            position = position + 1
            output(i,j,k) = recv_buffer(position)
          end do
        end do
      end do
    end do
    deallocate(send_counts, recv_counts, send_displs, recv_displs)
    deallocate(source_bounds, send_buffer, recv_buffer)
  end subroutine redistribute_decomp_x_to_local

  subroutine make_displacements(counts, displacements)
    integer, intent(in) :: counts(0:)
    integer, intent(out) :: displacements(0:)
    integer :: rank

    displacements(0) = 0
    do rank = 1, ubound(counts,1)
      displacements(rank) = displacements(rank-1) + counts(rank-1)
    end do
  end subroutine make_displacements

  pure integer function overlap_size(nx, ylo, yhi, zlo, zhi) result(count)
    integer, intent(in) :: nx, ylo, yhi, zlo, zhi

    if (yhi < ylo .or. zhi < zlo) then
      count = 0
    else
      count = nx * (yhi-ylo+1) * (zhi-zlo+1)
    end if
  end function overlap_size

  pure integer function signed_mode(global_index, n) result(mode)
    integer, intent(in) :: global_index, n

    mode = global_index - 1
    if (mode > n/2) mode = mode - n
  end function signed_mode

end module mod_nse_forcing
