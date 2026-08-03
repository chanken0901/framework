module mod_init_hit_spectral
  use, intrinsic :: iso_fortran_env, only : int64
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use module_mpi, only : my_rank, root, nprocs, ndiv_ny, ndiv_nz, &
    jjsta, jjend, kksta, kkend, mp_allmaxr8, mp_allsumr8, &
    MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION
  use decomp_2d, only : decomp_info, decomp_2d_init, decomp_2d_finalize, &
    alloc_x, alloc_z, xstart, xend, zstart, zend
  use decomp_2d_fft, only : decomp_2d_fft_init, decomp_2d_fft_finalize, &
    decomp_2d_fft_3d, decomp_2d_fft_get_ph
  use decomp_2d_constants, only : mytype, PHYSICAL_IN_X, &
    DECOMP_2D_FFT_BACKWARD
  implicit none
  private

  public :: initialize_hit_spectral

contains

  subroutine initialize_hit_spectral(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    type(decomp_info), pointer :: ph
    complex(mytype), allocatable :: work_x(:,:,:)
    complex(mytype), allocatable :: u_hat(:,:,:), v_hat(:,:,:), w_hat(:,:,:)
    real(dp) :: target_rms, turbulent_mach, measured_rms, scale
    real(dp) :: velocity_square_sum
    real(dp) :: spectral_divergence, imaginary_residual(3)
    real(dp) :: rho, pressure, u, v, w
    integer :: i, j, k

    if (nse%nv /= 5) then
      error stop 'Spectral HIT initialization requires five conserved variables'
    end if
    if (storage_size(0.0_mytype) /= storage_size(0.0_dp)) then
      error stop '2DECOMP&FFT precision must match mod_precision dp'
    end if
    if (sim%lx <= 0.0_dp .or. sim%ly <= 0.0_dp .or. sim%lz <= 0.0_dp) then
      error stop 'Spectral HIT initialization requires positive domain lengths'
    end if
    if (nse%hit_dealias_fraction <= 0.0_dp .or. &
        nse%hit_dealias_fraction > 1.0_dp) then
      error stop 'hit_dealias_fraction must be in the interval (0, 1]'
    end if

    target_rms = nse%hit_rms_velocity
    if (target_rms <= 0.0_dp) target_rms = nse%mach / sqrt(3.0_dp)
    if (target_rms <= 0.0_dp) then
      error stop 'hit_rms_velocity or mach must be positive for HIT initialization'
    end if
    turbulent_mach = sqrt(3.0_dp) * target_rms

    ! The NSE MPI layout is an x-pencil: x is local-complete while y and z
    ! are divided by ndiv_ny and ndiv_nz. Use the same process grid here.
    call decomp_2d_init(sim%nx, sim%ny, sim%nz, ndiv_ny, ndiv_nz, &
      complex_pool=.true.)
    call decomp_2d_fft_init(PHYSICAL_IN_X)
    ph => decomp_2d_fft_get_ph()

    call validate_decomp_x_extent(sim)
    call alloc_x(work_x, ph, .true.)
    call alloc_z(u_hat, ph, .true.)
    call alloc_z(v_hat, ph, .true.)
    call alloc_z(w_hat, ph, .true.)

    call generate_hit_fourier_field(u_hat, v_hat, w_hat, sim, nse, &
      spectral_divergence)

    call inverse_component_to_state(u_hat, work_x, q, sim, js, je, ks, &
      ke, 2, imaginary_residual(1))
    call inverse_component_to_state(v_hat, work_x, q, sim, js, je, ks, &
      ke, 3, imaginary_residual(2))
    call inverse_component_to_state(w_hat, work_x, q, sim, js, je, ks, &
      ke, 4, imaginary_residual(3))

    velocity_square_sum = 0.0_dp
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          velocity_square_sum = velocity_square_sum + q(i,j,k,2)**2 + &
            q(i,j,k,3)**2 + q(i,j,k,4)**2
        end do
      end do
    end do
    call mp_allsumr8(velocity_square_sum)
    measured_rms = sqrt(velocity_square_sum / &
      (3.0_dp * real(sim%nx,dp) * real(sim%ny,dp) * real(sim%nz,dp)))
    if (measured_rms <= tiny(1.0_dp)) then
      error stop 'Spectral HIT initialization produced a zero velocity field'
    end if
    scale = target_rms / measured_rms

    rho = nse%rho0
    pressure = nse%rho0 / nse%gamma
    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          u = scale * q(i,j,k,2)
          v = scale * q(i,j,k,3)
          w = scale * q(i,j,k,4)
          q(i,j,k,1) = rho
          q(i,j,k,2) = rho * u
          q(i,j,k,3) = rho * v
          q(i,j,k,4) = rho * w
          q(i,j,k,5) = pressure / (nse%gamma - 1.0_dp) + &
            0.5_dp * rho * (u*u + v*v + w*w)
        end do
      end do
    end do

    if (my_rank == root) then
      write(*,'(A)') '# Distributed spectral HIT initialization'
      write(*,'(A,A)') '# spectrum: ', trim(nse%hit_spectrum)
      write(*,'(A,I0,A,I0)') '# 2DECOMP process grid: ', ndiv_ny, ' x ', ndiv_nz
      write(*,'(A)') '# Fourier coefficients: direct transverse-mode construction'
      write(*,'(A,ES16.8)') '# max spectral divergence: ', spectral_divergence
      write(*,'(A,ES16.8)') '# max inverse FFT imaginary residual: ', &
        maxval(imaginary_residual)
      write(*,'(A,ES16.8)') '# unscaled component RMS: ', measured_rms
      write(*,'(A,ES16.8)') '# target component RMS:   ', target_rms
      write(*,'(A,ES16.8)') '# initial turbulent Mach number: ', &
        turbulent_mach
      if (turbulent_mach > 1.0_dp) then
        write(*,'(A)') '# WARNING: the requested HIT field is initially supersonic.'
        write(*,'(A)') '# hit_rms_velocity is nondimensional, not a velocity in m/s.'
      end if
    end if

    deallocate(work_x, u_hat, v_hat, w_hat)
    nullify(ph)
    call decomp_2d_fft_finalize
    call decomp_2d_finalize
  end subroutine initialize_hit_spectral

  subroutine validate_decomp_x_extent(sim)
    type(simulation_config), intent(in) :: sim

    if (xstart(1) /= 1 .or. xend(1) /= sim%nx) then
      error stop '2DECOMP physical layout is not an x-pencil'
    end if
  end subroutine validate_decomp_x_extent

  subroutine generate_hit_fourier_field(u_hat, v_hat, w_hat, sim, nse, &
      max_divergence)
    complex(mytype), intent(out) :: u_hat(zstart(1):, zstart(2):, zstart(3):)
    complex(mytype), intent(out) :: v_hat(zstart(1):, zstart(2):, zstart(3):)
    complex(mytype), intent(out) :: w_hat(zstart(1):, zstart(2):, zstart(3):)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    real(dp), intent(out) :: max_divergence
    complex(mytype) :: a_mode, b_mode, phase_a, phase_b
    complex(mytype) :: u_canonical, v_canonical, w_canonical
    real(dp) :: pi, kx, ky, kz, kxy, kmag, amplitude
    real(dp) :: e1x, e1y, e1z, e2x, e2y, e2z
    real(dp) :: theta_a, theta_b, theta_mix, local_divergence
    integer :: mx, my, mz, cmx, cmy, cmz, orientation, i, j, k
    integer :: cutoff_x, cutoff_y, cutoff_z

    pi = acos(-1.0_dp)
    cutoff_x = floor(nse%hit_dealias_fraction * real(sim%nx/2,dp))
    cutoff_y = floor(nse%hit_dealias_fraction * real(sim%ny/2,dp))
    cutoff_z = floor(nse%hit_dealias_fraction * real(sim%nz/2,dp))
    max_divergence = 0.0_dp

    do k = zstart(3), zend(3)
      mz = signed_mode(k, sim%nz)
      do j = zstart(2), zend(2)
        my = signed_mode(j, sim%ny)
        do i = zstart(1), zend(1)
          mx = signed_mode(i, sim%nx)

          if ((mx == 0 .and. my == 0 .and. mz == 0) .or. &
              abs(mx) > cutoff_x .or. &
              abs(my) > cutoff_y .or. abs(mz) > cutoff_z .or. &
              is_nyquist_mode(mx,sim%nx) .or. &
              is_nyquist_mode(my,sim%ny) .or. &
              is_nyquist_mode(mz,sim%nz)) then
            u_hat(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            v_hat(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            w_hat(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=mytype)
            cycle
          end if

          ! Generate one representative of each {k,-k} pair. The opposite
          ! coefficient is its complex conjugate, so the inverse FFT is real.
          orientation = canonical_orientation(mx, my, mz)
          cmx = orientation * mx
          cmy = orientation * my
          cmz = orientation * mz
          kx = 2.0_dp * pi * real(cmx,dp) / sim%lx
          ky = 2.0_dp * pi * real(cmy,dp) / sim%ly
          kz = 2.0_dp * pi * real(cmz,dp) / sim%lz
          kxy = sqrt(kx*kx + ky*ky)
          kmag = sqrt(kx*kx + ky*ky + kz*kz)

          ! These are the two transverse basis vectors used by the original
          ! HIT code. Both are perpendicular to k and mutually orthogonal.
          if (kxy > tiny(1.0_dp)) then
            e1x = ky / kxy
            e1y = -kx / kxy
            e1z = 0.0_dp
            e2x = kx * kz / (kxy * kmag)
            e2y = ky * kz / (kxy * kmag)
            e2z = -kxy / kmag
          else
            e1x = 1.0_dp
            e1y = 0.0_dp
            e1z = 0.0_dp
            e2x = 0.0_dp
            e2y = 1.0_dp
            e2z = 0.0_dp
          end if

          amplitude = spectrum_mode_amplitude(kmag, nse)
          theta_a = 2.0_dp*pi*deterministic_uniform(cmx, cmy, cmz, &
            nse%hit_seed, 1) + 0.5_dp*pi
          theta_b = 2.0_dp*pi*deterministic_uniform(cmx, cmy, cmz, &
            nse%hit_seed, 2)
          theta_mix = 2.0_dp*pi*deterministic_uniform(cmx, cmy, cmz, &
            nse%hit_seed, 3)
          phase_a = cmplx(cos(theta_a), sin(theta_a), kind=mytype)
          phase_b = cmplx(cos(theta_b), sin(theta_b), kind=mytype)
          a_mode = amplitude * cos(theta_mix) * phase_a
          b_mode = amplitude * sin(theta_mix) * phase_b

          u_canonical = e1x*a_mode + e2x*b_mode
          v_canonical = e1y*a_mode + e2y*b_mode
          w_canonical = e1z*a_mode + e2z*b_mode
          if (orientation > 0) then
            u_hat(i,j,k) = u_canonical
            v_hat(i,j,k) = v_canonical
            w_hat(i,j,k) = w_canonical
          else
            u_hat(i,j,k) = conjg(u_canonical)
            v_hat(i,j,k) = conjg(v_canonical)
            w_hat(i,j,k) = conjg(w_canonical)
          end if

          kx = 2.0_dp * pi * real(mx,dp) / sim%lx
          ky = 2.0_dp * pi * real(my,dp) / sim%ly
          kz = 2.0_dp * pi * real(mz,dp) / sim%lz
          local_divergence = abs(kx*u_hat(i,j,k) + ky*v_hat(i,j,k) + &
            kz*w_hat(i,j,k))
          max_divergence = max(max_divergence, local_divergence)
        end do
      end do
    end do
    call mp_allmaxr8(max_divergence)
  end subroutine generate_hit_fourier_field

  real(dp) function spectrum_mode_amplitude(kmag, nse) result(amplitude)
    real(dp), intent(in) :: kmag
    type(nse_config), intent(in) :: nse
    real(dp), parameter :: pope_c_l = 6.78_dp
    real(dp), parameter :: pope_c_eta = 0.40_dp
    real(dp), parameter :: pope_p0 = 2.0_dp
    real(dp), parameter :: pope_beta = 5.2_dp
    real(dp) :: energy_shape, kl, keta, f_l, f_eta, ratio

    select case (trim(adjustl(nse%hit_spectrum)))
    case ('johnsen', 'k4_gaussian')
      if (nse%hit_peak_wavenumber <= 0.0_dp) then
        error stop 'hit_peak_wavenumber must be positive for Johnsen spectrum'
      end if
      ratio = kmag / nse%hit_peak_wavenumber
      energy_shape = ratio**4 * exp(-2.0_dp * ratio**2)
    case ('pope')
      if (nse%hit_integral_length <= 0.0_dp .or. &
          nse%hit_kolmogorov_length <= 0.0_dp) then
        error stop 'Pope spectrum requires positive HIT length scales'
      end if
      kl = kmag * nse%hit_integral_length
      keta = kmag * nse%hit_kolmogorov_length
      f_l = (kl / sqrt(kl*kl + pope_c_l))**(5.0_dp/3.0_dp + pope_p0)
      f_eta = exp(-pope_beta * ((keta**4 + pope_c_eta**4)**0.25_dp - &
        pope_c_eta))
      energy_shape = kmag**(-5.0_dp/3.0_dp) * f_l * f_eta
    case default
      write(*,'(A,A)') 'ERROR: unsupported HIT spectrum: ', &
        trim(nse%hit_spectrum)
      error stop
    end select

    ! This is the original HIT coefficient magnitude
    ! sqrt(2 E(k) / (4 pi k^2)). A final global rescaling sets the requested
    ! component RMS without changing the spectral shape.
    amplitude = sqrt(2.0_dp * max(energy_shape, 0.0_dp) / &
      (4.0_dp * acos(-1.0_dp) * kmag * kmag))
  end function spectrum_mode_amplitude

  pure integer function canonical_orientation(mx, my, mz) result(orientation)
    integer, intent(in) :: mx, my, mz

    if (mx /= 0) then
      orientation = merge(1, -1, mx > 0)
    else if (my /= 0) then
      orientation = merge(1, -1, my > 0)
    else if (mz /= 0) then
      orientation = merge(1, -1, mz > 0)
    else
      orientation = 1
    end if
  end function canonical_orientation

  pure real(dp) function deterministic_uniform(mx, my, mz, seed, stream) &
      result(value)
    integer, intent(in) :: mx, my, mz, seed, stream
    integer(int64), parameter :: prime = 2147483647_int64
    integer(int64) :: hash

    hash = modulo(104729_int64 * int(mx,int64) + &
      130363_int64 * int(my,int64) + 15485863_int64 * int(mz,int64) + &
      32452843_int64 * int(stream,int64) + &
      49979687_int64 * int(seed,int64), prime)
    ! Squaring modulo the 31-bit prime breaks the linear phase correlation
    ! produced by a plain linear congruential hash. The product remains
    ! within the signed int64 range because (prime-1)^2 < huge(int64).
    hash = modulo(hash*hash + 104729_int64*int(stream,int64) + &
      12345_int64, prime)
    hash = modulo(hash*hash + 130363_int64*int(stream,int64) + &
      67891_int64, prime)
    hash = modulo(hash * 48271_int64 + 69621_int64, prime)
    value = (real(hash,dp) + 0.5_dp) / real(prime,dp)
  end function deterministic_uniform

  pure integer function signed_mode(global_index, n) result(mode)
    integer, intent(in) :: global_index, n

    mode = global_index - 1
    if (mode > n/2) mode = mode - n
  end function signed_mode

  pure logical function is_nyquist_mode(mode, n) result(is_nyquist)
    integer, intent(in) :: mode, n

    is_nyquist = (modulo(n,2) == 0 .and. abs(mode) == n/2)
  end function is_nyquist_mode

  subroutine inverse_component_to_state(field_hat, work_x, q, sim, &
      js, je, ks, ke, variable, imaginary_residual)
    complex(mytype), intent(inout) :: field_hat(zstart(1):, zstart(2):, &
      zstart(3):)
    complex(mytype), intent(inout) :: work_x(xstart(1):, xstart(2):, &
      xstart(3):)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: js, je, ks, ke, variable
    real(dp), intent(out) :: imaginary_residual
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp) :: fft_normalization

    call decomp_2d_fft_3d(field_hat, work_x, DECOMP_2D_FFT_BACKWARD)
    fft_normalization = 1.0_dp / &
      (real(sim%nx,dp) * real(sim%ny,dp) * real(sim%nz,dp))
    imaginary_residual = maxval(abs(aimag(work_x))) * fft_normalization
    call mp_allmaxr8(imaginary_residual)
    call redistribute_decomp_x_to_nse(work_x, q, sim, js, je, ks, ke, &
      variable, fft_normalization)
  end subroutine inverse_component_to_state

  subroutine redistribute_decomp_x_to_nse(work_x, q, sim, js, je, ks, &
      ke, variable, normalization)
    complex(mytype), intent(in) :: work_x(xstart(1):, xstart(2):, &
      xstart(3):)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: js, je, ks, ke, variable
    real(dp), intent(in) :: normalization
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    integer, allocatable :: send_counts(:), recv_counts(:)
    integer, allocatable :: send_displs(:), recv_displs(:)
    integer, allocatable :: source_bounds(:,:)
    real(dp), allocatable :: send_buffer(:), recv_buffer(:)
    integer :: local_bounds(4)
    integer :: dest, source, dest_j, dest_k
    integer :: ylo, yhi, zlo, zhi
    integer :: i, j, k, position, ierr
    integer :: send_total, recv_total, expected

    allocate(send_counts(0:nprocs-1), recv_counts(0:nprocs-1))
    allocate(send_displs(0:nprocs-1), recv_displs(0:nprocs-1))
    allocate(source_bounds(4,0:nprocs-1))

    local_bounds = [xstart(2), xend(2), xstart(3), xend(3)]
    call MPI_ALLGATHER(local_bounds, 4, MPI_INTEGER, source_bounds, 4, &
      MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) error stop 'MPI_Allgather failed during HIT redistribution'

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

    send_displs(0) = 0
    recv_displs(0) = 0
    do source = 1, nprocs - 1
      send_displs(source) = send_displs(source-1) + send_counts(source-1)
      recv_displs(source) = recv_displs(source-1) + recv_counts(source-1)
    end do
    send_total = sum(send_counts)
    recv_total = sum(recv_counts)
    expected = sim%nx * (je-js+1) * (ke-ks+1)
    if (send_total /= size(work_x) .or. recv_total /= expected) then
      error stop 'Invalid x-pencil overlap during HIT redistribution'
    end if

    allocate(send_buffer(max(1,send_total)))
    allocate(recv_buffer(max(1,recv_total)))
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
            send_buffer(position) = real(work_x(i,j,k),dp) * normalization
          end do
        end do
      end do
    end do

    call MPI_ALLTOALLV(send_buffer, send_counts, send_displs, &
      MPI_DOUBLE_PRECISION, recv_buffer, recv_counts, recv_displs, &
      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) error stop 'MPI_Alltoallv failed during HIT redistribution'

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
            q(i,j,k,variable) = recv_buffer(position)
          end do
        end do
      end do
    end do

    deallocate(send_counts, recv_counts, send_displs, recv_displs)
    deallocate(source_bounds, send_buffer, recv_buffer)
  end subroutine redistribute_decomp_x_to_nse

  pure integer function overlap_size(nx, ylo, yhi, zlo, zhi) result(count)
    integer, intent(in) :: nx, ylo, yhi, zlo, zhi

    if (yhi < ylo .or. zhi < zlo) then
      count = 0
    else
      count = nx * (yhi-ylo+1) * (zhi-zlo+1)
    end if
  end function overlap_size

end module mod_init_hit_spectral
