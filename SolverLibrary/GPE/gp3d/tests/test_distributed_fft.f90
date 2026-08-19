program test_distributed_fft
  use gp3d_types, only: dp, pi
  use gp3d_mpi, only: gp3d_mpi_t, gp3d_mpi_init, gp3d_mpi_finalize, &
    gp3d_mpi_is_root, gp3d_mpi_z_range, gp3d_mpi_sum_real
  use gp3d_fft, only: gp3d_fft_plan_t, gp3d_fft_init, gp3d_fft_forward, &
    gp3d_fft_inverse, gp3d_fft_finalize
  implicit none

  integer, parameter :: nx = 4, ny = 4, nz = 4
  integer, parameter :: mode_x = 1, mode_y = 1, mode_z = 1
  type(gp3d_mpi_t) :: mpi
  type(gp3d_fft_plan_t) :: plan
  complex(dp), allocatable :: input(:,:,:), spectrum(:,:,:), restored(:,:,:)
  complex(dp) :: expected
  real(dp) :: phase, local_spectral_error, spectral_error
  real(dp) :: local_roundtrip_error, roundtrip_error
  integer :: i, j, k, kg, k_start, k_end, local_nz

  call gp3d_mpi_init(mpi)
  if (mpi%nprocs > min(ny, nz)) error stop "test process count exceeds slab dimensions"
  call gp3d_mpi_z_range(nz, mpi%rank, mpi%nprocs, k_start, k_end)
  local_nz = k_end - k_start + 1
  allocate(input(nx, ny, local_nz), spectrum(nx, ny, local_nz), restored(nx, ny, local_nz))

  do k = 1, local_nz
    kg = k_start + k - 1
    do j = 1, ny
      do i = 1, nx
        phase = 2.0_dp * pi * (real(mode_x * (i - 1), dp) / real(nx, dp) + &
          real(mode_y * (j - 1), dp) / real(ny, dp) + &
          real(mode_z * (kg - 1), dp) / real(nz, dp))
        input(i,j,k) = cmplx(cos(phase), sin(phase), kind=dp)
      end do
    end do
  end do

  call gp3d_fft_init(plan, nx, ny, nz, mpi%comm, mpi%rank, mpi%nprocs)
  call gp3d_fft_forward(plan, input, spectrum)

  local_spectral_error = 0.0_dp
  do k = 1, local_nz
    kg = k_start + k - 1
    do j = 1, ny
      do i = 1, nx
        expected = (0.0_dp, 0.0_dp)
        if (i == mode_x + 1 .and. j == mode_y + 1 .and. kg == mode_z + 1) then
          expected = cmplx(real(nx * ny * nz, dp), 0.0_dp, kind=dp)
        end if
        local_spectral_error = local_spectral_error + abs(spectrum(i,j,k) - expected)**2
      end do
    end do
  end do
  call gp3d_mpi_sum_real(mpi, local_spectral_error, spectral_error)

  call gp3d_fft_inverse(plan, spectrum, restored)
  local_roundtrip_error = sum(abs(restored - input)**2)
  call gp3d_mpi_sum_real(mpi, local_roundtrip_error, roundtrip_error)

  if (gp3d_mpi_is_root(mpi)) then
    write(*,'(A,ES12.4)') "spectral_error=", sqrt(spectral_error)
    write(*,'(A,ES12.4)') "roundtrip_error=", sqrt(roundtrip_error)
  end if
  if (sqrt(spectral_error) > 1.0e-10_dp) error stop "distributed FFT spectral layout is incorrect"
  if (sqrt(roundtrip_error) > 1.0e-10_dp) error stop "distributed FFT round trip failed"

  call gp3d_fft_finalize(plan)
  call gp3d_mpi_finalize(mpi)
end program test_distributed_fft
