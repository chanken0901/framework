!> MPI分散FFTが各軸の1次元変換に使う、外部ライブラリ不要の参照DFT。
!> 分散転置の検証を主目的とし、大規模計算にはlocal FFTW版を使用する。
module gp3d_local_fft
  use gp3d_types, only: dp, pi
  implicit none
  private

  public :: gp3d_local_fft_plan_t
  public :: gp3d_local_fft_init
  public :: gp3d_local_fft_execute
  public :: gp3d_local_fft_finalize

  type :: gp3d_local_fft_plan_t
    integer :: n = 0
  end type gp3d_local_fft_plan_t

contains

  subroutine gp3d_local_fft_init(plan, n)
    type(gp3d_local_fft_plan_t), intent(out) :: plan
    integer, intent(in) :: n

    if (n <= 0) error stop "local FFT length must be positive"
    plan%n = n
  end subroutine gp3d_local_fft_init

  subroutine gp3d_local_fft_execute(plan, input, output, sign)
    type(gp3d_local_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%n)
    complex(dp), intent(out) :: output(plan%n)
    integer, intent(in) :: sign

    integer :: frequency, point
    real(dp) :: phase
    complex(dp) :: sum_value

    if (sign /= -1 .and. sign /= 1) error stop "local FFT sign must be -1 or 1"
    do frequency = 1, plan%n
      sum_value = (0.0_dp, 0.0_dp)
      do point = 1, plan%n
        phase = real(sign, dp) * 2.0_dp * pi * &
          real((frequency - 1) * (point - 1), dp) / real(plan%n, dp)
        sum_value = sum_value + input(point) * cmplx(cos(phase), sin(phase), kind=dp)
      end do
      output(frequency) = sum_value
    end do
  end subroutine gp3d_local_fft_execute

  subroutine gp3d_local_fft_finalize(plan)
    type(gp3d_local_fft_plan_t), intent(inout) :: plan

    plan%n = 0
  end subroutine gp3d_local_fft_finalize

end module gp3d_local_fft
