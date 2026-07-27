!> 外部FFTライブラリを使わない、逐次3次元DFTの参照実装。
!> 小規模テストと正解比較用であり、大規模計算にはFFTWまたはcuFFTを使用する。
module gp3d_fft
  use gp3d_types, only: dp, pi
  implicit none
  private

  public :: gp3d_fft_plan_t
  public :: gp3d_fft_init
  public :: gp3d_fft_forward
  public :: gp3d_fft_inverse
  public :: gp3d_fft_finalize

  type :: gp3d_fft_plan_t
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
  end type gp3d_fft_plan_t

contains

  subroutine gp3d_fft_init(plan, nx, ny, nz, comm, rank, nprocs)
    type(gp3d_fft_plan_t), intent(out) :: plan
    integer, intent(in) :: nx, ny, nz
    integer, intent(in), optional :: comm, rank, nprocs

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) error stop "FFT dimensions must be positive"
    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    if (present(comm)) continue
    if (present(rank)) continue
    if (present(nprocs)) continue
  end subroutine gp3d_fft_init

  subroutine gp3d_fft_forward(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%nz)

    call dft3(plan, input, output, -1.0_dp, .false.)
  end subroutine gp3d_fft_forward

  subroutine gp3d_fft_inverse(plan, input, output)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%nz)

    call dft3(plan, input, output, 1.0_dp, .true.)
  end subroutine gp3d_fft_inverse

  subroutine gp3d_fft_finalize(plan)
    type(gp3d_fft_plan_t), intent(inout) :: plan

    plan%nx = 0
    plan%ny = 0
    plan%nz = 0
  end subroutine gp3d_fft_finalize

  subroutine dft3(plan, input, output, sign, normalize)
    type(gp3d_fft_plan_t), intent(in) :: plan
    complex(dp), intent(in) :: input(plan%nx, plan%ny, plan%nz)
    complex(dp), intent(out) :: output(plan%nx, plan%ny, plan%nz)
    real(dp), intent(in) :: sign
    logical, intent(in) :: normalize

    integer :: i, j, k, p, q, r
    real(dp) :: phase
    complex(dp) :: factor, sum_value
    real(dp) :: scale

    scale = 1.0_dp
    if (normalize) scale = 1.0_dp / real(plan%nx * plan%ny * plan%nz, dp)

    do k = 1, plan%nz
      do j = 1, plan%ny
        do i = 1, plan%nx
          sum_value = (0.0_dp, 0.0_dp)
          do r = 1, plan%nz
            do q = 1, plan%ny
              do p = 1, plan%nx
                phase = sign * 2.0_dp * pi * ( &
                  real((i - 1) * (p - 1), dp) / real(plan%nx, dp) + &
                  real((j - 1) * (q - 1), dp) / real(plan%ny, dp) + &
                  real((k - 1) * (r - 1), dp) / real(plan%nz, dp))
                factor = cmplx(cos(phase), sin(phase), kind=dp)
                sum_value = sum_value + input(p, q, r) * factor
              end do
            end do
          end do
          output(i, j, k) = scale * sum_value
        end do
      end do
    end do
  end subroutine dft3

end module gp3d_fft
