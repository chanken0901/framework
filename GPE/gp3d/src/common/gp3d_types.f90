!> GP3D全体で共有する数値精度、格子、実行条件、GPE条件、波動関数の型定義。
!> このモジュールは計算処理を持たず、各モジュール間で受け渡すデータ構造だけを所有する。
module gp3d_types
  implicit none
  private

  integer, parameter, public :: dp = selected_real_kind(15, 307)

  real(dp), parameter, public :: pi = 4.0_dp * atan(1.0_dp)

  !> 全体格子とMPI rankが担当するzスラブ、および実空間・波数空間座標。
  type, public :: gp3d_grid_t
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: local_nz = 0
    integer :: k_start = 1
    integer :: k_end = 0
    integer :: rank = 0
    integer :: nprocs = 1
    real(dp) :: lx = 0.0_dp
    real(dp) :: ly = 0.0_dp
    real(dp) :: lz = 0.0_dp
    real(dp) :: dx = 0.0_dp
    real(dp) :: dy = 0.0_dp
    real(dp) :: dz = 0.0_dp
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: kx(:)
    real(dp), allocatable :: ky(:)
    real(dp), allocatable :: kz(:)
  end type gp3d_grid_t

  !> 時間積分ループで頻繁に参照する、導出済みの最小パラメータ集合。
  type, public :: gp3d_params_t
    real(dp) :: dt = 1.0e-3_dp
    real(dp) :: mass = 1.0_dp
    real(dp) :: hbar = 1.0_dp
    real(dp) :: g = 1.0_dp
    real(dp) :: norm = 1.0_dp
    integer :: nsteps = 100
    integer :: output_every = 10
    logical :: imaginary_time = .false.
  end type gp3d_params_t

  !> &simulation namelistから読み込む格子、時間、出力、実行バックエンド設定。
  type, public :: gp3d_run_config_t
    character(len=32) :: equation = "GPE"
    character(len=256) :: case_name = "gp3d"
    character(len=256) :: input_file = "input.nml"
    character(len=64) :: initial_condition = "gaussian"
    character(len=256) :: restart_file = ""
    character(len=256) :: output_dir = "output"
    character(len=32) :: output_format = "slf"
    character(len=32) :: precision_name = "float64"
    character(len=32) :: backend = "dft"
    integer :: nx = 16
    integer :: ny = 16
    integer :: nz = 16
    integer :: nghost = 0
    real(dp) :: x_min = -6.0_dp
    real(dp) :: x_max = 6.0_dp
    real(dp) :: y_min = -6.0_dp
    real(dp) :: y_max = 6.0_dp
    real(dp) :: z_min = -6.0_dp
    real(dp) :: z_max = 6.0_dp
    real(dp) :: lx = 12.0_dp
    real(dp) :: ly = 12.0_dp
    real(dp) :: lz = 12.0_dp
    real(dp) :: dx = 0.75_dp
    real(dp) :: dy = 0.75_dp
    real(dp) :: dz = 0.75_dp
    real(dp) :: dt = 2.5e-4_dp
    real(dp) :: t_max = 0.0_dp
    real(dp) :: cfl = 0.0_dp
    integer :: nsteps = 20
    integer :: output_frequency = 10
    logical :: use_fixed_dt = .true.
    logical :: write_initial = .true.
    logical :: write_meta = .true.
    logical :: timing_enabled = .false.
    logical :: use_mpi = .false.
    logical :: use_openmp = .false.
    logical :: use_cuda = .false.
    integer :: cuda_device = 0
    integer :: rank = 0
    integer :: nprocs = 1
  end type gp3d_run_config_t

  !> &gpe namelistから読み込む物理量と初期条件・ARGLE固有パラメータ。
  type, public :: gp3d_model_config_t
    logical :: use_dimensionless_parameters = .false.
    real(dp) :: alpha = 0.05_dp
    real(dp) :: beta = 40.0_dp
    real(dp) :: g = 1.0_dp
    real(dp) :: sigma0 = 1.0_dp
    real(dp) :: wx = 1.0_dp
    real(dp) :: wy = 1.0_dp
    real(dp) :: wz = 1.0_dp
    real(dp) :: hbar = 1.0_dp
    real(dp) :: mass = 1.0_dp
    real(dp) :: norm = 1.0_dp
    real(dp) :: mu = 8.0_dp
    real(dp) :: healing_length = 0.25_dp
    integer :: vortex_charge = 1
    real(dp) :: vortex_x0 = 0.0_dp
    real(dp) :: vortex_y0 = 0.0_dp
    real(dp) :: ring_radius = 2.0_dp
    real(dp) :: ring_z0 = 0.0_dp
    real(dp) :: phase_noise = 0.0_dp
    integer :: tangle_nlines = 8
    integer :: tangle_nrings = 8
    real(dp) :: ring_radius_min = 0.8_dp
    real(dp) :: ring_radius_max = 1.8_dp
    integer :: random_seed = 12345
    real(dp) :: tg_velocity_amplitude = 1.0_dp
    integer :: tg_winding = 1
    logical :: tg_auto_winding = .true.
    logical :: argle_enabled = .false.
    logical :: argle_write_seed = .true.
    integer :: argle_steps = 0
    integer :: argle_output_every = 100
    real(dp) :: argle_dtau = 1.0e-3_dp
    real(dp) :: argle_tolerance = 1.0e-8_dp
    logical :: imaginary_time = .true.
  end type gp3d_model_config_t

  !> 各rankが保持する複素波動関数psiと外部ポテンシャル。
  type, public :: gp3d_state_t
    complex(dp), allocatable :: psi(:,:,:)
    real(dp), allocatable :: potential(:,:,:)
  end type gp3d_state_t

end module gp3d_types
