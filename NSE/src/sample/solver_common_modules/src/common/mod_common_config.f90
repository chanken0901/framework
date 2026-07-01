module mod_common_config
  use mod_precision, only : dp
  implicit none
  private

  public :: simulation_config
  public :: init_simulation_config
  public :: update_derived_config
  public :: print_simulation_config
  public :: should_output

  type :: simulation_config
    ! --- model/case ---
    character(len=32)  :: equation = 'NES'      ! GPE, NSE, LES, Euler, ...
    character(len=256) :: case_name = 'case0001'
    character(len=256) :: input_file = 'input.dat'
    character(len=64)  :: initial_condition = 'default'

    ! --- grid ---
    integer :: nx = 64
    integer :: ny = 64
    integer :: nz = 64
    integer :: nghost = 3

    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 1.0_dp
    real(dp) :: z_min = 0.0_dp
    real(dp) :: z_max = 1.0_dp

    real(dp) :: lx = 1.0_dp
    real(dp) :: ly = 1.0_dp
    real(dp) :: lz = 1.0_dp
    real(dp) :: dx = 1.0_dp
    real(dp) :: dy = 1.0_dp
    real(dp) :: dz = 1.0_dp

    ! --- time ---
    real(dp) :: dt = 1.0e-4_dp
    real(dp) :: t_max = 1.0_dp
    integer  :: nsteps = 1000
    real(dp) :: cfl = 0.5_dp
    logical  :: use_fixed_dt = .true.

    ! --- output ---
    integer :: output_frequency = 100
    character(len=256) :: output_dir = 'output'
    character(len=32)  :: output_format = 'slf'  ! slf/bin/vtk etc.
    character(len=32)  :: precision_name = 'float64'
    logical :: write_initial = .true.
    logical :: write_meta = .true.

    ! --- execution ---
    character(len=32) :: backend = 'serial'
    logical :: use_mpi = .false.
    logical :: use_openmp = .false.
    integer :: rank = 0
    integer :: nprocs = 1
  end type simulation_config

contains

  subroutine init_simulation_config(cfg)
    type(simulation_config), intent(out) :: cfg
    cfg = simulation_config()
    call update_derived_config(cfg)
  end subroutine init_simulation_config

  subroutine update_derived_config(cfg)
    type(simulation_config), intent(inout) :: cfg

    cfg%lx = cfg%x_max - cfg%x_min
    cfg%ly = cfg%y_max - cfg%y_min
    cfg%lz = cfg%z_max - cfg%z_min

    if (cfg%nx > 0) cfg%dx = cfg%lx / real(cfg%nx, dp)
    if (cfg%ny > 0) cfg%dy = cfg%ly / real(cfg%ny, dp)
    if (cfg%nz > 0) cfg%dz = cfg%lz / real(cfg%nz, dp)
  end subroutine update_derived_config

  logical function should_output(cfg, step) result(flag)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step

    if (step == 0) then
      flag = cfg%write_initial
    else if (cfg%output_frequency > 0) then
      flag = (mod(step, cfg%output_frequency) == 0)
    else
      flag = .false.
    end if
  end function should_output

  subroutine print_simulation_config(cfg, unit)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in), optional :: unit
    integer :: u

    u = 6
    if (present(unit)) u = unit

    write(u,'(A)') '--- simulation_config ---'
    write(u,'(A,A)')    'equation          = ', trim(cfg%equation)
    write(u,'(A,A)')    'case_name         = ', trim(cfg%case_name)
    write(u,'(A,3I10)') 'nx, ny, nz        = ', cfg%nx, cfg%ny, cfg%nz
    write(u,'(A,3ES16.8)') 'lx, ly, lz     = ', cfg%lx, cfg%ly, cfg%lz
    write(u,'(A,3ES16.8)') 'dx, dy, dz     = ', cfg%dx, cfg%dy, cfg%dz
    write(u,'(A,ES16.8)') 'dt               = ', cfg%dt
    write(u,'(A,ES16.8)') 't_max            = ', cfg%t_max
    write(u,'(A,I10)')  'nsteps            = ', cfg%nsteps
    write(u,'(A,I10)')  'output_frequency  = ', cfg%output_frequency
    write(u,'(A,A)')    'output_dir        = ', trim(cfg%output_dir)
    write(u,'(A,A)')    'backend           = ', trim(cfg%backend)
  end subroutine print_simulation_config

end module mod_common_config
