!> GP3DのSLFバイナリ出力、meta.json、簡易密度スライスを担当する。
!> MPI時はrank別SLFと担当範囲を記録し、後処理ツールが全領域を再構成できるようにする。
module gp3d_io
  use, intrinsic :: iso_fortran_env, only: int32
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_state_t
  use gp3d_restart, only: gp3d_restart_info_t, gp3d_restart_load
  implicit none
  private

  public :: gp3d_output_config_t
  public :: gp3d_output_config_from_grid
  public :: gp3d_write_meta_json
  public :: gp3d_make_step_filename
  public :: gp3d_write_field_slf
  public :: gp3d_write_field_real3_slf
  public :: gp3d_write_field_real4_slf
  public :: gp3d_write_field_complex3_slf
  public :: gp3d_write_gpe_psi_slf
  public :: gp3d_read_gpe_psi_slf
  public :: gp3d_write_density_slice

  interface gp3d_write_field_slf
    module procedure gp3d_write_field_real3_slf
    module procedure gp3d_write_field_real4_slf
    module procedure gp3d_write_field_complex3_slf
  end interface gp3d_write_field_slf

  type :: gp3d_output_config_t
    character(len=256) :: output_dir = "output"
    character(len=64) :: equation = "GPE"
    character(len=128) :: case_name = "gp3d"
    character(len=32) :: precision_name = "float64"
    character(len=32) :: output_format = "slf"
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: nghost = 0
    real(dp) :: lx = 0.0_dp
    real(dp) :: ly = 0.0_dp
    real(dp) :: lz = 0.0_dp
    real(dp) :: dx = 0.0_dp
    real(dp) :: dy = 0.0_dp
    real(dp) :: dz = 0.0_dp
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 0.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 0.0_dp
    real(dp) :: z_min = 0.0_dp
    real(dp) :: z_max = 0.0_dp
    logical :: use_mpi = .false.
    logical :: use_openmp = .false.
    logical :: use_cuda = .false.
    integer :: cuda_device = 0
    integer :: mpi_rank = 0
    integer :: mpi_nprocs = 1
    integer :: local_nz = 0
    integer :: k_start = 1
    integer :: k_end = 0
  end type gp3d_output_config_t

contains

  subroutine gp3d_output_config_from_grid(cfg, grid, output_dir, case_name)
    type(gp3d_output_config_t), intent(out) :: cfg
    type(gp3d_grid_t), intent(in) :: grid
    character(len=*), intent(in), optional :: output_dir
    character(len=*), intent(in), optional :: case_name

    cfg%nx = grid%nx
    cfg%ny = grid%ny
    cfg%nz = grid%nz
    cfg%local_nz = grid%local_nz
    cfg%k_start = grid%k_start
    cfg%k_end = grid%k_end
    cfg%lx = grid%lx
    cfg%ly = grid%ly
    cfg%lz = grid%lz
    cfg%dx = grid%dx
    cfg%dy = grid%dy
    cfg%dz = grid%dz
    cfg%x_min = minval(grid%x)
    cfg%x_max = maxval(grid%x)
    cfg%y_min = minval(grid%y)
    cfg%y_max = maxval(grid%y)
    cfg%z_min = minval(grid%z)
    cfg%z_max = maxval(grid%z)

    if (present(output_dir)) cfg%output_dir = trim(output_dir)
    if (present(case_name)) cfg%case_name = trim(case_name)
  end subroutine gp3d_output_config_from_grid

  subroutine gp3d_write_meta_json(cfg)
    ! 格子・座標・変数・並列分割を後処理向けJSONへ一度だけ記録する。
    type(gp3d_output_config_t), intent(in) :: cfg

    integer :: unit, ios, rank, k_start, k_end
    character(len=512) :: filename

    call ensure_directory(cfg%output_dir)
    filename = trim(cfg%output_dir) // "/meta.json"

    open(newunit=unit, file=trim(filename), status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "ERROR: cannot write meta.json."

    write(unit,'(A)') "{"
    write(unit,'(A,A,A)') '  "equation": "', trim(cfg%equation), '",'
    write(unit,'(A,A,A)') '  "case_name": "', trim(cfg%case_name), '",'
    write(unit,'(A,I0,A,I0,A,I0,A)') '  "grid": [', cfg%nx, ', ', cfg%ny, ', ', cfg%nz, '],'
    write(unit,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') &
      '  "domain_length": [', cfg%lx, ', ', cfg%ly, ', ', cfg%lz, '],'
    write(unit,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') &
      '  "origin": [', cfg%x_min, ', ', cfg%y_min, ', ', cfg%z_min, '],'
    write(unit,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') &
      '  "spacing": [', cfg%dx, ', ', cfg%dy, ', ', cfg%dz, '],'
    write(unit,'(A,A,A)') '  "precision": "', trim(cfg%precision_name), '",'
    write(unit,'(A,A,A)') '  "format": "', trim(cfg%output_format), '",'
    write(unit,'(A)') '  "primary_variables_note": "GPE output stores psi_real and psi_imag only.",'
    write(unit,'(A)') '  "nvar": 2,'
    write(unit,'(A)') '  "variables": ["psi_real", "psi_imag"],'
    write(unit,'(A)') '  "parallel": {'
    write(unit,'(A,A,A)') '    "mpi_enabled": ', json_bool(cfg%use_mpi), ','
    write(unit,'(A,I0,A)') '    "mpi_nprocs": ', cfg%mpi_nprocs, ','
    write(unit,'(A,A,A)') '    "openmp_enabled": ', json_bool(cfg%use_openmp), ','
    write(unit,'(A)') '    "openmp_max_threads": 1,'
    write(unit,'(A,A,A)') '    "cuda_enabled": ', json_bool(cfg%use_cuda), ','
    write(unit,'(A,I0,A)') '    "cuda_device": ', cfg%cuda_device, ','
    if (cfg%use_mpi) then
      write(unit,'(A)') '    "decomposition": "z-slab-distributed-fft",'
    else
      write(unit,'(A)') '    "decomposition": "serial-global",'
    end if
    write(unit,'(A)') '    "rank_ranges": ['
    do rank = 0, cfg%mpi_nprocs - 1
      call block_range(cfg%nz, rank, cfg%mpi_nprocs, k_start, k_end)
      write(unit,'(A)') '      {'
      write(unit,'(A,I0,A)') '        "rank": ', rank, ','
      write(unit,'(A)') '        "i_start": 1,'
      write(unit,'(A,I0,A)') '        "i_end": ', cfg%nx, ','
      write(unit,'(A)') '        "j_start": 1,'
      write(unit,'(A,I0,A)') '        "j_end": ', cfg%ny, ','
      write(unit,'(A,I0,A)') '        "k_start": ', k_start, ','
      write(unit,'(A,I0)') '        "k_end": ', k_end
      if (rank < cfg%mpi_nprocs - 1) then
        write(unit,'(A)') '      },'
      else
        write(unit,'(A)') '      }'
      end if
    end do
    write(unit,'(A)') '    ]'
    write(unit,'(A)') '  }'
    write(unit,'(A)') "}"
    close(unit)
  end subroutine gp3d_write_meta_json

  subroutine gp3d_make_step_filename(cfg, step, rank, ext, filename)
    type(gp3d_output_config_t), intent(in) :: cfg
    integer, intent(in) :: step
    integer, intent(in), optional :: rank
    character(len=*), intent(in) :: ext
    character(len=*), intent(out) :: filename

    if (present(rank)) then
      write(filename,'(A,A,I0.6,A,I0.5,A,A)') trim(cfg%output_dir), "/field_", step, "_rank", rank, ".", trim(ext)
    else
      write(filename,'(A,A,I0.6,A,A)') trim(cfg%output_dir), "/field_", step, ".", trim(ext)
    end if
  end subroutine gp3d_make_step_filename

  subroutine gp3d_write_field_real4_slf(cfg, step, time, field, variable_names, rank)
    ! 複数実数変数を共通SLFヘッダーとFortran配列順でstream出力する。
    type(gp3d_output_config_t), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    real(dp), intent(in) :: field(:,:,:,:)
    character(len=*), intent(in), optional :: variable_names(:)
    integer, intent(in), optional :: rank

    integer :: unit, ios, nvar
    integer :: shape3(3)
    character(len=512) :: filename

    call ensure_directory(cfg%output_dir)
    call gp3d_make_step_filename(cfg, step, rank, "slf", filename)

    shape3 = [size(field, 1), size(field, 2), size(field, 3)]
    nvar = size(field, 4)

    open(newunit=unit, file=trim(filename), access="stream", form="unformatted", &
      status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "ERROR: cannot open SLF file."

    call write_header(unit, cfg, step, time, nvar, shape3, rank, variable_names)
    write(unit) field
    close(unit)
  end subroutine gp3d_write_field_real4_slf

  subroutine gp3d_write_field_real3_slf(cfg, step, time, field, variable_name, rank)
    type(gp3d_output_config_t), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    real(dp), intent(in) :: field(:,:,:)
    character(len=*), intent(in), optional :: variable_name
    integer, intent(in), optional :: rank

    real(dp), allocatable :: tmp(:,:,:,:)
    character(len=32) :: names(1)

    allocate(tmp(size(field, 1), size(field, 2), size(field, 3), 1))
    tmp(:,:,:,1) = field
    names(1) = "var1"
    if (present(variable_name)) names(1) = trim(variable_name)
    call gp3d_write_field_real4_slf(cfg, step, time, tmp, names, rank)
    deallocate(tmp)
  end subroutine gp3d_write_field_real3_slf

  subroutine gp3d_write_field_complex3_slf(cfg, step, time, psi, rank)
    type(gp3d_output_config_t), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    complex(dp), intent(in) :: psi(:,:,:)
    integer, intent(in), optional :: rank

    real(dp), allocatable :: tmp(:,:,:,:)
    character(len=32) :: names(2)

    allocate(tmp(size(psi, 1), size(psi, 2), size(psi, 3), 2))
    tmp(:,:,:,1) = real(psi, dp)
    tmp(:,:,:,2) = aimag(psi)
    names = [character(len=32) :: "psi_real", "psi_imag"]

    call gp3d_write_field_real4_slf(cfg, step, time, tmp, names, rank)
    deallocate(tmp)
  end subroutine gp3d_write_field_complex3_slf

  subroutine gp3d_write_gpe_psi_slf(cfg, step, time, state, rank)
    ! 複素psiを可搬なpsi_real/psi_imagの2変数へ分解して保存する。
    type(gp3d_output_config_t), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    type(gp3d_state_t), intent(in) :: state
    integer, intent(in), optional :: rank

    call gp3d_write_field_complex3_slf(cfg, step, time, state%psi, rank)
  end subroutine gp3d_write_gpe_psi_slf

  subroutine gp3d_read_gpe_psi_slf(filename, state, grid)
    character(len=*), intent(in) :: filename
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid

    type(gp3d_restart_info_t) :: info

    call gp3d_restart_load(filename, state, grid, info)
  end subroutine gp3d_read_gpe_psi_slf

  subroutine gp3d_write_density_slice(filename, state, grid, iz)
    character(len=*), intent(in) :: filename
    type(gp3d_state_t), intent(in) :: state
    type(gp3d_grid_t), intent(in) :: grid
    integer, intent(in) :: iz

    integer :: unit, i, j, local_iz

    if (iz < 1 .or. iz > grid%nz) error stop "slice index is out of range"
    if (iz < grid%k_start .or. iz > grid%k_end) error stop "slice is not owned by this MPI rank"
    local_iz = iz - grid%k_start + 1

    open(newunit=unit, file=filename, status="replace", action="write", form="formatted")
    write(unit, '(a)') "# x y density"
    do j = 1, grid%ny
      do i = 1, grid%nx
        write(unit, '(3(es24.16,1x))') grid%x(i), grid%y(j), abs(state%psi(i,j,local_iz))**2
      end do
      write(unit, *)
    end do
    close(unit)
  end subroutine gp3d_write_density_slice

  subroutine write_header(unit, cfg, step, time, nvar, shape3, rank, varnames)
    integer, intent(in) :: unit
    type(gp3d_output_config_t), intent(in) :: cfg
    integer, intent(in) :: step, nvar
    real(dp), intent(in) :: time
    integer, intent(in) :: shape3(3)
    integer, intent(in), optional :: rank
    character(len=*), intent(in), optional :: varnames(:)

    character(len=8) :: magic
    integer(int32) :: version, dtype_code, ndim, nvar_i4
    integer(int32) :: shp(4), meta(8)
    character(len=32) :: name32
    integer :: i

    magic = "SLF1" // char(0) // char(0) // char(0) // char(0)
    version = 1_int32
    dtype_code = 2_int32
    ndim = 4_int32
    nvar_i4 = int(nvar, int32)
    shp = [int(shape3(1), int32), int(shape3(2), int32), int(shape3(3), int32), int(nvar, int32)]

    meta = 0_int32
    meta(1) = int(step, int32)
    if (present(rank)) meta(2) = int(rank, int32)
    meta(3) = int(cfg%nx, int32)
    meta(4) = int(cfg%ny, int32)
    meta(5) = int(cfg%nz, int32)
    meta(6) = int(cfg%nghost, int32)
    meta(7) = int(cfg%mpi_nprocs, int32)
    meta(8) = int(cfg%k_start, int32)

    write(unit) magic
    write(unit) version
    write(unit) dtype_code
    write(unit) ndim
    write(unit) shp
    write(unit) meta
    write(unit) time
    write(unit) cfg%x_min, cfg%x_max, cfg%y_min, cfg%y_max, cfg%z_min, cfg%z_max
    write(unit) nvar_i4

    do i = 1, nvar
      name32 = ""
      if (present(varnames)) then
        if (i <= size(varnames)) name32 = trim(varnames(i))
      end if
      if (len_trim(name32) == 0) write(name32,'(A,I0)') "var", i
      write(unit) name32
    end do
  end subroutine write_header

  subroutine ensure_directory(dirname)
    character(len=*), intent(in) :: dirname
    character(len=768) :: command
    character(len=64) :: os_name
    integer :: status

    if (len_trim(dirname) == 0 .or. trim(dirname) == ".") return
    os_name = ""
    call get_environment_variable("OS", os_name)
    if (trim(os_name) == "Windows_NT") then
      command = 'if not exist "' // trim(dirname) // '" mkdir "' // trim(dirname) // '"'
    else
      command = 'mkdir -p "' // trim(dirname) // '"'
    end if
    call execute_command_line(trim(command), wait=.true., exitstat=status)
    if (status /= 0) error stop "ERROR: cannot create output directory."
  end subroutine ensure_directory

  pure function json_bool(flag) result(str)
    logical, intent(in) :: flag
    character(len=5) :: str

    if (flag) then
      str = "true "
    else
      str = "false"
    end if
  end function json_bool

  pure subroutine block_range(n, rank, nprocs, start_index, end_index)
    integer, intent(in) :: n, rank, nprocs
    integer, intent(out) :: start_index, end_index
    integer :: base, rest, count

    base = n / nprocs
    rest = mod(n, nprocs)
    if (rank < rest) then
      count = base + 1
      start_index = rank * count + 1
    else
      count = base
      start_index = rest * (base + 1) + (rank - rest) * base + 1
    end if
    end_index = start_index + count - 1
  end subroutine block_range

end module gp3d_io
