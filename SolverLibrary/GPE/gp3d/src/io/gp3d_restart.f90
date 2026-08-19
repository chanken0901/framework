!> 単一またはrank分割されたSLFから波動関数を復元し、途中計算を再開する。
!> 保存時と再開時のMPI分割数が異なっても、各rankの重なるz範囲だけを読み込む。
module gp3d_restart
  use, intrinsic :: iso_fortran_env, only: error_unit, int32, int64
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_run_config_t, gp3d_state_t
  implicit none
  private

  type, public :: gp3d_restart_info_t
    logical :: loaded = .false.
    integer :: step = 0
    real(dp) :: time = 0.0_dp
    integer :: source_nprocs = 1
    character(len=512) :: source_file = ""
  end type gp3d_restart_info_t

  type :: slf_header_t
    integer :: step = 0
    integer :: rank = 0
    integer :: source_nprocs = 1
    integer :: global_nx = 0
    integer :: global_ny = 0
    integer :: global_nz = 0
    integer :: file_nx = 0
    integer :: file_ny = 0
    integer :: file_nz = 0
    integer :: k_start = 1
    integer :: nvar = 0
    real(dp) :: time = 0.0_dp
    real(dp) :: bounds(6) = 0.0_dp
    integer(int64) :: data_pos = 0_int64
    character(len=32), allocatable :: names(:)
  end type slf_header_t

  public :: gp3d_restart_requested
  public :: gp3d_restart_load

contains

  pure logical function gp3d_restart_requested(run_cfg) result(requested)
    type(gp3d_run_config_t), intent(in) :: run_cfg

    requested = len_trim(run_cfg%restart_file) > 0
    if (requested) return
    select case (trim(run_cfg%initial_condition))
    case ("restart_slf", "restart", "slf")
      requested = .true.
    case default
      requested = .false.
    end select
  end function gp3d_restart_requested

  subroutine gp3d_restart_load(filename, state, grid, info)
    ! 入力がrankファミリーか単一SLFかを判定し、現在の局所スラブを完全に埋める。
    character(len=*), intent(in) :: filename
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    type(gp3d_restart_info_t), intent(out) :: info

    type(slf_header_t) :: first_header, header
    character(len=512) :: source_file, rank_prefix, rank_suffix, rank_file
    logical, allocatable :: covered(:)
    logical :: is_rank_family
    integer :: source_rank

    info = gp3d_restart_info_t()
    if (len_trim(filename) == 0) then
      call restart_error("restart_file must be specified")
    end if
    if (.not. allocated(state%psi)) then
      call restart_error("restart destination wave function is not allocated")
    end if
    if (size(state%psi, 1) /= grid%nx .or. size(state%psi, 2) /= grid%ny .or. &
        size(state%psi, 3) /= grid%local_nz) then
      call restart_error("restart destination shape does not match the local grid")
    end if

    call resolve_source_file(filename, source_file)
    call read_slf_header(source_file, first_header)
    call validate_header(first_header, grid)

    allocate(covered(grid%local_nz), source=.false.)
    is_rank_family = first_header%file_nz /= first_header%global_nz .or. &
      first_header%k_start /= 1
    if (is_rank_family) then
      call split_rank_filename(source_file, rank_prefix, rank_suffix, is_rank_family)
      if (.not. is_rank_family) then
        call restart_error("a slab SLF restart file must use the _rankNNNNN filename convention")
      end if
      do source_rank = 0, first_header%source_nprocs - 1
        call make_rank_filename(rank_prefix, rank_suffix, source_rank, rank_file)
        call read_slf_header(rank_file, header)
        call validate_header(header, grid)
        call validate_family_member(first_header, header, source_rank)
        call read_overlap(rank_file, header, state, grid, covered)
      end do
    else
      call read_overlap(source_file, first_header, state, grid, covered)
    end if

    if (.not. all(covered)) then
      call restart_error("restart SLF files do not cover the current rank's complete z slab")
    end if

    info%loaded = .true.
    info%step = first_header%step
    info%time = first_header%time
    info%source_nprocs = first_header%source_nprocs
    info%source_file = trim(source_file)
  end subroutine gp3d_restart_load

  subroutine resolve_source_file(requested_file, source_file)
    character(len=*), intent(in) :: requested_file
    character(len=*), intent(out) :: source_file

    character(len=512) :: rank_zero_file
    logical :: exists
    integer :: extension_pos, length

    source_file = trim(requested_file)
    inquire(file=trim(source_file), exist=exists)
    if (exists) return

    length = len_trim(requested_file)
    extension_pos = index(requested_file(1:length), ".slf", back=.true.)
    if (extension_pos > 0 .and. extension_pos + 3 == length) then
      rank_zero_file = requested_file(1:extension_pos - 1) // "_rank00000.slf"
    else
      rank_zero_file = trim(requested_file) // "_rank00000.slf"
    end if
    inquire(file=trim(rank_zero_file), exist=exists)
    if (.not. exists) then
      call restart_error("restart SLF file was not found: " // trim(requested_file))
    end if
    source_file = trim(rank_zero_file)
  end subroutine resolve_source_file

  subroutine read_slf_header(filename, header)
    character(len=*), intent(in) :: filename
    type(slf_header_t), intent(out) :: header

    character(len=8) :: magic
    integer(int32) :: version, dtype_code, ndim, shape(4), meta(8), nvar_i4
    integer :: unit, ios, i

    open(newunit=unit, file=trim(filename), access="stream", form="unformatted", &
      status="old", action="read", iostat=ios)
    if (ios /= 0) call restart_error("cannot open restart SLF file: " // trim(filename))

    read(unit, iostat=ios) magic
    if (ios == 0) read(unit, iostat=ios) version
    if (ios == 0) read(unit, iostat=ios) dtype_code
    if (ios == 0) read(unit, iostat=ios) ndim
    if (ios == 0) read(unit, iostat=ios) shape
    if (ios == 0) read(unit, iostat=ios) meta
    if (ios == 0) read(unit, iostat=ios) header%time
    if (ios == 0) read(unit, iostat=ios) header%bounds
    if (ios == 0) read(unit, iostat=ios) nvar_i4
    if (ios /= 0) then
      close(unit)
      call restart_error("cannot read restart SLF header: " // trim(filename))
    end if

    if (magic(1:4) /= "SLF1") call restart_error("restart file is not SLF1: " // trim(filename))
    if (version /= 1_int32) call restart_error("unsupported restart SLF version")
    if (dtype_code /= 2_int32) call restart_error("restart SLF must store float64 data")
    if (ndim /= 4_int32) call restart_error("restart SLF must contain rank-4 field data")
    if (nvar_i4 <= 0_int32 .or. nvar_i4 > 1024_int32) call restart_error("invalid restart SLF variable count")
    if (shape(4) /= nvar_i4) call restart_error("restart SLF shape and variable count disagree")

    header%step = int(meta(1))
    header%rank = int(meta(2))
    header%global_nx = int(meta(3))
    header%global_ny = int(meta(4))
    header%global_nz = int(meta(5))
    header%source_nprocs = max(1, int(meta(7)))
    header%k_start = max(1, int(meta(8)))
    header%file_nx = int(shape(1))
    header%file_ny = int(shape(2))
    header%file_nz = int(shape(3))
    header%nvar = int(nvar_i4)
    if (header%global_nx <= 0) header%global_nx = header%file_nx
    if (header%global_ny <= 0) header%global_ny = header%file_ny
    if (header%global_nz <= 0) header%global_nz = header%file_nz

    allocate(header%names(header%nvar))
    do i = 1, header%nvar
      read(unit, iostat=ios) header%names(i)
      if (ios /= 0) then
        close(unit)
        call restart_error("cannot read restart SLF variable names")
      end if
      header%names(i) = adjustl(header%names(i))
    end do
    inquire(unit=unit, pos=header%data_pos)
    close(unit)
  end subroutine read_slf_header

  subroutine validate_header(header, grid)
    ! 精度、変数名、全体格子など、誤った再開を防ぐための互換性を確認する。
    type(slf_header_t), intent(in) :: header
    type(gp3d_grid_t), intent(in) :: grid

    real(dp) :: expected_bounds(6), tolerance, scale

    if (header%global_nx /= grid%nx .or. header%global_ny /= grid%ny .or. &
        header%global_nz /= grid%nz) then
      call restart_error("restart SLF global grid does not match the configured grid")
    end if
    if (header%file_nx /= grid%nx .or. header%file_ny /= grid%ny) then
      call restart_error("restart SLF x-y shape does not match the configured grid")
    end if
    if (header%file_nz <= 0 .or. header%k_start + header%file_nz - 1 > grid%nz) then
      call restart_error("restart SLF z-slab metadata is invalid")
    end if
    if (find_name(header%names, "psi_real") <= 0 .or. &
        find_name(header%names, "psi_imag") <= 0) then
      call restart_error("restart SLF requires psi_real and psi_imag variables")
    end if

    expected_bounds = [grid%x(1), grid%x(grid%nx), grid%y(1), grid%y(grid%ny), &
      grid%z(1), grid%z(grid%nz)]
    scale = max(1.0_dp, max(maxval(abs(expected_bounds)), maxval(abs(header%bounds))))
    tolerance = 1.0e-10_dp * scale
    if (maxval(abs(expected_bounds - header%bounds)) > tolerance) then
      call restart_error("restart SLF domain bounds do not match the configured grid")
    end if
  end subroutine validate_header

  subroutine validate_family_member(first, member, expected_rank)
    type(slf_header_t), intent(in) :: first, member
    integer, intent(in) :: expected_rank

    real(dp) :: time_tolerance

    time_tolerance = 1.0e-12_dp * max(1.0_dp, abs(first%time))
    if (member%rank /= expected_rank) call restart_error("restart rank filename and SLF rank metadata disagree")
    if (member%step /= first%step .or. abs(member%time - first%time) > time_tolerance) then
      call restart_error("restart rank files do not have a common step and time")
    end if
    if (member%source_nprocs /= first%source_nprocs) then
      call restart_error("restart rank files disagree on the source process count")
    end if
    if (member%global_nx /= first%global_nx .or. member%global_ny /= first%global_ny .or. &
        member%global_nz /= first%global_nz) then
      call restart_error("restart rank files disagree on the global grid")
    end if
  end subroutine validate_family_member

  subroutine read_overlap(filename, header, state, grid, covered)
    ! ファイルのz範囲と現在rankのz範囲の共通部分だけをpsiへコピーする。
    character(len=*), intent(in) :: filename
    type(slf_header_t), intent(in) :: header
    type(gp3d_state_t), intent(inout) :: state
    type(gp3d_grid_t), intent(in) :: grid
    logical, intent(inout) :: covered(:)

    real(dp), allocatable :: real_part(:,:,:), imag_part(:,:,:)
    integer :: overlap_start, overlap_end, planes, destination_start
    integer :: real_index, imag_index, unit, ios
    integer(int64) :: values_per_variable, plane_offset, real_pos, imag_pos

    overlap_start = max(grid%k_start, header%k_start)
    overlap_end = min(grid%k_end, header%k_start + header%file_nz - 1)
    if (overlap_end < overlap_start) return

    planes = overlap_end - overlap_start + 1
    destination_start = overlap_start - grid%k_start + 1
    real_index = find_name(header%names, "psi_real")
    imag_index = find_name(header%names, "psi_imag")
    values_per_variable = int(header%file_nx, int64) * int(header%file_ny, int64) * &
      int(header%file_nz, int64)
    plane_offset = int(header%file_nx, int64) * int(header%file_ny, int64) * &
      int(overlap_start - header%k_start, int64)
    real_pos = header%data_pos + 8_int64 * &
      (int(real_index - 1, int64) * values_per_variable + plane_offset)
    imag_pos = header%data_pos + 8_int64 * &
      (int(imag_index - 1, int64) * values_per_variable + plane_offset)

    allocate(real_part(grid%nx, grid%ny, planes))
    allocate(imag_part(grid%nx, grid%ny, planes))
    open(newunit=unit, file=trim(filename), access="stream", form="unformatted", &
      status="old", action="read", iostat=ios)
    if (ios /= 0) call restart_error("cannot reopen restart SLF data: " // trim(filename))
    read(unit, pos=real_pos, iostat=ios) real_part
    if (ios == 0) read(unit, pos=imag_pos, iostat=ios) imag_part
    close(unit)
    if (ios /= 0) call restart_error("cannot read restart SLF field data: " // trim(filename))

    state%psi(:,:,destination_start:destination_start + planes - 1) = &
      cmplx(real_part, imag_part, kind=dp)
    covered(destination_start:destination_start + planes - 1) = .true.
  end subroutine read_overlap

  subroutine split_rank_filename(filename, prefix, suffix, matched)
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: prefix, suffix
    logical, intent(out) :: matched

    integer :: marker, i, length

    prefix = ""
    suffix = ""
    matched = .false.
    length = len_trim(filename)
    marker = index(filename(1:length), "_rank", back=.true.)
    if (marker <= 0 .or. marker + 9 > length) return
    do i = marker + 5, marker + 9
      if (filename(i:i) < "0" .or. filename(i:i) > "9") return
    end do
    prefix = filename(1:marker + 4)
    if (marker + 10 <= length) suffix = filename(marker + 10:length)
    matched = .true.
  end subroutine split_rank_filename

  subroutine make_rank_filename(prefix, suffix, rank, filename)
    character(len=*), intent(in) :: prefix, suffix
    integer, intent(in) :: rank
    character(len=*), intent(out) :: filename

    write(filename,'(A,I5.5,A)') trim(prefix), rank, trim(suffix)
  end subroutine make_rank_filename

  pure integer function find_name(names, target) result(index_value)
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: target
    integer :: i

    index_value = 0
    do i = 1, size(names)
      if (trim(names(i)) == trim(target)) then
        index_value = i
        return
      end if
    end do
  end function find_name

  subroutine restart_error(message)
    character(len=*), intent(in) :: message

    write(error_unit,'(A)') "ERROR: " // trim(message)
    error stop
  end subroutine restart_error

end module gp3d_restart
