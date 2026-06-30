module mod_slf_output
  use, intrinsic :: iso_fortran_env, only : int32
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  implicit none
  private

  public :: write_meta_json
  public :: write_field_slf
  public :: write_field_real3_slf
  public :: write_field_real4_slf
  public :: write_field_complex3_slf
  public :: make_step_filename

  interface write_field_slf
    module procedure write_field_real3_slf
    module procedure write_field_real4_slf
    module procedure write_field_complex3_slf
  end interface write_field_slf

contains

  subroutine ensure_directory(dirname)
    character(len=*), intent(in) :: dirname
    character(len=512) :: cmd
    if (len_trim(dirname) == 0) return
    write(cmd,'(A,A,A)') 'mkdir -p "', trim(dirname), '"'
    call execute_command_line(trim(cmd), wait=.true.)
  end subroutine ensure_directory

  subroutine make_step_filename(cfg, step, rank, ext, fname)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    integer, intent(in), optional :: rank
    character(len=*), intent(in) :: ext
    character(len=*), intent(out) :: fname

    if (present(rank)) then
      write(fname,'(A,A,I0.6,A,I0.5,A,A)') trim(cfg%output_dir), '/field_', step, '_rank', rank, '.', trim(ext)
    else
      write(fname,'(A,A,I0.6,A,A)') trim(cfg%output_dir), '/field_', step, '.', trim(ext)
    end if
  end subroutine make_step_filename

  subroutine write_meta_json(cfg, variable_names, filename)
    type(simulation_config), intent(in) :: cfg
    character(len=*), intent(in), optional :: variable_names(:)
    character(len=*), intent(in), optional :: filename

    integer :: u, i, nvar, ios
    character(len=512) :: fname

    call ensure_directory(cfg%output_dir)
    if (present(filename)) then
      fname = filename
    else
      fname = trim(cfg%output_dir)//'/meta.json'
    end if

    nvar = 0
    if (present(variable_names)) nvar = size(variable_names)

    open(newunit=u, file=trim(fname), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot write meta.json.'

    write(u,'(A)') '{'
    write(u,'(A,A,A)') '  "equation": "', trim(cfg%equation), '",'
    write(u,'(A,A,A)') '  "case_name": "', trim(cfg%case_name), '",'
    write(u,'(A,I0,A,I0,A,I0,A)') '  "grid": [', cfg%nx, ', ', cfg%ny, ', ', cfg%nz, '],'
    write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "domain_length": [', cfg%lx, ', ', cfg%ly, ', ', cfg%lz, '],'
    write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "origin": [', cfg%x_min, ', ', cfg%y_min, ', ', cfg%z_min, '],'
    write(u,'(A,A,A)') '  "precision": "', trim(cfg%precision_name), '",'
    write(u,'(A,A,A)') '  "format": "', trim(cfg%output_format), '",'
    write(u,'(A)', advance='no') '  "variables": ['
    if (nvar > 0) then
      do i = 1, nvar
        if (i > 1) write(u,'(A)', advance='no') ', '
        write(u,'(A,A,A)', advance='no') '"', trim(variable_names(i)), '"'
      end do
    end if
    write(u,'(A)') ']'
    write(u,'(A)') '}'
    close(u)
  end subroutine write_meta_json

  subroutine write_header(u, cfg, step, time, nvar, shape3, rank, varnames)
    integer, intent(in) :: u
    type(simulation_config), intent(in) :: cfg
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

    magic = 'SLF1'//char(0)//char(0)//char(0)//char(0)
    version = 1_int32
    dtype_code = 2_int32    ! 2 = float64
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
    meta(7) = int(cfg%nprocs, int32)

    write(u) magic
    write(u) version
    write(u) dtype_code
    write(u) ndim
    write(u) shp
    write(u) meta
    write(u) time
    write(u) cfg%x_min, cfg%x_max, cfg%y_min, cfg%y_max, cfg%z_min, cfg%z_max
    write(u) nvar_i4

    do i = 1, nvar
      name32 = ''
      if (present(varnames)) then
        if (i <= size(varnames)) name32 = trim(varnames(i))
      end if
      if (len_trim(name32) == 0) write(name32,'(A,I0)') 'var', i
      write(u) name32
    end do
  end subroutine write_header

  subroutine write_field_real4_slf(cfg, step, time, field, variable_names, rank)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    real(dp), intent(in) :: field(:,:,:,:)
    character(len=*), intent(in), optional :: variable_names(:)
    integer, intent(in), optional :: rank

    integer :: u, ios, nvar
    character(len=512) :: fname
    integer :: shape3(3)

    call ensure_directory(cfg%output_dir)
    if (present(rank)) then
      call make_step_filename(cfg, step, rank, 'slf', fname)
    else
      call make_step_filename(cfg, step, ext='slf', fname=fname)
    end if

    shape3 = [size(field,1), size(field,2), size(field,3)]
    nvar = size(field,4)

    open(newunit=u, file=trim(fname), access='stream', form='unformatted', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot open SLF file.'

    call write_header(u, cfg, step, time, nvar, shape3, rank, variable_names)
    write(u) field
    close(u)
  end subroutine write_field_real4_slf

  subroutine write_field_real3_slf(cfg, step, time, field, variable_name, rank)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    real(dp), intent(in) :: field(:,:,:)
    character(len=*), intent(in), optional :: variable_name
    integer, intent(in), optional :: rank

    real(dp), allocatable :: tmp(:,:,:,:)
    character(len=32) :: names(1)

    allocate(tmp(size(field,1), size(field,2), size(field,3), 1))
    tmp(:,:,:,1) = field
    names(1) = 'var1'
    if (present(variable_name)) names(1) = trim(variable_name)
    call write_field_real4_slf(cfg, step, time, tmp, names, rank)
    deallocate(tmp)
  end subroutine write_field_real3_slf

  subroutine write_field_complex3_slf(cfg, step, time, psi, rank)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    complex(dp), intent(in) :: psi(:,:,:)
    integer, intent(in), optional :: rank

    real(dp), allocatable :: tmp(:,:,:,:)
    character(len=32) :: names(4)

    allocate(tmp(size(psi,1), size(psi,2), size(psi,3), 4))
    tmp(:,:,:,1) = real(psi, dp)
    tmp(:,:,:,2) = aimag(psi)
    tmp(:,:,:,3) = real(psi*conjg(psi), dp)
    tmp(:,:,:,4) = atan2(aimag(psi), real(psi, dp))

    names = [character(len=32) :: 'psi_real', 'psi_imag', 'rho', 'phase']
    call write_field_real4_slf(cfg, step, time, tmp, names, rank)
    deallocate(tmp)
  end subroutine write_field_complex3_slf

end module mod_slf_output
