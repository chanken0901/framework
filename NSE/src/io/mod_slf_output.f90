module mod_slf_output
  use, intrinsic :: iso_fortran_env, only : int32
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use module_mpi
  implicit none
  private

  public :: write_meta_json
  public :: write_field_slf
  public :: write_nse_conserved_slf
  public :: write_gpe_psi_slf
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

    ! Works on Linux/macOS/Git Bash/MSYS2.  If this is not available on a
    ! target system, create the output directory before running the solver.
    write(cmd,'(A,A,A)') 'mkdir -p "', trim(dirname), '"'
    !call execute_command_line(trim(cmd), wait=.true.)
    call execute_command_line('if not exist "' // trim(dirname) // '" mkdir "' // trim(dirname) // '"')
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

  subroutine write_meta_json(cfg, variable_names, filename, primary_variables_note)
    type(simulation_config), intent(in) :: cfg
    character(len=*), intent(in), optional :: variable_names(:)
    character(len=*), intent(in), optional :: filename
    character(len=*), intent(in), optional :: primary_variables_note

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
    if (present(primary_variables_note)) then
      write(u,'(A,A,A)') '  "primary_variables_note": "', trim(primary_variables_note), '",'
    end if
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

  subroutine write_header(u, cfg, step, time, nvar, shape3, rank, varnames, dtype_code_in)
    integer, intent(in) :: u
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step, nvar
    real(dp), intent(in) :: time
    integer, intent(in) :: shape3(3)
    integer, intent(in), optional :: rank
    character(len=*), intent(in), optional :: varnames(:)
    integer, intent(in), optional :: dtype_code_in

    character(len=8) :: magic
    integer(int32) :: version, dtype_code, ndim, nvar_i4
    integer(int32) :: shp(4), meta(8)
    character(len=32) :: name32
    integer :: i

    magic = 'SLF1'//char(0)//char(0)//char(0)//char(0)
    version = 1_int32
    dtype_code = 2_int32       ! 2 = float64. Complex values are stored as two float64 variables.
    if (present(dtype_code_in)) dtype_code = int(dtype_code_in, int32)
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
    meta(7) = int(nprocs, int32)

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
    ! Generic writer for real primary variables stored as field(nx,ny,nz,nvar).
    ! For NSE, pass the conservative-variable array directly, e.g. Q(:,:,:,:).
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
    ! Generic writer for a single real scalar field.
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
    ! Generic complex-field writer.
    ! Important: this stores only the primary complex field as two real arrays.
    ! It does NOT compute rho or phase.  Those should be derived in post-processing.
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    complex(dp), intent(in) :: psi(:,:,:)
    integer, intent(in), optional :: rank

    real(dp), allocatable :: tmp(:,:,:,:)
    character(len=32) :: names(2)

    allocate(tmp(size(psi,1), size(psi,2), size(psi,3), 2))
    tmp(:,:,:,1) = real(psi, dp)
    tmp(:,:,:,2) = aimag(psi)

    names = [character(len=32) :: 'psi_real', 'psi_imag']
    call write_field_real4_slf(cfg, step, time, tmp, names, rank)
    deallocate(tmp)
  end subroutine write_field_complex3_slf

  subroutine write_nse_conserved_slf(cfg, step, time, q, rank)
    ! NSE standard output: conservative variables only.
    ! q(:,:,:,1:5) = [rho, rho_u, rho_v, rho_w, rho_E]
    ! Primitive variables such as u,v,w,p,T are intentionally not written here.
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    real(dp), intent(in) :: q(:,:,:,:)
    integer, intent(in), optional :: rank

    character(len=32), allocatable :: names(:)
    integer :: nvar, ivar

    nvar = size(q,4)
    allocate(names(nvar))

    if (nvar >= 5) then
      names(1:5) = [character(len=32) :: 'rho', 'rho_u', 'rho_v', 'rho_w', 'rho_E']
      do ivar = 6, nvar
        write(names(ivar),'(A,I0)') 'q', ivar
      end do
    else
      do ivar = 1, nvar
        write(names(ivar),'(A,I0)') 'q', ivar
      end do
    end if

    call write_field_real4_slf(cfg, step, time, q, names, rank)
    deallocate(names)
  end subroutine write_nse_conserved_slf

  subroutine write_gpe_psi_slf(cfg, step, time, psi, rank)
    ! GPE standard output: wave function only.
    ! psi is stored as [psi_real, psi_imag].  rho=|psi|^2 and phase=atan2(Im,Re)
    ! should be computed by the post-processing scripts.
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: step
    real(dp), intent(in) :: time
    complex(dp), intent(in) :: psi(:,:,:)
    integer, intent(in), optional :: rank

    call write_field_complex3_slf(cfg, step, time, psi, rank)
  end subroutine write_gpe_psi_slf

end module mod_slf_output
