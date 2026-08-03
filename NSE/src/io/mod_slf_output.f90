module mod_slf_output
  use, intrinsic :: iso_fortran_env, only : int32
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_openmp_runtime, only : nse_max_threads
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

  !subroutine write_meta_json(cfg)
  !  type(simulation_config), intent(in) :: cfg
  !
  !  integer :: u, i, nvar, ios
  subroutine write_meta_json(cfg, filename, primary_variables_note, is, ie, js, je, ks, ke, use_cuda)
    type(simulation_config), intent(in) :: cfg
    character(len=*), intent(in), optional :: filename
    character(len=*), intent(in), optional :: primary_variables_note
    character(len=32), allocatable :: variable_names(:)
    character(len=64) :: eq
    integer, intent(in), optional :: is, ie, js, je, ks, ke
    logical, intent(in), optional :: use_cuda
  
    integer :: u, i, nvar, ios, ierr, r
    integer :: local_range(6)
    integer, allocatable :: all_ranges(:,:)
    integer :: omp_threads
    logical :: cuda_enabled
    character(len=512) :: fname

    if (my_rank == root) call ensure_directory(cfg%output_dir)
    call mp_barrier
  
    fname = trim(cfg%output_dir)//'/meta.json'
    eq = adjustl(cfg%equation)

  ! local_range = [i_start, i_end, j_start, j_end, k_start, k_end]
    !local_range = [1, cfg%nx, 1, cfg%ny, 1, cfg%nz]
    !if (present(js)) local_range(3) = is
    !if (present(je)) local_range(4) = ie
    !if (present(js)) local_range(3) = js
    !if (present(je)) local_range(4) = je
    !if (present(ks)) local_range(5) = ks
    !if (present(ke)) local_range(6) = ke

local_range = [1, cfg%nx, 1, cfg%ny, 1, cfg%nz]

if (present(is)) local_range(1) = is
if (present(ie)) local_range(2) = ie
if (present(js)) local_range(3) = js
if (present(je)) local_range(4) = je
if (present(ks)) local_range(5) = ks
if (present(ke)) local_range(6) = ke
  
    ! MPI_Gather ignores the receive buffer on non-root ranks, but the
    ! Fortran actual argument must still be a valid allocated array.
    allocate(all_ranges(6,max(1,nprocs)))
    call MPI_Gather(local_range, 6, MPI_INTEGER, all_ranges, 6, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    if (my_rank /= root) return

    cuda_enabled = .false.
    if (present(use_cuda)) cuda_enabled = use_cuda
    omp_threads = 1
    if (cfg%use_openmp) omp_threads = nse_max_threads()

  
    select case (trim(eq))
    case ('NSE','nse')
      nvar = 5
      allocate(variable_names(nvar))
      variable_names = [character(len=32) :: &
        'rho', 'rho_u', 'rho_v', 'rho_w', 'rho_E']
  
    case ('GPE','gpe')
      nvar = 2
      allocate(variable_names(nvar))
      variable_names = [character(len=32) :: &
        'psi_real', 'psi_imag']
  
    case default
      nvar = 0
      allocate(variable_names(0))
    end select
  
    open(newunit=u, file=trim(fname), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'ERROR: cannot write meta.json.'
  
    write(u,'(A)') '{'
    write(u,'(A,A,A)') '  "equation": "', trim(cfg%equation), '",'
    write(u,'(A,A,A)') '  "case_name": "', trim(cfg%case_name), '",'
    write(u,'(A,I0,A,I0,A,I0,A)') '  "grid": [', cfg%nx, ', ', cfg%ny, ', ', cfg%nz, '],'
    write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "domain_length": [', cfg%lx, ', ', cfg%ly, ', ', cfg%lz, '],'
    write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "origin": [', cfg%x_min, ', ', cfg%y_min, ', ', cfg%z_min, '],'
    write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "spacing": [', cfg%dx, ', ', cfg%dy, ', ', cfg%dz, '],'
    write(u,'(A,A,A)') '  "precision": "', trim(cfg%precision_name), '",'
    write(u,'(A,A,A)') '  "format": "', trim(cfg%output_format), '",'
  
    select case (trim(eq))
    case ('NSE','nse')
      write(u,'(A)') '  "primary_variables_note": "NSE output stores conservative variables only.",'
    case ('GPE','gpe')
      write(u,'(A)') '  "primary_variables_note": "GPE output stores psi_real and psi_imag only.",'
    case default
      write(u,'(A)') '  "primary_variables_note": "Primary variables depend on the solver.",'
    end select
  
    write(u,'(A,I0,A)') '  "nvar": ', nvar, ','
    write(u,'(A)', advance='no') '  "variables": ['
  
    do i = 1, nvar
      if (i > 1) write(u,'(A)', advance='no') ', '
      write(u,'(A,A,A)', advance='no') '"', trim(variable_names(i)), '"'
    end do
  
    write(u,'(A)') '],'

    write(u,'(A)') '  "parallel": {'
    write(u,'(A,A,A)') '    "mpi_enabled": ', json_bool(cfg%use_mpi), ','
    write(u,'(A,I0,A)') '    "mpi_nprocs": ', nprocs, ','
    write(u,'(A,A,A)') '    "openmp_enabled": ', json_bool(cfg%use_openmp), ','
    write(u,'(A,I0,A)') '    "openmp_max_threads": ', omp_threads, ','
    write(u,'(A,A,A)') '    "cuda_enabled": ', json_bool(cuda_enabled), ','
    write(u,'(A)') '    "decomposition": "x-global_yz-block",'
    write(u,'(A)') '    "rank_ranges": ['

do r = 0, nprocs-1
  write(u,'(A)') '      {'
  write(u,'(A)') '        "rank": '    // trim(itoa(r)) // ','
  write(u,'(A)') '        "i_start": ' // trim(itoa(all_ranges(1,r+1))) // ','
  write(u,'(A)') '        "i_end": '   // trim(itoa(all_ranges(2,r+1))) // ','
  write(u,'(A)') '        "j_start": ' // trim(itoa(all_ranges(3,r+1))) // ','
  write(u,'(A)') '        "j_end": '   // trim(itoa(all_ranges(4,r+1))) // ','
  write(u,'(A)') '        "k_start": ' // trim(itoa(all_ranges(5,r+1))) // ','
  write(u,'(A)') '        "k_end": '   // trim(itoa(all_ranges(6,r+1)))

  if (r < nprocs-1) then
    write(u,'(A)') '      },'
  else
    write(u,'(A)') '      }'
  end if
end do

write(u,'(A)') '    ]'
write(u,'(A)') '  }'
write(u,'(A)') '}'
    close(u)

    if (allocated(variable_names)) deallocate(variable_names)

  end subroutine write_meta_json

  !subroutine write_meta_json(sim)
  !  use mod_model_config, only : nse_config, gpe_config
  !  type(simulation_config), intent(in) :: sim
  !  character(len=*), intent(in), optional :: filename
  !  character(len=*), intent(in), optional :: primary_variables_note

  !  integer :: u, i, nvar, ios
  !  character(len=512) :: fname

  !  call ensure_directory(cfg%output_dir)
  !  if (present(filename)) then
  !    fname = filename
  !  else
  !    fname = trim(cfg%output_dir)//'/meta.json'
  !  end if

  !  nvar = 0
  !  if (present(variable_names)) nvar = size(variable_names)

  !  open(newunit=u, file=trim(fname), status='replace', action='write', iostat=ios)
  !  if (ios /= 0) error stop 'ERROR: cannot write meta.json.'

  !  write(u,'(A)') '{'
  !  write(u,'(A,A,A)') '  "equation": "', trim(cfg%equation), '",'
  !  write(u,'(A,A,A)') '  "case_name": "', trim(cfg%case_name), '",'
  !  write(u,'(A,I0,A,I0,A,I0,A)') '  "grid": [', cfg%nx, ', ', cfg%ny, ', ', cfg%nz, '],'
  !  write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "domain_length": [', cfg%lx, ', ', cfg%ly, ', ', cfg%lz, '],'
  !  write(u,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "origin": [', cfg%x_min, ', ', cfg%y_min, ', ', cfg%z_min, '],'
  !  write(u,'(A,A,A)') '  "precision": "', trim(cfg%precision_name), '",'
  !  write(u,'(A,A,A)') '  "format": "', trim(cfg%output_format), '",'
  !  if (present(primary_variables_note)) then
  !    write(u,'(A,A,A)') '  "primary_variables_note": "', trim(primary_variables_note), '",'
  !  end if
  !  write(unit,'(A)') '  "variables": ['
  !  
  !  select case(trim(sim%equation))
  !  
  !  case("NSE")
  !  
  !      write(u,'(A)') '    "rho",'
  !      write(u,'(A)') '    "rho_u",'
  !      write(u,'(A)') '    "rho_v",'
  !      write(u,'(A)') '    "rho_w",'
  !      write(u,'(A)') '    "rho_E"'
  !  
  !  case("GPE")
  !  
  !      write(u,'(A)') '    "psi_real",'
  !      write(u,'(A)') '    "psi_imag"'
  !  
  !  end select
  !  
  !  write(u,'(A)') '  ],'    
  !  write(u,'(A)') '}'
  !  close(u)
  !end subroutine write_meta_json

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

  pure function itoa(i) result(str)
    integer, intent(in) :: i
    character(len=32) :: str
  
    write(str,'(I0)') i
    str = adjustl(str)
  end function itoa

  pure function json_bool(flag) result(str)
    logical, intent(in) :: flag
    character(len=5) :: str
  
    if (flag) then
      str = 'true '
    else
      str = 'false'
    end if
  end function json_bool

end module mod_slf_output
