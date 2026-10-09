module mod_init_imported_turbulence
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : int32, int64
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use mod_nse_checkpoint, only: load_checkpoint
  implicit none
  private

  integer, parameter :: nconserved = 5

  type :: slf_header
    integer :: step = 0
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: nvar = 0
    integer :: nghost = 0
    logical :: mapped = .false.
    real(dp) :: time = 0.0_dp
    real(dp) :: bounds(6) = 0.0_dp
    character(len=32), allocatable :: names(:)
    integer(int64) :: data_position = 0_int64
  end type slf_header

  public :: initialize_imported_turbulence
  public :: initialize_restart
  public :: imported_turbulence_weight
  public :: inherit_nse_parameters

contains

  subroutine inherit_nse_parameters(sim,nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(inout) :: nse
    type(slf_header) :: header
    character(len=1024) :: path
    character(len=8) :: magic
    character(len=32) :: viscous
    real(dp) :: values(5)
    integer :: unit,ios
    integer(int64) :: offset
    path=sim%restart_file
    if(len_trim(path)==0) then
      select case(trim(sim%initial_condition))
      case('imported_turbulence','shock_turbulence_interaction','shock_tube_turbulence_interaction')
        path=nse%imported_turbulence_file
      case default
        return
      end select
    end if
    open(newunit=unit,file=trim(path),access='stream',form='unformatted', &
      status='old',action='read',convert='little_endian',iostat=ios)
    if(ios/=0) error stop 'Cannot open source SLF for parameter inheritance'
    call read_slf_header(unit,header)
    offset=header%data_position+8_int64*header%nx*header%ny*header%nz*header%nvar
    read(unit,pos=offset,iostat=ios) magic
    if(ios/=0) error stop 'Source SLF lacks NSE parameters; reconvert with --source-input ORIGINAL/input.dat'
    if(magic/='NSEPAR1'//achar(0)) error stop 'Invalid NSE parameter trailer'
    read(unit,iostat=ios) values,viscous
    close(unit)
    if(ios/=0) error stop 'Truncated NSE parameter trailer'
    if(.not.all(ieee_is_finite(values))) error stop 'Nonfinite source NSE parameters'
    if(values(1)<=1.or.minval(values(3:4))<=0.or.values(5)<0.or.values(2)<0) &
      error stop 'Invalid source NSE parameters'
    if(trim(viscous)/='none'.and.values(2)<=0) error stop 'Source viscosity requires positive Reynolds number'
    if(trim(viscous)/='none'.and.trim(viscous)/='central6'.and.trim(viscous)/='fv2') &
      error stop 'Unsupported source viscosity model'
    if((trim(viscous)=='none').neqv.(trim(nse%viscous_scheme)=='none')) &
      error stop 'Source/destination disagree on viscous versus inviscid physics'
    nse%gamma=values(1);nse%reynolds=values(2);nse%prandtl=values(3)
    nse%rho0=values(4);nse%mach=values(5)
    if(sim%rank==0) write(*,'(A,5ES23.15)') '# inherited gamma,Re,Pr,rho0,mach: ',values
  end subroutine inherit_nse_parameters

  subroutine initialize_restart(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(inout) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    type(slf_header) :: header
    integer :: unit, ios, map(5), j, k, v, i
    integer(int64) :: offset
    real(dp) :: dx, dy, dz
    logical :: exists
    if (sim%rank == 0) then
      inquire(file=trim(sim%output_dir)//'/meta.json',exist=exists)
      if (exists) error stop 'restart requires a new output directory (existing meta.json found)'
    end if
    if(sim%grid_mapping/='uniform') then
      call load_checkpoint(q,sim,nse,js,je,ks,ke)
      return
    end if
    if (nse%nv /= 5 .or. sim%grid_mapping /= 'uniform') &
      error stop 'restart currently requires single-component uniform-grid NSE'
    open(newunit=unit,file=trim(sim%restart_file),access='stream',form='unformatted', &
      status='old',action='read',convert='little_endian',iostat=ios)
    if (ios /= 0) error stop 'cannot open restart SLF'
    call read_slf_header(unit,header)
    call validate_slf_header(header,sim,nse,map,dx,dy,dz)
    if (header%nx /= sim%nx .or. header%nvar /= 5 .or. &
        .not.nearly_equal(header%bounds(1),sim%x_min) .or. &
        .not.nearly_equal(header%bounds(2),sim%x_max)) error stop 'restart grid mismatch'
    if (.not.ieee_is_finite(header%time)) error stop 'nonfinite restart time'
    if (header%time < 0 .or. header%step < 0) error stop 'invalid restart time/step'
    if (sim%nsteps <= header%step .or. sim%t_max <= header%time) &
      error stop 'restart requires nsteps and t_max greater than saved step and time'
    do v=1,5
      do k=ks,ke
        do j=js,je
          offset = ((int(map(v)-1,int64)*header%nz+k-1)*header%ny+j-1)*header%nx
          read(unit,pos=header%data_position+8_int64*offset,iostat=ios) q(1:sim%nx,j,k,v)
          if (ios /= 0) error stop 'truncated restart SLF'
        end do
      end do
    end do
    close(unit)
    do k=ks,ke
      do j=js,je
        do i=1,sim%nx
          call require_admissible(q(i,j,k,1:5),nse,'restart state')
        end do
      end do
    end do
    sim%t=header%time
    sim%step=header%step
    if (sim%rank == 0) write(*,'(A,I0,A,ES24.16)') 'Restart loaded: step=',sim%step,', time=',sim%t
  end subroutine initialize_restart

  subroutine initialize_imported_turbulence(q, sim, nse, js, je, ks, ke, &
      embedded_first_i, embedded_last_i)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    integer, intent(out), optional :: embedded_first_i, embedded_last_i
    real(dp), intent(inout) :: q(1-sim%nghost:, &
      js-sim%nghost:, ks-sim%nghost:, :)

    type(slf_header) :: header
    character(len=32) :: mode
    real(dp), allocatable :: source_state(:,:)
    real(dp), allocatable :: source_line(:)
    real(dp) :: background(nconserved), imported(nconserved)
    real(dp) :: blended(nconserved)
    real(dp) :: background_rho, background_p, x_start, weight
    real(dp) :: source_dx, source_dy, source_dz
    integer :: unit, ios, i, j, k, ivar, iq, source_i
    integer :: first_i, last_i, offset_cells, region_cells
    integer :: variable_map(nconserved)
    integer(int64) :: element_offset, byte_position, bytes_per_real

    if (nse%nv /= nconserved) then
      error stop 'imported turbulence requires five conserved variables'
    end if
    if (len_trim(nse%imported_turbulence_file) == 0) then
      error stop 'imported turbulence requires imported_turbulence_file'
    end if

    open(newunit=unit, file=trim(nse%imported_turbulence_file), &
      access='stream', form='unformatted', status='old', action='read', &
      convert='little_endian', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open imported turbulence SLF: ', &
        trim(nse%imported_turbulence_file)
      error stop 'failed to open imported turbulence SLF'
    end if

    call read_slf_header(unit, header)
    call validate_slf_header(header, sim, nse, variable_map, &
      source_dx, source_dy, source_dz)

    mode = lowercase(trim(adjustl(nse%imported_turbulence_mode)))
    x_start = nse%imported_turbulence_x_start
    if (x_start < -1.0e250_dp) x_start = sim%x_min
    if (.not. ieee_is_finite(x_start)) error stop 'imported turbulence x_start must be finite'
    if (mode == 'periodic_embed') then
      if (x_start < sim%x_min .or. x_start >= sim%x_max) &
        error stop 'periodic_embed x_start lies outside the target x domain'
    end if
    offset_cells = nint((x_start - sim%x_min) / sim%dx)
    if (.not. nearly_equal(x_start, &
        sim%x_min + real(offset_cells, dp) * sim%dx)) then
      error stop 'imported turbulence x_start must lie on a target cell boundary'
    end if

    select case (mode)
    case ('periodic_embed')
      if (.not. ieee_is_finite(nse%imported_turbulence_x_length)) &
        error stop 'periodic_embed x_length must be finite'
      if (nse%imported_turbulence_x_length <= 0.0_dp .or. &
          nse%imported_turbulence_x_length > sim%x_max-sim%x_min) &
        error stop 'periodic_embed x_length must be positive and fit in the x domain'
      region_cells = nint(nse%imported_turbulence_x_length/sim%dx)
      if (region_cells < 1 .or. .not. nearly_equal(nse%imported_turbulence_x_length, &
          real(region_cells,dp)*sim%dx)) &
        error stop 'periodic_embed x_length must be an integer multiple of dx'
      first_i = offset_cells+1
      last_i = offset_cells+region_cells
      if (first_i < 1 .or. last_i > sim%nx) &
        error stop 'periodic_embed interval lies outside the target x domain'
      if (nse%imported_turbulence_blend_cells < 0 .or. &
          nse%imported_turbulence_blend_cells > region_cells/2) &
        error stop 'invalid periodic_embed blend cell count'
    case ('embed')
      first_i = offset_cells + 1
      last_i = first_i + header%nx - 1
      if (first_i < 1 .or. last_i > sim%nx) then
        error stop 'embedded turbulence block lies outside the target x domain'
      end if
      if (nse%imported_turbulence_blend_cells < 0 .or. &
          2*nse%imported_turbulence_blend_cells > header%nx) then
        error stop 'invalid imported turbulence blend cell count'
      end if
    case ('tile')
      first_i = 1
      last_i = sim%nx
      if (mod(sim%nx, header%nx) /= 0) then
        error stop 'tile mode requires target nx to be a multiple of source nx'
      end if
      if (nse%imported_turbulence_blend_cells /= 0) then
        error stop 'tile mode requires imported_turbulence_blend_cells=0'
      end if
    case default
      write(*,'(A,A)') 'ERROR: unsupported imported turbulence mode: ', &
        trim(mode)
      error stop 'unsupported imported turbulence mode'
    end select
    if (present(embedded_first_i)) embedded_first_i = first_i
    if (present(embedded_last_i)) embedded_last_i = last_i

    background_rho = nse%imported_turbulence_background_rho
    if (background_rho <= 0.0_dp) background_rho = nse%rho0
    background_p = nse%imported_turbulence_background_p
    if (background_p <= 0.0_dp) background_p = 1.0_dp / nse%gamma
    call primitive_to_conserved(background_rho, &
      nse%imported_turbulence_background_u, &
      nse%imported_turbulence_background_v, &
      nse%imported_turbulence_background_w, background_p, nse%gamma, &
      background)
    call require_admissible(background, nse, 'background state')

    do k = ks, ke
      do j = js, je
        do i = 1, sim%nx
          q(i,j,k,1:nconserved) = background
        end do
      end do
    end do

    allocate(source_state(header%nx,nconserved))
    allocate(source_line(header%nx))
    bytes_per_real = int(storage_size(0.0_dp) / 8, int64)

    do k = ks, ke
      do j = js, je
        do iq = 1, nconserved
          ivar = variable_map(iq)
          element_offset = int(ivar-1, int64) * &
            int(header%nx, int64) * int(header%ny, int64) * &
            int(header%nz, int64)
          element_offset = element_offset + int(k-1, int64) * &
            int(header%nx, int64) * int(header%ny, int64)
          element_offset = element_offset + int(j-1, int64) * &
            int(header%nx, int64)
          byte_position = header%data_position + &
            element_offset * bytes_per_real
          read(unit, pos=byte_position, iostat=ios) source_line
          if (ios /= 0) then
            error stop 'failed to read imported turbulence SLF data'
          end if
          source_state(:,iq) = source_line
        end do

        do i = first_i, last_i
          select case (mode)
          case ('periodic_embed')
            source_i = modulo(i-first_i,header%nx)+1
            weight = imported_turbulence_weight(i-first_i+1,region_cells, &
              nse%imported_turbulence_blend_cells)
          case ('embed')
            source_i = i - first_i + 1
            weight = imported_turbulence_weight(source_i, header%nx, &
              nse%imported_turbulence_blend_cells)
          case ('tile')
            source_i = modulo(i - 1 - offset_cells, header%nx) + 1
            weight = 1.0_dp
          end select

          imported = source_state(source_i,:)
          call require_admissible(imported, nse, 'source state')
          call add_velocity_offset(imported, nse)
          blended = (1.0_dp-weight) * background + &
            weight * imported
          call require_admissible(blended, nse, 'blended imported state')
          q(i,j,k,1:nconserved) = blended
        end do
      end do
    end do

    deallocate(source_state, source_line)
    close(unit)

    if (sim%rank == 0) then
      write(*,'(A,A,A,A,A,I0,A,I0,A,I0)') &
        'Imported turbulence initialized: file=', &
        trim(nse%imported_turbulence_file), ', mode=', trim(mode), &
        ', source_grid=', header%nx, 'x', header%ny, 'x', header%nz
    end if
  end subroutine initialize_imported_turbulence

  subroutine read_slf_header(unit, header)
    integer, intent(in) :: unit
    type(slf_header), intent(out) :: header

    character(len=8) :: magic
    integer(int32) :: version, dtype_code, ndim, shape4(4), meta(8)
    integer(int32) :: nvar_header
    integer :: ios, ivar

    read(unit, iostat=ios) magic
    if (ios /= 0 .or. magic(1:4) /= 'SLF1') then
      error stop 'imported turbulence file is not an SLF1 file'
    end if
    read(unit, iostat=ios) version
    if (ios /= 0 .or. version /= 1_int32) then
      error stop 'unsupported imported turbulence SLF version'
    end if
    read(unit, iostat=ios) dtype_code
    if (ios /= 0 .or. dtype_code /= 2_int32) then
      error stop 'imported turbulence SLF must contain float64 data'
    end if
    read(unit, iostat=ios) ndim
    if (ios /= 0 .or. ndim /= 4_int32) then
      error stop 'imported turbulence SLF must be a four-dimensional field'
    end if
    read(unit, iostat=ios) shape4
    if (ios /= 0) error stop 'failed to read imported turbulence SLF shape'
    read(unit, iostat=ios) meta
    if (ios /= 0) error stop 'failed to read imported turbulence SLF metadata'
    read(unit, iostat=ios) header%time
    if (ios /= 0) error stop 'failed to read imported turbulence SLF time'
    read(unit, iostat=ios) header%bounds
    if (ios /= 0) error stop 'failed to read imported turbulence SLF bounds'
    read(unit, iostat=ios) nvar_header
    if (ios /= 0 .or. nvar_header /= shape4(4)) then
      error stop 'invalid imported turbulence SLF variable count'
    end if

    header%nx = int(shape4(1))
    header%step = int(meta(1))
    header%ny = int(shape4(2))
    header%nz = int(shape4(3))
    header%nvar = int(shape4(4))
    header%nghost = int(meta(6))
    header%mapped = meta(8)==1
    allocate(header%names(header%nvar))
    do ivar = 1, header%nvar
      read(unit, iostat=ios) header%names(ivar)
      if (ios /= 0) then
        error stop 'failed to read imported turbulence SLF variable names'
      end if
    end do
    inquire(unit=unit, pos=header%data_position)
  end subroutine read_slf_header

  subroutine validate_slf_header(header, sim, nse, variable_map, &
      source_dx, source_dy, source_dz)
    type(slf_header), intent(in) :: header
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(out) :: variable_map(nconserved)
    real(dp), intent(out) :: source_dx, source_dy, source_dz

    character(len=32) :: name
    integer :: ivar

    if(header%mapped) error stop 'Mapped SLF requires nonuniform restart; uniform import is not supported'
    if (header%nx <= 0 .or. header%ny <= 0 .or. header%nz <= 0) then
      error stop 'imported turbulence SLF has an invalid grid shape'
    end if
    if (header%nvar < nconserved) then
      error stop 'imported turbulence SLF must contain five conserved fields'
    end if
    if (header%nghost /= 0) then
      error stop 'prepare a ghost-free SLF with nse_prepare_imported_turbulence.py'
    end if
    if (header%ny /= sim%ny .or. header%nz /= sim%nz) then
      error stop 'imported turbulence source and target ny/nz must match'
    end if

    source_dx = (header%bounds(2)-header%bounds(1)) / real(header%nx, dp)
    source_dy = (header%bounds(4)-header%bounds(3)) / real(header%ny, dp)
    source_dz = (header%bounds(6)-header%bounds(5)) / real(header%nz, dp)
    if (.not. nearly_equal(source_dx, sim%dx)) then
      error stop 'imported turbulence source and target dx must match'
    end if
    if (.not. nearly_equal(source_dy, sim%dy) .or. &
        .not. nearly_equal(source_dz, sim%dz)) then
      error stop 'imported turbulence source and target dy/dz must match'
    end if
    if (.not. nearly_equal(header%bounds(3), sim%y_min) .or. &
        .not. nearly_equal(header%bounds(4), sim%y_max) .or. &
        .not. nearly_equal(header%bounds(5), sim%z_min) .or. &
        .not. nearly_equal(header%bounds(6), sim%z_max)) then
      error stop 'imported turbulence source and target y/z bounds must match'
    end if
    if (nse%gamma <= 1.0_dp) then
      error stop 'imported turbulence requires gamma > 1'
    end if

    variable_map = 0
    do ivar = 1, header%nvar
      name = lowercase(trim(adjustl(header%names(ivar))))
      select case (name)
      case ('rho')
        variable_map(1) = ivar
      case ('rho_u')
        variable_map(2) = ivar
      case ('rho_v')
        variable_map(3) = ivar
      case ('rho_w')
        variable_map(4) = ivar
      case ('rho_e')
        variable_map(5) = ivar
      end select
    end do
    if (any(variable_map == 0)) then
      error stop 'imported turbulence SLF is missing a conserved variable'
    end if
  end subroutine validate_slf_header

  subroutine add_velocity_offset(state, nse)
    real(dp), intent(inout) :: state(nconserved)
    type(nse_config), intent(in) :: nse

    real(dp) :: rho, u, v, w, p

    rho = state(1)
    u = state(2) / rho
    v = state(3) / rho
    w = state(4) / rho
    p = (nse%gamma-1.0_dp) * (state(5) - &
      0.5_dp * rho * (u*u + v*v + w*w))
    u = u + nse%imported_turbulence_velocity_offset_x
    v = v + nse%imported_turbulence_velocity_offset_y
    w = w + nse%imported_turbulence_velocity_offset_z
    call primitive_to_conserved(rho, u, v, w, p, nse%gamma, state)
  end subroutine add_velocity_offset

  pure subroutine primitive_to_conserved(rho, u, v, w, p, gamma, state)
    real(dp), intent(in) :: rho, u, v, w, p, gamma
    real(dp), intent(out) :: state(nconserved)

    state(1) = rho
    state(2) = rho * u
    state(3) = rho * v
    state(4) = rho * w
    state(5) = p / (gamma-1.0_dp) + &
      0.5_dp * rho * (u*u + v*v + w*w)
  end subroutine primitive_to_conserved

  subroutine require_admissible(state, nse, label)
    real(dp), intent(in) :: state(nconserved)
    type(nse_config), intent(in) :: nse
    character(len=*), intent(in) :: label

    real(dp) :: pressure

    if (.not. all(ieee_is_finite(state))) then
      write(*,'(A,A)') 'ERROR: non-finite imported turbulence ', trim(label)
      error stop 'non-finite imported turbulence state'
    end if
    if (state(1) <= nse%small_rho) then
      write(*,'(A,A)') 'ERROR: non-positive imported turbulence ', trim(label)
      error stop 'imported turbulence density is below small_rho'
    end if
    pressure = (nse%gamma-1.0_dp) * (state(5) - &
      0.5_dp * sum(state(2:4)*state(2:4)) / state(1))
    if (.not. ieee_is_finite(pressure) .or. pressure <= nse%small_p) then
      write(*,'(A,A)') 'ERROR: non-positive imported turbulence ', trim(label)
      error stop 'imported turbulence pressure is below small_p'
    end if
  end subroutine require_admissible

  pure real(dp) function imported_turbulence_weight(source_i, source_nx, &
      blend_cells) result(weight)
    integer, intent(in) :: source_i, source_nx, blend_cells
    integer :: edge_distance
    real(dp) :: fraction, pi

    if (blend_cells <= 0) then
      weight = 1.0_dp
      return
    end if
    edge_distance = min(source_i, source_nx-source_i+1)
    if (edge_distance > blend_cells) then
      weight = 1.0_dp
      return
    end if
    fraction = real(edge_distance, dp) / real(blend_cells+1, dp)
    pi = acos(-1.0_dp)
    weight = 0.5_dp * (1.0_dp-cos(pi*fraction))
  end function imported_turbulence_weight

  pure logical function nearly_equal(left, right) result(equal)
    real(dp), intent(in) :: left, right
    real(dp) :: scale

    scale = max(1.0_dp, abs(left), abs(right))
    equal = abs(left-right) <= 1.0e-10_dp * scale
  end function nearly_equal

  pure function lowercase(value) result(lower)
    character(len=*), intent(in) :: value
    character(len=len(value)) :: lower
    integer :: code, index

    lower = value
    do index = 1, len(value)
      code = iachar(value(index:index))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(index:index) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase

end module mod_init_imported_turbulence
