module mod_mc_reactive_config
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mod_precision, only : dp
  implicit none
  private

  integer, parameter :: reactive_path_length = 256

  type, public :: mc_reactive_config
    character(len=32) :: splitting_scheme = 'strang'
    character(len=32) :: chemistry_integrator = 'ssprk3_subcycled'
    real(dp) :: chemistry_cfl = 0.1_dp
    integer :: maximum_chemistry_substeps = 10000
    logical :: write_snapshots = .false.
    logical :: write_history = .false.
    integer :: output_every = 1
    character(len=reactive_path_length) :: snapshot_prefix = &
      'multicomponent_reactive'
    character(len=reactive_path_length) :: history_file = &
      'multicomponent_reactive_history.csv'
  end type mc_reactive_config

  public :: read_mc_reactive_config
  public :: validate_mc_reactive_config

contains

  subroutine read_mc_reactive_config(path,config)
    character(len=*), intent(in) :: path
    type(mc_reactive_config), intent(out) :: config
    integer :: unit, ios
    character(len=512) :: message
    character(len=32) :: splitting_scheme, chemistry_integrator
    real(dp) :: chemistry_cfl
    integer :: maximum_chemistry_substeps, output_every
    logical :: write_snapshots, write_history
    character(len=reactive_path_length) :: snapshot_prefix, history_file
    namelist /reactive_navier_stokes/ splitting_scheme, &
      chemistry_integrator, chemistry_cfl, maximum_chemistry_substeps, &
      write_snapshots, write_history, output_every, snapshot_prefix, &
      history_file

    config = mc_reactive_config()
    splitting_scheme = config%splitting_scheme
    chemistry_integrator = config%chemistry_integrator
    chemistry_cfl = config%chemistry_cfl
    maximum_chemistry_substeps = config%maximum_chemistry_substeps
    write_snapshots = config%write_snapshots
    write_history = config%write_history
    output_every = config%output_every
    snapshot_prefix = config%snapshot_prefix
    history_file = config%history_file

    open(newunit=unit,file=trim(path),status='old',action='read', &
      iostat=ios,iomsg=message)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open reactive input: ', trim(path)
      write(*,'(A,A)') 'ERROR: ', trim(message)
      error stop 'failed to open reactive Navier-Stokes input'
    end if
    read(unit,nml=reactive_navier_stokes,iostat=ios,iomsg=message)
    close(unit)
    if (ios /= 0) then
      write(*,'(A,A)') &
        'ERROR: invalid reactive_navier_stokes namelist: ', trim(message)
      error stop 'failed to read reactive Navier-Stokes input'
    end if

    config%splitting_scheme = adjustl(splitting_scheme)
    config%chemistry_integrator = adjustl(chemistry_integrator)
    config%chemistry_cfl = chemistry_cfl
    config%maximum_chemistry_substeps = maximum_chemistry_substeps
    config%write_snapshots = write_snapshots
    config%write_history = write_history
    config%output_every = output_every
    config%snapshot_prefix = adjustl(snapshot_prefix)
    config%history_file = adjustl(history_file)
    call validate_mc_reactive_config(config)
  end subroutine read_mc_reactive_config

  subroutine validate_mc_reactive_config(config)
    type(mc_reactive_config), intent(in) :: config

    if (trim(config%splitting_scheme) /= 'strang') then
      error stop 'reactive solver supports splitting_scheme=strang'
    end if
    if (trim(config%chemistry_integrator) /= 'ssprk3_subcycled') then
      error stop 'reactive solver supports chemistry_integrator=ssprk3_subcycled'
    end if
    if (.not. ieee_is_finite(config%chemistry_cfl) .or. &
        config%chemistry_cfl <= 0.0_dp .or. &
        config%chemistry_cfl > 1.0_dp) then
      error stop 'reactive chemistry_cfl must be in (0,1]'
    end if
    if (config%maximum_chemistry_substeps < 1) then
      error stop 'reactive maximum chemistry substeps must be positive'
    end if
    if (config%output_every < 1) then
      error stop 'reactive output interval must be positive'
    end if
    if (config%write_snapshots .and. &
        len_trim(config%snapshot_prefix) == 0) then
      error stop 'reactive snapshot prefix must not be empty'
    end if
    if (config%write_history .and. len_trim(config%history_file) == 0) then
      error stop 'reactive history file must not be empty'
    end if
  end subroutine validate_mc_reactive_config

end module mod_mc_reactive_config
