module mod_input_reader
  use mod_precision, only: dp
  use mod_common_parameters
  use mod_nse_parameters

  implicit none
  private

  public :: read_input

contains

  subroutine read_input(filename)
    character(len=*), intent(in) :: filename

    integer :: unit_no
    integer :: ios
    integer :: pos

    character(len=1024) :: line
    character(len=256)  :: key
    character(len=768)  :: value

    open(newunit=unit_no, file=trim(filename), status="old", &
         action="read", iostat=ios)

    if (ios /= 0) then
      write(*,'(A,A)') "[ERROR] Cannot open input file: ", trim(filename)
      stop
    end if

    do
      read(unit_no, '(A)', iostat=ios) line

      if (ios < 0) exit

      if (ios > 0) then
        write(*,*) "[ERROR] Failed while reading input.dat"
        close(unit_no)
        stop
      end if

      call remove_comment(line)

      if (len_trim(line) == 0) cycle

      pos = index(line, "=")
      if (pos == 0) cycle

      key   = adjustl(trim(line(:pos-1)))
      value = adjustl(trim(line(pos+1:)))

      call assign_input_value(trim(key), trim(value))
    end do

    close(unit_no)

    write(*,'(A,A)') "[OK] Read input file: ", trim(filename)

  end subroutine read_input


  subroutine remove_comment(line)
    character(len=*), intent(inout) :: line
    integer :: pos

    pos = index(line, "#")

    if (pos > 0) then
      if (pos == 1) then
        line = ""
      else
        line = line(:pos-1)
      end if
    end if

  end subroutine remove_comment


  subroutine assign_input_value(key, value)
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value

    integer :: ios

    ios = 0

    select case (trim(key))

    !==========================
    ! common: case
    !==========================
    case ("case_id")
      case_id = trim(value)

    case ("case_label")
      case_label = trim(value)

    case ("physics_model")
      physics_model = trim(value)

    !==========================
    ! common: grid
    !==========================
    case ("nx")
      read(value, *, iostat=ios) nx

    case ("ny")
      read(value, *, iostat=ios) ny

    case ("nz")
      read(value, *, iostat=ios) nz

    case ("x_min")
      read(value, *, iostat=ios) x_min

    case ("x_max")
      read(value, *, iostat=ios) x_max

    case ("y_min")
      read(value, *, iostat=ios) y_min

    case ("y_max")
      read(value, *, iostat=ios) y_max

    case ("z_min")
      read(value, *, iostat=ios) z_min

    case ("z_max")
      read(value, *, iostat=ios) z_max

    !==========================
    ! common: time/output
    !==========================
    case ("dt")
      read(value, *, iostat=ios) dt

    case ("t_max")
      read(value, *, iostat=ios) t_max

    case ("output_frequency")
      read(value, *, iostat=ios) output_frequency

    case ("output_location")
      output_location = trim(value)

    case ("restart_location")
      restart_location = trim(value)

    !==========================
    ! NSE parameters
    !==========================
    case ("gamma")
      read(value, *, iostat=ios) gamma

    case ("rho_0")
      read(value, *, iostat=ios) rho_0

    case ("mach_number")
      read(value, *, iostat=ios) mach_number

    case ("cfl")
      read(value, *, iostat=ios) cfl

    case ("scheme")
      scheme = trim(value)

    case ("reconstruction")
      reconstruction = trim(value)

    case ("time_integration")
      time_integration = trim(value)

    case ("flux")
      flux = trim(value)

    case default
      ! Python側が出力するが、今のNSEコードでは使わない項目は無視
      ! 例: project_name, author, status, description, raw_data_location
      ios = 0

    end select

    if (ios /= 0) then
      write(*,'(A,A,A,A)') "[ERROR] Invalid value: ", &
                           trim(key), " = ", trim(value)
      stop
    end if

  end subroutine assign_input_value

end module mod_input_reader