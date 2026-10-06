program main_grid_preview
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_input_reader, only: read_common_input
  use mod_grid_axis, only: grid_axis,build_sinh_axis
  implicit none
  type(simulation_config) :: sim
  type(grid_axis) :: axis
  character(2048) :: input,output
  integer :: out,ios,d,i,n(3)
  real(dp) :: lower(3),upper(3),strength(3)
  if(command_argument_count()/=2) error stop 'Usage: nse_grid_preview input.dat new_grid.csv'
  call get_command_argument(1,input);call get_command_argument(2,output)
  call read_common_input(trim(input),sim,allow_grid_preview=.true.)
  n=[sim%nx,sim%ny,sim%nz]
  lower=[sim%x_min,sim%y_min,sim%z_min];upper=[sim%x_max,sim%y_max,sim%z_max]
  strength=sim%grid_stretch
  open(newunit=out,file=trim(output),status='new',action='write',iostat=ios)
  if(ios/=0) error stop 'Cannot create grid preview (existing file?)'
  write(out,'(a)') 'axis,cell,left_edge,right_edge,center,width'
  do d=1,3
    call build_sinh_axis(n(d),sim%nghost,lower(d),upper(d),strength(d),.false.,axis)
    do i=1,n(d)
      write(out,'(i0,",",i0,4(",",es24.16))') d,i,axis%edge(i-1),axis%edge(i),axis%center(i),axis%width(i)
    end do
  end do
  close(out,iostat=ios)
  if(ios/=0) error stop 'Failed to close grid preview'
  print *, '[OK] Geometry preview only; no NSE calculation performed.'
end program
