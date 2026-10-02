program main_rf_mesh2d
  use mod_rf_grid2d
  use mod_rf_thermo, only: require
  implicit none
  type(rf_grid2d) :: grid
  character(2048) :: profile,output,arg
  integer :: ny,ios
  call require(command_argument_count()==3,'Usage: rf_mesh2d profile.rf ny new_output.vtk')
  call get_command_argument(1,profile);call get_command_argument(2,arg);call get_command_argument(3,output)
  ny=0;read(arg,*,iostat=ios) ny
  call require(ios==0.and.ny>0,'ny must be a positive integer')
  call read_nozzle2d(trim(profile),ny,grid)
  call write_grid2d_vtk(trim(output),grid)
  write(*,'(a,2(1x,i0))') '[OK] Grid nx ny:',grid%nx,grid%ny
  write(*,'(a)') '[OK] Mesh written: '//trim(output)
end program
