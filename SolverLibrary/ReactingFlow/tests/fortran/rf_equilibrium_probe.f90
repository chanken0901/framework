program rf_equilibrium_probe
  use mod_rf_equilibrium
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: path
  real(dp) :: t,rho
  real(dp), allocatable :: y0(:),y(:)
  call get_command_argument(1,path)
  call read_mechanism(trim(path),m)
  allocate(y0(size(m%species)),y(size(m%species)))
  read(*,*) t,rho
  read(*,*) y0
  call equilibrium_tv(m,t,rho,y0,y)
  write(*,'(*(es25.16e3,1x))') y
end program
