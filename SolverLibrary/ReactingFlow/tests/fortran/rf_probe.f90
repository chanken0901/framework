program rf_probe
  use mod_rf_kinetics
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: path
  real(dp) :: t,p,cp,cv,h,e,r,s,rho,heat
  real(dp), allocatable :: y(:),qf(:),qr(:),net(:),omega(:)
  call get_command_argument(1,path)
  call read_mechanism(trim(path),m)
  allocate(y(size(m%species)),omega(size(m%species)))
  allocate(qf(size(m%reactions)),qr(size(m%reactions)),net(size(m%reactions)))
  read(*,*) t,p
  read(*,*) y
  call mixture(m,t,y,p,cp,cv,h,e,r,s)
  rho=p/(r*t)
  call rates(m,t,rho,y,qf,qr,net,omega,heat)
  write(*,'(*(es25.16e3,1x))') cp,cv,h,e,r,s,rho,temperature_from_energy(m,e,y),heat
  write(*,'(*(es25.16e3,1x))') qf
  write(*,'(*(es25.16e3,1x))') qr
  write(*,'(*(es25.16e3,1x))') net
  write(*,'(*(es25.16e3,1x))') omega
end program
