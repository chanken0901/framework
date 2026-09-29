program rf_transport_probe
  use mod_rf_transport
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  type(rf_transport) :: c
  character(2048) :: path
  real(dp) :: t,p,cp,cv,h,e,r,rho
  real(dp), allocatable :: y(:),gradient(:),flux(:)
  integer :: i,j,n
  call get_command_argument(1,path);call read_mechanism(trim(path),m)
  call get_command_argument(2,path);call read_transport_table(trim(path),m,c)
  call get_command_argument(3,path);call read_binary_diffusion_data(trim(path),m,c,.true.)
  n=size(m%species);allocate(y(n),gradient(n),flux(n+2))
  read(*,*) t,p
  read(*,*) y
  read(*,*) gradient
  call mixture(m,t,y,p,cp,cv,h,e,r);rho=p/(r*t)
  write(*,'(*(es25.16e3,1x))') mixture_viscosity(m,c,t,y),mixture_conductivity(m,c,t,y), &
    viscosity_bound(m,c,t),conductivity_bound(m,c,t,t)
  do i=1,n
    write(*,'(*(es25.16e3,1x))') (binary_coefficient(c,i,j,t,p),j=1,n)
  end do
  call fixed_diffusive_flux(m,c,rho,0._dp,t,y,0._dp,t,y-gradient,2._dp,1,flux)
  write(*,'(*(es25.16e3,1x))') flux
end program
