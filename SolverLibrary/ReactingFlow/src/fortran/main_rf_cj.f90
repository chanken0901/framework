program main_rf_cj
  use mod_rf_thermo
  use mod_rf_cj
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: mechanism_file,input_file,output_file
  real(dp) :: temperature=300,pressure=101325,ratio_min=1.2_dp,ratio_max=2.5_dp
  real(dp) :: speed,t,p,rho,residual
  real(dp), allocatable :: y(:),y0(:)
  integer :: unit,out,ios,i
  namelist /cj/ temperature,pressure,ratio_min,ratio_max
  call require(command_argument_count()==3,'Usage: rf_cj mechanism.rf cj.in output.csv')
  call get_command_argument(1,mechanism_file);call get_command_argument(2,input_file)
  call get_command_argument(3,output_file);call read_mechanism(trim(mechanism_file),m)
  allocate(y(size(m%species)),y0(size(m%species)))
  open(newunit=unit,file=trim(input_file),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open CJ input')
  read(unit,nml=cj,iostat=ios);call require(ios==0,'Invalid CJ namelist')
  read(unit,*,iostat=ios) y0;call require(ios==0,'Missing CJ composition');close(unit)
  call find_cj(m,temperature,pressure,y0,ratio_min,ratio_max,speed,t,p,rho,y,residual)
  open(newunit=out,file=trim(output_file),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create CJ output (must not already exist)')
  write(out,'(a)') '# model=ideal_gas_equilibrium_Hugoniot_minimum'
  write(out,'(a)') '# canonical_sha256='//m%canonical_hash
  write(out,'(a,*(es25.16e3,1x))') '# upstream_T_p=',temperature,pressure
  write(out,'(a,*(es25.16e3,1x))') '# upstream_Y=',y0
  write(out,'(a,es25.16e3)') '# energy_residual=',residual
  write(out,'(a)',advance='no') 'cj_speed,temperature,pressure,density'
  do i=1,size(y)
    write(out,'(a)',advance='no') ',Y_'//trim(m%species(i)%name)
  end do
  write(out,*)
  write(out,'(*(es25.16e3,:,","))') speed,t,p,rho,y
  write(out,'(a)') '# SUCCESS'
  close(out)
  write(*,'(a)') '[OK] Equilibrium CJ reference completed'
end program
