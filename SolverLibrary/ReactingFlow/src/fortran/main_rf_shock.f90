program main_rf_shock
  use mod_rf_thermo
  use mod_rf_shock
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: mechanism_file,input_file,output_file
  real(dp) :: temperature=300,pressure=101325,mach=2
  real(dp) :: t1,p1,rho1,u1,speed,residual(3)
  real(dp), allocatable :: y(:)
  integer :: unit,out,ios
  namelist /shock/ temperature,pressure,mach
  call require(command_argument_count()==3,'Usage: rf_normal_shock mechanism.rf shock.in output.csv')
  call get_command_argument(1,mechanism_file)
  call get_command_argument(2,input_file)
  call get_command_argument(3,output_file)
  call read_mechanism(trim(mechanism_file),m)
  allocate(y(size(m%species)))
  open(newunit=unit,file=trim(input_file),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open shock input')
  read(unit,nml=shock,iostat=ios);call require(ios==0,'Invalid shock namelist')
  read(unit,*,iostat=ios) y;call require(ios==0,'Missing shock mass fractions')
  close(unit)
  call frozen_normal_shock(m,temperature,pressure,y,mach,t1,p1,rho1,u1,speed,residual)
  open(newunit=out,file=trim(output_file),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create shock output (must not already exist)')
  write(out,'(a)') '# model=frozen_normal_shock_not_CJ_or_ZND'
  write(out,'(a)') '# canonical_sha256='//m%canonical_hash
  write(out,'(a,*(es25.16e3,1x))') '# upstream_T_p_Mach=',temperature,pressure,mach
  write(out,'(a,*(es25.16e3,1x))') '# frozen_Y=',y
  write(out,'(a,*(es25.16e3,1x))') '# mass_momentum_energy_residual=',residual
  write(out,'(a)') 'shock_speed,downstream_temperature,downstream_pressure,downstream_density,downstream_lab_velocity'
  write(out,'(*(es25.16e3,:,","))') speed,t1,p1,rho1,u1
  write(out,'(a)') '# SUCCESS'
  close(out)
  write(*,'(a)') '[OK] Frozen normal shock reference completed (not CJ/ZND)'
end program
