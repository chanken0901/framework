program main_rf_reactor
  use mod_rf_reactor
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: input,mechanism_file,output_file
  character(32) :: mode='constant_volume'
  real(dp) :: temperature=1000,pressure=101325,end_time=.001_dp,rtol=1.e-7_dp
  real(dp) :: atol_species=1.e-14_dp,atol_temperature=1.e-6_dp,max_step=.001_dp,ignition_rise=400
  integer :: u,out,ios,max_steps=100000
  real(dp), allocatable :: y(:)
  namelist /reactor/ mode,temperature,pressure,end_time,rtol,atol_species,atol_temperature,max_step,max_steps,ignition_rise
  call require(command_argument_count()==3,'Usage: rf_reactor mechanism.rf reactor.in output.csv')
  call get_command_argument(1,mechanism_file)
  call get_command_argument(2,input)
  call get_command_argument(3,output_file)
  call read_mechanism(trim(mechanism_file),m)
  allocate(y(size(m%species)))
  open(newunit=u,file=trim(input),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open reactor input')
  read(u,nml=reactor,iostat=ios)
  call require(ios==0,'Invalid reactor namelist')
  read(u,*,iostat=ios) y
  call require(ios==0,'Expected mass fractions in mechanism species order after namelist')
  close(u)
  call require(mode=='constant_volume'.or.mode=='constant_pressure','Invalid reactor mode')
  open(newunit=out,file=trim(output_file),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create output (existing files are never overwritten)')
  call run_reactor(m,temperature,pressure,y,end_time,mode=='constant_pressure',rtol,atol_species, &
                   atol_temperature,max_step,max_steps,ignition_rise,out)
  close(out)
  write(*,'(a)') '[OK] Fortran reactor completed: '//trim(output_file)
end program
