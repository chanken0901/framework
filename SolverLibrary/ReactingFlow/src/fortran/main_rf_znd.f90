program main_rf_znd
  use mod_rf_thermo
  use mod_rf_znd
  use mod_rf_cj
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: mechanism_file,input_file,output_file
  real(dp) :: temperature=300,pressure=101325,mach=5,end_time=1.e-6_dp
  real(dp) :: rtol=1.e-9_dp,atol_species=1.e-16_dp,atol_velocity=1.e-7_dp,atol_distance=1.e-12_dp
  real(dp) :: max_step=1.e-8_dp,sonic_margin=1.e-6_dp
  real(dp), allocatable :: y(:)
  real(dp), allocatable :: yeq(:)
  real(dp) :: overdrive=1.05_dp,ratio_min=1.2_dp,ratio_max=2.5_dp,equilibrium_tolerance=0
  real(dp) :: cj_speed,t,p,rho,residual,cp,cv,h,e,r
  character(32) :: speed_mode='mach'
  integer :: unit,out,ios,max_steps=100000
  namelist /znd/ temperature,pressure,mach,end_time,rtol,atol_species,atol_velocity,atol_distance, &
    max_step,max_steps,sonic_margin,speed_mode,overdrive,ratio_min,ratio_max,equilibrium_tolerance
  call require(command_argument_count()==3,'Usage: rf_znd mechanism.rf znd.in output.csv')
  call get_command_argument(1,mechanism_file);call get_command_argument(2,input_file)
  call get_command_argument(3,output_file);call read_mechanism(trim(mechanism_file),m)
  allocate(y(size(m%species)))
  open(newunit=unit,file=trim(input_file),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open ZND input')
  read(unit,nml=znd,iostat=ios);call require(ios==0,'Invalid ZND namelist')
  read(unit,*,iostat=ios) y;call require(ios==0,'Missing ZND mass fractions');close(unit)
  select case(trim(speed_mode))
  case('mach')
    ! Legacy direct specification retained.
  case('cj')
    call require(ieee_is_finite(overdrive).and.overdrive>=1,'CJ speed multiplier must be >= 1')
    allocate(yeq(size(y)))
    call find_cj(m,temperature,pressure,y,ratio_min,ratio_max,cj_speed,t,p,rho,yeq,residual)
    call mixture(m,temperature,y,pressure,cp,cv,h,e,r)
    mach=overdrive*cj_speed/sqrt(cp/cv*r*temperature)
  case default
    call require(.false.,'speed_mode must be mach or cj')
  end select
  open(newunit=out,file=trim(output_file),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create ZND output (must not already exist)')
  if(trim(speed_mode)=='cj') then
    write(out,'(a,es25.16e3)') '# cj_speed=',cj_speed
    write(out,'(a,es25.16e3)') '# speed_multiplier_D_over_Dcj=',overdrive
  end if
  call run_znd(m,temperature,pressure,y,mach,end_time,rtol,atol_species,atol_velocity, &
    atol_distance,max_step,max_steps,sonic_margin,out,equilibrium_tolerance)
  close(out)
  if(trim(speed_mode)=='cj') then
    write(*,'(a)') '[OK] ZND completed using equilibrium CJ speed and specified speed multiplier'
  else
    write(*,'(a)') '[OK] Prescribed-speed ZND reached requested residence time (not CJ search)'
  end if
end program
