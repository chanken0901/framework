program rf_unit
  use mod_rf_reactor
  use mod_rf_units
  implicit none
  type(rf_mechanism) :: m
  type(rf_reference_scales) :: scales
  real(dp) :: cp,cv,h,e,r,s,heat,qf(1),qr(1),net(1),omega(2),y(2),row(6)
  integer :: i,u,ios
  character(2048) :: line
  logical :: success
  scales=rf_reference_scales(2._dp,10._dp,.5_dp,300._dp)
  call require(abs(scales%to_si(2._dp,'time')-.1_dp)<1.e-14_dp,'Time scale test')
  call require(abs(scales%from_si(-200._dp,'energy')+2)<1.e-14_dp,'Energy scale test')
  m%ne=1; m%source_hash='analytic'; m%canonical_hash='analytic'
  allocate(m%species(2),m%reactions(1))
  do i=1,2
    m%species(i)%mass=.01_dp; m%species(i)%pref=101325; m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp]; m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp; m%species(i)%atoms=1
  end do
  m%species(1)%name='A'; m%species(2)%name='B'
  associate(a=>m%reactions(1))
    a%kind=1; a%reversible=0; a%high=[log(1.e6_dp),0._dp,0._dp,1._dp]
    a%reactants=[1._dp,0._dp]; a%products=[0._dp,1._dp]; a%orders=a%reactants
  end associate
  y=[1._dp,0._dp]
  call mixture(m,1000._dp,y,101325._dp,cp,cv,h,e,r,s)
  call require(abs(cp-3.5_dp*gas_r/.01_dp)<1.e-10_dp,'NASA cp test')
  call require(abs(temperature_from_energy(m,e,y)-1000)<1.e-7_dp,'Energy inversion test')
  call rates(m,1000._dp,1._dp,y,qf,qr,net,omega,heat)
  call require(abs(omega(1)/1.e6_dp+1)<1.e-12_dp.and.abs(sum(omega))<1.e-8_dp,'Rate test')
  open(newunit=u,status='scratch',action='readwrite')
  call run_reactor(m,1000._dp,101325._dp,y,5.e-6_dp,.false.,1.e-10_dp,1.e-17_dp,1.e-9_dp,5.e-6_dp, &
                   100000,400._dp,u)
  rewind(u); success=.false.
  do
    read(u,'(a)',iostat=ios) line
    if(ios/=0) exit
    if(trim(line)=='# SUCCESS') success=.true.
    if(line(1:1)=='#'.or.line(1:4)=='time') cycle
    read(line,*) row
  end do
  close(u)
  call require(success,'Missing reactor success marker')
  call require(abs(row(5)-exp(-5._dp))<1.e-9_dp,'Stiff analytic decay test')
  call require(abs(row(2)-1000)<1.e-7_dp,'Isothermal equal-enthalpy reaction test')
  write(*,'(a)') '[OK] Fortran thermodynamics, units, kinetics, stiff reactor tests'
end program
