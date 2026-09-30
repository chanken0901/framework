program rf_unit
  use mod_rf_reactor
  use mod_rf_units
  use mod_rf_shock
  use mod_rf_flow1d, only: primitive_to_conserved,advance_flow
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
  block
    real(dp) :: cp0,cv0,h0,e0,r0,tcheck
    call mixture(m,200._dp,y,101325._dp,cp0,cv0,h0,e0,r0)
    tcheck=temperature_from_energy(m,e0,y,ok=success)
    call require(success.and.tcheck==200._dp,'Exact NASA endpoint roundtrip')
    tcheck=temperature_from_energy(m,nearest(e0,-1._dp),y,ok=success)
    call require(success.and.abs(tcheck-200)<1.e-10_dp,'NASA endpoint roundoff accepted')
    tcheck=temperature_from_energy(m,e0-1._dp,y,ok=success)
    call require(.not.success,'Physical NASA undershoot still rejected')
  end block
  block
    real(dp) :: mach,t1,p1,rho1,u1,speed,residual(3),compression,pratio,rho0
    real(dp), parameter :: machs(4)=[1.001_dp,1.1_dp,2._dp,5._dp]
    rho0=101325._dp/(gas_r/.01_dp*300)
    do i=1,size(machs)
      mach=machs(i)
      call frozen_normal_shock(m,300._dp,101325._dp,y,mach,t1,p1,rho1,u1,speed,residual)
      compression=2.4_dp*mach**2/(.4_dp*mach**2+2)
      pratio=1+2*1.4_dp/2.4_dp*(mach**2-1)
      call require(abs(t1/(300*pratio/compression)-1)<1.e-8_dp,'Shock analytic temperature')
      call require(abs(p1/(101325*pratio)-1)<1.e-8_dp,'Shock analytic pressure')
      call require(abs(rho1/(rho0*compression)-1)<1.e-8_dp,'Shock analytic density')
      call require(abs(u1/speed-(1-1/compression))<1.e-8_dp,'Shock analytic velocity')
      call require(maxval(abs(residual))<1.e-10_dp,'Shock conservation')
    end do
  end block
  call mixture(m,1000._dp,y,101325._dp,cp,cv,h,e,r,s)
  call require(abs(cp-3.5_dp*gas_r/.01_dp)<1.e-10_dp,'NASA cp test')
  call require(abs(temperature_from_energy(m,e,y)-1000)<1.e-7_dp,'Energy inversion test')
  s=temperature_from_energy(m,e,y,ok=success)
  call require(success.and.abs(s-1000)<1.e-7_dp,'Recoverable energy inversion')
  s=temperature_from_energy(m,-huge(e),y,ok=success)
  call require(.not.success.and.s==0,'Recoverable NASA range rejection')
  block
    type(rf_mechanism) :: gap
    real(dp) :: gap_e,gap_t
    gap=m
    do i=1,2
      deallocate(gap%species(i)%bounds,gap%species(i)%coeff)
      allocate(gap%species(i)%bounds(3),gap%species(i)%coeff(9,2))
      gap%species(i)%bounds=[200._dp,1000._dp,4000._dp]
      gap%species(i)%coeff=0; gap%species(i)%coeff(1,:)=3.5_dp
      gap%species(i)%coeff(6,2)=1
    end do
    gap_e=e+gas_r/.01_dp/2
    gap_t=temperature_from_energy(gap,gap_e,y,ok=success)
    call require(.not.success.and.gap_t==0,'NASA polynomial gap rejected without stop or clipping')
  end block
  call rates(m,1000._dp,1._dp,y,qf,qr,net,omega,heat)
  block
    real(dp) :: trial_t,trial_y(2),field(4,2),initial(4,2),reference(4,2),boundary(4),delta,reference_dt
    integer :: rejected
    trial_t=1000;trial_y=y
    call advance_chemistry(m,1._dp,trial_t,trial_y,5.e-6_dp,1.e-9_dp,1.e-16_dp,1.e-8_dp,1,success)
    call require(.not.success,'Chemistry internal step cap must reject')
    call require(trial_t==1000.and.all(trial_y==y),'Rejected cell chemistry must not mutate T/Y')
    call primitive_to_conserved(m,1000._dp,101325._dp,0._dp,y,field(:,1))
    field(:,2)=field(:,1);initial=field;delta=5.e-6_dp
    call advance_flow(m,field,1._dp,delta,.4_dp,'periodic','periodic',.true., &
      1.e-9_dp,1.e-16_dp,1.e-8_dp,30,boundary,rejected_steps=rejected)
    call require(rejected>0.and.delta<5.e-6_dp,'Chemistry failure must halve complete flow timestep')
    reference=initial;reference_dt=delta
    call advance_flow(m,reference,1._dp,reference_dt,.4_dp,'periodic','periodic',.true., &
      1.e-9_dp,1.e-16_dp,1.e-8_dp,100000,boundary)
    call require(maxval(abs(field-reference))<1.e-12_dp,'Retried step must equal fresh step at accepted dt')
    call require(maxval(abs(field(3:,:)-initial(3:,:))/max(1._dp,abs(initial(3:,:))))<1.e-14_dp, &
      'Chemistry retry preserves momentum and total energy to roundoff')
  end block
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
