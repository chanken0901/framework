program rf_flow_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mod_rf_flow1d
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  real(dp) :: q(4,8),old(4,8),flux(4),expected(4),change(4),rho,u,t,p,a,y(2),dt,speed
  integer :: i
  integer :: rejected,k
  logical :: ok
  real(dp) :: result(4,8),ref(4,8),original_dt,reference_change(4),reference_dt
  character(16) :: method
  real(dp) :: fixed(4,2),fixed_copy(4,2),fixed_dt
  m%ne=1; allocate(m%species(2),m%reactions(0))
  do i=1,2
    m%species(i)%mass=.01_dp; m%species(i)%pref=101325; m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp];m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp;m%species(i)%atoms=1
  end do
  y=[.25_dp,.75_dp]
  block
    real(dp) :: inside(4),reservoir(4),ghost(4),rhoi,ui,ti,pi,ai,yi(2),sign_in,pressure_expected
    integer :: side
    do side=1,2
      sign_in=real(3-2*side,dp)
      call primitive_to_conserved(m,1100._dp,110000._dp,100*sign_in,[.25_dp,.75_dp],inside)
      call primitive_to_conserved(m,900._dp,101325._dp,105*sign_in,[.8_dp,.2_dp],reservoir)
      call conserved_to_primitive(m,inside,rhoi,ui,ti,pi,ai,yi)
      pressure_expected=pi+rhoi*ai*5
      call reacting_inlet(m,inside,reservoir,side,ghost)
      call conserved_to_primitive(m,ghost,rhoi,ui,ti,pi,ai,yi)
      call require(abs(pi/pressure_expected-1)<1.e-10_dp,'Inlet outgoing acoustic relation')
      call require(abs(ti-900)<1.e-7_dp.and.abs(ui-105*sign_in)<1.e-10_dp,'Inlet T and velocity')
      call require(maxval(abs(yi-[.8_dp,.2_dp]))<1.e-12_dp,'Inlet composition')
    end do
  end block
  do i=1,8
    call primitive_to_conserved(m,1000._dp,101325._dp,20._dp,y,q(:,i))
  end do
  call conserved_to_primitive(m,q(:,1),rho,u,t,p,a,y)
  call require(abs(t-1000)<1.e-7_dp.and.abs(u-20)<1.e-10_dp,'Primitive roundtrip')
  call physical_flux(m,q(:,1),flux,speed)
  expected=q(:,1)*20;expected(3)=expected(3)+101325;expected(4)=(q(4,1)+101325)*20
  call require(maxval(abs(flux-expected)/max(1._dp,abs(expected)))<1.e-10_dp,'Legacy Euler flux formula')
  call rusanov_flux(m,q(:,1),q(:,2),flux)
  call require(maxval(abs(flux-expected)/max(1._dp,abs(expected)))<1.e-10_dp,'Consistent numerical flux')
  fixed(:,1)=q(:,1); fixed(:,2)=q(:,8)
  call primitive_to_conserved(m,1000._dp,101325._dp,2000._dp,y,fixed(:,1))
  fixed_copy=fixed
  dt=flow_timestep(m,q,.1_dp,.4_dp)
  fixed_dt=flow_timestep(m,q,.1_dp,.4_dp,left_bc='dirichlet',right_bc='outflow',fixed_states=fixed)
  call require(fixed_dt<dt,'Reservoir wave speed restricts timestep')
  old=q
  call advance_flow(m,q,.1_dp,fixed_dt,.4_dp,'dirichlet','outflow',.false., &
    1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,change,reconstruction='muscl',fixed_states=fixed)
  call require(maxval(abs(.1_dp*sum(q-old,dim=2)-change))<1.e-8_dp,'Fixed inflow boundary conservation')
  call require(all(fixed==fixed_copy),'Reservoir is not evolved')
  q=old
  old=q;dt=flow_timestep(m,q,.1_dp,.4_dp)/2
  call advance_flow(m,q,.1_dp,dt,.4_dp,'periodic','periodic',.false.,1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,change)
  call require(maxval(abs(q-old)/max(1._dp,abs(old)))<1.e-12_dp,'Uniform periodic preservation')
  call require(maxval(abs(change))<1.e-12_dp,'Periodic boundary cancellation')
  do i=1,8
    call primitive_to_conserved(m,1000._dp,101325._dp,0._dp,y,q(:,i))
  end do
  old=q;dt=flow_timestep(m,q,.1_dp,.4_dp)/2
  call advance_flow(m,q,.1_dp,dt,.4_dp,'reflecting','reflecting',.true.,1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,change)
  call require(maxval(abs(q-old)/max(1._dp,abs(old)))<1.e-12_dp,'Inert chemistry and wall preservation')
  call require(admissible_flow(m,q),'Valid field accepted')
  old=q; q(1,4)=-1
  call require(.not.admissible_flow(m,q),'Negative species rejected without stop')
  q=old; q(4,4)=-1
  call require(.not.admissible_flow(m,q),'Energy outside NASA range rejected without stop')
  q=old; q(:,4)=0
  call require(.not.admissible_flow(m,q),'Vacuum rejected without stop')
  q=old; q(1,4)=ieee_value(0._dp,ieee_quiet_nan)
  call require(.not.admissible_flow(m,q),'NaN rejected without stop')
  do k=1,2
    method='first_order'
    if(k==2) method='muscl'
    do i=1,8
      call primitive_to_conserved(m,200.01_dp,101325._dp,10._dp*(i-4.5_dp),y,q(:,i))
    end do
    old=q
    dt=flow_timestep(m,q,.1_dp,.4_dp,reconstruction=method,left_bc='outflow',right_bc='outflow')
    original_dt=dt
    call transport_step(m,q,.1_dp,dt,.4_dp,'outflow','outflow',result,change,ok,reconstruction=method)
    call require(.not.ok,'Near NASA floor expansion rejects intermediate stage')
    call require(all(q==old).and.all(result==old).and.all(change==0),'Rejected step leaves no update or flux')
    call advance_flow(m,q,.1_dp,dt,.4_dp,'outflow','outflow',.false.,1.e-9_dp,1.e-16_dp,1.e-8_dp, &
      10000,change,reconstruction=method,rejected_steps=rejected)
    call require(rejected>0.and.dt<original_dt.and.admissible_flow(m,q),'Reduced step recovers valid state')
    call require(abs(dt-original_dt*2._dp**(-rejected))<1.e-20_dp,'Retry count matches accepted dt')
    ref=old; reference_dt=dt
    call advance_flow(m,ref,.1_dp,reference_dt,.4_dp,'outflow','outflow',.false.,1.e-9_dp,1.e-16_dp, &
      1.e-8_dp,10000,reference_change,reconstruction=method)
    call require(all(ref==q).and.all(change==reference_change),'Retry starts from original state')
    call require(maxval(abs(sum(q-old,dim=2)*.1_dp-change)/max(1._dp,abs(sum(old,dim=2)*.1_dp))) &
      <1.e-12_dp,'Rejected boundary fluxes excluded from conservation')
  end do
  write(*,'(a)') '[OK] invalid states rejected, reduced-step recovery and rollback for first-order/MUSCL'
  write(*,'(a)') '[OK] 1D flux, uniform state, periodic/wall conservation, inert chemistry'
  call check_characteristic()
contains
  subroutine check_characteristic()
    integer, parameter :: cells=64
    real(dp) :: base(4),refstates(4,2),field(4,cells),start(4),net(4),bc(4),ghost(4)
    real(dp) :: rho0,a0,step_dt,time,amp,temp,pressure,vel,rr,aa,yy(2),reflection(2),x,drho
    integer :: mode,j
    character(16) :: outlet
    yy=[.25_dp,.75_dp]
    call primitive_to_conserved(m,1000._dp,101325._dp,100._dp,yy,base)
    call conserved_to_primitive(m,base,rho0,vel,temp,pressure,a0,yy)
    refstates(:,1)=base;refstates(:,2)=base
    call characteristic_outlet(m,base,base,2,ghost)
    call require(maxval(abs(ghost-base)/max(1._dp,abs(base)))<1.e-10_dp,'Uniform characteristic outlet')
    call primitive_to_conserved(m,1000._dp,101325._dp,2*a0,yy,ghost)
    bc=ghost
    call characteristic_outlet(m,bc,base,2,ghost)
    call require(all(ghost==bc),'Supersonic outlet ignores reference')
    do mode=1,2
      outlet='dirichlet'
      if(mode==2) outlet='characteristic'
      do j=1,cells
        x=(real(j,dp)-.5_dp)/cells
        amp=10*exp(-((x-.3_dp)/.06_dp)**2)
        drho=amp/a0**2
        temp=(101325+amp)/(gas_r/.01_dp*(rho0+drho))
        call primitive_to_conserved(m,temp,101325+amp,100+amp/(rho0*a0),yy,field(:,j))
      end do
      start=sum(field,dim=2)/cells;net=0;time=0
      do while(time<.001_dp)
        step_dt=min(.001_dp-time,flow_timestep(m,field,1._dp/cells,.4_dp, &
          reconstruction='muscl',left_bc='dirichlet',right_bc=outlet,fixed_states=refstates))
        call advance_flow(m,field,1._dp/cells,step_dt,.4_dp,'dirichlet',outlet,.false., &
          1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,bc,reconstruction='muscl',fixed_states=refstates)
        net=net+bc;time=time+step_dt
      end do
      reflection(mode)=0
      do j=1,cells
        call conserved_to_primitive(m,field(:,j),rr,vel,temp,pressure,aa,yy)
        reflection(mode)=reflection(mode)+((pressure-101325)-rho0*a0*(vel-100))**2/cells
      end do
      reflection(mode)=sqrt(reflection(mode))
      call require(maxval(abs(sum(field,dim=2)/cells-start-net)/max(1._dp,abs(start)))<1.e-10_dp, &
        'Characteristic boundary flux conservation')
    end do
    write(*,'(a,2es16.8)') '[OK] reflected acoustic RMS: fixed / characteristic ',reflection
    ! A fixed exterior state with an upwind Riemann flux is also weakly reflecting here.
    ! Require small reflection; do not claim superiority over that existing boundary.
    call require(reflection(2)<1.e-3_dp,'Characteristic acoustic incoming RMS below 1e-4 of pulse amplitude')
  end subroutine
end program
