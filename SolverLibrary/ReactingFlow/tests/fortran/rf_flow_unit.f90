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
  m%ne=1; allocate(m%species(2),m%reactions(0))
  do i=1,2
    m%species(i)%mass=.01_dp; m%species(i)%pref=101325; m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp];m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp;m%species(i)%atoms=1
  end do
  y=[.25_dp,.75_dp]
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
end program
