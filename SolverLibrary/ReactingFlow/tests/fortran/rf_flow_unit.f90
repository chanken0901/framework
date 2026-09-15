program rf_flow_unit
  use mod_rf_flow1d
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  real(dp) :: q(4,8),old(4,8),flux(4),expected(4),change(4),rho,u,t,p,a,y(2),dt,speed
  integer :: i
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
  write(*,'(a)') '[OK] 1D flux, uniform state, periodic/wall conservation, inert chemistry'
end program
