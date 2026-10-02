program rf_flow2d_unit
  use mod_rf_flow2d
  use mod_rf_grid2d
  use mod_rf_finite_volume
  use mod_rf_initial2d
  use mod_rf_thermo
  implicit none
  type(rf_mechanism), target :: m
  type(rf_grid2d) :: grid
  type(rf_boundary2d) :: bc
  real(dp) :: q(5,4),old(5,4),dq(5,4),rate(5),change(5),dt,allowed,y(2)
  integer :: i,rejected,k
  logical :: ok
  m%ne=1;allocate(m%species(2),m%reactions(0))
  do i=1,2
    m%species(i)%mass=.01_dp;m%species(i)%pref=101325;m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp];m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp;m%species(i)%atoms=1
  end do
  call build_nozzle2d([0._dp,1._dp,2._dp],[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],2,grid)
  y=[.25_dp,.75_dp]
  block
    real(dp) :: bg(5),hot(5)
    call primitive_nd(m,1000._dp,101325._dp,[0._dp,0._dp],y,bg)
    call primitive_nd(m,1200._dp,202650._dp,[0._dp,0._dp],y,hot)
    call initialize_region2d(grid,bg,hot,'split',1,1._dp,[0._dp,0._dp,0._dp,0._dp],q)
    call require(all(q(:,1)==hot).and.all(q(:,3)==hot),'Split lower-x region')
    call require(all(q(:,2)==bg).and.all(q(:,4)==bg),'Split background region')
    call initialize_region2d(grid,bg,hot,'box',1,0._dp,[0._dp,1._dp,0._dp,.5_dp],q)
    call require(all(q(:,1)==hot).and.all(q(:,2)==bg).and.all(q(:,3)==bg).and.all(q(:,4)==bg), &
      'Box selects only cell centers inside half-open interval')
  end block
  do i=1,4
    call primitive_nd(m,1000._dp,101325._dp,[0._dp,0._dp],y,q(:,i))
  end do
  old=q
  do k=1,2
    q=old
    call flow2d_rhs(m,grid,bc,q,.3_dp,dq,rate,allowed,ok)
    call require(ok,'Valid initial 2D residual')
    dt=min(allowed/2,1.e-7_dp)
    call advance_flow2d(m,grid,bc,q,dt,.3_dp,k==2,1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,change,rejected)
    call require(maxval(abs(q-old)/max(1._dp,abs(old)))<1.e-12_dp,'Uniform inert 2D preservation')
    call require(rejected==0,'Uniform 2D state should not retry')
  end do
  call primitive_nd(m,1100._dp,120000._dp,[20._dp,10._dp],y,q(:,2))
  old=q
  call flow2d_rhs(m,grid,bc,q,.3_dp,dq,rate,allowed,ok)
  dt=allowed/4
  call advance_flow2d(m,grid,bc,q,dt,.3_dp,.false.,1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,change,rejected)
  rate=-change
  do i=1,4
    rate=rate+grid%mesh%volume(i)*(q(:,i)-old(:,i))
  end do
  call require(maxval(abs(rate)/max(1._dp,abs(change)))<1.e-9_dp,'SSPRK boundary flux accounting')
  write(*,'(a)') '[OK] 2D stepping, inert chemistry and boundary conservation'
end program
