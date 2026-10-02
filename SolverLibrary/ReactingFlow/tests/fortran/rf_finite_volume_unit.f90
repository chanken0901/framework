program rf_finite_volume_unit
  use mod_rf_finite_volume
  use mod_rf_thermo
  use mod_rf_flow1d, only: physical_flux
  implicit none
  type(rf_mechanism) :: m
  type(rf_face_mesh) :: mesh
  real(dp) :: q(5,2),ghost(5,7),dq(5,2),boundary(5),dt
  real(dp) :: y(2),rho,v(2),t,p,a,flux(5),reverse(5),speed,speed2
  real(dp) :: line(4),lineflux(4),wall(5),q3(6),f3(6),v3(3)
  integer :: i,f
  m%ne=1;allocate(m%species(2),m%reactions(0))
  do i=1,2
    m%species(i)%mass=.01_dp;m%species(i)%pref=101325;m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp];m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp;m%species(i)%atoms=1
  end do
  y=[.25_dp,.75_dp]
  call primitive_nd(m,1000._dp,101325._dp,[20._dp,0._dp],y,q(:,1))
  call conserved_nd(m,q(:,1),rho,v,t,p,a,y)
  call require(abs(t-1000)<1.e-7_dp.and.maxval(abs(v-[20._dp,0._dp]))<1.e-12_dp,'ND roundtrip')
  line=q([1,2,3,5],1)
  call physical_flux(m,line,lineflux,speed)
  call normal_flux(m,q(:,1),[1._dp,0._dp],flux,speed2)
  call require(maxval(abs(flux([1,2,3,5])-lineflux)/max(1._dp,abs(lineflux)))<1.e-12_dp, &
    'ND axis flux matches 1D')
  call require(abs(speed-speed2)<1.e-10_dp.and.flux(4)==0,'ND axis spectral speed')
  call primitive_nd(m,1100._dp,120000._dp,[-10._dp,30._dp],y,q(:,2))
  call rusanov_normal_flux(m,q(:,1),q(:,2),[3._dp,4._dp],flux,speed)
  call rusanov_normal_flux(m,q(:,2),q(:,1),[-3._dp,-4._dp],reverse,speed2)
  call require(all(flux==-reverse).and.speed==speed2,'Interior face orientation symmetry')
  call reflect_normal(m,q(:,2),[3._dp,4._dp],wall)
  call require(all(wall(:2)==q(:2,2)).and.wall(5)==q(5,2),'Wall preserves species and total energy')
  call require(abs(dot_product(wall(3:4)+q(3:4,2),[3._dp,4._dp]))<1.e-10_dp,'Wall normal reflection')
  call require(abs(dot_product(wall(3:4)-q(3:4,2),[-4._dp,3._dp]))<1.e-10_dp,'Wall tangential preservation')
  call primitive_nd(m,1000._dp,101325._dp,[20._dp,10._dp,5._dp],y,q3)
  call conserved_nd(m,q3,rho,v3,t,p,a,y)
  call require(maxval(abs(v3-[20._dp,10._dp,5._dp]))<1.e-12_dp,'3D kinetic energy roundtrip')
  call normal_flux(m,q3,[0._dp,0._dp,2._dp],f3,speed)
  call require(abs(f3(6)-(q3(6)+p)*10)<1.e-8_dp,'3D area-weighted energy flux')
  mesh%owner=[1,1,1,1,2,2,2];mesh%neighbor=[0,2,0,0,0,0,0]
  mesh%volume=[1._dp,1._dp]
  mesh%area_vector=reshape([-1._dp,0._dp,1._dp,0._dp,0._dp,-1._dp,0._dp,1._dp, &
    1._dp,0._dp,0._dp,-1._dp,0._dp,1._dp],[2,7])
  do f=1,7
    ghost(:,f)=q(:,mesh%owner(f))
  end do
  call finite_volume_rhs(m,mesh,q,ghost,.4_dp,dq,boundary,dt)
  call require(maxval(abs(sum(dq,dim=2)-boundary)/max(1._dp,abs(boundary)))<1.e-10_dp, &
    'Shared face flux cancels from volume integral')
  q(:,2)=q(:,1)
  do f=1,7
    ghost(:,f)=q(:,1)
  end do
  call finite_volume_rhs(m,mesh,q,ghost,.4_dp,dq,boundary,dt)
  call require(maxval(abs(dq))<1.e-8_dp,'Closed stationary geometry preserves uniform state')
  call require(dt>0,'Positive ND timestep')
  write(*,'(a)') '[OK] ND states, oriented flux, slip wall and finite-volume conservation'
end program
