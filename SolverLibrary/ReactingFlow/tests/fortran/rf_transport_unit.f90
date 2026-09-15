program rf_transport_unit
  use mod_rf_transport
  use mod_rf_flow1d
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  type(rf_transport) :: coeff,none
  type(rf_transport) :: separate
  real(dp) :: flux(4),yl(2),yr(2),y(2),cp,cv,h,e,r,t,p,rho,u,a,pi,expected,err32,err64
  real(dp) :: q(4,32),dq(4,32),boundary(4),dt,dt_none,sums(4),x
  integer :: i,kind
  m%ne=1; allocate(m%species(2),m%reactions(0))
  do i=1,2
    m%species(i)%mass=.01_dp; m%species(i)%pref=101325; m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp];m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp;m%species(i)%atoms=1
  end do
  m%species(1)%coeff(6,1)=1000 ! unequal formation enthalpies expose a missing species energy flux
  coeff=rf_transport(.75_dp,.25_dp,2._dp,.1_dp)
  yl=[.8_dp,.2_dp];yr=[.2_dp,.8_dp]
  call diffusive_flux(m,coeff,1._dp,2._dp,1000._dp,yl,1._dp,4._dp,1100._dp,yr,.5_dp,flux)
  call require(abs(flux(1)-.12_dp)<1.e-14_dp.and.abs(sum(flux(:2)))<1.e-14_dp,'Fick flux sign/mass closure')
  call require(abs(flux(3)+5)<1.e-14_dp,'Normal viscous stress including bulk viscosity')
  expected=-15-400+.12_dp*gas_r*1000/.01_dp
  call require(abs(flux(4)-expected)<1.e-8_dp,'Viscous work, Fourier heat, formation enthalpy transport')
  separate=coeff; separate%diffusivity=0; separate%species_diffusivity=[.1_dp,.3_dp]
  call diffusive_flux(m,separate,1._dp,2._dp,1000._dp,yl,1._dp,4._dp,1100._dp,yr,.5_dp,flux)
  call require(abs(flux(1)-.24_dp)<1.e-14_dp.and.abs(sum(flux(:2)))<1.e-14_dp, &
    'Unequal species diffusivity correction velocity')
  expected=-15-400+.24_dp*gas_r*1000/.01_dp
  call require(abs(flux(4)-expected)<1.e-8_dp,'Unequal diffusivity formation enthalpy flux')
  separate%species_diffusivity=.1_dp
  call diffusive_flux(m,separate,1._dp,2._dp,1000._dp,yl,1._dp,4._dp,1100._dp,yr,.5_dp,flux)
  call require(abs(flux(1)-.12_dp)<1.e-14_dp,'Equal species coefficients reproduce common D')
  call require(abs(diffusion_bound(separate)-.1_dp)<1.e-15_dp,'Equal species diffusion bound')
  m%species(1)%coeff(6,1)=0
  y=[.5_dp,.5_dp];pi=acos(-1._dp)
  call mixture(m,1000._dp,y,101325._dp,cp,cv,h,e,r)
  do i=1,32
    x=(i-.5_dp)/32
    call primitive_to_conserved(m,1000._dp,1000*r,sin(2*pi*x),y,q(:,i)) ! rho=1
  end do
  coeff=rf_transport(.75_dp,0._dp,0._dp,0._dp)
  call diffusion_rhs(m,q,1._dp/32,'periodic','periodic',coeff,dq,boundary)
  do i=1,32
    x=(i-.5_dp)/32
    expected=-(2*pi)**2*sin(2*pi*x)
    call require(abs(dq(3,i)-expected)<.13_dp,'Viscous momentum diffusion')
  end do
  call require(maxval(abs(sum(dq,dim=2)))<1.e-8_dp,'Periodic total energy conservation')
  call require(sum(dq(4,:)-q(3,:)*dq(3,:))>1,'Kinetic dissipation must heat internal energy')
  call diffusion_rhs(m,q,1._dp/32,'reflecting','reflecting',coeff,dq,boundary)
  call require(abs(boundary(4))<1.e-14_dp,'Stationary adiabatic wall has zero energy flux')
  call require(maxval(abs(sum(dq,dim=2)/32-boundary))<1.e-8_dp,'Wall traction conservation')

  ! A periodic heat/species Fourier mode is an analytic isolated diffusion test.
  ! It excludes Euler numerical diffusion, unlike a complete advecting flow run.
  do kind=1,2
    call sine_decay(32,kind,err32)
    call sine_decay(64,kind,err64)
    call require(err64<.3_dp*err32,'Second-order diffusion spatial convergence')
  end do
  dt_none=flow_timestep(m,q,.1_dp,.4_dp)
  dt=flow_timestep(m,q,.1_dp,.4_dp,none)
  call require(abs(dt-dt_none)<1.e-15_dp,'Disabled transport timestep compatibility')
  separate=rf_transport(); separate%species_diffusivity=[100._dp,200._dp]
  call require(transport_active(separate),'Species-only diffusion activates transport')
  dt=flow_timestep(m,q,.1_dp,.4_dp,separate)
  call require(dt<dt_none/2,'Species-only explicit timestep restriction')
  coeff=rf_transport(100._dp,100._dp,1.e6_dp,100._dp)
  dt=flow_timestep(m,q,.1_dp,.4_dp,coeff)
  call require(dt<dt_none/2,'Explicit diffusion timestep restriction')
  sums=sum(q,dim=2)*.1_dp
  call advance_flow(m,q,.1_dp,dt,.4_dp,'periodic','periodic',.false., &
                    1.e-9_dp,1.e-16_dp,1.e-8_dp,10000,boundary,coeff)
  call require(maxval(abs(sum(q,dim=2)*.1_dp-sums)/max(1._dp,abs(sums)))<1.e-12_dp,'Combined flow conservation')
  write(*,'(a)') '[OK] viscous/heat/species flux, enthalpy, boundaries, diffusion convergence, timestep'
  err32=unequal_rhs_error(32); err64=unequal_rhs_error(64)
  call require(err64<.3_dp*err32,'Unequal species diffusion operator second-order convergence')
  write(*,'(a,2es15.6)') '[OK] unequal species diffusion operator errors 32/64: ',err32,err64
contains
  real(dp) function unequal_rhs_error(nx) result(error)
    integer, intent(in) :: nx
    type(rf_transport) :: c
    real(dp) :: state(4,nx),rhs(4,nx),bnd(4),fractions(2),xx,exact,sn,cs
    integer :: j
    c%species_diffusivity=[.1_dp,.3_dp]
    do j=1,nx
      xx=(j-.5_dp)/nx
      fractions=[.5_dp+.1_dp*sin(2*pi*xx),.5_dp-.1_dp*sin(2*pi*xx)]
      call primitive_to_conserved(m,1000._dp,r*1000,0._dp,fractions,state(:,j))
    end do
    call diffusion_rhs(m,state,1._dp/nx,'periodic','periodic',c,rhs,bnd)
    error=0
    do j=1,nx
      xx=(j-.5_dp)/nx; sn=sin(2*pi*xx); cs=cos(2*pi*xx)
      ! J1=-(D1*(1-Y1)+D2*Y1)*grad(Y1), rho=1.
      exact=(2*pi)**2*(-(.2_dp+.02_dp*sn)*.1_dp*sn+.002_dp*cs*cs)
      error=error+abs(rhs(1,j)-exact)/nx
    end do
    call require(maxval(abs(sum(rhs,dim=2)))<1.e-9_dp,'Unequal diffusivity periodic conservation')
  end function
  subroutine sine_decay(nx,which,error)
    integer, intent(in) :: nx,which
    real(dp), intent(out) :: error
    type(rf_transport) :: c
    real(dp) :: state(4,nx),s1(4,nx),s2(4,nx),rhs(4,nx),bnd(4),step,time,spacing,value,exact
    real(dp) :: fractions(2),initial_sum(4)
    integer :: j
    c=rf_transport()
    if(which==1) c%conductivity=cv ! rho=1, alpha=k/(rho*cv)=1
    if(which==2) c%diffusivity=1
    spacing=1._dp/nx
    do j=1,nx
      x=(j-.5_dp)*spacing;fractions=[.5_dp,.5_dp];t=1000
      if(which==1) t=1000+20*sin(2*pi*x)
      if(which==2) fractions=[.5_dp+.1_dp*sin(2*pi*x),.5_dp-.1_dp*sin(2*pi*x)]
      call primitive_to_conserved(m,t,r*t,0._dp,fractions,state(:,j))
    end do
    initial_sum=sum(state,dim=2);time=0
    do while(time<.002_dp)
      step=min(.2_dp*spacing**2,.002_dp-time)
      call diffusion_rhs(m,state,spacing,'periodic','periodic',c,rhs,bnd)
      s1=state+step*rhs
      call diffusion_rhs(m,s1,spacing,'periodic','periodic',c,rhs,bnd)
      s2=.75_dp*state+.25_dp*(s1+step*rhs)
      call diffusion_rhs(m,s2,spacing,'periodic','periodic',c,rhs,bnd)
      state=state/3+2._dp/3*(s2+step*rhs)
      time=time+step
    end do
    call require(maxval(abs(sum(state,dim=2)-initial_sum)/max(1._dp,abs(initial_sum)))<1.e-12_dp, &
                 'Isolated diffusion global conservation')
    error=0
    do j=1,nx
      call conserved_to_primitive(m,state(:,j),rho,u,t,p,a,fractions)
      x=(j-.5_dp)*spacing
      if(which==1) then
        value=t;exact=1000+20*sin(2*pi*x)*exp(-(2*pi)**2*time)
      else
        value=fractions(1);exact=.5_dp+.1_dp*sin(2*pi*x)*exp(-(2*pi)**2*time)
      end if
      error=error+abs(value-exact)/nx
    end do
  end subroutine
end program
