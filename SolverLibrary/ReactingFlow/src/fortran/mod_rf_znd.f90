module mod_rf_znd
  use mod_rf_reactor, only: dvode_t,unpack
  use mod_rf_kinetics
  use mod_rf_shock
  use mod_rf_equilibrium
  implicit none
  private
  public :: run_znd
  type, extends(dvode_t) :: znd_context
    type(rf_mechanism), pointer :: mechanism=>null()
    real(dp) :: flux,b,sonic_margin
    integer :: dependent
  end type
contains
  subroutine znd_state(me,state,y,t,p,rho,cp,cv,h,r,eta)
    type(znd_context), intent(in) :: me
    real(dp), intent(in) :: state(:),y(:)
    real(dp), intent(out) :: t,p,rho,cp,cv,h,r,eta
    real(dp) :: e,u
    integer :: i
    u=state(1)
    call require(ieee_is_finite(u).and.u>0.and.u<me%b,'Invalid ZND velocity/pressure')
    r=0
    do i=1,size(y)
      r=r+gas_r*y(i)/me%mechanism%species(i)%mass
    end do
    rho=me%flux/u;p=me%flux*(me%b-u);t=u*(me%b-u)/r
    call mixture(me%mechanism,t,y,p,cp,cv,h,e,r)
    eta=1-u*u/(cp/cv*r*t)
    call require(eta>me%sonic_margin,'ZND sonic limit reached; not a completed CJ solution')
  end subroutine

  subroutine znd_rhs(me,neq,time,state,derivative)
    class(dvode_t), intent(inout) :: me
    integer :: neq
    real(dp) :: time,state(neq),derivative(neq)
    real(dp) :: y(neq-1),t,p,rho,cp,cv,h,r,eta,hdot,rdot,cpi,hi,si,heat
    integer :: i,j
    select type(me)
    type is(znd_context)
      associate(m=>me%mechanism)
        block
          real(dp) :: qf(size(m%reactions)),qr(size(m%reactions)),net(size(m%reactions)),omega(size(y))
          call require(all(ieee_is_finite(state)),'Nonfinite ZND Newton trial')
          ! unpack uses its first slot as a dummy; x is that slot here.
          call unpack(state(2:),me%dependent,y)
          ! RHS extension only, identical policy to the homogeneous reactor.
          y=max(y,0._dp);call require(sum(y)>0,'Invalid ZND trial composition');y=y/sum(y)
          call znd_state(me,state,y,t,p,rho,cp,cv,h,r,eta)
          call rates(m,t,rho,y,qf,qr,net,omega,heat)
          hdot=0;rdot=0;j=2
          do i=1,size(y)
            call species_thermo(m%species(i),t,cpi,hi,si)
            hdot=hdot+hi/m%species(i)%mass*omega(i)/rho
            rdot=rdot+gas_r/m%species(i)%mass*omega(i)/rho
            if(i==me%dependent) cycle
            j=j+1;derivative(j)=omega(i)/rho
          end do
          derivative(1)=state(1)*(rdot/r-hdot/(cp*t))/eta
          derivative(2)=state(1)
        end block
      end associate
    class default
      call require(.false.,'Invalid ZND context')
    end select
  end subroutine

  subroutine run_znd(m,t0,p0,y0,mach,tend,rtol,atoly,atolu,atolx,maxstep,maxsteps,sonic_margin,unit,equilibrium_tolerance)
    type(rf_mechanism), intent(in), target :: m
    real(dp), intent(in) :: t0,p0,y0(:),mach,tend,rtol,atoly,atolu,atolx,maxstep,sonic_margin
    integer, intent(in) :: maxsteps,unit
    real(dp), optional, intent(in) :: equilibrium_tolerance
    type(znd_context) :: solver
    real(dp) :: t,p,rho,u,speed,residual(3),cp,cv,h,e,r,eta,h0,scale,time,previous,err,elemerr
    real(dp) :: state(size(y0)+1),atol(size(y0)+1),y(size(y0)),elem0(m%ne),elem(m%ne)
    real(dp) :: rw(22+9*(size(y0)+1)+2*(size(y0)+1)**2)
    integer :: iw(31+size(y0)),n,i,j,istate,steps
    real(dp) :: eqtol,yeq(size(y0)),eqerr,heat,peak_heat,peak_time,peak_distance
    real(dp) :: qf(size(m%reactions)),qr(size(m%reactions)),net(size(m%reactions)),omega(size(y0))
    eqtol=0
    if(present(equilibrium_tolerance)) eqtol=equilibrium_tolerance
    call require(ieee_is_finite(eqtol).and.eqtol>=0.and.eqtol<1,'Invalid equilibrium tolerance')
    peak_heat=-huge(peak_heat);peak_time=0;peak_distance=0
    call require(all(ieee_is_finite([tend,rtol,atoly,atolu,atolx,maxstep,sonic_margin])), &
      'Nonfinite ZND controls')
    call require(min(tend,atoly,atolu,atolx,maxstep)>0.and.maxsteps>0,'Invalid ZND time/tolerance')
    call require(rtol>=1.e-12_dp.and.rtol<=1.e-2_dp,'Invalid ZND rtol')
    call require(sonic_margin>=1.e-8_dp.and.sonic_margin<1,'Invalid ZND sonic margin')
    call frozen_normal_shock(m,t0,p0,y0,mach,t,p,rho,u,speed,residual)
    call mixture(m,t0,y0,p0,cp,cv,h,e,r)
    h0=h+.5_dp*speed**2;scale=max(abs(h0),cp*t0,speed**2,1._dp)
    solver%mechanism=>m;solver%flux=p0/(r*t0)*speed
    solver%b=p0/solver%flux+speed;solver%sonic_margin=sonic_margin
    solver%dependent=maxloc(y0,dim=1)
    n=size(state);state(1)=speed-u;state(2)=0;j=2
    elem0=0
    do i=1,size(y0)
      elem0=elem0+y0(i)/m%species(i)%mass*m%species(i)%atoms
      if(i==solver%dependent) cycle
      j=j+1;state(j)=y0(i)
    end do
    atol=atoly;atol(1)=atolu;atol(2)=atolx
    rw=0;iw=0;rw(1)=tend;rw(6)=maxstep;iw(6)=maxsteps
    time=0;istate=1;steps=0;err=0;elemerr=0
    call solver%initialize(f=znd_rhs)
    write(unit,'(a)') '# model=prescribed_speed_planar_inviscid_ZND_not_CJ_search'
    write(unit,'(a)') '# canonical_sha256='//m%canonical_hash
    write(unit,'(a,es25.16e3)') '# shock_speed=',speed
    write(unit,'(a,*(es25.16e3,1x))') '# upstream_T_p_Mach=',t0,p0,mach
    write(unit,'(a,*(es25.16e3,1x))') '# rtol_atolY_atolU_atolX_maxstep_sonic=',rtol,atoly,atolu,atolx,maxstep,sonic_margin
    write(unit,'(a)',advance='no') 'residence_time,distance,temperature,pressure,density,shock_frame_velocity,sonic_eta'
    do i=1,size(y0)
      write(unit,'(a)',advance='no') ',Y_'//trim(m%species(i)%name)
    end do
    write(unit,*)
    do
      call unpack(state(2:),solver%dependent,y);call check_y(m,y)
      call znd_state(solver,state,y,t,p,rho,cp,cv,h,r,eta)
      err=max(err,abs(h+.5_dp*state(1)**2-h0)/scale)
      elem=0
      do i=1,size(y)
        elem=elem+y(i)/m%species(i)%mass*m%species(i)%atoms
      end do
      elemerr=max(elemerr,maxval(abs(elem-elem0)/max(1._dp,abs(elem0))))
      call require(err<=max(100*rtol,1.e-7_dp).and.elemerr<=1.e-8_dp,'ZND conservation failure')
      call rates(m,t,rho,y,qf,qr,net,omega,heat)
      if(heat>peak_heat) then
        peak_heat=heat;peak_time=time;peak_distance=state(2)
      end if
      write(unit,'(*(es25.16e3,:,","))') time,state(2),t,p,rho,state(1),eta,y
      if(time>=tend) exit
      call require(steps<maxsteps,'ZND max_steps reached; partial output is not success')
      previous=time
      call solver%solve(n,state,time,tend,2,[rtol],atol,5,istate,1,rw,size(rw),iw,size(iw),22)
      call require(istate>=0.and.ieee_is_finite(time).and.time>previous,'ZND DVODE failure')
      call require(all(ieee_is_finite(state)),'Nonfinite ZND accepted state')
      steps=steps+1
    end do
    write(unit,'(a,es25.16e3)') '# energy_relative_error=',err
    write(unit,'(a,es25.16e3)') '# element_relative_error=',elemerr
    write(unit,'(a,i0)') '# accepted_steps=',steps
    write(unit,'(a,es25.16e3)') '# sampled_peak_heat_release_W_m3=',peak_heat
    write(unit,'(a,es25.16e3)') '# sampled_peak_time_s=',peak_time
    write(unit,'(a,es25.16e3)') '# sampled_peak_distance_m=',peak_distance
    if(eqtol>0) then
      call equilibrium_tv(m,t,rho,y0,yeq)
      eqerr=maxval(abs(y-yeq))
      write(unit,'(a,es25.16e3)') '# equilibrium_max_abs_Y_error=',eqerr
      call require(eqerr<=eqtol,'ZND end state not equilibrated; extend residence time or inspect mechanism')
      write(unit,'(a)') '# SUCCESS (requested residence time and equilibrium composition tolerance reached)'
    else
      write(unit,'(a)') '# SUCCESS (requested residence time reached; equilibrium not asserted)'
    end if
  end subroutine
end module
