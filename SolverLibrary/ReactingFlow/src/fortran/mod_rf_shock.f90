module mod_rf_shock
  use mod_rf_thermo
  implicit none
  private
  public :: frozen_normal_shock
contains
  subroutine frozen_normal_shock(m,t0,p0,y,mach,t1,p1,rho1,u1,speed,residual)
    ! Stationary upstream, right-running normal shock; composition is frozen.
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t0,p0,y(:),mach
    real(dp), intent(out) :: t1,p1,rho1,u1,speed,residual(3)
    real(dp) :: cp,cv,h0,e,r,rho0,a0,lo,hi,mid,f,previous_f,ratio,limit,h1,s0,s1
    real(dp) :: tlo,thi,pressure,temp,scale,downstream_sound
    logical :: bracketed,valid
    integer :: i
    call require(ieee_is_finite(mach).and.mach>=1.001_dp,'Frozen shock requires upstream Mach >= 1.001')
    call mixture(m,t0,y,p0,cp,cv,h0,e,r,s0)
    rho0=p0/(r*t0);a0=sqrt(cp/cv*r*t0);speed=mach*a0
    tlo=0;thi=huge(thi)
    do i=1,size(m%species)
      tlo=max(tlo,m%species(i)%bounds(1));thi=min(thi,m%species(i)%bounds(size(m%species(i)%bounds)))
    end do
    ! At this compression T returns to T0; use the nontrivial compressive root.
    limit=speed**2/(r*t0);scale=max(abs(h0),speed**2,cp*t0,1._dp)
    lo=1+1.e-7_dp;call evaluate(lo,previous_f,valid)
    call require(valid.and.previous_f>0,'Cannot bracket weak frozen shock within thermodynamic range')
    bracketed=.false.
    do i=1,2000
      hi=exp(log(limit)*real(i,dp)/2000)
      if(hi<=lo) cycle
      call evaluate(hi,f,valid)
      ! Never extrapolate NASA data or bridge an out-of-range interval.
      call require(valid,'Frozen shock search exceeds NASA range; extend validated thermodynamic data')
      if(f<=0) then
        bracketed=.true.;exit
      end if
      lo=hi
    end do
    call require(bracketed,'No compressive frozen shock root found')
    do i=1,100
      mid=(lo+hi)/2;call evaluate(mid,f,valid)
      call require(valid,'Frozen shock root leaves NASA range')
      if(abs(f)<=1.e-12_dp*scale) exit
      if(f>0) then
        lo=mid
      else
        hi=mid
      end if
    end do
    ratio=mid
    rho1=rho0*ratio;p1=p0+rho0*speed**2*(1-1/ratio);t1=p1/(rho1*r)
    u1=speed*(1-1/ratio)
    call mixture(m,t1,y,p1,cp,cv,h1,e,r,s1)
    downstream_sound=sqrt(cp/cv*r*t1)
    residual(1)=(rho1*(speed-u1)-rho0*speed)/(rho0*speed)
    residual(2)=(p1+rho1*(speed-u1)**2-p0-rho0*speed**2)/(p0+rho0*speed**2)
    residual(3)=(h1+.5_dp*(speed-u1)**2-h0-.5_dp*speed**2)/scale
    call require(maxval(abs(residual))<1.e-9_dp,'Frozen shock conservation residual too large')
    call require(ratio>1.and.p1>p0.and.s1>=s0-1.e-8_dp*max(1._dp,abs(s0)), &
      'Frozen shock violates compression/entropy condition')
    call require(speed-u1<downstream_sound,'Frozen shock downstream is not subsonic in shock frame')
  contains
    subroutine evaluate(compression,value,ok)
      real(dp), intent(in) :: compression
      real(dp), intent(out) :: value
      logical, intent(out) :: ok
      real(dp) :: cpv,cvv,hv,ev,rv
      pressure=p0+rho0*speed**2*(1-1/compression)
      temp=pressure/(rho0*compression*r)
      ok=ieee_is_finite(temp).and.temp>=tlo.and.temp<=thi
      value=0
      if(.not.ok) return
      call mixture(m,temp,y,pressure,cpv,cvv,hv,ev,rv)
      value=hv+.5_dp*(speed/compression)**2-h0-.5_dp*speed**2
    end subroutine
  end subroutine
end module
