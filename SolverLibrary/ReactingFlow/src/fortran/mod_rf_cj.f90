module mod_rf_cj
  use mod_rf_equilibrium
  use mod_rf_thermo
  implicit none
  private
  public :: find_cj
contains
  subroutine find_cj(m,t0,p0,y0,rmin,rmax,speed,t,p,rho,y,residual)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t0,p0,y0(:),rmin,rmax
    real(dp), intent(out) :: speed,t,p,rho,y(:),residual
    real(dp) :: cp,cv,h,e0,r,rho0,tlo,thi,left,right,c,d,fc,fd,ratio,f,fa,fb,es,hs
    real(dp) :: scan(21),ratios(21),s0,s1,aeq,a0
    integer :: i,kmin
    call require(all(ieee_is_finite([rmin,rmax])).and.rmin>1.and.rmax>rmin,'Invalid CJ density bounds')
    call mixture(m,t0,y0,p0,cp,cv,h,e0,r,s0)
    a0=sqrt(cp/cv*r*t0)
    rho0=p0/(r*t0);es=max(abs(e0),cp*t0,1._dp)
    tlo=0;thi=huge(thi)
    do i=1,size(y)
      tlo=max(tlo,m%species(i)%bounds(1));thi=min(thi,m%species(i)%bounds(size(m%species(i)%bounds)))
    end do
    left=rmin;right=rmax
    call wave(left,fa);call wave(right,fb)
    do i=1,21
      ratios(i)=rmin+(rmax-rmin)*real(i-1,dp)/20
      call wave(ratios(i),scan(i))
    end do
    kmin=minloc(scan,dim=1)
    call require(kmin>1.and.kmin<21,'CJ minimum not bracketed by supplied density interval')
    call require(all(scan(2:kmin)<scan(:kmin-1)).and.all(scan(kmin+1:)>scan(kmin:20)), &
      'CJ interval is not unimodal; narrow interval and inspect equilibrium branches')
    left=ratios(kmin-1);right=ratios(kmin+1)
    c=right-(right-left)*.6180339887498949_dp
    d=left+(right-left)*.6180339887498949_dp
    call wave(c,fc);call wave(d,fd)
    do i=1,80
      if(right-left<1.e-7_dp) exit
      if(fc<fd) then
        right=d;d=c;fd=fc;c=right-(right-left)*.6180339887498949_dp;call wave(c,fc)
      else
        left=c;c=d;fc=fd;d=left+(right-left)*.6180339887498949_dp;call wave(d,fd)
      end if
    end do
    ratio=(left+right)/2
    call wave(ratio,f);speed=sqrt(f)
    call require(ratio-rmin>1.e-4_dp.and.rmax-ratio>1.e-4_dp.and.f<min(fa,fb), &
      'CJ minimum not bracketed; no interior detonation solution')
    call mixture(m,t,y,p,cp,cv,h,hs,r,s1)
    residual=(h+.5_dp*(speed/ratio)**2-e0-p0/rho0-.5_dp*speed**2)/max(es,speed**2)
    call require(abs(residual)<1.e-8_dp,'CJ energy conservation failure')
    call require(s1>=s0.and.speed>a0,'CJ entropy/admissibility failure')
    aeq=equilibrium_sound(m,t,rho,y0)
    call require(abs(speed/ratio/aeq-1)<1.e-4_dp,'CJ equilibrium sonic condition failure')
  contains
    subroutine wave(compression,d2)
      real(dp), intent(in) :: compression
      real(dp), intent(out) :: d2
      real(dp) :: lo,hi,fl,fh,fmid
      integer :: k
      rho=rho0*compression;lo=tlo;hi=thi
      call hugoniot(lo,fl);call hugoniot(hi,fh)
      call require(fl<0.and.fh>0,'CJ Hugoniot root outside NASA range or unsupported branch')
      do k=1,80
        t=(lo+hi)/2;call hugoniot(t,fmid)
        if(abs(fmid)<1.e-11_dp*es) exit
        if(fmid>0) then
          hi=t
        else
          lo=t
        end if
      end do
      call require(p>p0,'CJ requires compressive high-pressure branch')
      d2=(p-p0)/(rho0*(1-1/compression))
    end subroutine
    subroutine hugoniot(temp,value)
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: value
      real(dp) :: cpv,cvv,hv,ev,rv
      call equilibrium_tv(m,temp,rho,y0,y)
      call mixture(m,temp,y,p0,cpv,cvv,hv,ev,rv)
      p=rho*rv*temp
      value=ev-e0-.5_dp*(p+p0)*(1/rho0-1/rho)
    end subroutine
  end subroutine
end module
