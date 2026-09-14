module mod_rf_kinetics
  use mod_rf_thermo
  implicit none
contains
  real(dp) function arr_log(a,t) result(v)
    real(dp), intent(in) :: a(4),t
    v=log_zero
    if(a(1)>log_zero/2) v=a(1)+a(2)*log(t)-a(3)/t
  end function
  real(dp) function exp_safe(x) result(v)
    real(dp), intent(in) :: x
    call require(x<=log(huge(x)).and.ieee_is_finite(x),'Rate overflow/nonfinite logarithm')
    v=0
    if(x>log(tiny(x))) v=exp(x)
  end function
  real(dp) function expm1_local(x) result(v)
    real(dp), intent(in) :: x
    if(abs(x)<1.e-5_dp) then
      v=x*(1+x*(0.5_dp+x*(1._dp/6+x*(1._dp/24+x/120))))
    else
      v=exp(x)-1
    end if
  end function
  real(dp) function group_log(g,t) result(v)
    type(rf_plog), intent(in) :: g
    real(dp), intent(in) :: t
    real(dp) :: terms(size(g%arr,2)),peak,total
    integer :: i
    do i=1,size(terms)
      terms(i)=arr_log(g%arr(:,i),t)
    end do
    peak=maxval(terms)
    call require(peak>log_zero/2,'Zero PLOG group')
    total=sum(g%arr(4,:)*exp(terms-peak))
    call require(total>0,'Nonpositive PLOG group')
    v=peak+log(total)
  end function
  function cheb_sequence(x,n) result(v)
    real(dp), intent(in) :: x
    integer, intent(in) :: n
    real(dp) :: v(n)
    integer :: i
    v(1)=1
    if(n>1) v(2)=x
    do i=3,n
      v(i)=2*x*v(i-1)-v(i-2)
    end do
  end function
  real(dp) function coefficient(r,t,p,c) result(v)
    type(rf_reaction), intent(in) :: r
    real(dp), intent(in) :: t,p,c(:)
    real(dp) :: f,collider,lh,ll,pr,blend,correction,fc,x,d,base,tt,pp,tx,px
    integer :: i,j,ng
    if(r%kind==1) then
      v=arr_log(r%high,t)
      return
    end if
    if(r%kind==6) then
      ng=size(r%groups)
      if(p<=r%groups(1)%pressure) then
        v=group_log(r%groups(1),t)
      else if(p>=r%groups(ng)%pressure) then
        v=group_log(r%groups(ng),t)
      else
        do i=1,ng-1
          if(p<r%groups(i+1)%pressure) exit
        end do
        f=log(p/r%groups(i)%pressure)/log(r%groups(i+1)%pressure/r%groups(i)%pressure)
        v=(1-f)*group_log(r%groups(i),t)+f*group_log(r%groups(i+1),t)
      end if
      return
    end if
    if(r%kind==7) then
      tt=t; pp=p
      do i=1,2
        if(abs(tt-r%bounds(i))<=2.e-14_dp*r%bounds(i)) tt=r%bounds(i)
        if(abs(pp-r%bounds(i+2))<=2.e-14_dp*r%bounds(i+2)) pp=r%bounds(i+2)
      end do
      call require(tt>=r%bounds(1).and.tt<=r%bounds(2).and.pp>=r%bounds(3).and.pp<=r%bounds(4), &
                   'Chebyshev outside fitted domain')
      tx=(2/tt-1/r%bounds(1)-1/r%bounds(2))/(1/r%bounds(2)-1/r%bounds(1))
      px=(2*log(pp)-log(r%bounds(3))-log(r%bounds(4)))/log(r%bounds(4)/r%bounds(3))
      block
        real(dp) :: ts(size(r%cheb,1)),ps(size(r%cheb,2))
        ts=cheb_sequence(tx,size(ts)); ps=cheb_sequence(px,size(ps)); v=0
        do j=1,size(ps)
          do i=1,size(ts)
            v=v+r%cheb(i,j)*ts(i)*ps(j)
          end do
        end do
        v=v*log(10._dp)
      end block
      return
    end if
    collider=sum(r%efficiencies*c); v=log_zero
    if(collider==0) return
    lh=arr_log(r%high,t)
    if(r%kind==2) then
      v=lh+log(collider)
      return
    end if
    ll=arr_log(r%low,t)
    if(min(ll,lh)<=log_zero/2) return
    pr=ll+log(collider)-lh
    blend=min(pr,0._dp)-log(1+exp(-abs(pr)))
    correction=0
    if(r%kind==4) then
      fc=(1-r%params(1))*exp(-t/r%params(2))+r%params(1)*exp(-t/r%params(3))
      if(size(r%params)==4) fc=fc+exp(-r%params(4)/t)
      call require(fc>0,'Invalid Troe Fcent')
      fc=log10(fc); x=pr/log(10._dp)-.4_dp-.67_dp*fc
      d=.75_dp-1.27_dp*fc-.14_dp*x
      correction=log(10._dp)*fc*d*d/(d*d+x*x)
    else if(r%kind==5) then
      base=r%params(1)*exp(-r%params(2)/t)+exp(-t/r%params(3))
      call require(base>0,'Invalid SRI base')
      correction=log(r%params(4))+log(base)/(1+(pr/log(10._dp))**2)+r%params(5)*log(t)
    end if
    v=lh+blend+correction
  end function
  real(dp) function mass_action(c,orders) result(v)
    real(dp), intent(in) :: c(:),orders(:)
    integer :: i
    v=0
    do i=1,size(c)
      if(orders(i)==0) cycle
      if(c(i)==0) then
        v=log_zero
        return
      end if
      v=v+orders(i)*log(c(i))
    end do
  end function
  subroutine rates(m,t,rho,y,qf,qr,net,omega,heat)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t,rho,y(:)
    real(dp), intent(out) :: qf(:),qr(:),net(:),omega(:),heat
    real(dp) :: c(size(y)),chemical(size(y)),hm(size(y)),cp,s,p,lk,lf,lr,kc
    integer :: i
    call check_y(m,y)
    call require(ieee_is_finite(rho).and.rho>0,'Invalid density')
    do i=1,size(y)
      c(i)=rho*y(i)/m%species(i)%mass
      call species_thermo(m%species(i),t,cp,hm(i),s)
      chemical(i)=-(hm(i)-t*s)/(gas_r*t)+log(m%species(i)%pref/(gas_r*t))
    end do
    p=gas_r*t*sum(c); omega=0
    do i=1,size(m%reactions)
      associate(r=>m%reactions(i))
        lk=coefficient(r,t,p,c); lf=lk+mass_action(c,r%orders); lr=log_zero
        if(r%reversible==1) then
          kc=sum((r%products-r%reactants)*chemical)
          lr=lk-kc+mass_action(c,r%products)
        end if
        qf(i)=exp_safe(lf); qr(i)=exp_safe(lr)
        if(lf==lr) then
          net(i)=0
        else if(lf>lr) then
          net(i)=-qf(i)*expm1_local(lr-lf)
        else
          net(i)=qr(i)*expm1_local(lf-lr)
        end if
        omega=omega+(r%products-r%reactants)*net(i)
      end associate
    end do
    heat=-sum(hm*omega)
    do i=1,size(y)
      omega(i)=omega(i)*m%species(i)%mass
    end do
    call require(all(ieee_is_finite(omega)).and.ieee_is_finite(heat),'Nonfinite production rates')
  end subroutine
end module
