module mod_rf_thermo
  use mod_rf_mechanism
  implicit none
contains
  logical function temperature_supported(m,t) result(valid)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t
    integer :: i
    valid=.false.
    if(.not.ieee_is_finite(t)) return
    do i=1,size(m%species)
      if(t<m%species(i)%bounds(1).or.t>m%species(i)%bounds(size(m%species(i)%bounds))) return
    end do
    valid=.true.
  end function

  subroutine species_thermo(s,t,cp,h,entropy)
    type(rf_species), intent(in) :: s
    real(dp), intent(in) :: t
    real(dp), intent(out) :: cp,h,entropy
    real(dp) :: a(9)
    integer :: j
    call require(ieee_is_finite(t).and.t>=s%bounds(1).and.t<=s%bounds(size(s%bounds)), 'T outside NASA range')
    j=1
    do while(j<size(s%coeff,2))
      if(t<=s%bounds(j+1)) exit
      j=j+1
    end do
    a=s%coeff(:,j)
    if(s%model==7) then
      cp=a(1)+a(2)*t+a(3)*t**2+a(4)*t**3+a(5)*t**4
      h=a(1)+a(2)*t/2+a(3)*t**2/3+a(4)*t**3/4+a(5)*t**4/5+a(6)/t
      entropy=a(1)*log(t)+a(2)*t+a(3)*t**2/2+a(4)*t**3/3+a(5)*t**4/4+a(7)
    else
      cp=a(1)/t**2+a(2)/t+a(3)+a(4)*t+a(5)*t**2+a(6)*t**3+a(7)*t**4
      h=-a(1)/t**2+a(2)*log(t)/t+a(3)+a(4)*t/2+a(5)*t**2/3+a(6)*t**3/4+a(7)*t**4/5+a(8)/t
      entropy=-a(1)/(2*t**2)-a(2)/t+a(3)*log(t)+a(4)*t+a(5)*t**2/2+a(6)*t**3/3+a(7)*t**4/4+a(9)
    end if
    call require(cp>1.and.all(ieee_is_finite([cp,h,entropy])),'Invalid NASA properties')
    cp=cp*gas_r; h=h*gas_r*t; entropy=entropy*gas_r
  end subroutine

  subroutine check_y(m,y)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: y(:)
    call require(size(y)==size(m%species),'Invalid composition size')
    call require(all(ieee_is_finite(y)).and.all(y>=0),'Negative/nonfinite mass fractions; no clipping')
    call require(abs(sum(y)-1._dp)<=1.e-12_dp,'Mass fractions must sum to one')
  end subroutine

  subroutine mixture(m,t,y,p,cp,cv,h,e,rmix,entropy)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t,y(:),p
    real(dp), intent(out) :: cp,cv,h,e,rmix
    real(dp), optional, intent(out) :: entropy
    real(dp) :: cpi,hi,si,amount,total,smix
    integer :: i
    call check_y(m,y)
    call require(p>0.and.ieee_is_finite(p),'Invalid pressure')
    total=0
    do i=1,size(y)
      total=total+y(i)/m%species(i)%mass
    end do
    cp=0; h=0; smix=0
    do i=1,size(y)
      call species_thermo(m%species(i),t,cpi,hi,si)
      amount=y(i)/m%species(i)%mass
      cp=cp+amount*cpi; h=h+amount*hi
      if(amount>0) smix=smix+amount*(si-gas_r*log(amount/total*p/m%species(i)%pref))
    end do
    rmix=gas_r*total; cv=cp-rmix; e=h-rmix*t
    if(present(entropy)) entropy=smix
  end subroutine

  function temperature_from_energy(m,value,y,enthalpy,ok) result(t)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: value,y(:)
    logical, optional, intent(in) :: enthalpy
    logical, optional, intent(out) :: ok
    real(dp) :: t,lo,hi,cp,cv,h,e,r,res,low,high,roundoff,trial,derivative
    logical :: use_h
    integer :: i
    t=0
    if(present(ok)) ok=.false.
    use_h=.false.
    if(present(enthalpy)) use_h=enthalpy
    lo=0; hi=huge(hi)
    do i=1,size(m%species)
      lo=max(lo,m%species(i)%bounds(1))
      hi=min(hi,m%species(i)%bounds(size(m%species(i)%bounds)))
    end do
    call mixture(m,lo,y,101325._dp,cp,cv,h,e,r)
    low=merge(h,e,use_h)
    call mixture(m,hi,y,101325._dp,cp,cv,h,e,r)
    high=merge(h,e,use_h)
    ! Conservative -> internal energy subtracts kinetic energy. Permit only a
    ! machine-roundoff-sized endpoint discrepancy; never modify conserved state.
    roundoff=64*epsilon(value)*max(1._dp,abs(low),abs(high))
    if(ieee_is_finite(value)) then
      if(abs(value-low)<=roundoff) then
        t=lo
        if(present(ok)) ok=.true.
        return
      end if
      if(abs(value-high)<=roundoff) then
        t=hi
        if(present(ok)) ok=.true.
        return
      end if
    end if
    if(.not.(ieee_is_finite(value).and.value>=low.and.value<=high)) then
      if(present(ok)) return
      call require(.false.,'Energy outside NASA range')
    end if
    ! Safeguarded Newton: keep the NASA interval bracket. This avoids dozens
    ! of property evaluations per CFD state while retaining bisection fallback.
    t=lo+(hi-lo)*(value-low)/(high-low)
    do i=1,100
      call mixture(m,t,y,101325._dp,cp,cv,h,e,r)
      res=merge(h,e,use_h)-value
      ! The inversion must be more accurate than admissibility at NASA endpoints.
      ! The former 1e-11 relative residual injected energy errors into MUSCL's
      ! primitive -> conservative face conversion larger than that endpoint band.
      if(abs(res)<=8*epsilon(value)*max(1._dp,abs(value),abs(h),abs(e),abs(cp*t))) then
        if(present(ok)) ok=.true.
        return
      end if
      if(res>0) then
        hi=t
      else
        lo=t
      end if
      derivative=merge(cp,cv,use_h)
      trial=t-res/derivative
      if(.not.ieee_is_finite(trial).or.trial<=lo.or.trial>=hi.or.trial==t) trial=(lo+hi)/2
      t=trial
    end do
    t=0
    if(present(ok)) return
    call require(.false.,'Temperature inversion failed; polynomial discontinuity?')
  end function
end module
