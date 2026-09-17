module mod_rf_transport
  use mod_rf_thermo
  implicit none
  private
  public :: rf_transport,validate_transport,diffusive_flux,transport_active
  public :: diffusion_bound,transport_scale
  public :: fixed_diffusive_flux
  public :: mixture_viscosity,viscosity_bound
  public :: mixture_conductivity,conductivity_bound
  type :: rf_transport
    ! Constant SI coefficients. The common species D is not a detailed mixture-averaged model.
    real(dp) :: viscosity=0,bulk_viscosity=0,conductivity=0,diffusivity=0
    real(dp), allocatable :: species_diffusivity(:)
    real(dp) :: reference_temperature=1,temperature_exponent=0
    real(dp), allocatable :: species_viscosity(:),sutherland_temperature(:)
    real(dp) :: viscosity_reference_temperature=300
    logical :: eucken_wms=.false.
  end type
contains
  function pure_viscosities(c,t) result(mu)
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: t
    real(dp) :: mu(size(c%species_viscosity))
    call require(ieee_is_finite(t).and.t>0,'Invalid viscosity temperature')
    mu=c%species_viscosity*(t/c%viscosity_reference_temperature)**1.5_dp * &
      (c%viscosity_reference_temperature+c%sutherland_temperature)/(t+c%sutherland_temperature)
    call require(all(ieee_is_finite(mu)).and.all(mu>0),'Nonfinite/zero Sutherland viscosity')
  end function

  real(dp) function mixture_viscosity(m,c,t,y) result(mu)
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: t,y(:)
    real(dp) :: pure(size(y))
    mu=c%viscosity
    if(.not.allocated(c%species_viscosity)) return
    call validate_transport(c,size(m%species)); call check_y(m,y)
    pure=pure_viscosities(c,t)
    mu=wilke_property(m,pure,y,pure)
    call require(ieee_is_finite(mu).and.mu>0,'Invalid Wilke viscosity')
  end function

  real(dp) function wilke_property(m,pure,y,property) result(value)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: pure(:),y(:),property(:)
    real(dp) :: x(size(y)),mass(size(y)),phi,denom
    integer :: i,j
    mass=m%species%mass;x=y/mass;x=x/sum(x);value=0
    do i=1,size(y)
      if(x(i)==0) cycle
      denom=0
      do j=1,size(y)
        if(x(j)==0) cycle
        phi=(1+sqrt(pure(i)/pure(j))*(mass(j)/mass(i))**.25_dp)**2 / &
          sqrt(8*(1+mass(i)/mass(j)))
        denom=denom+x(j)*phi
      end do
      value=value+x(i)*property(i)/denom
    end do
  end function

  real(dp) function mixture_conductivity(m,c,t,y) result(kappa)
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: t,y(:)
    real(dp) :: pure(size(y)),ki(size(y)),cp,h,s
    integer :: i
    kappa=c%conductivity
    if(.not.c%eucken_wms) return
    call validate_transport(c,size(m%species));call check_y(m,y)
    pure=pure_viscosities(c,t)
    do i=1,size(y)
      call species_thermo(m%species(i),t,cp,h,s)
      ki(i)=pure(i)*(cp+1.25_dp*gas_r)/m%species(i)%mass
    end do
    kappa=wilke_property(m,pure,y,ki)
    call require(ieee_is_finite(kappa).and.kappa>0,'Invalid Eucken/WMS conductivity')
  end function

  real(dp) function conductivity_bound(m,c,tmin,tmax) result(kappa)
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: tmin,tmax
    real(dp) :: upper(size(m%species)),cp_upper,piece,lo,hi,mass(size(m%species))
    integer :: i,j,n,power
    kappa=c%conductivity
    if(.not.c%eucken_wms) return
    call validate_transport(c,size(m%species))
    call require(tmin>0.and.tmax>=tmin,'Invalid conductivity bound interval')
    mass=m%species%mass
    do i=1,size(m%species)
      call require(tmin>=m%species(i)%bounds(1).and. &
        tmax<=m%species(i)%bounds(size(m%species(i)%bounds)),'Conductivity interval outside NASA range')
      cp_upper=0
      do j=1,size(m%species(i)%coeff,2)
        lo=max(tmin,m%species(i)%bounds(j));hi=min(tmax,m%species(i)%bounds(j+1))
        if(lo>hi) cycle
        piece=0
        do n=1,merge(5,7,m%species(i)%model==7)
          power=n-1
          if(m%species(i)%model==9) power=n-3
          piece=piece+abs(m%species(i)%coeff(n,j))*max(lo**power,hi**power)
        end do
        cp_upper=max(cp_upper,piece*gas_r)
      end do
      upper(i)=(cp_upper+1.25_dp*gas_r)/mass(i)
    end do
    kappa=maxval(pure_viscosities(c,tmax)*upper*sqrt(8*(1+mass/minval(mass))))
    call require(ieee_is_finite(kappa).and.kappa>0,'Invalid conductivity upper bound')
  end function

  real(dp) function viscosity_bound(m,c,tmax) result(mu)
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: tmax
    real(dp) :: mass(size(m%species))
    mu=c%viscosity
    if(.not.allocated(c%species_viscosity)) return
    call validate_transport(c,size(m%species))
    mass=m%species%mass
    ! phi_ij >= 1/sqrt(8*(1+Mi/Mj)); sum(x)=1. Bound holds for every composition.
    ! S>=0 makes each pure viscosity monotone in T, so Tmax bounds all diffusion faces.
    mu=maxval(pure_viscosities(c,tmax)*sqrt(8*(1+mass/minval(mass))))
    call require(ieee_is_finite(mu).and.mu>0,'Invalid Wilke viscosity bound')
  end function

  real(dp) function transport_scale(c,t) result(scale)
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: t
    call require(ieee_is_finite(t).and.t>0,'Invalid transport temperature')
    scale=1
    if(c%temperature_exponent==0) return
    scale=(t/c%reference_temperature)**c%temperature_exponent
    call require(ieee_is_finite(scale).and.scale>0,'Transport temperature scaling overflow/underflow')
  end function

  subroutine fixed_diffusive_flux(m,c,rho,u,t,y,ui,ti,yi,dx,side,flux)
    ! Prescribed physical boundary state, cell-center distance dx/2.
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: rho,u,t,y(:),ui,ti,yi(:),dx
    integer, intent(in) :: side ! -1: left; +1: right
    real(dp), intent(out) :: flux(:)
    real(dp) :: jmass(size(y)),tau,cp,h,s,gradt,gradu
    integer :: i,ns
    ns=size(m%species)
    call validate_transport(c,ns)
    call check_y(m,y); call check_y(m,yi)
    call require(all(ieee_is_finite([rho,u,t,ui,ti,dx])),'Nonfinite fixed boundary state')
    call require(abs(side)==1.and.dx>0.and.rho>0,'Invalid fixed boundary geometry/density')
    call require(size(flux)==ns+2,'Invalid fixed boundary flux size')
    gradu=side*(u-ui)/(dx/2); gradt=side*(t-ti)/(dx/2)
    jmass=-rho*c%diffusivity*side*(y-yi)/(dx/2)
    if(allocated(c%species_diffusivity)) jmass=-rho*c%species_diffusivity*side*(y-yi)/(dx/2)
    jmass=jmass-y*sum(jmass); jmass(ns)=-sum(jmass(:ns-1))
    tau=(4._dp/3*mixture_viscosity(m,c,t,y)+c%bulk_viscosity)*gradu
    flux(:ns)=jmass; flux(ns+1)=-tau; flux(ns+2)=-u*tau-mixture_conductivity(m,c,t,y)*gradt
    do i=1,ns
      call species_thermo(m%species(i),t,cp,h,s)
      flux(ns+2)=flux(ns+2)+h/m%species(i)%mass*jmass(i)
    end do
    flux=flux*transport_scale(c,t)
    call require(all(ieee_is_finite(flux)),'Nonfinite fixed boundary flux')
  end subroutine

  subroutine validate_transport(c,ns)
    type(rf_transport), intent(in) :: c
    integer, intent(in), optional :: ns
    real(dp) :: v(4)
    v=[c%viscosity,c%bulk_viscosity,c%conductivity,c%diffusivity]
    call require(all(ieee_is_finite(v)).and.all(v>=0),'Transport coefficients must be finite and nonnegative')
    call require(ieee_is_finite(c%reference_temperature).and.c%reference_temperature>0, &
      'Transport reference temperature must be positive and finite')
    call require(ieee_is_finite(c%temperature_exponent).and.c%temperature_exponent>=0, &
      'Transport exponent must be finite and nonnegative')
    call require(allocated(c%species_viscosity).eqv.allocated(c%sutherland_temperature), &
      'Sutherland viscosity requires both species arrays')
    if(c%eucken_wms) then
      call require(allocated(c%species_viscosity),'Eucken/WMS requires species viscosities')
      call require(c%conductivity==0,'Eucken/WMS cannot combine with constant conductivity')
    end if
    if(allocated(c%species_viscosity)) then
      call require(size(c%species_viscosity)>0,'Empty species viscosity array')
      if(present(ns)) call require(size(c%species_viscosity)==ns,'Species viscosity count mismatch')
      call require(size(c%sutherland_temperature)==size(c%species_viscosity),'Sutherland array size mismatch')
      call require(all(ieee_is_finite(c%species_viscosity)).and.all(c%species_viscosity>0), &
        'Species reference viscosities must be positive and finite')
      call require(all(ieee_is_finite(c%sutherland_temperature)).and.all(c%sutherland_temperature>=0), &
        'Sutherland temperatures must be nonnegative and finite')
      call require(ieee_is_finite(c%viscosity_reference_temperature).and.c%viscosity_reference_temperature>0, &
        'Invalid viscosity reference temperature')
      call require(c%viscosity==0.and.c%temperature_exponent==0, &
        'Wilke cannot combine with scalar viscosity or common temperature scaling')
    end if
    if(allocated(c%species_diffusivity)) then
      call require(size(c%species_diffusivity)>0,'Empty species diffusivity array')
      if(present(ns)) call require(size(c%species_diffusivity)==ns,'Species diffusivity count mismatch')
      call require(all(ieee_is_finite(c%species_diffusivity)).and.all(c%species_diffusivity>=0), &
        'Species diffusivities must be finite and nonnegative')
      call require(c%diffusivity==0,'Common and species diffusivities cannot be combined')
    end if
  end subroutine

  real(dp) function diffusion_bound(c) result(d)
    type(rf_transport), intent(in) :: c
    d=c%diffusivity
    ! Include a conservative estimate of the correction-velocity coupling.
    if(allocated(c%species_diffusivity)) d=2*maxval(c%species_diffusivity)-minval(c%species_diffusivity)
  end function

  logical function transport_active(c) result(active)
    type(rf_transport), intent(in) :: c
    active=max(c%viscosity,c%bulk_viscosity,c%conductivity,diffusion_bound(c))>0
    active=active.or.allocated(c%species_viscosity)
  end function

  subroutine diffusive_flux(m,c,rhol,ul,tl,yl,rhor,ur,tr,yr,dx,flux)
    ! Additional physical flux: dQ/dt = -d(F_Euler + F_diffusive)/dx.
    ! J=-rho*D*grad(Y), tau=(4mu/3+zeta)*grad(u), q_heat=-kappa*grad(T).
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: rhol,ul,tl,yl(:),rhor,ur,tr,yr(:),dx
    real(dp), intent(out) :: flux(:)
    real(dp) :: yf(size(yl)),jmass(size(yl)),tau,cp,h,s,temp,velocity
    integer :: i,ns
    ns=size(m%species)
    call validate_transport(c,ns)
    call require(size(flux)==ns+2,'Invalid diffusive flux size')
    call require(all(ieee_is_finite([rhol,ul,tl,rhor,ur,tr,dx])), 'Nonfinite transport face state')
    call require(min(rhol,rhor,dx)>0,'Invalid transport density/spacing')
    call check_y(m,yl); call check_y(m,yr)
    yf=(yl+yr)/2; temp=(tl+tr)/2; velocity=(ul+ur)/2
    jmass=-(rhol+rhor)/2*c%diffusivity*(yr-yl)/dx
    if(allocated(c%species_diffusivity)) jmass=-(rhol+rhor)/2*c%species_diffusivity*(yr-yl)/dx
    ! Correction velocity, followed by a roundoff-only closure on the final species.
    jmass=jmass-yf*sum(jmass)
    jmass(ns)=-sum(jmass(:ns-1))
    tau=(4._dp/3*mixture_viscosity(m,c,temp,yf)+c%bulk_viscosity)*(ur-ul)/dx
    flux(:ns)=jmass; flux(ns+1)=-tau
    flux(ns+2)=-velocity*tau-mixture_conductivity(m,c,temp,yf)*(tr-tl)/dx
    do i=1,ns
      call species_thermo(m%species(i),temp,cp,h,s)
      flux(ns+2)=flux(ns+2)+h/m%species(i)%mass*jmass(i)
    end do
    flux=flux*transport_scale(c,temp)
    call require(all(ieee_is_finite(flux)),'Nonfinite diffusive flux')
  end subroutine
end module
