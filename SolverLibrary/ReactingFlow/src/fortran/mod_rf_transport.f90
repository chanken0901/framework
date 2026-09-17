module mod_rf_transport
  use mod_rf_thermo
  implicit none
  private
  public :: rf_transport,validate_transport,diffusive_flux,transport_active
  public :: diffusion_bound,transport_scale
  public :: fixed_diffusive_flux
  public :: mixture_viscosity,viscosity_bound
  public :: mixture_conductivity,conductivity_bound
  public :: read_transport_data
  type :: rf_transport
    ! Constant SI coefficients. The common species D is not a detailed mixture-averaged model.
    real(dp) :: viscosity=0,bulk_viscosity=0,conductivity=0,diffusivity=0
    real(dp), allocatable :: species_diffusivity(:)
    real(dp) :: reference_temperature=1,temperature_exponent=0
    real(dp), allocatable :: species_viscosity(:),sutherland_temperature(:)
    real(dp) :: viscosity_reference_temperature=300
    logical :: eucken_wms=.false.
    real(dp), allocatable :: binary_diffusivity(:,:)
  end type
contains
  subroutine read_transport_data(path,m,c)
    ! Versioned SI data, matched by exact species name rather than row order.
    character(*), intent(in) :: path
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(inout) :: c
    character(256) :: fields(4)
    integer :: unit,ios,n,ns,i,j,index
    logical :: found,seen(size(m%species))
    real(dp) :: mass
    ns=size(m%species)
    call require(.not.allocated(c%species_viscosity).and..not.allocated(c%sutherland_temperature), &
      'Transport file cannot overwrite existing species data')
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'Cannot open transport file: '//path)
    call transport_record(unit,fields(:1),found)
    call require(found,'Empty transport file')
    call require(trim(fields(1))=='RF_TRANSPORT_V1','Unsupported transport file version')
    call transport_record(unit,fields(:2),found)
    call require(found,'Missing transport count/reference temperature')
    call require(verify(trim(fields(1)),'0123456789')==0,'Invalid transport species count')
    read(fields(1),*,iostat=ios) n
    call require(ios==0,'Invalid transport species count')
    call require(n==ns,'Transport file must contain every mechanism species exactly once')
    c%viscosity_reference_temperature=transport_number(fields(2))
    allocate(c%species_viscosity(ns),c%sutherland_temperature(ns));seen=.false.
    do i=1,ns
      call transport_record(unit,fields,found)
      call require(found,'Missing transport species record')
      index=0
      do j=1,ns
        if(trim(fields(1))==trim(m%species(j)%name)) index=j
      end do
      call require(index>0,'Unknown transport species: '//trim(fields(1)))
      call require(.not.seen(index),'Duplicate transport species: '//trim(fields(1)))
      seen(index)=.true.;mass=transport_number(fields(2))
      call require(mass>0,'Transport molar mass must be positive (kg/mol)')
      call require(abs(mass-m%species(index)%mass)<=1.e-8_dp*m%species(index)%mass, &
        'Transport molar mass differs from mechanism (require kg/mol): '//trim(fields(1)))
      c%species_viscosity(index)=transport_number(fields(3))
      c%sutherland_temperature(index)=transport_number(fields(4))
    end do
    call transport_record(unit,fields,found)
    call require(.not.found,'Extra transport file records')
    close(unit)
    call validate_transport(c,ns)
  end subroutine

  subroutine transport_record(unit,fields,found)
    integer, intent(in) :: unit
    character(*), intent(out) :: fields(:)
    logical, intent(out) :: found
    character(4096) :: line
    integer :: ios,i,j,n,last
    found=.false.;fields=''
    do
      read(unit,'(a)',advance='no',iostat=ios) line
      if(is_iostat_end(ios)) return
      call require(is_iostat_eor(ios),'Transport line too long or read failure (limit 4095 characters)')
      i=scan(line,'#!');if(i>0) line(i:)=''
      if(len_trim(line)==0) cycle
      exit
    end do
    last=len_trim(line);i=1;n=0
    do while(i<=last)
      if(line(i:i)==' '.or.line(i:i)==achar(9)) then
        i=i+1;cycle
      end if
      j=i
      do while(j<=last)
        if(line(j:j)==' '.or.line(j:j)==achar(9)) exit
        j=j+1
      end do
      n=n+1
      call require(n<=size(fields),'Too many fields in transport record')
      call require(j-i<=len(fields),'Transport token too long')
      fields(n)=line(i:j-1);i=j
    end do
    call require(n==size(fields),'Wrong number of transport record fields')
    found=.true.
  end subroutine

  real(dp) function transport_number(token) result(value)
    character(*), intent(in) :: token
    integer :: ios
    call require(len_trim(token)>0.and.verify(trim(token),'0123456789+-.eEdD')==0, &
      'Invalid transport numeric token: '//trim(token))
    read(token,*,iostat=ios) value
    call require(ios==0,'Invalid transport number: '//trim(token))
    call require(ieee_is_finite(value),'Nonfinite transport value')
  end function

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
    if(allocated(c%binary_diffusivity)) jmass=binary_mass_flux(m,c,rho,y,side*(y-yi)/(dx/2))
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
    integer :: i,j,n
    v=[c%viscosity,c%bulk_viscosity,c%conductivity,c%diffusivity]
    if(allocated(c%binary_diffusivity)) then
      n=size(c%binary_diffusivity,1)
      call require(n>0.and.size(c%binary_diffusivity,2)==n,'Binary diffusion matrix must be square')
      if(present(ns)) call require(n==ns,'Binary diffusion species count mismatch')
      call require(c%diffusivity==0.and..not.allocated(c%species_diffusivity), &
        'Binary diffusion cannot combine with common/species diffusivity')
      call require(all(ieee_is_finite(c%binary_diffusivity)),'Nonfinite binary diffusivity')
      do i=1,n
        call require(c%binary_diffusivity(i,i)==0,'Binary diffusion diagonal must be zero')
        do j=i+1,n
          call require(c%binary_diffusivity(i,j)>0,'Binary diffusivities must be positive')
          call require(c%binary_diffusivity(i,j)==c%binary_diffusivity(j,i),'Binary diffusion must be symmetric')
        end do
      end do
    end if
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

  real(dp) function diffusion_bound(c,m) result(d)
    type(rf_transport), intent(in) :: c
    type(rf_mechanism), intent(in), optional :: m
    real(dp) :: ratio
    d=c%diffusivity
    ! Include a conservative estimate of the correction-velocity coupling.
    if(allocated(c%species_diffusivity)) d=2*maxval(c%species_diffusivity)-minval(c%species_diffusivity)
    if(allocated(c%binary_diffusivity)) then
      call require(present(m),'Binary diffusion bound requires mechanism')
      ratio=maxval(m%species%mass)/minval(m%species%mass)
      ! Bound the mole-gradient Jacobian and correction velocity in induced 1-norm.
      d=2*maxval(c%binary_diffusivity)*ratio*(1+ratio)
    end if
  end function

  function binary_mass_flux(m,c,rho,y,gradient) result(jmass)
    type(rf_mechanism), intent(in) :: m
    type(rf_transport), intent(in) :: c
    real(dp), intent(in) :: rho,y(:),gradient(:)
    real(dp) :: jmass(size(y)),x(size(y)),mass(size(y)),w,denom,numerator,d,weighted_gradient
    integer :: i,j
    mass=m%species%mass;w=1/sum(y/mass);x=y*w/mass
    weighted_gradient=w*sum(gradient/mass)
    do i=1,size(y)
      denom=0;numerator=0
      do j=1,size(y)
        if(j==i) cycle
        denom=denom+x(j)/c%binary_diffusivity(i,j)
        numerator=numerator+y(j)
      end do
      d=0
      if(denom>0) d=numerator/denom
      ! (M_i/W)*grad(X_i), evaluated by the chain rule at the face.
      ! At a pure-species face the undefined dominant raw flux is removed by correction.
      jmass(i)=-rho*d*(gradient(i)-y(i)*weighted_gradient)
    end do
  end function

  logical function transport_active(c) result(active)
    type(rf_transport), intent(in) :: c
    active=max(c%viscosity,c%bulk_viscosity,c%conductivity,c%diffusivity)>0
    if(allocated(c%species_diffusivity)) active=active.or.maxval(c%species_diffusivity)>0
    if(allocated(c%binary_diffusivity)) active=active.or.maxval(c%binary_diffusivity)>0
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
    if(allocated(c%binary_diffusivity)) jmass=binary_mass_flux(m,c,(rhol+rhor)/2,yf,(yr-yl)/dx)
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
