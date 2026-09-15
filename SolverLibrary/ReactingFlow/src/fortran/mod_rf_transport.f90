module mod_rf_transport
  use mod_rf_thermo
  implicit none
  private
  public :: rf_transport,validate_transport,diffusive_flux,transport_active
  public :: diffusion_bound
  type :: rf_transport
    ! Constant SI coefficients. The common species D is not a detailed mixture-averaged model.
    real(dp) :: viscosity=0,bulk_viscosity=0,conductivity=0,diffusivity=0
    real(dp), allocatable :: species_diffusivity(:)
  end type
contains
  subroutine validate_transport(c,ns)
    type(rf_transport), intent(in) :: c
    integer, intent(in), optional :: ns
    real(dp) :: v(4)
    v=[c%viscosity,c%bulk_viscosity,c%conductivity,c%diffusivity]
    call require(all(ieee_is_finite(v)).and.all(v>=0),'Transport coefficients must be finite and nonnegative')
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
    tau=(4._dp/3*c%viscosity+c%bulk_viscosity)*(ur-ul)/dx
    flux(:ns)=jmass; flux(ns+1)=-tau
    flux(ns+2)=-velocity*tau-c%conductivity*(tr-tl)/dx
    do i=1,ns
      call species_thermo(m%species(i),temp,cp,h,s)
      flux(ns+2)=flux(ns+2)+h/m%species(i)%mass*jmass(i)
    end do
    call require(all(ieee_is_finite(flux)),'Nonfinite diffusive flux')
  end subroutine
end module
