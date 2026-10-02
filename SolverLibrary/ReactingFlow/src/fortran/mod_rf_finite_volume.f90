module mod_rf_finite_volume
  ! Dimension-independent, stationary-mesh inviscid finite-volume kernels.
  ! No mesh generator, time integrator, chemistry splitting or parallel runtime here.
  use mod_rf_thermo
  use mod_rf_flow1d, only: primitive_1d => primitive_to_conserved, &
    conserved_1d => conserved_to_primitive
  implicit none
  private
  public :: rf_face_mesh,validate_face_mesh,primitive_nd,conserved_nd
  public :: normal_flux,rusanov_normal_flux,reflect_normal,finite_volume_rhs

  type :: rf_face_mesh
    ! Area vector points out of owner, into neighbor. neighbor=0 is a boundary.
    integer, allocatable :: owner(:),neighbor(:)
    real(dp), allocatable :: area_vector(:,:),volume(:)
  end type
contains
  subroutine primitive_nd(m,t,p,velocity,y,q)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t,p,velocity(:),y(:)
    real(dp), intent(out) :: q(:)
    real(dp) :: line(size(m%species)+2),rho
    integer :: ns,nd
    ns=size(m%species);nd=size(velocity)
    call require(nd>=1.and.nd<=3.and.size(q)==ns+nd+1,'Invalid ND state layout')
    call require(all(ieee_is_finite(velocity)),'Nonfinite ND velocity')
    call primitive_1d(m,t,p,norm2(velocity),y,line)
    rho=sum(line(:ns))
    q(:ns)=line(:ns);q(ns+1:ns+nd)=rho*velocity;q(ns+nd+1)=line(ns+2)
  end subroutine

  subroutine conserved_nd(m,q,rho,velocity,t,p,sound,y,ok)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:)
    real(dp), intent(out) :: rho,velocity(:),t,p,sound,y(:)
    logical, optional, intent(out) :: ok
    real(dp) :: line(size(m%species)+2),speed
    logical :: valid
    integer :: ns,nd
    ns=size(m%species);nd=size(velocity)
    call require(nd>=1.and.nd<=3.and.size(q)==ns+nd+1.and.size(y)==ns,'Invalid ND state layout')
    rho=0;velocity=0;t=0;p=0;sound=0;y=0
    if(present(ok)) ok=.false.
    valid=all(ieee_is_finite(q)).and.all(q(:ns)>=0)
    if(valid) then
      rho=sum(q(:ns));valid=ieee_is_finite(rho).and.rho>0
    end if
    if(valid) then
      line(:ns)=q(:ns);line(ns+1)=norm2(q(ns+1:ns+nd));line(ns+2)=q(ns+nd+1)
      ! Reuse the same NASA inversion/admissibility policy as the 1D solver.
      call conserved_1d(m,line,rho,speed,t,p,sound,y,valid)
      if(valid) velocity=q(ns+1:ns+nd)/rho
    end if
    if(present(ok)) then
      ok=valid
    else
      call require(valid,'Invalid ND conservative state')
    end if
  end subroutine

  subroutine normal_flux(m,q,area_vector,flux,spectral_area)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:),area_vector(:)
    real(dp), intent(out) :: flux(:),spectral_area
    real(dp) :: rho,v(size(area_vector)),t,p,a,y(size(m%species)),vs,area
    integer :: ns,nd
    ns=size(m%species);nd=size(area_vector)
    call require(size(flux)==size(q),'Invalid ND flux shape')
    area=norm2(area_vector)
    call require(all(ieee_is_finite(area_vector)).and.ieee_is_finite(area).and.area>0, &
      'Invalid face area vector')
    call conserved_nd(m,q,rho,v,t,p,a,y)
    vs=dot_product(v,area_vector)
    flux=q*vs
    flux(ns+1:ns+nd)=flux(ns+1:ns+nd)+p*area_vector
    flux(ns+nd+1)=(q(ns+nd+1)+p)*vs
    spectral_area=abs(vs)+a*area
  end subroutine

  subroutine rusanov_normal_flux(m,left,right,area_vector,flux,spectral_area)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: left(:),right(:),area_vector(:)
    real(dp), intent(out) :: flux(:),spectral_area
    real(dp) :: fl(size(left)),fr(size(left)),sl,sr
    call require(size(left)==size(right).and.size(flux)==size(left),'Invalid face state shape')
    call normal_flux(m,left,area_vector,fl,sl)
    call normal_flux(m,right,area_vector,fr,sr)
    spectral_area=max(sl,sr)
    flux=(fl+fr-spectral_area*(right-left))/2
  end subroutine

  subroutine reflect_normal(m,q,area_vector,ghost)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:),area_vector(:)
    real(dp), intent(out) :: ghost(:)
    real(dp) :: n(size(area_vector)),area
    integer :: ns,nd
    ns=size(m%species);nd=size(area_vector)
    call require(nd>=1.and.nd<=3.and.size(q)==ns+nd+1.and.size(ghost)==size(q),'Invalid wall state shape')
    area=norm2(area_vector)
    call require(all(ieee_is_finite(area_vector)).and.ieee_is_finite(area).and.area>0,'Invalid wall normal')
    n=area_vector/area
    ghost=q
    ghost(ns+1:ns+nd)=q(ns+1:ns+nd)-2*dot_product(q(ns+1:ns+nd),n)*n
  end subroutine

  subroutine validate_face_mesh(mesh)
    type(rf_face_mesh), intent(in) :: mesh
    real(dp), allocatable :: closure(:,:),surface(:)
    real(dp) :: area
    integer :: nc,nf,nd,f,l,r
    call require(allocated(mesh%owner).and.allocated(mesh%neighbor).and. &
      allocated(mesh%area_vector).and.allocated(mesh%volume),'Unallocated face mesh')
    nc=size(mesh%volume);nf=size(mesh%owner);nd=size(mesh%area_vector,1)
    call require(nc>0.and.nf>0.and.nd>=1.and.nd<=3,'Empty or invalid face mesh')
    call require(size(mesh%neighbor)==nf.and.size(mesh%area_vector,2)==nf,'Face mesh shape mismatch')
    call require(all(ieee_is_finite(mesh%volume)).and.all(mesh%volume>0),'Nonpositive/nonfinite cell volume')
    call require(all(ieee_is_finite(mesh%area_vector)),'Nonfinite area vectors')
    allocate(closure(nd,nc),surface(nc));closure=0;surface=0
    do f=1,nf
      l=mesh%owner(f);r=mesh%neighbor(f)
      call require(l>=1.and.l<=nc.and.r>=0.and.r<=nc.and.l/=r,'Invalid face connectivity')
      area=norm2(mesh%area_vector(:,f))
      call require(ieee_is_finite(area).and.area>0,'Zero/overflow face area')
      closure(:,l)=closure(:,l)+mesh%area_vector(:,f);surface(l)=surface(l)+area
      if(r>0) then
        closure(:,r)=closure(:,r)-mesh%area_vector(:,f);surface(r)=surface(r)+area
      end if
    end do
    call require(all(ieee_is_finite(surface)).and.all(surface>0),'Disconnected cell/overflow surface')
    do l=1,nc
      call require(norm2(closure(:,l))/surface(l)<=1.e-12_dp,'Cell area vectors do not close')
    end do
  end subroutine

  subroutine finite_volume_rhs(m,mesh,q,ghost,cfl,dq,boundary_rate,dt)
    ! First-order inviscid residual. One shared flux per face ensures conservation.
    ! ghost(:,f) is required only for boundary faces; interior columns are ignored.
    type(rf_mechanism), intent(in) :: m
    type(rf_face_mesh), intent(in) :: mesh
    real(dp), intent(in) :: q(:,:),ghost(:,:),cfl
    real(dp), intent(out) :: dq(:,:),boundary_rate(:),dt
    real(dp) :: flux(size(q,1)),spectral_area,spectral_sum(size(q,2))
    integer :: ns,nd,nf,nc,f,l,r
    call validate_face_mesh(mesh)
    ns=size(m%species);nd=size(mesh%area_vector,1);nf=size(mesh%owner);nc=size(mesh%volume)
    call require(size(q,1)==ns+nd+1.and.size(q,2)==nc,'Invalid finite-volume field shape')
    call require(all(shape(dq)==shape(q)).and.size(boundary_rate)==size(q,1),'Invalid residual shape')
    call require(size(ghost,1)==size(q,1).and.size(ghost,2)==nf,'Invalid boundary ghost shape')
    call require(ieee_is_finite(cfl).and.cfl>0.and.cfl<=1,'Require CFL in (0,1]')
    dq=0;boundary_rate=0;spectral_sum=0
    do f=1,nf
      l=mesh%owner(f);r=mesh%neighbor(f)
      if(r>0) then
        call rusanov_normal_flux(m,q(:,l),q(:,r),mesh%area_vector(:,f),flux,spectral_area)
        dq(:,r)=dq(:,r)+flux
        spectral_sum(r)=spectral_sum(r)+spectral_area
      else
        call rusanov_normal_flux(m,q(:,l),ghost(:,f),mesh%area_vector(:,f),flux,spectral_area)
        boundary_rate=boundary_rate-flux
      end if
      dq(:,l)=dq(:,l)-flux
      spectral_sum(l)=spectral_sum(l)+spectral_area
    end do
    call require(all(ieee_is_finite(spectral_sum)).and.all(spectral_sum>0),'Invalid face spectral sum')
    do l=1,nc
      dq(:,l)=dq(:,l)/mesh%volume(l)
    end do
    dt=cfl*minval(mesh%volume/spectral_sum)
    call require(ieee_is_finite(dt).and.dt>0,'Invalid finite-volume timestep')
  end subroutine
end module
