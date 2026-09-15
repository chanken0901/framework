module mod_rf_flow1d
  use mod_rf_reactor, only: advance_chemistry
  use mod_rf_thermo
  use mod_rf_transport
  implicit none
  private
  public :: primitive_to_conserved,conserved_to_primitive,physical_flux,rusanov_flux
  public :: flow_timestep,advance_flow,transport_step,diffusion_rhs
  public :: reconstruct_faces,validate_reconstruction
  public :: admissible_flow
contains
  subroutine validate_reconstruction(method)
    character(*), intent(in), optional :: method
    if(present(method)) call require(method=='first_order'.or.method=='muscl','Unknown flow reconstruction')
  end subroutine

  pure elemental real(dp) function mc_slope(dl,dr) result(s)
    real(dp), intent(in) :: dl,dr
    s=0
    if((dl>0.and.dr>0).or.(dl<0.and.dr<0)) s=sign(min(2*abs(dl),2*abs(dr),abs((dl+dr)/2)),dl)
  end function

  subroutine reconstruct_faces(m,q,left_bc,right_bc,qminus,qplus,method)
    ! Primitive MC-MUSCL: [T,p,u,Y]. Only face states change; cell averages remain conservative.
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:)
    character(*), intent(in) :: left_bc,right_bc
    character(*), intent(in), optional :: method
    real(dp), intent(out) :: qminus(:,:),qplus(:,:)
    real(dp) :: v(size(m%species)+3,0:size(q,2)+1),s(size(m%species)+3)
    real(dp) :: lo(size(m%species)),hi(size(m%species)),rho,a,theta,bound
    integer :: i,j,ns,nx
    call validate_reconstruction(method)
    ns=size(m%species); nx=size(q,2)
    call require(nx>=2.and.size(q,1)==ns+2,'Invalid reconstruction field')
    call require(all(shape(qminus)==shape(q)).and.all(shape(qplus)==shape(q)),'Invalid face output shape')
    call require((left_bc=='periodic').eqv.(right_bc=='periodic'),'Periodic boundary must be paired')
    qminus=q; qplus=q
    if(.not.present(method)) return
    if(method=='first_order') return
    do i=1,nx
      call conserved_to_primitive(m,q(:,i),rho,v(3,i),v(1,i),v(2,i),a,v(4:,i))
    end do
    if(left_bc=='periodic') then
      v(:,0)=v(:,nx); v(:,nx+1)=v(:,1)
    else
      v(:,0)=v(:,1); v(:,nx+1)=v(:,nx)
      call require(left_bc=='outflow'.or.left_bc=='reflecting','Unknown reconstruction left boundary')
      call require(right_bc=='outflow'.or.right_bc=='reflecting','Unknown reconstruction right boundary')
      if(left_bc=='reflecting') v(3,0)=-v(3,1)
      if(right_bc=='reflecting') v(3,nx+1)=-v(3,nx)
    end if
    do i=1,nx
      s=mc_slope(v(:,i)-v(:,i-1),v(:,i+1)-v(:,i))
      ! Close the largest species slope so sum(dY)=0. Do not perturb constant/absent species.
      j=maxloc(v(4:,i),dim=1)+3
      s(j)=0; s(j)=-sum(s(4:))
      lo=min(v(4:,i-1),v(4:,i),v(4:,i+1))
      hi=max(v(4:,i-1),v(4:,i),v(4:,i+1))
      theta=1
      do j=1,ns
        if(abs(s(j+3))>0) then
          bound=2*min(v(j+3,i)-lo(j),hi(j)-v(j+3,i))/abs(s(j+3))
          theta=min(theta,max(0._dp,bound))
        end if
      end do
      s(4:)=s(4:)*theta*(1-16*epsilon(theta))
      call primitive_to_conserved(m,v(1,i)-s(1)/2,v(2,i)-s(2)/2,v(3,i)-s(3)/2,v(4:,i)-s(4:)/2,qminus(:,i))
      call primitive_to_conserved(m,v(1,i)+s(1)/2,v(2,i)+s(2)/2,v(3,i)+s(3)/2,v(4:,i)+s(4:)/2,qplus(:,i))
    end do
  end subroutine

  subroutine primitive_to_conserved(m,t,p,u,y,q)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t,p,u,y(:)
    real(dp), intent(out) :: q(:)
    real(dp) :: cp,cv,h,e,r,rho
    integer :: ns
    ns=size(m%species)
    call require(size(q)==ns+2,'1D conserved state dimension mismatch')
    call require(ieee_is_finite(u),'Nonfinite velocity')
    call mixture(m,t,y,p,cp,cv,h,e,r)
    rho=p/(r*t)
    q(:ns)=rho*y; q(ns+1)=rho*u; q(ns+2)=rho*(e+u*u/2)
  end subroutine

  subroutine conserved_to_primitive(m,q,rho,u,t,p,sound,y,ok)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:)
    real(dp), intent(out) :: rho,u,t,p,sound,y(:)
    logical, optional, intent(out) :: ok
    logical :: valid
    real(dp) :: cp,cv,h,e,r
    integer :: ns
    ns=size(m%species)
    if(present(ok)) ok=.false.
    rho=0; u=0; t=0; p=0; sound=0; y=0
    call require(size(q)==ns+2.and.size(y)==ns,'1D state dimension mismatch')
    if(.not.all(ieee_is_finite(q)).or.any(q(:ns)<0)) then
      if(present(ok)) return
      call require(.false.,'Invalid conserved species; no clipping')
    end if
    rho=sum(q(:ns))
    if(.not.ieee_is_finite(rho).or.rho<=0) then
      if(present(ok)) return
      call require(.false.,'Nonpositive/nonfinite total density')
    end if
    y=q(:ns)/rho; u=q(ns+1)/rho; e=q(ns+2)/rho-u*u/2
    if(present(ok)) then
      t=temperature_from_energy(m,e,y,ok=valid)
      if(.not.valid) return
    else
      t=temperature_from_energy(m,e,y)
    end if
    call mixture(m,t,y,101325._dp,cp,cv,h,e,r)
    p=rho*r*t; sound=sqrt(cp/cv*r*t)
    valid=p>0.and.all(ieee_is_finite([u,p,sound]))
    if(present(ok)) then
      ok=valid
    else
      call require(valid,'Invalid pressure/sound speed')
    end if
  end subroutine

  logical function admissible_flow(m,q) result(valid)
    ! No modification, clipping, or renormalization of the candidate conserved field.
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:)
    real(dp) :: rho,u,t,p,a,y(size(m%species))
    integer :: i
    valid=.false.
    call require(size(q,1)==size(m%species)+2.and.size(q,2)>=2,'Invalid admissibility field shape')
    do i=1,size(q,2)
      call conserved_to_primitive(m,q(:,i),rho,u,t,p,a,y,ok=valid)
      if(.not.valid) return
    end do
  end function

  subroutine physical_flux(m,q,f,speed)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:)
    real(dp), intent(out) :: f(:),speed
    real(dp) :: rho,u,t,p,a,y(size(m%species))
    integer :: ns
    ns=size(m%species)
    call conserved_to_primitive(m,q,rho,u,t,p,a,y)
    f=q*u; f(ns+1)=f(ns+1)+p; f(ns+2)=(q(ns+2)+p)*u
    speed=abs(u)+a
  end subroutine

  subroutine rusanov_flux(m,left,right,f)
    ! Same physical Rusanov formula as legacy mod_mc_euler_flux; independent EOS/layout.
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: left(:),right(:)
    real(dp), intent(out) :: f(:)
    real(dp) :: fl(size(left)),fr(size(left)),sl,sr
    call physical_flux(m,left,fl,sl)
    call physical_flux(m,right,fr,sr)
    f=(fl+fr-max(sl,sr)*(right-left))/2
  end subroutine

  real(dp) function flow_timestep(m,q,dx,cfl,transport,reconstruction,left_bc,right_bc) result(dt)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx,cfl
    type(rf_transport), intent(in), optional :: transport
    character(*), intent(in), optional :: reconstruction,left_bc,right_bc
    character(16) :: lb,rb
    real(dp) :: qm(size(q,1),size(q,2)),qp(size(q,1),size(q,2))
    real(dp) :: f(size(q,1)),speed,maxspeed
    real(dp) :: rho,u,t,p,a,y(size(m%species)),cp,cv,h,e,r,rhomin,rhomax,rhocvmin,diff
    integer :: i
    call require(dx>0.and.cfl>0.and.cfl<=.5_dp,'Require dx>0 and 0<CFL<=0.5')
    call validate_reconstruction(reconstruction)
    maxspeed=0
    do i=1,size(q,2)
      call physical_flux(m,q(:,i),f,speed)
      maxspeed=max(maxspeed,speed)
    end do
    if(present(reconstruction)) then
      if(reconstruction=='muscl') then
        lb='outflow'; rb='outflow'
        if(present(left_bc)) lb=left_bc
        if(present(right_bc)) rb=right_bc
        call reconstruct_faces(m,q,lb,rb,qm,qp,reconstruction)
        do i=1,size(q,2)
          call physical_flux(m,qm(:,i),f,speed)
          maxspeed=max(maxspeed,speed)
          call physical_flux(m,qp(:,i),f,speed)
          maxspeed=max(maxspeed,speed)
        end do
      end if
    end if
    dt=cfl*dx/maxspeed
    if(present(transport)) then
      call validate_transport(transport)
      if(transport_active(transport)) then
        rhomin=huge(rhomin); rhomax=0; rhocvmin=huge(rhocvmin)
        do i=1,size(q,2)
          call conserved_to_primitive(m,q(:,i),rho,u,t,p,a,y)
          call mixture(m,t,y,p,cp,cv,h,e,r)
          rhomin=min(rhomin,rho); rhomax=max(rhomax,rho); rhocvmin=min(rhocvmin,rho*cv)
        end do
        ! Conservative explicit convection/diffusion estimate (cv, not cp, for compressible energy).
        diff=(4._dp/3*transport%viscosity+transport%bulk_viscosity)/rhomin &
              +transport%conductivity/rhocvmin+transport%diffusivity*rhomax/rhomin
        dt=cfl/(maxspeed/dx+2*diff/dx**2)
      end if
    end if
  end function

  subroutine boundary_state(q,kind,ghost)
    real(dp), intent(in) :: q(:)
    character(*), intent(in) :: kind
    real(dp), intent(out) :: ghost(:)
    ghost=q
    select case(kind)
    case('outflow') ! zero-gradient extrapolation, NOT a nonreflecting boundary
    case('reflecting')
      ghost(size(q)-1)=-ghost(size(q)-1)
    case default
      call require(.false.,'Unknown 1D boundary: '//kind)
    end select
  end subroutine

  subroutine diffusion_rhs(m,q,dx,left_bc,right_bc,transport,dq,net_boundary)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx
    character(*), intent(in) :: left_bc,right_bc
    type(rf_transport), intent(in) :: transport
    real(dp), intent(out) :: dq(:,:),net_boundary(:)
    real(dp) :: face(size(q,1),0:size(q,2)),rho(size(q,2)),u(size(q,2)),t(size(q,2))
    real(dp) :: y(size(m%species),size(q,2)),p,a
    integer :: i,nx
    nx=size(q,2)
    call validate_transport(transport)
    call require(nx>=2.and.size(q,1)==size(m%species)+2,'Invalid diffusion field shape')
    call require(all(shape(dq)==shape(q)).and.size(net_boundary)==size(q,1),'Invalid diffusion output shape')
    call require(dx>0.and.ieee_is_finite(dx),'Invalid diffusion spacing')
    call require((left_bc=='periodic').eqv.(right_bc=='periodic'),'Periodic boundary must be paired')
    dq=0; net_boundary=0
    if(.not.transport_active(transport)) return
    do i=1,nx
      call conserved_to_primitive(m,q(:,i),rho(i),u(i),t(i),p,a,y(:,i))
    end do
    do i=1,nx-1
      call diffusive_flux(m,transport,rho(i),u(i),t(i),y(:,i),rho(i+1),u(i+1),t(i+1),y(:,i+1),dx,face(:,i))
    end do
    face(:,0)=0; face(:,nx)=0
    if(left_bc=='periodic') then
      call diffusive_flux(m,transport,rho(nx),u(nx),t(nx),y(:,nx),rho(1),u(1),t(1),y(:,1),dx,face(:,0))
      face(:,nx)=face(:,0)
    else
      ! Reflecting: u_wall=0, zero temperature/species normal gradients.
      ! Outflow: zero primitive gradients, hence zero diffusive flux.
      select case(left_bc)
      case('reflecting')
        call diffusive_flux(m,transport,rho(1),-u(1),t(1),y(:,1),rho(1),u(1),t(1),y(:,1),dx,face(:,0))
      case('outflow')
      case default
        call require(.false.,'Unknown diffusion left boundary')
      end select
      select case(right_bc)
      case('reflecting')
        call diffusive_flux(m,transport,rho(nx),u(nx),t(nx),y(:,nx),rho(nx),-u(nx),t(nx),y(:,nx),dx,face(:,nx))
      case('outflow')
      case default
        call require(.false.,'Unknown diffusion right boundary')
      end select
    end if
    do i=1,nx
      dq(:,i)=-(face(:,i)-face(:,i-1))/dx
    end do
    net_boundary=face(:,0)-face(:,nx)
  end subroutine

  subroutine rhs(m,q,dx,left_bc,right_bc,dq,net_boundary,transport,reconstruction)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx
    character(*), intent(in) :: left_bc,right_bc
    real(dp), intent(out) :: dq(:,:),net_boundary(:)
    type(rf_transport), intent(in), optional :: transport
    real(dp) :: faces(size(q,1),0:size(q,2)),ghost(size(q,1))
    real(dp) :: diffusion(size(q,1),size(q,2)),diff_boundary(size(q,1))
    character(*), intent(in), optional :: reconstruction
    real(dp) :: qm(size(q,1),size(q,2)),qp(size(q,1),size(q,2))
    integer :: i,nx
    nx=size(q,2)
    call reconstruct_faces(m,q,left_bc,right_bc,qm,qp,reconstruction)
    call require((left_bc=='periodic').eqv.(right_bc=='periodic'),'Periodic boundary must be paired')
    if(left_bc=='periodic') then
      call rusanov_flux(m,qp(:,nx),qm(:,1),faces(:,0))
      faces(:,nx)=faces(:,0)
    else
      call boundary_state(qm(:,1),left_bc,ghost)
      call rusanov_flux(m,ghost,qm(:,1),faces(:,0))
      call boundary_state(qp(:,nx),right_bc,ghost)
      call rusanov_flux(m,qp(:,nx),ghost,faces(:,nx))
    end if
    do i=1,nx-1
      call rusanov_flux(m,qp(:,i),qm(:,i+1),faces(:,i))
    end do
    do i=1,nx
      dq(:,i)=-(faces(:,i)-faces(:,i-1))/dx
    end do
    net_boundary=faces(:,0)-faces(:,nx)
    if(present(transport)) then
      call diffusion_rhs(m,q,dx,left_bc,right_bc,transport,diffusion,diff_boundary)
      dq=dq+diffusion; net_boundary=net_boundary+diff_boundary
    end if
  end subroutine

  subroutine transport_step(m,q,dx,dt,cfl,left_bc,right_bc,result,boundary_change,ok,transport,reconstruction)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx,dt,cfl
    character(*), intent(in) :: left_bc,right_bc
    real(dp), intent(out) :: result(:,:),boundary_change(:)
    logical, intent(out) :: ok
    character(*), intent(in), optional :: reconstruction
    type(rf_transport), intent(in), optional :: transport
    real(dp) :: a(size(q,1),size(q,2)),b(size(q,1),size(q,2)),dq(size(q,1),size(q,2))
    real(dp) :: f1(size(q,1)),f2(size(q,1)),f3(size(q,1))
    ok=.false.; boundary_change=0
    result=q
    if(.not.admissible_flow(m,q)) return
    if(dt>flow_timestep(m,q,dx,cfl,transport,reconstruction,left_bc,right_bc)*(1+1.e-12_dp)) return
    call rhs(m,q,dx,left_bc,right_bc,dq,f1,transport,reconstruction)
    a=q+dt*dq
    if(.not.admissible_flow(m,a)) return
    if(dt>flow_timestep(m,a,dx,cfl,transport,reconstruction,left_bc,right_bc)*(1+1.e-12_dp)) return
    call rhs(m,a,dx,left_bc,right_bc,dq,f2,transport,reconstruction)
    b=.75_dp*q+.25_dp*(a+dt*dq)
    if(.not.admissible_flow(m,b)) return
    if(dt>flow_timestep(m,b,dx,cfl,transport,reconstruction,left_bc,right_bc)*(1+1.e-12_dp)) return
    call rhs(m,b,dx,left_bc,right_bc,dq,f3,transport,reconstruction)
    result=q/3+2._dp/3*(b+dt*dq)
    if(.not.admissible_flow(m,result)) then
      result=q
      return
    end if
    boundary_change=dt*(f1/6+f2/6+2._dp/3*f3)
    ok=.true.
  end subroutine

  subroutine chemistry_cells(m,q,dt,rtol,atoly,atolt,maxsteps)
    type(rf_mechanism), intent(in), target :: m
    real(dp), intent(inout) :: q(:,:)
    real(dp), intent(in) :: dt,rtol,atoly,atolt
    integer, intent(in) :: maxsteps
    real(dp) :: rho,u,t,p,a,y(size(m%species)),elem0(m%ne),elem1(m%ne),target_e,restored
    integer :: i,j,ns
    ns=size(m%species)
    do i=1,size(q,2)
      call conserved_to_primitive(m,q(:,i),rho,u,t,p,a,y)
      target_e=q(ns+2,i)/rho-u*u/2
      elem0=0
      do j=1,ns
        elem0=elem0+y(j)/m%species(j)%mass*m%species(j)%atoms
      end do
      call advance_chemistry(m,rho,t,y,dt,rtol,atoly,atolt,maxsteps)
      elem1=0
      do j=1,ns
        elem1=elem1+y(j)/m%species(j)%mass*m%species(j)%atoms
      end do
      call require(maxval(abs(elem1-elem0)/max(1._dp,abs(elem0)))<=1.e-8_dp,'Cell chemistry element drift')
      restored=temperature_from_energy(m,target_e,y)
      call require(abs(restored-t)<=max(1.e-3_dp,100*rtol*t),'Chemistry energy drift too large')
      ! No heat-release source: formation energy is already part of rho*E.
      q(:ns,i)=rho*y
    end do
  end subroutine

  subroutine advance_flow(m,q,dx,dt,cfl,left_bc,right_bc,chemistry,rtol,atoly,atolt,maxsteps,boundary_change, &
                          transport,reconstruction,rejected_steps)
    type(rf_mechanism), intent(in), target :: m
    real(dp), intent(inout) :: q(:,:),dt
    real(dp), intent(in) :: dx,cfl,rtol,atoly,atolt
    character(*), intent(in) :: left_bc,right_bc
    logical, intent(in) :: chemistry
    integer, intent(in) :: maxsteps
    real(dp), intent(out) :: boundary_change(:)
    character(*), intent(in), optional :: reconstruction
    type(rf_transport), intent(in), optional :: transport
    integer, intent(out), optional :: rejected_steps
    real(dp) :: old(size(q,1),size(q,2)),stage(size(q,1),size(q,2)),newq(size(q,1),size(q,2)),check_dt
    logical :: ok
    integer :: retry
    call require(ieee_is_finite(dt).and.dt>0,'Require positive finite flow dt')
    call require(size(q,1)==size(m%species)+2.and.size(q,2)>=2,'Invalid 1D field shape')
    call require(admissible_flow(m,q),'Invalid input flow state; cannot recover by reducing dt')
    if(present(rejected_steps)) rejected_steps=0
    old=q
    do retry=1,30
      stage=old
      if(chemistry) call chemistry_cells(m,stage,dt/2,rtol,atoly,atolt,maxsteps)
      call transport_step(m,stage,dx,dt,cfl,left_bc,right_bc,newq,boundary_change,ok,transport,reconstruction)
      if(ok) then
        if(chemistry) call chemistry_cells(m,newq,dt/2,rtol,atoly,atolt,maxsteps)
        check_dt=flow_timestep(m,newq,dx,cfl,transport,reconstruction,left_bc,right_bc) ! validate complete step
        call require(check_dt>0,'Invalid final flow state')
        q=newq
        if(present(rejected_steps)) rejected_steps=retry-1
        return
      end if
      dt=dt/2
      call require(dt>0.and.ieee_is_finite(dt),'Flow retry timestep underflow')
    end do
    call require(.false.,'Flow step CFL/admissibility retry limit exceeded (30 attempts)')
  end subroutine
end module
