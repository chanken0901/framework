module mod_rf_flow1d
  use mod_rf_reactor, only: advance_chemistry
  use mod_rf_thermo
  implicit none
  private
  public :: primitive_to_conserved,conserved_to_primitive,physical_flux,rusanov_flux
  public :: flow_timestep,advance_flow,transport_step
contains
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

  subroutine conserved_to_primitive(m,q,rho,u,t,p,sound,y)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:)
    real(dp), intent(out) :: rho,u,t,p,sound,y(:)
    real(dp) :: cp,cv,h,e,r
    integer :: ns
    ns=size(m%species)
    call require(size(q)==ns+2.and.size(y)==ns,'1D state dimension mismatch')
    call require(all(ieee_is_finite(q)).and.all(q(:ns)>=0),'Invalid conserved species; no clipping')
    rho=sum(q(:ns)); call require(rho>0,'Nonpositive total density')
    y=q(:ns)/rho; u=q(ns+1)/rho; e=q(ns+2)/rho-u*u/2
    t=temperature_from_energy(m,e,y)
    call mixture(m,t,y,101325._dp,cp,cv,h,e,r)
    p=rho*r*t; sound=sqrt(cp/cv*r*t)
    call require(p>0.and.ieee_is_finite(sound),'Invalid pressure/sound speed')
  end subroutine

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

  real(dp) function flow_timestep(m,q,dx,cfl) result(dt)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx,cfl
    real(dp) :: f(size(q,1)),speed,maxspeed
    integer :: i
    call require(dx>0.and.cfl>0.and.cfl<=.5_dp,'Require dx>0 and 0<CFL<=0.5')
    maxspeed=0
    do i=1,size(q,2)
      call physical_flux(m,q(:,i),f,speed)
      maxspeed=max(maxspeed,speed)
    end do
    dt=cfl*dx/maxspeed
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

  subroutine rhs(m,q,dx,left_bc,right_bc,dq,net_boundary)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx
    character(*), intent(in) :: left_bc,right_bc
    real(dp), intent(out) :: dq(:,:),net_boundary(:)
    real(dp) :: faces(size(q,1),0:size(q,2)),ghost(size(q,1))
    integer :: i,nx
    nx=size(q,2)
    call require((left_bc=='periodic').eqv.(right_bc=='periodic'),'Periodic boundary must be paired')
    if(left_bc=='periodic') then
      call rusanov_flux(m,q(:,nx),q(:,1),faces(:,0))
      faces(:,nx)=faces(:,0)
    else
      call boundary_state(q(:,1),left_bc,ghost)
      call rusanov_flux(m,ghost,q(:,1),faces(:,0))
      call boundary_state(q(:,nx),right_bc,ghost)
      call rusanov_flux(m,q(:,nx),ghost,faces(:,nx))
    end if
    do i=1,nx-1
      call rusanov_flux(m,q(:,i),q(:,i+1),faces(:,i))
    end do
    do i=1,nx
      dq(:,i)=-(faces(:,i)-faces(:,i-1))/dx
    end do
    net_boundary=faces(:,0)-faces(:,nx)
  end subroutine

  subroutine transport_step(m,q,dx,dt,cfl,left_bc,right_bc,result,boundary_change,ok)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:),dx,dt,cfl
    character(*), intent(in) :: left_bc,right_bc
    real(dp), intent(out) :: result(:,:),boundary_change(:)
    logical, intent(out) :: ok
    real(dp) :: a(size(q,1),size(q,2)),b(size(q,1),size(q,2)),dq(size(q,1),size(q,2))
    real(dp) :: f1(size(q,1)),f2(size(q,1)),f3(size(q,1))
    ok=.false.; boundary_change=0
    if(dt>flow_timestep(m,q,dx,cfl)*(1+1.e-12_dp)) return
    call rhs(m,q,dx,left_bc,right_bc,dq,f1)
    a=q+dt*dq
    if(dt>flow_timestep(m,a,dx,cfl)*(1+1.e-12_dp)) return
    call rhs(m,a,dx,left_bc,right_bc,dq,f2)
    b=.75_dp*q+.25_dp*(a+dt*dq)
    if(dt>flow_timestep(m,b,dx,cfl)*(1+1.e-12_dp)) return
    call rhs(m,b,dx,left_bc,right_bc,dq,f3)
    result=q/3+2._dp/3*(b+dt*dq)
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

  subroutine advance_flow(m,q,dx,dt,cfl,left_bc,right_bc,chemistry,rtol,atoly,atolt,maxsteps,boundary_change)
    type(rf_mechanism), intent(in), target :: m
    real(dp), intent(inout) :: q(:,:),dt
    real(dp), intent(in) :: dx,cfl,rtol,atoly,atolt
    character(*), intent(in) :: left_bc,right_bc
    logical, intent(in) :: chemistry
    integer, intent(in) :: maxsteps
    real(dp), intent(out) :: boundary_change(:)
    real(dp) :: old(size(q,1),size(q,2)),stage(size(q,1),size(q,2)),newq(size(q,1),size(q,2)),check_dt
    logical :: ok
    integer :: retry
    call require(ieee_is_finite(dt).and.dt>0,'Require positive finite flow dt')
    call require(size(q,1)==size(m%species)+2.and.size(q,2)>=2,'Invalid 1D field shape')
    old=q
    do retry=1,30
      stage=old
      if(chemistry) call chemistry_cells(m,stage,dt/2,rtol,atoly,atolt,maxsteps)
      call transport_step(m,stage,dx,dt,cfl,left_bc,right_bc,newq,boundary_change,ok)
      if(ok) then
        if(chemistry) call chemistry_cells(m,newq,dt/2,rtol,atoly,atolt,maxsteps)
        check_dt=flow_timestep(m,newq,dx,cfl) ! validate the complete step
        call require(check_dt>0,'Invalid final flow state')
        q=newq
        return
      end if
      dt=dt/2
    end do
    call require(.false.,'Flow step CFL retry limit exceeded')
  end subroutine
end module
