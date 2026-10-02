module mod_rf_flow2d
  use mod_rf_thermo
  use mod_rf_grid2d
  use mod_rf_finite_volume
  use mod_rf_flow1d, only: chemistry_cells
  implicit none
  private
  public :: rf_boundary2d,flow2d_rhs,advance_flow2d
  type :: rf_boundary2d
    character(16) :: kind(4)=[character(16)::'outflow','outflow','reflecting','reflecting']
    real(dp), allocatable :: fixed(:,:)
  end type
contains
  logical function admissible(m,q) result(ok)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: q(:,:)
    real(dp) :: rho,v(2),t,p,a,y(size(m%species))
    integer :: c
    ok=.false.
    do c=1,size(q,2)
      call conserved_nd(m,q(:,c),rho,v,t,p,a,y,ok)
      if(.not.ok) return
    end do
  end function

  subroutine flow2d_rhs(m,grid,bc,q,cfl,dq,rate,dt,ok)
    type(rf_mechanism), intent(in) :: m
    type(rf_grid2d), intent(in) :: grid
    type(rf_boundary2d), intent(in) :: bc
    real(dp), intent(in) :: q(:,:),cfl
    real(dp), intent(out) :: dq(:,:),rate(:),dt
    logical, intent(out) :: ok
    real(dp), allocatable :: ghost(:,:)
    real(dp) :: rho,v(2),t,p,a,y(size(m%species))
    integer :: f,b,l,nv
    nv=size(m%species)+3
    call validate_face_mesh(grid%mesh)
    call require(size(grid%mesh%area_vector,1)==2,'flow2d requires a 2D mesh')
    call require(size(q,1)==nv.and.size(q,2)==size(grid%mesh%volume),'Invalid 2D field shape')
    call require(allocated(grid%boundary),'Missing boundary tags')
    call require(size(grid%boundary)==size(grid%mesh%owner),'Invalid boundary tags')
    call require(all(shape(dq)==shape(q)).and.size(rate)==nv,'Invalid 2D residual shape')
    do b=1,4
      call require(bc%kind(b)=='outflow'.or.bc%kind(b)=='reflecting'.or.bc%kind(b)=='dirichlet', &
        'Unsupported 2D boundary (periodic/characteristic not implemented)')
      if(bc%kind(b)=='dirichlet') then
        call require(allocated(bc%fixed),'Missing fixed boundary states')
        call require(size(bc%fixed,1)==nv.and.size(bc%fixed,2)==4,'Invalid fixed boundary shape')
        call conserved_nd(m,bc%fixed(:,b),rho,v,t,p,a,y)
      end if
    end do
    dq=0;rate=0;dt=0;ok=admissible(m,q)
    if(.not.ok) return
    allocate(ghost(nv,size(grid%boundary)));ghost=0
    do f=1,size(grid%boundary)
      b=grid%boundary(f)
      if(grid%mesh%neighbor(f)>0) then
        call require(b==0,'Internal face has a boundary tag')
        cycle
      end if
      call require(b>=1.and.b<=4,'Boundary face lacks a valid tag')
      l=grid%mesh%owner(f)
      select case(bc%kind(b))
      case('outflow')
        ghost(:,f)=q(:,l)
      case('reflecting')
        call reflect_normal(m,q(:,l),grid%mesh%area_vector(:,f),ghost(:,f))
      case('dirichlet')
        ghost(:,f)=bc%fixed(:,b)
      end select
      call conserved_nd(m,ghost(:,f),rho,v,t,p,a,y,ok)
      if(.not.ok) return
    end do
    call finite_volume_rhs(m,grid%mesh,q,ghost,cfl,dq,rate,dt)
    ok=all(ieee_is_finite(dq)).and.all(ieee_is_finite(rate))
  end subroutine

  subroutine chemistry2d(m,q,dt,rtol,atoly,atolt,maxsteps,ok)
    type(rf_mechanism), intent(in), target :: m
    real(dp), intent(inout) :: q(:,:)
    real(dp), intent(in) :: dt,rtol,atoly,atolt
    integer, intent(in) :: maxsteps
    logical, intent(out) :: ok
    real(dp) :: line(size(m%species)+2,size(q,2))
    integer :: ns,c
    ns=size(m%species)
    line(:ns,:)=q(:ns,:);line(ns+2,:)=q(ns+3,:)
    do c=1,size(q,2)
      line(ns+1,c)=norm2(q(ns+1:ns+2,c))
    end do
    ! Same total kinetic energy and constant-volume chemistry as 1D.
    ! Momentum components and total (formation-inclusive) energy are not modified.
    call chemistry_cells(m,line,dt,rtol,atoly,atolt,maxsteps,ok)
    if(ok) q(:ns,:)=line(:ns,:)
  end subroutine

  subroutine advance_flow2d(m,grid,bc,q,dt,cfl,chemistry,rtol,atoly,atolt,maxsteps,change,rejected)
    type(rf_mechanism), intent(in), target :: m
    type(rf_grid2d), intent(in) :: grid
    type(rf_boundary2d), intent(in) :: bc
    real(dp), intent(inout) :: q(:,:),dt
    real(dp), intent(in) :: cfl,rtol,atoly,atolt
    logical, intent(in) :: chemistry
    integer, intent(in) :: maxsteps
    real(dp), intent(out) :: change(:)
    integer, intent(out) :: rejected
    real(dp) :: base(size(q,1),size(q,2)),a(size(q,1),size(q,2)),b(size(q,1),size(q,2))
    real(dp) :: candidate(size(q,1),size(q,2)),dq(size(q,1),size(q,2))
    real(dp) :: f1(size(q,1)),f2(size(q,1)),f3(size(q,1)),unused(size(q,1)),allowed
    logical :: ok
    integer :: attempt
    call require(ieee_is_finite(dt).and.dt>0,'Invalid 2D timestep')
    call flow2d_rhs(m,grid,bc,q,cfl,dq,unused,allowed,ok)
    call require(ok,'Invalid initial 2D operator state')
    call require(size(change)==size(q,1),'Invalid boundary change shape')
    change=0;rejected=0
    do attempt=1,30
      base=q;ok=.true.
      if(chemistry) call chemistry2d(m,base,dt/2,rtol,atoly,atolt,maxsteps,ok)
      if(ok) then
        call flow2d_rhs(m,grid,bc,base,cfl,dq,f1,allowed,ok)
        if(ok) ok=dt<=allowed*(1+1.e-12_dp)
      end if
      if(ok) then
        a=base+dt*dq
        call flow2d_rhs(m,grid,bc,a,cfl,dq,f2,allowed,ok)
        if(ok) ok=dt<=allowed*(1+1.e-12_dp)
      end if
      if(ok) then
        b=base+.25_dp*((a-base)+dt*dq)
        call flow2d_rhs(m,grid,bc,b,cfl,dq,f3,allowed,ok)
        if(ok) ok=dt<=allowed*(1+1.e-12_dp)
      end if
      if(ok) then
        candidate=base+2._dp/3*((b-base)+dt*dq)
        ok=admissible(m,candidate)
      end if
      if(ok.and.chemistry) call chemistry2d(m,candidate,dt/2,rtol,atoly,atolt,maxsteps,ok)
      if(ok) call flow2d_rhs(m,grid,bc,candidate,cfl,dq,unused,allowed,ok)
      if(ok) then
        q=candidate;change=dt*(f1/6+f2/6+2._dp/3*f3);rejected=attempt-1
        return
      end if
      dt=dt/2
      call require(ieee_is_finite(dt).and.dt>0,'2D retry timestep underflow')
    end do
    call require(.false.,'2D transport/chemistry retry limit exceeded')
  end subroutine
end module
