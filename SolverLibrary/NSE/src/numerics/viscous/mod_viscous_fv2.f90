module mod_viscous_fv2
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_model_config, only: nse_config
  use mod_grid_fvm, only: axis_x,axis_y,axis_z
  implicit none
  private
  public :: add_viscous_fv2_rhs
contains
  pure real(dp) function coordinate(axis,index) result(x)
    integer, intent(in) :: axis,index
    select case(axis)
    case(1);x=axis_x%center(index)
    case(2);x=axis_y%center(index)
    case(3);x=axis_z%center(index)
    case default;error stop 'Invalid FV2 coordinate direction'
    end select
  end function

  pure subroutine primitive_at(q,p,sim,nse,js,ks,v)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: p(3),js,ks
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(out) :: v(4)
    real(dp) :: rho,pressure
    rho=max(q(p(1),p(2),p(3),1),nse%small_rho)
    v(1:3)=q(p(1),p(2),p(3),2:4)/rho
    pressure=max((nse%gamma-1)*(q(p(1),p(2),p(3),5)-.5_dp*rho*sum(v(1:3)**2)),nse%small_p)
    v(4)=pressure/rho
  end subroutine

  pure subroutine gradient_at(q,p,axis,sim,nse,js,ks,g)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: p(3),axis,js,ks
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(out) :: g(4)
    real(dp) :: vm(4),v0(4),vp(4),dl,dr
    integer :: pm(3),pp(3)
    pm=p;pp=p;pm(axis)=pm(axis)-1;pp(axis)=pp(axis)+1
    call primitive_at(q,pm,sim,nse,js,ks,vm)
    call primitive_at(q,p,sim,nse,js,ks,v0)
    call primitive_at(q,pp,sim,nse,js,ks,vp)
    dl=coordinate(axis,p(axis))-coordinate(axis,pm(axis))
    dr=coordinate(axis,pp(axis))-coordinate(axis,p(axis))
    g=(dr*(v0-vm)/dl+dl*(vp-v0)/dr)/(dl+dr)
  end subroutine

  pure subroutine face_flux(q,p,axis,sim,nse,js,ks,flux)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: p(3),axis,js,ks
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(out) :: flux(5)
    real(dp) :: vl(4),vr(4),vf(4),gl(4),gr(4),grad(4,3),tau(3),theta,xf,distance,divu
    integer :: right(3),a
    right=p;right(axis)=right(axis)+1
    call primitive_at(q,p,sim,nse,js,ks,vl)
    call primitive_at(q,right,sim,nse,js,ks,vr)
    distance=coordinate(axis,right(axis))-coordinate(axis,p(axis))
    select case(axis)
    case(1);xf=axis_x%edge(p(axis))
    case(2);xf=axis_y%edge(p(axis))
    case(3);xf=axis_z%edge(p(axis))
    case default;error stop 'Invalid FV2 flux direction'
    end select
    theta=(xf-coordinate(axis,p(axis)))/distance
    vf=(1-theta)*vl+theta*vr
    do a=1,3
      if(a==axis) then
        grad(:,a)=(vr-vl)/distance
      else
        call gradient_at(q,p,a,sim,nse,js,ks,gl)
        call gradient_at(q,right,a,sim,nse,js,ks,gr)
        grad(:,a)=(1-theta)*gl+theta*gr
      end if
    end do
    divu=grad(1,1)+grad(2,2)+grad(3,3)
    do a=1,3
      tau(a)=grad(a,axis)+grad(axis,a)
    end do
    tau(axis)=tau(axis)-2._dp/3*divu
    flux(1)=0;flux(2:4)=tau/nse%reynolds
    flux(5)=(dot_product(vf(1:3),tau)+nse%gamma/((nse%gamma-1)*nse%prandtl)*grad(4,axis))/nse%reynolds
  end subroutine

  subroutine add_viscous_fv2_rhs(q,rhs,sim,nse,js,je,ks,ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js,je,ks,ke
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(inout) :: rhs(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp) :: plus(5),minus(5),width(3)
    integer :: i,j,k,a,p(3),pm(3)
    if(sim%grid_mapping/='sinh') error stop 'fv2 currently requires the stretched grid geometry'
    if(.not.allocated(axis_x%center)) error stop 'fv2 geometry is not initialized'
    if(nse%nv/=5) error stop 'fv2 requires five conserved variables'
    !$OMP DO collapse(2) schedule(static)
    do k=ks,ke
      do j=js,je
        do i=1,sim%nx
          p=[i,j,k];width=[axis_x%width(i),axis_y%width(j),axis_z%width(k)]
          ! Use the same diagonal integration norm as mapped KEEP6 convection.
          ! The face gradient remains FV2 (not sixth-order diffusion).
          if(sim%mapped_keep6) width=[axis_x%keep6_metric(i),axis_y%keep6_metric(j),axis_z%keep6_metric(k)]
          do a=1,3
            pm=p;pm(a)=pm(a)-1
            call face_flux(q,p,a,sim,nse,js,ks,plus)
            call face_flux(q,pm,a,sim,nse,js,ks,minus)
            rhs(i,j,k,1:5)=rhs(i,j,k,1:5)+(plus-minus)/width(a)
          end do
        end do
      end do
    end do
    !$OMP END DO
  end subroutine
end module
