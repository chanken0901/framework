! Optional LLNS baseline. D+ gradients and -D+^T divergences form a matched
! dissipative/noise pair. This changes transport ONLY when explicitly enabled.
module mod_nse_fluctuating
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_model_config, only: nse_config
  use mod_fh_random, only: fh_normal
  implicit none
  private
  public :: validate_fh, add_fh_transport, add_fh_noise, fh_sample
  real(dp), allocatable, save :: flux(:,:,:,:,:)
contains
  subroutine validate_fh(sim,nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    if (.not. nse%fh_enabled) return
    if (.not. sim%use_fixed_dt) error stop 'fluctuating hydrodynamics requires fixed dt'
    if (.not. ieee_is_finite(sim%dt) .or. sim%dt<=0) error stop 'FH requires finite positive dt'
    if (any(nse%boundary_face_type/='periodic') .or. nse%boundary_condition/='periodic') &
      error stop 'fluctuating hydrodynamics requires all-periodic boundaries'
    if (nse%nv/=5 .or. sim%nghost<2) error stop 'FH requires five variables and at least two ghosts'
    if (min(sim%nx,sim%ny,sim%nz)<2 .or. min(sim%dx,sim%dy,sim%dz)<=0 .or. &
        .not. all(ieee_is_finite([sim%dx,sim%dy,sim%dz]))) error stop 'FH requires a positive 3D grid'
    if (nse%viscous_scheme/='central6' .or. nse%convective_scheme/='keep6') &
      error stop 'FH requires central6 build and KEEP6; transport is replaced by matched FH operator'
    if (.not. ieee_is_finite(nse%reynolds) .or. nse%reynolds<=0 .or. &
        .not. ieee_is_finite(nse%prandtl) .or. nse%prandtl<=0 .or. &
        .not. ieee_is_finite(nse%gamma) .or. nse%gamma<=1) error stop 'FH requires positive transport'
    if (.not. ieee_is_finite(nse%fh_boltzmann_number) .or. nse%fh_boltzmann_number<0 .or. &
        nse%fh_seed<0) error stop 'invalid FH Boltzmann number or seed'
  end subroutine

  subroutine workspace(sim,js,je,ks,ke)
    type(simulation_config), intent(in) :: sim
    integer, intent(in) :: js,je,ks,ke
    if (allocated(flux)) then
      if (all(shape(flux)==[sim%nx+2,je-js+3,ke-ks+3,4,3]) .and. &
          lbound(flux,2)==js-1 .and. lbound(flux,3)==ks-1) return
      deallocate(flux)
    end if
    allocate(flux(0:sim%nx+1,js-1:je+1,ks-1:ke+1,4,3))
  end subroutine

  pure subroutine primitive(s,gamma,u,t)
    real(dp), intent(in) :: s(:),gamma
    real(dp), intent(out) :: u(3),t
    u=s(2:4)/s(1)
    t=(gamma-1)*(s(5)/s(1)-0.5_dp*sum(u*u))
  end subroutine

  pure subroutine fh_sample(i,j,k,sim,nse,t,dt,stress,heat)
    integer, intent(in) :: i,j,k
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    real(dp), intent(in) :: t,dt
    real(dp), intent(out) :: stress(3,3),heat(3)
    real(dp) :: z(9),a,b,mu,kappa
    integer :: c,ii,jj,kk
    ii=modulo(i-1,sim%nx); jj=modulo(j-1,sim%ny); kk=modulo(k-1,sim%nz)
    do c=1,9
      z(c)=fh_normal(ii,jj,kk,sim%step,nse%fh_seed,c)
    end do
    mu=1.0_dp/nse%reynolds
    kappa=mu*nse%gamma/((nse%gamma-1)*nse%prandtl)
    a=sqrt(2*nse%fh_boltzmann_number*mu*t/(sim%dx*sim%dy*sim%dz*dt))
    b=sum(z(1:3))/3
    stress=0
    do c=1,3
      stress(c,c)=sqrt(2.0_dp)*a*(z(c)-b)
    end do
    stress(1,2)=a*z(4); stress(2,1)=stress(1,2)
    stress(1,3)=a*z(5); stress(3,1)=stress(1,3)
    stress(2,3)=a*z(6); stress(3,2)=stress(2,3)
    heat=sqrt(2*nse%fh_boltzmann_number*kappa*t*t/(sim%dx*sim%dy*sim%dz*dt))*z(7:9)
  end subroutine

  subroutine add_fh_transport(q,rhs,sim,nse,js,je,ks,ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js,je,ks,ke
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(inout) :: rhs(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    call add_flux(q,rhs,sim,nse,js,je,ks,ke,.false.,1.0_dp)
  end subroutine

  subroutine add_fh_noise(q,rhs,dt,sim,nse,js,je,ks,ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js,je,ks,ke
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(inout) :: rhs(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    if (.not. nse%fh_enabled .or. nse%fh_boltzmann_number==0) return
    call add_flux(q,rhs,sim,nse,js,je,ks,ke,.true.,dt)
  end subroutine

  subroutine add_flux(q,rhs,sim,nse,js,je,ks,ke,stochastic,dt)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js,je,ks,ke
    logical, intent(in) :: stochastic
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    real(dp), intent(inout) :: rhs(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    integer :: i,j,k,d,offset(3),p(3)
    real(dp) :: u(3),up(3),t,tp,h(3),grad(3,3),heat(3),stress(3,3),vel(3,3),mu,kappa,divu
    if (.not. nse%fh_enabled) return
    h=[sim%dx,sim%dy,sim%dz]; mu=1/nse%reynolds
    kappa=mu*nse%gamma/((nse%gamma-1)*nse%prandtl)
    !$OMP MASKED
    call workspace(sim,js,je,ks,ke)
    !$OMP END MASKED
    !$OMP BARRIER
    !$OMP DO collapse(2) schedule(static)
    do k=ks-1,ke+1
      do j=js-1,je+1
        do i=0,sim%nx+1
          call primitive(q(i,j,k,:),nse%gamma,u,t)
          do d=1,3
            offset=0; offset(d)=1; p=[i,j,k]+offset
            call primitive(q(p(1),p(2),p(3),:),nse%gamma,up,tp)
            grad(:,d)=(up-u)/h(d)
            heat(d)=kappa*(tp-t)/h(d)
            vel(:,d)=0.5_dp*(up+u)
          end do
          if (stochastic) then
            call fh_sample(i,j,k,sim,nse,t,dt,stress,heat)
          else
            divu=grad(1,1)+grad(2,2)+grad(3,3)
            stress=mu*(grad+transpose(grad))
            do d=1,3
              stress(d,d)=stress(d,d)-2*mu*divu/3
            end do
          end if
          do d=1,3
            flux(i,j,k,1:3,d)=stress(:,d)
            flux(i,j,k,4,d)=dot_product(vel(:,d),stress(:,d))+heat(d)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP DO collapse(2) schedule(static)
    do k=ks,ke
      do j=js,je
        do i=1,sim%nx
          rhs(i,j,k,2:5)=rhs(i,j,k,2:5)+(flux(i,j,k,:,1)-flux(i-1,j,k,:,1))/h(1) &
            +(flux(i,j,k,:,2)-flux(i,j-1,k,:,2))/h(2)+(flux(i,j,k,:,3)-flux(i,j,k-1,:,3))/h(3)
        end do
      end do
    end do
    !$OMP END DO
  end subroutine
end module
