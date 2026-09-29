module mod_rf_equilibrium
  use mod_rf_thermo
  implicit none
  private
  public :: equilibrium_tv,equilibrium_sound
contains
  function equilibrium_sound(m,t,rho,y0) result(sound)
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t,rho,y0(:)
    real(dp) :: sound,pp(4),ss(4),tt,rr,y(size(y0)),cp,cv,h,e,r,pt,pr,st,sr
    real(dp), parameter :: delta=1.e-4_dp
    integer :: i
    do i=1,4
      tt=t;rr=rho
      select case(i)
      case(1);tt=t*(1-delta)
      case(2);tt=t*(1+delta)
      case(3);rr=rho*(1-delta)
      case(4);rr=rho*(1+delta)
      end select
      call equilibrium_tv(m,tt,rr,y0,y)
      call mixture(m,tt,y,101325._dp,cp,cv,h,e,r)
      pp(i)=rr*r*tt
      call mixture(m,tt,y,pp(i),cp,cv,h,e,r,ss(i))
    end do
    pt=(pp(2)-pp(1))/(2*delta*t);pr=(pp(4)-pp(3))/(2*delta*rho)
    st=(ss(2)-ss(1))/(2*delta*t);sr=(ss(4)-ss(3))/(2*delta*rho)
    call require(st>0,'Invalid equilibrium entropy derivative')
    sound=pr-pt*sr/st
    call require(sound>0.and.ieee_is_finite(sound),'Invalid equilibrium sound speed')
    sound=sqrt(sound)
  end function

  subroutine linear_solve(a,b,x)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(:)
    real(dp) :: c(size(b),size(b)),v(size(b)),row(size(b)),q
    integer :: i,j,k,n
    n=size(b);c=a;v=b
    do i=1,n
      k=i-1+maxloc(abs(c(i:,i)),dim=1)
      call require(abs(c(k,i))>1.e-14_dp*maxval(abs(a)),'Singular equilibrium element system')
      row=c(i,:);c(i,:)=c(k,:);c(k,:)=row;q=v(i);v(i)=v(k);v(k)=q
      do j=i+1,n
        q=c(j,i)/c(i,i);c(j,i:)=c(j,i:)-q*c(i,i:);v(j)=v(j)-q*v(i)
      end do
    end do
    do i=n,1,-1
      x(i)=(v(i)-dot_product(c(i,i+1:),x(i+1:)))/c(i,i)
    end do
    call require(all(ieee_is_finite(x)),'Nonfinite equilibrium Newton step')
  end subroutine

  subroutine equilibrium_tv(m,t,rho,y0,y)
    ! Ideal neutral gas equilibrium: element potentials at fixed T and density.
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: t,rho,y0(:)
    real(dp), intent(out) :: y(:)
    real(dp) :: bfull(m%ne),a(m%ne,size(y)),b(m%ne),lambda(m%ne),step(m%ne)
    real(dp) :: jac(m%ne,m%ne),res(m%ne),trial(m%ne),n(size(y)),nt(size(y)),base(size(y))
    real(dp) :: cp,h,s,alpha,err,newerr,scale,normal(m%ne,m%ne),rhs(m%ne)
    integer :: elements(m%ne),species(size(y)),ne,ns,i,j,k,it,ls
    call check_y(m,y0)
    call require(ieee_is_finite(rho).and.rho>0,'Invalid equilibrium density')
    bfull=0
    do i=1,size(y)
      bfull=bfull+y0(i)/m%species(i)%mass*m%species(i)%atoms
    end do
    ne=0
    do i=1,m%ne
      if(bfull(i)<=0) cycle
      ne=ne+1;elements(ne)=i;b(ne)=bfull(i)
    end do
    ns=0
    do i=1,size(y)
      if(any(m%species(i)%atoms>0.and.bfull==0)) cycle
      ns=ns+1;species(ns)=i
      a(:ne,ns)=m%species(i)%atoms(elements(:ne))
      call species_thermo(m%species(i),t,cp,h,s)
      base(ns)=log(m%species(i)%pref/(rho*gas_r*t))-(h-t*s)/(gas_r*t)
    end do
    call require(ne>0.and.ns>=ne,'Insufficient equilibrium species/element rank')
    ! Fit element potentials to moderate initial mole amounts, not to zero trace Y.
    normal(:ne,:ne)=matmul(a(:ne,:ns),transpose(a(:ne,:ns)))
    rhs(:ne)=matmul(a(:ne,:ns),log(sum(b(:ne))/ns)-base(:ns))
    call linear_solve(normal(:ne,:ne),rhs(:ne),lambda(:ne))
    do it=1,500
      call amounts(lambda(:ne),n(:ns))
      res(:ne)=matmul(a(:ne,:ns),n(:ns))-b(:ne)
      err=maxval(abs(res(:ne))/b(:ne))
      if(err<1.e-11_dp) exit
      do i=1,ne
        do j=1,ne
          jac(i,j)=sum(a(i,:ns)*a(j,:ns)*n(:ns))
        end do
      end do
      ! Row/column scaling by elemental abundance.
      do i=1,ne
        rhs(i)=-res(i)/sqrt(b(i))
        do j=1,ne
          normal(i,j)=jac(i,j)/sqrt(b(i)*b(j))
        end do
      end do
      ! Trace species can make the Hessian numerically rank deficient at low T.
      ! Regularize the Newton direction only; convergence uses unmodified balances.
      scale=1.e-12_dp*max(1._dp,maxval(abs(normal(:ne,:ne))))
      do i=1,ne
        normal(i,i)=normal(i,i)+scale
      end do
      call linear_solve(normal(:ne,:ne),rhs(:ne),step(:ne))
      step(:ne)=step(:ne)/sqrt(b(:ne))
      scale=maxval(abs(matmul(transpose(a(:ne,:ns)),step(:ne))))
      alpha=min(1._dp,2._dp/max(scale,tiny(scale)))
      do ls=1,60
        trial(:ne)=lambda(:ne)+alpha*step(:ne)
        call amounts(trial(:ne),nt(:ns))
        newerr=maxval(abs(matmul(a(:ne,:ns),nt(:ns))-b(:ne))/b(:ne))
        if(newerr<err) exit
        alpha=alpha*.5_dp
      end do
      call require(ls<=60,'Equilibrium line search failed')
      lambda(:ne)=trial(:ne)
    end do
    call require(it<=500,'Equilibrium failed to converge')
    y=0
    do k=1,ns
      y(species(k))=n(k)*m%species(species(k))%mass
    end do
    call require(abs(sum(y)-1)<1.e-9_dp,'Equilibrium mass closure failed')
    y=y/sum(y)
  contains
    subroutine amounts(potential,values)
      real(dp), intent(in) :: potential(:)
      real(dp), intent(out) :: values(:)
      real(dp) :: ex(size(values))
      ex=base(:ns)+matmul(transpose(a(:ne,:ns)),potential)
      call require(maxval(ex)<600,'Equilibrium exponential overflow')
      values=exp(max(ex,-700._dp))
    end subroutine
  end subroutine
end module
