module mod_reconstruction_nonuniform
  ! Cell-average polynomial reconstruction, distinct from point differentiation.
  use mod_precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: weno_face_geometry,build_weno_geometry,reconstruct_nonuniform_left
  type :: weno_face_geometry
    real(dp) :: polynomial(0:2,0:2,0:2)=0 ! degree, stencil cell, candidate
    real(dp) :: optimal(0:2)=0
  end type
contains
  subroutine build_weno_geometry(edges,geometry)
    ! Five cell averages; reconstruction at edges(3), from the left.
    real(dp), intent(in) :: edges(0:5)
    type(weno_face_geometry), intent(out) :: geometry
    real(dp) :: x(0:5),h,moment(0:4,0:4),inverse(0:4,0:4)
    real(dp) :: local(0:2,0:2),local_inverse(0:2,0:2),combined(0:4),full(0:4)
    integer :: r,i,k
    if(.not.all(ieee_is_finite(edges))) error stop 'Nonfinite WENO geometry'
    if(any(edges(1:)<=edges(:4))) error stop 'WENO edges must increase'
    h=edges(3)-edges(2);x=(edges-edges(3))/h
    do i=0,4
      do k=0,4
        moment(i,k)=(x(i+1)**(k+1)-x(i)**(k+1))/((k+1)*(x(i+1)-x(i)))
      end do
    end do
    call invert(moment,inverse);full=inverse(0,:)
    do r=0,2
      local=moment(r:r+2,0:2)
      call invert(local,local_inverse)
      geometry%polynomial(:,:,r)=local_inverse
    end do
    geometry%optimal(0)=full(0)/geometry%polynomial(0,0,0)
    geometry%optimal(2)=full(4)/geometry%polynomial(0,2,2)
    geometry%optimal(1)=1-geometry%optimal(0)-geometry%optimal(2)
    if(.not.all(ieee_is_finite(geometry%optimal)).or.any(geometry%optimal<=0)) &
      error stop 'WENO geometry requires positive optimal weights; reduce mesh stretching'
    combined=0
    do r=0,2
      combined(r:r+2)=combined(r:r+2)+geometry%optimal(r)*geometry%polynomial(0,:,r)
    end do
    if(maxval(abs(combined-full))>1.e-10_dp*max(1._dp,maxval(abs(full)))) &
      error stop 'Inconsistent nonuniform WENO optimal weights'
  end subroutine

  subroutine invert(a,inverse)
    real(dp), intent(in) :: a(0:,0:)
    real(dp), intent(out) :: inverse(0:,0:)
    real(dp) :: work(0:size(a,1)-1,0:2*size(a,1)-1),row(0:2*size(a,1)-1),pivot
    integer :: n,i,j,p
    n=size(a,1);work=0;work(:,0:n-1)=a
    if(.not.all(ieee_is_finite(a))) error stop 'Nonfinite WENO moments'
    do i=0,n-1
      work(i,n+i)=1
    end do
    do i=0,n-1
      p=i+maxloc(abs(work(i:n-1,i)),dim=1)-1
      pivot=work(p,i)
      if(abs(pivot)<=tiny(pivot)) error stop 'Singular WENO moment matrix'
      row=work(i,:);work(i,:)=work(p,:);work(p,:)=row
      work(i,:)=work(i,:)/pivot
      do j=0,n-1
        if(j/=i) work(j,:)=work(j,:)-work(j,i)*work(i,:)
      end do
    end do
    inverse=work(:,n:2*n-1)
    if(.not.all(ieee_is_finite(inverse))) error stop 'Nonfinite WENO inverse moments'
  end subroutine

  pure real(dp) function reconstruct_nonuniform_left(value,geometry,epsilon) result(face)
    real(dp), intent(in) :: value(5)
    type(weno_face_geometry), intent(in) :: geometry
    real(dp), optional, intent(in) :: epsilon
    real(dp) :: v(0:4),coeff(0:2),candidate(0:2),beta(0:2),alpha(0:2),logalpha(0:2)
    real(dp) :: scale,eps,tau,den,hi,lo
    integer :: r
    v=value-value(3);scale=maxval(abs(v))
    if(scale==0) then
      face=value(3);return
    end if
    v=v/scale;eps=1.e-20_dp
    if(present(epsilon)) eps=max(epsilon,tiny(eps))
    eps=exp(max(log(tiny(eps)),min(log(huge(eps))-2,log(eps)-2*log(scale))))
    do r=0,2
      coeff=matmul(geometry%polynomial(:,:,r),v(r:r+2))
      candidate(r)=coeff(0)
      ! Integral over the left adjacent cell x in [-1,0].
      beta(r)=(coeff(1)-coeff(2))**2+13._dp/3*coeff(2)**2
    end do
    tau=abs(beta(0)-beta(2))
    do r=0,2
      den=beta(r)+eps;hi=max(den,tau);lo=min(den,tau)
      logalpha(r)=log(geometry%optimal(r))+2*(log(hi)-log(den))+log(1+(lo/hi)**2)
    end do
    alpha=exp(logalpha-maxval(logalpha))
    face=value(3)+scale*dot_product(alpha,candidate)/sum(alpha)
  end function
end module
