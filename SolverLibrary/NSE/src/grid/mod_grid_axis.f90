module mod_grid_axis
  ! Backend-neutral axis geometry. Not yet wired to production NSE operators.
  use mod_precision, only: dp
  use mod_reconstruction_nonuniform, only: weno_face_geometry,build_weno_geometry
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: grid_axis,build_axis,build_sinh_axis,axis_derivative_weights
  public :: prepare_axis_weno
  type :: grid_axis
    integer :: n=0,ng=0
    logical :: periodic=.false.
    real(dp) :: minimum_width=0
    real(dp), allocatable :: edge(:),center(:),width(:)
    ! Offset is contiguous: directly packable for GPU, global cell indexing.
    real(dp), allocatable :: d1(:,:),d2(:,:)
    type(weno_face_geometry), allocatable :: weno_left(:),weno_right(:)
  end type
contains
  subroutine prepare_axis_weno(axis)
    type(grid_axis), intent(inout) :: axis
    integer :: f
    if(allocated(axis%weno_left)) deallocate(axis%weno_left,axis%weno_right)
    allocate(axis%weno_left(0:axis%n),axis%weno_right(0:axis%n))
    do f=0,axis%n
      call build_weno_geometry(axis%edge(f-3:f+2),axis%weno_left(f))
      call build_weno_geometry(-axis%edge(f+3:f-2:-1),axis%weno_right(f))
    end do
  end subroutine
  subroutine check(valid,message)
    logical, intent(in) :: valid
    character(*), intent(in) :: message
    if(.not.valid) then
      write(*,'(a)') 'ERROR: '//message
      error stop 'Invalid grid axis'
    end if
  end subroutine

  subroutine build_sinh_axis(n,ng,lower,upper,strength,periodic,axis)
    integer, intent(in) :: n,ng
    real(dp), intent(in) :: lower,upper,strength
    logical, intent(in) :: periodic
    type(grid_axis), intent(out) :: axis
    real(dp), allocatable :: edges(:)
    real(dp) :: s,fraction
    integer :: i
    call check(n>=3.and.ng>=3.and.ng<=n,'Axis requires n>=ng>=3')
    call check(all(ieee_is_finite([lower,upper,strength])),'Nonfinite axis bounds/strength')
    call check(upper>lower.and.strength>=0.and.strength<=20,'Invalid axis bounds/strength')
    allocate(edges(0:n))
    do i=0,n
      s=real(i,dp)/n
      if(strength<sqrt(epsilon(strength))) then
        fraction=s
      else
        fraction=.5_dp*(1+sinh(strength*(2*s-1))/sinh(strength))
      end if
      edges(i)=lower+(upper-lower)*fraction
    end do
    edges(0)=lower;edges(n)=upper
    call build_axis(edges,ng,periodic,axis)
  end subroutine

  subroutine build_axis(edges,ng,periodic,axis)
    real(dp), intent(in) :: edges(0:)
    integer, intent(in) :: ng
    logical, intent(in) :: periodic
    type(grid_axis), intent(out) :: axis
    real(dp) :: weights(0:6,0:2)
    integer :: n,i
    n=size(edges)-1
    call check(n>=3.and.ng>=3.and.ng<=n,'Axis requires n>=ng>=3')
    call check(all(ieee_is_finite(edges)),'Nonfinite edge coordinates')
    call check(all(edges(1:n)>edges(0:n-1)),'Edges must strictly increase')
    axis%n=n;axis%ng=ng;axis%periodic=periodic
    allocate(axis%edge(-ng:n+ng),axis%center(1-ng:n+ng),axis%width(1-ng:n+ng))
    allocate(axis%d1(-3:3,n),axis%d2(-3:3,n))
    axis%edge(0:n)=edges
    do i=1,ng
      if(periodic) then
        axis%edge(-i)=axis%edge(1-i)-(edges(n-i+1)-edges(n-i))
        axis%edge(n+i)=axis%edge(n+i-1)+(edges(i)-edges(i-1))
      else
        axis%edge(-i)=axis%edge(1-i)-(edges(i)-edges(i-1))
        axis%edge(n+i)=axis%edge(n+i-1)+(edges(n-i+1)-edges(n-i))
      end if
    end do
    call check(all(ieee_is_finite(axis%edge)),'Ghost coordinate overflow')
    do i=1-ng,n+ng
      axis%width(i)=axis%edge(i)-axis%edge(i-1)
      axis%center(i)=axis%edge(i-1)+axis%width(i)/2
    end do
    call check(all(ieee_is_finite(axis%width)).and.all(axis%width>0),'Invalid axis widths')
    call check(all(ieee_is_finite(axis%center)),'Nonfinite centers')
    axis%minimum_width=minval(axis%width(1:n))
    do i=1,n
      call axis_derivative_weights(axis%center(i-3:i+3),axis%center(i),weights)
      axis%d1(:,i)=weights(:,1);axis%d2(:,i)=weights(:,2)
    end do
  end subroutine

  subroutine axis_derivative_weights(points,target,weights)
    ! Polynomial differentiation on seven point values, with scaled coordinates.
    ! NOT coefficients for finite-volume reconstruction or conservative fluxes.
    real(dp), intent(in) :: points(0:6),target
    real(dp), intent(out) :: weights(0:6,0:2)
    real(dp) :: x(0:6),scale,c1,c2,c3,c4,c5
    integer :: i,j,k,mn
    call check(all(ieee_is_finite(points)).and.ieee_is_finite(target),'Nonfinite stencil')
    call check(all(points(1:)>points(:5)),'Stencil points must strictly increase')
    scale=maxval(abs(points-target))
    call check(ieee_is_finite(scale).and.scale>0,'Invalid stencil scale')
    x=(points-target)/scale
    weights=0;weights(0,0)=1;c1=1;c4=x(0)
    do i=1,6
      mn=min(i,2);c2=1;c5=c4;c4=x(i)
      do j=0,i-1
        c3=x(i)-x(j);c2=c2*c3
        call check(c3>0.and.c2>tiny(c2),'Degenerate derivative stencil')
        if(j==i-1) then
          do k=mn,1,-1
            weights(i,k)=c1*(k*weights(i-1,k-1)-c5*weights(i-1,k))/c2
          end do
          weights(i,0)=-c1*c5*weights(i-1,0)/c2
        end if
        do k=mn,1,-1
          weights(j,k)=(c4*weights(j,k)-k*weights(j,k-1))/c3
        end do
        weights(j,0)=c4*weights(j,0)/c3
      end do
      c1=c2
    end do
    weights(:,1)=weights(:,1)/scale
    weights(:,2)=(weights(:,2)/scale)/scale
    call check(all(ieee_is_finite(weights)),'Derivative coefficient overflow')
  end subroutine
end module
