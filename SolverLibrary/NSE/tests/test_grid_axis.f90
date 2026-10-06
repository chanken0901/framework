program test_grid_axis
  use mod_precision, only: dp
  use mod_grid_axis
  implicit none
  type(grid_axis) :: axis,periodic
  real(dp) :: edges(0:12),w(0:6,0:2),x(0:6),target,value,exact,scale
  integer :: i,k,d
  character(32) :: mode
  if(command_argument_count()>0) then
    call get_command_argument(1,mode)
    edges=[(real(i,dp),i=0,12)];edges(4)=edges(3)
    call build_axis(edges,3,.false.,axis)
    stop
  end if
  call build_sinh_axis(12,3,-1._dp,1._dp,0._dp,.false.,axis)
  if(maxval(abs(axis%width(1:12)-1._dp/6))>1.e-14_dp) error stop 'Uniform widths'
  if(maxval(abs(axis%d1(:,6)*axis%width(6)- &
    [-1._dp/60,3._dp/20,-3._dp/4,0._dp,3._dp/4,-3._dp/20,1._dp/60]))>1.e-12_dp) &
    error stop 'Uniform first derivative coefficients'
  call build_sinh_axis(12,3,-1._dp,1._dp,2._dp,.false.,axis)
  if(axis%width(6)>=axis%width(1)) error stop 'Center refinement'
  do i=0,12
    edges(i)=(real(i,dp)/12)**1.5_dp
  end do
  call build_axis(edges,3,.true.,periodic)
  if(abs(periodic%width(0)-periodic%width(12))>1.e-14_dp) error stop 'Periodic left ghost'
  if(abs(periodic%width(13)-periodic%width(1))>1.e-14_dp) error stop 'Periodic right ghost'
  do i=1,12
    x=axis%center(i-3:i+3);target=axis%center(i)
    scale=maxval(abs(x-target));x=(x-target)/scale
    call axis_derivative_weights(x,0._dp,w)
    do d=1,2
      do k=0,6
        value=sum(w(:,d)*x**k);exact=0
        if(k==d) exact=real(d,dp)
        if(abs(value-exact)>1.e-10_dp) error stop 'Polynomial derivative exactness'
      end do
    end do
  end do
  print *, '[OK] nonuniform axis geometry and polynomial derivative coefficients'
end program
