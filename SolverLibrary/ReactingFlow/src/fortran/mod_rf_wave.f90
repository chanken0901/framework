module mod_rf_wave
  use mod_rf_thermo, only: dp,require,ieee_is_finite
  implicit none
  private
  public :: wave_history,sample_wave
  type :: wave_history
    real(dp) :: time=-huge(1._dp),position=0
    logical :: detected=.false.
  end type
contains
  subroutine sample_wave(pressure,dx,time,xmin,xmax,min_jump,history,x,speed,jump,detected,speed_valid)
    real(dp), intent(in) :: pressure(:),dx,time,xmin,xmax,min_jump
    type(wave_history), intent(inout) :: history
    real(dp), intent(out) :: x,speed,jump
    logical, intent(out) :: detected,speed_valid
    integer :: i,index
    real(dp) :: candidate
    call require(size(pressure)>=2,'Wave diagnostic needs at least two cells')
    call require(all(ieee_is_finite(pressure)).and.all(pressure>0),'Invalid wave diagnostic pressure')
    call require(all(ieee_is_finite([dx,time,xmin,xmax,min_jump])),'Nonfinite wave diagnostic control')
    call require(dx>0.and.xmin>=0.and.xmax>xmin.and.xmax<=size(pressure)*dx, &
                 'Invalid wave search interval')
    call require(min_jump>0,'Wave minimum pressure jump must be positive')
    call require(time>history%time,'Wave samples must have strictly increasing time')
    x=0;speed=0;jump=0;index=0
    do i=1,size(pressure)-1
      if(i*dx<xmin.or.i*dx>xmax) cycle
      candidate=abs(pressure(i+1)-pressure(i))
      ! Equal maxima deterministically select the leftmost face.
      if(candidate>jump) then
        jump=candidate;index=i
      end if
    end do
    detected=index>0.and.jump>=min_jump
    speed_valid=detected.and.history%detected
    if(detected) x=index*dx
    if(speed_valid) speed=(x-history%position)/(time-history%time)
    history%time=time;history%position=x;history%detected=detected
  end subroutine
end module
