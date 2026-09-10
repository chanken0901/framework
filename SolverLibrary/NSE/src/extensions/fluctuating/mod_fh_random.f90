! Philox4x32-10 (Random123), evaluated without signed integer overflow.
module mod_fh_random
  use iso_fortran_env, only: int64
  use mod_precision, only: dp
  implicit none
  private
  public :: fh_normal, philox4x32
  integer(int64), parameter :: mask32=4294967295_int64, base=4294967296_int64
contains
  pure subroutine multiply32(a,b,hi,lo)
    integer(int64), intent(in) :: a,b
    integer(int64), intent(out) :: hi,lo
    integer(int64) :: p0,p1,p2,p3,t
    p0=iand(a,65535_int64)*iand(b,65535_int64)
    p1=ishft(a,-16)*iand(b,65535_int64)
    p2=iand(a,65535_int64)*ishft(b,-16)
    p3=ishft(a,-16)*ishft(b,-16)
    t=p0+ishft(iand(p1,65535_int64)+iand(p2,65535_int64),16)
    lo=iand(t,mask32)
    hi=iand(p3+ishft(p1,-16)+ishft(p2,-16)+ishft(t,-32),mask32)
  end subroutine

  pure function philox4x32(counter,key) result(c)
    integer(int64), intent(in) :: counter(4),key(2)
    integer(int64) :: c(4),k(2),h0,h1,l0,l1
    integer :: r
    c=counter; k=key
    do r=1,10
      call multiply32(3528531795_int64,c(1),h0,l0)
      call multiply32(3449720151_int64,c(3),h1,l1)
      c=[ieor(ieor(h1,c(2)),k(1)),l1,ieor(ieor(h0,c(4)),k(2)),l0]
      k=iand(k+[2654435769_int64,3144134277_int64],mask32)
    end do
  end function

  pure real(dp) function fh_normal(i,j,k,step,seed,stream) result(z)
    integer, intent(in) :: i,j,k,step,seed,stream
    integer(int64) :: c(4)
    real(dp) :: u,v
    c=philox4x32(int([i,j,k,step],int64),int([seed,stream],int64))
    u=(real(c(1),dp)+0.5_dp)/real(base,dp)
    v=(real(c(2),dp)+0.5_dp)/real(base,dp)
    z=sqrt(-2.0_dp*log(u))*cos(2.0_dp*acos(-1.0_dp)*v)
  end function
end module
