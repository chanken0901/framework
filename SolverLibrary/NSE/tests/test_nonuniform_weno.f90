program test_nonuniform_weno
  use mod_precision, only: dp
  use mod_reconstruction_nonuniform
  use mod_reconstruction_weno5z, only: reconstruct_weno5z_left
  implicit none
  type(weno_face_geometry) :: geometry
  real(dp) :: edges(0:5),values(5),actual,expected,combined(0:4),h
  integer :: i,r,p
  edges=[-3._dp,-2._dp,-1._dp,0._dp,1._dp,2._dp]
  call build_weno_geometry(edges,geometry)
  if(maxval(abs(geometry%optimal-[.1_dp,.6_dp,.3_dp]))>1.e-12_dp) error stop 'Uniform optimal weights'
  values=[1._dp,.8_dp,1.2_dp,2._dp,1.5_dp]
  actual=reconstruct_nonuniform_left(values,geometry)
  expected=reconstruct_weno5z_left(values)
  if(abs(actual-expected)>1.e-12_dp) error stop 'Uniform WENO-Z equivalence'
  do i=1,3
    edges=[-4._dp,-2.5_dp,-1._dp,0._dp,1.2_dp,3._dp]*real(i,dp)
    call build_weno_geometry(edges,geometry)
    combined=0
    do r=0,2
      combined(r:r+2)=combined(r:r+2)+geometry%optimal(r)*geometry%polynomial(0,:,r)
    end do
    h=edges(3)-edges(2)
    do p=0,4
      do r=1,5
        values(r)=((edges(r)/h)**(p+1)-(edges(r-1)/h)**(p+1))/ &
          ((p+1)*(edges(r)-edges(r-1))/h)
      end do
      expected=0
      if(p==0) expected=1
      if(abs(dot_product(combined,values)-expected)>1.e-11_dp) error stop 'Cell-average polynomial exactness'
    end do
    values=3.7_dp
    if(reconstruct_nonuniform_left(values,geometry)/=3.7_dp) error stop 'Constant preservation'
  end do
  print *, '[OK] nonuniform WENO geometry, uniform limit and polynomial reconstruction'
end program
