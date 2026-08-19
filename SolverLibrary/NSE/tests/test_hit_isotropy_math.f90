program test_hit_isotropy_math
  use mod_precision, only : dp
  use mod_hit_isotropy_math, only : isotropy_error_3x3, &
    symmetric_inverse_sqrt_3x3
  implicit none

  real(dp) :: reynolds(3,3), inverse_sqrt(3,3), transformed(3,3)
  real(dp) :: rotation(3,3), diagonal(3,3), target
  logical :: success

  rotation = reshape([ &
    0.8_dp, -0.6_dp, 0.0_dp, &
    0.6_dp,  0.8_dp, 0.0_dp, &
    0.0_dp,  0.0_dp, 1.0_dp ], [3,3])
  diagonal = 0.0_dp
  diagonal(1,1) = 1.0_dp
  diagonal(2,2) = 2.0_dp
  diagonal(3,3) = 4.0_dp
  reynolds = matmul(rotation, matmul(diagonal, transpose(rotation)))

  call symmetric_inverse_sqrt_3x3(reynolds, inverse_sqrt, success)
  if (.not. success) error stop 'inverse square root unexpectedly failed'
  target = (reynolds(1,1)+reynolds(2,2)+reynolds(3,3)) / 3.0_dp
  transformed = target * matmul(inverse_sqrt, &
    matmul(reynolds, transpose(inverse_sqrt)))

  if (isotropy_error_3x3(transformed) > 1.0e-12_dp) then
    error stop 'inverse square root did not isotropize the test tensor'
  end if
  if (abs(transformed(1,1)-target) > 1.0e-12_dp .or. &
      abs(transformed(2,2)-target) > 1.0e-12_dp .or. &
      abs(transformed(3,3)-target) > 1.0e-12_dp) then
    error stop 'isotropic tensor energy was not preserved'
  end if

  write(*,'(A)') 'HIT isotropy tensor test passed'
end program test_hit_isotropy_math
