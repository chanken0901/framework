program test_hit_isotropy_math
  use mod_precision, only : dp
  use mod_hit_isotropy_math, only : isotropy_error_3x3, &
    symmetric_inverse_sqrt_3x3, velocity_gradient_contraction_3x3, &
    pressure_poisson_multiplier
  implicit none

  real(dp) :: reynolds(3,3), inverse_sqrt(3,3), transformed(3,3)
  real(dp) :: rotation(3,3), diagonal(3,3), target
  real(dp) :: gradient(3,3), x, y, expected_contraction
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

  ! Taylor-Green velocity u=sin(x)cos(y), v=-cos(x)sin(y) has
  ! (du_i/dx_j)(du_j/dx_i)=cos(2x)+cos(2y).  This also fixes the
  ! sign convention used by the spectral pressure Poisson solve.
  x = 0.37_dp
  y = 0.91_dp
  gradient = 0.0_dp
  gradient(1,1) = cos(x)*cos(y)
  gradient(1,2) = -sin(x)*sin(y)
  gradient(2,1) = sin(x)*sin(y)
  gradient(2,2) = -cos(x)*cos(y)
  expected_contraction = cos(2.0_dp*x) + cos(2.0_dp*y)
  if (abs(velocity_gradient_contraction_3x3(gradient) - &
      expected_contraction) > 1.0e-12_dp) then
    error stop 'HIT pressure source contraction has the wrong sign'
  end if
  if (abs(pressure_poisson_multiplier(4.0_dp,2.0_dp)-0.5_dp) > &
      1.0e-15_dp) then
    error stop 'HIT pressure Poisson multiplier is incorrect'
  end if
  if (abs(pressure_poisson_multiplier(0.0_dp,2.0_dp)) > tiny(1.0_dp)) then
    error stop 'HIT pressure zero mode must not be divided by k squared'
  end if

  write(*,'(A)') 'HIT isotropy tensor test passed'
end program test_hit_isotropy_math
