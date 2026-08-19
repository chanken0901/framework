module mod_hit_isotropy_math
  use mod_precision, only : dp
  implicit none
  private

  public :: isotropy_error_3x3
  public :: symmetric_inverse_sqrt_3x3

contains

  pure real(dp) function isotropy_error_3x3(reynolds) result(error_value)
    real(dp), intent(in) :: reynolds(3,3)
    real(dp) :: trace_value
    integer :: i, j

    trace_value = reynolds(1,1) + reynolds(2,2) + reynolds(3,3)
    if (trace_value <= tiny(1.0_dp)) then
      error_value = 0.0_dp
      return
    end if

    error_value = 0.0_dp
    do i = 1, 3
      error_value = max(error_value, &
        abs(3.0_dp*reynolds(i,i)/trace_value - 1.0_dp))
      do j = i + 1, 3
        error_value = max(error_value, &
          abs(3.0_dp*reynolds(i,j)/trace_value))
      end do
    end do
  end function isotropy_error_3x3

  subroutine symmetric_inverse_sqrt_3x3(matrix, inverse_sqrt, success)
    real(dp), intent(in) :: matrix(3,3)
    real(dp), intent(out) :: inverse_sqrt(3,3)
    logical, intent(out) :: success
    real(dp) :: a(3,3), eigenvectors(3,3), eigenvalues(3)
    real(dp) :: app, aqq, apq, arp, arq, vrp, vrq
    real(dp) :: tau, tangent, cosine, sine, scale, floor_value
    real(dp) :: largest_off_diagonal
    integer :: i, j, p, q, sweep

    a = 0.5_dp * (matrix + transpose(matrix))
    eigenvectors = 0.0_dp
    do i = 1, 3
      eigenvectors(i,i) = 1.0_dp
    end do

    do sweep = 1, 32
      p = 1
      q = 2
      largest_off_diagonal = abs(a(1,2))
      if (abs(a(1,3)) > largest_off_diagonal) then
        p = 1
        q = 3
        largest_off_diagonal = abs(a(1,3))
      end if
      if (abs(a(2,3)) > largest_off_diagonal) then
        p = 2
        q = 3
        largest_off_diagonal = abs(a(2,3))
      end if

      scale = max(1.0_dp, maxval(abs(a)))
      if (largest_off_diagonal <= 64.0_dp*epsilon(1.0_dp)*scale) exit

      app = a(p,p)
      aqq = a(q,q)
      apq = a(p,q)
      tau = (aqq-app) / (2.0_dp*apq)
      if (tau >= 0.0_dp) then
        tangent = 1.0_dp / (tau + sqrt(1.0_dp+tau*tau))
      else
        tangent = -1.0_dp / (-tau + sqrt(1.0_dp+tau*tau))
      end if
      cosine = 1.0_dp / sqrt(1.0_dp+tangent*tangent)
      sine = tangent * cosine

      do i = 1, 3
        if (i == p .or. i == q) cycle
        arp = a(i,p)
        arq = a(i,q)
        a(i,p) = cosine*arp - sine*arq
        a(p,i) = a(i,p)
        a(i,q) = sine*arp + cosine*arq
        a(q,i) = a(i,q)
      end do
      a(p,p) = cosine*cosine*app - 2.0_dp*sine*cosine*apq + &
        sine*sine*aqq
      a(q,q) = sine*sine*app + 2.0_dp*sine*cosine*apq + &
        cosine*cosine*aqq
      a(p,q) = 0.0_dp
      a(q,p) = 0.0_dp

      do i = 1, 3
        vrp = eigenvectors(i,p)
        vrq = eigenvectors(i,q)
        eigenvectors(i,p) = cosine*vrp - sine*vrq
        eigenvectors(i,q) = sine*vrp + cosine*vrq
      end do
    end do

    eigenvalues = [a(1,1), a(2,2), a(3,3)]
    scale = sum(eigenvalues) / 3.0_dp
    success = scale > tiny(1.0_dp)
    if (.not. success) then
      inverse_sqrt = 0.0_dp
      return
    end if

    if (minval(eigenvalues) < -sqrt(epsilon(1.0_dp))*scale) then
      success = .false.
      inverse_sqrt = 0.0_dp
      return
    end if

    floor_value = max(scale*1.0e-12_dp, tiny(1.0_dp))
    inverse_sqrt = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        inverse_sqrt(i,j) = sum(eigenvectors(i,:) * eigenvectors(j,:) / &
          sqrt(max(eigenvalues, floor_value)))
      end do
    end do
  end subroutine symmetric_inverse_sqrt_3x3

end module mod_hit_isotropy_math
