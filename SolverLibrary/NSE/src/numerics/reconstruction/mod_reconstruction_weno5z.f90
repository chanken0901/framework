module mod_reconstruction_weno5z
  use mod_precision, only : dp
  implicit none
  private

  real(dp), parameter, public :: weno5z_default_epsilon = 1.0e-20_dp

  public :: reconstruct_weno5z_left
  public :: reconstruct_weno5z_right

contains

  pure real(dp) function reconstruct_weno5z_left(value, epsilon) result(face)
    real(dp), intent(in) :: value(5)
    real(dp), intent(in), optional :: epsilon
    real(dp) :: candidate(3), beta(3), alpha(3), weight_sum
    real(dp) :: tau5, eps

    eps = weno5z_default_epsilon
    if (present(epsilon)) eps = max(epsilon, tiny(1.0_dp))

    candidate(1) = (2.0_dp*value(1) - 7.0_dp*value(2) + &
      11.0_dp*value(3)) / 6.0_dp
    candidate(2) = (-value(2) + 5.0_dp*value(3) + &
      2.0_dp*value(4)) / 6.0_dp
    candidate(3) = (2.0_dp*value(3) + 5.0_dp*value(4) - &
      value(5)) / 6.0_dp

    call smoothness_indicators(value, beta)
    tau5 = abs(beta(1) - beta(3))
    alpha(1) = 0.1_dp * (1.0_dp + (tau5/(beta(1)+eps))**2)
    alpha(2) = 0.6_dp * (1.0_dp + (tau5/(beta(2)+eps))**2)
    alpha(3) = 0.3_dp * (1.0_dp + (tau5/(beta(3)+eps))**2)
    weight_sum = sum(alpha)
    face = dot_product(alpha, candidate) / weight_sum
  end function reconstruct_weno5z_left

  pure real(dp) function reconstruct_weno5z_right(value, epsilon) result(face)
    real(dp), intent(in) :: value(5)
    real(dp), intent(in), optional :: epsilon
    real(dp) :: candidate(3), beta(3), alpha(3), weight_sum
    real(dp) :: tau5, eps

    eps = weno5z_default_epsilon
    if (present(epsilon)) eps = max(epsilon, tiny(1.0_dp))

    candidate(1) = (-value(1) + 5.0_dp*value(2) + &
      2.0_dp*value(3)) / 6.0_dp
    candidate(2) = (2.0_dp*value(2) + 5.0_dp*value(3) - &
      value(4)) / 6.0_dp
    candidate(3) = (11.0_dp*value(3) - 7.0_dp*value(4) + &
      2.0_dp*value(5)) / 6.0_dp

    call smoothness_indicators(value, beta)
    tau5 = abs(beta(1) - beta(3))

    ! The candidates above are ordered from the leftmost to the rightmost
    ! stencil. Their right-biased optimal weights are the mirror image of
    ! the left reconstruction weights.
    alpha(1) = 0.3_dp * (1.0_dp + (tau5/(beta(1)+eps))**2)
    alpha(2) = 0.6_dp * (1.0_dp + (tau5/(beta(2)+eps))**2)
    alpha(3) = 0.1_dp * (1.0_dp + (tau5/(beta(3)+eps))**2)
    weight_sum = sum(alpha)
    face = dot_product(alpha, candidate) / weight_sum
  end function reconstruct_weno5z_right

  pure subroutine smoothness_indicators(value, beta)
    real(dp), intent(in) :: value(5)
    real(dp), intent(out) :: beta(3)

    beta(1) = (13.0_dp/12.0_dp) * &
      (value(1)-2.0_dp*value(2)+value(3))**2 + 0.25_dp * &
      (value(1)-4.0_dp*value(2)+3.0_dp*value(3))**2
    beta(2) = (13.0_dp/12.0_dp) * &
      (value(2)-2.0_dp*value(3)+value(4))**2 + 0.25_dp * &
      (value(2)-value(4))**2
    beta(3) = (13.0_dp/12.0_dp) * &
      (value(3)-2.0_dp*value(4)+value(5))**2 + 0.25_dp * &
      (3.0_dp*value(3)-4.0_dp*value(4)+value(5))**2
  end subroutine smoothness_indicators

end module mod_reconstruction_weno5z
