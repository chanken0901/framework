module mod_constants
    use mod_precision
    implicit none
    private

    public :: pi, twopi, halfpi
    public :: c1_2, c1_3, c1_4 

    real(dp), parameter :: pi     = acos(-1.0_dp)
    real(dp), parameter :: twopi  = 2.0_dp*pi
    real(dp), parameter :: halfpi = 0.5_dp*pi

    real(dp), parameter :: c1_2 = 1.0_dp/2.0_dp
    real(dp), parameter :: c1_3 = 1.0_dp/3.0_dp
    real(dp), parameter :: c1_4 = 1.0_dp/4.0_dp

end module mod_constants