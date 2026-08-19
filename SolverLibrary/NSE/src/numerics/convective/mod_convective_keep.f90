module mod_convective_keep
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  implicit none
  private

  ! Positive-side coefficients of the sixth-order central first derivative.
  real(dp), parameter :: central6_coefficient(3) = [ &
    3.0_dp/4.0_dp, -3.0_dp/20.0_dp, 1.0_dp/60.0_dp ]

  public :: compute_keep_flux
  public :: compute_keep_face_flux
  public :: validate_keep_scheme
  public :: keep_required_ghost_cells

contains

  subroutine compute_keep_flux(q, fface, direction, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: direction, js, je, ks, ke
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(out) :: fface(0:, js-1:, ks-1:, :)
    integer :: i, j, k
    integer :: selected_order

    if (nse%nv /= 5) error stop 'KEEP flux requires five conserved variables'
    if (sim%nghost < keep_required_ghost_cells()) then
      error stop 'selectable-order KEEP flux requires three ghost cells'
    end if
    selected_order = requested_keep_order(nse)
    if (selected_order /= 2 .and. selected_order /= 6) then
      error stop 'KEEP scheme must be keep2 or keep6'
    end if
    ! A conservative face flux whose adjacent difference is
    !
    !   2 sum_s d_s [f#(q_i,q_{i+s}) - f#(q_{i-s},q_i)].
    !
    ! Here f# is the symmetric two-point KEEP flux.  The selected positive
    ! derivative coefficients are d=[1/2] for second order or
    ! d=[3/4,-3/20,1/60] for sixth order.  The telescoping face construction
    ! retains the conservative KEEP pairwise split form for either order.
    select case (direction)
    case (1)
      !$OMP DO collapse(2) schedule(static)
      do k = ks-1, ke
        do j = js-1, je
          do i = 0, sim%nx
            call compute_keep_face_flux(q, i, j, k, direction, &
              selected_order, sim, nse, js, ks, fface(i,j,k,1:5))
          end do
        end do
      end do
      !$OMP END DO

    case (2)
      !$OMP DO collapse(2) schedule(static)
      do k = ks-1, ke
        do j = js-1, je
          do i = 0, sim%nx
            call compute_keep_face_flux(q, i, j, k, direction, &
              selected_order, sim, nse, js, ks, fface(i,j,k,1:5))
          end do
        end do
      end do
      !$OMP END DO

    case (3)
      !$OMP DO collapse(2) schedule(static)
      do k = ks-1, ke
        do j = js-1, je
          do i = 0, sim%nx
            call compute_keep_face_flux(q, i, j, k, direction, &
              selected_order, sim, nse, js, ks, fface(i,j,k,1:5))
          end do
        end do
      end do
      !$OMP END DO

    case default
      error stop 'convective flux direction must be 1, 2, or 3'
    end select
  end subroutine compute_keep_flux

  pure subroutine compute_keep_face_flux(q, i, j, k, direction, order, &
      sim, nse, js, ks, flux)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: i, j, k, direction, order, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, &
      ks-sim%nghost:, :)
    real(dp), intent(out) :: flux(5)
    real(dp) :: pair_flux(5), weight, derivative_coefficient(3)
    integer :: separation, offset, maximum_separation
    integer :: im, jm, km, ip, jp, kp

    call select_derivative_coefficients(order, derivative_coefficient, &
      maximum_separation)
    flux = 0.0_dp
    do separation = 1, maximum_separation
      weight = 2.0_dp * derivative_coefficient(separation)
      do offset = 0, separation-1
        im = i
        jm = j
        km = k
        ip = i
        jp = j
        kp = k
        select case (direction)
        case (1)
          im = i-offset
          ip = i-offset+separation
        case (2)
          jm = j-offset
          jp = j-offset+separation
        case (3)
          km = k-offset
          kp = k-offset+separation
        case default
          return
        end select
        call keep_two_point_flux(q, im, jm, km, ip, jp, kp, direction, &
          sim, nse, js, ks, pair_flux)
        flux = flux + weight*pair_flux
      end do
    end do
  end subroutine compute_keep_face_flux

  pure subroutine select_derivative_coefficients(order, coefficient, &
      maximum_separation)
    integer, intent(in) :: order
    real(dp), intent(out) :: coefficient(3)
    integer, intent(out) :: maximum_separation

    coefficient = 0.0_dp
    select case (order)
    case (2)
      coefficient(1) = 0.5_dp
      maximum_separation = 1
    case (6)
      coefficient = central6_coefficient
      maximum_separation = 3
    case default
      coefficient = 0.0_dp
      maximum_separation = 0
    end select
  end subroutine select_derivative_coefficients

  pure subroutine keep_two_point_flux(q, im, jm, km, ip, jp, kp, &
      direction, sim, nse, js, ks, flux)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: im, jm, km, ip, jp, kp, direction, js, ks
    real(dp), intent(in) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), intent(out) :: flux(5)
    real(dp) :: rm, um, vm, wm, pm
    real(dp) :: rp, up, vp, wp, pp
    real(dp) :: rmp, ump, vmp, wmp, lmp
    real(dp) :: normal_velocity, minus_normal, plus_normal
    real(dp) :: ck, mxk, myk, mzk, kk, lk, gk, pk

    rm = max(q(im,jm,km,1), nse%small_rho)
    um = q(im,jm,km,2) / rm
    vm = q(im,jm,km,3) / rm
    wm = q(im,jm,km,4) / rm
    pm = max((nse%gamma-1.0_dp) * &
      (q(im,jm,km,5)-0.5_dp*rm*(um*um+vm*vm+wm*wm)), nse%small_p)

    rp = max(q(ip,jp,kp,1), nse%small_rho)
    up = q(ip,jp,kp,2) / rp
    vp = q(ip,jp,kp,3) / rp
    wp = q(ip,jp,kp,4) / rp
    pp = max((nse%gamma-1.0_dp) * &
      (q(ip,jp,kp,5)-0.5_dp*rp*(up*up+vp*vp+wp*wp)), nse%small_p)

    rmp = 0.5_dp * (rm + rp)
    ump = 0.5_dp * (um + up)
    vmp = 0.5_dp * (vm + vp)
    wmp = 0.5_dp * (wm + wp)
    lmp = 0.5_dp * (pm/rm + pp/rp) / (nse%gamma-1.0_dp)

    select case (direction)
    case (1)
      normal_velocity = ump
      minus_normal = um
      plus_normal = up
    case (2)
      normal_velocity = vmp
      minus_normal = vm
      plus_normal = vp
    case default
      normal_velocity = wmp
      minus_normal = wm
      plus_normal = wp
    end select

    ck = rmp * normal_velocity
    mxk = ck * ump
    myk = ck * vmp
    mzk = ck * wmp
    kk = ck * 0.5_dp * (um*up + vm*vp + wm*wp)
    lk = ck * lmp
    gk = 0.5_dp * (pm + pp)
    pk = 0.5_dp * (plus_normal*pm + minus_normal*pp)

    flux(1) = ck
    flux(2) = mxk
    flux(3) = myk
    flux(4) = mzk
    flux(1+direction) = flux(1+direction) + gk
    flux(5) = kk + lk + pk
  end subroutine keep_two_point_flux

  subroutine validate_keep_scheme(nse)
    type(nse_config), intent(in) :: nse
    integer :: selected_order

    selected_order = requested_keep_order(nse)
    if (selected_order /= 2 .and. selected_order /= 6) then
      write(*,'(A,A,A)') 'ERROR: unsupported KEEP scheme "', &
        trim(adjustl(nse%convective_scheme)), '"; use keep2 or keep6'
      error stop
    end if
  end subroutine validate_keep_scheme

  pure integer function requested_keep_order(nse) result(order)
    type(nse_config), intent(in) :: nse

    select case (trim(adjustl(nse%convective_scheme)))
    case ('keep2')
      order = 2
    case ('keep6')
      order = 6
    case default
      order = 0
    end select
  end function requested_keep_order

  integer function keep_required_ghost_cells() result(nghost)
    nghost = 3
  end function keep_required_ghost_cells

end module mod_convective_keep
