module mod_mc_geometry
  use mod_precision, only: dp
  use mod_mc_euler_config, only: mc_euler_config
  implicit none
  private
  public :: mc_cell_volume, mc_cell_center, mc_face_area, mc_spacing, mc_metric_inverse
contains
  pure function mc_spacing(c) result(h)
    type(mc_euler_config), intent(in) :: c
    real(dp) :: h(3)
    h=[c%x_max-c%x_min,c%y_max-c%y_min,c%z_max-c%z_min]/real([c%nx,c%ny,c%nz],dp)
  end function
  pure real(dp) function height(c,x) result(h)
    type(mc_euler_config), intent(in) :: c
    real(dp), intent(in) :: x
    h=1.0_dp
    if(c%geometry /= 'planar_nozzle') return
    if(x <= c%nozzle_throat_x) then
      h=c%nozzle_throat_half_height+(c%nozzle_inlet_half_height-c%nozzle_throat_half_height)* &
        ((x-c%nozzle_throat_x)/(c%x_min-c%nozzle_throat_x))**2
    else
      h=c%nozzle_throat_half_height+(c%nozzle_exit_half_height-c%nozzle_throat_half_height)* &
        ((x-c%nozzle_throat_x)/(c%x_max-c%nozzle_throat_x))**2
    end if
  end function
  pure subroutine cell_heights(c,i,hl,hr)
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::i
    real(dp),intent(out)::hl,hr
    real(dp)::h(3),x
    h=mc_spacing(c)
    x=c%x_min+(i-1)*h(1)
    hl=height(c,x)
    hr=height(c,x+h(1))
  end subroutine
  pure real(dp) function mc_cell_volume(c,i,j,k) result(v)
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::i,j,k
    real(dp)::h(3),hl,hr
    h=mc_spacing(c)
    call cell_heights(c,i,hl,hr)
    v=product(h)*0.5_dp*(hl+hr)
  end function
  pure function mc_cell_center(c,i,j,k) result(x)
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::i,j,k
    real(dp)::x(3),h(3),hl,hr
    h=mc_spacing(c)
    call cell_heights(c,i,hl,hr)
    x=[c%x_min,c%y_min,c%z_min]+(real([i,j,k],dp)-0.5_dp)*h
    x(2)=x(2)*0.5_dp*(hl+hr)
  end function
  ! Positive-coordinate face area vector. p(d)=0..N identifies a face;
  ! transverse p values are cell indices. Shared faces have identical metrics.
  pure function mc_face_area(c,p,d) result(a)
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::p(3),d
    real(dp)::a(3),h(3),hl,hr,eta,x
    h=mc_spacing(c)
    a=0
    select case(d)
    case(1)
      x=c%x_min+p(1)*h(1)
      a(1)=height(c,x)*h(2)*h(3)
    case(2)
      call cell_heights(c,p(1),hl,hr)
      eta=c%y_min+p(2)*h(2)
      a=[-eta*(hr-hl)*h(3),h(1)*h(3),0.0_dp]
    case(3)
      call cell_heights(c,p(1),hl,hr)
      a(3)=h(1)*h(2)*0.5_dp*(hl+hr)
    end select
  end function
  ! Physical gradients = computational gradients * inverse map Jacobian.
  pure function mc_metric_inverse(c,i,j,k) result(m)
    type(mc_euler_config),intent(in)::c
    integer,intent(in)::i,j,k
    real(dp)::m(3,3),h(3),hl,hr,eta,hbar
    h=mc_spacing(c)
    call cell_heights(c,i,hl,hr)
    hbar=0.5_dp*(hl+hr)
    eta=c%y_min+(j-0.5_dp)*h(2)
    m=0
    m(1,1)=1
    m(2,1)=-eta*(hr-hl)/(h(1)*hbar)
    m(2,2)=1/hbar
    m(3,3)=1
  end function
end module mod_mc_geometry

