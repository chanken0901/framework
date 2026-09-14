module mod_rf_units
  use mod_rf_mechanism, only: dp,require,ieee_is_finite
  implicit none
  type :: rf_reference_scales
    real(dp) :: density=1,velocity=1,length=1,temperature=1
  contains
    procedure :: scale=>reference_scale
    procedure :: to_si
    procedure :: from_si
  end type
contains
  real(dp) function reference_scale(self,quantity) result(v)
    class(rf_reference_scales), intent(in) :: self
    character(*), intent(in) :: quantity
    call require(all(ieee_is_finite([self%density,self%velocity,self%length,self%temperature])), 'Nonfinite scales')
    call require(min(self%density,self%velocity,self%length,self%temperature)>0,'Reference scales must be positive')
    select case(quantity)
    case('density'); v=self%density
    case('velocity'); v=self%velocity
    case('length'); v=self%length
    case('temperature'); v=self%temperature
    case('time'); v=self%length/self%velocity
    case('pressure','energy_density'); v=self%density*self%velocity**2
    case('energy'); v=self%velocity**2
    case default
      call require(.false.,'Unknown reference quantity: '//quantity)
      v=0
    end select
    call require(ieee_is_finite(v).and.v>0,'Invalid derived reference scale')
  end function
  real(dp) function to_si(self,value,quantity) result(v)
    class(rf_reference_scales), intent(in) :: self
    real(dp), intent(in) :: value
    character(*), intent(in) :: quantity
    call require(ieee_is_finite(value),'Nonfinite dimensionless value')
    v=value*self%scale(quantity)
    call require(ieee_is_finite(v),'SI conversion overflow')
  end function
  real(dp) function from_si(self,value,quantity) result(v)
    class(rf_reference_scales), intent(in) :: self
    real(dp), intent(in) :: value
    character(*), intent(in) :: quantity
    call require(ieee_is_finite(value),'Nonfinite SI value')
    v=value/self%scale(quantity)
    call require(ieee_is_finite(v),'Dimensionless conversion overflow')
  end function
end module
