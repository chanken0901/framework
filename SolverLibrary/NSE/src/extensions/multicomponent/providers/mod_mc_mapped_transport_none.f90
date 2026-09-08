module mod_mc_mapped_transport_provider
  use mod_precision, only: dp
  use mod_mc_state_layout, only: mc_state_layout
  implicit none
  private
  public :: transport_flux
contains
  subroutine transport_flux(v,rho,grad,l,normal,flux)
    real(dp),intent(in)::v(:),rho,grad(:,:),normal(3)
    type(mc_state_layout),intent(in)::l
    real(dp),intent(out)::flux(:)
    error stop 'mapped transport requires mixture_averaged transport provider'
  end subroutine
end module

