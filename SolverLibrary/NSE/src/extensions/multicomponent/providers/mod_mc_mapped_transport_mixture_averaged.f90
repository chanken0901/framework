module mod_mc_mapped_transport_provider
  use mod_precision, only: dp
  use mod_mc_state_layout, only: mc_state_layout
  use mod_mc_thermodynamics_provider, only: mc_species_enthalpies, mc_mixture_cp
  use mod_mc_transport_provider, only: mc_dynamic_viscosity, mc_thermal_conductivity, mc_species_diffusivities
  implicit none
  private
  public :: transport_flux
contains
  subroutine transport_flux(v,rho,grad,l,normal,flux)
    real(dp),intent(in)::v(:),rho,grad(:,:),normal(3)
    type(mc_state_layout),intent(in)::l
    real(dp),intent(out)::flux(:)
    real(dp)::mu,kappa,cp,diff(l%nspecies),enthalpy(l%nspecies)
    real(dp)::js(l%nspecies),tau(3,3),stress(3),div
    integer::a
    mu=mc_dynamic_viscosity(v(5:),l,v(4))
    cp=mc_mixture_cp(v(5:),l,v(4))
    kappa=mc_thermal_conductivity(v(5:),l,v(4),cp)
    call mc_species_diffusivities(v(5:),l,v(4),diff)
    call mc_species_enthalpies(l,v(4),enthalpy)
    js=-rho*diff*matmul(grad(5:,:),normal)
    js=js-v(5:)*sum(js)
    js(l%nspecies)=js(l%nspecies)-sum(js)
    div=grad(1,1)+grad(2,2)+grad(3,3)
    tau=mu*(grad(1:3,:)+transpose(grad(1:3,:)))
    do a=1,3
      tau(a,a)=tau(a,a)-2*mu*div/3
    end do
    stress=matmul(tau,normal)
    flux=0
    flux(l%first_species:l%last_species)=-js
    flux(l%momentum)=stress
    flux(l%total_energy)=dot_product(stress,v(1:3))+kappa*dot_product(grad(4,:),normal)- &
      dot_product(enthalpy,js)
  end subroutine
end module

