module mod_nse_field
  use mod_precision,     only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config,  only : nse_config
  implicit none

  private

  public :: allocate_nse_fields
  public :: deallocate_nse_fields
  public :: Q, Q0, RHS, F, Qw

  real(dp), allocatable :: Q(:,:,:,:)
  real(dp), allocatable :: Q0(:,:,:,:)
  real(dp), allocatable :: RHS(:,:,:,:)
  real(dp), allocatable :: F(:,:,:,:)
  real(dp), allocatable :: Qw(:,:,:,:)

  real(dp), allocatable :: QL (:,:,:,:), QR(:,:,:,:) ! reconstructed face states
  real(dp), allocatable :: Q_vis (:,:,:,:)              ! [1-sim%nghost:sim%nx+sim%nghost, 1-sim%nghost:sim%ny+sim%nghost, 1-sim%nghost:sim%nz+sim%nghost, nse%nv]


contains

  subroutine allocate_nse_fields(sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config),        intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke

    allocate(Q  (1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))
    allocate(Q0 (1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))
    allocate(RHS(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))

    allocate(F  (0:sim%nx, js-1:je, ks-1:ke, nse%nv))
    allocate(Qw (1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))

  end subroutine allocate_nse_fields


  subroutine deallocate_nse_fields()
    if (allocated(Q))   deallocate(Q)
    if (allocated(Q0))  deallocate(Q0)
    if (allocated(RHS)) deallocate(RHS)
    if (allocated(F))   deallocate(F)
    if (allocated(Qw))  deallocate(Qw)
  end subroutine deallocate_nse_fields

end module mod_nse_field