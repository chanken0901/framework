program test_multicomponent_mpi_pencil
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config,read_mc_config
  use mod_mc_state_layout, only : mc_state_layout,initialize_mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config,read_mc_euler_config
  use mod_mc_reactive_config, only : mc_reactive_config,read_mc_reactive_config
  use mod_mc_euler_field, only : set_mc_euler_conservative_state
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_rhs
  use mod_mc_reactive_solver, only : advance_mc_reactive_strang
  use mod_mc_thermodynamics_provider, only : configure_mc_thermodynamics,mc_mixture_gas_constant
  use mod_mc_transport_provider, only : configure_mc_transport
  use mod_mc_chemistry_provider, only : configure_mc_chemistry
  use mod_mc_mpi_pencil
  implicit none
  include 'mpif.h'
  type(mc_config) :: model
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: config
  type(mc_reactive_config) :: reactive
  type(mc_pencil_domain) :: domain
  real(dp),allocatable :: reference(:,:,:,:),q(:,:,:,:),rhs(:,:,:,:),q0(:,:,:,:),gathered(:,:,:,:)
  real(dp),allocatable :: ref_rhs(:,:,:,:),ref_q0(:,:,:,:)
  real(dp) :: phase,rho,t,p,y(3),u(3)
  integer :: mode,i,j,k,fj,fk,nj,nk,ierr,steps,v,substeps
  character(len=512) :: path
  character(len=16),parameter :: kinds(4) = &
    [character(len=16) :: 'periodic','reflective','dirichlet','non_reflecting']
  call mc_mpi_start()
  call get_command_argument(1,path)
  call read_mc_config(trim(path),model)
  call initialize_mc_state_layout(layout,model%nspecies)
  call configure_mc_thermodynamics(trim(path),model%nspecies,model%species_names)
  call configure_mc_transport(trim(path),model%nspecies,model%species_names)
  call configure_mc_chemistry(trim(path),model%nspecies,model%species_names)
  call read_mc_euler_config(trim(path),model%nspecies,config)
  call read_mc_reactive_config(trim(path),reactive)
  config%nx=5
  config%ny=7
  config%nz=5
  config%boundary_reference_densities=1.0_dp
  config%boundary_reference_pressures=300000.0_dp
  config%boundary_reference_velocities=0.0_dp
  config%boundary_reference_mass_fractions=0.0_dp
  config%boundary_reference_mass_fractions(:,1)=0.3_dp
  config%boundary_reference_mass_fractions(:,2)=0.4_dp
  config%boundary_reference_mass_fractions(:,3)=0.3_dp
  allocate(reference(config%nx,config%ny,config%nz,layout%nvariables))
  allocate(ref_rhs,mold=reference)
  allocate(ref_q0,mold=reference)
  do mode=1,6
    config%boundary_face_types=kinds(min(mode,4))
    if (mode == 5) config%boundary_face_types=[character(len=32) :: &
      'reflective','non_reflecting','dirichlet','non_reflecting','reflective','dirichlet']
    if(mode==6) then
      config%geometry='planar_nozzle'
      config%y_min=-1
      config%y_max=1
      config%boundary_face_types=[character(len=32):: &
        'dirichlet','non_reflecting','reflective','reflective','periodic','periodic']
    end if
    do k=1,config%nz
      do j=1,config%ny
        do i=1,config%nx
          phase=real(i+2*j+3*k,dp)
          rho=1.0_dp+0.02_dp*sin(phase)
          t=950.0_dp+100.0_dp*cos(phase)
          y=[0.3_dp+0.01_dp*cos(phase),0.4_dp,0.3_dp-0.01_dp*cos(phase)]
          u=[sin(phase),cos(phase),2.0_dp*sin(phase)]
          p=rho*mc_mixture_gas_constant(y,layout)*t
          call set_mc_euler_conservative_state(reference(i,j,k,:),layout,config%gamma,rho,u,p,y)
        end do
      end do
    end do
    call initialize_mc_pencil(domain,config,layout,[0,0])
    fj=domain%first(1)
    fk=domain%first(2)
    nj=domain%count(1)
    nk=domain%count(2)
    allocate(q(config%nx,nj,nk,layout%nvariables))
    allocate(q0,mold=q)
    allocate(rhs,mold=q)
    q=reference(:,fj:fj+nj-1,fk:fk+nk-1,:)
    call compute_mc_pencil_rhs(q,rhs,layout,domain)
    call gather_mc_pencil(rhs,gathered,config,domain)
    if (mc_rank == 0) then
      call compute_mc_navier_stokes_rhs(reference,ref_rhs,layout,config)
      call compare(gathered,ref_rhs,2.0e-10_dp)
      deallocate(gathered)
    end if
    do steps=1,2
      call advance_mc_pencil_strang(q,q0,rhs,1.0e-7_dp,layout,reactive,domain,substeps)
      if (mc_rank == 0) call advance_mc_reactive_strang(reference,ref_q0,ref_rhs, &
        1.0e-7_dp,layout,config,reactive)
    end do
    call gather_mc_pencil(q,gathered,config,domain)
    if (mc_rank == 0) then
      call compare(gathered,reference,2.0e-10_dp)
      deallocate(gathered)
    end if
    call MPI_Comm_free(domain%comm,ierr)
    deallocate(q,q0,rhs)
  end do
  if (mc_rank == 0) write(*,'(A)') 'MPI pencil transverse/corner and Strang regression passed.'
  call mc_mpi_finish()
contains
  subroutine compare(actual,expected,tolerance)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp),intent(in) :: actual(:,:,:,:),expected(:,:,:,:),tolerance
    do v=1,size(actual,4)
      if (.not. all(ieee_is_finite(actual(:,:,:,v))) .or. &
          maxval(abs(actual(:,:,:,v)-expected(:,:,:,v))) > &
          tolerance*max(1.0_dp,maxval(abs(expected(:,:,:,v))))) then
        write(*,*) 'Mismatch mode, variable, error:',mode,v,maxval(abs(actual(:,:,:,v)-expected(:,:,:,v)))
        call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      end if
    end do
  end subroutine
end program
