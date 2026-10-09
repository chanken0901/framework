! Exercise variable Ns, including the maximum, through the normal providers.
program test_multicomponent_cuda_species
  use iso_c_binding,only:c_ptr
  use mod_precision,only:dp
  use mod_mc_state_layout,only:mc_state_layout,initialize_mc_state_layout
  use mod_mc_euler_config,only:mc_euler_config
  use mod_mc_reactive_config,only:mc_reactive_config
  use mod_mc_euler_field,only:set_mc_euler_conservative_state
  use mod_mc_reactive_solver,only:advance_mc_reactive_strang
  use mod_mc_thermodynamics_provider,only:configure_mc_thermodynamics
  use mod_mc_transport_provider,only:configure_mc_transport
  use mod_mc_chemistry_provider,only:configure_mc_chemistry
  use mod_mc_cuda
  implicit none
  type(mc_state_layout)::l
  type(mc_euler_config)::c
  type(mc_reactive_config)::r
  type(c_ptr)::gpu
  character(len=32)::thermo_species_names(64),transport_species_names(64),chemistry_species_names(64)
  real(dp)::molecular_weights(64),temperature_midpoints(64),nasa_low_coefficients(7,64),nasa_high_coefficients(7,64)
  real(dp)::reference_dynamic_viscosity,prandtl_number,species_diffusivities(64)
  real(dp)::reactant_stoich(64),product_stoich(64),reaction_orders(64)
  real(dp)::pre_exponential_factor,temperature_exponent,activation_temperature
  real(dp),allocatable::q(:,:,:,:),q0(:,:,:,:),rhs(:,:,:,:),actual(:,:,:,:),y(:)
  real(dp)::dt,err
  integer::n,trial,i,s,unit,stage,steps
  integer,parameter::sizes(3)=[2,4,64]
  namelist /thermally_perfect/ thermo_species_names,molecular_weights,temperature_midpoints, &
    nasa_low_coefficients,nasa_high_coefficients
  namelist /mixture_averaged_transport/ transport_species_names,reference_dynamic_viscosity, &
    prandtl_number,species_diffusivities
  namelist /one_step_arrhenius/ chemistry_species_names,reactant_stoich,product_stoich,reaction_orders, &
    pre_exponential_factor,temperature_exponent,activation_temperature
  do trial=1,3
    n=sizes(trial)
    do s=1,64
      write(thermo_species_names(s),'(A,I0)')'species',s
    end do
    transport_species_names=thermo_species_names;chemistry_species_names=thermo_species_names
    molecular_weights=28;temperature_midpoints=1000
    nasa_low_coefficients=0;nasa_low_coefficients(1,:)=3.5_dp;nasa_low_coefficients(6,2)=-5000
    nasa_high_coefficients=nasa_low_coefficients
    reference_dynamic_viscosity=1.8e-5_dp;prandtl_number=.72_dp;species_diffusivities=2e-5_dp
    reactant_stoich=0;reactant_stoich(1)=1;product_stoich=0;product_stoich(2)=1
    reaction_orders=reactant_stoich;pre_exponential_factor=1e3_dp;temperature_exponent=0;activation_temperature=2000
    open(newunit=unit,file='cuda_species_test.dat',status='replace')
    write(unit,nml=thermally_perfect)
    write(unit,nml=mixture_averaged_transport)
    write(unit,nml=one_step_arrhenius)
    close(unit)
    call configure_mc_thermodynamics('cuda_species_test.dat',n,thermo_species_names)
    call configure_mc_transport('cuda_species_test.dat',n,thermo_species_names)
    call configure_mc_chemistry('cuda_species_test.dat',n,thermo_species_names)
    open(newunit=unit,file='cuda_species_test.dat',status='old');close(unit,status='delete')
    call initialize_mc_state_layout(l,n)
    c=mc_euler_config();c%nx=3;c%ny=2;c%nz=1
    allocate(q(3,2,1,l%nvariables),q0(3,2,1,l%nvariables),rhs(3,2,1,l%nvariables), &
      actual(3,2,1,l%nvariables),y(n))
    y=1.0_dp/n
    do i=1,3
      y(1)=(1+.01_dp*sin(real(i,dp)))/n;y(2)=2.0_dp/n-y(1)
      do s=1,2
        call set_mc_euler_conservative_state(q(i,s,1,:),l,c%gamma,1.0_dp, &
          [real(i,dp),.2_dp,.1_dp],300000.0_dp,y)
      end do
    end do
    call mc_gpu_create(gpu,l,c,r,0);call mc_gpu_upload(gpu,q)
    dt=1e-7_dp
    call advance_mc_reactive_strang(q,q0,rhs,dt,l,c,r)
    call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,steps)
    call mc_gpu_begin(gpu)
    do stage=1,3
      call mc_gpu_stage(gpu,dt,stage)
    end do
    call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,steps)
    call mc_gpu_download(gpu,actual)
    err=maxval(abs(q-actual)/max(1.0_dp,abs(q)))
    print *, 'CUDA variable species count/error ',n,err
    if(err>2e-8_dp)error stop 'CUDA variable species comparison failed'
    call mc_gpu_destroy(gpu)
    deallocate(q,q0,rhs,actual,y)
  end do
  print *, 'CUDA variable species comparison passed'
end program
