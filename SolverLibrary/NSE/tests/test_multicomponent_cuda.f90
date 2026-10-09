program test_multicomponent_cuda
  use iso_c_binding,only:c_ptr
  use mod_precision,only:dp
  use mod_mc_config,only:mc_config,read_mc_config
  use mod_mc_state_layout,only:mc_state_layout,initialize_mc_state_layout
  use mod_mc_euler_config,only:mc_euler_config,read_mc_euler_config
  use mod_mc_reactive_config,only:mc_reactive_config,read_mc_reactive_config
  use mod_mc_thermodynamics_provider,only:configure_mc_thermodynamics
  use mod_mc_transport_provider,only:configure_mc_transport
  use mod_mc_chemistry_provider,only:configure_mc_chemistry
  use mod_mc_euler_field,only:set_mc_euler_conservative_state,compute_mc_euler_totals
  use mod_mc_viscous_flux,only:compute_mc_navier_stokes_rhs
  use mod_mc_reactive_solver,only:advance_mc_reactive_strang,compute_mc_reactive_timestep,apply_mc_field_chemistry
  use mod_mc_cuda
  implicit none
  type(mc_config)::model
  type(mc_state_layout)::l
  type(mc_euler_config)::c
  type(mc_reactive_config)::r
  type(c_ptr)::gpu
  real(dp),allocatable::q(:,:,:,:),q0(:,:,:,:),rhs(:,:,:,:),gq(:,:,:,:),grhs(:,:,:,:)
  real(dp)::dt,gdt,err,scale,y(3),u(3),totals(7),initial(7),theta
  integer::mode,i,j,k,s,stage,step,substeps,counts(2),cpu_substeps
  character(len=512)::path
  character(len=32)::failure_case
  call get_command_argument(1,path)
  call get_command_argument(2,failure_case)
  call read_mc_config(trim(path),model)
  call initialize_mc_state_layout(l,model%nspecies)
  call configure_mc_thermodynamics(trim(path),model%nspecies,model%species_names)
  call configure_mc_transport(trim(path),model%nspecies,model%species_names)
  call configure_mc_chemistry(trim(path),model%nspecies,model%species_names)
  call read_mc_reactive_config(trim(path),r)
  do mode=1,6
    call read_mc_euler_config(trim(path),model%nspecies,c)
    c%nx=7;c%ny=5;c%nz=3;c%dt=0
    c%boundary_face_types='periodic'
    if(mode==2)c%boundary_face_types='reflective'
    if(mode==3)c%boundary_face_types='dirichlet'
    if(mode==4)c%boundary_face_types='non_reflecting'
    if(mode>=5)then
      c%geometry='planar_nozzle';c%y_min=-1;c%y_max=1
      c%nozzle_inlet_half_height=.4_dp;c%nozzle_throat_half_height=.2_dp;c%nozzle_exit_half_height=.6_dp
      c%boundary_face_types='reflective';c%boundary_face_types(5:6)='periodic'
      if(mode==6)then
        c%boundary_face_types(1)='dirichlet';c%boundary_face_types(2)='non_reflecting'
      end if
    end if
    ! All six references must be valid even if the fixture only specifies x faces.
    do i=1,6
      c%boundary_reference_densities(i)=1
      c%boundary_reference_pressures(i)=300000
      c%boundary_reference_velocities(:,i)=[1.0_dp,.2_dp,.3_dp]
      c%boundary_reference_mass_fractions(i,1:3)=[.3_dp,.4_dp,.3_dp]
    end do
    allocate(q(c%nx,c%ny,c%nz,l%nvariables),q0(c%nx,c%ny,c%nz,l%nvariables), &
      rhs(c%nx,c%ny,c%nz,l%nvariables),gq(c%nx,c%ny,c%nz,l%nvariables),grhs(c%nx,c%ny,c%nz,l%nvariables))
    do k=1,c%nz
      do j=1,c%ny
        do i=1,c%nx
          theta=real(i+2*j+3*k,dp)
          y=[.3_dp+.01_dp*sin(theta),.4_dp,.3_dp-.01_dp*sin(theta)]
          u=[sin(theta),.2_dp*cos(theta),.3_dp*sin(theta/2)]
          call set_mc_euler_conservative_state(q(i,j,k,:),l,c%gamma,1+.01_dp*cos(theta), &
            u,300000.0_dp*(1+.01_dp*sin(theta)),y)
        end do
      end do
    end do
    if(failure_case=='bad_device')then
      call mc_gpu_create(gpu,l,c,r,999999)
      stop 0
    end if
    if(failure_case=='bad_state')q(1,1,1,1)=-1.0_dp
    if(failure_case=='bad_temperature')q(1,1,1,l%nvariables)=1e20_dp
    call mc_gpu_create(gpu,l,c,r,0);call mc_gpu_upload(gpu,q)
    if(failure_case=='substep_limit')then
      call mc_gpu_chemistry(gpu,1e6_dp,1,substeps)
      stop 0
    end if
    if(failure_case/='')then
      gdt=mc_gpu_timestep(gpu,.true.)
      stop 0
    end if
    dt=compute_mc_reactive_timestep(q,l,c,r);gdt=mc_gpu_timestep(gpu,.true.)
    if(abs(dt-gdt)/dt>1e-8_dp)error stop 'CUDA timestep differs from CPU'
    call compute_mc_navier_stokes_rhs(q,rhs,l,c)
    call mc_gpu_stage(gpu,0.0_dp,0);call mc_gpu_download(gpu,grhs,.true.)
    err=0
    do s=1,l%nvariables
      scale=max(1.0_dp,maxval(abs(rhs(:,:,:,s))))
      err=max(err,maxval(abs(grhs(:,:,:,s)-rhs(:,:,:,s)))/scale)
    end do
    print *, 'CUDA RHS mode/error ',mode,err
    if(err>2e-8_dp)error stop 'CUDA RHS differs from CPU'
    call compute_mc_euler_totals(q,c,initial)
    dt=min(dt,1e-7_dp)
    do step=1,4
      call advance_mc_reactive_strang(q,q0,rhs,dt,l,c,r)
      call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,substeps)
      call mc_gpu_begin(gpu)
      do stage=1,3
        call mc_gpu_stage(gpu,dt,stage)
      end do
      call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,substeps)
    end do
    call mc_gpu_transfer_counts(gpu,counts)
    if(any(counts/=[1,1]))error stop 'full-field transfer occurred during resident integration'
    call mc_gpu_download(gpu,gq)
    err=0
    do s=1,l%nvariables
      err=max(err,maxval(abs(gq(:,:,:,s)-q(:,:,:,s)))/max(1.0_dp,maxval(abs(q(:,:,:,s)))))
    end do
    print *, 'CUDA Strang mode/error ',mode,err
    if(err>2e-8_dp)error stop 'CUDA Strang update differs from CPU'
    if(mode==1.or.mode==2.or.mode==5)then
      call compute_mc_euler_totals(gq,c,totals)
      if(abs(sum(totals(1:3)-initial(1:3)))>1e-11_dp)error stop 'CUDA total mass conservation'
      if(abs(totals(7)-initial(7))/abs(initial(7))>1e-11_dp)error stop 'CUDA total energy conservation'
    end if
    if(mode==1)then
      call apply_mc_field_chemistry(q,1.0_dp,l,c,r,cpu_substeps)
      call mc_gpu_chemistry(gpu,1.0_dp,r%maximum_chemistry_substeps,substeps)
      if(substeps<2.or.cpu_substeps<2)error stop 'subcycling regression did not subcycle'
      call mc_gpu_transfer_counts(gpu,counts)
      if(any(counts/=[1,2]))error stop 'full-field transfer occurred inside chemistry'
      call mc_gpu_download(gpu,gq)
      err=maxval(abs(q-gq)/max(1.0_dp,abs(q)))
      print *, 'CUDA chemistry subcycling error/steps ',err,substeps
      if(err>2e-8_dp)error stop 'CUDA subcycling differs from CPU'
    end if
    call mc_gpu_destroy(gpu)
    deallocate(q,q0,rhs,gq,grhs)
  end do
  print *, 'Resident CUDA CPU comparison passed'
end program
