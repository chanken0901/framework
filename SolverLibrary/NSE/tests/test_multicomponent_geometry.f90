program test_multicomponent_geometry
  use mod_precision, only: dp
  use mod_mc_config, only: mc_config,read_mc_config
  use mod_mc_state_layout, only: mc_state_layout,initialize_mc_state_layout
  use mod_mc_euler_config, only: mc_euler_config,read_mc_euler_config
  use mod_mc_geometry
  use mod_mc_euler_field, only: set_mc_euler_conservative_state,compute_mc_euler_totals,validate_mc_euler_state
  use mod_mc_mapped_flux, only: mc_mapped_rhs
  use mod_mc_viscous_flux, only: compute_mc_navier_stokes_rhs
  use mod_mc_reactive_config, only: mc_reactive_config,read_mc_reactive_config
  use mod_mc_reactive_solver, only: advance_mc_reactive_strang
  use mod_mc_thermodynamics_provider, only: configure_mc_thermodynamics
  use mod_mc_transport_provider, only: configure_mc_transport
  use mod_mc_chemistry_provider, only: configure_mc_chemistry
  implicit none
  type(mc_config)::model
  type(mc_state_layout)::l
  type(mc_euler_config)::c
  type(mc_reactive_config)::reactive
  real(dp),allocatable::q(:,:,:,:),rhs(:,:,:,:),q0(:,:,:,:)
  real(dp)::area(3),initial(7),final(7),y(3),u(3),position(3)
  integer::i,j,k,d,idx(3)
  character(len=512)::path
  call get_command_argument(1,path)
  call read_mc_config(trim(path),model)
  call initialize_mc_state_layout(l,model%nspecies)
  call configure_mc_thermodynamics(trim(path),model%nspecies,model%species_names)
  call configure_mc_transport(trim(path),model%nspecies,model%species_names)
  call configure_mc_chemistry(trim(path),model%nspecies,model%species_names)
  call read_mc_euler_config(trim(path),model%nspecies,c)
  call read_mc_reactive_config(trim(path),reactive)
  c%nx=16
  c%ny=8
  c%nz=4
  c%y_min=-1
  c%y_max=1
  c%geometry='planar_nozzle'
  c%nozzle_inlet_half_height=0.4_dp
  c%nozzle_throat_half_height=0.2_dp
  c%nozzle_exit_half_height=0.6_dp
  c%boundary_face_types=[character(len=32)::'reflective','reflective','reflective','reflective','periodic','periodic']
  allocate(q(c%nx,c%ny,c%nz,l%nvariables),rhs(c%nx,c%ny,c%nz,l%nvariables),q0(c%nx,c%ny,c%nz,l%nvariables))
  y=[0.3_dp,0.4_dp,0.3_dp]
  do k=1,c%nz
    do j=1,c%ny
      do i=1,c%nx
        area=0
        do d=1,3
          idx=[i,j,k]
          area=area+mc_face_area(c,idx,d)
          idx(d)=idx(d)-1
          area=area-mc_face_area(c,idx,d)
        end do
        if(maxval(abs(area))>1e-14_dp) error stop 'discrete metric identity failed'
        if(mc_cell_volume(c,i,j,k)<=0) error stop 'invalid cell volume'
        call set_mc_euler_conservative_state(q(i,j,k,:),l,c%gamma,1.0_dp, &
          [0.0_dp,0.0_dp,1.0_dp],300000.0_dp,y)
      end do
    end do
  end do
  call mc_mapped_rhs(q,rhs,l,c,.false.)
  if(maxval(abs(rhs))>1e-7_dp) error stop 'mapped free stream preservation failed'
  call mc_mapped_rhs(q,rhs,l,c,.true.)
  if(maxval(abs(rhs))>1e-7_dp) error stop 'mapped uniform transport RHS failed'
  do k=1,c%nz
    do j=1,c%ny
      do i=1,c%nx
        position=mc_cell_center(c,i,j,k)
        u=[sin(position(1)),0.2_dp*cos(position(2)),0.1_dp*sin(position(3))]
        call set_mc_euler_conservative_state(q(i,j,k,:),l,c%gamma, &
          1.0_dp+0.01_dp*sin(position(1)),u,300000.0_dp,y)
      end do
    end do
  end do
  call compute_mc_navier_stokes_rhs(q,rhs,l,c)
  call compute_mc_euler_totals(rhs,c,final)
  if(abs(sum(final(1:3)))>1e-10_dp .or. abs(final(7))>1e-7_dp) &
    error stop 'closed mapped domain RHS does not conserve mass/energy'
  call compute_mc_euler_totals(q,c,initial)
  do i=1,3
    call advance_mc_reactive_strang(q,q0,rhs,1e-7_dp,l,c,reactive)
  end do
  call compute_mc_euler_totals(q,c,final)
  if(abs(sum(final(1:3)-initial(1:3)))>1e-11_dp .or. &
      abs(final(7)-initial(7))/abs(initial(7))>1e-12_dp) error stop 'mapped Strang conservation failed'
  call validate_mc_euler_state(q,l,c)
  print *, 'Mapped nozzle geometry, free-stream, transport and reactive conservation passed.'
end program
