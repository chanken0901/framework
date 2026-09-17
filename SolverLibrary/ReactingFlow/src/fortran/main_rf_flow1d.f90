program main_rf_flow1d
  use mod_rf_flow1d
  use mod_rf_thermo
  use mod_rf_transport
  implicit none
  type(rf_mechanism) :: m
  type(rf_transport) :: transport
  character(16) :: transport_model='none'
  character(16) :: viscosity_model='constant'
  character(16) :: conductivity_model='constant'
  real(dp) :: viscosity_reference_temperature=-1
  real(dp), allocatable :: species_reference_viscosities(:),species_sutherland_temperatures(:)
  character(16) :: transport_temperature_model='constant'
  real(dp) :: transport_reference_temperature=-1,transport_temperature_exponent=-1
  character(16) :: reconstruction='first_order'
  real(dp) :: viscosity=0,bulk_viscosity=0,thermal_conductivity=0,mass_diffusivity=0
  real(dp), allocatable :: species_diffusivities(:)
  real(dp) :: left_boundary_temperature=-1,left_boundary_pressure=-1,left_boundary_velocity=0
  real(dp) :: right_boundary_temperature=-1,right_boundary_pressure=-1,right_boundary_velocity=0
  real(dp), allocatable :: left_boundary_y(:),right_boundary_y(:),fixed_states(:,:)
  character(2048) :: mechanism_file,input_file,output_file
  character(512) :: input_error=''
  character(2048) :: transport_file='',resolved_transport_file=''
  character(16) :: left_bc='outflow',right_bc='outflow'
  integer :: nx=100,max_steps=100000,chemistry_max_steps=100000,write_every=50
  real(dp) :: length=1,interface_x=.5_dp,end_time=.001_dp,cfl=.4_dp,max_dt=.00001_dp
  real(dp) :: left_temperature=1100,left_pressure=101325,left_velocity=0
  real(dp) :: right_temperature=1100,right_pressure=101325,right_velocity=0
  real(dp) :: chemistry_rtol=1.e-9_dp,chemistry_atol_species=1.e-16_dp,chemistry_atol_temperature=1.e-8_dp
  logical :: chemistry=.false.
  real(dp), allocatable :: q(:,:),yl(:),yr(:),initial(:),boundary(:),change(:)
  real(dp) :: time,dt,dx,rho,u,t,p,a,mass_error,momentum_error,energy_error,element_error,x
  real(dp), allocatable :: y(:),delta(:),elements0(:),elements(:)
  integer :: ns,io,out,ios,i,j,step
  integer :: rejected_steps,total_rejected=0
  namelist /flow1d/ nx,length,interface_x,end_time,cfl,max_dt,max_steps,write_every,left_bc,right_bc, &
    left_temperature,left_pressure,left_velocity,right_temperature,right_pressure,right_velocity,chemistry, &
    chemistry_rtol,chemistry_atol_species,chemistry_atol_temperature,chemistry_max_steps, &
    transport_model,viscosity,bulk_viscosity,thermal_conductivity,mass_diffusivity,reconstruction,species_diffusivities, &
    left_boundary_temperature,left_boundary_pressure,left_boundary_velocity,left_boundary_y, &
    right_boundary_temperature,right_boundary_pressure,right_boundary_velocity,right_boundary_y, &
    transport_temperature_model,transport_reference_temperature,transport_temperature_exponent, &
    viscosity_model,viscosity_reference_temperature,species_reference_viscosities,species_sutherland_temperatures, &
    conductivity_model,transport_file
  call require(command_argument_count()==3,'Usage: rf_flow1d mechanism.rf flow.in output.csv')
  call get_command_argument(1,mechanism_file)
  call get_command_argument(2,input_file)
  call get_command_argument(3,output_file)
  call read_mechanism(trim(mechanism_file),m)
  ns=size(m%species)
  allocate(species_diffusivities(ns)); species_diffusivities=-1
  allocate(species_reference_viscosities(ns),species_sutherland_temperatures(ns))
  species_reference_viscosities=-1;species_sutherland_temperatures=-1
  allocate(left_boundary_y(ns),right_boundary_y(ns),fixed_states(ns+2,2))
  left_boundary_y=-1; right_boundary_y=-1; fixed_states=0
  open(newunit=io,file=trim(input_file),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open 1D input')
  read(io,nml=flow1d,iostat=ios,iomsg=input_error)
  call require(ios==0,'Invalid flow1d namelist: '//trim(input_error))
  call validate_reconstruction(reconstruction)
  transport=rf_transport(viscosity,bulk_viscosity,thermal_conductivity,mass_diffusivity)
  call require(conductivity_model=='constant'.or.conductivity_model=='eucken_wms','Unknown conductivity model')
  transport%eucken_wms=conductivity_model=='eucken_wms'
  select case(viscosity_model)
  case('constant')
    call require(len_trim(transport_file)==0,'transport_file requires sutherland_wilke')
    call require(viscosity_reference_temperature==-1.and.all(species_reference_viscosities==-1).and. &
      all(species_sutherland_temperatures==-1),'Species viscosity parameters require sutherland_wilke')
  case('sutherland_wilke')
    call require(transport_model/='none'.and.transport_temperature_model=='constant', &
      'sutherland_wilke requires transport and cannot combine with common power_law')
    if(len_trim(transport_file)>0) then
      call require(viscosity_reference_temperature==-1.and.all(species_reference_viscosities==-1).and. &
        all(species_sutherland_temperatures==-1),'Transport file cannot combine with inline viscosity data')
      resolved_transport_file=resolve_transport_path(trim(input_file),trim(transport_file))
      call read_transport_data(trim(resolved_transport_file),m,transport)
      viscosity_reference_temperature=transport%viscosity_reference_temperature
      species_reference_viscosities=transport%species_viscosity
      species_sutherland_temperatures=transport%sutherland_temperature
    else
      transport%species_viscosity=species_reference_viscosities
      transport%sutherland_temperature=species_sutherland_temperatures
      transport%viscosity_reference_temperature=viscosity_reference_temperature
    end if
  case default
    call require(.false.,'Unknown viscosity model')
  end select
  select case(transport_temperature_model)
  case('constant')
    call require(transport_reference_temperature==-1.and.transport_temperature_exponent==-1, &
      'Temperature parameters require power_law transport')
  case('power_law')
    call require(transport_model/='none','power_law requires active transport model')
    transport%reference_temperature=transport_reference_temperature
    transport%temperature_exponent=transport_temperature_exponent
  case default
    call require(.false.,'Unknown transport temperature model')
  end select
  call require(transport_model=='none'.or.transport_model=='constant'.or. &
    transport_model=='species_constant','Unknown transport model')
  if(transport_model=='species_constant') then
    ! -1 sentinels make incomplete lists fail rather than silently supplying zero diffusivity.
    transport%species_diffusivity=species_diffusivities
  else
    call require(all(species_diffusivities==-1),'species_diffusivities requires transport_model=species_constant')
  end if
  call validate_transport(transport,ns)
  if(transport_model=='none') call require(.not.transport_active(transport), &
    'Nonzero transport coefficients require constant or species_constant transport')
  call require(nx>=2.and.max_steps>0.and.write_every>0,'Invalid grid/step/output count')
  call require(all(ieee_is_finite([length,interface_x,end_time,cfl,max_dt])), 'Nonfinite flow control')
  call require(min(length,end_time,max_dt,cfl)>0.and.cfl<=.5_dp,'Require positive controls, CFL<=0.5')
  call require(interface_x>=0.and.interface_x<=length,'Interface outside domain')
  call require(left_bc=='periodic'.or.left_bc=='outflow'.or.left_bc=='reflecting'.or. &
    left_bc=='dirichlet','Invalid left boundary')
  call require(right_bc=='periodic'.or.right_bc=='outflow'.or.right_bc=='reflecting'.or. &
    right_bc=='dirichlet','Invalid right boundary')
  call require((left_bc=='periodic').eqv.(right_bc=='periodic'),'Periodic boundary must be paired')
  allocate(yl(ns),yr(ns),y(ns),q(ns+2,nx),initial(ns+2),boundary(ns+2),change(ns+2),delta(ns+2))
  allocate(elements0(m%ne),elements(m%ne))
  read(io,*,iostat=ios) yl
  call require(ios==0,'Expected full left mass-fraction row')
  read(io,*,iostat=ios) yr
  call require(ios==0,'Expected full right mass-fraction row')
  close(io)
  call initialize_boundary(left_bc,left_boundary_temperature,left_boundary_pressure,left_boundary_velocity, &
    left_boundary_y,fixed_states(:,1))
  call initialize_boundary(right_bc,right_boundary_temperature,right_boundary_pressure,right_boundary_velocity, &
    right_boundary_y,fixed_states(:,2))
  dx=length/nx
  do i=1,nx
    x=(i-.5_dp)*dx
    if(x<interface_x) then
      call primitive_to_conserved(m,left_temperature,left_pressure,left_velocity,yl,q(:,i))
    else
      call primitive_to_conserved(m,right_temperature,right_pressure,right_velocity,yr,q(:,i))
    end if
  end do
  initial=sum(q,dim=2)*dx; boundary=0; elements0=0
  do j=1,ns
    elements0=elements0+initial(j)/m%species(j)%mass*m%species(j)%atoms
  end do
  open(newunit=out,file=trim(output_file),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create output; existing files are never overwritten')
  write(out,'(a)') '# 1D conservative flow + optional transport + Strang/DVODE chemistry, SI units'
  write(out,'(a)') '# transport_model='//trim(transport_model)
  if(len_trim(resolved_transport_file)>0) write(out,'(a)') '# transport_file='//trim(resolved_transport_file)
  write(out,'(a)') '# viscosity_model='//trim(viscosity_model)
  write(out,'(a)') '# conductivity_model='//trim(conductivity_model)
  if(viscosity_model=='sutherland_wilke') then
    write(out,'(a,es25.16e3)') '# viscosity_reference_temperature=',viscosity_reference_temperature
    do j=1,ns
      write(out,'(a,es25.16e3)') '# reference_viscosity_'//trim(m%species(j)%name)//'=',species_reference_viscosities(j)
      write(out,'(a,es25.16e3)') '# sutherland_temperature_'//trim(m%species(j)%name)//'=',species_sutherland_temperatures(j)
    end do
  end if
  write(out,'(a)') '# transport_temperature_model='//trim(transport_temperature_model)
  if(transport_temperature_model=='power_law') then
    write(out,'(a,es25.16e3)') '# transport_reference_temperature=',transport%reference_temperature
    write(out,'(a,es25.16e3)') '# transport_temperature_exponent=',transport%temperature_exponent
  end if
  write(out,'(a)') '# reconstruction='//trim(reconstruction)
  write(out,'(a,es25.16e3)') '# viscosity=',viscosity
  write(out,'(a,es25.16e3)') '# bulk_viscosity=',bulk_viscosity
  write(out,'(a,es25.16e3)') '# thermal_conductivity=',thermal_conductivity
  write(out,'(a,es25.16e3)') '# mass_diffusivity=',mass_diffusivity
  if(allocated(transport%species_diffusivity)) then
    do j=1,ns
      write(out,'(a,es25.16e3)') '# diffusivity_'//trim(m%species(j)%name)//'=',species_diffusivities(j)
    end do
  end if
  write(out,'(a)') '# source_sha256='//m%source_hash
  write(out,'(a)') '# canonical_sha256='//m%canonical_hash
  write(out,'(a,l1)') '# chemistry=',chemistry
  write(out,'(a,es25.16e3)') '# cfl=',cfl
  write(out,'(a,es25.16e3)') '# max_dt=',max_dt
  write(out,'(a,es25.16e3)') '# chemistry_rtol=',chemistry_rtol
  write(out,'(a,es25.16e3)') '# chemistry_atol_species=',chemistry_atol_species
  write(out,'(a,es25.16e3)') '# chemistry_atol_temperature=',chemistry_atol_temperature
  write(out,'(a)') '# left_bc='//trim(left_bc)//' right_bc='//trim(right_bc)
  if(left_bc=='dirichlet') write(out,'(a,*(es25.16e3,1x))') '# fixed_state_left=',fixed_states(:,1)
  if(right_bc=='dirichlet') write(out,'(a,*(es25.16e3,1x))') '# fixed_state_right=',fixed_states(:,2)
  write(out,'(a)',advance='no') 'step,time,x,density,velocity,temperature,pressure'
  do j=1,ns
    write(out,'(a)',advance='no') ',Y_'//trim(m%species(j)%name)
  end do
  write(out,*)
  step=0; time=0
  mass_error=0; momentum_error=0; energy_error=0; element_error=0
  do
    delta=sum(q,dim=2)*dx-initial-boundary
    mass_error=max(mass_error,abs(sum(delta(:ns)))/max(1._dp,abs(sum(initial(:ns)))))
    momentum_error=max(momentum_error,abs(delta(ns+1))/max(1._dp,abs(initial(ns+1))))
    energy_error=max(energy_error,abs(delta(ns+2))/max(1._dp,abs(initial(ns+2))))
    elements=0
    do j=1,ns
      elements=elements+delta(j)/m%species(j)%mass*m%species(j)%atoms
    end do
    element_error=max(element_error,maxval(abs(elements)/max(1._dp,abs(elements0))))
    call require(max(mass_error,momentum_error,energy_error,element_error)<1.e-7_dp, &
                 'Boundary-corrected global conservation failed')
    if(mod(step,write_every)==0.or.time>=end_time) then
      do i=1,nx
        call conserved_to_primitive(m,q(:,i),rho,u,t,p,a,y)
        write(out,'(i0,",",*(es25.16e3,:,","))') step,time,(i-.5_dp)*dx,rho,u,t,p,y
      end do
      flush(out)
    end if
    if(time>=end_time) exit
    call require(step<max_steps,'Flow exceeded max_steps; output is incomplete')
    dt=min(max_dt,end_time-time,flow_timestep(m,q,dx,cfl,transport,reconstruction,left_bc,right_bc,fixed_states))
    call require(time+dt>time,'Flow timestep underflow')
    call advance_flow(m,q,dx,dt,cfl,left_bc,right_bc,chemistry,chemistry_rtol,chemistry_atol_species, &
                      chemistry_atol_temperature,chemistry_max_steps,change,transport,reconstruction,rejected_steps,fixed_states)
    call require(time+dt>time,'Accepted flow timestep cannot advance time')
    total_rejected=total_rejected+rejected_steps
    boundary=boundary+change
    time=time+dt; step=step+1
  end do
  write(out,'(a,es25.16e3)') '# mass_error=',mass_error
  write(out,'(a,i0)') '# rejected_steps=',total_rejected
  write(out,'(a,es25.16e3)') '# momentum_error=',momentum_error
  write(out,'(a,es25.16e3)') '# energy_error=',energy_error
  write(out,'(a,es25.16e3)') '# element_error=',element_error
  write(out,'(a)') '# SUCCESS'
  close(out)
  write(*,'(a)') '[OK] Fortran 1D flow completed: '//trim(output_file)
contains
  function resolve_transport_path(input_path,data_path) result(path)
    character(*), intent(in) :: input_path,data_path
    character(:), allocatable :: path
    integer :: i,last
    path=data_path
    if(data_path(1:1)=='/'.or.data_path(1:1)==achar(92)) return
    if(len(data_path)>=2) then
      if(data_path(2:2)==':') then
        call require(len(data_path)>=3,'Invalid drive path')
        call require(data_path(3:3)=='/'.or.data_path(3:3)==achar(92),'Drive-relative transport path is ambiguous')
        return
      end if
    end if
    last=0
    do i=1,len(input_path)
      if(input_path(i:i)=='/'.or.input_path(i:i)==achar(92)) last=i
    end do
    path=input_path(:last)//data_path
    call require(len(path)<=len(resolved_transport_file),'Transport path too long')
  end function

  subroutine initialize_boundary(kind,bt,bp,bu,by,state)
    character(*), intent(in) :: kind
    real(dp), intent(in) :: bt,bp,bu,by(:)
    real(dp), intent(out) :: state(:)
    state=0
    if(kind=='dirichlet') then
      call require(bt>0.and.bp>0,'Dirichlet requires positive boundary temperature and pressure')
      call primitive_to_conserved(m,bt,bp,bu,by,state)
    else
      call require(bt==-1.and.bp==-1.and.bu==0.and.all(by==-1),'Boundary state supplied for non-Dirichlet face')
    end if
  end subroutine
end program
