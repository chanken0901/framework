program main_rf_flow1d
  use mod_rf_flow1d
  use mod_rf_thermo
  use mod_rf_transport
  implicit none
  type(rf_mechanism) :: m
  type(rf_transport) :: transport
  character(16) :: transport_model='none'
  character(16) :: reconstruction='first_order'
  real(dp) :: viscosity=0,bulk_viscosity=0,thermal_conductivity=0,mass_diffusivity=0
  character(2048) :: mechanism_file,input_file,output_file
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
  namelist /flow1d/ nx,length,interface_x,end_time,cfl,max_dt,max_steps,write_every,left_bc,right_bc, &
    left_temperature,left_pressure,left_velocity,right_temperature,right_pressure,right_velocity,chemistry, &
    chemistry_rtol,chemistry_atol_species,chemistry_atol_temperature,chemistry_max_steps, &
    transport_model,viscosity,bulk_viscosity,thermal_conductivity,mass_diffusivity,reconstruction
  call require(command_argument_count()==3,'Usage: rf_flow1d mechanism.rf flow.in output.csv')
  call get_command_argument(1,mechanism_file)
  call get_command_argument(2,input_file)
  call get_command_argument(3,output_file)
  call read_mechanism(trim(mechanism_file),m)
  ns=size(m%species)
  open(newunit=io,file=trim(input_file),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open 1D input')
  read(io,nml=flow1d,iostat=ios)
  call require(ios==0,'Invalid flow1d namelist')
  call validate_reconstruction(reconstruction)
  transport=rf_transport(viscosity,bulk_viscosity,thermal_conductivity,mass_diffusivity)
  call validate_transport(transport)
  call require(transport_model=='none'.or.transport_model=='constant','Unknown transport model')
  if(transport_model=='none') call require(.not.transport_active(transport), &
    'Nonzero transport coefficients require transport_model=constant')
  call require(nx>=2.and.max_steps>0.and.write_every>0,'Invalid grid/step/output count')
  call require(all(ieee_is_finite([length,interface_x,end_time,cfl,max_dt])), 'Nonfinite flow control')
  call require(min(length,end_time,max_dt,cfl)>0.and.cfl<=.5_dp,'Require positive controls, CFL<=0.5')
  call require(interface_x>=0.and.interface_x<=length,'Interface outside domain')
  call require(left_bc=='periodic'.or.left_bc=='outflow'.or.left_bc=='reflecting','Invalid left boundary')
  call require(right_bc=='periodic'.or.right_bc=='outflow'.or.right_bc=='reflecting','Invalid right boundary')
  call require((left_bc=='periodic').eqv.(right_bc=='periodic'),'Periodic boundary must be paired')
  allocate(yl(ns),yr(ns),y(ns),q(ns+2,nx),initial(ns+2),boundary(ns+2),change(ns+2),delta(ns+2))
  allocate(elements0(m%ne),elements(m%ne))
  read(io,*,iostat=ios) yl
  call require(ios==0,'Expected full left mass-fraction row')
  read(io,*,iostat=ios) yr
  call require(ios==0,'Expected full right mass-fraction row')
  close(io)
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
  write(out,'(a)') '# reconstruction='//trim(reconstruction)
  write(out,'(a,es25.16e3)') '# viscosity=',viscosity
  write(out,'(a,es25.16e3)') '# bulk_viscosity=',bulk_viscosity
  write(out,'(a,es25.16e3)') '# thermal_conductivity=',thermal_conductivity
  write(out,'(a,es25.16e3)') '# mass_diffusivity=',mass_diffusivity
  write(out,'(a)') '# source_sha256='//m%source_hash
  write(out,'(a)') '# canonical_sha256='//m%canonical_hash
  write(out,'(a,l1)') '# chemistry=',chemistry
  write(out,'(a,es25.16e3)') '# cfl=',cfl
  write(out,'(a,es25.16e3)') '# max_dt=',max_dt
  write(out,'(a,es25.16e3)') '# chemistry_rtol=',chemistry_rtol
  write(out,'(a,es25.16e3)') '# chemistry_atol_species=',chemistry_atol_species
  write(out,'(a,es25.16e3)') '# chemistry_atol_temperature=',chemistry_atol_temperature
  write(out,'(a)') '# left_bc='//trim(left_bc)//' right_bc='//trim(right_bc)
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
    dt=min(max_dt,end_time-time,flow_timestep(m,q,dx,cfl,transport,reconstruction,left_bc,right_bc))
    call require(time+dt>time,'Flow timestep underflow')
    call advance_flow(m,q,dx,dt,cfl,left_bc,right_bc,chemistry,chemistry_rtol,chemistry_atol_species, &
                      chemistry_atol_temperature,chemistry_max_steps,change,transport,reconstruction)
    boundary=boundary+change
    time=time+dt; step=step+1
  end do
  write(out,'(a,es25.16e3)') '# mass_error=',mass_error
  write(out,'(a,es25.16e3)') '# momentum_error=',momentum_error
  write(out,'(a,es25.16e3)') '# energy_error=',energy_error
  write(out,'(a,es25.16e3)') '# element_error=',element_error
  write(out,'(a)') '# SUCCESS'
  close(out)
  write(*,'(a)') '[OK] Fortran 1D flow completed: '//trim(output_file)
end program
