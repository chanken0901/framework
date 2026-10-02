program main_rf_flow2d
  use mod_rf_flow2d
  use mod_rf_grid2d
  use mod_rf_finite_volume
  use mod_rf_thermo
  implicit none
  type(rf_mechanism), target :: m
  type(rf_grid2d) :: grid
  type(rf_boundary2d) :: bc
  character(2048) :: mechanism,input,output,mesh_profile='',mesh_path
  character(512) :: message
  character(16) :: boundary_kind(4)=[character(16)::'outflow','outflow','reflecting','reflecting']
  integer :: ny=16,max_steps=100000,write_every=50,chemistry_max_steps=100000
  real(dp) :: temperature=300,pressure=101325,velocity(2)=0
  real(dp) :: boundary_temperature(4)=300,boundary_pressure(4)=101325,boundary_velocity(2,4)=0
  real(dp) :: end_time=1.e-6_dp,max_dt=1.e-9_dp,cfl=.4_dp
  real(dp) :: chemistry_rtol=1.e-9_dp,chemistry_atol_species=1.e-16_dp,chemistry_atol_temperature=1.e-8_dp
  logical :: chemistry=.false.,ok
  real(dp), allocatable :: mass_fractions(:),boundary_y(:,:),q(:,:),dq(:,:),change(:),rate(:),total0(:),balance(:),y(:)
  real(dp) :: time,dt,allowed,rho,v(2),t,p,sound
  integer :: unit,out,ios,ns,nv,nc,c,b,step,rejected,total_rejected,i,slash
  namelist /flow2d/ mesh_profile,ny,temperature,pressure,velocity,mass_fractions,boundary_kind, &
    boundary_temperature,boundary_pressure,boundary_velocity,boundary_y,end_time,max_dt,cfl,max_steps, &
    write_every,chemistry,chemistry_rtol,chemistry_atol_species,chemistry_atol_temperature,chemistry_max_steps
  call require(command_argument_count()==3,'Usage: rf_flow2d mechanism.rf flow2d.in new_output.csv')
  call get_command_argument(1,mechanism);call get_command_argument(2,input);call get_command_argument(3,output)
  call read_mechanism(trim(mechanism),m)
  ns=size(m%species);nv=ns+3
  allocate(mass_fractions(ns),boundary_y(ns,4),y(ns));mass_fractions=-1;boundary_y=-1
  open(newunit=unit,file=trim(input),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open flow2d input')
  read(unit,nml=flow2d,iostat=ios,iomsg=message);close(unit)
  call require(ios==0,'Invalid flow2d namelist: '//trim(message))
  call require(len_trim(mesh_profile)>0,'mesh_profile is required')
  call require(ieee_is_finite(end_time).and.end_time>=0.and.ieee_is_finite(max_dt).and.max_dt>0,'Invalid time limits')
  call require(max_steps>0.and.write_every>0.and.chemistry_max_steps>0,'Invalid step limits')
  call require(all(ieee_is_finite([chemistry_rtol,chemistry_atol_species,chemistry_atol_temperature])).and. &
    min(chemistry_rtol,chemistry_atol_species,chemistry_atol_temperature)>0,'Invalid chemistry tolerances')
  call validate_input_y(mass_fractions)
  ! Mesh path is relative to the input file, not the shell working directory.
  mesh_path=trim(mesh_profile)
  if(mesh_profile(1:1)/='/'.and.mesh_profile(1:1)/=achar(92).and.index(mesh_profile,':')==0) then
    slash=0
    do i=1,len_trim(input)
      if(input(i:i)=='/'.or.input(i:i)==achar(92)) slash=i
    end do
    call require(slash+len_trim(mesh_profile)<=len(mesh_path),'Mesh path too long')
    mesh_path=input(:slash)//trim(mesh_profile)
  end if
  call read_nozzle2d(trim(mesh_path),ny,grid);nc=size(grid%mesh%volume)
  allocate(q(nv,nc),dq(nv,nc),change(nv),rate(nv),total0(nv),balance(nv),bc%fixed(nv,4))
  bc%kind=boundary_kind;bc%fixed=0
  do b=1,4
    if(bc%kind(b)/='dirichlet') cycle
    call validate_input_y(boundary_y(:,b))
    call primitive_nd(m,boundary_temperature(b),boundary_pressure(b),boundary_velocity(:,b),boundary_y(:,b),bc%fixed(:,b))
  end do
  do c=1,nc
    call primitive_nd(m,temperature,pressure,velocity,mass_fractions,q(:,c))
  end do
  call flow2d_rhs(m,grid,bc,q,cfl,dq,rate,allowed,ok)
  call require(ok,'Invalid initial 2D flow')
  open(newunit=out,file=trim(output),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create new flow2d output (existing path?)')
  write(out,'(a)') '# Experimental CPU serial planar reactive Euler; first-order Rusanov / SSPRK3 / Strang'
  write(out,'(a)') '# mechanism_hash='//trim(m%canonical_hash)
  write(out,'(a)',advance='no') 'step,time,cell,x,y,area,density,u,v,temperature,pressure'
  do i=1,ns
    write(out,'(a)',advance='no') ',Y_'//trim(m%species(i)%name)
  end do
  write(out,*)
  total0=0;balance=0;time=0;step=0;total_rejected=0
  do c=1,nc
    total0=total0+grid%mesh%volume(c)*q(:,c)
  end do
  call snapshot()
  do while(time<end_time)
    call require(step<max_steps,'2D maximum step count reached before end_time')
    call flow2d_rhs(m,grid,bc,q,cfl,dq,rate,allowed,ok)
    call require(ok,'Invalid accepted 2D state')
    dt=min(max_dt,allowed,end_time-time)
    call advance_flow2d(m,grid,bc,q,dt,cfl,chemistry,chemistry_rtol,chemistry_atol_species, &
      chemistry_atol_temperature,chemistry_max_steps,change,rejected)
    call require(time+dt>time,'2D time increment below machine resolution')
    time=time+dt;step=step+1;balance=balance+change;total_rejected=total_rejected+rejected
    if(mod(step,write_every)==0.or.time>=end_time) call snapshot()
  end do
  rate=-total0-balance
  do c=1,nc
    rate=rate+grid%mesh%volume(c)*q(:,c)
  end do
  write(out,'(a,es24.16)') '# mass_balance_residual=',sum(rate(:ns))
  write(out,'(a,2(es24.16,1x))') '# momentum_balance_residual=',rate(ns+1:ns+2)
  write(out,'(a,es24.16)') '# energy_balance_residual=',rate(nv)
  write(out,'(a,i0)') '# rejected_steps=',total_rejected
  write(out,'(a)') '# SUCCESS'
  close(out,iostat=ios);call require(ios==0,'Failed closing flow2d output')
  write(*,'(a)') '[OK] Final flow2d snapshot written; calculation finished normally.'
contains
  subroutine validate_input_y(values)
    real(dp), intent(in) :: values(:)
    call require(all(ieee_is_finite(values)).and.all(values>=0),'Missing/invalid mass fractions')
    call require(abs(sum(values)-1)<1.e-12_dp,'Mass fractions must sum to one; no normalization')
  end subroutine
  subroutine snapshot()
    integer :: cell
    do cell=1,nc
      call conserved_nd(m,q(:,cell),rho,v,t,p,sound,y)
      write(out,'(i0,",",es24.16,",",i0,*(",",es24.16))') step,time,cell,grid%cell_center(:,cell), &
        grid%mesh%volume(cell),rho,v,t,p,y
    end do
    flush(out)
  end subroutine
end program
