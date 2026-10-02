program rf_legacy_unit
  use mod_rf_flow1d
  use mod_rf_thermo
  use mod_mc_state_layout, only: mc_state_layout,initialize_mc_state_layout
  use mod_mc_euler_config, only: mc_euler_config,initialize_mc_euler_config
  use mod_mc_thermodynamics_provider, only: configure_mc_thermodynamics
  use mod_mc_chemistry_provider, only: configure_mc_chemistry
  use mod_mc_chemistry_integrator, only: advance_mc_chemistry_interval
  use mod_mc_euler_flux, only: advance_mc_euler_ssprk3
  implicit none
  integer, parameter :: nx=24
  type(rf_mechanism) :: m
  type(mc_state_layout) :: layout
  type(mc_euler_config) :: config
  real(dp) :: q(4,nx),old(nx,1,1,6),q0(nx,1,1,6),rhs(nx,1,1,6)
  real(dp) :: bc(4),y(2),dt,err,x,pressure,temperature,velocity,initial(4),netbc(4)
  integer :: i,step,mode,wall
  character(2048) :: path
  character(16) :: boundary
  logical :: chemistry
  call require(command_argument_count()==1,'Usage: rf_legacy_unit legacy_one_step.in')
  call get_command_argument(1,path)
  call initialize_mc_state_layout(layout,2)
  call initialize_mc_euler_config(config,2)
  call configure_mc_thermodynamics(trim(path),2,['A','B'])
  call configure_mc_chemistry(trim(path),2,['A','B'])
  config%nx=nx;config%ny=1;config%nz=1
  m%ne=1;allocate(m%species(2),m%reactions(1))
  do i=1,2
    m%species(i)%mass=.01_dp;m%species(i)%pref=101325;m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp];m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp;m%species(i)%atoms=1
  end do
  m%species(2)%coeff(6,1)=-100
  associate(a=>m%reactions(1))
    a%kind=1;a%reversible=0;a%high=[log(1000._dp),0._dp,0._dp,1._dp]
    a%reactants=[1._dp,0._dp];a%products=[0._dp,1._dp];a%orders=a%reactants
  end associate
  do mode=1,2
    chemistry=mode==2
    do wall=1,2
      config%boundary_face_types='periodic';boundary='periodic'
      if(wall==2) then
        config%boundary_face_types(1:2)='reflective';boundary='reflecting'
      end if
      do i=1,nx
        x=(i-.5_dp)/nx
        y=[.5_dp+.2_dp*sin(2*acos(-1._dp)*x),.5_dp-.2_dp*sin(2*acos(-1._dp)*x)]
        pressure=101325;temperature=1100;velocity=20
        if(wall==2) then
          pressure=merge(202650._dp,101325._dp,i<=nx/2);velocity=0
        end if
        call primitive_to_conserved(m,temperature,pressure,velocity,y,q(:,i))
        old(i,1,1,:)=0;old(i,1,1,1:3)=q(1:3,i);old(i,1,1,6)=q(4,i)
      end do
      initial=sum(q,dim=2)/nx;netbc=0
      do step=1,50
        dt=1.e-7_dp
        if(chemistry) then
          do i=1,nx
            call advance_mc_chemistry_interval(old(i,1,1,:),dt/2,layout,config%gamma,.1_dp,10000)
          end do
        end if
        call advance_mc_euler_ssprk3(old,q0,rhs,dt,layout,config)
        if(chemistry) then
          do i=1,nx
            call advance_mc_chemistry_interval(old(i,1,1,:),dt/2,layout,config%gamma,.1_dp,10000)
          end do
        end if
        call advance_flow(m,q,1._dp/nx,dt,.4_dp,boundary,boundary,chemistry, &
          1.e-10_dp,1.e-17_dp,1.e-9_dp,10000,bc,reconstruction='first_order')
        call require(dt==1.e-7_dp,'Legacy comparison requires the same accepted timestep')
        netbc=netbc+bc
      end do
      err=0
      do i=1,nx
        err=max(err,maxval(abs(q(1:3,i)-old(i,1,1,1:3))/max(1._dp,abs(q(1:3,i)))))
        err=max(err,abs(q(4,i)-old(i,1,1,6))/max(1._dp,abs(q(4,i))))
      end do
      write(*,'(a,l1,a,a,a,es16.8)') 'Legacy CFD chemistry=',chemistry,' boundary=',trim(boundary),' error=',err
      call require(err<1.e-7_dp,'Legacy/new CFD fields must agree at matched EOS, flux, dt and chemistry')
      call require(abs(sum(q(:2,:))/nx-sum(initial(:2))-sum(netbc(:2)))<1.e-10_dp,'Migration mass conservation')
      call require(abs(sum(q(4,:))/nx-initial(4)-netbc(4))/abs(initial(4))<1.e-10_dp,'Migration energy conservation')
    end do
  end do
  write(*,'(a)') '[OK] actual legacy Euler/one-step chemistry vs independent ReactingFlow'
end program
