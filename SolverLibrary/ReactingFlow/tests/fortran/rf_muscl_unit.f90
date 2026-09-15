program rf_muscl_unit
  use mod_rf_flow1d
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  real(dp) :: q(5,8),qm(5,8),qp(5,8),rho,u,t,p,a,y(3),err32,err64,err128,errfirst
  real(dp) :: change(5),old(5,8),dt
  integer :: i,j,k
  m%ne=1; allocate(m%species(3),m%reactions(0))
  do i=1,3
    m%species(i)%mass=.01_dp; m%species(i)%pref=101325; m%species(i)%model=7
    allocate(m%species(i)%bounds(2),m%species(i)%coeff(9,1),m%species(i)%atoms(1))
    m%species(i)%bounds=[200._dp,4000._dp]; m%species(i)%coeff=0
    m%species(i)%coeff(1,1)=3.5_dp; m%species(i)%atoms=1
  end do
  do i=1,8
    y=[real(i-1,dp)/7,real(8-i,dp)/7,0._dp]
    call primitive_to_conserved(m,800._dp+50*i,101325._dp+100*i,20._dp,y,q(:,i))
  end do
  call reconstruct_faces(m,q,'periodic','periodic',qm,qp,'first_order')
  call require(all(qm==q).and.all(qp==q),'First order face identity')
  call reconstruct_faces(m,q,'periodic','periodic',qm,qp,'muscl')
  do i=1,8
    do j=1,2
      if(j==1) then
        call conserved_to_primitive(m,qm(:,i),rho,u,t,p,a,y)
      else
        call conserved_to_primitive(m,qp(:,i),rho,u,t,p,a,y)
      end if
      call require(minval(y)>=0.and.maxval(y)<=1.and.abs(sum(y)-1)<1.e-14_dp,'Bounded face composition')
      call require(y(3)==0,'Absent species remains absent')
      call require(t>=850-1.e-7_dp.and.t<=1200+1.e-7_dp.and.p>0,'Bounded face temperature/pressure')
    end do
  end do
  old=q
  do k=1,2
    q=old
    if(k==1) then
      dt=flow_timestep(m,q,.1_dp,.3_dp,reconstruction='muscl',left_bc='periodic',right_bc='periodic')/2
      call advance_flow(m,q,.1_dp,dt,.3_dp,'periodic','periodic',.false.,1.e-9_dp,1.e-16_dp, &
        1.e-8_dp,10000,change,reconstruction='muscl')
    else
      dt=flow_timestep(m,q,.1_dp,.3_dp,reconstruction='muscl',left_bc='reflecting',right_bc='reflecting')/2
      call advance_flow(m,q,.1_dp,dt,.3_dp,'reflecting','reflecting',.false.,1.e-9_dp,1.e-16_dp, &
        1.e-8_dp,10000,change,reconstruction='muscl')
      call require(maxval(abs(change(:3)))<1.e-12_dp.and.abs(change(5))<1.e-12_dp,'Impermeable adiabatic wall')
    end if
    call require(maxval(abs(sum(q-old,dim=2)*.1_dp-change)/max(1._dp,abs(sum(old,dim=2)*.1_dp))) &
      <1.e-12_dp,'MUSCL boundary-corrected conservation')
  end do
  err32=advection_error(32,'muscl'); err64=advection_error(64,'muscl'); err128=advection_error(128,'muscl')
  errfirst=advection_error(64,'first_order')
  write(*,'(a,4es15.6)') 'Advection L1 errors: MUSCL 32/64/128, first-order 64: ',err32,err64,err128,errfirst
  call require(err64<.32_dp*err32.and.err128<.32_dp*err64,'Second-order smooth advection convergence')
  call require(err64<.2_dp*errfirst,'MUSCL reduces smooth advection error')
  err32=advection_error(32,'muscl',.true.); err64=advection_error(64,'muscl',.true.)
  err128=advection_error(128,'muscl',.true.)
  write(*,'(a,3es15.6)') 'Density wave relative L1 errors: MUSCL 32/64/128: ',err32,err64,err128
  call require(err64<.32_dp*err32.and.err128<.32_dp*err64,'Second-order density/energy advection convergence')
  write(*,'(a)') '[OK] MUSCL composition bounds, wall/periodic conservation, second-order advection'
contains
  real(dp) function advection_error(nx,method,density_wave) result(err)
    integer, intent(in) :: nx
    character(*), intent(in) :: method
    logical, intent(in), optional :: density_wave
    logical :: density_test
    real(dp) :: state(5,nx),start(5),bc(5),yy(3),dx,time,step_dt,x,pi,amp,exact,rr,uu,tt,pp,aa
    integer :: cell
    density_test=.false.
    if(present(density_wave)) density_test=density_wave
    pi=acos(-1._dp); dx=1._dp/nx; amp=.2_dp*sin(pi*dx)/(pi*dx)
    do cell=1,nx
      x=(cell-.5_dp)*dx
      yy(1)=.4_dp+amp*sin(2*pi*x); yy(2)=.5_dp-amp*sin(2*pi*x); yy(3)=.1_dp
      tt=1100
      if(density_test) then
        yy=[.4_dp,.5_dp,.1_dp]
        rr=.1_dp*(1+amp*sin(2*pi*x))
        tt=101325/(rr*gas_r/.01_dp)
      end if
      call primitive_to_conserved(m,tt,101325._dp,1000._dp,yy,state(:,cell))
    end do
    start=sum(state,dim=2); time=0
    do while(time<.001_dp)
      step_dt=min(.001_dp-time,flow_timestep(m,state,dx,.3_dp, &
        reconstruction=method,left_bc='periodic',right_bc='periodic'))
      call advance_flow(m,state,dx,step_dt,.3_dp,'periodic','periodic',.false.,1.e-9_dp,1.e-16_dp, &
        1.e-8_dp,10000,bc,reconstruction=method)
      time=time+step_dt
    end do
    call require(maxval(abs(sum(state,dim=2)-start)/max(1._dp,abs(start)))<1.e-11_dp,'Advection conservation')
    err=0
    do cell=1,nx
      call conserved_to_primitive(m,state(:,cell),rr,uu,tt,pp,aa,yy)
      exact=.4_dp+amp*sin(2*pi*((cell-.5_dp)*dx-1000*time))
      if(density_test) then
        exact=1+amp*sin(2*pi*((cell-.5_dp)*dx-1000*time))
        err=err+abs(rr/.1_dp-exact)/nx
      else
        err=err+abs(yy(1)-exact)/nx
      end if
    end do
  end function
end program
