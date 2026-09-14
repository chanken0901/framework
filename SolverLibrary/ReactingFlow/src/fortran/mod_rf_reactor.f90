module mod_rf_reactor
  use mod_rf_kinetics
  use dvode_module, only: dvode_t
  implicit none
  type, extends(dvode_t) :: reactor_context
    type(rf_mechanism), pointer :: mechanism=>null()
    real(dp) :: rho0,p0
    logical :: constant_pressure
    integer :: negative_trials=0
    integer :: dependent_species=1
  end type
contains
  subroutine reactor_rhs(me,neq,time,state,derivative)
    class(dvode_t), intent(inout) :: me
    integer :: neq
    real(dp) :: time,state(neq),derivative(neq)
    real(dp) :: y(neq),cp,cv,h,e,r,rho,cpi,hi,si,heat,capacity
    integer :: i,j
    select type(me)
    type is(reactor_context)
      associate(m=>me%mechanism)
        block
          real(dp) :: qf(size(m%reactions)),qr(size(m%reactions)),net(size(m%reactions)),omega(neq)
          call require(all(ieee_is_finite(state)),'Nonfinite Newton trial')
          call unpack(state,me%dependent_species,y)
          if(any(y<0)) me%negative_trials=me%negative_trials+1
          ! Numerical extension for Newton trials only. Never change the solution vector.
          y=max(y,0._dp)
          call require(sum(y)>0,'No positive species in Newton trial')
          y=y/sum(y)
          call mixture(m,state(1),y,me%p0,cp,cv,h,e,r)
          rho=me%rho0; capacity=cv
          if(me%constant_pressure) then
            rho=me%p0/(r*state(1)); capacity=cp
          end if
          call rates(m,state(1),rho,y,qf,qr,net,omega,heat)
          derivative(1)=0; j=1
          do i=1,neq
            if(i/=me%dependent_species) then
              j=j+1
              derivative(j)=omega(i)/rho
            end if
            call species_thermo(m%species(i),state(1),cpi,hi,si)
            if(.not.me%constant_pressure) hi=hi-gas_r*state(1)
            derivative(1)=derivative(1)-hi/m%species(i)%mass*omega(i)/rho/capacity
          end do
        end block
      end associate
    class default
      call require(.false.,'Invalid reactor context')
    end select
  end subroutine

  subroutine unpack(state,dependent,y)
    real(dp), intent(in) :: state(:)
    integer, intent(in) :: dependent
    real(dp), intent(out) :: y(:)
    integer :: i,j
    j=1; y=0
    do i=1,size(y)
      if(i==dependent) cycle
      j=j+1; y(i)=state(j)
    end do
    y(dependent)=1-sum(y)
  end subroutine

  subroutine run_reactor(m,t0,p0,y0,tend,constant_pressure,rtol,atoly,atolt,maxstep,maxsteps,rise,unit)
    type(rf_mechanism), intent(in), target :: m
    real(dp), intent(in) :: t0,p0,y0(:),tend,rtol,atoly,atolt,maxstep,rise
    logical, intent(in) :: constant_pressure
    integer, intent(in) :: unit,maxsteps
    type(reactor_context) :: solver
    real(dp), allocatable :: state(:),atol(:),rw(:),interp(:)
    integer, allocatable :: iw(:)
    real(dp) :: time,previous,cp,cv,h,e,r,energy0,scale,err,elemerr,masserr,rho,p,ignition,lo,hi,mid
    real(dp) :: elem0(m%ne),elem(m%ne),y(size(y0))
    integer :: n,i,j,istate,step,iflag
    call require(all(ieee_is_finite([t0,p0,tend,rtol,atoly,atolt,maxstep,rise])),'Nonfinite reactor input')
    call require(min(t0,p0,tend,atoly,atolt,maxstep,rise)>0,'Reactor inputs must be positive')
    call require(rtol>=1.e-12_dp.and.rtol<=1.e-2_dp.and.maxsteps>0,'Invalid tolerances/step limit')
    call mixture(m,t0,y0,p0,cp,cv,h,e,r)
    solver%mechanism=>m; solver%rho0=p0/(r*t0); solver%p0=p0; solver%constant_pressure=constant_pressure
    energy0=merge(h,e,constant_pressure); scale=max(1._dp,abs(energy0),cp*t0)
    elem0=0
    do i=1,size(y0)
      elem0=elem0+y0(i)/m%species(i)%mass*m%species(i)%atoms
    end do
    n=size(y0)
    solver%dependent_species=maxloc(y0,dim=1)
    allocate(state(n),interp(n),atol(n),rw(22+9*n+2*n*n),iw(30+n))
    state(1)=t0; j=1
    do i=1,n
      if(i==solver%dependent_species) cycle
      j=j+1; state(j)=y0(i)
    end do
    atol=atoly; atol(1)=atolt
    rw=0; iw=0; rw(1)=tend; rw(6)=maxstep; iw(6)=maxsteps
    call solver%initialize(f=reactor_rhs)
    time=0; istate=1; step=0; ignition=-1; err=0; elemerr=0; masserr=0
    write(unit,'(a)') '# source_sha256='//m%source_hash
    write(unit,'(a)') '# canonical_sha256='//m%canonical_hash
    write(unit,'(a,l1)') '# constant_pressure=',constant_pressure
    write(unit,'(a,es25.16e3)') '# rtol=',rtol
    write(unit,'(a,es25.16e3)') '# atol_species=',atoly
    write(unit,'(a,es25.16e3)') '# atol_temperature=',atolt
    write(unit,'(a,es25.16e3)') '# max_step=',maxstep
    write(unit,'(a,es25.16e3)') '# ignition_temperature=',t0+rise
    write(unit,'(a)',advance='no') 'time,temperature,pressure,density'
    do i=1,size(y0)
      write(unit,'(a)',advance='no') ',Y_'//trim(m%species(i)%name)
    end do
    write(unit,*)
    do
      call unpack(state,solver%dependent_species,y)
      call mixture(m,state(1),y,p0,cp,cv,h,e,r)
      rho=solver%rho0; p=rho*r*state(1)
      if(constant_pressure) then
        p=p0; rho=p/(r*state(1))
      end if
      err=max(err,abs(merge(h,e,constant_pressure)-energy0)/scale)
      masserr=max(masserr,abs(sum(y)-1))
      elem=0
      do i=1,size(y0)
        elem=elem+y(i)/m%species(i)%mass*m%species(i)%atoms
      end do
      elemerr=max(elemerr,maxval(abs(elem-elem0)/max(1._dp,abs(elem0))))
      call require(err<=max(100*rtol,1.e-7_dp).and.elemerr<=1.e-8_dp,'Reactor conservation check failed')
      write(unit,'(*(es25.16e3,:,","))') time,state(1),p,rho,y
      if(time>=tend) exit
      call require(step<maxsteps,'Reactor exceeded max_steps (partial output is not success)')
      previous=time
      ! ITASK=5: one internal step, never beyond TCRIT=tend. MF=22: BDF, numerical dense Jacobian.
      call solver%solve(n,state,time,tend,2,[rtol],atol,5,istate,1,rw,size(rw),iw,size(iw),22)
      call require(istate>=0,'DVODE failed (partial output is not success)')
      call unpack(state,solver%dependent_species,y)
      call check_y(m,y)
      if(ignition<0.and.state(1)>=t0+rise) then
        lo=previous; hi=time
        do j=1,60
          mid=(lo+hi)/2
          call solver%dvindy(mid,0,rw(21:),n,interp,iflag)
          call require(iflag==0,'DVODE ignition interpolation failed')
          if(interp(1)>=t0+rise) then
            hi=mid
          else
            lo=mid
          end if
        end do
        ignition=(lo+hi)/2
      end if
      step=step+1
    end do
    write(unit,'(a,es25.16e3)') '# ignition_delay_s=',ignition
    write(unit,'(a,es25.16e3)') '# mass_sum_error=',masserr
    write(unit,'(a,es25.16e3)') '# element_relative_error=',elemerr
    write(unit,'(a,es25.16e3)') '# energy_relative_error=',err
    write(unit,'(a,i0)') '# negative_trial_rhs_calls=',solver%negative_trials
    write(unit,'(a,i0)') '# accepted_steps=',step
    write(unit,'(a,i0)') '# rhs_calls=',iw(12)
    write(unit,'(a,i0)') '# jacobian_calls=',iw(13)
    write(unit,'(a)') '# SUCCESS'
  end subroutine
end module
