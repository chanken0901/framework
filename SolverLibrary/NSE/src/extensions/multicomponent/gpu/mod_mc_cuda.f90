module mod_mc_cuda
  use iso_c_binding
  use mod_precision, only: dp
  use mod_mc_config, only: mc_config
  use mod_mc_state_layout, only: mc_state_layout
  use mod_mc_euler_config, only: mc_euler_config
  use mod_mc_reactive_config, only: mc_reactive_config
  use mod_mc_thermodynamics_provider, only: mc_export_thermodynamics
  use mod_mc_transport_provider, only: mc_export_transport
  use mod_mc_chemistry_provider, only: mc_export_chemistry
  use mod_mc_euler_field, only: set_mc_euler_conservative_state,initialize_mc_euler_state
  use mod_mc_euler_solver, only: write_mc_euler_csv
  use mod_mc_reactive_solver, only: write_reactive_snapshot,write_reactive_history_header,write_reactive_history_row
  implicit none
  private
  public::mc_gpu_create,mc_gpu_destroy,mc_gpu_upload,mc_gpu_download,mc_gpu_timestep
  public::mc_gpu_chemistry,mc_gpu_begin,mc_gpu_stage,mc_gpu_check,run_mc_cuda
  public::mc_gpu_halo,mc_gpu_boundary,mc_gpu_failure
  public::mc_gpu_statistics,mc_gpu_transfer_counts,mc_gpu_check_totals
  abstract interface
    subroutine failure_interface()
    end subroutine
  end interface
  procedure(failure_interface),pointer::mc_gpu_failure=>null()
  type,bind(c)::device_model
    integer(c_int)::n(3),ns,mapped,bc(6)
    real(c_double)::origin(3),h(3),xrange(2),nozzle(4),lengths(3),relaxation,length
    real(c_double)::cfl,diffusion_cfl,chemistry_cfl
    real(c_double)::thermo(4),species(16,64),transport(2),diff(64),chem(3),reaction(3,64),reference(68,6)
  end type
  interface
    integer(c_int) function statistics_c(handle,stats) bind(c,name='mc_cuda_statistics')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      real(c_double),intent(out)::stats(*)
    end function
    integer(c_int) function counts_c(handle,counts) bind(c,name='mc_cuda_transfer_counts')
      import c_ptr,c_int
      type(c_ptr),value::handle
      integer(c_int),intent(out)::counts(2)
    end function
    integer(c_int) function boundary_c(handle) bind(c,name='mc_cuda_local_boundary')
      import c_ptr,c_int
      type(c_ptr),value::handle
    end function
    integer(c_int) function halo_c(handle,d,side,unpack,data) bind(c,name='mc_cuda_halo')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      integer(c_int),value::d,side,unpack
      real(c_double)::data(*)
    end function
    integer(c_int) function create_c(m,device,handle) bind(c,name='mc_cuda_create')
      import device_model,c_int,c_ptr
      type(device_model),intent(in)::m
      integer(c_int),value::device
      type(c_ptr),intent(out)::handle
    end function
    subroutine mc_gpu_destroy(handle) bind(c,name='mc_cuda_destroy')
      import c_ptr
      type(c_ptr),value::handle
    end subroutine
    integer(c_int) function upload_c(handle,q) bind(c,name='mc_cuda_upload')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      real(c_double),intent(in)::q(*)
    end function
    integer(c_int) function download_c(handle,q,rhs) bind(c,name='mc_cuda_download')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      real(c_double),intent(out)::q(*)
      integer(c_int),value::rhs
    end function
    integer(c_int) function dt_c(handle,chem,dt) bind(c,name='mc_cuda_timestep')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      integer(c_int),value::chem
      real(c_double),intent(out)::dt
    end function
    integer(c_int) function chemistry_c(handle,duration,limit,steps) bind(c,name='mc_cuda_chemistry')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      real(c_double),value::duration
      integer(c_int),value::limit
      integer(c_int),intent(out)::steps
    end function
    integer(c_int) function begin_c(handle) bind(c,name='mc_cuda_begin')
      import c_ptr,c_int
      type(c_ptr),value::handle
    end function
    integer(c_int) function stage_c(handle,dt,stage) bind(c,name='mc_cuda_stage')
      import c_ptr,c_int,c_double
      type(c_ptr),value::handle
      real(c_double),value::dt
      integer(c_int),value::stage
    end function
  end interface
contains
  subroutine mc_gpu_statistics(h,stats)
    type(c_ptr),intent(in)::h
    real(dp),contiguous,intent(out)::stats(:)
    call mc_gpu_check(statistics_c(h,stats))
  end subroutine
  subroutine mc_gpu_transfer_counts(h,counts)
    type(c_ptr),intent(in)::h
    integer,intent(out)::counts(2)
    integer(c_int)::values(2)
    call mc_gpu_check(counts_c(h,values));counts=values
  end subroutine
  subroutine mc_gpu_check_totals(initial,final,l,c)
    real(dp),intent(in)::initial(:),final(:)
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    real(dp)::error
    if(all(c%boundary_face_types=='periodic'))then
      error=abs(sum(final(1:l%nspecies)-initial(1:l%nspecies)))/max(1.0_dp,abs(sum(initial(1:l%nspecies))))
      error=max(error,maxval(abs(final(l%nspecies+1:l%nvariables)-initial(l%nspecies+1:l%nvariables))/ &
        max(1.0_dp,abs(initial(l%nspecies+1:l%nvariables)))))
      if(error>1e-10_dp)then
        print *, 'ERROR: CUDA periodic conservation check failed: ',error
        call mc_gpu_check(1)
      end if
    end if
    if(final(l%nvariables+1)<-1e-12_dp.or.any(final(l%nvariables+2:l%nvariables+4)<=0))call mc_gpu_check(1)
  end subroutine
  subroutine mc_gpu_check(status)
    integer,intent(in)::status
    if(status/=0)then
      if(associated(mc_gpu_failure))call mc_gpu_failure()
      error stop 'multicomponent CUDA backend failure (see preceding diagnostic)'
    end if
  end subroutine
  subroutine mc_gpu_boundary(h)
    type(c_ptr),intent(in)::h
    call mc_gpu_check(boundary_c(h))
  end subroutine
  subroutine mc_gpu_halo(h,d,side,unpack,data)
    type(c_ptr),intent(in)::h
    integer,intent(in)::d,side
    logical,intent(in)::unpack
    real(dp),contiguous,intent(inout)::data(:)
    call mc_gpu_check(halo_c(h,int(d-1,c_int),int(side,c_int),merge(1_c_int,0_c_int,unpack),data))
  end subroutine
  subroutine mc_gpu_create(handle,l,c,r,device,remote)
    type(c_ptr),intent(out)::handle
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    type(mc_reactive_config),intent(in)::r
    integer,intent(in)::device
    logical,intent(in),optional::remote(6)
    type(device_model)::m
    integer::f
    m%n=[c%nx,c%ny,c%nz];m%ns=l%nspecies;m%mapped=0
    if(c%geometry=='planar_nozzle')m%mapped=1
    if(c%geometry/='cartesian'.and.c%geometry/='planar_nozzle') error stop 'unsupported CUDA geometry'
    m%origin=[c%x_min,c%y_min,c%z_min]
    m%lengths=[c%x_max-c%x_min,c%y_max-c%y_min,c%z_max-c%z_min]
    m%h=m%lengths/real(m%n,dp)
    where(c%global_boundary_lengths>0)m%lengths=c%global_boundary_lengths
    m%xrange=[c%x_min,c%x_max]
    m%nozzle=[c%nozzle_inlet_half_height,c%nozzle_throat_half_height,c%nozzle_exit_half_height,c%nozzle_throat_x]
    m%relaxation=c%boundary_relaxation_strength;m%length=c%boundary_length_scale
    m%cfl=c%cfl;m%diffusion_cfl=c%diffusion_cfl;m%chemistry_cfl=r%chemistry_cfl
    m%species=0;m%diff=0;m%reaction=0;m%reference=0
    call mc_export_thermodynamics(l%nspecies,m%thermo,m%species(:,1:l%nspecies))
    call mc_export_transport(l%nspecies,m%transport,m%diff(1:l%nspecies))
    call mc_export_chemistry(l%nspecies,m%chem,m%reaction(:,1:l%nspecies))
    do f=1,6
      if(present(remote))then
        if(remote(f))then
          m%bc(f)=4
          cycle
        end if
      end if
      select case(c%boundary_face_types(f))
      case('periodic');m%bc(f)=0
      case('reflective');m%bc(f)=1
      case('dirichlet');m%bc(f)=2
      case('non_reflecting');m%bc(f)=3
      case default;error stop 'unsupported CUDA boundary'
      end select
      if(m%bc(f)<2)cycle
      call set_mc_euler_conservative_state(m%reference(1:l%nvariables,f),l,c%gamma, &
        c%boundary_reference_densities(f),c%boundary_reference_velocities(:,f), &
        c%boundary_reference_pressures(f),c%boundary_reference_mass_fractions(f,1:l%nspecies))
    end do
    call mc_gpu_check(create_c(m,int(device,c_int),handle))
  end subroutine
  subroutine mc_gpu_upload(h,q)
    type(c_ptr),intent(in)::h
    real(dp),contiguous,intent(in)::q(:,:,:,:)
    call mc_gpu_check(upload_c(h,q))
  end subroutine
  subroutine mc_gpu_download(h,q,rhs)
    type(c_ptr),intent(in)::h
    real(dp),contiguous,intent(out)::q(:,:,:,:)
    logical,intent(in),optional::rhs
    integer(c_int)::mode
    mode=0
    if(present(rhs))then
      if(rhs)mode=1
    end if
    call mc_gpu_check(download_c(h,q,mode))
  end subroutine
  real(dp) function mc_gpu_timestep(h,chem,fixed) result(dt)
    type(c_ptr),intent(in)::h
    logical,intent(in)::chem
    logical,intent(in),optional::fixed
    integer(c_int)::mode
    mode=merge(1_c_int,0_c_int,chem)
    if(present(fixed))then
      if(fixed)mode=2
    end if
    call mc_gpu_check(dt_c(h,mode,dt))
  end function
  subroutine mc_gpu_chemistry(h,duration,limit,steps)
    type(c_ptr),intent(in)::h
    real(dp),intent(in)::duration
    integer,intent(in)::limit
    integer,intent(out)::steps
    integer(c_int)::used
    call mc_gpu_check(chemistry_c(h,duration,int(limit,c_int),used))
    steps=used
  end subroutine
  subroutine mc_gpu_begin(h)
    type(c_ptr),intent(in)::h
    call mc_gpu_check(begin_c(h))
  end subroutine
  subroutine mc_gpu_stage(h,dt,stage)
    type(c_ptr),intent(in)::h
    real(dp),intent(in)::dt
    integer,intent(in)::stage
    call mc_gpu_check(stage_c(h,dt,int(stage,c_int)))
  end subroutine
  subroutine run_mc_cuda(model,l,c,r)
    type(mc_config),intent(in)::model
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::c
    type(mc_reactive_config),intent(in)::r
    type(c_ptr)::gpu
    real(dp),allocatable::q(:,:,:,:)
    real(dp)::time,dt,stable,started,elapsed
    real(dp)::initial(l%nvariables+4),final(l%nvariables+4)
    integer::step,stage,substeps,maximum_steps,unit,ios,count,rate,device
    character(len=32)::device_text
    logical::output
    device=0
    call get_environment_variable('NSE_MC_CUDA_DEVICE',device_text,status=ios)
    if(ios==0)then
      read(device_text,*,iostat=ios)device
      if(ios/=0.or.device<0)error stop 'invalid NSE_MC_CUDA_DEVICE'
    end if
    allocate(q(c%nx,c%ny,c%nz,l%nvariables))
    call initialize_mc_euler_state(q,l,c)
    call mc_gpu_create(gpu,l,c,r,device)
    call mc_gpu_upload(gpu,q)
    call mc_gpu_statistics(gpu,initial)
    time=0;maximum_steps=0
    if(r%write_history)then
      open(newunit=unit,file=trim(r%history_file),status='replace',iostat=ios)
      if(ios/=0)error stop 'cannot open CUDA reactive history'
      call write_reactive_history_header(unit,model,l)
      call write_reactive_history_row(unit,0,time,0.0_dp,q,l,c)
    end if
    if(r%write_snapshots)call write_reactive_snapshot(0,q,model,l,c,r%snapshot_prefix)
    call system_clock(count,rate);started=real(count,dp)/rate
    write(*,'(A)') 'Resident CUDA reactive multicomponent solver: q/RK/flux/gradients/chemistry on GPU'
    write(*,'(A)') '# step time dt (full-state downloads only for output)'
    do step=1,c%nsteps
      stable=mc_gpu_timestep(gpu,c%dt<=0,c%dt>0)
      dt=stable
      if(c%dt>0)then
        if(c%dt>stable*(1+100*epsilon(dt)))error stop 'fixed CUDA dt violates fluid CFL/diffusion limit'
        dt=c%dt
      end if
      call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,substeps)
      maximum_steps=max(maximum_steps,substeps)
      stable=mc_gpu_timestep(gpu,.false.,.true.)
      if(dt>stable*(1+100*epsilon(dt)))error stop 'CUDA dt violates post-chemistry fluid limit'
      call mc_gpu_begin(gpu)
      do stage=1,3
        call mc_gpu_stage(gpu,dt,stage)
      end do
      call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,substeps)
      maximum_steps=max(maximum_steps,substeps);time=time+dt
      output=mod(step,r%output_every)==0.or.step==c%nsteps
      if(output.and.(r%write_history.or.r%write_snapshots))then
        call mc_gpu_download(gpu,q)
        if(r%write_history)call write_reactive_history_row(unit,step,time,dt,q,l,c)
        if(r%write_snapshots)call write_reactive_snapshot(step,q,model,l,c,r%snapshot_prefix)
      end if
      if(output)write(*,'(I10,2ES24.16)')step,time,dt
    end do
    call mc_gpu_statistics(gpu,final)
    call mc_gpu_check_totals(initial,final,l,c)
    if(c%write_final)then
      call mc_gpu_download(gpu,q)
      call write_mc_euler_csv(trim(c%output_file),q,model,l,c)
    end if
    if(r%write_history)close(unit)
    call mc_gpu_destroy(gpu)
    call system_clock(count);elapsed=real(count,dp)/rate-started
    write(*,'(A,I0)') 'maximum chemistry substeps per half-step = ',maximum_steps
    write(*,'(A,ES14.6)') 'elapsed wall seconds (including output) = ',elapsed
    write(*,'(A)') 'Reactive multicomponent Navier-Stokes calculation completed successfully.'
  end subroutine
end module
