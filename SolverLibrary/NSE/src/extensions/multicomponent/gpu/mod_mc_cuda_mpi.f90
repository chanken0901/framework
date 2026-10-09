module mod_mc_cuda_mpi
  use iso_c_binding,only:c_ptr
  use mod_precision,only:dp
  use mod_mc_cuda
  use mod_mc_config,only:mc_config
  use mod_mc_state_layout,only:mc_state_layout
  use mod_mc_euler_config,only:mc_euler_config
  use mod_mc_reactive_config,only:mc_reactive_config
  use mod_mc_euler_field,only:initialize_mc_euler_state
  use mod_mc_euler_solver,only:write_mc_euler_csv
  use mod_mc_reactive_solver,only:write_reactive_snapshot,write_reactive_history_header,write_reactive_history_row
  use mod_mc_mpi_pencil,only:mc_pencil_domain,initialize_mc_pencil,gather_mc_pencil,mc_rank
  implicit none
  include 'mpif.h'
  private
  public::run_mc_cuda_mpi,mc_cuda_exchange,mc_cuda_remote,mc_cuda_abort
contains
  subroutine mc_cuda_abort()
    integer::ierr
    call MPI_Abort(MPI_COMM_WORLD,10,ierr)
  end subroutine
  function mc_cuda_remote(domain) result(remote)
    type(mc_pencil_domain),intent(in)::domain
    logical::remote(6)
    remote=.false.
    remote(3)=domain%lower(1)>0;remote(4)=domain%upper(1)>0
    remote(5)=domain%lower(2)>0;remote(6)=domain%upper(2)>0
  end function
  subroutine mc_cuda_exchange(gpu,domain,send,recv)
    type(c_ptr),intent(in)::gpu
    type(mc_pencil_domain),intent(in)::domain
    real(dp),intent(inout),contiguous::send(:),recv(:)
    integer::d,count,n(3),status(MPI_STATUS_SIZE),ierr
    call mc_gpu_boundary(gpu)
    n=[domain%local%nx,domain%local%ny,domain%local%nz]+4
    do d=1,2
      if(domain%dims(d)==1)cycle
      ! Buffer size is supplied per variable by the caller via total capacity;
      ! all padded transverse cells (including corners) are sent, y then z.
      count=size(send)/(2*n(1)*max(n(2),n(3)))
      count=count*2*n(1)*n(4-d)
      call mc_gpu_halo(gpu,d+1,0,.false.,send(1:count))
      call MPI_Sendrecv(send,count,MPI_DOUBLE_PRECISION,domain%neighbor_minus(d),810+d, &
        recv,count,MPI_DOUBLE_PRECISION,domain%neighbor_plus(d),810+d,domain%comm,status,ierr)
      if(ierr/=MPI_SUCCESS)call mc_cuda_abort()
      if(domain%upper(d)>0)call mc_gpu_halo(gpu,d+1,1,.true.,recv(1:count))
      call mc_gpu_halo(gpu,d+1,1,.false.,send(1:count))
      call MPI_Sendrecv(send,count,MPI_DOUBLE_PRECISION,domain%neighbor_plus(d),820+d, &
        recv,count,MPI_DOUBLE_PRECISION,domain%neighbor_minus(d),820+d,domain%comm,status,ierr)
      if(ierr/=MPI_SUCCESS)call mc_cuda_abort()
      if(domain%lower(d)>0)call mc_gpu_halo(gpu,d+1,0,.true.,recv(1:count))
    end do
  end subroutine
  subroutine run_mc_cuda_mpi(model,l,global,r,path)
    type(mc_config),intent(in)::model
    type(mc_state_layout),intent(in)::l
    type(mc_euler_config),intent(in)::global
    type(mc_reactive_config),intent(in)::r
    character(len=*),intent(in)::path
    type(mc_pencil_domain)::domain
    type(c_ptr)::gpu
    real(dp),allocatable::q(:,:,:,:),full(:,:,:,:),send(:),recv(:)
    real(dp)::time,dt,stable,local_dt,started,elapsed,total_elapsed
    real(dp)::local_stats(l%nvariables+4),initial(l%nvariables+4),final(l%nvariables+4)
    integer::unit,ios,ierr,grid(2),shared,local_rank,device,n(3),count
    integer::step,stage,substeps,maximum_steps,global_steps,history
    character(len=32)::decomposition,device_text
    character(len=512)::line
    logical::output
    integer::process_grid(2)
    namelist /multicomponent_parallel/ decomposition,process_grid
    ! MPI is initialized by the selected main, using THREAD_FUNNELED.
    mc_gpu_failure=>mc_cuda_abort
    process_grid=0;decomposition='pencil'
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)call mc_cuda_abort()
    do
      read(unit,'(A)',iostat=ios)line
      if(ios/=0)exit
      if(index(adjustl(line),'&multicomponent_parallel')==1)then
        backspace(unit);read(unit,nml=multicomponent_parallel,iostat=ios)
        if(ios/=0)call mc_cuda_abort()
        exit
      end if
    end do
    close(unit)
    if(decomposition/='pencil')then
      print *, 'ERROR: CUDA MPI requires pencil decomposition'
      call mc_cuda_abort()
    end if
    grid=process_grid
    call initialize_mc_pencil(domain,global,l,grid,.false.)
    call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,mc_rank,MPI_INFO_NULL,shared,ierr)
    call MPI_Comm_rank(shared,local_rank,ierr)
    call MPI_Comm_free(shared,ierr)
    device=local_rank
    call get_environment_variable('NSE_MC_CUDA_DEVICE',device_text,status=ios)
    if(ios==0)then
      read(device_text,*,iostat=ios)device
      if(ios/=0.or.device<0)call mc_cuda_abort()
    end if
    call mc_gpu_create(gpu,l,domain%local,r,device,mc_cuda_remote(domain))
    n=[domain%local%nx,domain%local%ny,domain%local%nz]
    allocate(q(n(1),n(2),n(3),l%nvariables))
    count=2*(n(1)+4)*max(n(2)+4,n(3)+4)*l%nvariables
    allocate(send(count),recv(count))
    call initialize_mc_euler_state(q,l,domain%local,global)
    call mc_gpu_upload(gpu,q)
    call mc_gpu_statistics(gpu,local_stats)
    call MPI_Allreduce(local_stats,initial,l%nvariables,MPI_DOUBLE_PRECISION,MPI_SUM,domain%comm,ierr)
    time=0;maximum_steps=0;started=MPI_Wtime()
    if(mc_rank==0)then
      print *, 'Resident CUDA MPI reactive flow; CPU-staged halo only; pencil grid ',domain%dims
      if(r%write_history)then
        open(newunit=history,file=trim(r%history_file),status='replace',iostat=ios)
        if(ios/=0)call mc_cuda_abort()
        call write_reactive_history_header(history,model,l)
      end if
    end if
    if(r%write_history.or.r%write_snapshots)then
      call gather_mc_pencil(q,full,global,domain)
      if(mc_rank==0)then
        if(r%write_history)call write_reactive_history_row(history,0,time,0.0_dp,full,l,global)
        if(r%write_snapshots)call write_reactive_snapshot(0,full,model,l,global,r%snapshot_prefix)
        deallocate(full)
      end if
    end if
    do step=1,global%nsteps
      local_dt=mc_gpu_timestep(gpu,global%dt<=0,global%dt>0)
      call MPI_Allreduce(local_dt,stable,1,MPI_DOUBLE_PRECISION,MPI_MIN,domain%comm,ierr)
      dt=stable
      if(global%dt>0)then
        if(global%dt>stable*(1+100*epsilon(dt)))then
          print *, 'ERROR: CUDA MPI fixed dt violates fluid limit'
          call mc_cuda_abort()
        end if
        dt=global%dt
      end if
      call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,substeps)
      maximum_steps=max(maximum_steps,substeps)
      local_dt=mc_gpu_timestep(gpu,.false.,.true.)
      call MPI_Allreduce(local_dt,stable,1,MPI_DOUBLE_PRECISION,MPI_MIN,domain%comm,ierr)
      if(dt>stable*(1+100*epsilon(dt)))then
        print *, 'ERROR: CUDA MPI post-chemistry fluid limit'
        call mc_cuda_abort()
      end if
      call mc_gpu_begin(gpu)
      do stage=1,3
        call mc_cuda_exchange(gpu,domain,send,recv)
        call mc_gpu_stage(gpu,dt,stage)
      end do
      call mc_gpu_chemistry(gpu,dt/2,r%maximum_chemistry_substeps,substeps)
      maximum_steps=max(maximum_steps,substeps);time=time+dt
      output=mod(step,r%output_every)==0.or.step==global%nsteps
      if(output.and.(r%write_history.or.r%write_snapshots))then
        call mc_gpu_download(gpu,q)
        call gather_mc_pencil(q,full,global,domain)
        if(mc_rank==0)then
          if(r%write_history)call write_reactive_history_row(history,step,time,dt,full,l,global)
          if(r%write_snapshots)call write_reactive_snapshot(step,full,model,l,global,r%snapshot_prefix)
          deallocate(full)
        end if
      end if
      if(mc_rank==0.and.output)write(*,'(I10,2ES24.16)')step,time,dt
    end do
    call mc_gpu_statistics(gpu,local_stats)
    call MPI_Allreduce(local_stats,final,l%nvariables,MPI_DOUBLE_PRECISION,MPI_SUM,domain%comm,ierr)
    call MPI_Allreduce(local_stats(l%nvariables+1),final(l%nvariables+1),4, &
      MPI_DOUBLE_PRECISION,MPI_MIN,domain%comm,ierr)
    call mc_gpu_check_totals(initial,final,l,global)
    if(global%write_final)then
      call mc_gpu_download(gpu,q)
      call gather_mc_pencil(q,full,global,domain)
      if(mc_rank==0)call write_mc_euler_csv(trim(global%output_file),full,model,l,global)
    end if
    call mc_gpu_destroy(gpu)
    call MPI_Allreduce(maximum_steps,global_steps,1,MPI_INTEGER,MPI_MAX,domain%comm,ierr)
    elapsed=MPI_Wtime()-started
    call MPI_Allreduce(elapsed,total_elapsed,1,MPI_DOUBLE_PRECISION,MPI_MAX,domain%comm,ierr)
    if(mc_rank==0)then
      if(r%write_history)close(history)
      print *, 'maximum chemistry substeps per half-step = ',global_steps
      print *, 'elapsed wall seconds (including output) = ',total_elapsed
      print *, 'Reactive multicomponent Navier-Stokes calculation completed successfully.'
    end if
    nullify(mc_gpu_failure)
    call MPI_Comm_free(domain%comm,ierr)
  end subroutine
end module
