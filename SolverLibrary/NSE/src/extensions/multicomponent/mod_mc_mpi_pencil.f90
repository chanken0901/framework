module mod_mc_mpi_pencil
  use mod_precision, only : dp
  use mod_mc_config, only : mc_config
  use mod_mc_state_layout, only : mc_state_layout
  use mod_mc_euler_config, only : mc_euler_config
  use mod_mc_euler_field, only : initialize_mc_euler_state, compute_mc_euler_totals, &
    compute_mc_euler_minima
  use mod_mc_viscous_flux, only : compute_mc_navier_stokes_rhs, &
    compute_mc_navier_stokes_timestep, mc_navier_stokes_workspace
  use mod_mc_reactive_config, only : mc_reactive_config
  use mod_mc_reactive_solver, only : apply_mc_field_chemistry, compute_mc_reactive_timestep, &
    write_reactive_snapshot, write_reactive_history_header, write_reactive_history_row
  use mod_mc_euler_solver, only : write_mc_euler_csv
  use mod_mc_boundary, only : mc_all_boundaries_periodic
  implicit none
  include 'mpif.h'
  private
  integer, public :: mc_rank=0, mc_processes=1
  type, public :: mc_pencil_domain
    integer :: comm=MPI_COMM_NULL, dims(2)=0, coords(2)=0
    integer :: first(2)=1, count(2)=0, lower(2)=0, upper(2)=0
    integer :: neighbor_minus(2),neighbor_plus(2)
    type(mc_euler_config) :: local, padded
    real(dp), allocatable :: halo(:,:,:,:), halo_rhs(:,:,:,:)
    real(dp), allocatable :: send_buffer(:),recv_buffer(:)
    type(mc_navier_stokes_workspace) :: work
  end type mc_pencil_domain
  public :: mc_mpi_start,mc_mpi_finish,initialize_mc_pencil
  public :: compute_mc_pencil_rhs,advance_mc_pencil_strang,gather_mc_pencil
  public :: run_mc_reactive_mpi
contains
  subroutine mc_mpi_start()
    integer :: ierr,provided
    call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,mc_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,mc_processes,ierr)
    if (provided < MPI_THREAD_FUNNELED) call fail('MPI_THREAD_FUNNELED is required')
  end subroutine

  subroutine mc_mpi_finish()
    integer :: ierr
    call MPI_Finalize(ierr)
  end subroutine

  subroutine fail(message)
    character(len=*), intent(in) :: message
    integer :: ierr
    if (mc_rank == 0) write(*,'(A)') 'ERROR: '//trim(message)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    error stop 'MPI multicomponent failure'
  end subroutine

  subroutine partition(n,p,c,first,count)
    integer,intent(in) :: n,p,c
    integer,intent(out) :: first,count
    count=n/p+merge(1,0,c < mod(n,p))
    first=c*(n/p)+min(c,mod(n,p))+1
  end subroutine

  subroutine initialize_mc_pencil(domain,global,layout,process_grid,allocate_host_halo)
    type(mc_pencil_domain),intent(out) :: domain
    type(mc_euler_config),intent(in) :: global
    type(mc_state_layout),intent(in) :: layout
    integer,intent(in) :: process_grid(2)
    logical,intent(in),optional :: allocate_host_halo
    logical :: periods(2)
    integer :: ierr,d,extent(2),rank,n(3),buffer_size
    real(dp) :: spacing(2),origin(2),lo(2),hi(2)

    domain%dims=process_grid
    if (any(domain%dims < 0)) call fail('process_grid must be non-negative')
    if (all(domain%dims > 0)) then
      if (product(domain%dims) /= mc_processes) call fail('process_grid product must equal MPI size')
    else
      do d=1,2
        if (domain%dims(d) > 0) then
          if (mod(mc_processes,domain%dims(d)) /= 0) call fail('invalid process_grid factor')
        end if
      end do
    end if
    call MPI_Dims_create(mc_processes,2,domain%dims,ierr)
    if (ierr /= MPI_SUCCESS) call fail('MPI_Dims_create failed')
    extent=[global%ny,global%nz]
    do d=1,2
      if (domain%dims(d) > 1 .and. extent(d)/domain%dims(d) < 2) &
        call fail('each split pencil direction needs at least two owned cells per rank')
    end do
    periods=[global%boundary_face_types(3) == 'periodic',global%boundary_face_types(5) == 'periodic']
    call MPI_Cart_create(MPI_COMM_WORLD,2,domain%dims,periods,.false.,domain%comm,ierr)
    if (ierr /= MPI_SUCCESS) call fail('MPI_Cart_create failed')
    call MPI_Comm_rank(domain%comm,rank,ierr)
    call MPI_Cart_coords(domain%comm,rank,2,domain%coords,ierr)
    do d=1,2
      call partition(extent(d),domain%dims(d),domain%coords(d),domain%first(d),domain%count(d))
      call MPI_Cart_shift(domain%comm,d-1,1,domain%neighbor_minus(d),domain%neighbor_plus(d),ierr)
      if (domain%dims(d) == 1) then
        domain%neighbor_minus(d)=MPI_PROC_NULL
        domain%neighbor_plus(d)=MPI_PROC_NULL
      end if
      domain%lower(d)=merge(2,0,domain%neighbor_minus(d) /= MPI_PROC_NULL)
      domain%upper(d)=merge(2,0,domain%neighbor_plus(d) /= MPI_PROC_NULL)
    end do
    origin=[global%y_min,global%z_min]
    spacing=[global%y_max-global%y_min,global%z_max-global%z_min]/real(extent,dp)
    lo=origin+real(domain%first-1,dp)*spacing
    hi=lo+real(domain%count,dp)*spacing
    domain%local=global
    domain%local%ny=domain%count(1)
    domain%local%nz=domain%count(2)
    domain%local%y_min=lo(1)
    domain%local%y_max=hi(1)
    domain%local%z_min=lo(2)
    domain%local%z_max=hi(2)
    domain%local%global_boundary_lengths=[global%x_max-global%x_min, &
      global%y_max-global%y_min,global%z_max-global%z_min]
    domain%padded=domain%local
    domain%padded%ny=domain%count(1)+domain%lower(1)+domain%upper(1)
    domain%padded%nz=domain%count(2)+domain%lower(2)+domain%upper(2)
    domain%padded%y_min=lo(1)-real(domain%lower(1),dp)*spacing(1)
    domain%padded%y_max=hi(1)+real(domain%upper(1),dp)*spacing(1)
    domain%padded%z_min=lo(2)-real(domain%lower(2),dp)*spacing(2)
    domain%padded%z_max=hi(2)+real(domain%upper(2),dp)*spacing(2)
    do d=1,2
      ! Artificial outer edges are two layers beyond the owned domain.
      if (domain%lower(d) > 0) domain%padded%boundary_face_types(2*d+1)='reflective'
      if (domain%upper(d) > 0) domain%padded%boundary_face_types(2*d+2)='reflective'
    end do
    n=[global%nx,domain%padded%ny,domain%padded%nz]
    if(present(allocate_host_halo))then
      if(.not.allocate_host_halo)return
    end if
    allocate(domain%halo(n(1),n(2),n(3),layout%nvariables))
    allocate(domain%halo_rhs,mold=domain%halo)
    buffer_size=2*n(1)*max(n(2),n(3))*layout%nvariables
    allocate(domain%send_buffer(buffer_size),domain%recv_buffer(buffer_size))
  end subroutine

  subroutine exchange_halo(q,domain)
    real(dp),intent(in) :: q(:,:,:,:)
    type(mc_pencil_domain),intent(inout) :: domain
    integer :: j,k,j0,k0,ny,nz,nx,nv,count,ierr,status(MPI_STATUS_SIZE)
    nx=size(q,1)
    ny=size(q,2)
    nz=size(q,3)
    nv=size(q,4)
    j0=domain%lower(1)
    k0=domain%lower(2)
    ! Seed corners; the y exchange followed by full-width z exchange replaces
    ! every inter-rank corner with its true diagonal neighbor.
    do k=1,size(domain%halo,3)
      do j=1,size(domain%halo,2)
        domain%halo(:,j,k,:)=q(:,min(ny,max(1,j-j0)),min(nz,max(1,k-k0)),:)
      end do
    end do
    if (domain%dims(1) > 1) then
      count=nx*2*size(domain%halo,3)*nv
      domain%send_buffer(1:count)=reshape(domain%halo(:,j0+1:j0+2,:,:),[count])
      call MPI_Sendrecv(domain%send_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_minus(1),101, &
        domain%recv_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_plus(1),101,domain%comm,status,ierr)
      if (domain%upper(1) > 0) domain%halo(:,j0+ny+1:j0+ny+2,:,:)= &
        reshape(domain%recv_buffer(1:count),[nx,2,size(domain%halo,3),nv])
      domain%send_buffer(1:count)=reshape(domain%halo(:,j0+ny-1:j0+ny,:,:),[count])
      call MPI_Sendrecv(domain%send_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_plus(1),102, &
        domain%recv_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_minus(1),102,domain%comm,status,ierr)
      if (j0 > 0) domain%halo(:,1:2,:,:)=reshape(domain%recv_buffer(1:count),[nx,2,size(domain%halo,3),nv])
    end if
    if (domain%dims(2) > 1) then
      count=nx*size(domain%halo,2)*2*nv
      domain%send_buffer(1:count)=reshape(domain%halo(:,:,k0+1:k0+2,:),[count])
      call MPI_Sendrecv(domain%send_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_minus(2),103, &
        domain%recv_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_plus(2),103,domain%comm,status,ierr)
      if (domain%upper(2) > 0) domain%halo(:,:,k0+nz+1:k0+nz+2,:)= &
        reshape(domain%recv_buffer(1:count),[nx,size(domain%halo,2),2,nv])
      domain%send_buffer(1:count)=reshape(domain%halo(:,:,k0+nz-1:k0+nz,:),[count])
      call MPI_Sendrecv(domain%send_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_plus(2),104, &
        domain%recv_buffer,count,MPI_DOUBLE_PRECISION,domain%neighbor_minus(2),104,domain%comm,status,ierr)
      if (k0 > 0) domain%halo(:,:,1:2,:)=reshape(domain%recv_buffer(1:count),[nx,size(domain%halo,2),2,nv])
    end if
  end subroutine

  subroutine compute_mc_pencil_rhs(q,rhs,layout,domain)
    real(dp),intent(in) :: q(:,:,:,:)
    real(dp),intent(out) :: rhs(:,:,:,:)
    type(mc_state_layout),intent(in) :: layout
    type(mc_pencil_domain),intent(inout) :: domain
    integer :: j,k
    call exchange_halo(q,domain)
    call compute_mc_navier_stokes_rhs(domain%halo,domain%halo_rhs,layout,domain%padded,domain%work)
    j=domain%lower(1)
    k=domain%lower(2)
    rhs=domain%halo_rhs(:,j+1:j+size(q,2),k+1:k+size(q,3),:)
  end subroutine

  subroutine advance_mc_pencil_strang(q,q0,rhs,dt,layout,reactive,domain,substeps)
    real(dp),intent(inout) :: q(:,:,:,:),q0(:,:,:,:),rhs(:,:,:,:)
    real(dp),intent(in) :: dt
    type(mc_state_layout),intent(in) :: layout
    type(mc_reactive_config),intent(in) :: reactive
    type(mc_pencil_domain),intent(inout) :: domain
    integer,intent(out) :: substeps
    integer :: first,second,ierr
    real(dp) :: checked
    type(mc_euler_config) :: fixed
    call apply_mc_field_chemistry(q,dt/2,layout,domain%local,reactive,first)
    fixed=domain%local
    fixed%dt=dt
    checked=compute_mc_navier_stokes_timestep(q,layout,fixed)
    q0=q
    call compute_mc_pencil_rhs(q,rhs,layout,domain)
    q=q0+checked*rhs
    call compute_mc_pencil_rhs(q,rhs,layout,domain)
    q=0.75_dp*q0+0.25_dp*(q+checked*rhs)
    call compute_mc_pencil_rhs(q,rhs,layout,domain)
    q=q0/3.0_dp+(2.0_dp/3.0_dp)*(q+checked*rhs)
    call apply_mc_field_chemistry(q,dt/2,layout,domain%local,reactive,second)
    first=max(first,second)
    call MPI_Allreduce(first,substeps,1,MPI_INTEGER,MPI_MAX,domain%comm,ierr)
  end subroutine

  subroutine gather_mc_pencil(q,global_q,global,domain)
    real(dp),intent(in) :: q(:,:,:,:)
    real(dp),allocatable,intent(out) :: global_q(:,:,:,:)
    type(mc_euler_config),intent(in) :: global
    type(mc_pencil_domain),intent(in) :: domain
    real(dp),allocatable :: received(:)
    integer :: counts(mc_processes),offsets(mc_processes),r,c(2),f(2),n(2),ierr,nv
    nv=size(q,4)
    do r=0,mc_processes-1
      call MPI_Cart_coords(domain%comm,r,2,c,ierr)
      call partition(global%ny,domain%dims(1),c(1),f(1),n(1))
      call partition(global%nz,domain%dims(2),c(2),f(2),n(2))
      counts(r+1)=global%nx*product(n)*nv
    end do
    offsets(1)=0
    do r=2,mc_processes
      offsets(r)=offsets(r-1)+counts(r-1)
    end do
    allocate(received(merge(sum(counts),1,mc_rank == 0)))
    call MPI_Gatherv(q,size(q),MPI_DOUBLE_PRECISION,received,counts,offsets, &
      MPI_DOUBLE_PRECISION,0,domain%comm,ierr)
    if (mc_rank == 0) then
      allocate(global_q(global%nx,global%ny,global%nz,nv))
      do r=0,mc_processes-1
        call MPI_Cart_coords(domain%comm,r,2,c,ierr)
        call partition(global%ny,domain%dims(1),c(1),f(1),n(1))
        call partition(global%nz,domain%dims(2),c(2),f(2),n(2))
        global_q(:,f(1):f(1)+n(1)-1,f(2):f(2)+n(2)-1,:)= &
          reshape(received(offsets(r+1)+1:offsets(r+1)+counts(r+1)),[global%nx,n(1),n(2),nv])
      end do
    end if
  end subroutine

  subroutine run_mc_reactive_mpi(model,layout,global,reactive,path)
    type(mc_config),intent(in) :: model
    type(mc_state_layout),intent(in) :: layout
    type(mc_euler_config),intent(in) :: global
    type(mc_reactive_config),intent(in) :: reactive
    character(len=*),intent(in) :: path
    type(mc_pencil_domain) :: domain
    real(dp),allocatable :: q(:,:,:,:),q0(:,:,:,:),rhs(:,:,:,:),full(:,:,:,:)
    real(dp) :: initial(layout%nvariables),totals(layout%nvariables),local_totals(layout%nvariables)
    real(dp) :: dt,local_dt,time,minimum_species,minimum_density,minimum_pressure,minimum_temperature
    real(dp) :: error,started,elapsed
    real(dp) :: local_minima(4),global_minima(4)
    integer :: process_grid(2),unit,ios,ierr,step,substeps,maximum_substeps,history_unit
    character(len=32) :: decomposition
    character(len=512) :: line
    logical :: write_step
    namelist /multicomponent_parallel/ decomposition,process_grid

    decomposition='pencil'
    process_grid=0
    open(newunit=unit,file=path,status='old',action='read')
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (index(adjustl(line),'&multicomponent_parallel') == 1) then
        backspace(unit)
        read(unit,nml=multicomponent_parallel,iostat=ios)
        if (ios /= 0) call fail('invalid multicomponent_parallel namelist')
        exit
      end if
    end do
    close(unit)
    if (decomposition /= 'pencil') call fail('Stage 8 requires decomposition=pencil')
    call initialize_mc_pencil(domain,global,layout,process_grid)
    allocate(q(global%nx,domain%count(1),domain%count(2),layout%nvariables))
    allocate(q0,mold=q)
    allocate(rhs,mold=q)
    call initialize_mc_euler_state(q,layout,domain%local,global)
    call compute_mc_euler_totals(q,domain%local,local_totals)
    call MPI_Allreduce(local_totals,initial,layout%nvariables,MPI_DOUBLE_PRECISION,MPI_SUM,domain%comm,ierr)
    time=0.0_dp
    maximum_substeps=0
    history_unit=-1
    started=MPI_Wtime()
    if (mc_rank == 0) then
      write(*,'(A,I0,A,2I5)') 'MPI ranks = ',mc_processes,'; y-z pencil process grid = ',domain%dims
      if (reactive%write_history) then
        open(newunit=history_unit,file=trim(reactive%history_file),status='replace',iostat=ios)
        if (ios /= 0) call fail('cannot open reactive history output')
        call write_reactive_history_header(history_unit,model,layout)
      end if
    end if
    if (reactive%write_history .or. reactive%write_snapshots) then
      call gather_mc_pencil(q,full,global,domain)
      if (mc_rank == 0) then
        if (reactive%write_history) call write_reactive_history_row(history_unit,0,time,0.0_dp,full,layout,global)
        if (reactive%write_snapshots) call write_reactive_snapshot(0,full,model,layout,global,reactive%snapshot_prefix)
        deallocate(full)
      end if
    end if
    do step=1,global%nsteps
      local_dt=compute_mc_reactive_timestep(q,layout,domain%local,reactive)
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,domain%comm,ierr)
      call advance_mc_pencil_strang(q,q0,rhs,dt,layout,reactive,domain,substeps)
      maximum_substeps=max(maximum_substeps,substeps)
      time=time+dt
      write_step=mod(step,reactive%output_every) == 0 .or. step == global%nsteps
      if (write_step .and. (reactive%write_history .or. reactive%write_snapshots)) then
        call gather_mc_pencil(q,full,global,domain)
        if (mc_rank == 0) then
          if (reactive%write_history) call write_reactive_history_row(history_unit,step,time,dt,full,layout,global)
          if (reactive%write_snapshots) &
            call write_reactive_snapshot(step,full,model,layout,global,reactive%snapshot_prefix)
          deallocate(full)
        end if
      end if
    end do
    call compute_mc_euler_totals(q,domain%local,local_totals)
    call MPI_Allreduce(local_totals,totals,layout%nvariables,MPI_DOUBLE_PRECISION,MPI_SUM,domain%comm,ierr)
    if (mc_all_boundaries_periodic(global)) then
      error=abs(sum(totals(1:layout%nspecies))-sum(initial(1:layout%nspecies)))/ &
        max(1.0_dp,abs(sum(initial(1:layout%nspecies))))
      error=max(error,maxval(abs(totals(layout%momentum)-initial(layout%momentum))/ &
        max(1.0_dp,abs(initial(layout%momentum)))))
      error=max(error,abs(totals(layout%total_energy)-initial(layout%total_energy))/ &
        max(1.0_dp,abs(initial(layout%total_energy))))
      if (error > 1.0e-10_dp) call fail('MPI reactive conservation check failed')
    end if
    call compute_mc_euler_minima(q,layout,domain%local,minimum_species,minimum_density, &
      minimum_pressure,minimum_temperature)
    local_minima=[minimum_species,minimum_density,minimum_pressure,minimum_temperature]
    call MPI_Allreduce(local_minima,global_minima,4,MPI_DOUBLE_PRECISION,MPI_MIN,domain%comm,ierr)
    if (minimum_species < -1.0e-12_dp .or. minimum_density <= 0.0_dp .or. &
        minimum_pressure <= 0.0_dp .or. minimum_temperature <= 0.0_dp) call fail('MPI positivity check failed')
    if (global%write_final) then
      call gather_mc_pencil(q,full,global,domain)
      if (mc_rank == 0) call write_mc_euler_csv(trim(global%output_file),full,model,layout,global)
    end if
    local_dt=MPI_Wtime()-started
    call MPI_Allreduce(local_dt,elapsed,1,MPI_DOUBLE_PRECISION,MPI_MAX,domain%comm,ierr)
    if (mc_rank == 0) then
      if (reactive%write_history) close(history_unit)
      write(*,'(A,ES14.6)') 'final time = ',time
      write(*,'(A,I0)') 'maximum chemistry substeps per half-step = ',maximum_substeps
      if (mc_all_boundaries_periodic(global)) then
        write(*,'(A,ES14.6)') 'maximum relative conservation error = ',error
      end if
      write(*,'(A,ES14.6)') 'minimum species partial density = ',global_minima(1)
      write(*,'(A,ES14.6)') 'minimum mixture density = ',global_minima(2)
      write(*,'(A,ES14.6)') 'minimum pressure = ',global_minima(3)
      write(*,'(A,ES14.6)') 'minimum temperature = ',global_minima(4)
      write(*,'(A,ES14.6)') 'elapsed wall seconds (including output) = ',elapsed
      write(*,'(A)') 'Reactive multicomponent Navier-Stokes calculation completed successfully.'
    end if
    call MPI_Comm_free(domain%comm,ierr)
  end subroutine
end module mod_mc_mpi_pencil
