program test_cuda_fluctuating_mpi
  use module_mpi
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config, update_derived_config
  use mod_model_config, only: nse_config
  use mod_nse_gpu
  use mod_nse_gpu_mpi
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  type(nse_gpu_context) :: full,local
  type(nse_gpu_mpi_halo) :: halo
  real(dp), allocatable :: q(:,:,:,:),part(:,:,:,:)
  real(dp) :: pi,x,y,z,u(3),rho,p,err,totals(5),dt
  integer :: ierr,i,j,k,g,js,je,ks,ke,step,stage,v,status,global_status,device
  character(len=32) :: policy
  character(len=1024) :: visible
  integer :: local_comm, local_rank
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
  sim%nx=8; sim%ny=8; sim%nz=8; sim%nghost=3
  pi=acos(-1.0_dp); sim%x_max=2*pi; sim%y_max=2*pi; sim%z_max=2*pi
  call update_derived_config(sim)
  nse%viscous_scheme='central6'; nse%reynolds=25
  nse%fh_enabled=.true.; nse%fh_boltzmann_number=1.e-4_dp; nse%fh_seed=987
  call mp_setup_division(sim%nx,sim%ny,sim%nz)
  js=j_sta; je=j_end; ks=k_sta; ke=k_end; g=sim%nghost
  allocate(q(1-g:sim%nx+g,1-g:sim%ny+g,1-g:sim%nz+g,5))
  allocate(part(1-g:sim%nx+g,js-g:je+g,ks-g:ke+g,5))
  do k=1-g,sim%nz+g
    do j=1-g,sim%ny+g
      do i=1-g,sim%nx+g
        x=2*pi*modulo(i-1,sim%nx)/sim%nx
        y=2*pi*modulo(j-1,sim%ny)/sim%ny
        z=2*pi*modulo(k-1,sim%nz)/sim%nz
        rho=1+0.03_dp*sin(x+y-z); p=1+0.04_dp*cos(x-y+z)
        u=[0.12_dp*cos(x+y),-0.08_dp*sin(x-z),0.05_dp*cos(y+z)]
        q(i,j,k,1)=rho; q(i,j,k,2:4)=rho*u
        q(i,j,k,5)=p/(nse%gamma-1)+0.5_dp*rho*sum(u*u)
      end do
    end do
  end do
  do v=1,5
    totals(v)=sum(q(1:sim%nx,1:sim%ny,1:sim%nz,v))
  end do
  part=q(:,js-g:je+g,ks-g:ke+g,:)
  device=0; call get_environment_variable('NSE_CUDA_DEVICE_POLICY',policy,status=status)
  if(status/=0) policy='local_rank'
  if(trim(policy)/='fixed') then
    call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,my_rank,MPI_INFO_NULL,local_comm,ierr)
    call MPI_Comm_rank(local_comm,local_rank,ierr)
    call MPI_Comm_free(local_comm,ierr)
    call get_environment_variable('CUDA_VISIBLE_DEVICES',visible,status=status)
    if(status/=0 .or. len_trim(visible)==0 .or. index(visible,',')>0) device=local_rank
  end if
  call nse_gpu_initialize(full,sim,nse,device=device)
  call nse_gpu_initialize(local,sim,nse,local_ny=je-js+1,local_nz=ke-ks+1, &
    distributed_y=ndiv_ny>1,distributed_z=ndiv_nz>1,device=device,global_y_start=js-1,global_z_start=ks-1)
  call nse_gpu_upload(full,q); call nse_gpu_upload(local,part)
  call nse_gpu_mpi_halo_initialize(halo,local,nse)
  dt=1.e-5_dp
  do step=1,4
    call nse_gpu_advance_ssprk3(full,dt)
    call nse_gpu_begin_ssprk3(local,dt)
    do stage=1,3
      call nse_gpu_mpi_exchange(halo,local)
      call nse_gpu_advance_ssprk3_stage(local,dt,stage,status)
      call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(global_status/=0) call MPI_Abort(MPI_COMM_WORLD,21,ierr)
    end do
    call nse_gpu_download(full,q); call nse_gpu_download(local,part)
    ! Compare physical cells only: halo contents depend on the last RK exchange.
    err=maxval(abs(q(1:sim%nx,js:je,ks:ke,:)-part(1:sim%nx,js:je,ks:ke,:)))
    if(err>3.e-12_dp) call MPI_Abort(MPI_COMM_WORLD,22,ierr)
    do v=1,5
      if(abs(sum(q(1:sim%nx,1:sim%ny,1:sim%nz,v))-totals(v))>1.e-10_dp) &
        call MPI_Abort(MPI_COMM_WORLD,23,ierr)
    end do
  end do
  call nse_gpu_mpi_halo_finalize(halo)
  call nse_gpu_finalize(full); call nse_gpu_finalize(local)
  if(my_rank==0) print *, 'MPI resident LLNS matches single-domain CUDA; conservation passed'
  call MPI_Finalize(ierr)
end program
