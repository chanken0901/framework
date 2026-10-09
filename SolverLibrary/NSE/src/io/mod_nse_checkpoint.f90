module mod_nse_checkpoint
  ! Portable, little-endian, rank-independent mapped-grid restart snapshots.
  use iso_fortran_env, only: int32,int64
  use ieee_arithmetic, only: ieee_is_finite
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_model_config, only: nse_config
  implicit none
  private
  public :: append_checkpoint,load_checkpoint
  character(8), parameter :: magic='NSEGRID1',ending='NSREND1'//achar(0)
  integer(int64), parameter :: data_pos=289_int64
  type :: piece_header
    ! version,step,nx,ny,nz,js,je,ks,ke,mapped_norm,viscosity,mapping,periodicity(3)
    integer(int32) :: h(15)
    ! time,gamma,Re,Pr,rho0,mach,bounds(6),stretch(3)
    real(dp) :: r(15)
    integer :: shape4(4),ghost,nparts,rank
  end type
contains
  integer function visc_code(nse) result(code)
    type(nse_config), intent(in) :: nse
    select case(trim(nse%viscous_scheme))
    case('none');code=0
    case('fv2');code=1
    case('central6');code=2
    case default;error stop 'Unsupported checkpoint viscosity'
    end select
  end function

  subroutine append_checkpoint(path,sim,nse,step,time,js,je,ks,ke)
    character(*), intent(in) :: path
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: step,js,je,ks,ke
    real(dp), intent(in) :: time
    type(piece_header) :: p
    integer :: u,ios
    if(sim%grid_mapping=='uniform') return
    p%h=[1,step,sim%nx,sim%ny,sim%nz,js,je,ks,ke,merge(1,0,sim%mapped_keep6),visc_code(nse),1, &
      merge(1,0,nse%boundary_face_type([1,3,5])=='periodic')]
    p%r=[time,nse%gamma,nse%reynolds,nse%prandtl,nse%rho0,nse%mach, &
      sim%x_min,sim%x_max,sim%y_min,sim%y_max,sim%z_min,sim%z_max,sim%grid_stretch]
    open(newunit=u,file=trim(path),status='old',access='stream',form='unformatted', &
      convert='little_endian',action='write',position='append',iostat=ios)
    if(ios/=0) error stop 'Cannot create NSE checkpoint piece'
    write(u,iostat=ios) magic,p%h,p%r,ending
    if(ios/=0) error stop 'Cannot write NSE checkpoint header'
    close(u,iostat=ios)
    if(ios/=0) error stop 'Cannot close NSE checkpoint piece'
  end subroutine

  subroutine read_header(path,p)
    character(*), intent(in) :: path
    type(piece_header), intent(out) :: p
    character(8) :: mark
    integer :: u,ios,i
    integer(int32) :: version,dtype,ndim,shp(4),meta(8),nv
    real(dp) :: time,bounds(6),physics(5)
    character(32) :: name,viscous
    character(32), parameter :: names(5)=[character(32)::'rho','rho_u','rho_v','rho_w','rho_E']
    integer(int64) :: bytes,expected,cells,trailer
    open(newunit=u,file=trim(path),status='old',access='stream',form='unformatted', &
      convert='little_endian',action='read',iostat=ios)
    if(ios/=0) error stop 'Missing NSE checkpoint piece'
    read(u,iostat=ios) mark,version,dtype,ndim,shp,meta,time,bounds,nv
    if(ios/=0) error stop 'Truncated restart SLF header'
    if(mark/='SLF1'//repeat(achar(0),4).or.version/=1.or.dtype/=2.or.ndim/=4.or.nv/=5.or.shp(4)/=5) &
      error stop 'Restart requires five-variable float64 SLF1'
    if(minval(shp)<1.or.meta(6)<0.or.meta(7)<1.or.meta(7)>100000.or.meta(8)/=1) &
      error stop 'Invalid mapped restart SLF header'
    do i=1,5
      read(u,iostat=ios) name
      if(ios/=0.or.name/=names(i)) error stop 'Restart SLF variable order mismatch'
    end do
    cells=int(shp(1),int64)*shp(2)
    if(cells>(huge(cells)-data_pos-300)/40/shp(3)) error stop 'SLF field size overflow'
    trailer=data_pos+40*cells*shp(3)
    read(u,pos=trailer,iostat=ios) mark,physics,viscous
    if(ios/=0.or.mark/='NSEPAR1'//achar(0)) error stop 'Restart SLF lacks physical metadata'
    read(u,iostat=ios) mark,p%h,p%r
    if(ios/=0) error stop 'Truncated NSE checkpoint header'
    if(mark/=magic.or.p%h(1)/=1) error stop 'Invalid NSE checkpoint format/version'
    if(p%h(2)<0.or.minval(p%h(3:5))<1.or.p%h(6)<1.or.p%h(8)<1.or. &
       p%h(7)<p%h(6).or.p%h(9)<p%h(8).or.p%h(7)>p%h(4).or.p%h(9)>p%h(5)) &
      error stop 'Invalid NSE checkpoint extent'
    if(p%h(10)<0.or.p%h(10)>1.or.p%h(11)<0.or.p%h(11)>2.or.p%h(12)/=1.or. &
       any(p%h(13:15)<0).or.any(p%h(13:15)>1)) error stop 'Invalid checkpoint discretization'
    if(.not.all(ieee_is_finite(p%r))) error stop 'Nonfinite checkpoint metadata'
    if(p%r(1)<0.or.p%r(2)<=1.or.p%r(3)<0.or.minval(p%r(4:5))<=0.or.p%r(6)<0.or. &
       any(p%r([8,10,12])<=p%r([7,9,11])).or.any(p%r(13:15)<0).or.any(p%r(13:15)>20)) &
      error stop 'Invalid checkpoint physical parameters'
    if(p%h(11)/=0.and.p%r(3)<=0) error stop 'Invalid checkpoint Reynolds number'
    p%shape4=shp;p%ghost=meta(6);p%nparts=meta(7);p%rank=meta(2)
    if(p%ghost>minval(p%h(3:5))) error stop 'Invalid SLF ghost count'
    if(any(shp(1:3)/=([p%h(3),p%h(7)-p%h(6)+1,p%h(9)-p%h(8)+1]+2*p%ghost))) &
      error stop 'SLF local shape disagrees with grid trailer'
    if(any(meta(3:5)/=p%h(3:5)).or.meta(1)/=p%h(2).or.time/=p%r(1).or. &
       any(bounds/=p%r(7:12)).or.any(physics/=p%r(2:6))) error stop 'SLF header/trailer metadata mismatch'
    select case(p%h(11))
    case(0);if(viscous/='none') error stop 'SLF viscosity mismatch'
    case(1);if(viscous/='fv2') error stop 'SLF viscosity mismatch'
    case(2);if(viscous/='central6') error stop 'SLF viscosity mismatch'
    end select
    expected=trailer-1+80+196
    inquire(unit=u,size=bytes)
    if(bytes/=expected) error stop 'Checkpoint size mismatch or truncated field'
    read(u,pos=expected-7,iostat=ios) mark
    close(u)
    if(ios/=0.or.mark/=ending) error stop 'Incomplete NSE checkpoint piece'
  end subroutine

  subroutine read_index(path,paths,headers)
    character(*), intent(in) :: path
    character(1024), allocatable, intent(out) :: paths(:)
    type(piece_header), allocatable, intent(out) :: headers(:)
    character(1024) :: base
    character(5) :: suffix
    type(piece_header) :: first
    integer :: n,i,j,slash
    integer(int64) :: area
    call read_header(path,first)
    n=first%nparts
    slash=index(trim(path),'_rank',back=.true.)
    if(n>1.and.slash==0) error stop 'Multi-rank restart requires field_STEP_rankNNNNN.slf names'
    base=''
    if(slash>0) base=path(:slash+4)
    allocate(paths(n),headers(n))
    area=0
    do i=1,n
      paths(i)=path
      if(n>1) then
        write(suffix,'(I5.5)') i-1
        paths(i)=trim(base)//suffix//'.slf'
      end if
      call read_header(paths(i),headers(i))
      if(headers(i)%nparts/=n.or.headers(i)%rank/=i-1) error stop 'SLF rank metadata mismatch'
      if(i>1) then
        if(any(headers(i)%h([1,2,3,4,5,10,11,12,13,14,15])/= &
               headers(1)%h([1,2,3,4,5,10,11,12,13,14,15])).or. &
           any(headers(i)%r/=headers(1)%r)) error stop 'Checkpoint pieces belong to different snapshots'
      end if
      do j=1,i-1
        if(max(headers(i)%h(6),headers(j)%h(6))<=min(headers(i)%h(7),headers(j)%h(7)).and. &
           max(headers(i)%h(8),headers(j)%h(8))<=min(headers(i)%h(9),headers(j)%h(9))) &
          error stop 'Overlapping checkpoint pieces'
      end do
      area=area+int(headers(i)%h(7)-headers(i)%h(6)+1,int64)*int(headers(i)%h(9)-headers(i)%h(8)+1,int64)
    end do
    if(area/=int(headers(1)%h(4),int64)*headers(1)%h(5)) error stop 'Missing checkpoint domain coverage'
  end subroutine

  subroutine load_checkpoint(q,sim,nse,js,je,ks,ke)
    type(simulation_config), intent(inout) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js,je,ks,ke
    real(dp), intent(inout) :: q(1-sim%nghost:,js-sim%nghost:,ks-sim%nghost:,:)
    character(1024), allocatable :: paths(:)
    type(piece_header), allocatable :: p(:)
    real(dp) :: geom(9),rho,pressure
    logical :: output_exists
    integer :: a,u,ios,i,j,k,v
    integer(int64) :: offset
    inquire(file=trim(sim%output_dir)//'/.',exist=output_exists)
    if(output_exists) error stop 'Mapped restart requires a new output directory'
    call read_index(sim%restart_file,paths,p)
    if(sim%grid_mapping/='sinh'.or.nse%nv/=5) error stop 'Mapped SLF restart requires single-component NSE'
    if(any(p(1)%h(3:5)/=[sim%nx,sim%ny,sim%nz])) error stop 'Restart grid size mismatch'
    if(p(1)%h(10)/=merge(1,0,sim%mapped_keep6)) error stop 'Restart integration norm mismatch'
    if(any(p(1)%h(13:15)/=merge(1,0,nse%boundary_face_type([1,3,5])=='periodic'))) &
      error stop 'Restart grid periodicity mismatch'
    geom=[sim%x_min,sim%x_max,sim%y_min,sim%y_max,sim%z_min,sim%z_max,sim%grid_stretch]
    if(any(abs(geom-p(1)%r(7:15))>64*epsilon(1._dp)*max(1._dp,abs(geom),abs(p(1)%r(7:15))))) &
      error stop 'Restart grid bounds/stretch mismatch'
    if(sim%nsteps<=p(1)%h(2).or.sim%t_max<=p(1)%r(1)) &
      error stop 'restart requires nsteps and t_max greater than saved step and time'
    do a=1,size(p)
      if(p(a)%h(7)<js.or.p(a)%h(6)>je.or.p(a)%h(9)<ks.or.p(a)%h(8)>ke) cycle
      open(newunit=u,file=trim(paths(a)),status='old',access='stream',form='unformatted', &
        convert='little_endian',action='read',iostat=ios)
      if(ios/=0) error stop 'Cannot open checkpoint field'
      do v=1,5
        do k=max(ks,p(a)%h(8)),min(ke,p(a)%h(9))
          do j=max(js,p(a)%h(6)),min(je,p(a)%h(7))
            offset=((int(v-1,int64)*p(a)%shape4(3)+k-p(a)%h(8)+p(a)%ghost)* &
              p(a)%shape4(2)+j-p(a)%h(6)+p(a)%ghost)*p(a)%shape4(1)+p(a)%ghost
            read(u,pos=data_pos+8*offset,iostat=ios) q(1:sim%nx,j,k,v)
            if(ios/=0) error stop 'Cannot read checkpoint field row'
          end do
        end do
      end do
      close(u)
    end do
    do k=ks,ke
      do j=js,je
        do i=1,sim%nx
          if(.not.all(ieee_is_finite(q(i,j,k,1:5)))) error stop 'Nonfinite checkpoint state'
          rho=q(i,j,k,1)
          if(rho<=nse%small_rho) error stop 'Nonpositive checkpoint density'
          pressure=(nse%gamma-1)*(q(i,j,k,5)-.5_dp*sum(q(i,j,k,2:4)**2)/rho)
          if(.not.ieee_is_finite(pressure).or.pressure<=nse%small_p) error stop 'Nonpositive checkpoint pressure'
        end do
      end do
    end do
    sim%step=p(1)%h(2);sim%t=p(1)%r(1)
    if(sim%rank==0) write(*,'(A,I0,A,ES24.16)') 'Restart loaded: step=',sim%step,', time=',sim%t
  end subroutine
end module
