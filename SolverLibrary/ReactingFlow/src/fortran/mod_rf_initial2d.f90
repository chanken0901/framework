module mod_rf_initial2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mod_rf_thermo
  use mod_rf_grid2d
  use mod_rf_finite_volume
  implicit none
  private
  public :: initialize_region2d,read_initial2d,write_initial2d
contains
  subroutine initialize_region2d(grid,background,region,mode,axis,position,box,q)
    type(rf_grid2d), intent(in) :: grid
    real(dp), intent(in) :: background(:),region(:),position,box(4)
    character(*), intent(in) :: mode
    integer, intent(in) :: axis
    real(dp), intent(out) :: q(:,:)
    logical :: selected
    integer :: c
    call require(size(background)==size(region).and.size(q,1)==size(background),'Initial state size mismatch')
    call require(size(q,2)==size(grid%cell_center,2),'Initial field size mismatch')
    call require(mode=='split'.or.mode=='box','Unknown region initialization')
    if(mode=='split') then
      call require(axis==1.or.axis==2,'Split axis must be 1 or 2')
      call require(ieee_is_finite(position),'Nonfinite split position')
    else
      call require(all(ieee_is_finite(box)).and.box(2)>box(1).and.box(4)>box(3),'Invalid initial box')
    end if
    do c=1,size(q,2)
      if(mode=='split') then
        selected=grid%cell_center(axis,c)<position
      else
        selected=grid%cell_center(1,c)>=box(1).and.grid%cell_center(1,c)<box(2).and. &
          grid%cell_center(2,c)>=box(3).and.grid%cell_center(2,c)<box(4)
      end if
      q(:,c)=background
      if(selected) q(:,c)=region
    end do
  end subroutine

  subroutine read_initial2d(path,m,grid,q)
    character(*), intent(in) :: path
    type(rf_mechanism), intent(in) :: m
    type(rf_grid2d), intent(in) :: grid
    real(dp), intent(out) :: q(:,:)
    character(256) :: line
    real(dp) :: point(2),state(size(m%species)+4),tol(2),extent(2)
    integer :: unit,ios,ns,nx,ny,i,j,c,index
    call require(size(q,1)==size(m%species)+3.and.size(q,2)==grid%nx*grid%ny,'Initial 2D shape mismatch')
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'Cannot open initial 2D profile')
    read(unit,'(a)',iostat=ios) line
    call require(ios==0,'Missing initial 2D header')
    call require(trim(line)=='RF_INITIAL2D_V1','Invalid initial 2D version')
    read(unit,'(a)',iostat=ios) line
    call require(ios==0,'Missing initial 2D mechanism hash')
    call require(trim(line)==m%canonical_hash,'Initial 2D mechanism hash mismatch')
    ns=0;nx=0;ny=0;read(unit,*,iostat=ios) ns,nx,ny
    call require(ios==0.and.ns==size(m%species).and.nx==grid%nx.and.ny==grid%ny,'Initial 2D dimensions mismatch')
    do i=1,2
      extent(i)=maxval(grid%nodes(i,:,:))-minval(grid%nodes(i,:,:))
      tol(i)=1.e-12_dp*extent(i)+64*epsilon(1._dp)*maxval(abs(grid%nodes(i,:,:)))
    end do
    do j=0,ny
      do i=0,nx
        point=ieee_value(0._dp,ieee_quiet_nan)
        read(unit,*,iostat=ios) point
        call require(ios==0.and.all(ieee_is_finite(point)),'Invalid initial 2D node')
        call require(all(abs(point-grid%nodes(:,i,j))<=tol),'Initial 2D geometry mismatch')
      end do
    end do
    do c=1,size(q,2)
      index=0;state=ieee_value(0._dp,ieee_quiet_nan)
      read(unit,*,iostat=ios) index,state
      call require(ios==0.and.index==c.and.all(ieee_is_finite(state)),'Invalid/missing initial 2D cell')
      call require(all(state(5:)>=0).and.abs(sum(state(5:))-1)<1.e-12_dp,'Invalid initial 2D composition')
      call primitive_nd(m,state(1),state(2),state(3:4),state(5:),q(:,c))
    end do
    do
      read(unit,'(a)',iostat=ios) line
      if(ios<0) exit
      call require(ios==0,'Initial 2D profile read error')
      call require(len_trim(line)==0,'Extra initial 2D records')
    end do
    close(unit)
  end subroutine

  subroutine write_initial2d(path,m,grid,q)
    character(*), intent(in) :: path
    type(rf_mechanism), intent(in) :: m
    type(rf_grid2d), intent(in) :: grid
    real(dp), intent(in) :: q(:,:)
    real(dp) :: rho,v(2),t,p,a,y(size(m%species))
    integer :: unit,ios,i,j,c
    call require(size(q,1)==size(m%species)+3.and.size(q,2)==grid%nx*grid%ny,'Initial 2D export shape mismatch')
    open(newunit=unit,file=path,status='new',action='write',iostat=ios)
    call require(ios==0,'Cannot create initial 2D profile (existing path?)')
    write(unit,'(a)') 'RF_INITIAL2D_V1',trim(m%canonical_hash)
    write(unit,'(3(i0,1x))') size(m%species),grid%nx,grid%ny
    do j=0,grid%ny
      do i=0,grid%nx
        write(unit,'(2(es25.17,1x))') grid%nodes(:,i,j)
      end do
    end do
    do c=1,size(q,2)
      call conserved_nd(m,q(:,c),rho,v,t,p,a,y)
      write(unit,'(i0,*(1x,es25.17))') c,t,p,v,y
    end do
    close(unit,iostat=ios)
    call require(ios==0,'Failed closing initial 2D output')
  end subroutine
end module
