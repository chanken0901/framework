module mod_rf_grid2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mod_rf_thermo, only: dp,require,ieee_is_finite
  use mod_rf_finite_volume, only: rf_face_mesh,validate_face_mesh
  implicit none
  private
  public :: rf_grid2d,build_grid2d,build_nozzle2d,read_nozzle2d,write_grid2d_vtk
  public :: rf_imin,rf_imax,rf_jmin,rf_jmax
  integer, parameter :: rf_imin=1,rf_imax=2,rf_jmin=3,rf_jmax=4
  type :: rf_grid2d
    integer :: nx=0,ny=0
    type(rf_face_mesh) :: mesh
    real(dp), allocatable :: nodes(:,:,:),cell_center(:,:),face_center(:,:)
    integer, allocatable :: boundary(:)
  end type
contains
  pure real(dp) function cross(a,b)
    real(dp), intent(in) :: a(2),b(2)
    cross=a(1)*b(2)-a(2)*b(1)
  end function

  subroutine build_grid2d(nodes,grid)
    real(dp), intent(in) :: nodes(:,0:,0:)
    type(rf_grid2d), intent(out) :: grid
    real(dp) :: v(2,4),edge(2,4),s(2),origin(2),a1,a2,scale
    integer :: nx,ny,nc,nf,i,j,c,k,next,f,l,r,tag
    nx=size(nodes,2)-1;ny=size(nodes,3)-1
    call require(size(nodes,1)==2.and.nx>=1.and.ny>=1,'Invalid 2D node shape')
    call require(all(ieee_is_finite(nodes)),'Nonfinite grid nodes')
    ! Check integer products before computing cell/face counts.
    call require(real(nx,dp)*real(ny,dp)<real(huge(nc),dp)/4,'Grid too large for integer indexing')
    nc=nx*ny;nf=(nx+1)*ny+nx*(ny+1)
    grid%nx=nx;grid%ny=ny
    allocate(grid%nodes(2,0:nx,0:ny));grid%nodes=nodes
    allocate(grid%cell_center(2,nc),grid%face_center(2,nf),grid%boundary(nf))
    allocate(grid%mesh%owner(nf),grid%mesh%neighbor(nf),grid%mesh%area_vector(2,nf),grid%mesh%volume(nc))
    do j=1,ny
      do i=1,nx
        c=i+nx*(j-1);origin=nodes(:,i-1,j-1)
        ! Translation reduces cancellation for small cells at large coordinates.
        v(:,1)=0;v(:,2)=nodes(:,i,j-1)-origin
        v(:,3)=nodes(:,i,j)-origin;v(:,4)=nodes(:,i-1,j)-origin
        do k=1,4
          next=mod(k,4)+1;edge(:,k)=v(:,next)-v(:,k)
        end do
        scale=maxval(abs(edge))
        call require(ieee_is_finite(scale).and.scale>0,'Degenerate grid cell')
        do k=1,4
          next=mod(k,4)+1
          call require(cross(edge(:,k)/scale,edge(:,next)/scale)>64*epsilon(scale), &
            'Cell must be counterclockwise, strictly convex and nondegenerate')
        end do
        a1=cross(v(:,2),v(:,3))/2;a2=cross(v(:,3),v(:,4))/2
        grid%mesh%volume(c)=a1+a2
        call require(ieee_is_finite(a1+a2).and.a1+a2>0,'Invalid grid cell area')
        grid%cell_center(:,c)=origin+(a1/(a1+a2)*(v(:,2)+v(:,3))+ &
          a2/(a1+a2)*(v(:,3)+v(:,4)))/3
      end do
    end do
    f=0
    ! i-faces: tangent from j-1 to j; right normal points toward increasing i.
    do j=1,ny
      do i=0,nx
        s=nodes(:,i,j)-nodes(:,i,j-1);s=[s(2),-s(1)]
        tag=0
        if(i==0) then
          l=1+nx*(j-1);r=0;s=-s;tag=rf_imin
        else if(i==nx) then
          l=nx*j;r=0;tag=rf_imax
        else
          l=i+nx*(j-1);r=l+1
        end if
        call append_face(nodes(:,i,j-1),nodes(:,i,j))
      end do
    end do
    ! j-faces: tangent from i-1 to i; left normal points toward increasing j.
    do j=0,ny
      do i=1,nx
        s=nodes(:,i,j)-nodes(:,i-1,j);s=[-s(2),s(1)]
        tag=0
        if(j==0) then
          l=i;r=0;s=-s;tag=rf_jmin
        else if(j==ny) then
          l=i+nx*(ny-1);r=0;tag=rf_jmax
        else
          l=i+nx*(j-1);r=l+nx
        end if
        call append_face(nodes(:,i-1,j),nodes(:,i,j))
      end do
    end do
    call require(all(ieee_is_finite(grid%cell_center)).and.all(ieee_is_finite(grid%face_center)), &
      'Nonfinite grid centers')
    call validate_face_mesh(grid%mesh)
  contains
    subroutine append_face(first,last)
      real(dp), intent(in) :: first(2),last(2)
      f=f+1
      grid%mesh%owner(f)=l;grid%mesh%neighbor(f)=r;grid%mesh%area_vector(:,f)=s
      grid%boundary(f)=tag;grid%face_center(:,f)=first+(last-first)/2
    end subroutine
  end subroutine

  subroutine build_nozzle2d(x,lower,upper,ny,grid)
    real(dp), intent(in) :: x(:),lower(:),upper(:)
    integer, intent(in) :: ny
    type(rf_grid2d), intent(out) :: grid
    real(dp), allocatable :: nodes(:,:,:)
    real(dp) :: eta
    integer :: nx,i,j
    nx=size(x)-1
    call require(nx>=1.and.ny>=1,'Nozzle needs >=2 stations and >=1 transverse cell')
    call require(real(nx,dp)*real(ny,dp)<real(huge(nx),dp)/4,'Nozzle grid too large')
    call require(size(lower)==nx+1.and.size(upper)==nx+1,'Nozzle wall shape mismatch')
    call require(all(ieee_is_finite(x)).and.all(ieee_is_finite(lower)).and.all(ieee_is_finite(upper)), &
      'Nonfinite nozzle profile')
    call require(all(x(2:)>x(:nx)).and.all(upper>lower),'Nozzle x must increase and height must be positive')
    allocate(nodes(2,0:nx,0:ny))
    do j=0,ny
      eta=real(j,dp)/ny
      do i=0,nx
        nodes(1,i,j)=x(i+1)
        nodes(2,i,j)=(1-eta)*lower(i+1)+eta*upper(i+1)
      end do
    end do
    call build_grid2d(nodes,grid)
  end subroutine

  subroutine read_nozzle2d(path,ny,grid)
    character(*), intent(in) :: path
    integer, intent(in) :: ny
    type(rf_grid2d), intent(out) :: grid
    real(dp), allocatable :: x(:),lower(:),upper(:)
    integer :: unit,ios,n,i
    character(256) :: line
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'Cannot open nozzle profile')
    read(unit,'(a)',iostat=ios) line
    call require(ios==0,'Missing nozzle profile header')
    call require(trim(line)=='RF_NOZZLE_PROFILE_V1','Invalid nozzle profile version')
    n=0;read(unit,*,iostat=ios) n
    call require(ios==0.and.n>=2,'Invalid nozzle station count')
    call require(ny>=1,'Invalid transverse cell count')
    call require(real(n-1,dp)*real(ny,dp)<real(huge(n),dp)/4,'Nozzle input too large')
    allocate(x(n),lower(n),upper(n))
    x=ieee_value(0._dp,ieee_quiet_nan);lower=x;upper=x
    do i=1,n
      read(unit,*,iostat=ios) x(i),lower(i),upper(i)
      call require(ios==0,'Incomplete nozzle profile row')
    end do
    ! Reject additional nonempty records instead of silently ignoring another profile.
    do
      read(unit,'(a)',iostat=ios) line
      if(ios<0) exit
      call require(ios==0,'Nozzle profile read error')
      call require(len_trim(line)==0,'Unexpected trailing nozzle profile data')
    end do
    close(unit)
    call build_nozzle2d(x,lower,upper,ny,grid)
  end subroutine

  subroutine write_grid2d_vtk(path,grid)
    character(*), intent(in) :: path
    type(rf_grid2d), intent(in) :: grid
    integer :: out,ios,i,j,c
    ! Refuse overwriting an existing artifact; caller chooses a new path.
    call validate_face_mesh(grid%mesh)
    call require(allocated(grid%nodes),'Missing grid nodes')
    call require(grid%nx>=1.and.grid%ny>=1,'Invalid stored grid dimensions')
    call require(size(grid%nodes,1)==2.and.size(grid%nodes,2)==grid%nx+1.and. &
      size(grid%nodes,3)==grid%ny+1,'Stored grid node shape mismatch')
    call require(real(grid%nx,dp)*real(grid%ny,dp)<real(huge(i),dp)/4,'Stored grid too large')
    call require(size(grid%mesh%volume)==grid%nx*grid%ny,'Stored grid cell shape mismatch')
    open(newunit=out,file=path,status='new',action='write',iostat=ios)
    call require(ios==0,'Cannot create new grid VTK file (existing path?)')
    write(out,'(a)') '# vtk DataFile Version 3.0','ReactingFlow stationary planar grid','ASCII','DATASET STRUCTURED_GRID'
    write(out,'(a,3(1x,i0))') 'DIMENSIONS',grid%nx+1,grid%ny+1,1
    write(out,'(a,1x,i0,1x,a)') 'POINTS',(grid%nx+1)*(grid%ny+1),'double'
    do j=0,grid%ny
      do i=0,grid%nx
        write(out,'(3(es24.16,1x))') grid%nodes(:,i,j),0._dp
      end do
    end do
    write(out,'(a,1x,i0)') 'CELL_DATA',grid%nx*grid%ny
    write(out,'(a)') 'SCALARS cell_area double 1','LOOKUP_TABLE default'
    do c=1,size(grid%mesh%volume)
      write(out,'(es24.16)') grid%mesh%volume(c)
    end do
    close(out,iostat=ios)
    call require(ios==0,'Failed closing grid VTK output')
  end subroutine
end module
