program rf_grid2d_unit
  use mod_rf_grid2d
  use mod_rf_thermo, only: dp,require
  implicit none
  type(rf_grid2d) :: grid
  real(dp) :: nodes(2,0:2,0:1),x(3),lower(3),upper(3)
  character(32) :: mode
  integer :: i,j,f
  x=[0._dp,1._dp,2._dp];lower=0;upper=1
  if(command_argument_count()>0) then
    call get_command_argument(1,mode)
    select case(trim(mode))
    case('reversed')
      x=[0._dp,2._dp,1._dp]
      call build_nozzle2d(x,lower,upper,2,grid)
    case('closed')
      upper(2)=0
      call build_nozzle2d(x,lower,upper,2,grid)
    case('concave')
      nodes(:,0,0)=[0._dp,0._dp];nodes(:,1,0)=[1._dp,0._dp]
      nodes(:,1,1)=[.2_dp,.2_dp];nodes(:,0,1)=[0._dp,1._dp]
      call build_grid2d(nodes(:,0:1,:),grid)
    case default
      error stop 'Unknown negative test'
    end select
    ! A missing rejection returns success and therefore fails a WILL_FAIL test.
    stop
  end if
  call build_nozzle2d(x,lower,upper,2,grid)
  call require(grid%nx==2.and.grid%ny==2,'Structured dimensions')
  call require(size(grid%mesh%owner)==12,'Unique face count')
  call require(maxval(abs(grid%mesh%volume-.5_dp))<1.e-14_dp,'Rectangle areas')
  call require(maxval(abs(grid%cell_center(:,1)-[.5_dp,.25_dp]))<1.e-14_dp,'Rectangle centroid')
  call require(count(grid%boundary==0)==4,'Interior face classification')
  call require(count(grid%boundary==rf_imin)==2.and.count(grid%boundary==rf_imax)==2,'i boundary counts')
  call require(count(grid%boundary==rf_jmin)==2.and.count(grid%boundary==rf_jmax)==2,'j boundary counts')
  do f=1,size(grid%boundary)
    call require((grid%boundary(f)==0).eqv.(grid%mesh%neighbor(f)>0),'Boundary connectivity')
  end do
  ! Linearly tapered channel: exact polygon area and non-vertex-average centroid.
  call build_nozzle2d([0._dp,2._dp],[0._dp,0._dp],[1._dp,2._dp],1,grid)
  call require(abs(grid%mesh%volume(1)-3)<1.e-14_dp,'Trapezoid area')
  call require(maxval(abs(grid%cell_center(:,1)-[10._dp/9,7._dp/9]))<1.e-14_dp,'Trapezoid centroid')
  ! Generic skew nodes also support affine shears, not just vertical sections.
  do j=0,1
    do i=0,2
      nodes(:,i,j)=[real(i,dp)+.2_dp*j,real(j,dp)+.1_dp*i]
    end do
  end do
  call build_grid2d(nodes,grid)
  call require(maxval(abs(grid%mesh%volume-.98_dp))<1.e-14_dp,'Skew quadrilateral area')
  write(*,'(a)') '[OK] structured/nozzle geometry, orientation, centroids and boundary tags'
end program
