program test_nonuniform_operators
  use mod_precision, only: dp
  use mod_common_config, only: simulation_config
  use mod_model_config, only: nse_config
  use mod_grid_axis, only: build_sinh_axis,prepare_axis_weno,pack_axis_geometry
  use mod_grid_fvm, only: axis_x,axis_y,axis_z
  use mod_convective_hybrid, only: hybrid_face_weight
  use mod_viscous_fv2, only: add_viscous_fv2_rhs
  implicit none
  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp) :: q(-2:11,-2:11,-2:11,5),rhs(-2:11,-2:11,-2:11,5)
  real(dp), allocatable :: centers(:),widths(:),coefficients(:,:,:)
  real(dp) :: x,u,p,alpha,expected
  integer :: i,j,k,mode
  sim%nx=8;sim%ny=8;sim%nz=8;sim%nghost=3;sim%grid_mapping='sinh'
  call build_sinh_axis(8,3,0._dp,1._dp,2._dp,.false.,axis_x)
  call build_sinh_axis(8,3,0._dp,1._dp,1._dp,.false.,axis_y)
  call build_sinh_axis(8,3,0._dp,1._dp,3._dp,.false.,axis_z)
  sim%dx=axis_x%minimum_width;sim%dy=axis_y%minimum_width;sim%dz=axis_z%minimum_width
  call prepare_axis_weno(axis_x)
  call pack_axis_geometry(axis_x,2,3,centers,widths,coefficients)
  if(size(centers)/=9.or.size(coefficients)/=240) error stop 'Packed local geometry shape'
  if(maxval(abs(centers-axis_x%center(0:8)))>1.e-14_dp) error stop 'Packed global offset'
  nse%reynolds=100; nse%prandtl=.72_dp
  nse%hybrid_sensor='ducros_pressure'
  nse%hybrid_sensor_onset=0; nse%hybrid_sensor_full=.1_dp
  do mode=1,3
    do k=-2,11
      do j=-2,11
        do i=-2,11
          x=axis_x%center(i);u=-x;p=2._dp
          if(mode==1) p=p+.1_dp*x
          if(mode==3) u=.25_dp
          q(i,j,k,:)=[1._dp,u,0._dp,0._dp,p/(nse%gamma-1)+.5_dp*u*u]
        end do
      end do
    end do
    if(mode==1) then
      do i=0,8
        alpha=hybrid_face_weight(q,i,4,4,1,sim,nse,1,1)
        if(alpha>1.e-20_dp) error stop 'Linear pressure incorrectly detected as shock'
      end do
      do k=-2,11
        do j=-2,11
          do i=-2,11
            p=1._dp
            if(i>4) p=4._dp
            q(i,j,k,5)=p/(nse%gamma-1)+.5_dp*q(i,j,k,2)**2
          end do
        end do
      end do
      alpha=hybrid_face_weight(q,4,4,4,1,sim,nse,1,1)
      if(alpha<.9_dp) error stop 'Stretched pressure jump not detected'
    else
      rhs=0
      !$OMP PARALLEL DEFAULT(shared)
      call add_viscous_fv2_rhs(q,rhs,sim,nse,1,8,1,8)
      !$OMP END PARALLEL
      if(maxval(abs(rhs(1:8,1:8,1:8,1:4)))>1.e-11_dp) error stop 'FV2 momentum/continuity'
      expected=0
      if(mode==2) expected=4._dp/(3*nse%reynolds)
      if(maxval(abs(rhs(1:8,1:8,1:8,5)-expected))>1.e-11_dp) error stop 'FV2 viscous work'
    end if
  end do
  print *, '[OK] stretched sensor, FV2 constant/linear flow, local GPU geometry packing'
end program
