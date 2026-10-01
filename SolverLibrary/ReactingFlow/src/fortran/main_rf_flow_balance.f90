program main_rf_flow_balance
  use mod_rf_flow1d
  use mod_rf_profile
  use mod_rf_thermo
  implicit none
  type(rf_mechanism) :: m
  character(2048) :: mechanism_path,profile_path,output_path
  character(256) :: line
  character(32) :: reconstruction
  character(128) :: name
  integer :: unit,out,ios,ns,nx,i,j
  real(dp) :: length,dx,rho,u,t,p,a
  real(dp), allocatable :: q(:,:),adv(:,:),chem(:,:),boundary(:),flux(:),y(:)
  call require(command_argument_count()==4,'Usage: rf_flow_balance mechanism.rf profile.rf muscl|first_order output.csv')
  call get_command_argument(1,mechanism_path);call get_command_argument(2,profile_path)
  call get_command_argument(3,reconstruction);call get_command_argument(4,output_path)
  call validate_reconstruction(trim(reconstruction))
  call read_mechanism(trim(mechanism_path),m)
  open(newunit=unit,file=trim(profile_path),status='old',action='read',iostat=ios)
  call require(ios==0,'Cannot open balance profile')
  read(unit,'(a)',iostat=ios) line
  call require(ios==0.and.trim(line)=='RF_FLOW_PROFILE_V1','Invalid balance profile version')
  read(unit,'(a)',iostat=ios) line
  call require(ios==0.and.trim(line)==m%canonical_hash,'Balance mechanism hash mismatch')
  ns=0;nx=0;length=0
  read(unit,*,iostat=ios) ns,nx,length;close(unit)
  call require(ios==0.and.ns==size(m%species).and.nx>=2,'Invalid balance profile dimensions')
  call require(ieee_is_finite(length).and.length>0,'Invalid balance domain length')
  dx=length/nx
  allocate(q(ns+2,nx),adv(ns+2,nx),chem(ns+2,nx),boundary(ns+2),flux(ns+2),y(ns))
  call read_flow_profile(trim(profile_path),m,length,q)
  call inviscid_balance(m,q,dx,'outflow','outflow',adv,chem,boundary,trim(reconstruction))
  open(newunit=out,file=trim(output_path),status='new',action='write',iostat=ios)
  call require(ios==0,'Cannot create balance output; file must not already exist')
  write(out,'(a)') '# semidiscrete_inviscid_balance; boundary=outflow; exclude endpoint-adjacent cells in interior comparisons'
  write(out,'(a)') '# canonical_sha256='//m%canonical_hash
  write(out,'(a)') '# reconstruction='//trim(reconstruction)
  write(out,'(a)',advance='no') 'x,density,velocity,temperature,pressure,mass_flux,momentum_flux,energy_flux'
  do j=1,ns+2
    if(j<=ns) then
      name=m%species(j)%name
    else if(j==ns+1) then
      name='momentum'
    else
      name='energy'
    end if
    write(out,'(a)',advance='no') ',advection_'//trim(name)//',chemistry_'//trim(name)//',residual_'//trim(name)
  end do
  write(out,*)
  do i=1,nx
    call conserved_to_primitive(m,q(:,i),rho,u,t,p,a,y)
    call physical_flux(m,q(:,i),flux,a)
    write(out,'(es25.16e3,7(",",es25.16e3))',advance='no') &
      (i-.5_dp)*dx,rho,u,t,p,sum(flux(:ns)),flux(ns+1),flux(ns+2)
    do j=1,ns+2
      write(out,'(3(",",es25.16e3))',advance='no') adv(j,i),chem(j,i),adv(j,i)+chem(j,i)
    end do
    write(out,*)
  end do
  write(out,'(a)') '# SUCCESS'
  close(out)
  write(*,*) '[OK] Fortran inviscid balance diagnostic completed'
end program
