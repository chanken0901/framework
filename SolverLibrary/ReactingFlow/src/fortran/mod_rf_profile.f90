module mod_rf_profile
  use mod_rf_flow1d, only: primitive_to_conserved
  use mod_rf_thermo
  implicit none
  private
  public :: read_flow_profile
contains
  subroutine read_flow_profile(path,m,length,q)
    character(*), intent(in) :: path
    type(rf_mechanism), intent(in) :: m
    real(dp), intent(in) :: length
    real(dp), intent(out) :: q(:,:)
    character(256) :: magic,hash,line
    real(dp) :: domain,t,p,u,y(size(m%species))
    integer :: unit,ios,ns,nx,i
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'Cannot open initial profile')
    read(unit,'(a)',iostat=ios) magic
    call require(ios==0.and.trim(magic)=='RF_FLOW_PROFILE_V1','Invalid initial profile version')
    read(unit,'(a)',iostat=ios) hash
    call require(ios==0.and.trim(hash)==m%canonical_hash,'Initial profile mechanism hash mismatch')
    ns=0;nx=0;domain=0
    read(unit,*,iostat=ios) ns,nx,domain
    call require(ios==0.and.ns==size(m%species).and.nx==size(q,2),'Initial profile dimensions mismatch')
    call require(ieee_is_finite(domain).and.abs(domain-length)<=1.e-12_dp*length,'Initial profile length mismatch')
    do i=1,nx
      read(unit,*,iostat=ios) t,p,u,y
      call require(ios==0,'Missing initial profile cell')
      call primitive_to_conserved(m,t,p,u,y,q(:,i))
    end do
    do
      read(unit,'(a)',iostat=ios) line
      if(ios<0) exit
      call require(ios==0.and.len_trim(line)==0,'Extra initial profile records')
    end do
    close(unit)
  end subroutine
end module
