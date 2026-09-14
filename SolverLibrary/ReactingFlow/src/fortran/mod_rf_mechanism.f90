module mod_rf_mechanism
  use iso_fortran_env, only: dp=>real64, error_unit
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(dp), parameter :: gas_r=8.31446261815324_dp, log_zero=-1.e300_dp
  type :: rf_species
    character(128) :: name
    real(dp) :: mass, pref
    integer :: model
    real(dp), allocatable :: bounds(:), coeff(:,:), atoms(:)
  end type
  type :: rf_plog
    real(dp) :: pressure
    real(dp), allocatable :: arr(:,:)
  end type
  type :: rf_reaction
    integer :: kind, reversible
    real(dp) :: high(4), low(4), bounds(4)
    real(dp), allocatable :: reactants(:), products(:), orders(:), efficiencies(:), params(:), cheb(:,:)
    type(rf_plog), allocatable :: groups(:)
  end type
  type :: rf_mechanism
    character(64) :: source_hash, canonical_hash
    integer :: ne
    type(rf_species), allocatable :: species(:)
    type(rf_reaction), allocatable :: reactions(:)
  end type
contains
  subroutine require(ok, message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not.ok) then
      write(error_unit,'(a)') '[ERROR] '//message
      error stop 1
    end if
  end subroutine

  subroutine read_mechanism(path,m)
    character(*), intent(in) :: path
    type(rf_mechanism), intent(out) :: m
    integer :: u,ios,ns,nr,i,j,k,nreg,npar,ng,nt,np,nterm
    character(64) :: magic
    real(dp) :: residual, scale
    open(newunit=u,file=path,status='old',action='read',iostat=ios)
    call require(ios==0,'Cannot open mechanism: '//path)
    read(u,*,iostat=ios) magic
    call require(ios==0.and.magic=='RFMECH1','Unknown mechanism format')
    read(u,*) m%source_hash
    read(u,*) m%canonical_hash
    read(u,*) ns,m%ne,nr
    call require(ns>0.and.m%ne>0.and.nr>=0,'Invalid mechanism dimensions')
    allocate(m%species(ns),m%reactions(nr))
    do i=1,ns
      associate(s=>m%species(i))
        read(u,*) s%name,s%mass,s%model,nreg,s%pref
        call require(all(ieee_is_finite([s%mass,s%pref])),'Nonfinite species constants')
        call require(s%mass>0.and.s%pref>0.and.nreg>0,'Invalid species data')
        call require(s%model==7.or.s%model==9,'Only NASA7/9 supported')
        allocate(s%bounds(nreg+1),s%coeff(9,nreg),s%atoms(m%ne))
        read(u,*) s%bounds
        do j=1,nreg
          read(u,*) s%coeff(:,j)
        end do
        read(u,*) s%atoms
        call require(all(ieee_is_finite(s%coeff)).and.all(ieee_is_finite(s%bounds)), 'Nonfinite NASA data')
        call require(all(s%bounds>0).and.all(s%bounds(2:)>s%bounds(:nreg)),'Invalid NASA intervals')
        call require(all(s%atoms>=0).and.sum(s%atoms)>0,'Invalid element counts')
        do j=1,i-1
          call require(s%name/=m%species(j)%name,'Duplicate species name')
        end do
      end associate
    end do
    do i=1,nr
      associate(r=>m%reactions(i))
        read(u,*) r%kind,r%reversible,npar,ng,nt,np
        call require(r%reversible==0.or.r%reversible==1,'Invalid reversible flag')
        call require(r%kind>=1.and.r%kind<=7,'Unsupported reaction kind')
        call require(npar>=0.and.ng>=0.and.nt>=0.and.np>=0,'Invalid rate dimensions')
        allocate(r%reactants(ns),r%products(ns),r%orders(ns),r%efficiencies(ns),r%params(npar))
        allocate(r%groups(ng),r%cheb(nt,np))
        read(u,*) r%high
        read(u,*) r%low
        read(u,*) r%reactants
        read(u,*) r%products
        read(u,*) r%orders
        read(u,*) r%efficiencies
        if(npar>0) read(u,*) r%params
        do j=1,ng
          read(u,*) r%groups(j)%pressure,nterm
          call require(nterm>0.and.r%groups(j)%pressure>0,'Invalid PLOG group')
          allocate(r%groups(j)%arr(4,nterm))
          do k=1,nterm
            read(u,*) r%groups(j)%arr(:,k)
          end do
        end do
        if(nt>0) then
          read(u,*) r%bounds
          do j=1,nt
            read(u,*) r%cheb(j,:)
          end do
        end if
        call require(all(ieee_is_finite(r%high)).and.all(ieee_is_finite(r%low)), 'Nonfinite Arrhenius data')
        call require(all(ieee_is_finite(r%reactants)).and.all(ieee_is_finite(r%products)), 'Nonfinite stoichiometry')
        call require(all(ieee_is_finite(r%orders)).and.all(ieee_is_finite(r%efficiencies)), 'Nonfinite rate data')
        call require(all(ieee_is_finite(r%params)).and.all(ieee_is_finite(r%cheb)), 'Nonfinite rate parameters')
        do j=1,ng
          call require(ieee_is_finite(r%groups(j)%pressure), 'Nonfinite PLOG pressure')
          call require(all(ieee_is_finite(r%groups(j)%arr)), 'Nonfinite PLOG Arrhenius')
          call require(all(abs(r%groups(j)%arr(4,:))==1),'Invalid PLOG sign')
          if(j>1) call require(r%groups(j)%pressure>r%groups(j-1)%pressure,'PLOG pressures must increase')
        end do
        call require(all(r%reactants>=0).and.all(r%products>=0).and.all(r%orders>=0), 'Negative stoichiometry/order')
        call require(all(r%efficiencies>=0),'Negative collider efficiency')
        call require(r%kind/=4.or.npar==3.or.npar==4,'Invalid Troe dimensions')
        call require(r%kind/=5.or.npar==5,'Invalid SRI dimensions')
        call require(r%kind/=6.or.ng>0,'Empty PLOG')
        call require(r%kind/=7.or.(nt>0.and.np>0),'Empty Chebyshev')
        if(r%kind==4) then
          call require(r%params(1)>=0.and.r%params(1)<=1.and.all(r%params(2:3)>0),'Invalid Troe parameters')
          if(npar==4) call require(r%params(4)>=0,'Invalid Troe T2')
        else if(r%kind==5) then
          call require(r%params(1)>=0.and.r%params(3)>0.and.r%params(4)>0,'Invalid SRI parameters')
        else if(r%kind==7) then
          call require(all(ieee_is_finite(r%bounds)).and.all(r%bounds>0),'Invalid Chebyshev bounds')
          call require(r%bounds(2)>r%bounds(1).and.r%bounds(4)>r%bounds(3),'Invalid Chebyshev intervals')
        end if
        do j=1,m%ne
          residual=0; scale=0
          do k=1,ns
            residual=residual+(r%products(k)-r%reactants(k))*m%species(k)%atoms(j)
            scale=scale+(r%products(k)+r%reactants(k))*m%species(k)%atoms(j)
          end do
          call require(abs(residual)<=1.e-12_dp*max(1._dp,scale),'Unbalanced reaction elements')
        end do
        residual=0; scale=0
        do k=1,ns
          residual=residual+(r%products(k)-r%reactants(k))*m%species(k)%mass
          scale=scale+(r%products(k)+r%reactants(k))*m%species(k)%mass
        end do
        call require(abs(residual)<=1.e-12_dp*max(1._dp,scale),'Unbalanced reaction mass')
      end associate
    end do
    close(u)
  end subroutine
end module
