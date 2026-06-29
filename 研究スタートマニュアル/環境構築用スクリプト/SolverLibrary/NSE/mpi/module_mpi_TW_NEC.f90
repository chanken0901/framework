!-----------------------------------------------------------------------------
! MPI routines
!-----------------------------------------------------------------------------
module module_mpi
  !use mpi
  implicit none
  include 'mpif.h'

!  nprocs   = number of processes
!  my_rank  = process rank number
!  root     = root process = rank 0
!  nblock   = minimum cell number 
!  ndiv_nx  = # division in i direction
!  ndiv_ny  = # division in j direction
!  i_sta    = start index of do loop of i
!  i_end    = end index of do loop of i
!  i_down   = neighbor rank for i-down direction
!  i_up     = neighbor rank for i-up direction
!  i_myrank = 
!  j_sta    = start index of do loop of j
!  j_end    = end index of do loop of j
!  j_down   = neighbor rank for j-down direction
!  j_up     = neighbor rank for j-up direction
!  j_myrank = 
!  k_sta    = start index of do loop of k
!  k_end    = end index of do loop of k
!  k_down   = neighbor rank for k-down direction
!  k_up     = neighbor rank for k-up direction
!  k_myrank = 
!  i_len    = loop length of i ( i_end-i_sta+1 )
!  j_len    = loop length of j ( j_end-j_sta+1 )
!  k_len    = loop length of k ( k_end-k_sta+1 )
!  iimax    = maximum of i_len
!  jjmax    = maximum of j_len
!  kkmax    = maximum of k_len
!  itable   = communication table
!  iista    = 
!  iiend    = 
!  jjsta    = 
!  jjend    = 
!  kksta    = 
!  kkend    = 
!  i_allocate   = allocation flag
!  i_deallocate = deallocation flag

  integer :: nprocs
  integer :: my_rank
  integer,parameter :: root = 0
  integer,parameter :: nblock = 2
  integer :: ndiv_nx , ndiv_ny , ndiv_nz
  integer :: i_sta , i_end , i_up , i_down , i_myrank
  integer :: j_sta , j_end , j_up , j_down , j_myrank
  integer :: k_sta , k_end , k_up , k_down , k_myrank
  integer :: i_len , j_len , k_len , iimax , jjmax , kkmax
  integer,allocatable :: itable(:,:)
  integer,allocatable :: iista(:),iiend(:)
  integer,allocatable :: jjsta(:),jjend(:)
  integer,allocatable :: kksta(:),kkend(:)
  integer,parameter :: i_allocate=0 , i_deallocate=1
!NECs
  integer :: inum
  integer :: len_ij1,len_ik1
  integer :: len_ij3,len_ik3
  real,allocatable,dimension(:) :: s_jp1,s_jm1,s_kp1,s_km1
  real,allocatable,dimension(:) :: r_jp1,r_jm1,r_kp1,r_km1
  real,allocatable,dimension(:) :: s_jp3,s_jm3,s_kp3,s_km3
  real,allocatable,dimension(:) :: r_jp3,r_jm3,r_kp3,r_km3
  integer :: ireq1(8),ireq3(8)
!NECe
!KTs
  integer :: len_ij4,len_ik4
  real,allocatable,dimension(:) :: s_jp4,s_jm4,s_kp4,s_km4
  real,allocatable,dimension(:) :: r_jp4,r_jm4,r_kp4,r_km4
  integer :: ireq4(8)
!KTe
  integer :: len_ij1_r8,len_ik1_r8
  real*8,allocatable,dimension(:) :: s_jp1_r8,s_jm1_r8,s_kp1_r8,s_km1_r8
  real*8,allocatable,dimension(:) :: r_jp1_r8,r_jm1_r8,r_kp1_r8,r_km1_r8
  integer :: ireq1_r8(8)

  integer :: len_ij3_r8,len_ik3_r8
  real*8,allocatable,dimension(:) :: s_jp3_r8,s_jm3_r8,s_kp3_r8,s_km3_r8
  real*8,allocatable,dimension(:) :: r_jp3_r8,r_jm3_r8,r_kp3_r8,r_km3_r8
  integer :: ireq3_r8(8)

  contains

!-----------------------------------------------------------------------------
  subroutine mp_setup_division( nx , ny , nz )
    implicit none
    integer,intent(in)    :: nx , ny , nz
    integer :: i , j , k
    integer :: irank , ierr
    integer :: nx_h , ny_h , nz_h , ii_sta , ii_end , jj_sta , jj_end , kk_sta , kk_end

  !--- parallel division
  
    ierr = 0
    call mp_process_check( ny, nz, ierr )
    if ( ierr .ne. 0 ) then
       if ( ierr .eq. 1 ) &
          call error ( 0 , '#CPU should be the power of 2.' , 8 )
       if ( ierr .eq. 2 ) &
          call error ( 0 , 'Parallel job only.' , 8 )
       if ( ierr .eq. 3 ) &
          call error ( 0 , '#CPU should be larger than 4.' , 8 )
       call mp_stop( 10 )
    endif
  
  !--- communication table
  
    allocate( itable(-1:ndiv_ny,-1:ndiv_nz) )
    do k = -1 , ndiv_nz
       do j = -1 , ndiv_ny
!    do j = -1 , ndiv_ny
!       do i = -1 , ndiv_nx
         itable(j,k) = MPI_PROC_NULL
       end do
    end do
    irank = 0
    do k = 0 , ndiv_nz-1
       do j = 0 , ndiv_ny-1
         itable(j,k) = irank
         if( my_rank == irank )then
           j_myrank = j
           k_myrank = k
         end if
         irank = irank + 1
       end do
    end do
  
  !--- calculate start/end for loop i,j
    allocate( jjsta(0:ndiv_ny) )
    allocate( jjend(0:ndiv_ny) )
    allocate( kksta(0:ndiv_nz) )
    allocate( kkend(0:ndiv_nz) )
    ny_h = ( ny - 1 ) / nblock + 1
    nz_h = ( nz - 1 ) / nblock + 1
  
    iimax=nx
    jjmax=-999
    kkmax=-999

    do i = 0 , ndiv_nz-1
       call para_range ( 1 , nz_h , ndiv_nz , i , kk_sta , kk_end )
       kksta(i) = 1 + ( kk_sta - 1 )*nblock
       kkend(i) = min( 1 + kk_end * nblock - 1 , nz )
       if( kkend(i)-kksta(i) + 1 < nblock ) then
         write(*,*) kksta(i),kkend(i)
         call error ( 0 , '*wincol* error at set kksta(i),kkend(i).' , 8 )
         call mp_stop( 0 )
       end if
       kkmax = max( kkmax , kkend(i)-kksta(i)+1 )
    end do
  
    do i = 0 , ndiv_ny-1
       call para_range ( 1 , ny_h , ndiv_ny , i , jj_sta , jj_end )
       jjsta(i) = 1 + ( jj_sta - 1 )*nblock
       jjend(i) = min( 1 + jj_end * nblock - 1 , ny )
       if( jjend(i)-jjsta(i) + 1 < nblock ) then
         write(*,*) jjsta(i),jjend(i)
         call error ( 0 , '*wincol* error at set jjsta(i),jjend(i).' , 8 )
         call mp_stop( 0 )
       end if
       jjmax = max( jjmax , jjend(i)-jjsta(i)+1 )
    end do
    i_sta = 1 ! iista( i_myrank )
    i_end = nx ! iiend( i_myrank )
    j_sta = jjsta( j_myrank )
    j_end = jjend( j_myrank )
    k_sta = kksta(k_myrank)
    k_end = kkend(k_myrank)
  
    j_up   = itable( j_myrank+1   , k_myrank     )
    j_down = itable( j_myrank-1   , k_myrank     )
    k_up   = itable( j_myrank     , k_myrank+1   )
    k_down = itable( j_myrank     , k_myrank-1   )

    i_len = i_end - i_sta + 1
    j_len = j_end - j_sta + 1
    k_len = k_end - k_sta + 1
  
!NECs
    len_ik1 = (Nx+6)*((k_end+1)-(k_sta-1)+1)
    len_ij1 = (Nx+6)*((j_end+1)-(j_sta-1)+1)
    allocate( s_jp1( len_ik1 ) ) 
    allocate( s_jm1( len_ik1 ) ) 
    allocate( s_kp1( len_ij1 ) ) 
    allocate( s_km1( len_ij1 ) ) 
    allocate( r_jp1( len_ik1 ) ) 
    allocate( r_jm1( len_ik1 ) ) 
    allocate( r_kp1( len_ij1 ) ) 
    allocate( r_km1( len_ij1 ) ) 

    len_ik3 = (Nx+6)*((k_end+3)-(k_sta-3)+1)*3
    len_ij3 = (Nx+6)*((j_end+3)-(j_sta-3)+1)*3
    allocate( s_jp3( len_ik3 ) ) 
    allocate( s_jm3( len_ik3 ) ) 
    allocate( s_kp3( len_ij3 ) ) 
    allocate( s_km3( len_ij3 ) ) 
    allocate( r_jp3( len_ik3 ) ) 
    allocate( r_jm3( len_ik3 ) ) 
    allocate( r_kp3( len_ij3 ) ) 
    allocate( r_km3( len_ij3 ) ) 

    call mpi_send_init(s_jm1,len_ik1,MPI_REAL,j_down,1,MPI_COMM_WORLD,ireq1(1),ierr)
    call mpi_send_init(s_jp1,len_ik1,MPI_REAL,j_up  ,2,MPI_COMM_WORLD,ireq1(2),ierr)
    call mpi_send_init(s_km1,len_ij1,MPI_REAL,k_down,3,MPI_COMM_WORLD,ireq1(3),ierr)
    call mpi_send_init(s_kp1,len_ij1,MPI_REAL,k_up  ,4,MPI_COMM_WORLD,ireq1(4),ierr)

    call mpi_recv_init(r_jp1,len_ik1,MPI_REAL,j_up  ,1,MPI_COMM_WORLD,ireq1(5),ierr)
    call mpi_recv_init(r_jm1,len_ik1,MPI_REAL,j_down,2,MPI_COMM_WORLD,ireq1(6),ierr)
    call mpi_recv_init(r_kp1,len_ij1,MPI_REAL,k_up  ,3,MPI_COMM_WORLD,ireq1(7),ierr)
    call mpi_recv_init(r_km1,len_ij1,MPI_REAL,k_down,4,MPI_COMM_WORLD,ireq1(8),ierr)

    call mpi_send_init(s_jm3,len_ik3,MPI_REAL,j_down,1,MPI_COMM_WORLD,ireq3(1),ierr)
    call mpi_send_init(s_jp3,len_ik3,MPI_REAL,j_up  ,2,MPI_COMM_WORLD,ireq3(2),ierr)
    call mpi_send_init(s_km3,len_ij3,MPI_REAL,k_down,3,MPI_COMM_WORLD,ireq3(3),ierr)
    call mpi_send_init(s_kp3,len_ij3,MPI_REAL,k_up  ,4,MPI_COMM_WORLD,ireq3(4),ierr)

    call mpi_recv_init(r_jp3,len_ik3,MPI_REAL,j_up  ,1,MPI_COMM_WORLD,ireq3(5),ierr)
    call mpi_recv_init(r_jm3,len_ik3,MPI_REAL,j_down,2,MPI_COMM_WORLD,ireq3(6),ierr)
    call mpi_recv_init(r_kp3,len_ij3,MPI_REAL,k_up  ,3,MPI_COMM_WORLD,ireq3(7),ierr)
    call mpi_recv_init(r_km3,len_ij3,MPI_REAL,k_down,4,MPI_COMM_WORLD,ireq3(8),ierr)
!NECe
!KTs
    len_ik4 = (Nx+10)*((k_end+4)-(k_sta-4)+1)*4
    len_ij4 = (Nx+10)*((j_end+4)-(j_sta-4)+1)*4
    allocate( s_jp4( len_ik4 ) ) 
    allocate( s_jm4( len_ik4 ) ) 
    allocate( s_kp4( len_ij4 ) ) 
    allocate( s_km4( len_ij4 ) ) 
    allocate( r_jp4( len_ik4 ) ) 
    allocate( r_jm4( len_ik4 ) ) 
    allocate( r_kp4( len_ij4 ) ) 
    allocate( r_km4( len_ij4 ) ) 

    call mpi_send_init(s_jm4,len_ik4,MPI_REAL,j_down,1,MPI_COMM_WORLD,ireq4(1),ierr)
    call mpi_send_init(s_jp4,len_ik4,MPI_REAL,j_up  ,2,MPI_COMM_WORLD,ireq4(2),ierr)
    call mpi_send_init(s_km4,len_ij4,MPI_REAL,k_down,3,MPI_COMM_WORLD,ireq4(3),ierr)
    call mpi_send_init(s_kp4,len_ij4,MPI_REAL,k_up  ,4,MPI_COMM_WORLD,ireq4(4),ierr)

    call mpi_recv_init(r_jp4,len_ik4,MPI_REAL,j_up  ,1,MPI_COMM_WORLD,ireq4(5),ierr)
    call mpi_recv_init(r_jm4,len_ik4,MPI_REAL,j_down,2,MPI_COMM_WORLD,ireq4(6),ierr)
    call mpi_recv_init(r_kp4,len_ij4,MPI_REAL,k_up  ,3,MPI_COMM_WORLD,ireq4(7),ierr)
    call mpi_recv_init(r_km4,len_ij4,MPI_REAL,k_down,4,MPI_COMM_WORLD,ireq4(8),ierr)
!KTe






    len_ik1_r8 = (nx)*((k_end+1)-(k_sta-1)+1)*1
    len_ij1_r8 = (nx)*((j_end+1)-(j_sta-1)+1)*1
    allocate( s_jp1_r8( len_ik1_r8 ) ) 
    allocate( s_jm1_r8( len_ik1_r8 ) ) 
    allocate( s_kp1_r8( len_ij1_r8 ) ) 
    allocate( s_km1_r8( len_ij1_r8 ) ) 
    allocate( r_jp1_r8( len_ik1_r8 ) ) 
    allocate( r_jm1_r8( len_ik1_r8 ) ) 
    allocate( r_kp1_r8( len_ij1_r8 ) ) 
    allocate( r_km1_r8( len_ij1_r8 ) ) 

    call mpi_send_init(s_jm1_r8,len_ik1_r8,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,ireq1_r8(1),ierr)
    call mpi_send_init(s_jp1_r8,len_ik1_r8,MPI_DOUBLE_PRECISION,j_up  ,2,MPI_COMM_WORLD,ireq1_r8(2),ierr)
    call mpi_send_init(s_km1_r8,len_ij1_r8,MPI_DOUBLE_PRECISION,k_down,3,MPI_COMM_WORLD,ireq1_r8(3),ierr)
    call mpi_send_init(s_kp1_r8,len_ij1_r8,MPI_DOUBLE_PRECISION,k_up  ,4,MPI_COMM_WORLD,ireq1_r8(4),ierr)

    call mpi_recv_init(r_jp1_r8,len_ik1_r8,MPI_DOUBLE_PRECISION,j_up  ,1,MPI_COMM_WORLD,ireq1_r8(5),ierr)
    call mpi_recv_init(r_jm1_r8,len_ik1_r8,MPI_DOUBLE_PRECISION,j_down,2,MPI_COMM_WORLD,ireq1_r8(6),ierr)
    call mpi_recv_init(r_kp1_r8,len_ij1_r8,MPI_DOUBLE_PRECISION,k_up  ,3,MPI_COMM_WORLD,ireq1_r8(7),ierr)
    call mpi_recv_init(r_km1_r8,len_ij1_r8,MPI_DOUBLE_PRECISION,k_down,4,MPI_COMM_WORLD,ireq1_r8(8),ierr)




    len_ik3_r8 = (nx+6)*((k_end+3)-(k_sta-3)+1)*3*5
    len_ij3_r8 = (nx+6)*((j_end+3)-(j_sta-3)+1)*3*5
    allocate( s_jp3_r8( len_ik3_r8 ) ) 
    allocate( s_jm3_r8( len_ik3_r8 ) ) 
    allocate( s_kp3_r8( len_ij3_r8 ) ) 
    allocate( s_km3_r8( len_ij3_r8 ) ) 
    allocate( r_jp3_r8( len_ik3_r8 ) ) 
    allocate( r_jm3_r8( len_ik3_r8 ) ) 
    allocate( r_kp3_r8( len_ij3_r8 ) ) 
    allocate( r_km3_r8( len_ij3_r8 ) ) 

    call mpi_send_init(s_jm3_r8,len_ik3_r8,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,ireq3_r8(1),ierr)
    call mpi_send_init(s_jp3_r8,len_ik3_r8,MPI_DOUBLE_PRECISION,j_up  ,2,MPI_COMM_WORLD,ireq3_r8(2),ierr)
    call mpi_send_init(s_km3_r8,len_ij3_r8,MPI_DOUBLE_PRECISION,k_down,3,MPI_COMM_WORLD,ireq3_r8(3),ierr)
    call mpi_send_init(s_kp3_r8,len_ij3_r8,MPI_DOUBLE_PRECISION,k_up  ,4,MPI_COMM_WORLD,ireq3_r8(4),ierr)

    call mpi_recv_init(r_jp3_r8,len_ik3_r8,MPI_DOUBLE_PRECISION,j_up  ,1,MPI_COMM_WORLD,ireq3_r8(5),ierr)
    call mpi_recv_init(r_jm3_r8,len_ik3_r8,MPI_DOUBLE_PRECISION,j_down,2,MPI_COMM_WORLD,ireq3_r8(6),ierr)
    call mpi_recv_init(r_kp3_r8,len_ij3_r8,MPI_DOUBLE_PRECISION,k_up  ,3,MPI_COMM_WORLD,ireq3_r8(7),ierr)
    call mpi_recv_init(r_km3_r8,len_ij3_r8,MPI_DOUBLE_PRECISION,k_down,4,MPI_COMM_WORLD,ireq3_r8(8),ierr)













!!KTs
!    len_ik4 = (Nx+10)*((k_end+4)-(k_sta-4)+1)*4
!    len_ij4 = (Nx+10)*((j_end+4)-(j_sta-4)+1)*4
!    allocate( s_jp4( len_ik4 ) ) 
!    allocate( s_jm4( len_ik4 ) ) 
!    allocate( s_kp4( len_ij4 ) ) 
!    allocate( s_km4( len_ij4 ) ) 
!    allocate( r_jp4( len_ik4 ) ) 
!    allocate( r_jm4( len_ik4 ) ) 
!    allocate( r_kp4( len_ij4 ) ) 
!    allocate( r_km4( len_ij4 ) ) 
!
!    call mpi_send_init(s_jm4,len_ik4,MPI_REAL,j_down,1,MPI_COMM_WORLD,ireq4(1),ierr)
!    call mpi_send_init(s_jp4,len_ik4,MPI_REAL,j_up  ,2,MPI_COMM_WORLD,ireq4(2),ierr)
!    call mpi_send_init(s_km4,len_ij4,MPI_REAL,k_down,3,MPI_COMM_WORLD,ireq4(3),ierr)
!    call mpi_send_init(s_kp4,len_ij4,MPI_REAL,k_up  ,4,MPI_COMM_WORLD,ireq4(4),ierr)
!
!    call mpi_recv_init(r_jp4,len_ik4,MPI_REAL,j_up  ,1,MPI_COMM_WORLD,ireq4(5),ierr)
!    call mpi_recv_init(r_jm4,len_ik4,MPI_REAL,j_down,2,MPI_COMM_WORLD,ireq4(6),ierr)
!    call mpi_recv_init(r_kp4,len_ij4,MPI_REAL,k_up  ,3,MPI_COMM_WORLD,ireq4(7),ierr)
!    call mpi_recv_init(r_km4,len_ij4,MPI_REAL,k_down,4,MPI_COMM_WORLD,ireq4(8),ierr)
!!KTe


















  end subroutine mp_setup_division
!-----------------------------------------------------------------------------
  subroutine mp_process_check( ny , nz , ierr )
    implicit none
    integer,intent(in)    :: ny, nz
    integer,intent(inout) :: ierr
    integer :: i , j , k , jj , idum , jdum1 , jdum2 , kdum

    ierr = 0
    if( nprocs .le. 1 ) then ! parallel job
      ierr = 2
      return
    end if
    if( nprocs .le. 3 ) then ! at least, 2 x 2 = 4 processes
      ierr = 3
      return
    end if
!    i = int( log(dble(nprocs)) / log(2.0d0) )
!    j = int( (i+1) / 2.0 )
!    k = int( i / 2.0)
!    if( (2**j)*(2**k) .ne. nprocs ) then
!      ierr = 1
!      return
!    end if

    if(nprocs==4) then
    j=2
    k=2
    else
    idum=int(dsqrt(dble(nprocs)))
    kdum=nprocs
    do i=1,idum
     jdum1=i
     jdum2=nprocs/i
    if(jdum1*jdum2/=nprocs) then
      cycle
    else
     if(jdum1+jdum2<kdum) then
      kdum=jdum1+jdum2
      j=jdum1
      k=jdum2
     endif
    endif
    enddo
    endif

!    if( ny .lt. nz ) then
!      ndiv_nz = j
!      ndiv_ny = k
!    else
      ndiv_nz = k
      ndiv_ny = j
!    end if

    print*,'mp_prc_chk',j,k

    return



  end subroutine mp_process_check
!-----------------------------------------------------------------------------
  subroutine mp_stop( x )
    implicit none
    integer,intent(in) :: x
    integer :: ecode , ierr

    ecode = x
    if( ecode>9 ) then
       write(*,*) 'ERROR CODE=',x
       ecode = ecode - 10
       call MPI_ABORT( MPI_COMM_WORLD , ecode , ierr )
    end if

    call MPI_FINALIZE(ierr)

    select case(ecode)
       case(1)
          stop 1
       case(2)
          stop 2
       case default
          stop
    end select

  end subroutine mp_stop
!-----------------------------------------------------------------------------
  subroutine mp_barrier
    implicit none
    integer :: ierr

    call MPI_BARRIER( MPI_COMM_WORLD , ierr )

  end subroutine mp_barrier
!-----------------------------------------------------------------------------
  subroutine mp_bcast_i( x )
    implicit none
    integer,intent(inout) :: x
    integer :: ierr
    integer :: ibuff

    ibuff = x
    call MPI_BCAST( ibuff , 1 , MPI_INTEGER , &
                     root , MPI_COMM_WORLD , ierr )
    x = ibuff

  end subroutine mp_bcast_i
!-----------------------------------------------------------------------------
  subroutine mp_bcast_i1( n , x )
    implicit none
    integer,intent(in)    :: n
    integer,intent(inout) :: x(n)
    integer :: ierr
    integer :: ibuff(n)

    ibuff = x
    call MPI_BCAST( ibuff , n , MPI_INTEGER , &
                     root , MPI_COMM_WORLD , ierr )
    x = ibuff

  end subroutine mp_bcast_i1
!-----------------------------------------------------------------------------
  subroutine mp_bcast_i12( n1 , n2 , x )
    implicit none
    integer,intent(in)    :: n1,n2
    integer,intent(inout) :: x(n1:n2)
    integer :: nn , ierr
    integer :: ibuff(n1:n2)

    nn = (n2-n1+1)
    ibuff = x
    call MPI_BCAST( ibuff(n1) , nn , MPI_INTEGER , &
                         root , MPI_COMM_WORLD , ierr )
    x = ibuff

  end subroutine mp_bcast_i12
!-----------------------------------------------------------------------------
  subroutine mp_bcast_i2( n1, n2, x )
    implicit none
    integer,intent(in)    :: n1 , n2
    integer,intent(inout) :: x(n1,n2)
    integer :: ierr
    integer :: ibuff(n1,n2)

    ibuff = x
    call MPI_BCAST( ibuff , n1*n2 , MPI_INTEGER , &
                     root , MPI_COMM_WORLD , ierr )
    x = ibuff

  end subroutine mp_bcast_i2
!-----------------------------------------------------------------------------
  subroutine mp_bcast_i32( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    integer,intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    integer :: ibuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    ibuff = x
    call MPI_BCAST( ibuff(n11,n21,n31) , nnn , MPI_INTEGER , &
                                  root , MPI_COMM_WORLD , ierr )
    x = ibuff

  end subroutine mp_bcast_i32
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r( x )
    implicit none
    real,intent(inout) :: x
    integer :: ierr
    real :: rbuff

    rbuff = x
    call MPI_BCAST( rbuff , 1 , MPI_REAL , &
                     root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8( x )
    implicit none
    real(8),intent(inout) :: x
    integer :: ierr
    real(8) :: rbuff

    rbuff = x
    call MPI_BCAST( rbuff , 1 , MPI_DOUBLE_PRECISION , &
                     root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r8
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8_1( n , x )
    implicit none
    integer,intent(in) :: n
    real(8),intent(inout) :: x(n)
    integer :: ierr
    real(8) :: rbuff(n)

    rbuff = x
    call MPI_BCAST( rbuff , n , MPI_DOUBLE_PRECISION , &
                     root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r8_1
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8_12( n1 , n2 , x )
    implicit none
    integer,intent(in) :: n1,n2
    real(8),intent(inout) :: x(n1:n2)
    integer :: nn , ierr
    real(8) :: rbuff(n1:n2)

    nn = n2 - n1 + 1
    rbuff = x
    call MPI_BCAST( rbuff(n1) , nn , MPI_DOUBLE_PRECISION , &
                         root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r8_12
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8_2( n1 , n2 , x )
    implicit none
    integer,intent(in) :: n1 , n2
    real(8),intent(inout) :: x(n1,n2)
    integer :: ierr
    real(8) :: rbuff(n1,n2)

    rbuff = x
    call MPI_BCAST( rbuff , n1*n2 , MPI_DOUBLE_PRECISION , &
                     root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r8_2
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r22_myr( n1 , n2 , n3 , n4 , x )
    implicit none
    integer,intent(in) :: n1 , n2 , n3 , n4
    real(4),intent(inout) :: x(n1:n2,n3:n4)
    integer :: ierr , nnn
    real(4) :: rbuff(n1:n2,n3:n4)

    nnn   = (n2-n1+1)*(n4-n3+1)
    rbuff = x
    call MPI_BCAST( rbuff   , nnn , MPI_REAL , &
                    my_rank , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r22_myr
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8_3( n1 , n2 , n3 , x )
    implicit none
    integer,intent(in) :: n1 , n2 , n3
    real(8),intent(inout) :: x(n1,n2,n3)
    integer :: ierr
    real(8) :: rbuff(n1,n2,n3)

    rbuff = x
    call MPI_BCAST( rbuff , n1*n2*n3 , MPI_DOUBLE_PRECISION , &
                     root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r8_3
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r32( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(4),intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    real(4) :: rbuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    rbuff = x

!    call MPI_BCAST( rbuff(n11,n21,n31) , nnn , MPI_DOUBLE_PRECISION , &
    call MPI_BCAST( rbuff , nnn , MPI_REAL , &
                                  root , MPI_COMM_WORLD , ierr )

    if(ierr/=0) then
    print*,my_rank,nnn,n11,n21,n31
    call MPI_FINALIZE(ierr)
    endif

    x = rbuff

  end subroutine mp_bcast_r32
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r32_myr( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(4),intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    real(4) :: rbuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    rbuff = x
!    print*,'PPPP',my_rank,nnn
!    print*,'RRRR',my_rank,n12-n11,n22-n21,n32-n31
!    print*,'SSSS',my_rank,n22,n21
!    call MPI_BCAST( rbuff(n11,n21,n31) , nnn , MPI_DOUBLE_PRECISION , &
    call MPI_BCAST( rbuff , nnn , MPI_REAL , &
                                  my_rank , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r32_myr
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8_4( n1 , n2 , n3 , n4 , x )
    implicit none
    integer,intent(in) :: n1 , n2 , n3 , n4
    real(8),intent(inout) :: x(n1,n2,n3,n4)
    integer :: ierr
    real(8) :: rbuff(n1,n2,n3,n4)

    rbuff = x
    call MPI_BCAST( rbuff , n1*n2*n3*n4 , MPI_DOUBLE_PRECISION , &
                     root , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_bcast_r8_4
!-----------------------------------------------------------------------------
  subroutine mp_bcast_chr( n , x )
    implicit none
    integer,intent(in) :: n
    character(n),intent(inout) :: x
    integer :: ierr
    character(n) :: cbuff

    cbuff = x
    call MPI_BCAST( cbuff , n , MPI_CHARACTER , &
                     root , MPI_COMM_WORLD , ierr )
    x = cbuff

  end subroutine mp_bcast_chr
!-----------------------------------------------------------------------------
  subroutine mp_allmaxi( x )
    implicit none
    integer,intent(inout) :: x
    integer :: ierr
    integer :: io

    io = 0
    call MPI_ALLREDUCE( x , io , 1 , MPI_INTEGER , &
                  MPI_MAX , MPI_COMM_WORLD , ierr )
    x = io

  end subroutine mp_allmaxi
!-----------------------------------------------------------------------------
  subroutine mp_allmaxr( x )
    implicit none
    real(4),intent(inout) :: x
    integer :: ierr
    real(4) :: io

    io = 0.0d0
!    call MPI_ALLREDUCE( x , io , 1 , MPI_DOUBLE_PRECISION , &
    call MPI_ALLREDUCE( x , io , 1 , MPI_REAL , &
                  MPI_MAX , MPI_COMM_WORLD , ierr )
    x = io

  end subroutine mp_allmaxr
!-----------------------------------------------------------------------------
  subroutine mp_allmaxr8( x )
    implicit none
    real(8),intent(inout) :: x
    integer :: ierr
    real(8) :: io

    io = 0.0d0
    call MPI_ALLREDUCE( x , io , 1 , MPI_DOUBLE_PRECISION , &
!    call MPI_ALLREDUCE( x , io , 1 , MPI_REAL , &
                  MPI_MAX , MPI_COMM_WORLD , ierr )
    x = io

  end subroutine mp_allmaxr8
!-----------------------------------------------------------------------------
  subroutine mp_allminr( x )
    implicit none
    real(4),intent(inout) :: x
    integer :: ierr
    real(4) :: io

    io = 0.0d0
    call MPI_ALLREDUCE( x , io , 1 , MPI_REAL , &
                  MPI_MIN , MPI_COMM_WORLD , ierr )
    x = io

  end subroutine mp_allminr
!-----------------------------------------------------------------------------
  subroutine mp_allminr8( x )
    implicit none
    real(8),intent(inout) :: x
    integer :: ierr
    real(8) :: io

    io = 0.0d0
    call MPI_ALLREDUCE( x , io , 1 , MPI_DOUBLE_PRECISION , &
                  MPI_MIN , MPI_COMM_WORLD , ierr )
    x = io

  end subroutine mp_allminr8
!-----------------------------------------------------------------------------
  subroutine mp_allsumi( x )
    implicit none
    integer,intent(inout) :: x
    integer :: ierr
    integer :: ro

    ro = 0
    call MPI_ALLREDUCE( x , ro , 1 , MPI_INTEGER , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = ro

  end subroutine mp_allsumi
!-----------------------------------------------------------------------------
  subroutine mp_allsumr( x )
    implicit none
    real(4),intent(inout) :: x
    integer :: ierr
    real(4) :: ro

    ro = 0.0
    call MPI_ALLREDUCE( x , ro , 1 , MPI_REAL , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = ro

  end subroutine mp_allsumr
!-----------------------------------------------------------------------------
  subroutine mp_allsumr8( x )
    implicit none
    real(8),intent(inout) :: x
    integer :: ierr
    real(8) :: ro

    ro = 0.0
    call MPI_ALLREDUCE( x , ro , 1 , MPI_DOUBLE_PRECISION , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = ro

  end subroutine mp_allsumr8
!NECs
  subroutine mp_allsumr8_2( x,y )
    implicit none
    real(8),intent(inout) :: x,y
    integer :: ierr
    real(8) :: ro(2)

    ro(1) = x
    ro(2) = y
    call MPI_ALLREDUCE( MPI_IN_PLACE , ro, 2 , MPI_DOUBLE_PRECISION , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = ro(1)
    y = ro(2)

  end subroutine mp_allsumr8_2
!NECe

!KTs
  subroutine mp_allsumr8_3( x1,x2,x3 )
    implicit none
    real(8),intent(inout) :: x1,x2,x3
    integer :: ierr
    real(8) :: ro(3)

    ro(1) = x1
    ro(2) = x2
    ro(3) = x3
    call MPI_ALLREDUCE( MPI_IN_PLACE , ro, 3 , MPI_DOUBLE_PRECISION , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x1 = ro(1)
    x2 = ro(2)
    x3 = ro(3)

  end subroutine mp_allsumr8_3
!KTe
!-----------------------------------------------------------------------------
  subroutine mp_allsum_r32( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(4),intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    real(4) :: rbuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    rbuff = 0.0
    call MPI_ALLREDUCE(x , rbuff , nnn , MPI_REAL , &
               MPI_SUM , MPI_COMM_WORLD , ierr )

    x = rbuff

  end subroutine mp_allsum_r32
!NECs
  subroutine mp_allsum_r32_2( n11 , n12 , n21 , n22 , n31 , n32 , x , m)
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(4),intent(inout) :: x(n11:n12,n21:n22,n31:n32,m)
    integer,intent(in) :: m
    integer :: nnn , ierr

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)*m
    call MPI_ALLREDUCE(MPI_IN_PLACE , x , nnn , MPI_REAL , &
               MPI_SUM , MPI_COMM_WORLD , ierr )

  end subroutine mp_allsum_r32_2
!NECe
!-----------------------------------------------------------------------------
  subroutine mp_allor( x )
    implicit none
    logical,intent(inout) :: x
    integer :: ierr
    logical :: lo

    lo = .false.
    call MPI_ALLREDUCE( x , lo , 1 , MPI_LOGICAL , &
                  MPI_LOR , MPI_COMM_WORLD , ierr )
    x = lo

  end subroutine mp_allor
!-----------------------------------------------------------------------------
subroutine para_range ( nn1 , nn2 , ncpu , irnk , ista , iend )
  implicit none

  integer,intent(in)  :: nn1 , nn2 , ncpu , irnk
  integer,intent(out) :: ista, iend

  integer :: i1,i2

  i1 = ( nn2 - nn1 + 1 ) / ncpu
  i2 = mod( nn2 - nn1 + 1 , ncpu )
  ista = irnk * i1 + nn1 + min( irnk , i2 )
  iend = ista + i1 - 1
  if ( i2 > irnk ) iend = iend + 1

end subroutine para_range 
!-----------------------------------------------------------------------------

!=== TW ======================================================================
!-----------------------------------------------------------------------------
  subroutine mp_allsum_r12( n1 , n2 , x )
    implicit none
    integer,intent(in) :: n1,n2
    real(4),intent(inout) :: x(n1:n2)
    integer :: nn , ierr
    real(4) :: rbuff(n1:n2)

    nn = n2 - n1 + 1
    rbuff = 0.0
    call MPI_ALLREDUCE(x , rbuff , nn , MPI_REAL , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_allsum_r12
!=== TW ======================================================================
!-----------------------------------------------------------------------------
  subroutine mp_allsum_r8_12( n1 , n2 , x )
    implicit none
    integer,intent(in) :: n1,n2
    real(8),intent(inout) :: x(n1:n2)
    integer :: nn , ierr
    real(8) :: rbuff(n1:n2)

    nn = n2 - n1 + 1
    rbuff = 0.0d0
    call MPI_ALLREDUCE(x , rbuff , nn , MPI_DOUBLE_PRECISION , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_allsum_r8_12
!=== TW ======================================================================
!-----------------------------------------------------------------------------
  subroutine mp_allsum_i12( n1 , n2 , x )
    implicit none
    integer,intent(in) :: n1,n2
    integer,intent(inout) :: x(n1:n2)
    integer :: nn , ierr
    integer :: rbuff(n1:n2)

    nn = n2 - n1 + 1
    rbuff = 0
    call MPI_ALLREDUCE(x , rbuff , nn , MPI_INTEGER , &
                  MPI_SUM , MPI_COMM_WORLD , ierr )
    x = rbuff

  end subroutine mp_allsum_i12
!-----------------------------------------------------------------------------
  subroutine mp_sum_r22( n11 , n12 , n21 , n22 , x )
    implicit none
    integer,intent(in) :: n11 , n21
    integer,intent(in) :: n12 , n22
    real(4),intent(inout) :: x(n11:n12,n21:n22)
    integer :: nnn , ierr
    real(4) :: rbuff(n11:n12,n21:n22)

    nnn = (n12-n11+1)*(n22-n21+1)
    rbuff = 0.0
    call MPI_REDUCE(x , rbuff , nnn , MPI_REAL , &
               MPI_SUM , root, MPI_COMM_WORLD , ierr )

    If(my_rank==root) x = rbuff

  end subroutine mp_sum_r22
!-----------------------------------------------------------------------------
  subroutine mp_sum_r8_12( n11 , n12 , n21 , n22 , x )
    implicit none
    integer,intent(in) :: n11 , n21
    integer,intent(in) :: n12 , n22
    real(8),intent(inout) :: x(n11:n12,n21:n22)
    integer :: nnn , ierr
    real(8) :: rbuff(n11:n12,n21:n22)

    nnn = (n12-n11+1)*(n22-n21+1)
    rbuff = 0.0
    call MPI_REDUCE(x , rbuff , nnn , MPI_DOUBLE_PRECISION , &
               MPI_SUM , root, MPI_COMM_WORLD , ierr )

    If(my_rank==root) x = rbuff

  end subroutine mp_sum_r8_12
!-----------------------------------------------------------------------------
  subroutine mp_sum_r32( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(4),intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    real(4) :: rbuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    rbuff = 0.0
    call MPI_REDUCE(x , rbuff , nnn , MPI_REAL , &
               MPI_SUM , root, MPI_COMM_WORLD , ierr )

    If(my_rank==root) x = rbuff

  end subroutine mp_sum_r32
!-----------------------------------------------------------------------------
  subroutine mp_sum_r8_32( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(8),intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    real(8) :: rbuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    rbuff = 0.0
    call MPI_REDUCE(x , rbuff , nnn , MPI_DOUBLE_PRECISION , &
               MPI_SUM , root, MPI_COMM_WORLD , ierr )

    If(my_rank==root) x = rbuff

  end subroutine mp_sum_r8_32
!-----------------------------------------------------------------------------
  subroutine mp_bcast_r8_32( n11 , n12 , n21 , n22 , n31 , n32 , x )
    implicit none
    integer,intent(in) :: n11 , n21 , n31
    integer,intent(in) :: n12 , n22 , n32
    real(8),intent(inout) :: x(n11:n12,n21:n22,n31:n32)
    integer :: nnn , ierr
    real(8) :: rbuff(n11:n12,n21:n22,n31:n32)

    nnn = (n12-n11+1)*(n22-n21+1)*(n32-n31+1)
    rbuff = x

    call MPI_BCAST( rbuff , nnn , MPI_DOUBLE_PRECISION , &
                                  root , MPI_COMM_WORLD , ierr )

    if(ierr/=0) then
    print*,my_rank,nnn,n11,n21,n31
    call MPI_FINALIZE(ierr)
    endif

    x = rbuff

  end subroutine mp_bcast_r8_32
!-----------------------------------------------------------------------------
!
!--- MPI communication, sending and reciving message ---
!

subroutine mp_send_recv_pre(AA,imode,nx1,nx2,ny1,ny2,nz1,nz2)

  implicit none

  integer,intent(in) :: nx1,nx2,ny1,ny2,nz1,nz2
  real,intent(inout) :: AA(nx1:nx2,ny1:ny2,nz1:nz2)
  real               :: A1(nx1:nx2,j_sta-1:j_end+1,k_sta-1:k_end+1)
  real               :: A2(nx1:nx2,j_sta-2:j_end+2,k_sta-2:k_end+2)
  real               :: A3(nx1:nx2,j_sta-3:j_end+3,k_sta-3:k_end+3)
  real               :: A4(nx1:nx2,j_sta-4:j_end+4,k_sta-4:k_end+4)
  real               :: A5(nx1:nx2,j_sta-5:j_end+5,k_sta-5:k_end+5)
  integer            :: i,j,k,imode,kmax,kmin

   !----------------------------
   if(imode==1) then

!NECs
    inum=0
    do k=k_sta-1,k_end+1; do i=nx1,nx2
      inum=inum+1
      s_jm1(inum)=AA(i,j_sta,k)
      s_jp1(inum)=AA(i,j_end,k)
    end do; end do
    inum=0
    do j=j_sta-1,j_end+1; do i=nx1,nx2
      inum=inum+1
      s_km1(inum)=AA(i,j,k_sta)
      s_kp1(inum)=AA(i,j,k_end)
    end do; end do
    call mp_send_recv_1_r4_2
    if(j_up/=MPI_PROC_NULL) then
      inum=0
      do k=k_sta-1,k_end+1; do i=nx1,nx2
        inum=inum+1
        AA(i,j_end+1,k)=r_jp1(inum)
      end do; end do
    end if
    if(j_down/=MPI_PROC_NULL) then
      inum=0
      do k=k_sta-1,k_end+1; do i=nx1,nx2
        inum=inum+1
        AA(i,j_sta-1,k)=r_jm1(inum)
      end do; end do
    end if
    if(k_up/=MPI_PROC_NULL) then
      inum=0
      do j=j_sta-1,j_end+1; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k_end+1)=r_kp1(inum)
      end do; end do
    end if
    if(k_down/=MPI_PROC_NULL) then
      inum=0
      do j=j_sta-1,j_end+1; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k_sta-1)=r_km1(inum)
      end do; end do
    end if
!NECe
    !do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
    ! A1(i,j,k)=AA(i,j,k)
    !enddo; enddo; enddo
    !call mp_send_recv_1_r(A1,nx1,nx2,j_sta-1,j_end+1,k_sta-1,k_end+1)
    !do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
    ! AA(i,j,k)=A1(i,j,k)
    !enddo; enddo; enddo

   !----------------------------
   elseif(imode==2) then

    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     A2(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_2_r(A2,nx1,nx2,j_sta-2,j_end+2,k_sta-2,k_end+2)
    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     AA(i,j,k)=A2(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==3) then

!NECs
    inum=0
    do j=j_sta,j_sta+2; do k=k_sta-3,k_end+3; do i=nx1,nx2
      inum=inum+1
      s_jm3(inum)=AA(i,j,k)
    end do; end do; end do
    inum=0
    do j=j_end,j_end-2,-1; do k=k_sta-3,k_end+3; do i=nx1,nx2
      inum=inum+1
      s_jp3(inum)=AA(i,j,k)
    end do; end do; end do
    inum=0
    do k=k_sta,k_sta+2; do j=j_sta-3,j_end+3; do i=nx1,nx2
      inum=inum+1
      s_km3(inum)=AA(i,j,k)
    end do; end do; end do
    inum=0
    do k=k_end,k_end-2,-1; do j=j_sta-3,j_end+3; do i=nx1,nx2
      inum=inum+1
      s_kp3(inum)=AA(i,j,k)
    end do; end do; end do
    call mp_send_recv_3_r4_2
    if(j_up/=MPI_PROC_NULL) then
      inum=0
      do j=j_end+1,j_end+3; do k=k_sta-3,k_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_jp3(inum)
      end do; end do; end do
    end if
    if(j_down/=MPI_PROC_NULL) then
      inum=0
      do j=j_sta-1,j_sta-3,-1; do k=k_sta-3,k_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_jm3(inum)
      end do; end do; end do
    end if
    if(k_up/=MPI_PROC_NULL) then
      inum=0
      do k=k_end+1,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_kp3(inum)
      end do; end do; end do
    end if
    if(k_down/=MPI_PROC_NULL) then
      inum=0
      do k=k_sta-1,k_sta-3,-1; do j=j_sta-3,j_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_km3(inum)
      end do; end do; end do
    end if
!NECe
   ! do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
   !  A3(i,j,k)=AA(i,j,k)
   ! enddo; enddo; enddo
   ! call mp_send_recv_3_r(A3,nx1,nx2,j_sta-3,j_end+3,k_sta-3,k_end+3)
   ! do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
   !  AA(i,j,k)=A3(i,j,k)
   ! enddo; enddo; enddo

   !----------------------------
   elseif(imode==4) then

!KTs
    inum=0
    do j=j_sta,j_sta+3; do k=k_sta-4,k_end+4; do i=nx1,nx2
      inum=inum+1
      s_jm4(inum)=AA(i,j,k)
    end do; end do; end do
    inum=0
    do j=j_end,j_end-3,-1; do k=k_sta-4,k_end+4; do i=nx1,nx2
      inum=inum+1
      s_jp4(inum)=AA(i,j,k)
    end do; end do; end do
    inum=0
    do k=k_sta,k_sta+3; do j=j_sta-4,j_end+4; do i=nx1,nx2
      inum=inum+1
      s_km4(inum)=AA(i,j,k)
    end do; end do; end do
    inum=0
    do k=k_end,k_end-3,-1; do j=j_sta-4,j_end+4; do i=nx1,nx2
      inum=inum+1
      s_kp4(inum)=AA(i,j,k)
    end do; end do; end do
    call mp_send_recv_4_r4_2
    if(j_up/=MPI_PROC_NULL) then
      inum=0
      do j=j_end+1,j_end+4; do k=k_sta-4,k_end+4; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_jp4(inum)
      end do; end do; end do
    end if
    if(j_down/=MPI_PROC_NULL) then
      inum=0
      do j=j_sta-1,j_sta-4,-1; do k=k_sta-4,k_end+4; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_jm4(inum)
      end do; end do; end do
    end if
    if(k_up/=MPI_PROC_NULL) then
      inum=0
      do k=k_end+1,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_kp4(inum)
      end do; end do; end do
    end if
    if(k_down/=MPI_PROC_NULL) then
      inum=0
      do k=k_sta-1,k_sta-4,-1; do j=j_sta-4,j_end+4; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k)=r_km4(inum)
      end do; end do; end do
    end if
!KTe

    !do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
    ! A4(i,j,k)=AA(i,j,k)
    !enddo; enddo; enddo
    !call mp_send_recv_4_r(A4,nx1,nx2,j_sta-4,j_end+4,k_sta-4,k_end+4)
    !do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
    ! AA(i,j,k)=A4(i,j,k)
    !enddo; enddo; enddo

   !----------------------------
   elseif(imode==5) then

    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     A5(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_5_r(A5,nx1,nx2,j_sta-5,j_end+5,k_sta-5,k_end+5)
    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     AA(i,j,k)=A5(i,j,k)
    enddo; enddo; enddo

   endif

end subroutine mp_send_recv_pre

!
!-----------------------------------------------------------------------------
!
subroutine mp_send_recv_1_r(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype,dum_len,idum1
  real                  B1((nx2-nx1+1)*(ke-ks+1)*1),B2((nx2-nx1+1)*(ke-ks+1)*1)
  real                  B3((nx2-nx1+1)*(je-js+1)*1),B4((nx2-nx1+1)*(je-js+1)*1)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,js+1,j)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
if(j_up/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,je  ,j)=B2(idum1)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,je-1,j)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
if(j_down/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,js  ,j)=B2(idum1)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

! do i=1,Nlenx-1; do j=js,je; do k=ks,ke
!    A(i,j,k)=AA(i,j,k)
! enddo;enddo;enddo

!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ks+1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ke  )=B4(idum1)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ke-1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ks  )=B4(idum1)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_1_r
!NECs
subroutine mp_send_recv_1_r4_2

  implicit none
  integer iicom(8),ierr,istat(MPI_STATUS_SIZE,8)

  call mpi_startall(8,ireq1,ierr)
  call mpi_waitall(8,ireq1,istat,ierr)

end subroutine
!NECe
!
!-----------------------------------------------------------------------------
!
subroutine mp_send_recv_2_r(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype,dum_len,dum_len2,idum1,idum2
  real                  B1((nx2-nx1+1)*(ke-ks+1)*2),B2((nx2-nx1+1)*(ke-ks+1)*2)
  real                  B3((nx2-nx1+1)*(je-js+1)*2),B4((nx2-nx1+1)*(je-js+1)*2)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,js+2,j)
   B1(idum2)=A(i,js+3,j)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_up/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,je-1,j)=B2(idum1)
   AA(i,je  ,j)=B2(idum2)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,je-2,j)
   B1(idum2)=A(i,je-3,j)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_down/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,js+1,j)=B2(idum1)
   AA(i,js  ,j)=B2(idum2)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

! do i=1,Nlenx-1; do j=js,je; do k=ks,ke
!    A(i,j,k)=AA(i,j,k)
! enddo;enddo;enddo

!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ks+2)
   B3(idum2)=A(i,j,ks+3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ke-1)=B4(idum1)
   AA(i,j,ke  )=B4(idum2)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ke-2)
   B3(idum2)=A(i,j,ke-3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ks+1)=B4(idum1)
   AA(i,j,ks  )=B4(idum2)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_2_r
!
!-----------------------------------------------------------------------------
!
subroutine mp_send_recv_3_r(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3
  integer               idum1,idum2,idum3
  real                  B1((nx2-nx1+1)*(ke-ks+1)*3),B2((nx2-nx1+1)*(ke-ks+1)*3)
  real                  B3((nx2-nx1+1)*(je-js+1)*3),B4((nx2-nx1+1)*(je-js+1)*3)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,js+3,k)
   B1(idum2)=A(i,js+4,k)
   B1(idum3)=A(i,js+5,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,je-2,k)=B2(idum1)
   AA(i,je-1,k)=B2(idum2)
   AA(i,je  ,k)=B2(idum3)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,je-3,k)
   B1(idum2)=A(i,je-4,k)
   B1(idum3)=A(i,je-5,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,js+2,k)=B2(idum1)
   AA(i,js+1,k)=B2(idum2)
   AA(i,js  ,k)=B2(idum3)
  enddo
enddo
endif


!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ks+3)
   B3(idum2)=A(i,j,ks+4)
   B3(idum3)=A(i,j,ks+5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ke-2)=B4(idum1)
   AA(i,j,ke-1)=B4(idum2)
   AA(i,j,ke  )=B4(idum3)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ke-3)
   B3(idum2)=A(i,j,ke-4)
   B3(idum3)=A(i,j,ke-5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ks+2)=B4(idum1)
   AA(i,j,ks+1)=B4(idum2)
   AA(i,j,ks  )=B4(idum3)
  enddo
enddo
endif

return
end subroutine mp_send_recv_3_r
!NECs
subroutine mp_send_recv_3_r4_2

  implicit none
  integer iicom(8),ierr,istat(MPI_STATUS_SIZE,8)

  call mpi_startall(8,ireq3,ierr)
  call mpi_waitall(8,ireq3,istat,ierr)

end subroutine
!NECe
!NECs
subroutine mp_send_recv_3_r8_2

  implicit none
  integer iicom(8),ierr,istat(MPI_STATUS_SIZE,8)

  call mpi_startall(8,ireq3_r8,ierr)
  call mpi_waitall (8,ireq3_r8,istat,ierr)

end subroutine
!NECe
!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_4_r(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3,dum_len4
  integer               idum1,idum2,idum3,idum4
  real                  B1((nx2-nx1+1)*(ke-ks+1)*4),B2((nx2-nx1+1)*(ke-ks+1)*4)
  real                  B3((nx2-nx1+1)*(je-js+1)*4),B4((nx2-nx1+1)*(je-js+1)*4)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,js+4,k)
   B1(idum2)=A(i,js+5,k)
   B1(idum3)=A(i,js+6,k)
   B1(idum4)=A(i,js+7,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,je-3,k)=B2(idum1)
   AA(i,je-2,k)=B2(idum2)
   AA(i,je-1,k)=B2(idum3)
   AA(i,je  ,k)=B2(idum4)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,je-4,k)
   B1(idum2)=A(i,je-5,k)
   B1(idum3)=A(i,je-6,k)
   B1(idum4)=A(i,je-7,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,js+3,k)=B2(idum1)
   AA(i,js+2,k)=B2(idum2)
   AA(i,js+1,k)=B2(idum3)
   AA(i,js  ,k)=B2(idum4)
  enddo
enddo
endif


!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ks+4)
   B3(idum2)=A(i,j,ks+5)
   B3(idum3)=A(i,j,ks+6)
   B3(idum4)=A(i,j,ks+7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ke-3)=B4(idum1)
   AA(i,j,ke-2)=B4(idum2)
   AA(i,j,ke-1)=B4(idum3)
   AA(i,j,ke  )=B4(idum4)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ke-4)
   B3(idum2)=A(i,j,ke-5)
   B3(idum3)=A(i,j,ke-6)
   B3(idum4)=A(i,j,ke-7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ks+3)=B4(idum1)
   AA(i,j,ks+2)=B4(idum2)
   AA(i,j,ks+1)=B4(idum3)
   AA(i,j,ks  )=B4(idum4)
  enddo
enddo
endif

return
end subroutine mp_send_recv_4_r
!KTs
subroutine mp_send_recv_4_r4_2

  implicit none
  integer iicom(8),ierr,istat(MPI_STATUS_SIZE,8)

  call mpi_startall(8,ireq4,ierr)
  call mpi_waitall(8,ireq4,istat,ierr)

end subroutine
!KTe

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_5_r(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3,dum_len4,dum_len5
  integer               idum1,idum2,idum3,idum4,idum5
  real                  B1((nx2-nx1+1)*(ke-ks+1)*5),B2((nx2-nx1+1)*(ke-ks+1)*5)
  real                  B3((nx2-nx1+1)*(je-js+1)*5),B4((nx2-nx1+1)*(je-js+1)*5)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,js+5,k)
   B1(idum2)=A(i,js+6,k)
   B1(idum3)=A(i,js+7,k)
   B1(idum4)=A(i,js+8,k)
   B1(idum5)=A(i,js+9,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,je-4,k)=B2(idum1)
   AA(i,je-3,k)=B2(idum2)
   AA(i,je-2,k)=B2(idum3)
   AA(i,je-1,k)=B2(idum4)
   AA(i,je  ,k)=B2(idum5)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,je-5,k)
   B1(idum2)=A(i,je-6,k)
   B1(idum3)=A(i,je-7,k)
   B1(idum4)=A(i,je-8,k)
   B1(idum5)=A(i,je-9,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,js+4,k)=B2(idum1)
   AA(i,js+3,k)=B2(idum2)
   AA(i,js+2,k)=B2(idum3)
   AA(i,js+1,k)=B2(idum4)
   AA(i,js  ,k)=B2(idum5)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

! do i=1,Nlenx-1; do j=js,je; do k=ks,ke
!    A(i,j,k)=AA(i,j,k)
! enddo;enddo;enddo

!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ks+5)
   B3(idum2)=A(i,j,ks+6)
   B3(idum3)=A(i,j,ks+7)
   B3(idum4)=A(i,j,ks+8)
   B3(idum5)=A(i,j,ks+9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ke-4)=B4(idum1)
   AA(i,j,ke-3)=B4(idum2)
   AA(i,j,ke-2)=B4(idum3)
   AA(i,j,ke-1)=B4(idum4)
   AA(i,j,ke  )=B4(idum5)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ke-5)
   B3(idum2)=A(i,j,ke-6)
   B3(idum3)=A(i,j,ke-7)
   B3(idum4)=A(i,j,ke-8)
   B3(idum5)=A(i,j,ke-9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ks+4)=B4(idum1)
   AA(i,j,ks+3)=B4(idum2)
   AA(i,j,ks+2)=B4(idum3)
   AA(i,j,ks+1)=B4(idum4)
   AA(i,j,ks  )=B4(idum5)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_5_r

!=============================================================================
!=============================================================================
!=============================================================================

!
!--- MPI communication, sending and reciving message ---
!

subroutine mp_send_recv_pre_jk(AA,imode,idir,nx1,nx2,ny1,ny2,nz1,nz2)

  implicit none

  integer,intent(in) :: nx1,nx2,ny1,ny2,nz1,nz2
  real,intent(inout) :: AA(nx1:nx2,ny1:ny2,nz1:nz2)
  real               :: A1(nx1:nx2,j_sta-1:j_end+1,k_sta-1:k_end+1)
  real               :: A2(nx1:nx2,j_sta-2:j_end+2,k_sta-2:k_end+2)
  real               :: A3(nx1:nx2,j_sta-3:j_end+3,k_sta-3:k_end+3)
  real               :: A4(nx1:nx2,j_sta-4:j_end+4,k_sta-4:k_end+4)
  real               :: A5(nx1:nx2,j_sta-5:j_end+5,k_sta-5:k_end+5)
  integer            :: i,j,k,imode,idir,kmax,kmin

   !----------------------------
   if(imode==1) then

    do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
     A1(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_1_r(A1,idir,nx1,nx2,j_sta-1,j_end+1,k_sta-1,k_end+1)
    do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
     AA(i,j,k)=A1(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==2) then

    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     A2(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_2_r(A2,idir,nx1,nx2,j_sta-2,j_end+2,k_sta-2,k_end+2)
    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     AA(i,j,k)=A2(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==3) then

    do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
     A3(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_3_r(A3,idir,nx1,nx2,j_sta-3,j_end+3,k_sta-3,k_end+3)
    do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
     AA(i,j,k)=A3(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==4) then

    do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
     A4(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_4_r(A4,idir,nx1,nx2,j_sta-4,j_end+4,k_sta-4,k_end+4)
    do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
     AA(i,j,k)=A4(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==5) then

    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     A5(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_5_r(A5,idir,nx1,nx2,j_sta-5,j_end+5,k_sta-5,k_end+5)
    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     AA(i,j,k)=A5(i,j,k)
    enddo; enddo; enddo

   endif

end subroutine mp_send_recv_pre_jk

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_1_r(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len
  integer               idum1
  real                  B1((nx2-nx1+1)*(ke-ks+1)*1),B2((nx2-nx1+1)*(ke-ks+1)*1)
  real                  B3((nx2-nx1+1)*(je-js+1)*1),B4((nx2-nx1+1)*(je-js+1)*1)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,js+1,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,je  ,k)=B2(idum1)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,je-1,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,js  ,k)=B2(idum1)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ks+1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ke  )=B4(idum1)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ke-1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ks  )=B4(idum1)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_1_r

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_2_r(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2
  integer               idum1,idum2
  real                  B1((nx2-nx1+1)*(ke-ks+1)*2),B2((nx2-nx1+1)*(ke-ks+1)*2)
  real                  B3((nx2-nx1+1)*(je-js+1)*2),B4((nx2-nx1+1)*(je-js+1)*2)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,js+2,k)
   B1(idum2)=A(i,js+3,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,je-1,k)=B2(idum1)
   AA(i,je  ,k)=B2(idum2)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,je-2,k)
   B1(idum2)=A(i,je-3,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,js+1,k)=B2(idum1)
   AA(i,js  ,k)=B2(idum2)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ks+2)
   B3(idum2)=A(i,j,ks+3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ke-1)=B4(idum1)
   AA(i,j,ke  )=B4(idum2)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ke-2)
   B3(idum2)=A(i,j,ke-3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ks+1)=B4(idum1)
   AA(i,j,ks  )=B4(idum2)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_2_r

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_3_r(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3
  integer               idum1,idum2,idum3
  real                  B1((nx2-nx1+1)*(ke-ks+1)*3),B2((nx2-nx1+1)*(ke-ks+1)*3)
  real                  B3((nx2-nx1+1)*(je-js+1)*3),B4((nx2-nx1+1)*(je-js+1)*3)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,js+3,k)
   B1(idum2)=A(i,js+4,k)
   B1(idum3)=A(i,js+5,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,je-2,k)=B2(idum1)
   AA(i,je-1,k)=B2(idum2)
   AA(i,je  ,k)=B2(idum3)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,je-3,k)
   B1(idum2)=A(i,je-4,k)
   B1(idum3)=A(i,je-5,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,js+2,k)=B2(idum1)
   AA(i,js+1,k)=B2(idum2)
   AA(i,js  ,k)=B2(idum3)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ks+3)
   B3(idum2)=A(i,j,ks+4)
   B3(idum3)=A(i,j,ks+5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ke-2)=B4(idum1)
   AA(i,j,ke-1)=B4(idum2)
   AA(i,j,ke  )=B4(idum3)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ke-3)
   B3(idum2)=A(i,j,ke-4)
   B3(idum3)=A(i,j,ke-5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ks+2)=B4(idum1)
   AA(i,j,ks+1)=B4(idum2)
   AA(i,j,ks  )=B4(idum3)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_3_r

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_4_r(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3,dum_len4
  integer               idum1,idum2,idum3,idum4
  real                  B1((nx2-nx1+1)*(ke-ks+1)*4),B2((nx2-nx1+1)*(ke-ks+1)*4)
  real                  B3((nx2-nx1+1)*(je-js+1)*4),B4((nx2-nx1+1)*(je-js+1)*4)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,js+4,k)
   B1(idum2)=A(i,js+5,k)
   B1(idum3)=A(i,js+6,k)
   B1(idum4)=A(i,js+7,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,je-3,k)=B2(idum1)
   AA(i,je-2,k)=B2(idum2)
   AA(i,je-1,k)=B2(idum3)
   AA(i,je  ,k)=B2(idum4)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,je-4,k)
   B1(idum2)=A(i,je-5,k)
   B1(idum3)=A(i,je-6,k)
   B1(idum4)=A(i,je-7,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,js+3,k)=B2(idum1)
   AA(i,js+2,k)=B2(idum2)
   AA(i,js+1,k)=B2(idum3)
   AA(i,js  ,k)=B2(idum4)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ks+4)
   B3(idum2)=A(i,j,ks+5)
   B3(idum3)=A(i,j,ks+6)
   B3(idum4)=A(i,j,ks+7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ke-3)=B4(idum1)
   AA(i,j,ke-2)=B4(idum2)
   AA(i,j,ke-1)=B4(idum3)
   AA(i,j,ke  )=B4(idum4)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ke-4)
   B3(idum2)=A(i,j,ke-5)
   B3(idum3)=A(i,j,ke-6)
   B3(idum4)=A(i,j,ke-7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ks+3)=B4(idum1)
   AA(i,j,ks+2)=B4(idum2)
   AA(i,j,ks+1)=B4(idum3)
   AA(i,j,ks  )=B4(idum4)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_4_r

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_5_r(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3,dum_len4,dum_len5
  integer               idum1,idum2,idum3,idum4,idum5
  real                  B1((nx2-nx1+1)*(ke-ks+1)*5),B2((nx2-nx1+1)*(ke-ks+1)*5)
  real                  B3((nx2-nx1+1)*(je-js+1)*5),B4((nx2-nx1+1)*(je-js+1)*5)
  real                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,js+5,k)
   B1(idum2)=A(i,js+6,k)
   B1(idum3)=A(i,js+7,k)
   B1(idum4)=A(i,js+8,k)
   B1(idum5)=A(i,js+9,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,je-4,k)=B2(idum1)
   AA(i,je-3,k)=B2(idum2)
   AA(i,je-2,k)=B2(idum3)
   AA(i,je-1,k)=B2(idum4)
   AA(i,je  ,k)=B2(idum5)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,je-5,k)
   B1(idum2)=A(i,je-6,k)
   B1(idum3)=A(i,je-7,k)
   B1(idum4)=A(i,je-8,k)
   B1(idum5)=A(i,je-9,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_REAL,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_REAL,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,js+4,k)=B2(idum1)
   AA(i,js+3,k)=B2(idum2)
   AA(i,js+2,k)=B2(idum3)
   AA(i,js+1,k)=B2(idum4)
   AA(i,js  ,k)=B2(idum5)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ks+5)
   B3(idum2)=A(i,j,ks+6)
   B3(idum3)=A(i,j,ks+7)
   B3(idum4)=A(i,j,ks+8)
   B3(idum5)=A(i,j,ks+9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ke-4)=B4(idum1)
   AA(i,j,ke-3)=B4(idum2)
   AA(i,j,ke-2)=B4(idum3)
   AA(i,j,ke-1)=B4(idum4)
   AA(i,j,ke  )=B4(idum5)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ke-5)
   B3(idum2)=A(i,j,ke-6)
   B3(idum3)=A(i,j,ke-7)
   B3(idum4)=A(i,j,ke-8)
   B3(idum5)=A(i,j,ke-9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_REAL,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_REAL,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ks+4)=B4(idum1)
   AA(i,j,ks+3)=B4(idum2)
   AA(i,j,ks+2)=B4(idum3)
   AA(i,j,ks+1)=B4(idum4)
   AA(i,j,ks  )=B4(idum5)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_5_r

!
!*****************************************************************************
!
!=============================================================================
!=============================================================================
!=============================================================================



! --- Real*8 ---
!
!--- MPI communication, sending and reciving message ---
!

subroutine mp_send_recv_pre_r8(AA,imode,nx1,nx2,ny1,ny2,nz1,nz2)

  implicit none

  integer,intent(in)   :: nx1,nx2,ny1,ny2,nz1,nz2
  real*8,intent(inout) :: AA(nx1:nx2,ny1:ny2,nz1:nz2)
  real*8               :: A1(nx1:nx2,j_sta-1:j_end+1,k_sta-1:k_end+1)
  real*8               :: A2(nx1:nx2,j_sta-2:j_end+2,k_sta-2:k_end+2)
  real*8               :: A3(nx1:nx2,j_sta-3:j_end+3,k_sta-3:k_end+3)
  real*8               :: A4(nx1:nx2,j_sta-4:j_end+4,k_sta-4:k_end+4)
  real*8               :: A5(nx1:nx2,j_sta-5:j_end+5,k_sta-5:k_end+5)
  integer              :: i,j,k,imode,kmax,kmin

   !----------------------------
   if(imode==1) then
    do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
     A1(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_1_r8(A1,nx1,nx2,j_sta-1,j_end+1,k_sta-1,k_end+1)
    do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
     AA(i,j,k)=A1(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==2) then

    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     A2(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_2_r8(A2,nx1,nx2,j_sta-2,j_end+2,k_sta-2,k_end+2)
    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     AA(i,j,k)=A2(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==3) then

    do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
     A3(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_3_r8(A3,nx1,nx2,j_sta-3,j_end+3,k_sta-3,k_end+3)
    do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
     AA(i,j,k)=A3(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==4) then

    do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
     A4(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_4_r8(A4,nx1,nx2,j_sta-4,j_end+4,k_sta-4,k_end+4)
    do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
     AA(i,j,k)=A4(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==5) then

    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     A5(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_5_r8(A5,nx1,nx2,j_sta-5,j_end+5,k_sta-5,k_end+5)
    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     AA(i,j,k)=A5(i,j,k)
    enddo; enddo; enddo

   endif

end subroutine mp_send_recv_pre_r8

subroutine mp_send_recv_pre_r8_Vec(AA,imode,nx1,nx2,ny1,ny2,nz1,nz2)

  implicit none

  integer,intent(in)   :: nx1,nx2,ny1,ny2,nz1,nz2
  real*8,intent(inout) :: AA(nx1:nx2,ny1:ny2,nz1:nz2,5)
  real*8               :: A3(nx1:nx2,j_sta-3:j_end+3,k_sta-3:k_end+3,5)
  integer              :: i,j,k,l,imode,kmax,kmin,inum

    inum=0
    do l=1,5
    do j=j_sta,j_sta+2; do k=k_sta-3,k_end+3; do i=nx1,nx2
      inum=inum+1
      s_jm3_r8(inum)=AA(i,j,k,l)
    end do; end do; end do; end do

    inum=0
    do l=1,5
    do j=j_end,j_end-2,-1; do k=k_sta-3,k_end+3; do i=nx1,nx2
      inum=inum+1
      s_jp3_r8(inum)=AA(i,j,k,l)
    end do; end do; end do; end do

    inum=0
    do l=1,5
    do k=k_sta,k_sta+2; do j=j_sta-3,j_end+3; do i=nx1,nx2
      inum=inum+1
      s_km3_r8(inum)=AA(i,j,k,l)
    end do; end do; end do; end do

    inum=0
    do l=1,5
    do k=k_end,k_end-2,-1; do j=j_sta-3,j_end+3; do i=nx1,nx2
      inum=inum+1
      s_kp3_r8(inum)=AA(i,j,k,l)
    end do; end do; end do; end do

    call mp_send_recv_3_r8_2

    if(j_up/=MPI_PROC_NULL) then
      inum=0
      do l=1,5
      do j=j_end+1,j_end+3; do k=k_sta-3,k_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k,l)=r_jp3_r8(inum)
      end do; end do; end do; end do
    end if

    if(j_down/=MPI_PROC_NULL) then
      inum=0
      do l=1,5
      do j=j_sta-1,j_sta-3,-1; do k=k_sta-3,k_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k,l)=r_jm3_r8(inum)
      end do; end do; end do; end do
    end if

    if(k_up/=MPI_PROC_NULL) then
      inum=0
      do l=1,5
      do k=k_end+1,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k,l)=r_kp3_r8(inum)
      end do; end do; end do; end do
    end if

    if(k_down/=MPI_PROC_NULL) then
      inum=0
      do l=1,5
      do k=k_sta-1,k_sta-3,-1; do j=j_sta-3,j_end+3; do i=nx1,nx2
        inum=inum+1
        AA(i,j,k,l)=r_km3_r8(inum)
      end do; end do; end do; end do
    end if


end subroutine mp_send_recv_pre_r8_Vec


!
!-----------------------------------------------------------------------------
!
subroutine mp_send_recv_1_r8(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in)   :: nx1,nx2,js,je,ks,ke
  real*8,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer                 iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer                 jsize(2),jjsize(2),jstart(2),inewtype,dum_len,idum1
  real*8                  B1((nx2-nx1+1)*(ke-ks+1)*1),B2((nx2-nx1+1)*(ke-ks+1)*1)
  real*8                  B3((nx2-nx1+1)*(je-js+1)*1),B4((nx2-nx1+1)*(je-js+1)*1)
  real*8                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,js+1,j)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
if(j_up/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,je  ,j)=B2(idum1)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,je-1,j)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
if(j_down/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,js  ,j)=B2(idum1)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

! do i=1,Nlenx-1; do j=js,je; do k=ks,ke
!    A(i,j,k)=AA(i,j,k)
! enddo;enddo;enddo

!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ks+1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ke  )=B4(idum1)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ke-1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ks  )=B4(idum1)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_1_r8
!
!-----------------------------------------------------------------------------
!
subroutine mp_send_recv_2_r8(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in)   :: nx1,nx2,js,je,ks,ke
  real*8,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer                 iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer                 jsize(2),jjsize(2),jstart(2),inewtype,dum_len,dum_len2,idum1,idum2
  real*8                  B1((nx2-nx1+1)*(ke-ks+1)*2),B2((nx2-nx1+1)*(ke-ks+1)*2)
  real*8                  B3((nx2-nx1+1)*(je-js+1)*2),B4((nx2-nx1+1)*(je-js+1)*2)
  real*8                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,js+2,j)
   B1(idum2)=A(i,js+3,j)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_up/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,je-1,j)=B2(idum1)
   AA(i,je  ,j)=B2(idum2)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,je-2,j)
   B1(idum2)=A(i,je-3,j)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_down/=MPI_PROC_NULL) then
do j=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,js+1,j)=B2(idum1)
   AA(i,js  ,j)=B2(idum2)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

! do i=1,Nlenx-1; do j=js,je; do k=ks,ke
!    A(i,j,k)=AA(i,j,k)
! enddo;enddo;enddo

!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ks+2)
   B3(idum2)=A(i,j,ks+3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ke-1)=B4(idum1)
   AA(i,j,ke  )=B4(idum2)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ke-2)
   B3(idum2)=A(i,j,ke-3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ks+1)=B4(idum1)
   AA(i,j,ks  )=B4(idum2)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_2_r8
!
!-----------------------------------------------------------------------------
!
subroutine mp_send_recv_3_r8(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in)   :: nx1,nx2,js,je,ks,ke
  real*8,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer                 iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer                 jsize(2),jjsize(2),jstart(2),inewtype
  integer                 dum_len,dum_len2,dum_len3
  integer                 idum1,idum2,idum3
  real*8                  B1((nx2-nx1+1)*(ke-ks+1)*3),B2((nx2-nx1+1)*(ke-ks+1)*3)
  real*8                  B3((nx2-nx1+1)*(je-js+1)*3),B4((nx2-nx1+1)*(je-js+1)*3)
  real*8                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,js+3,k)
   B1(idum2)=A(i,js+4,k)
   B1(idum3)=A(i,js+5,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,je-2,k)=B2(idum1)
   AA(i,je-1,k)=B2(idum2)
   AA(i,je  ,k)=B2(idum3)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,je-3,k)
   B1(idum2)=A(i,je-4,k)
   B1(idum3)=A(i,je-5,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,js+2,k)=B2(idum1)
   AA(i,js+1,k)=B2(idum2)
   AA(i,js  ,k)=B2(idum3)
  enddo
enddo
endif


!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ks+3)
   B3(idum2)=A(i,j,ks+4)
   B3(idum3)=A(i,j,ks+5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ke-2)=B4(idum1)
   AA(i,j,ke-1)=B4(idum2)
   AA(i,j,ke  )=B4(idum3)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ke-3)
   B3(idum2)=A(i,j,ke-4)
   B3(idum3)=A(i,j,ke-5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ks+2)=B4(idum1)
   AA(i,j,ks+1)=B4(idum2)
   AA(i,j,ks  )=B4(idum3)
  enddo
enddo
endif

return
end subroutine mp_send_recv_3_r8
!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_4_r8(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in)   :: nx1,nx2,js,je,ks,ke
  real*8,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer                 iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer                 jsize(2),jjsize(2),jstart(2),inewtype
  integer                 dum_len,dum_len2,dum_len3,dum_len4
  integer                 idum1,idum2,idum3,idum4
  real*8                  B1((nx2-nx1+1)*(ke-ks+1)*4),B2((nx2-nx1+1)*(ke-ks+1)*4)
  real*8                  B3((nx2-nx1+1)*(je-js+1)*4),B4((nx2-nx1+1)*(je-js+1)*4)
  real*8                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,js+4,k)
   B1(idum2)=A(i,js+5,k)
   B1(idum3)=A(i,js+6,k)
   B1(idum4)=A(i,js+7,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,je-3,k)=B2(idum1)
   AA(i,je-2,k)=B2(idum2)
   AA(i,je-1,k)=B2(idum3)
   AA(i,je  ,k)=B2(idum4)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,je-4,k)
   B1(idum2)=A(i,je-5,k)
   B1(idum3)=A(i,je-6,k)
   B1(idum4)=A(i,je-7,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,js+3,k)=B2(idum1)
   AA(i,js+2,k)=B2(idum2)
   AA(i,js+1,k)=B2(idum3)
   AA(i,js  ,k)=B2(idum4)
  enddo
enddo
endif


!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ks+4)
   B3(idum2)=A(i,j,ks+5)
   B3(idum3)=A(i,j,ks+6)
   B3(idum4)=A(i,j,ks+7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ke-3)=B4(idum1)
   AA(i,j,ke-2)=B4(idum2)
   AA(i,j,ke-1)=B4(idum3)
   AA(i,j,ke  )=B4(idum4)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ke-4)
   B3(idum2)=A(i,j,ke-5)
   B3(idum3)=A(i,j,ke-6)
   B3(idum4)=A(i,j,ke-7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ks+3)=B4(idum1)
   AA(i,j,ks+2)=B4(idum2)
   AA(i,j,ks+1)=B4(idum3)
   AA(i,j,ks  )=B4(idum4)
  enddo
enddo
endif

return
end subroutine mp_send_recv_4_r8
!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_5_r8(AA,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in)   :: nx1,nx2,js,je,ks,ke
  real*8,intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer                 iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer                 jsize(2),jjsize(2),jstart(2),inewtype
  integer                 dum_len,dum_len2,dum_len3,dum_len4,dum_len5
  integer                 idum1,idum2,idum3,idum4,idum5
  real*8                  B1((nx2-nx1+1)*(ke-ks+1)*5),B2((nx2-nx1+1)*(ke-ks+1)*5)
  real*8                  B3((nx2-nx1+1)*(je-js+1)*5),B4((nx2-nx1+1)*(je-js+1)*5)
  real*8                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,js+5,k)
   B1(idum2)=A(i,js+6,k)
   B1(idum3)=A(i,js+7,k)
   B1(idum4)=A(i,js+8,k)
   B1(idum5)=A(i,js+9,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,je-4,k)=B2(idum1)
   AA(i,je-3,k)=B2(idum2)
   AA(i,je-2,k)=B2(idum3)
   AA(i,je-1,k)=B2(idum4)
   AA(i,je  ,k)=B2(idum5)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,je-5,k)
   B1(idum2)=A(i,je-6,k)
   B1(idum3)=A(i,je-7,k)
   B1(idum4)=A(i,je-8,k)
   B1(idum5)=A(i,je-9,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,js+4,k)=B2(idum1)
   AA(i,js+3,k)=B2(idum2)
   AA(i,js+2,k)=B2(idum3)
   AA(i,js+1,k)=B2(idum4)
   AA(i,js  ,k)=B2(idum5)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

! do i=1,Nlenx-1; do j=js,je; do k=ks,ke
!    A(i,j,k)=AA(i,j,k)
! enddo;enddo;enddo

!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ks+5)
   B3(idum2)=A(i,j,ks+6)
   B3(idum3)=A(i,j,ks+7)
   B3(idum4)=A(i,j,ks+8)
   B3(idum5)=A(i,j,ks+9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ke-4)=B4(idum1)
   AA(i,j,ke-3)=B4(idum2)
   AA(i,j,ke-2)=B4(idum3)
   AA(i,j,ke-1)=B4(idum4)
   AA(i,j,ke  )=B4(idum5)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ke-5)
   B3(idum2)=A(i,j,ke-6)
   B3(idum3)=A(i,j,ke-7)
   B3(idum4)=A(i,j,ke-8)
   B3(idum5)=A(i,j,ke-9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ks+4)=B4(idum1)
   AA(i,j,ks+3)=B4(idum2)
   AA(i,j,ks+2)=B4(idum3)
   AA(i,j,ks+1)=B4(idum4)
   AA(i,j,ks  )=B4(idum5)
  enddo
enddo
endif

! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_5_r8

!
!*****************************************************************************
!
!
!--- MPI communication, sending and reciving message ---
!

subroutine mp_send_recv_pre_r8_jk(AA,imode,idir,nx1,nx2,ny1,ny2,nz1,nz2)

  implicit none

  integer,intent(in) :: nx1,nx2,ny1,ny2,nz1,nz2
  real(8),intent(inout) :: AA(nx1:nx2,ny1:ny2,nz1:nz2)
  real(8)               :: A1(nx1:nx2,j_sta-1:j_end+1,k_sta-1:k_end+1)
  real(8)               :: A2(nx1:nx2,j_sta-2:j_end+2,k_sta-2:k_end+2)
  real(8)               :: A3(nx1:nx2,j_sta-3:j_end+3,k_sta-3:k_end+3)
  real(8)               :: A4(nx1:nx2,j_sta-4:j_end+4,k_sta-4:k_end+4)
  real(8)               :: A5(nx1:nx2,j_sta-5:j_end+5,k_sta-5:k_end+5)
  integer            :: i,j,k,imode,idir,kmax,kmin

   !----------------------------
   if(imode==1) then

    do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
     A1(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_1_r8(A1,idir,nx1,nx2,j_sta-1,j_end+1,k_sta-1,k_end+1)
    do k=k_sta-1,k_end+1; do j=j_sta-1,j_end+1; do i=nx1,nx2
     AA(i,j,k)=A1(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==2) then

    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     A2(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_2_r8(A2,idir,nx1,nx2,j_sta-2,j_end+2,k_sta-2,k_end+2)
    do k=k_sta-2,k_end+2; do j=j_sta-2,j_end+2; do i=nx1,nx2
     AA(i,j,k)=A2(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==3) then

    do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
     A3(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_3_r8(A3,idir,nx1,nx2,j_sta-3,j_end+3,k_sta-3,k_end+3)
    do k=k_sta-3,k_end+3; do j=j_sta-3,j_end+3; do i=nx1,nx2
     AA(i,j,k)=A3(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==4) then

    do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
     A4(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_4_r8(A4,idir,nx1,nx2,j_sta-4,j_end+4,k_sta-4,k_end+4)
    do k=k_sta-4,k_end+4; do j=j_sta-4,j_end+4; do i=nx1,nx2
     AA(i,j,k)=A4(i,j,k)
    enddo; enddo; enddo

   !----------------------------
   elseif(imode==5) then

    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     A5(i,j,k)=AA(i,j,k)
    enddo; enddo; enddo
    call mp_send_recv_jk_5_r8(A5,idir,nx1,nx2,j_sta-5,j_end+5,k_sta-5,k_end+5)
    do k=k_sta-5,k_end+5; do j=j_sta-5,j_end+5; do i=nx1,nx2
     AA(i,j,k)=A5(i,j,k)
    enddo; enddo; enddo

   endif

end subroutine mp_send_recv_pre_r8_jk

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_1_r8(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real(8),intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len
  integer               idum1
  real(8)               B1((nx2-nx1+1)*(ke-ks+1)*1),B2((nx2-nx1+1)*(ke-ks+1)*1)
  real(8)               B3((nx2-nx1+1)*(je-js+1)*1),B4((nx2-nx1+1)*(je-js+1)*1)
  real(8)               A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,js+1,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,je  ,k)=B2(idum1)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
    idum1=0
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   B1(idum1)=A(i,je-1,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,js  ,k)=B2(idum1)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ks+1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ke  )=B4(idum1)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
    idum1=0
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   B3(idum1)=A(i,j,ke-1)
  enddo
enddo
  call mpi_isend(B3(1),dum_len,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   AA(i,j,ks  )=B4(idum1)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_1_r8

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_2_r8(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real(8),intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2
  integer               idum1,idum2
  real(8)                  B1((nx2-nx1+1)*(ke-ks+1)*2),B2((nx2-nx1+1)*(ke-ks+1)*2)
  real(8)                  B3((nx2-nx1+1)*(je-js+1)*2),B4((nx2-nx1+1)*(je-js+1)*2)
  real(8)                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,js+2,k)
   B1(idum2)=A(i,js+3,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,je-1,k)=B2(idum1)
   AA(i,je  ,k)=B2(idum2)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B1(idum1)=A(i,je-2,k)
   B1(idum2)=A(i,je-3,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len2,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len2,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,js+1,k)=B2(idum1)
   AA(i,js  ,k)=B2(idum2)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ks+2)
   B3(idum2)=A(i,j,ks+3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ke-1)=B4(idum1)
   AA(i,j,ke  )=B4(idum2)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
    idum1=0
    idum2=dum_len
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   B3(idum1)=A(i,j,ke-2)
   B3(idum2)=A(i,j,ke-3)
  enddo
enddo
  call mpi_isend(B3(1),dum_len2,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len2,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   AA(i,j,ks+1)=B4(idum1)
   AA(i,j,ks  )=B4(idum2)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_2_r8

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_3_r8(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real(8),intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3
  integer               idum1,idum2,idum3
  real(8)                  B1((nx2-nx1+1)*(ke-ks+1)*3),B2((nx2-nx1+1)*(ke-ks+1)*3)
  real(8)                  B3((nx2-nx1+1)*(je-js+1)*3),B4((nx2-nx1+1)*(je-js+1)*3)
  real(8)                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,js+3,k)
   B1(idum2)=A(i,js+4,k)
   B1(idum3)=A(i,js+5,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,je-2,k)=B2(idum1)
   AA(i,je-1,k)=B2(idum2)
   AA(i,je  ,k)=B2(idum3)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B1(idum1)=A(i,je-3,k)
   B1(idum2)=A(i,je-4,k)
   B1(idum3)=A(i,je-5,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len3,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len3,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,js+2,k)=B2(idum1)
   AA(i,js+1,k)=B2(idum2)
   AA(i,js  ,k)=B2(idum3)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ks+3)
   B3(idum2)=A(i,j,ks+4)
   B3(idum3)=A(i,j,ks+5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ke-2)=B4(idum1)
   AA(i,j,ke-1)=B4(idum2)
   AA(i,j,ke  )=B4(idum3)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
    idum1=0
    idum2=dum_len
    idum3=dum_len2
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   B3(idum1)=A(i,j,ke-3)
   B3(idum2)=A(i,j,ke-4)
   B3(idum3)=A(i,j,ke-5)
  enddo
enddo
  call mpi_isend(B3(1),dum_len3,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len3,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   AA(i,j,ks+2)=B4(idum1)
   AA(i,j,ks+1)=B4(idum2)
   AA(i,j,ks  )=B4(idum3)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_3_r8

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_4_r8(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real(8),intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3,dum_len4
  integer               idum1,idum2,idum3,idum4
  real(8)                  B1((nx2-nx1+1)*(ke-ks+1)*4),B2((nx2-nx1+1)*(ke-ks+1)*4)
  real(8)                  B3((nx2-nx1+1)*(je-js+1)*4),B4((nx2-nx1+1)*(je-js+1)*4)
  real(8)                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,js+4,k)
   B1(idum2)=A(i,js+5,k)
   B1(idum3)=A(i,js+6,k)
   B1(idum4)=A(i,js+7,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,je-3,k)=B2(idum1)
   AA(i,je-2,k)=B2(idum2)
   AA(i,je-1,k)=B2(idum3)
   AA(i,je  ,k)=B2(idum4)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B1(idum1)=A(i,je-4,k)
   B1(idum2)=A(i,je-5,k)
   B1(idum3)=A(i,je-6,k)
   B1(idum4)=A(i,je-7,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len4,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len4,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,js+3,k)=B2(idum1)
   AA(i,js+2,k)=B2(idum2)
   AA(i,js+1,k)=B2(idum3)
   AA(i,js  ,k)=B2(idum4)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ks+4)
   B3(idum2)=A(i,j,ks+5)
   B3(idum3)=A(i,j,ks+6)
   B3(idum4)=A(i,j,ks+7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ke-3)=B4(idum1)
   AA(i,j,ke-2)=B4(idum2)
   AA(i,j,ke-1)=B4(idum3)
   AA(i,j,ke  )=B4(idum4)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   B3(idum1)=A(i,j,ke-4)
   B3(idum2)=A(i,j,ke-5)
   B3(idum3)=A(i,j,ke-6)
   B3(idum4)=A(i,j,ke-7)
  enddo
enddo
  call mpi_isend(B3(1),dum_len4,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len4,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   AA(i,j,ks+3)=B4(idum1)
   AA(i,j,ks+2)=B4(idum2)
   AA(i,j,ks+1)=B4(idum3)
   AA(i,j,ks  )=B4(idum4)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_4_r8

!
!-----------------------------------------------------------------------------
!

subroutine mp_send_recv_jk_5_r8(AA,idir,nx1,nx2,js,je,ks,ke)

  implicit none

  integer,intent(in) :: nx1,nx2,js,je,ks,ke,idir
  real(8),intent(inout) :: AA(nx1:nx2,js:je,ks:ke)
  integer               iicom(8),ierr,istatus(MPI_STATUS_SIZE),i,j,k
  integer               jsize(2),jjsize(2),jstart(2),inewtype
  integer               dum_len,dum_len2,dum_len3,dum_len4,dum_len5
  integer               idum1,idum2,idum3,idum4,idum5
  real(8)                  B1((nx2-nx1+1)*(ke-ks+1)*5),B2((nx2-nx1+1)*(ke-ks+1)*5)
  real(8)                  B3((nx2-nx1+1)*(je-js+1)*5),B4((nx2-nx1+1)*(je-js+1)*5)
  real(8)                  A(nx1:nx2,js:je,ks:ke)

 do k=ks,ke; do j=js,je; do i=nx1,nx2
    A(i,j,k)=AA(i,j,k)
 enddo;enddo;enddo

if(idir==2) then
!-- comm. j-1 dir.

! call MPI_Barrier(MPI_COMM_WORLD,ierr)

   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,js+5,k)
   B1(idum2)=A(i,js+6,k)
   B1(idum3)=A(i,js+7,k)
   B1(idum4)=A(i,js+8,k)
   B1(idum5)=A(i,js+9,k)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(1),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(3),ierr)
  call mpi_wait(iicom(1),istatus,ierr)
  call mpi_wait(iicom(3),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_up/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,je-4,k)=B2(idum1)
   AA(i,je-3,k)=B2(idum2)
   AA(i,je-2,k)=B2(idum3)
   AA(i,je-1,k)=B2(idum4)
   AA(i,je  ,k)=B2(idum5)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_down,idum1,B1(idum1)
  enddo
enddo
endif

!-- comm. j+1 dir.
   dum_len=(nx2-nx1+1)*(ke-ks+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B1(idum1)=A(i,je-5,k)
   B1(idum2)=A(i,je-6,k)
   B1(idum3)=A(i,je-7,k)
   B1(idum4)=A(i,je-8,k)
   B1(idum5)=A(i,je-9,k)
!   if(idum1==8290) print*,'LLLLLL',my_rank,j_up,idum1,B1(idum1)
  enddo
enddo
  call mpi_isend(B1(1),dum_len5,MPI_DOUBLE_PRECISION,j_up,1,MPI_COMM_WORLD,iicom(2),ierr)
  call mpi_irecv(B2(1),dum_len5,MPI_DOUBLE_PRECISION,j_down,1,MPI_COMM_WORLD,iicom(4),ierr)
  call mpi_wait(iicom(2),istatus,ierr)
  call mpi_wait(iicom(4),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(j_down/=MPI_PROC_NULL) then
do k=ks,ke
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,js+4,k)=B2(idum1)
   AA(i,js+3,k)=B2(idum2)
   AA(i,js+2,k)=B2(idum3)
   AA(i,js+1,k)=B2(idum4)
   AA(i,js  ,k)=B2(idum5)
  enddo
enddo
endif
endif

if(idir==3) then
!-- comm. k-1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ks+5)
   B3(idum2)=A(i,j,ks+6)
   B3(idum3)=A(i,j,ks+7)
   B3(idum4)=A(i,j,ks+8)
   B3(idum5)=A(i,j,ks+9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(5),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(8),ierr)
  call mpi_wait(iicom(5),istatus,ierr)
  call mpi_wait(iicom(8),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_up/=MPI_PROC_NULL) then
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ke-4)=B4(idum1)
   AA(i,j,ke-3)=B4(idum2)
   AA(i,j,ke-2)=B4(idum3)
   AA(i,j,ke-1)=B4(idum4)
   AA(i,j,ke  )=B4(idum5)
  enddo
enddo
 endif

!-- comm. k+1 dir.
   dum_len=(nx2-nx1+1)*(je-js+1)
   dum_len2=dum_len*2
   dum_len3=dum_len*3
   dum_len4=dum_len*4
   dum_len5=dum_len*5
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
do j=js,je
  do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   B3(idum1)=A(i,j,ke-5)
   B3(idum2)=A(i,j,ke-6)
   B3(idum3)=A(i,j,ke-7)
   B3(idum4)=A(i,j,ke-8)
   B3(idum5)=A(i,j,ke-9)
  enddo
enddo
  call mpi_isend(B3(1),dum_len5,MPI_DOUBLE_PRECISION,k_up,1,MPI_COMM_WORLD,iicom(7),ierr)
  call mpi_irecv(B4(1),dum_len5,MPI_DOUBLE_PRECISION,k_down,1,MPI_COMM_WORLD,iicom(6),ierr)
  call mpi_wait(iicom(7),istatus,ierr)
  call mpi_wait(iicom(6),istatus,ierr)
    idum1=0
    idum2=dum_len
    idum3=dum_len2
    idum4=dum_len3
    idum5=dum_len4
if(k_down/=MPI_PROC_NULL) then
do j=js,je
 do i=nx1,nx2
   idum1=idum1+1
   idum2=idum2+1
   idum3=idum3+1
   idum4=idum4+1
   idum5=idum5+1
   AA(i,j,ks+4)=B4(idum1)
   AA(i,j,ks+3)=B4(idum2)
   AA(i,j,ks+2)=B4(idum3)
   AA(i,j,ks+1)=B4(idum4)
   AA(i,j,ks  )=B4(idum5)
  enddo
enddo
endif
endif
! do i=1,Nlenx-1; do j=js-2,je+2; do k=ks-2,ke+2
!    AA(i,j,k)=A(i,j,k)
! enddo;enddo;enddo

return
end subroutine mp_send_recv_jk_5_r8

!
!*****************************************************************************
!

!-----------------------------------------------------------------------------

end module module_mpi
!-----------------------------------------------------------------------------
!=============================================================================
!==== [  ctime  ] ============================================================
!=============================================================================
! common ctime stores the present date and time.                              
!   use VMS option to link this routine.                                      
!-----------------------------------------------------------------------------
      Module ctime

        Implicit none

        Character(8)  :: adate
        Character(10) :: atime
        Real :: cput1 , cput2

        !Common / ctime0 / adate , atime
      End Module ctime
!=============================================================================
!==== [  cunit  ] ============================================================
!=============================================================================
! common for unit number.                                                     
!-----------------------------------------------------------------------------
      Module cunit

        Implicit none

!     inpunt = unit for input file.
!     iwkunt = work file for input.
!              input data are copied to this unit.
!     iotunt = unit for output file.
!     imsunt = unit for message file.
!     iplunt = unit for plot file.
!     idmunt = unit for dump file.
!     irsunt = unit for restart file.
!     igmunt = unit for geometry file
!     ivpunt = unit for initial data file

        Integer,parameter :: inpunt = 5
        Integer,parameter :: iwkunt = 55
        Integer,parameter :: iotunt = 9
        Integer,parameter :: imsunt = 0
        Integer,parameter :: iplunt = 10
        Integer,parameter :: idmunt = 11
        Integer,parameter :: irsunt = 12
        Integer,parameter :: igmunt = 21
        Integer,parameter :: ivpunt = 22

!      common / cunit  / inpunt , iwkunt , iotunt , imsunt , iplunt ,
!     .                  idmunt , irsunt , igmunt , ivpunt
       End module cunit

!-----------------------------------------------------------------------------
! subroutine error prints error message.
!-----------------------------------------------------------------------------
subroutine error ( leng , mess , icode  )

  use module_mpi
  use ctime,only : cput1,cput2
  use cunit,only : iotunt,imsunt

  implicit none

  integer,intent(in)      :: leng , icode
  character(*),intent(in) :: mess

!    leng   = character length of string mess.
!    mess   = error message
!    icode  = error code
!             0      : no error
!             4      : warnig. continue execution.
!             8      : error. stop execution.

  integer :: length

  length = leng
  if ( length .le. 0 ) length = len(mess)

  if ( icode .eq. 0 ) then
     if( my_rank .eq. root ) then
        write(iotunt,6900) mess(1:length)
        write(imsunt,6900) mess(1:length)
     end if
  else if ( icode .eq. 4 ) then
     if( my_rank .eq. root ) then
        write(iotunt,6100 )
        write(iotunt,6900 ) mess(1:length)
        write(iotunt,6200 )
        write(imsunt,6100 )
        write(imsunt,6900 ) mess(1:length)
        write(imsunt,6200 )
     end if
  else if ( icode .eq. 8 ) then
     write(iotunt,6300 )
     write(iotunt,6900 ) mess(1:length)
     write(iotunt,6400 )
     write(imsunt,6300 )
     write(imsunt,6900 ) mess(1:length)
     write(imsunt,6400 )
     call cpu_time ( cput2 )
     write(iotunt,7000) cput2 - cput1
     write(imsunt,7000) cput2 - cput1
     call mp_stop( 11 )
  else
     call cpu_time ( cput2 )
     write(iotunt,9000) icode
     write(iotunt,7000) cput2 - cput1
     write(imsunt,9000) icode
     write(imsunt,7000) cput2 - cput1
     call mp_stop( 12 )
  end if

  return

!--- format

6100 format(/10x,'***********'/&
             10x,'* warning *'/&
             10x,'***********')
6200 format(/60x,'job continues')
6300 format(/10x,'***************'/&
             10x,'* fatal error *'/&
             10x,'***************')
6400 format(/60x,'job aborts')
6900 format(/10x,a)
7000 format(10x,'cputime = ' , f15.6 , ' sec ')
9000 format(1x ,'*error* icode not defined. icode =' , i5 )

end subroutine error

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
