program main
  use module_mpi
  use mod_precision, only : dp
  use mod_constants, only : pi
  implicit none
  !integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nv = 5

  ! ------------------------------
  ! Simulation parameters
  ! ------------------------------
  integer, parameter :: nx =  64
  integer, parameter :: ny =  64
  integer, parameter :: nz =  64
  integer, parameter :: nghost = 3         ! WENO5 needs 3 ghosts
  !integer, parameter :: mode = 0           ! mode=0:constant temperature,  mode=1:constant density

  !real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: x_min = 0.0_dp
  real(dp), parameter :: x_max = 2.0_dp*pi
  real(dp), parameter :: x_center = (x_max-x_min)/2.0_dp
  real(dp), parameter :: y_min = 0.0_dp
  real(dp), parameter :: y_max = 2.0_dp*pi
  real(dp), parameter :: y_center = (y_max-y_min)/2.0_dp
  real(dp), parameter :: z_min = 0.0_dp
  real(dp), parameter :: z_max = 2.0_dp*pi
  real(dp), parameter :: z_center = (z_max-z_min)/2.0_dp

  real(dp), parameter :: t_max = 100.0_dp
  real(dp), parameter :: cfl   = 0.50_dp   ! a bit safer for sharper capture
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: small_rho = 1.0d-12
  real(dp), parameter :: small_p   = 1.0d-12
  integer, parameter :: output_frequency =  10 ! Output every 20 steps for Paraview

  ! ------------------------------
  ! Grid / geometry (3D)
  ! ------------------------------
  real(dp), allocatable :: x_edge(:,:,:)      ! faces: [1:2,0..nx, 0..ny]
  real(dp), allocatable :: x_cell(:,:,:)      ! centers: 1..nx
  real(dp), allocatable :: y_edge(:,:,:)      ! faces: [1:2,0..nx, 0..ny]
  real(dp), allocatable :: y_cell(:,:,:)      ! centers: 1..nx
  real(dp), allocatable :: z_edge(:,:,:)      ! faces: [1:2,0..nx, 0..ny]
  real(dp), allocatable :: z_cell(:,:,:)      ! centers: 1..nx
  real(dp), allocatable :: vol   (:,:,:)      ! cell volumes (length in 1D)
  real(dp), allocatable :: area_x(:,:,:)      ! face areas: 0..nx
  real(dp), allocatable :: area_y(:,:,:)      ! face areas: 0..nx
  real(dp), allocatable :: area_z(:,:,:)      ! face areas: 0..nx

  ! ------------------------------
  ! Conserved variables & work
  ! ------------------------------
  real(dp), allocatable :: Q  (:,:,:,:)              ! [1-nghost:nx+nghost, 1-nghost:ny+nghost, 1-nghost:nz+nghost, nv]
  real(dp), allocatable :: Q0 (:,:,:,:)              ! RK3 snapshot
  real(dp), allocatable :: RHS(:,:,:,:)              ! residual
  real(dp), allocatable :: F  (:,:,:,:)              ! x-direction face flux [0:nx, 0:ny, nv]
  real(dp), allocatable :: QL (:,:,:,:), QR(:,:,:,:) ! reconstructed face states
  real(dp), allocatable :: Qw(:,:,:,:)

  real(dp), allocatable :: Q_vis (:,:,:,:)              ! [1-nghost:nx+nghost, 1-nghost:ny+nghost, 1-nghost:nz+nghost, nv]

  ! time
  real(dp) :: t, dt
  real(dp) :: dx_min,dy_min,dz_min
  integer  :: step
  integer  :: i,j,k,l

  ! realtime
  real(dp) :: t1,t2,ttotal


  !--- MPI ---
  integer js,je,ks,ke,ierror,ierr,ierf
  !-----------

  !--- MPI ---
  Call mpi_init(ierror)
  Call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  Call mpi_comm_rank(mpi_comm_world,my_rank,ierr)
  Call mp_setup_division(nx, ny, nz)

  js=j_sta
  je=j_end
  ks=k_sta
  ke=k_end
  !-----------

  call build_grid_uniform(x_min, x_max, y_min, y_max, z_min, z_max)
  if(my_rank == root) write(*,*) "Complete build grid"
  call allocate_fields()
  if(my_rank == root) write(*,*) "Complete allocate"
  call initialize_Sym()
  if(my_rank == root) write(*,*) "Complete initialize"
  t = 0.0_dp; step = 0; ttotal = 0.0_dp

  !$OMP parallel default(none)         &
  !$OMP & shared(Q,Q0,QL,QR,Qw,RHS,F,  &
  !$OMP &        ks,ke,js,je,          &
  !$OMP &        t,dt,my_rank,step,    &
  !$OMP &        t2,t1,ttotal          )

  do while (t < t_max)
    !$OMP master
    if(my_rank == root) write(*,*) step,t,dt,t2-t1,ttotal
    t1 = MPI_Wtime()
    !$OMP end master
    !$OMP barrier

    !$OMP single
    call compute_dt(Q, dt)
    !$OMP end single

    !$OMP master
    if (t + dt > t_max) dt = t_max - t
    !$OMP end master
    !$OMP barrier

    call step_rk3(Q, dt)

    !$OMP master
    t = t + dt; step = step + 1
    if (mod(step, output_frequency) == 0) then
      !call write_vtk_data(step, my_rank, Q(1:nx,js:je,ks:ke,:), t) ! 3D VTKデータを出力
      call write_bin_data(step, my_rank, Q(1:nx,js:je,ks:ke,:), t) ! 3D VTKデータを出力
    end if
    t2 = MPI_Wtime()
    ttotal = ttotal+(t2-t1)
    call mp_barrier
    !$OMP end master
    !$OMP barrier


  end do

  !$OMP end parallel

contains

  subroutine build_grid_uniform(a, b, c, d, e, f)
    real(dp), intent(in) :: a, b, c, d, e, f 
    real(dp) :: dx,dy,dz
    integer :: i,j,k
    allocate(x_edge(-1:nx,js-2:je,ks-2:ke), x_cell(-2:nx,js-2:je,ks-2:ke))
    allocate(y_edge(-1:nx,js-2:je,ks-2:ke), y_cell(-2:nx,js-2:je,ks-2:ke))
    allocate(z_edge(-1:nx,js-2:je,ks-2:ke), z_cell(-2:nx,js-2:je,ks-2:ke))
    allocate(vol   (0:nx,js-1:je,ks-1:ke) )
    allocate(area_x(0:nx,js-1:je,ks-1:ke), area_y(0:nx,js-1:je,ks-1:ke), area_z(0:nx,js-1:je,ks-1:ke) )

    dx_min = x_max; dy_min = y_max; dz_min = z_max

    do k = ks-2, ke
    do j = js-2, je
    do i = -1, nx
      x_edge(i,j,k) = a + (b-a) * real(i,dp) / real(nx,dp)
      y_edge(i,j,k) = c + (d-c) * real(j,dp) / real(ny,dp)
      z_edge(i,j,k) = e + (f-e) * real(k,dp) / real(nz,dp)
    end do; end do; end do

    do k = ks-1, ke
    do j = js-1, je
    do i = 0, nx
      x_cell(i,j,k) = 0.5_dp*(x_edge(i-1,j  ,k  )+x_edge(i,j,k))
      y_cell(i,j,k) = 0.5_dp*(y_edge(i  ,j-1,k  )+y_edge(i,j,k))
      z_cell(i,j,k) = 0.5_dp*(z_edge(i  ,j  ,k-1)+z_edge(i,j,k))

      dx = x_edge(i  ,j  ,k  )-x_edge(i-1,j  ,k  )
      dy = y_edge(i  ,j  ,k  )-y_edge(i  ,j-1,k  )
      dz = z_edge(i  ,j  ,k  )-z_edge(i  ,j  ,k-1)

      dx_min = min(dx_min,dx)
      dy_min = min(dy_min,dy)
      dz_min = min(dz_min,dz)

      area_x(i,j,k) = dy*dz       !Sx = deltay*deltaz
      area_y(i,j,k) = dz*dx       !Sy = deltaz*deltax
      area_z(i,j,k) = dx*dy       !Sz = deltax*deltay

      vol(i,j,k)    = dx*dy*dz
    end do; end do; end do

  end subroutine build_grid_uniform

  subroutine allocate_fields()
    allocate(Q  (1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv))
    allocate(Q0 (1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv))
    allocate(RHS(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv))
    allocate(F (0:nx, js-1:je, ks-1:ke, nv))
    allocate(QL(0:nx, js-1:je, ks-1:ke, nv), QR(0:nx, js-1:je, ks-1:ke, nv))
    allocate(Q_vis(1:nx,1:ny,1:nz,nv))
    allocate(Qw(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv))
  end subroutine allocate_fields

  subroutine initialize_Sym()
    integer :: i,j,k
    real(dp) :: rho, u, v, w, p
    real(dp) :: C_1,C_2
    real(dp) :: x2,y2,z2

    real(dp),parameter :: rho_0 = 1.0_dp
    real(dp),parameter :: Mach  = 0.5_dp

      C_1=1.0_dp/gamma
      C_2=(rho_0*Mach*Mach)/16.0_dp

      do k = ks, ke
      do j = js, je
      do i = 1, nx
        !if (y_cell(i,j,k) < y_center) then
        !  rho = 1.0_dp; u = 0.0_dp; v = 0.0_dp; w = 0.0_dp; p = 1.0_dp
        !else
        !  rho = 0.1_dp; u = 0.0_dp; v = 0.0_dp; w = 0.0_dp; p = 0.1_dp
        !end if

        x2 = 2.0_dp*x_cell(i,j,k); y2 = 2.0_dp*y_cell(i,j,k); z2 = 2.0_dp*z_cell(i,j,k)

        rho=rho_0
        u= Mach*dsin(x_cell(i,j,k))*dcos(y_cell(i,j,k))*dcos(z_cell(i,j,k))
        v=-Mach*dcos(x_cell(i,j,k))*dsin(y_cell(i,j,k))*dcos(z_cell(i,j,k))
        w=0.0d0
        p=C_1+C_2*( dcos(x2)+dcos(y2) )*(dcos(z2)+2.0_Dp)

        Q(i,j,k,1) = rho
        Q(i,j,k,2) = rho*u
        Q(i,j,k,3) = rho*v
        Q(i,j,k,4) = rho*w
        Q(i,j,k,5) = p/(gamma-1.0_dp) + 0.5_dp*rho*(u*u+v*v+w*w)
      end do; end do; end do
    call apply_bc(Q)

    call write_bin_data(0, my_rank, Q(1:nx,js:je,ks:ke,:), 0.0_dp) ! 3D VTKデータを出力

  end subroutine initialize_Sym


  subroutine apply_bc(A)
    real(dp), intent(inout) :: A(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer :: i,j,k,g

    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do g = 1, nghost
      A( 1-g,j,k,:)  = A((nx+1)-g,j,k,:)
      A(nx+g,j,k,:)  = A(     0+g,j,k,:)
    end do; end do; end do
    !$OMP END DO

    !$OMP master
    call BC_Periodic_y_dir_MPI_R8(A)
    !$OMP end master
    !$OMP barrier

    !$OMP master
    call BC_Periodic_z_dir_MPI_R8(A)
    !$OMP end master
    !$OMP barrier

    !$OMP master
    call mp_send_recv_pre_r8_Vec(A,nghost,1-nghost,nx+nghost,js-nghost,je+nghost,ks-nghost,ke+nghost)
    !$OMP end master
    !$OMP barrier


  end subroutine apply_bc


  Subroutine BC_Periodic_y_dir_MPI_R8( F1 )

    Use module_mpi
    Implicit none
    real(dp),intent(inout) :: F1(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    !real(dp), dimension(:,:,:,:), allocatable :: dumcomzx_r8,dumcomzx_s8
    real(dp), allocatable :: sendbuf(:,:,:,:), recv_hi(:,:,:,:), recv_lo(:,:,:,:)
    Integer i, k, g ! 変数
    Integer icom

    !--- MPI ---
    integer ierror,ierr
    integer dum_len, partner, tag
    !integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
    integer isend(2),irecv(2)
    integer :: istatus(MPI_STATUS_SIZE)
    !-----------

    If(js==1.or.je==ny) then
      !Allocate(dumcomzx_r8(1:nx,nghost,ks:ke,nv))
      !Allocate(dumcomzx_s8(1:nx,nghost,ks:ke,nv))
      allocate(sendbuf(1:nx,nghost,ks:ke,nv))
      allocate(recv_hi(1:nx,nghost,ks:ke,nv))
      allocate(recv_lo(1:nx,nghost,ks:ke,nv))
      dum_len=nx*nghost*(ke-ks+1)*nv
      tag=1
      Do icom=0,Ndiv_Nz-1
        If(ks==kksta(icom)) then
          If(je==ny) then
           partner = itable(0, icom)
            Do k = ks,ke
            Do i = 1, nx, 1
              Do g = 1, nghost, 1
              sendbuf(i,g,k,:) = F1(i,ny+g-nghost,k,:)
              End Do
            End Do; End Do
            call MPI_Sendrecv( sendbuf(1,1,ks,1), dum_len, MPI_DOUBLE_PRECISION, partner, tag, &
                               recv_hi(1,1,ks,1), dum_len, MPI_DOUBLE_PRECISION, partner, tag, &
                               MPI_COMM_WORLD, istatus, ierr )
            !Call mpi_isend(dumcomzx_s8(1,1,ks,1),dum_len,MPI_DOUBLE_PRECISION,                        &
            !                                  itable(0,icom),1,MPI_COMM_WORLD,isend(1),ierr)
            !Call mpi_irecv(dumcomzx_r8(1,1,ks,1),dum_len,MPI_DOUBLE_PRECISION,                        &
            !                                  itable(0,icom),1,MPI_COMM_WORLD,isend(2),ierr)
          End If

          If(js==1) then
            partner = itable(Ndiv_Ny-1, icom)
            Do k = ks,ke
            Do i = 1, nx, 1
              Do g = 1, nghost, 1
              sendbuf(i,g,k,:) = F1(i,g,k,:)
              End Do
            End Do; End Do
            call MPI_Sendrecv( sendbuf(1,1,ks,1), dum_len, MPI_DOUBLE_PRECISION, partner, tag, &
                               recv_lo(1,1,ks,1), dum_len, MPI_DOUBLE_PRECISION, partner, tag, &
                               MPI_COMM_WORLD, istatus, ierr )
            !dum_len=nx*nghost*(ke-ks+1)*nv
            !Call mpi_isend(dumcomzx_s8(1,1,ks,1),dum_len,MPI_DOUBLE_PRECISION,                        &
            !                                  itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(2),ierr)
            !Call mpi_irecv(dumcomzx_r8(1,1,ks,1),dum_len,MPI_DOUBLE_PRECISION,                        &
            !                                  itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(1),ierr)
          End if

          !Call mpi_wait(isend(1),istatus,ierr)
          !Call mpi_wait(isend(2),istatus,ierr)

        End if
      End Do

      If(je==ny) then
        Do k = ks,ke
        Do i = 1, nx, 1
          Do g = 1, nghost, 1
          F1(i,ny+g,k,:)=recv_hi(i,g,k,:)
          End Do
        End Do; End Do
      End If

      If(js==1) then
        Do k = ks,ke
        Do i = 1, nx, 1
          Do g = 1, nghost, 1
          F1(i,1-g,k,:)=recv_lo(i,nghost-g+1,k,:)
          End Do
        End Do; End Do
      End if

      Deallocate(sendbuf, recv_hi, recv_lo)
    End If

    Return
  End Subroutine






  Subroutine BC_Periodic_z_dir_MPI_R8( F1 )

    Use module_mpi
    Implicit none
    real(dp),intent(inout) :: F1(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    real(dp), dimension(:,:,:,:), allocatable :: dumcomxy_s8,dumcomxy_r8

    Integer i, j, k, g ! 変数
    Integer icom

    !--- MPI ---
    integer ierror,ierr
    integer dum_len
    integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
    !-----------

    If(ks==1.or.ke==nz) then
      Allocate(dumcomxy_s8(1:nx,js:je,nghost,1:nv))
      Allocate(dumcomxy_r8(1:nx,js:je,nghost,1:nv))
      Do icom=0,Ndiv_Ny-1
        If(j_sta==jjsta(icom)) then

          If(ke==nz) then
            Do j = js,je
            Do i = 1, nx, 1
              Do g = 1, nghost, 1
              dumcomxy_s8(i,j,g,:) = F1(i,j,nz+g-nghost,:)
              End Do
            End Do; End Do
            dum_len=nx*nghost*(ke-ks+1)*nv
            Call mpi_isend(dumcomxy_s8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,0),1,MPI_COMM_WORLD,isend(1),ierr)
            Call mpi_irecv(dumcomxy_r8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,0),1,MPI_COMM_WORLD,isend(2),ierr)
          End If

          If(ks==1) then
            Do j = js,je
            Do i = 1, nx, 1
              Do g = 1, nghost, 1
              dumcomxy_s8(i,j,g,:) = F1(i,j,g,:)
              End Do
            End Do; End Do
            dum_len=nx*nghost*(ke-ks+1)*nv
            Call mpi_isend(dumcomxy_s8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(2),ierr)
            Call mpi_irecv(dumcomxy_r8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(1),ierr)
          End if

          Call mpi_wait(isend(1),istatus,ierr)
          Call mpi_wait(isend(2),istatus,ierr)

        End if
      End Do


      If(ke==nz) then
        Do j = js,je
        Do i = 1, nx, 1
          Do g = 1, nghost, 1
          F1(i,j,nz+g,:) = dumcomxy_r8(i,j,g,:)
          End Do
        End Do; End Do
      End If

      If(ks==1) then
        Do j = js,je
        Do i = 1, nx, 1
          Do g = 1, nghost, 1
          F1(i,j,1-g,:) = dumcomxy_r8(i,j,nghost-g+1,:)
          End Do
        End Do; End Do
      End if

      Deallocate(dumcomxy_r8)
      Deallocate(dumcomxy_s8)
    End If

    Return
  End Subroutine


  subroutine compute_dt(Qin, dt)
    real(dp), intent(in)  :: Qin(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    real(dp), intent(out) :: dt
    integer :: i,j,k
    real(dp) :: rho,u,v,w,p,a,maxs

    maxs = 0.0_dp

    do k = ks, ke
    do j = js, je
    do i = 1, nx
      rho = max(Qin(i,j,k,1), small_rho)
      u   = Qin(i,j,k,2)/rho
      v   = Qin(i,j,k,3)/rho
      w   = Qin(i,j,k,4)/rho
      p   = max( (gamma-1.0_dp)*(Qin(i,j,k,5) - 0.5_dp*rho*(u*u+v*v+w*w)), small_p )
      a   = sqrt(gamma*p/rho)
      maxs = max(maxs, abs(u)+a, abs(v)+a, abs(w)+a)
    end do;end do;end do

    dt = cfl * min( dx_min/maxs, dy_min/maxs , dz_min/maxs )

    call mp_barrier
    call mp_allminr8(dt)

  end subroutine compute_dt

  subroutine step_rk3(Qinout, dt)
    real(dp), intent(inout) :: Qinout(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    real(dp), intent(in)    :: dt
    integer :: i,j,k

    !$OMP DO collapse(2) schedule(static)
    do k = ks-nghost,ke+nghost
    do j = js-nghost,je+nghost
    do i = 1-nghost,nx+nghost
      Q0(i,j,k,:) = Qinout(i,j,k,:)
    end do; end do; end do
    !$OMP end do

    ! --- 1st step --- !
    call compute_rhs(Qinout, RHS)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-nghost,ke+nghost
    do j = js-nghost,je+nghost
    do i = 1-nghost,nx+nghost
      Qinout(i,j,k,:) = Q0(i,j,k,:)+dt*RHS(i,j,k,:)
    end do; end do; end do
    !$OMP end do

    ! --- 2nd step --- !
    call compute_rhs(Qinout, RHS)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-nghost,ke+nghost
    do j = js-nghost,je+nghost
    do i = 1-nghost,nx+nghost
      Qinout(i,j,k,:) = 0.75_dp*Q0(i,j,k,:) + 0.25_dp*(Qinout(i,j,k,:) + dt*RHS(i,j,k,:))
    end do; end do; end do
    !$OMP end do

    ! --- 3rd step --- !
    call compute_rhs(Qinout, RHS)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-nghost,ke+nghost
    do j = js-nghost,je+nghost
    do i = 1-nghost,nx+nghost
      Qinout(i,j,k,:) = (1.0_dp/3.0_dp)*Q0(i,j,k,:) + (2.0_dp/3.0_dp)*(Qinout(i,j,k,:) + dt*RHS(i,j,k,:))
    end do; end do; end do
    !$OMP end do

  end subroutine step_rk3

  subroutine compute_rhs(Qin, R)
    real(dp), intent(inout)  :: Qin(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    real(dp), intent(out) :: R  (1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer :: i,j,k,g
    integer :: dir
    real(dp) :: dR(nv)

    !$OMP DO collapse(2) schedule(static)
    do k = ks-nghost,ke+nghost
    do j = js-nghost,je+nghost
    do i = 1-nghost,nx+nghost
      R (i,j,k,:) = 0.0_dp
    end do; end do; end do
    !$OMP end do
    call apply_bc(Qin)

    dir = 1
    call flux_KEEP(Qin, F, dir)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do i = 1, nx
      dR(:) = ( area_x(i,j,k)*F(i,j,k,:) - area_x(i-1,j,k)*F(i-1,j,k,:) ) / vol(i,j,k)
      R(i,j,k,1) = R(i,j,k,1)-dR(1)
      R(i,j,k,2) = R(i,j,k,2)-dR(2)
      R(i,j,k,3) = R(i,j,k,3)-dR(3)
      R(i,j,k,4) = R(i,j,k,4)-dR(4)
      R(i,j,k,5) = R(i,j,k,5)-dR(5)
    end do; end do; end do
    !$OMP end do

    dir = 2
    call flux_KEEP(Qin, F, dir)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do i = 1, nx
      dR(:) = ( area_y(i,j,k)*F(i,j,k,:) - area_y(i,j-1,k)*F(i,j-1,k,:) ) / vol(i,j,k)
      R(i,j,k,1) = R(i,j,k,1)-dR(1)
      R(i,j,k,2) = R(i,j,k,2)-dR(2)
      R(i,j,k,3) = R(i,j,k,3)-dR(3)
      R(i,j,k,4) = R(i,j,k,4)-dR(4)
      R(i,j,k,5) = R(i,j,k,5)-dR(5)
    end do; end do; end do
    !$OMP end do

    dir = 3
    call flux_KEEP(Qin, F, dir)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do i = 1, nx
      dR(:) = ( area_z(i,j,k)*F(i,j,k,:) - area_z(i,j,k-1)*F(i,j,k-1,:) ) / vol(i,j,k)
      R(i,j,k,1) = R(i,j,k,1)-dR(1)
      R(i,j,k,2) = R(i,j,k,2)-dR(2)
      R(i,j,k,3) = R(i,j,k,3)-dR(3)
      R(i,j,k,4) = R(i,j,k,4)-dR(4)
      R(i,j,k,5) = R(i,j,k,5)-dR(5)
    end do; end do; end do
    !$OMP end do

  end subroutine compute_rhs

  subroutine rotation(Qin, Qout, direction)
    real(dp), intent(in)  :: Qin(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer , intent(in)  :: direction
    real(dp), intent(out) :: Qout(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer :: i, j, k
    if (direction == 1) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-nghost, ke+nghost
      do j = js-nghost, je+nghost
      do i = 1-nghost, nx+nghost
        Qout(i,j,k,1) =  Qin(i,j,k,1)
        Qout(i,j,k,2) =  Qin(i,j,k,2)
        Qout(i,j,k,3) =  Qin(i,j,k,3)
        Qout(i,j,k,4) =  Qin(i,j,k,4)
        Qout(i,j,k,5) =  Qin(i,j,k,5)
      end do; end do; end do
      !$OMP end do

    else if (direction == 2) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-nghost, ke+nghost
      do j = js-nghost, je+nghost
      do i = 1-nghost, nx+nghost
        Qout(i,j,k,1) =  Qin(i,j,k,1)
        Qout(i,j,k,2) =  Qin(i,j,k,3)
        Qout(i,j,k,3) = -Qin(i,j,k,2)
        Qout(i,j,k,4) =  Qin(i,j,k,4)
        Qout(i,j,k,5) =  Qin(i,j,k,5)
      end do; end do; end do
      !$OMP end do

    else if (direction == 3) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-nghost, ke+nghost
      do j = js-nghost, je+nghost
      do i = 1-nghost, nx+nghost
        Qout(i,j,k,1) =  Qin(i,j,k,1)
        Qout(i,j,k,2) =  Qin(i,j,k,4)
        Qout(i,j,k,3) =  Qin(i,j,k,3)
        Qout(i,j,k,4) = -Qin(i,j,k,2)
        Qout(i,j,k,5) =  Qin(i,j,k,5)
      end do; end do; end do
      !$OMP end do
    end if

  end subroutine rotation

  ! ---------------- Roe flux (same as before, with HH2 entropy fix) ---------
  subroutine flux_KEEP(Qin, Fface, direction)
    real(dp), intent(in)  :: Qin(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer , intent(in)  :: direction
    real(dp), intent(out) :: Fface(0:nx,js-1:je,ks-1:ke,nv)
    integer :: i,j,k
    real(dp) :: rp1,up1,vp1,wp1,pp1,Hp1
    real(dp) :: rm1,um1,vm1,wm1,pm1,Hm1
    real(dp) :: Qp1,Qp2,Qp3,Qp4,Qp5
    real(dp) :: Qm1,Qm2,Qm3,Qm4,Qm5
    real(dp) :: rmp,ump,vmp,wmp,Lmp
    real(dp) :: Ck,Mxk,Myk,Mzk,Kk,Lk,Gk,Pk

    if(direction == 1) then
       !$OMP DO collapse(2) schedule(static)
       do k = ks-1, ke
       do j = js-1, je
       do i = 0, nx
         ! Build Roe eigenvectors from a central pair around face j
         !Qm1 = Qin(i  , j, k, :)
         !Qp1 = Qin(i+1, j, k, :)

         Qm1 = Qin(i  ,j,k,1) 
         Qm2 = Qin(i  ,j,k,2) 
         Qm3 = Qin(i  ,j,k,3) 
         Qm4 = Qin(i  ,j,k,4) 
         Qm5 = Qin(i  ,j,k,5) 

         Qp1 = Qin(i+1,j,k,1) 
         Qp2 = Qin(i+1,j,k,2) 
         Qp3 = Qin(i+1,j,k,3) 
         Qp4 = Qin(i+1,j,k,4) 
         Qp5 = Qin(i+1,j,k,5) 

         rm1 = max(Qm1, small_rho)
         um1 = Qm2/rm1
         vm1 = Qm3/rm1
         wm1 = Qm4/rm1
         pm1 = max((gamma-1.0_dp)*(Qm5 - 0.5_dp*rm1*(um1*um1+vm1*vm1+wm1*wm1)), small_p)
         Hm1 = (Qm5 + pm1) / rm1

         rp1 = max(Qp1, small_rho)
         up1 = Qp2/rp1
         vp1 = Qp3/rp1
         wp1 = Qp4/rp1
         pp1 = max((gamma-1.0_dp)*(Qp5 - 0.5_dp*rp1*(up1*up1+vp1*vp1+wp1*wp1)), small_p)
         Hp1 = (Qp5 + pp1) / rp1

         rmp=0.5_dp*(rm1+rp1)
         ump=0.5_dp*(um1+up1)
         vmp=0.5_dp*(vm1+vp1)
         wmp=0.5_dp*(wm1+wp1)
         Lmp=0.5_dp*((pm1/rm1)+(pp1/rp1))/(gamma-1.0_dp)

         Ck  = rmp*ump
         Mxk = Ck*ump
         Myk = Ck*vmp
         Mzk = Ck*wmp
         Kk  = Ck*0.5_dp*(um1*up1+vm1*vp1+wm1*wp1)
         Lk  = Ck*Lmp
         Gk  = 0.5_dp*(pm1+pp1)
         Pk  = 0.5_dp*(up1*pm1+um1*pp1)

         Fface(i,j,k,1) = Ck
         Fface(i,j,k,2) = Mxk+Gk
         Fface(i,j,k,3) = Myk
         Fface(i,j,k,4) = Mzk
         Fface(i,j,k,5) = Kk+Lk+Pk

       end do; end do; end do
       !$OMP end do

    else if(direction == 2) then
       !$OMP DO collapse(2) schedule(static)
       do k = ks-1, ke
       do j = js-1, je
       do i = 0, nx
         ! Build Roe eigenvectors from a central pair around face j
         Qm1 = Qin(i,j  ,k,1) 
         Qm2 = Qin(i,j  ,k,2) 
         Qm3 = Qin(i,j  ,k,3) 
         Qm4 = Qin(i,j  ,k,4) 
         Qm5 = Qin(i,j  ,k,5) 

         Qp1 = Qin(i,j+1,k,1) 
         Qp2 = Qin(i,j+1,k,2) 
         Qp3 = Qin(i,j+1,k,3) 
         Qp4 = Qin(i,j+1,k,4) 
         Qp5 = Qin(i,j+1,k,5) 

         rm1 = max(Qm1, small_rho)
         um1 = Qm2/rm1
         vm1 = Qm3/rm1
         wm1 = Qm4/rm1
         pm1 = max((gamma-1.0_dp)*(Qm5 - 0.5_dp*rm1*(um1*um1+vm1*vm1+wm1*wm1)), small_p)
         Hm1 = (Qm5 + pm1) / rm1

         rp1 = max(Qp1, small_rho)
         up1 = Qp2/rp1
         vp1 = Qp3/rp1
         wp1 = Qp4/rp1
         pp1 = max((gamma-1.0_dp)*(Qp5 - 0.5_dp*rp1*(up1*up1+vp1*vp1+wp1*wp1)), small_p)
         Hp1 = (Qp5 + pp1) / rp1

         rmp=0.5_dp*(rm1+rp1)
         ump=0.5_dp*(um1+up1)
         vmp=0.5_dp*(vm1+vp1)
         wmp=0.5_dp*(wm1+wp1)
         Lmp=0.5_dp*((pm1/rm1)+(pp1/rp1))/(gamma-1.0_dp)

         Ck  = rmp*vmp
         Mxk = Ck*ump
         Myk = Ck*vmp
         Mzk = Ck*wmp
         Kk  = Ck*0.5_dp*(um1*up1+vm1*vp1+wm1*wp1)
         Lk  = Ck*Lmp
         Gk  = 0.5_dp*(pm1+pp1)
         Pk  = 0.5_dp*(vp1*pm1+vm1*pp1)

         Fface(i,j,k,1) = Ck
         Fface(i,j,k,2) = Mxk
         Fface(i,j,k,3) = Myk+Gk
         Fface(i,j,k,4) = Mzk
         Fface(i,j,k,5) = Kk+Lk+Pk

       end do; end do; end do
       !$OMP end do

    else if(direction == 3) then
       !$OMP DO collapse(2) schedule(static)
       do k = ks-1, ke
       do j = js-1, je
       do i = 0, nx
         ! Build Roe eigenvectors from a central pair around face j
         Qm1 = Qin(i,j,k  ,1) 
         Qm2 = Qin(i,j,k  ,2) 
         Qm3 = Qin(i,j,k  ,3) 
         Qm4 = Qin(i,j,k  ,4) 
         Qm5 = Qin(i,j,k  ,5) 

         Qp1 = Qin(i,j,k+1,1) 
         Qp2 = Qin(i,j,k+1,2) 
         Qp3 = Qin(i,j,k+1,3) 
         Qp4 = Qin(i,j,k+1,4) 
         Qp5 = Qin(i,j,k+1,5) 

         rm1 = max(Qm1, small_rho)
         um1 = Qm2/rm1
         vm1 = Qm3/rm1
         wm1 = Qm4/rm1
         pm1 = max((gamma-1.0_dp)*(Qm5 - 0.5_dp*rm1*(um1*um1+vm1*vm1+wm1*wm1)), small_p)
         Hm1 = (Qm5 + pm1) / rm1

         rp1 = max(Qp1, small_rho)
         up1 = Qp2/rp1
         vp1 = Qp3/rp1
         wp1 = Qp4/rp1
         pp1 = max((gamma-1.0_dp)*(Qp5 - 0.5_dp*rp1*(up1*up1+vp1*vp1+wp1*wp1)), small_p)
         Hp1 = (Qp5 + pp1) / rp1

         rmp=0.5_dp*(rm1+rp1)
         ump=0.5_dp*(um1+up1)
         vmp=0.5_dp*(vm1+vp1)
         wmp=0.5_dp*(wm1+wp1)
         Lmp=0.5_dp*((pm1/rm1)+(pp1/rp1))/(gamma-1.0_dp)

         Ck  = rmp*wmp
         Mxk = Ck*ump
         Myk = Ck*vmp
         Mzk = Ck*wmp
         Kk  = Ck*0.5_dp*(um1*up1+vm1*vp1+wm1*wp1)
         Lk  = Ck*Lmp
         Gk  = 0.5_dp*(pm1+pp1)
         Pk  = 0.5_dp*(wp1*pm1+wm1*pp1)

         Fface(i,j,k,1) = Ck
         Fface(i,j,k,2) = Mxk
         Fface(i,j,k,3) = Myk
         Fface(i,j,k,4) = Mzk+Gk
         Fface(i,j,k,5) = Kk+Lk+Pk

       end do; end do; end do
       !$OMP end do
    end if

  end subroutine flux_KEEP









































  subroutine write_vtk_data(step, my_rank, Qin, t)
    integer, intent(in) :: step, my_rank
    real(dp), intent(in) :: Qin(1:nx, js:je, ks:ke, nv)
    real(dp), intent(in) :: t
    integer :: i, j, k
    real(dp) :: rho, u, v, w, p
    real(dp) :: xm,xp,ym,yp,zm,zp
    real(dp) :: xc1,xc2,yc1,yc2,zc1,zc2
    real(dp) :: xc,yc,zc
    character(len=256) :: fname
    character(len=*), parameter :: vtk_fmt = '(ES22.12E3)'

    ! VTKファイルはoutput_frequencyごとにのみ書き出す
    write(fname, '(A,I0.5,A,I0.5,A)') '3d_result_step_', step, '_rank', my_rank, '.vtk'
    open(unit=30, file=fname, status='replace')

    ! VTK Header
    write(30, '(A)') '# vtk DataFile Version 2.0'
    write(30, '(A, F12.6)') 'Time = ', t
    write(30, '(A)') 'ASCII'
    write(30, '(A)') 'DATASET STRUCTURED_GRID'
    write(30, '(A,I6,I6,I6)') 'DIMENSIONS ', nx, (je-js)+1, (ke-ks)+1
    write(30, '(A,I12,A)') 'POINTS ', nx*((je-js)+1)*((ke-ks)+1), ' double'

    xm = x_min
    xp = x_max
    ym = y_min
    yp = y_max
    zm = z_min
    zp = z_max

    ! Write grid points
    do k = ks, ke
      do j = js, je
        do i = 1, nx

          !xc1 = xm + (xp-xm) * real(i-1,dp) / real(nx,dp)
          !xc2 = xm + (xp-xm) * real(i  ,dp) / real(nx,dp)
          !yc1 = ym + (yp-ym) * real(j-1,dp) / real(ny,dp)
          !yc2 = ym + (yp-ym) * real(j  ,dp) / real(ny,dp)
          !zc1 = zm + (zp-zm) * real(k-1,dp) / real(nz,dp)
          !zc2 = zm + (zp-zm) * real(k  ,dp) / real(nz,dp)

          !xc = 0.5_dp*(xc1+xc2)
          !yc = 0.5_dp*(yc1+yc2)
          !zc = 0.5_dp*(zc1+zc2)

          !write(30, vtk_fmt) xc, yc, zc
          write(30, vtk_fmt) x_cell(i,j,k), y_cell(i,j,k), z_cell(i,j,k)
        end do
      end do
    end do

    ! Write data
    write(30, '(A,I12)') 'POINT_DATA ', nx*((je-js)+1)*((ke-ks)+1)

    ! --- Density (rho) ---
    write(30, '(A)') 'SCALARS rho double 1'
    write(30, '(A)') 'LOOKUP_TABLE default'
    do k = ks, ke
      do j = js, je
        do i = 1, nx
          rho = Qin(i,j,k,1)
          write(30, vtk_fmt) rho
        end do
      end do
    end do

    ! --- X-Velocity (u) ---
    write(30, '(A)') 'SCALARS u double 1'
    write(30, '(A)') 'LOOKUP_TABLE default'
    do k = ks, ke
      do j = js, je
        do i = 1, nx
          rho = max(Qin(i,j,k,1), small_rho)
          u   = Qin(i,j,k,2) / rho
          write(30, vtk_fmt) u
        end do
      end do
    end do

    ! --- Y-Velocity (v) ---
    write(30, '(A)') 'SCALARS v double 1'
    write(30, '(A)') 'LOOKUP_TABLE default'
    do k = ks, ke
      do j = js, je
        do i = 1, nx
          rho = max(Qin(i,j,k,1), small_rho)
          v   = Qin(i,j,k,3) / rho
          write(30, vtk_fmt) v
        end do
      end do
    end do

    ! --- Z-Velocity (w) ---
    write(30, '(A)') 'SCALARS w double 1'
    write(30, '(A)') 'LOOKUP_TABLE default'
    do k = ks, ke
      do j = js, je
        do i = 1, nx
          rho = max(Qin(i,j,k,1), small_rho)
          w   = Qin(i,j,k,4) / rho
          write(30, vtk_fmt) w
        end do
      end do
    end do

    ! --- Pressure (p) ---
    write(30, '(A)') 'SCALARS p double 1'
    write(30, '(A)') 'LOOKUP_TABLE default'
    do k = ks, ke
      do j = js, je
        do i = 1, nx
          rho = max(Qin(i,j,k,1), small_rho)
          u   = Qin(i,j,k,2) / rho
          v   = Qin(i,j,k,3) / rho
          w   = Qin(i,j,k,4) / rho
          p   = (gamma-1.0_dp) * (Qin(i,j,k,5) - 0.5_dp*rho*(u*u + v*v + w*w))
          p   = max(p, small_p)
          write(30, vtk_fmt) p
        end do
      end do
    end do

    close(30)
    write(*,'(A,A)') 'Wrote VTK: ', trim(fname)
  end subroutine write_vtk_data


  subroutine write_bin_data(step, my_rank, Qin, t)
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none
    integer, intent(in) :: step, my_rank
    real(dp), intent(in) :: Qin(1:nx, js:je, ks:ke, nv)
    real(dp), intent(in) :: t

    integer :: u, ios
    character(len=256) :: fname
    character(len=8)   :: magic
    integer(int32) :: ndim, dtype_code
    integer(int32) :: shp(4)
    integer(int32) :: meta(6)
    real(dp) :: t_write

    ! ---- file name (rankごと) ----
    write(fname, '(A,I0.5,A,I0.5,A)') 'output/3d_result_step_', step, '_rank', my_rank, '.fbn'

    ! ---- header ----
    magic = 'FBN1' // char(0) // char(0) // char(0) // char(0)
    ndim  = 4_int32
    ! shape: (nx, ny_local, nz_local, nv)
    shp   = [ int(nx, int32), int((je-js)+1, int32), int((ke-ks)+1, int32), int(nv, int32) ]
    dtype_code = 2_int32   ! 1=float32, 2=float64(dp)

    ! 追加メタ情報（任意だが解析で便利）:
    ! meta = [js, je, ks, ke, step, rank]
    meta = [ int(js,int32), int(je,int32), int(ks,int32), int(ke,int32), int(step,int32), int(my_rank,int32) ]
    t_write = t

    open(newunit=u, file=fname, access='stream', form='unformatted', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'ERROR: cannot open binary file: ', trim(fname)
      error stop
    end if

    ! ---- write header ----
    write(u) magic
    write(u) ndim
    write(u) shp
    write(u) dtype_code
    write(u) meta
    write(u) t_write

    ! ---- write data (Fortran配列順のまま) ----
    write(u) Qin

    close(u)
    write(*,'(A,A)') 'Wrote BIN: ', trim(fname)
  end subroutine write_bin_data


end program 