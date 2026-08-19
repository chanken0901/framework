program shock_tube_fv_weno_roe_charZ
  Use module_mpi
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nv = 5

  ! ------------------------------
  ! Simulation parameters
  ! ------------------------------
  integer, parameter :: nx =  256
  integer, parameter :: ny =  256
  integer, parameter :: nz =  256
  integer, parameter :: nghost = 3         ! WENO5 needs 3 ghosts
  !integer, parameter :: mode = 0           ! mode=0:constant temperature,  mode=1:constant density

  real(dp), parameter :: x_min = 0.0_dp
  real(dp), parameter :: x_max = 5.0_dp
  real(dp), parameter :: x_center = (x_max-x_min)/2.0_dp
  real(dp), parameter :: y_min = 0.0_dp
  real(dp), parameter :: y_max = 5.0_dp
  real(dp), parameter :: y_center = (y_max-y_min)/2.0_dp
  real(dp), parameter :: z_min = 0.0_dp
  real(dp), parameter :: z_max = 5.0_dp
  real(dp), parameter :: z_center = (z_max-z_min)/2.0_dp

  real(dp), parameter :: t_max = 20.0_dp
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

    call compute_dt(Q, dt)

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
      do k = ks, ke
      do j = js, je
      do i = 1, nx
        if (y_cell(i,j,k) < y_center) then
          rho = 1.0_dp; u = 0.0_dp; v = 0.0_dp; w = 0.0_dp; p = 1.0_dp
        else
          rho = 0.1_dp; u = 0.0_dp; v = 0.0_dp; w = 0.0_dp; p = 0.1_dp
        end if
        Q(i,j,k,1) = rho
        Q(i,j,k,2) = rho*u
        Q(i,j,k,3) = rho*v
        Q(i,j,k,4) = rho*w
        Q(i,j,k,5) = p/(gamma-1.0_dp) + 0.5_dp*rho*(u*u+v*v+w*w)
      end do; end do; end do
    call apply_bc(Q)

    call write_bin_data(0, my_rank, Q(1:nx,js:je,ks:ke,:), t) ! 3D VTKデータを出力

  end subroutine initialize_Sym


  subroutine apply_bc(A)
    real(dp), intent(inout) :: A(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer :: i,j,k,g

    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do g = 1, nghost
      A( 1-g,j,k,:)  = A( 1,j,k,:)
      A(nx+g,j,k,:)  = A(nx,j,k,:)
    end do; end do; end do
    !$OMP END DO

    if(js==1) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks, ke
      do i = 1-nghost,nx+nghost
      do g = 1, nghost
        A(i,1-g,k,:)  = A(i,1,k,:)
      end do; end do; end do
      !$OMP END DO

      !do g = 1, nghost
      !  A(:,1-g,ks:ke, :)  = A(:, 1,ks:ke, :)
      !end do
    end if
    if(je==ny) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-1, ke
      do i = 1-nghost,nx+nghost
      do g = 1, nghost
        A(i,ny+g,k,:)  = A(i,ny,k,:)
      end do; end do; end do
      !$OMP END DO

      !do g = 1, nghost
      !  A(:,ny+g,ks:ke,:)  = A(:,ny,ks:ke,:)
      !end do
    end if

    if(ks==1) then
      !$OMP DO collapse(2) schedule(static)
      do j = js-nghost,je+nghost
      do i = 1-nghost,nx+nghost
      do g = 1,nghost
        A(i,j,1-g,:)  = A(i,j,1,:)
      end do; end do; end do
      !$OMP END DO

      !do g = 1, nghost
      !  A(:,:, 1-g, :)  = A(:,:,1, :)
      !end do
    end if
    if(ke==nz) then
      !$OMP DO collapse(2) schedule(static)
      do j = js-nghost,je+nghost
      do i = 1-nghost,nx+nghost
      do g = 1,nghost
        A(i,j,nz+g,:)  = A(i,j,nz,:)
      end do; end do; end do
      !$OMP END DO

      !do g = 1, nghost
      !  A(:,:,nz+g,:)  = A(:,:,nz,:)
      !end do
    end if

    !$OMP master
    call mp_send_recv_pre_r8_Vec(A,nghost,1-nghost,nx+nghost,js-nghost,je+nghost,ks-nghost,ke+nghost)
    !$OMP end master
    !$OMP barrier


  end subroutine apply_bc

  subroutine compute_dt(Qin, dt)
    real(dp), intent(in)  :: Qin(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    real(dp), intent(out) :: dt
    integer :: i,j,k
    real(dp) :: rho,u,v,w,p,a,maxs
    maxs = 0.0_dp
    !$OMP DO collapse(2) schedule(static)
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
    !$OMP end do

    dt = cfl * min( dx_min/maxs, dy_min/maxs , dz_min/maxs )

    !$OMP master
    call mp_barrier
    call mp_allminr8(dt)
    !$OMP end master
    !$OMP barrier

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
    call rotation(Qin,Qw,dir)
    !call apply_bc(Qw)
    call reconstruct_weno5_charZ(Qw, QL, QR, dir)
    call flux_roe(QL, QR, F)
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
    call rotation(Qin,Qw,dir)
    !call apply_bc(Qw)
    call reconstruct_weno5_charZ(Qw, QL, QR, dir)
    call flux_roe(QL, QR, F)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do i = 1, nx
      dR(:) = ( area_y(i,j,k)*F(i,j,k,:) - area_y(i,j-1,k)*F(i,j-1,k,:) ) / vol(i,j,k)
      R(i,j,k,1) = R(i,j,k,1)-dR(1)
      R(i,j,k,2) = R(i,j,k,2)+dR(3)
      R(i,j,k,3) = R(i,j,k,3)-dR(2)
      R(i,j,k,4) = R(i,j,k,4)-dR(4)
      R(i,j,k,5) = R(i,j,k,5)-dR(5)
    end do; end do; end do
    !$OMP end do

    dir = 3
    call rotation(Qin,Qw,dir)
    !call apply_bc(Qw)
    call reconstruct_weno5_charZ(Qw, QL, QR, dir)
    call flux_roe(QL, QR, F)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do i = 1, nx
      dR(:) = ( area_z(i,j,k)*F(i,j,k,:) - area_z(i,j,k-1)*F(i,j,k-1,:) ) / vol(i,j,k)
      R(i,j,k,1) = R(i,j,k,1)-dR(1)
      R(i,j,k,2) = R(i,j,k,2)+dR(4)
      R(i,j,k,3) = R(i,j,k,3)-dR(3)
      R(i,j,k,4) = R(i,j,k,4)-dR(2)
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

  ! ---------------- Characteristic-wise WENO-Z reconstruction ----------------
  subroutine reconstruct_weno5_charZ(Qin, qL, qR, direction)
    real(dp), intent(in)  :: Qin(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    integer , intent(in)  :: direction
    real(dp), intent(out) :: qL(0:nx,js-1:je,ks-1:ke,nv), qR(0:nx,js-1:je,ks-1:ke,nv)
    integer :: i, j, k, l
    real(dp) :: Rm(5,5), Lm(5,5), lam(5)
    real(dp) :: weL(5), weR(5)
    real(dp) :: v1, v2, v3, v4, v5
    real(dp) :: Hroe, uroe, vroe, wroe, aroe, rroe
    real(dp) :: rL,uL,vL,wL,pL,HL, rR,uR,vR,wR,pR,HR
    real(dp) :: QcL(5), QcR(5)

    if(direction == 1) then
       !$OMP DO collapse(2) schedule(static)
       do k = ks-1, ke
       do j = js-1, je
       do i = 0, nx
         ! Build Roe eigenvectors from a central pair around face j
         QcL = Qin(i  , j, k, :)
         QcR = Qin(i+1, j, k, :)

         call primitives(QcL, rL, uL, vL, wL, pL, HL)
         call primitives(QcR, rR, uR, vR, wR, pR, HR)

         rroe = sqrt(max(rL* rR, small_rho*small_rho))
         uroe = (sqrt(rL)*uL + sqrt(rR)*uR) / (sqrt(rL)+sqrt(rR))
         vroe = (sqrt(rL)*vL + sqrt(rR)*vR) / (sqrt(rL)+sqrt(rR))
         wroe = (sqrt(rL)*wL + sqrt(rR)*wR) / (sqrt(rL)+sqrt(rR))
         Hroe = (sqrt(rL)*HL  + sqrt(rR)*HR)  / (sqrt(rL)+sqrt(rR))
         aroe = sqrt(max((gamma-1.0_dp)*(Hroe - 0.5_dp*(uroe*uroe+vroe*vroe+wroe*wroe)), small_p/rroe))

         call eigen_matrices(uroe, vroe, wroe, Hroe, aroe, Rm, Lm, lam)

         ! Reconstruct each characteristic component with WENO-Z
         do l = 1, 5
           ! left trace at face j from cells j-2..j+2
           v1 = dot_product(Lm(l,:), Qin(i-2,j,k,:))
           v2 = dot_product(Lm(l,:), Qin(i-1,j,k,:))
           v3 = dot_product(Lm(l,:), Qin(i  ,j,k,:))
           v4 = dot_product(Lm(l,:), Qin(i+1,j,k,:))
           v5 = dot_product(Lm(l,:), Qin(i+2,j,k,:))
           weL(l) = weno5z_left(v1,v2,v3,v4,v5)

           ! right trace at face j from cells j-1..j+3
           v1 = dot_product(Lm(l,:), Qin(i-1,j,k,:))
           v2 = dot_product(Lm(l,:), Qin(i  ,j,k,:))
           v3 = dot_product(Lm(l,:), Qin(i+1,j,k,:))
           v4 = dot_product(Lm(l,:), Qin(i+2,j,k,:))
           v5 = dot_product(Lm(l,:), Qin(i+3,j,k,:))
           weR(l) = weno5z_right(v1,v2,v3,v4,v5)
         end do

         qL(i,j,k,:) = matmul(Rm, weL)
         qR(i,j,k,:) = matmul(Rm, weR)
       end do; end do; end do
      !$OMP end do

    else if(direction == 2) then
       !$OMP DO collapse(2) schedule(static)
       do k = ks-1, ke
       do j = js-1, je
       do i = 0, nx
         ! Build Roe eigenvectors from a central pair around face j
         !QcL = Qin(max( 1-nghost, i  ), max( 1-nghost, j), :)
         !QcR = Qin(min(nx+nghost, i+1), min(ny+nghost, j), :)
         QcL = Qin(i, j  , k, :)
         QcR = Qin(i, j+1, k, :)

         call primitives(QcL, rL, uL, vL, wL, pL, HL)
         call primitives(QcR, rR, uR, vR, wR, pR, HR)

         rroe = sqrt(max(rL* rR, small_rho*small_rho))
         uroe = (sqrt(rL)*uL + sqrt(rR)*uR) / (sqrt(rL)+sqrt(rR))
         vroe = (sqrt(rL)*vL + sqrt(rR)*vR) / (sqrt(rL)+sqrt(rR))
         wroe = (sqrt(rL)*wL + sqrt(rR)*wR) / (sqrt(rL)+sqrt(rR))
         Hroe = (sqrt(rL)*HL  + sqrt(rR)*HR)  / (sqrt(rL)+sqrt(rR))
         aroe = sqrt(max((gamma-1.0_dp)*(Hroe - 0.5_dp*(uroe*uroe+vroe*vroe+wroe*wroe)), small_p/rroe))

         call eigen_matrices(uroe, vroe, wroe, Hroe, aroe, Rm, Lm, lam)

         ! Reconstruct each characteristic component with WENO-Z
         do l = 1, 5
           ! left trace at face j from cells j-2..j+2
           v1 = dot_product(Lm(l,:), Qin(i,j-2,k,:))
           v2 = dot_product(Lm(l,:), Qin(i,j-1,k,:))
           v3 = dot_product(Lm(l,:), Qin(i,j  ,k,:))
           v4 = dot_product(Lm(l,:), Qin(i,j+1,k,:))
           v5 = dot_product(Lm(l,:), Qin(i,j+2,k,:))
           weL(l) = weno5z_left(v1,v2,v3,v4,v5)

           ! right trace at face j from cells j-1..j+3
           v1 = dot_product(Lm(l,:), Qin(i,j-1,k,:))
           v2 = dot_product(Lm(l,:), Qin(i,j  ,k,:))
           v3 = dot_product(Lm(l,:), Qin(i,j+1,k,:))
           v4 = dot_product(Lm(l,:), Qin(i,j+2,k,:))
           v5 = dot_product(Lm(l,:), Qin(i,j+3,k,:))
           weR(l) = weno5z_right(v1,v2,v3,v4,v5)
         end do

         qL(i,j,k,:) = matmul(Rm, weL)
         qR(i,j,k,:) = matmul(Rm, weR)
       end do; end do; end do
      !$OMP end do

    else if(direction == 3) then
       !$OMP DO collapse(2) schedule(static)
       do k = ks-1, ke
       do j = js-1, je
       do i = 0, nx
         ! Build Roe eigenvectors from a central pair around face j
         !QcL = Qin(max( 1-nghost, i  ), max( 1-nghost, j), :)
         !QcR = Qin(min(nx+nghost, i+1), min(ny+nghost, j), :)
         QcL = Qin(i, j  , k  ,:)
         QcR = Qin(i, j  , k+1,:)

         call primitives(QcL, rL, uL, vL, wL, pL, HL)
         call primitives(QcR, rR, uR, vR, wR, pR, HR)

         rroe = sqrt(max(rL* rR, small_rho*small_rho))
         uroe = (sqrt(rL)*uL + sqrt(rR)*uR) / (sqrt(rL)+sqrt(rR))
         vroe = (sqrt(rL)*vL + sqrt(rR)*vR) / (sqrt(rL)+sqrt(rR))
         wroe = (sqrt(rL)*wL + sqrt(rR)*wR) / (sqrt(rL)+sqrt(rR))
         Hroe = (sqrt(rL)*HL  + sqrt(rR)*HR)  / (sqrt(rL)+sqrt(rR))
         aroe = sqrt(max((gamma-1.0_dp)*(Hroe - 0.5_dp*(uroe*uroe+vroe*vroe+wroe*wroe)), small_p/rroe))

         call eigen_matrices(uroe, vroe, wroe, Hroe, aroe, Rm, Lm, lam)

         ! Reconstruct each characteristic component with WENO-Z
         do l = 1, 5
           ! left trace at face j from cells j-2..j+2
           v1 = dot_product(Lm(l,:), Qin(i,j,k-2,:))
           v2 = dot_product(Lm(l,:), Qin(i,j,k-1,:))
           v3 = dot_product(Lm(l,:), Qin(i,j,k  ,:))
           v4 = dot_product(Lm(l,:), Qin(i,j,k+1,:))
           v5 = dot_product(Lm(l,:), Qin(i,j,k+2,:))
           weL(l) = weno5z_left(v1,v2,v3,v4,v5)

           ! right trace at face j from cells j-1..j+3
           v1 = dot_product(Lm(l,:), Qin(i,j,k-1,:))
           v2 = dot_product(Lm(l,:), Qin(i,j,k  ,:))
           v3 = dot_product(Lm(l,:), Qin(i,j,k+1,:))
           v4 = dot_product(Lm(l,:), Qin(i,j,k+2,:))
           v5 = dot_product(Lm(l,:), Qin(i,j,k+3,:))
           weR(l) = weno5z_right(v1,v2,v3,v4,v5)
         end do

         qL(i,j,k,:) = matmul(Rm, weL)
         qR(i,j,k,:) = matmul(Rm, weR)
       end do; end do; end do
      !$OMP end do
    end if

  end subroutine reconstruct_weno5_charZ

  subroutine primitives(Qv, rho, u, v, w, p, H)
    real(dp), intent(in)  :: Qv(5)
    real(dp), intent(out) :: rho, u, v, w, p, H
    rho = max(Qv(1), small_rho)
    u   = Qv(2)/rho
    v   = Qv(3)/rho
    w   = Qv(4)/rho
    p   = max((gamma-1.0_dp)*(Qv(5) - 0.5_dp*rho*(u*u+v*v+w*w)), small_p)
    H   = (Qv(5) + p) / rho
  end subroutine primitives

  ! ---------------- Roe flux (same as before, with HH2 entropy fix) ---------
  subroutine flux_roe(qL, qR, Fface)
    real(dp), intent(in)  :: qL(0:nx,js-1:je,ks-1:ke,nv), qR(0:nx,js-1:je,ks-1:ke,nv)
    real(dp), intent(out) :: Fface(0:nx,js-1:je,ks-1:ke,nv)
    integer :: i,j,k
    real(dp) :: FL(nv), FR(nv), QLv(nv), QRv(nv)
    real(dp) :: rL, uL, vL, wL, pL, HL
    real(dp) :: rR, uR, vR, wR, pR, HR
    real(dp) :: rroe, uroe, vroe, wroe, Hroe, aroe
    real(dp) :: Rmat(5,5), Lmat(5,5), lam(5), abslam(5), alpha(5), dQ(5)
    real(dp) :: aL, aR, dl1, dl2, dl3, dl4, dl5

    !$OMP DO collapse(2) schedule(static)
    do k = ks-1, ke
    do j = js-1, je
    do i = 0, nx
      QLv = qL(i,j,k,:)
      QRv = qR(i,j,k,:)
      !write(*,*) QLv(1),QRv(1)

      ! primitives L
      call primitives(QLv, rL, uL, vL, wL, pL, HL)
      ! primitives R
      call primitives(QRv, rR, uR, vR, wR, pR, HR)

      ! physical fluxes
      FL(1) = rL*uL
      FL(2) = rL*uL*uL + pL
      FL(3) = rL*uL*vL     
      FL(4) = rL*uL*wL     
      FL(5) = uL*(QLv(5) + pL)

      FR(1) = rR*uR
      FR(2) = rR*uR*uR + pR
      FR(3) = rR*uR*vR
      FR(4) = rR*uR*wR
      FR(5) = uR*(QRv(5) + pR)

      ! Roe averages
      rroe = sqrt(rL*rR)
      uroe = (sqrt(rL)*uL + sqrt(rR)*uR) / (sqrt(rL)+sqrt(rR))
      vroe = (sqrt(rL)*vL + sqrt(rR)*vR) / (sqrt(rL)+sqrt(rR))
      wroe = (sqrt(rL)*wL + sqrt(rR)*wR) / (sqrt(rL)+sqrt(rR))
      Hroe = (sqrt(rL)*HL  + sqrt(rR)*HR)  / (sqrt(rL)+sqrt(rR))
      aroe = sqrt(max((gamma-1.0_dp)*(Hroe - 0.5_dp*(uroe*uroe+vroe*vroe+wroe*wroe)), small_p/rroe))

      call eigen_matrices(uroe, vroe, wroe, Hroe, aroe, Rmat, Lmat, lam)

      ! HH2 entropy fix
      aL = sqrt(gamma*pL/rL); aR = sqrt(gamma*pR/rR)
      dl1 = max(0.0_dp, (uR-aR) - (uL-aL))
      dl2 = max(0.0_dp,  uR     -  uL   )
      dl3 = max(0.0_dp,  uR     -  uL   )
      dl4 = max(0.0_dp,  uR     -  uL   )
      dl5 = max(0.0_dp, (uR+aR) - (uL+aL))
      abslam(1) = hh_fix(lam(1), dl1)
      abslam(2) = hh_fix(lam(2), dl2)
      abslam(3) = hh_fix(lam(3), dl3)
      abslam(4) = hh_fix(lam(4), dl4)
      abslam(5) = hh_fix(lam(5), dl5)

      dQ    = QRv - QLv
      alpha = matmul(Lmat, dQ)
      Fface(i,j,k,:) = 0.5_dp*(FL + FR) - 0.5_dp*   &
                     matmul(Rmat, (/abslam(1)*alpha(1), abslam(2)*alpha(2), abslam(3)*alpha(3), &
                                                        abslam(4)*alpha(4), abslam(5)*alpha(5)/))

    end do;end do;end do
    !$OMP end do

  end subroutine flux_roe

  function hh_fix(lambda, delta) result(val)
    real(dp), intent(in) :: lambda, delta
    real(dp) :: val
    if (abs(lambda) < delta) then
      val = 0.5_dp*(lambda*lambda/delta + delta)
    else
      val = abs(lambda)
    end if
  end function hh_fix

  subroutine eigen_matrices(u, v, w, H, a, R, L, lam)
    real(dp), intent(in)  :: u, v, w, H, a
    real(dp), intent(out) :: R(5,5), L(5,5), lam(5)
    real(dp) :: q2, a2, b1, b2, g1
    
    g1 = gamma - 1.0_dp
    q2 = u*u + v*v + w*w
    a2 = a*a

    lam(1) = u - a
    lam(2) = u
    lam(3) = u
    lam(4) = u
    lam(5) = u + a

    ! --- Right eigenvectors ---
    R = 0.0_dp
    ! acoustic (u-a)
    R(1,1) = 1.0_dp
    R(2,1) = u - a
    R(3,1) = v
    R(4,1) = w
    R(5,1) = H - u*a
    ! entropy
    R(1,2) = 1.0_dp
    R(2,2) = u
    R(3,2) = v
    R(4,2) = w
    R(5,2) = 0.5_dp*q2
    ! shear-y
    R(1,3) = 0.0_dp
    R(2,3) = 0.0_dp
    R(3,3) = 1.0_dp
    R(4,3) = 0.0_dp
    R(5,3) = v
    ! shear-z
    R(1,4) = 0.0_dp
    R(2,4) = 0.0_dp
    R(3,4) = 0.0_dp
    R(4,4) = 1.0_dp
    R(5,4) = w
    ! acoustic (u+a)
    R(1,5) = 1.0_dp
    R(2,5) = u + a
    R(3,5) = v
    R(4,5) = w
    R(5,5) = H + u*a

    ! --- Left eigenvectors (L = R^-1) ---
    L = 0.0_dp
    b1 = g1*q2/(2.0_dp*a2)
    b2 = g1/a2

    L(1,1) = 0.5_dp*(b1 + u/a)
    L(1,2) = -0.5_dp*(b2*u + 1.0_dp/a)
    L(1,3) = -0.5_dp*b2*v
    L(1,4) = -0.5_dp*b2*w
    L(1,5) = 0.5_dp*b2

    L(2,1) = 1.0_dp - b1
    L(2,2) = b2*u
    L(2,3) = b2*v
    L(2,4) = b2*w
    L(2,5) = -b2

    L(3,1) = -v
    L(3,3) = 1.0_dp

    L(4,1) = -w
    L(4,4) = 1.0_dp

    L(5,1) = 0.5_dp*(b1 - u/a)
    L(5,2) = 0.5_dp*(-b2*u + 1.0_dp/a)
    L(5,3) = -0.5_dp*b2*v
    L(5,4) = -0.5_dp*b2*w
    L(5,5) = 0.5_dp*b2
  end subroutine eigen_matrices

  ! ------------------------- WENO-Z weights ---------------------------
  pure function weno5z_left(v1,v2,v3,v4,v5) result(vf)
    real(dp), intent(in) :: v1,v2,v3,v4,v5
    real(dp) :: vf, p0,p1,p2, b0,b1,b2, a0,a1,a2, s, tau5
    real(dp), parameter :: eps = 1.0d-20, d0=1.0d0/10.0d0, d1=6.0d0/10.0d0, d2=3.0d0/10.0d0
    p0 = ( 2.0_dp*v1 - 7.0_dp*v2 + 11.0_dp*v3)/6.0_dp
    p1 = (-1.0_dp*v2 + 5.0_dp*v3 +  2.0_dp*v4)/6.0_dp
    p2 = ( 2.0_dp*v3 + 5.0_dp*v4 -  1.0_dp*v5)/6.0_dp
    b0 = (13.0_dp/12.0_dp)*(v1 - 2.0_dp*v2 + v3)**2 + 0.25_dp*(v1 - 4.0_dp*v2 + 3.0_dp*v3)**2
    b1 = (13.0_dp/12.0_dp)*(v2 - 2.0_dp*v3 + v4)**2 + 0.25_dp*(v2 - v4)**2
    b2 = (13.0_dp/12.0_dp)*(v3 - 2.0_dp*v4 + v5)**2 + 0.25_dp*(3.0_dp*v3 - 4.0_dp*v4 + v5)**2
    tau5 = abs(b0 - b2)
    a0 = d0 * (1.0_dp + (tau5/(b0+eps))**2)
    a1 = d1 * (1.0_dp + (tau5/(b1+eps))**2)
    a2 = d2 * (1.0_dp + (tau5/(b2+eps))**2)
    s  = a0 + a1 + a2
    vf = (a0*p0 + a1*p1 + a2*p2) / s
  end function weno5z_left

  pure function weno5z_right(v1,v2,v3,v4,v5) result(vf)
    real(dp), intent(in) :: v1,v2,v3,v4,v5
    real(dp) :: vf, p0,p1,p2, b0,b1,b2, a0,a1,a2, s, tau5
    real(dp), parameter :: eps = 1.0d-20, d0=1.0d0/10.0d0, d1=6.0d0/10.0d0, d2=3.0d0/10.0d0
    p0 = (-1.0_dp*v1 + 5.0_dp*v2 +  2.0_dp*v3)/6.0_dp
    p1 = ( 2.0_dp*v2 + 5.0_dp*v3 -  1.0_dp*v4)/6.0_dp
    p2 = (11.0_dp*v3 - 7.0_dp*v4 +  2.0_dp*v5)/6.0_dp
    b0 = (13.0_dp/12.0_dp)*(v1 - 2.0_dp*v2 + v3)**2 + 0.25_dp*(v1 - 4.0_dp*v2 + 3.0_dp*v3)**2
    b1 = (13.0_dp/12.0_dp)*(v2 - 2.0_dp*v3 + v4)**2 + 0.25_dp*(v2 - v4)**2
    b2 = (13.0_dp/12.0_dp)*(v3 - 2.0_dp*v4 + v5)**2 + 0.25_dp*(3.0_dp*v3 - 4.0_dp*v4 + v5)**2
    tau5 = abs(b0 - b2)
    a0 = d0 * (1.0_dp + (tau5/(b0+eps))**2)
    a1 = d1 * (1.0_dp + (tau5/(b1+eps))**2)
    a2 = d2 * (1.0_dp + (tau5/(b2+eps))**2)
    s  = a0 + a1 + a2
    vf = (a0*p0 + a1*p1 + a2*p2) / s
  end function weno5z_right


  subroutine compute_flux(Fface,R,dir)
    real(dp), intent(in)  :: Fface(0:nx,js-1:je,ks-1:ke,nv)
    integer , intent(in)  :: dir
    real(dp), intent(out) :: R(1-nghost:nx+nghost,js-nghost:je+nghost,ks-nghost:ke+nghost,nv)
    real(dp) :: dR(5)
    integer :: i,j,k

    If(dir==1) then
      !$OMP DO schedule(static)
      do k = ks, ke
      do j = js, je
      do i = 1, nx
        dR(:) = ( area_x(i,j,k)*Fface(i,j,k,:) - area_x(i-1,j,k)*Fface(i-1,j,k,:) ) / vol(i,j,k)
        R(i,j,k,1) = R(i,j,k,1)-dR(1)
        R(i,j,k,2) = R(i,j,k,2)-dR(2)
        R(i,j,k,3) = R(i,j,k,3)-dR(3)
        R(i,j,k,4) = R(i,j,k,4)-dR(4)
        R(i,j,k,5) = R(i,j,k,5)-dR(5)
      end do; end do; end do
      !$OMP end do

    else if(dir==2) then
      !$OMP DO schedule(static)
      do k = ks, ke
      do j = js, je
      do i = 1, nx
        dR(:) = ( area_y(i,j,k)*Fface(i,j,k,:) - area_y(i,j-1,k)*Fface(i,j-1,k,:) ) / vol(i,j,k)
        R(i,j,k,1) = R(i,j,k,1)-dR(1)
        R(i,j,k,2) = R(i,j,k,2)+dR(3)
        R(i,j,k,3) = R(i,j,k,3)-dR(2)
        R(i,j,k,4) = R(i,j,k,4)-dR(4)
        R(i,j,k,5) = R(i,j,k,5)-dR(5)
      end do; end do; end do
      !$OMP end do

    else if(dir==3) then
      !$OMP DO schedule(static)
      do k = ks, ke
      do j = js, je
      do i = 1, nx
        dR(:) = ( area_z(i,j,k)*Fface(i,j,k,:) - area_z(i,j,k-1)*Fface(i,j,k-1,:) ) / vol(i,j,k)
        R(i,j,k,1) = R(i,j,k,1)-dR(1)
        R(i,j,k,2) = R(i,j,k,2)+dR(4)
        R(i,j,k,3) = R(i,j,k,3)-dR(3)
        R(i,j,k,4) = R(i,j,k,4)-dR(2)
        R(i,j,k,5) = R(i,j,k,5)-dR(5)
      end do; end do; end do
      !$OMP end do
    end if

  end subroutine compute_flux










































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


end program shock_tube_fv_weno_roe_charZ