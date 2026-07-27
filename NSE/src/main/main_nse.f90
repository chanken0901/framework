program main
  use module_mpi
  use mod_precision, only : dp
  use mod_constants
  use mod_common_config, only : simulation_config, init_simulation_config, print_simulation_config
  use mod_model_config, only : nse_config, init_nse_config, print_nse_config
  use mod_input_reader, only : read_all_inputs
  use mod_grid_fvm, only : build_uniform_grid
  use mod_nse_field, only : allocate_nse_fields, deallocate_nse_fields, &
                            Q, Q0, RHS, F, Qw
  use mod_slf_output, only : write_nse_conserved_slf,write_meta_json

  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse

  !--- MPI ---
  integer js,je,ks,ke,ierror,ierr,ierf
  !-----------

  !--- MPI ---
  Call mpi_init(ierror)
  Call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  Call mpi_comm_rank(mpi_comm_world,my_rank,ierr)

  call init_simulation_config(sim)
  call init_nse_config(nse)
  call read_all_inputs('input.dat', sim, nse=nse)

  call print_simulation_config(sim)
  call print_nse_config(nse)

  Call mp_setup_division(sim%nx, sim%ny, sim%nz)

  js=j_sta
  je=j_end
  ks=k_sta
  ke=k_end
!write(*,*) my_rank,js,je,ks,ke
!stop
  !-----------


  call build_uniform_grid(sim, js, je, ks, ke)
  !call build_grid_uniform(sim%x_min, sim%x_max, sim%y_min, sim%y_max, sim%z_min, sim%z_max)
  if(my_rank == root) write(*,*) "Complete build grid"
  call allocate_nse_fields(sim, nse, js, je, ks, ke)
  if(my_rank == root) write(*,*) "Complete allocate"
  call initialize_Sym()
  if(my_rank == root) write(*,*) "Complete initialize"
  sim%t = 0.0_dp; sim%step = 0; sim%ttotal = 0.0_dp

  !$OMP parallel default(none)         &
  !$OMP & shared(Q,Q0,QL,QR,Qw,RHS,F,  &
  !$OMP &        ks,ke,js,je,          &
  !$OMP &        my_rank,              &
  !$OMP &        sim,nse               )

  do while (sim%t < sim%t_max)
    !$OMP masked
    if(my_rank == root) write(*,*) sim%step,sim%t,sim%dt,sim%t2-sim%t1,sim%ttotal
    sim%t1 = MPI_Wtime()
    !$OMP end masked
    !$OMP barrier

    !$OMP single
    call compute_dt(Q, sim%dt)
    !$OMP end single

    !$OMP masked
    if (sim%t + sim%dt > sim%t_max) sim%dt = sim%t_max - sim%t
    !$OMP end masked
    !$OMP barrier

    call step_rk3(Q, sim%dt)

    !$OMP masked
    sim%t = sim%t + sim%dt; sim%step = sim%step + 1
    if (mod(sim%step, sim%output_frequency) == 0) then
      !call write_vtk_data(step, my_rank, Q(1:sim%nx,js:je,ks:ke,:), t) ! 3D VTKデータを出力
      !call write_bin_data(step, my_rank, Q(1:sim%nx,js:je,ks:ke,:), t) ! 3D VTKデータを出力
      call write_nse_conserved_slf(sim, sim%step, sim%t, Q, rank=my_rank)
    end if
    sim%t2 = MPI_Wtime()
    sim%ttotal = sim%ttotal+(sim%t2-sim%t1)
    call mp_barrier
    !$OMP end masked
    !$OMP barrier


  end do

  !$OMP end parallel

  call mp_stop(ierr)

contains


  subroutine allocate_fields()
    allocate(Q  (1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))
    allocate(Q0 (1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))
    allocate(RHS(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))
    allocate(F (0:sim%nx, js-1:je, ks-1:ke, nse%nv))
    !allocate(QL(0:sim%nx, js-1:je, ks-1:ke, nse%nv), QR(0:sim%nx, js-1:je, ks-1:ke, nse%nv))
    !allocate(Q_vis(1:sim%nx,1:sim%ny,1:sim%nz,nse%nv))
    allocate(Qw(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv))
  end subroutine allocate_fields

  subroutine initialize_Sym()
  use mod_grid_fvm!, only : build_uniform_grid

    integer :: i,j,k
    real(dp) :: rho, u, v, w, p
    real(dp) :: C_1,C_2
    real(dp) :: x2,y2,z2

    real(dp),parameter :: rho_0 = 1.0_dp
    real(dp),parameter :: Mach  = 0.5_dp

      C_1=1.0_dp/nse%gamma
      C_2=(rho_0*Mach*Mach)/16.0_dp

      do k = ks, ke
      do j = js, je
      do i = 1, sim%nx
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
        Q(i,j,k,5) = p/(nse%gamma-1.0_dp) + 0.5_dp*rho*(u*u+v*v+w*w)
      end do; end do; end do
    call apply_bc(Q)

    !call write_bin_data(0, my_rank, Q(1:sim%nx,js:je,ks:ke,:), 0.0_dp) ! 3D VTKデータを出力
    call write_meta_json(sim, is=1, ie=sim%nx, js=js, je=je, ks=ks, ke=ke, use_cuda=.false.)
    call write_nse_conserved_slf(sim, 0, 0.0_dp, Q, rank=my_rank)

  end subroutine initialize_Sym


  subroutine apply_bc(A)
    real(dp), intent(inout) :: A(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    integer :: i,j,k,g

    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do g = 1, sim%nghost
      A( 1-g,j,k,:)  = A((sim%nx+1)-g,j,k,:)
      A(sim%nx+g,j,k,:)  = A(     0+g,j,k,:)
    end do; end do; end do
    !$OMP END DO

    !$OMP masked
    call BC_Periodic_y_dir_MPI_R8(A)
    !$OMP end masked
    !$OMP barrier

    !$OMP masked
    call BC_Periodic_z_dir_MPI_R8(A)
    !$OMP end masked
    !$OMP barrier

    !$OMP masked
    call mp_send_recv_pre_r8_Vec(A,sim%nghost,1-sim%nghost,sim%nx+sim%nghost,js-sim%nghost,je+sim%nghost,ks-sim%nghost,ke+sim%nghost)
    !$OMP end masked
    !$OMP barrier


  end subroutine apply_bc


  Subroutine BC_Periodic_y_dir_MPI_R8( F1 )

    Use module_mpi
    Implicit none
    real(dp),intent(inout) :: F1(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
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

    If(js==1.or.je==sim%ny) then
      !Allocate(dumcomzx_r8(1:sim%nx,sim%nghost,ks:ke,nse%nv))
      !Allocate(dumcomzx_s8(1:sim%nx,sim%nghost,ks:ke,nse%nv))
      allocate(sendbuf(1:sim%nx,sim%nghost,ks:ke,nse%nv))
      allocate(recv_hi(1:sim%nx,sim%nghost,ks:ke,nse%nv))
      allocate(recv_lo(1:sim%nx,sim%nghost,ks:ke,nse%nv))
      dum_len=sim%nx*sim%nghost*(ke-ks+1)*nse%nv
      tag=1
      Do icom=0,Ndiv_Nz-1
        If(ks==kksta(icom)) then
          If(je==sim%ny) then
           partner = itable(0, icom)
            Do k = ks,ke
            Do i = 1, sim%nx, 1
              Do g = 1, sim%nghost, 1
              sendbuf(i,g,k,:) = F1(i,sim%ny+g-sim%nghost,k,:)
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
            Do i = 1, sim%nx, 1
              Do g = 1, sim%nghost, 1
              sendbuf(i,g,k,:) = F1(i,g,k,:)
              End Do
            End Do; End Do
            call MPI_Sendrecv( sendbuf(1,1,ks,1), dum_len, MPI_DOUBLE_PRECISION, partner, tag, &
                               recv_lo(1,1,ks,1), dum_len, MPI_DOUBLE_PRECISION, partner, tag, &
                               MPI_COMM_WORLD, istatus, ierr )
            !dum_len=sim%nx*sim%nghost*(ke-ks+1)*nse%nv
            !Call mpi_isend(dumcomzx_s8(1,1,ks,1),dum_len,MPI_DOUBLE_PRECISION,                        &
            !                                  itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(2),ierr)
            !Call mpi_irecv(dumcomzx_r8(1,1,ks,1),dum_len,MPI_DOUBLE_PRECISION,                        &
            !                                  itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(1),ierr)
          End if

          !Call mpi_wait(isend(1),istatus,ierr)
          !Call mpi_wait(isend(2),istatus,ierr)

        End if
      End Do

      If(je==sim%ny) then
        Do k = ks,ke
        Do i = 1, sim%nx, 1
          Do g = 1, sim%nghost, 1
          F1(i,sim%ny+g,k,:)=recv_hi(i,g,k,:)
          End Do
        End Do; End Do
      End If

      If(js==1) then
        Do k = ks,ke
        Do i = 1, sim%nx, 1
          Do g = 1, sim%nghost, 1
          F1(i,1-g,k,:)=recv_lo(i,sim%nghost-g+1,k,:)
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
    real(dp),intent(inout) :: F1(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    real(dp), dimension(:,:,:,:), allocatable :: dumcomxy_s8,dumcomxy_r8

    Integer i, j, k, g ! 変数
    Integer icom

    !--- MPI ---
    integer ierror,ierr
    integer dum_len
    integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
    !-----------

    If(ks==1.or.ke==sim%nz) then
      Allocate(dumcomxy_s8(1:sim%nx,js:je,sim%nghost,1:nse%nv))
      Allocate(dumcomxy_r8(1:sim%nx,js:je,sim%nghost,1:nse%nv))
      Do icom=0,Ndiv_Ny-1
        If(j_sta==jjsta(icom)) then

          If(ke==sim%nz) then
            Do j = js,je
            Do i = 1, sim%nx, 1
              Do g = 1, sim%nghost, 1
              dumcomxy_s8(i,j,g,:) = F1(i,j,sim%nz+g-sim%nghost,:)
              End Do
            End Do; End Do
            !dum_len=sim%nx*sim%nghost*(ke-ks+1)*nse%nv
            dum_len = sim%nx * (je-js+1) * sim%nghost * nse%nv
            Call mpi_isend(dumcomxy_s8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,0),1,MPI_COMM_WORLD,isend(1),ierr)
            Call mpi_irecv(dumcomxy_r8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,0),1,MPI_COMM_WORLD,isend(2),ierr)
          End If

          If(ks==1) then
            Do j = js,je
            Do i = 1, sim%nx, 1
              Do g = 1, sim%nghost, 1
              dumcomxy_s8(i,j,g,:) = F1(i,j,g,:)
              End Do
            End Do; End Do
            !dum_len=sim%nx*sim%nghost*(ke-ks+1)*nse%nv
            dum_len = sim%nx * (je-js+1) * sim%nghost * nse%nv
            Call mpi_isend(dumcomxy_s8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(2),ierr)
            Call mpi_irecv(dumcomxy_r8(1,js,1,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                              itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(1),ierr)
          End if

          Call mpi_wait(isend(1),istatus,ierr)
          Call mpi_wait(isend(2),istatus,ierr)

        End if
      End Do


      If(ke==sim%nz) then
        Do j = js,je
        Do i = 1, sim%nx, 1
          Do g = 1, sim%nghost, 1
          F1(i,j,sim%nz+g,:) = dumcomxy_r8(i,j,g,:)
          End Do
        End Do; End Do
      End If

      If(ks==1) then
        Do j = js,je
        Do i = 1, sim%nx, 1
          Do g = 1, sim%nghost, 1
          F1(i,j,1-g,:) = dumcomxy_r8(i,j,sim%nghost-g+1,:)
          End Do
        End Do; End Do
      End if

      Deallocate(dumcomxy_r8)
      Deallocate(dumcomxy_s8)
    End If

    Return
  End Subroutine


  subroutine compute_dt(Qin, dt)
    real(dp), intent(in)  :: Qin(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    real(dp), intent(out) :: dt
    integer :: i,j,k
    real(dp) :: rho,u,v,w,p,a,maxs

    maxs = 0.0_dp

    do k = ks, ke
    do j = js, je
    do i = 1, sim%nx
      rho = max(Qin(i,j,k,1), nse%small_rho)
      u   = Qin(i,j,k,2)/rho
      v   = Qin(i,j,k,3)/rho
      w   = Qin(i,j,k,4)/rho
      p   = max( (nse%gamma-1.0_dp)*(Qin(i,j,k,5) - 0.5_dp*rho*(u*u+v*v+w*w)), nse%small_p )
      a   = sqrt(nse%gamma*p/rho)
      maxs = max(maxs, abs(u)+a, abs(v)+a, abs(w)+a)
    end do;end do;end do

    !dt = nse%cfl * min( dx_min/maxs, dy_min/maxs , dz_min/maxs )
    dt = nse%cfl * min( sim%dx/maxs, sim%dy/maxs , sim%dz/maxs )

    call mp_barrier
    call mp_allminr8(dt)

  end subroutine compute_dt

  subroutine step_rk3(Qinout, dt)
    real(dp), intent(inout) :: Qinout(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    real(dp), intent(in)    :: dt
    integer :: i,j,k

    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost,ke+sim%nghost
    do j = js-sim%nghost,je+sim%nghost
    do i = 1-sim%nghost,sim%nx+sim%nghost
      Q0(i,j,k,:) = Qinout(i,j,k,:)
    end do; end do; end do
    !$OMP end do

    ! --- 1st step --- !
    call compute_rhs(Qinout, RHS)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost,ke+sim%nghost
    do j = js-sim%nghost,je+sim%nghost
    do i = 1-sim%nghost,sim%nx+sim%nghost
      Qinout(i,j,k,:) = Q0(i,j,k,:)+dt*RHS(i,j,k,:)
    end do; end do; end do
    !$OMP end do

    ! --- 2nd step --- !
    call compute_rhs(Qinout, RHS)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost,ke+sim%nghost
    do j = js-sim%nghost,je+sim%nghost
    do i = 1-sim%nghost,sim%nx+sim%nghost
      Qinout(i,j,k,:) = 0.75_dp*Q0(i,j,k,:) + 0.25_dp*(Qinout(i,j,k,:) + dt*RHS(i,j,k,:))
    end do; end do; end do
    !$OMP end do

    ! --- 3rd step --- !
    call compute_rhs(Qinout, RHS)
    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost,ke+sim%nghost
    do j = js-sim%nghost,je+sim%nghost
    do i = 1-sim%nghost,sim%nx+sim%nghost
      Qinout(i,j,k,:) = (1.0_dp/3.0_dp)*Q0(i,j,k,:) + (2.0_dp/3.0_dp)*(Qinout(i,j,k,:) + dt*RHS(i,j,k,:))
    end do; end do; end do
    !$OMP end do

  end subroutine step_rk3

  subroutine compute_rhs(Qin, R)
    use mod_grid_fvm!, only : build_uniform_grid
    real(dp), intent(inout)  :: Qin(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    real(dp), intent(out) :: R  (1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    integer :: i,j,k,g
    integer :: dir
    real(dp) :: dR(nse%nv)

    !$OMP DO collapse(2) schedule(static)
    do k = ks-sim%nghost,ke+sim%nghost
    do j = js-sim%nghost,je+sim%nghost
    do i = 1-sim%nghost,sim%nx+sim%nghost
      R (i,j,k,:) = 0.0_dp
    end do; end do; end do
    !$OMP end do
    call apply_bc(Qin)

    dir = 1
    call flux_KEEP(Qin, F, dir)
    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
    do j = js, je
    do i = 1, sim%nx
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
    do i = 1, sim%nx
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
    do i = 1, sim%nx
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
    real(dp), intent(in)  :: Qin(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    integer , intent(in)  :: direction
    real(dp), intent(out) :: Qout(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    integer :: i, j, k
    if (direction == 1) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
      do i = 1-sim%nghost, sim%nx+sim%nghost
        Qout(i,j,k,1) =  Qin(i,j,k,1)
        Qout(i,j,k,2) =  Qin(i,j,k,2)
        Qout(i,j,k,3) =  Qin(i,j,k,3)
        Qout(i,j,k,4) =  Qin(i,j,k,4)
        Qout(i,j,k,5) =  Qin(i,j,k,5)
      end do; end do; end do
      !$OMP end do

    else if (direction == 2) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
      do i = 1-sim%nghost, sim%nx+sim%nghost
        Qout(i,j,k,1) =  Qin(i,j,k,1)
        Qout(i,j,k,2) =  Qin(i,j,k,3)
        Qout(i,j,k,3) = -Qin(i,j,k,2)
        Qout(i,j,k,4) =  Qin(i,j,k,4)
        Qout(i,j,k,5) =  Qin(i,j,k,5)
      end do; end do; end do
      !$OMP end do

    else if (direction == 3) then
      !$OMP DO collapse(2) schedule(static)
      do k = ks-sim%nghost, ke+sim%nghost
      do j = js-sim%nghost, je+sim%nghost
      do i = 1-sim%nghost, sim%nx+sim%nghost
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
    real(dp), intent(in)  :: Qin(1-sim%nghost:sim%nx+sim%nghost,js-sim%nghost:je+sim%nghost,ks-sim%nghost:ke+sim%nghost,nse%nv)
    integer , intent(in)  :: direction
    real(dp), intent(out) :: Fface(0:sim%nx,js-1:je,ks-1:ke,nse%nv)
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
       do i = 0, sim%nx
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

         rm1 = max(Qm1, nse%small_rho)
         um1 = Qm2/rm1
         vm1 = Qm3/rm1
         wm1 = Qm4/rm1
         pm1 = max((nse%gamma-1.0_dp)*(Qm5 - 0.5_dp*rm1*(um1*um1+vm1*vm1+wm1*wm1)), nse%small_p)
         Hm1 = (Qm5 + pm1) / rm1

         rp1 = max(Qp1, nse%small_rho)
         up1 = Qp2/rp1
         vp1 = Qp3/rp1
         wp1 = Qp4/rp1
         pp1 = max((nse%gamma-1.0_dp)*(Qp5 - 0.5_dp*rp1*(up1*up1+vp1*vp1+wp1*wp1)), nse%small_p)
         Hp1 = (Qp5 + pp1) / rp1

         rmp=0.5_dp*(rm1+rp1)
         ump=0.5_dp*(um1+up1)
         vmp=0.5_dp*(vm1+vp1)
         wmp=0.5_dp*(wm1+wp1)
         Lmp=0.5_dp*((pm1/rm1)+(pp1/rp1))/(nse%gamma-1.0_dp)

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
       do i = 0, sim%nx
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

         rm1 = max(Qm1, nse%small_rho)
         um1 = Qm2/rm1
         vm1 = Qm3/rm1
         wm1 = Qm4/rm1
         pm1 = max((nse%gamma-1.0_dp)*(Qm5 - 0.5_dp*rm1*(um1*um1+vm1*vm1+wm1*wm1)), nse%small_p)
         Hm1 = (Qm5 + pm1) / rm1

         rp1 = max(Qp1, nse%small_rho)
         up1 = Qp2/rp1
         vp1 = Qp3/rp1
         wp1 = Qp4/rp1
         pp1 = max((nse%gamma-1.0_dp)*(Qp5 - 0.5_dp*rp1*(up1*up1+vp1*vp1+wp1*wp1)), nse%small_p)
         Hp1 = (Qp5 + pp1) / rp1

         rmp=0.5_dp*(rm1+rp1)
         ump=0.5_dp*(um1+up1)
         vmp=0.5_dp*(vm1+vp1)
         wmp=0.5_dp*(wm1+wp1)
         Lmp=0.5_dp*((pm1/rm1)+(pp1/rp1))/(nse%gamma-1.0_dp)

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
       do i = 0, sim%nx
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

         rm1 = max(Qm1, nse%small_rho)
         um1 = Qm2/rm1
         vm1 = Qm3/rm1
         wm1 = Qm4/rm1
         pm1 = max((nse%gamma-1.0_dp)*(Qm5 - 0.5_dp*rm1*(um1*um1+vm1*vm1+wm1*wm1)), nse%small_p)
         Hm1 = (Qm5 + pm1) / rm1

         rp1 = max(Qp1, nse%small_rho)
         up1 = Qp2/rp1
         vp1 = Qp3/rp1
         wp1 = Qp4/rp1
         pp1 = max((nse%gamma-1.0_dp)*(Qp5 - 0.5_dp*rp1*(up1*up1+vp1*vp1+wp1*wp1)), nse%small_p)
         Hp1 = (Qp5 + pp1) / rp1

         rmp=0.5_dp*(rm1+rp1)
         ump=0.5_dp*(um1+up1)
         vmp=0.5_dp*(vm1+vp1)
         wmp=0.5_dp*(wm1+wp1)
         Lmp=0.5_dp*((pm1/rm1)+(pp1/rp1))/(nse%gamma-1.0_dp)

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









































  !subroutine write_vtk_data(step, my_rank, Qin, t)
  !  integer, intent(in) :: step, my_rank
  !  real(dp), intent(in) :: Qin(1:sim%nx, js:je, ks:ke, nse%nv)
  !  real(dp), intent(in) :: t
  !  integer :: i, j, k
  !  real(dp) :: rho, u, v, w, p
  !  real(dp) :: xm,xp,ym,yp,zm,zp
  !  real(dp) :: xc1,xc2,yc1,yc2,zc1,zc2
  !  real(dp) :: xc,yc,zc
  !  character(len=256) :: fname
  !  character(len=*), parameter :: vtk_fmt = '(ES22.12E3)'

  !  ! VTKファイルはsim%output_frequencyごとにのみ書き出す
  !  write(fname, '(A,I0.5,A,I0.5,A)') '3d_result_step_', step, '_rank', my_rank, '.vtk'
  !  open(unit=30, file=fname, status='replace')

  !  ! VTK Header
  !  write(30, '(A)') '# vtk DataFile Version 2.0'
  !  write(30, '(A, F12.6)') 'Time = ', t
  !  write(30, '(A)') 'ASCII'
  !  write(30, '(A)') 'DATASET STRUCTURED_GRID'
  !  write(30, '(A,I6,I6,I6)') 'DIMENSIONS ', sim%nx, (je-js)+1, (ke-ks)+1
  !  write(30, '(A,I12,A)') 'POINTS ', sim%nx*((je-js)+1)*((ke-ks)+1), ' double'

  !  xm = sim%x_min
  !  xp = sim%x_max
  !  ym = sim%y_min
  !  yp = sim%y_max
  !  zm = sim%z_min
  !  zp = sim%z_max

  !  ! Write grid points
  !  do k = ks, ke
  !    do j = js, je
  !      do i = 1, sim%nx

  !        !xc1 = xm + (xp-xm) * real(i-1,dp) / real(sim%nx,dp)
  !        !xc2 = xm + (xp-xm) * real(i  ,dp) / real(sim%nx,dp)
  !        !yc1 = ym + (yp-ym) * real(j-1,dp) / real(sim%ny,dp)
  !        !yc2 = ym + (yp-ym) * real(j  ,dp) / real(sim%ny,dp)
  !        !zc1 = zm + (zp-zm) * real(k-1,dp) / real(sim%nz,dp)
  !        !zc2 = zm + (zp-zm) * real(k  ,dp) / real(sim%nz,dp)

  !        !xc = 0.5_dp*(xc1+xc2)
  !        !yc = 0.5_dp*(yc1+yc2)
  !        !zc = 0.5_dp*(zc1+zc2)

  !        !write(30, vtk_fmt) xc, yc, zc
  !        write(30, vtk_fmt) x_cell(i,j,k), y_cell(i,j,k), z_cell(i,j,k)
  !      end do
  !    end do
  !  end do

  !  ! Write data
  !  write(30, '(A,I12)') 'POINT_DATA ', sim%nx*((je-js)+1)*((ke-ks)+1)

  !  ! --- Density (rho) ---
  !  write(30, '(A)') 'SCALARS rho double 1'
  !  write(30, '(A)') 'LOOKUP_TABLE default'
  !  do k = ks, ke
  !    do j = js, je
  !      do i = 1, sim%nx
  !        rho = Qin(i,j,k,1)
  !        write(30, vtk_fmt) rho
  !      end do
  !    end do
  !  end do

  !  ! --- X-Velocity (u) ---
  !  write(30, '(A)') 'SCALARS u double 1'
  !  write(30, '(A)') 'LOOKUP_TABLE default'
  !  do k = ks, ke
  !    do j = js, je
  !      do i = 1, sim%nx
  !        rho = max(Qin(i,j,k,1), nse%small_rho)
  !        u   = Qin(i,j,k,2) / rho
  !        write(30, vtk_fmt) u
  !      end do
  !    end do
  !  end do

  !  ! --- Y-Velocity (v) ---
  !  write(30, '(A)') 'SCALARS v double 1'
  !  write(30, '(A)') 'LOOKUP_TABLE default'
  !  do k = ks, ke
  !    do j = js, je
  !      do i = 1, sim%nx
  !        rho = max(Qin(i,j,k,1), nse%small_rho)
  !        v   = Qin(i,j,k,3) / rho
  !        write(30, vtk_fmt) v
  !      end do
  !    end do
  !  end do

  !  ! --- Z-Velocity (w) ---
  !  write(30, '(A)') 'SCALARS w double 1'
  !  write(30, '(A)') 'LOOKUP_TABLE default'
  !  do k = ks, ke
  !    do j = js, je
  !      do i = 1, sim%nx
  !        rho = max(Qin(i,j,k,1), nse%small_rho)
  !        w   = Qin(i,j,k,4) / rho
  !        write(30, vtk_fmt) w
  !      end do
  !    end do
  !  end do

  !  ! --- Pressure (p) ---
  !  write(30, '(A)') 'SCALARS p double 1'
  !  write(30, '(A)') 'LOOKUP_TABLE default'
  !  do k = ks, ke
  !    do j = js, je
  !      do i = 1, sim%nx
  !        rho = max(Qin(i,j,k,1), nse%small_rho)
  !        u   = Qin(i,j,k,2) / rho
  !        v   = Qin(i,j,k,3) / rho
  !        w   = Qin(i,j,k,4) / rho
  !        p   = (nse%gamma-1.0_dp) * (Qin(i,j,k,5) - 0.5_dp*rho*(u*u + v*v + w*w))
  !        p   = max(p, nse%small_p)
  !        write(30, vtk_fmt) p
  !      end do
  !    end do
  !  end do

  !  close(30)
  !  write(*,'(A,A)') 'Wrote VTK: ', trim(fname)
  !end subroutine write_vtk_data


  !subroutine write_bin_data(step, my_rank, Qin, t)
  !  use, intrinsic :: iso_fortran_env, only: int32
  !  implicit none
  !  integer, intent(in) :: step, my_rank
  !  real(dp), intent(in) :: Qin(1:sim%nx, js:je, ks:ke, nse%nv)
  !  real(dp), intent(in) :: t

  !  integer :: u, ios
  !  character(len=256) :: fname
  !  character(len=8)   :: magic
  !  integer(int32) :: ndim, dtype_code
  !  integer(int32) :: shp(4)
  !  integer(int32) :: meta(6)
  !  real(dp) :: t_write

  !  ! ---- file name (rankごと) ----
  !  write(fname, '(A,I0.5,A,I0.5,A)') 'output/3d_result_step_', step, '_rank', my_rank, '.fbn'

  !  ! ---- header ----
  !  magic = 'FBN1' // char(0) // char(0) // char(0) // char(0)
  !  ndim  = 4_int32
  !  ! shape: (sim%nx, ny_local, nz_local, nse%nv)
  !  shp   = [ int(sim%nx, int32), int((je-js)+1, int32), int((ke-ks)+1, int32), int(nse%nv, int32) ]
  !  dtype_code = 2_int32   ! 1=float32, 2=float64(dp)

  !  ! 追加メタ情報（任意だが解析で便利）:
  !  ! meta = [js, je, ks, ke, step, rank]
  !  meta = [ int(js,int32), int(je,int32), int(ks,int32), int(ke,int32), int(step,int32), int(my_rank,int32) ]
  !  t_write = t

  !  open(newunit=u, file=fname, access='stream', form='unformatted', &
  !       status='replace', action='write', iostat=ios)
  !  if (ios /= 0) then
  !    write(*,'(A,A)') 'ERROR: cannot open binary file: ', trim(fname)
  !    error stop
  !  end if

  !  ! ---- write header ----
  !  write(u) magic
  !  write(u) ndim
  !  write(u) shp
  !  write(u) dtype_code
  !  write(u) meta
  !  write(u) t_write

  !  ! ---- write data (Fortran配列順のまま) ----
  !  write(u) Qin

  !  close(u)
  !  write(*,'(A,A)') 'Wrote BIN: ', trim(fname)
  !end subroutine write_bin_data


end program 