module mod_nse_boundary
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config
  use mod_model_config, only : nse_config
  use module_mpi, only : ndiv_ny, ndiv_nz, itable, jjsta, kksta, j_sta, &
    mp_send_recv_pre_r8_Vec, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, &
    MPI_STATUS_SIZE
  implicit none
  private

  public :: apply_nse_boundary
  public :: validate_boundary_scheme
  public :: boundary_required_ghost_cells
  public :: boundary_scheme_name

contains

  subroutine apply_nse_boundary(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    integer :: g, j, k

    !$OMP DO collapse(2) schedule(static)
    do k = ks, ke
      do j = js, je
        do g = 1, sim%nghost
          q(1-g,j,k,:) = q(sim%nx+1-g,j,k,:)
          q(sim%nx+g,j,k,:) = q(g,j,k,:)
        end do
      end do
    end do
    !$OMP END DO

    !$OMP MASKED
    call apply_periodic_y(q, sim, nse, js, je, ks, ke)
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP MASKED
    call apply_periodic_z(q, sim, nse, js, je, ks, ke)
    !$OMP END MASKED
    !$OMP BARRIER

    !$OMP MASKED
    call mp_send_recv_pre_r8_Vec(q, sim%nghost, &
      1-sim%nghost, sim%nx+sim%nghost, &
      js-sim%nghost, je+sim%nghost, ks-sim%nghost, ke+sim%nghost)
    ! Mixed sixth-order derivatives need edge and corner halo values.  A
    ! second staged exchange propagates the face halos received above into
    ! the transverse ghost layers without enlarging the legacy MPI buffers.
    if (trim(adjustl(nse%viscous_scheme)) == 'central6') then
      call mp_send_recv_pre_r8_Vec(q, sim%nghost, &
        1-sim%nghost, sim%nx+sim%nghost, &
        js-sim%nghost, je+sim%nghost, ks-sim%nghost, ke+sim%nghost)
    end if
    !$OMP END MASKED
    !$OMP BARRIER
  end subroutine apply_nse_boundary

  subroutine apply_periodic_y(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), allocatable :: sendbuf(:,:,:,:), recv_high(:,:,:,:), recv_low(:,:,:,:)
    integer :: g, i, k, column, partner, count, ierr
    integer :: ilo, ihi
    integer :: status(MPI_STATUS_SIZE)

    ! X is not decomposed, so its periodic ghost cells are already available.
    ! Include them in the Y update to populate X-Y edges as well as Y faces.
    ilo = 1 - sim%nghost
    ihi = sim%nx + sim%nghost

    if (ndiv_ny == 1) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,sim%ny+g,k,:) = q(i,g,k,:)
            q(i,1-g,k,:) = q(i,sim%ny+1-g,k,:)
          end do
        end do
      end do
      return
    end if

    if (js /= 1 .and. je /= sim%ny) return

    allocate(sendbuf(ilo:ihi,1:sim%nghost,ks:ke,1:nse%nv))
    allocate(recv_high(ilo:ihi,1:sim%nghost,ks:ke,1:nse%nv))
    allocate(recv_low(ilo:ihi,1:sim%nghost,ks:ke,1:nse%nv))
    count = (ihi-ilo+1) * sim%nghost * (ke-ks+1) * nse%nv

    do column = 0, ndiv_nz-1
      if (ks /= kksta(column)) cycle

      if (je == sim%ny) then
        partner = itable(0,column)
        do k = ks, ke
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,g,k,:) = q(i,sim%ny+g-sim%nghost,k,:)
            end do
          end do
        end do
        call MPI_Sendrecv(sendbuf(ilo,1,ks,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          recv_high(ilo,1,ks,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic y-high exchange failed'
      end if

      if (js == 1) then
        partner = itable(ndiv_ny-1,column)
        do k = ks, ke
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,g,k,:) = q(i,g,k,:)
            end do
          end do
        end do
        call MPI_Sendrecv(sendbuf(ilo,1,ks,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          recv_low(ilo,1,ks,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic y-low exchange failed'
      end if
    end do

    if (je == sim%ny) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,sim%ny+g,k,:) = recv_high(i,g,k,:)
          end do
        end do
      end do
    end if

    if (js == 1) then
      do k = ks, ke
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,1-g,k,:) = recv_low(i,sim%nghost-g+1,k,:)
          end do
        end do
      end do
    end if

    deallocate(sendbuf, recv_high, recv_low)
  end subroutine apply_periodic_y

  subroutine apply_periodic_z(q, sim, nse, js, je, ks, ke)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse
    integer, intent(in) :: js, je, ks, ke
    real(dp), intent(inout) :: q(1-sim%nghost:, js-sim%nghost:, ks-sim%nghost:, :)
    real(dp), allocatable :: sendbuf(:,:,:,:), recv_high(:,:,:,:), recv_low(:,:,:,:)
    integer :: g, i, j, row, partner, count, ierr
    integer :: ilo, ihi, jlo, jhi
    integer :: status(MPI_STATUS_SIZE)

    ! Y is updated before Z.  Carry both the X and Y ghost layers through
    ! the Z update so all periodic edges and corners receive valid values.
    ilo = 1 - sim%nghost
    ihi = sim%nx + sim%nghost
    jlo = js - sim%nghost
    jhi = je + sim%nghost

    if (ndiv_nz == 1) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,j,sim%nz+g,:) = q(i,j,g,:)
            q(i,j,1-g,:) = q(i,j,sim%nz+1-g,:)
          end do
        end do
      end do
      return
    end if

    if (ks /= 1 .and. ke /= sim%nz) return

    allocate(sendbuf(ilo:ihi,jlo:jhi,1:sim%nghost,1:nse%nv))
    allocate(recv_high(ilo:ihi,jlo:jhi,1:sim%nghost,1:nse%nv))
    allocate(recv_low(ilo:ihi,jlo:jhi,1:sim%nghost,1:nse%nv))
    count = (ihi-ilo+1) * (jhi-jlo+1) * sim%nghost * nse%nv

    do row = 0, ndiv_ny-1
      if (j_sta /= jjsta(row)) cycle

      if (ke == sim%nz) then
        partner = itable(row,0)
        do j = jlo, jhi
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,j,g,:) = q(i,j,sim%nz+g-sim%nghost,:)
            end do
          end do
        end do
        call MPI_Sendrecv(sendbuf(ilo,jlo,1,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          recv_high(ilo,jlo,1,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic z-high exchange failed'
      end if

      if (ks == 1) then
        partner = itable(row,ndiv_nz-1)
        do j = jlo, jhi
          do i = ilo, ihi
            do g = 1, sim%nghost
              sendbuf(i,j,g,:) = q(i,j,g,:)
            end do
          end do
        end do
        call MPI_Sendrecv(sendbuf(ilo,jlo,1,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          recv_low(ilo,jlo,1,1), count, MPI_DOUBLE_PRECISION, partner, 1, &
          MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) error stop 'periodic z-low exchange failed'
      end if
    end do

    if (ke == sim%nz) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,j,sim%nz+g,:) = recv_high(i,j,g,:)
          end do
        end do
      end do
    end if

    if (ks == 1) then
      do j = jlo, jhi
        do i = ilo, ihi
          do g = 1, sim%nghost
            q(i,j,1-g,:) = recv_low(i,j,sim%nghost-g+1,:)
          end do
        end do
      end do
    end if

    deallocate(sendbuf, recv_high, recv_low)
  end subroutine apply_periodic_z

  subroutine validate_boundary_scheme(sim, nse)
    type(simulation_config), intent(in) :: sim
    type(nse_config), intent(in) :: nse

    if (trim(adjustl(nse%boundary_condition)) /= boundary_scheme_name()) then
      write(*,'(A,A,A,A)') 'ERROR: executable contains boundary scheme "', &
        boundary_scheme_name(), '", but input requested "', &
        trim(adjustl(nse%boundary_condition)) // '"'
      error stop
    end if
    if (sim%nghost /= boundary_required_ghost_cells()) then
      error stop 'legacy MPI halo exchange currently requires exactly three ghost cells'
    end if
    if (nse%nv /= 5) then
      error stop 'legacy MPI halo exchange currently requires five conserved variables'
    end if
  end subroutine validate_boundary_scheme

  integer function boundary_required_ghost_cells() result(nghost)
    nghost = 3
  end function boundary_required_ghost_cells

  pure function boundary_scheme_name() result(name)
    character(len=32) :: name
    name = 'periodic'
  end function boundary_scheme_name

end module mod_nse_boundary
