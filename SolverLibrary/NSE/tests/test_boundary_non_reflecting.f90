program test_boundary_non_reflecting
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use module_mpi
  use mod_precision, only : dp
  use mod_grid_axis, only : build_sinh_axis
  use mod_grid_fvm, only : axis_x,axis_y,axis_z
  use mod_common_config, only : simulation_config, init_simulation_config
  use mod_model_config, only : nse_config, init_nse_config, &
    nse_boundary_face_count, nse_face_x_min, nse_face_x_max
  use mod_nse_boundary, only : apply_nse_boundary, validate_boundary_scheme
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: q(:,:,:,:)
  real(dp) :: reference_state(5), maximum_error, pressure, rho
  integer :: face, i, j, k, js, je, ks, ke, ierr
  character(len=32) :: test_mode
  call get_command_argument(1,test_mode)

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  call init_simulation_config(sim)
  call init_nse_config(nse)
  sim%nx = 12
  sim%ny = 12
  sim%nz = 12
  sim%nghost = 3
  sim%x_min = 0.0_dp
  sim%x_max = 1.0_dp
  sim%y_min = 0.0_dp
  sim%y_max = 1.0_dp
  sim%z_min = 0.0_dp
  sim%z_max = 1.0_dp
  nse%viscous_scheme = 'central6'
  nse%boundary_face_type = 'periodic'
  nse%boundary_face_type(nse_face_x_min) = 'non_reflecting'
  nse%boundary_face_type(nse_face_x_max) = 'non_reflecting'
  nse%boundary_condition = 'mixed'
  nse%boundary_reference_rho(nse_face_x_min) = 1.0_dp
  nse%boundary_reference_rho(nse_face_x_max) = 1.0_dp
  nse%boundary_reference_velocity(:,nse_face_x_min) = &
    [0.25_dp, 0.02_dp, -0.01_dp]
  nse%boundary_reference_velocity(:,nse_face_x_max) = &
    [0.25_dp, 0.02_dp, -0.01_dp]
  nse%boundary_reference_p(nse_face_x_min) = 1.0_dp/nse%gamma
  nse%boundary_reference_p(nse_face_x_max) = 1.0_dp/nse%gamma
  nse%boundary_relaxation_strength = 0.1_dp
  nse%boundary_length_scale = -1.0_dp
  call prepare_geometry()

  call validate_boundary_scheme(sim, nse)
  call mp_setup_division(sim%nx, sim%ny, sim%nz)
  js = j_sta
  je = j_end
  ks = k_sta
  ke = k_end
  allocate(q(1-sim%nghost:sim%nx+sim%nghost, &
    js-sim%nghost:je+sim%nghost, &
    ks-sim%nghost:ke+sim%nghost, nse%nv))

  call make_conserved(1.0_dp, [0.25_dp, 0.02_dp, -0.01_dp], &
    1.0_dp/nse%gamma, nse%gamma, reference_state)
  q = -huge(1.0_dp)
  do k = ks, ke
    do j = js, je
      do i = 1, sim%nx
        q(i,j,k,:) = reference_state
      end do
    end do
  end do

  !$OMP PARALLEL DEFAULT(SHARED)
  call apply_nse_boundary(q, sim, nse, js, je, ks, ke)
  !$OMP END PARALLEL

  maximum_error = 0.0_dp
  do k = ks-sim%nghost, ke+sim%nghost
    do j = js-sim%nghost, je+sim%nghost
      do i = 1-sim%nghost, sim%nx+sim%nghost
        maximum_error = max(maximum_error, maxval(abs(q(i,j,k,:)-reference_state)))
      end do
    end do
  end do
  if (maximum_error > 1.0e-12_dp) then
    write(*,'(A,I0,A,ES16.8)') 'rank ', my_rank, &
      ' mixed-boundary constant-state error = ', maximum_error
    error stop 'mixed non-reflecting boundary failed to preserve a constant state'
  end if

  q = -huge(1.0_dp)
  do k = ks, ke
    do j = js, je
      do i = 1, sim%nx
        rho = 1.0_dp + 0.01_dp*real(i,dp)
        call make_conserved(rho, [0.25_dp, 0.02_dp, -0.01_dp], &
          1.0_dp/nse%gamma, nse%gamma, q(i,j,k,:))
      end do
    end do
  end do
  !$OMP PARALLEL DEFAULT(SHARED)
  call apply_nse_boundary(q, sim, nse, js, je, ks, ke)
  !$OMP END PARALLEL

  if (abs(q(0,js,ks,1)-q(sim%nx,js,ks,1)) < 1.0e-6_dp) then
    error stop 'x non-reflecting boundary was replaced by periodic wrapping'
  end if
  do k = ks-sim%nghost, ke+sim%nghost
    do j = js-sim%nghost, je+sim%nghost
      do i = 1-sim%nghost, sim%nx+sim%nghost
        rho = q(i,j,k,1)
        pressure = conserved_pressure(q(i,j,k,:), nse%gamma)
        if (.not. ieee_is_finite(rho) .or. rho <= 0.0_dp .or. &
            .not. ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
          write(*,'(A,4I6,2ES16.8)') 'invalid mixed-boundary state: ', &
            my_rank, i, j, k, rho, pressure
          error stop 'mixed non-reflecting boundary produced an invalid state'
        end if
      end do
    end do
  end do

  ! Exercise non-reflecting physical faces in the decomposed y and z
  ! directions as well as their edges and corners.  Internal transverse halos
  ! initially contain sentinels, so this also detects use-before-exchange.
  nse%boundary_face_type = 'non_reflecting'
  call prepare_geometry()
  nse%boundary_reference_rho = 1.0_dp
  nse%boundary_reference_p = 1.0_dp/nse%gamma
  do face = 1, nse_boundary_face_count
    nse%boundary_reference_velocity(:,face) = &
      [0.25_dp, 0.02_dp, -0.01_dp]
  end do
  call validate_boundary_scheme(sim, nse)
  q = -huge(1.0_dp)
  do k = ks, ke
    do j = js, je
      do i = 1, sim%nx
        q(i,j,k,:) = reference_state
      end do
    end do
  end do
  !$OMP PARALLEL DEFAULT(SHARED)
  call apply_nse_boundary(q, sim, nse, js, je, ks, ke)
  !$OMP END PARALLEL

  maximum_error = 0.0_dp
  do k = ks-sim%nghost, ke+sim%nghost
    do j = js-sim%nghost, je+sim%nghost
      do i = 1-sim%nghost, sim%nx+sim%nghost
        maximum_error = max(maximum_error, maxval(abs(q(i,j,k,:)-reference_state)))
      end do
    end do
  end do
  if (maximum_error > 1.0e-12_dp) then
    write(*,'(A,I0,A,ES16.8)') 'rank ', my_rank, &
      ' all-non-reflecting constant-state error = ', maximum_error
    error stop 'decomposed non-reflecting faces failed constant-state test'
  end if

  call check_characteristic_relaxation(-1._dp)
  call check_characteristic_relaxation(.37_dp)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (my_rank == root) write(*,'(A)') &
    'non-reflecting mixed boundary tests passed'
  deallocate(q)
  call MPI_Finalize(ierr)

contains

  subroutine check_characteristic_relaxation(length_override)
    real(dp), intent(in) :: length_override
    real(dp) :: inside(5),expected(5),velocity(3),distance,beta,alpha,ci,cr,cb,rhob,pb,length
    integer :: f,a,g,probe(3)
    nse%boundary_reference_velocity=0
    nse%boundary_length_scale=length_override
    nse%boundary_relaxation_strength=.7_dp
    call make_conserved(1._dp,[0._dp,0._dp,0._dp],1.1_dp/nse%gamma,nse%gamma,inside)
    do k=ks,ke
      do j=js,je
        do i=1,sim%nx
          q(i,j,k,:)=inside
        end do
      end do
    end do
    !$OMP PARALLEL DEFAULT(SHARED)
    call apply_nse_boundary(q,sim,nse,js,je,ks,ke)
    !$OMP END PARALLEL
    do f=1,6
      if(f==3.and.js/=1) cycle
      if(f==4.and.je/=sim%ny) cycle
      if(f==5.and.ks/=1) cycle
      if(f==6.and.ke/=sim%nz) cycle
      a=(f+1)/2
      select case(a)
      case(1);beta=2._dp
      case(2);beta=1._dp
      case(3);beta=1.5_dp
      end select
      do g=1,3
        distance=real(g,dp)/12
        if(test_mode=='stretched') then
          ! Analytic sinh edge coordinates, independent of boundary helper.
          distance=.5_dp*(edge_at(1,beta)+edge_at(g-1,beta)+edge_at(g,beta))
        end if
        length=1
        if(length_override>0) length=length_override
        alpha=1-exp(-.7_dp*distance/length)
        ci=sqrt(1.1_dp);cr=1
        cb=ci+.5_dp*alpha*(cr-ci)
        rhob=(cb*cb/1.1_dp)**(1/(nse%gamma-1))
        pb=(1.1_dp/nse%gamma)*rhob**nse%gamma
        velocity=0
        velocity(a)=alpha*(ci-cr)/(nse%gamma-1)
        if(mod(f,2)==1) velocity(a)=-velocity(a)
        call make_conserved(rhob,velocity,pb,nse%gamma,expected)
        probe=[1,js,ks]
        probe(a)=1-g
        if(mod(f,2)==0) probe(a)=12+g
        if(maxval(abs(q(probe(1),probe(2),probe(3),:)-expected))>2.e-12_dp) &
          error stop 'Characteristic relaxation differs from analytic ghost-distance result'
      end do
    end do
  end subroutine

  pure real(dp) function edge_at(index,beta) result(value)
    integer, intent(in) :: index
    real(dp), intent(in) :: beta
    value=.5_dp*(1+sinh(beta*(2*real(index,dp)/12-1))/sinh(beta))
  end function

  subroutine prepare_geometry()
    if(test_mode/='stretched') return
    sim%grid_mapping='sinh';nse%viscous_scheme='fv2'
    call build_sinh_axis(sim%nx,3,sim%x_min,sim%x_max,2._dp, &
      nse%boundary_face_type(1)=='periodic',axis_x)
    call build_sinh_axis(sim%ny,3,sim%y_min,sim%y_max,1._dp, &
      nse%boundary_face_type(3)=='periodic',axis_y)
    call build_sinh_axis(sim%nz,3,sim%z_min,sim%z_max,1.5_dp, &
      nse%boundary_face_type(5)=='periodic',axis_z)
    sim%dx=axis_x%minimum_width;sim%dy=axis_y%minimum_width;sim%dz=axis_z%minimum_width
  end subroutine

  subroutine make_conserved(rho_value, velocity, pressure_value, gamma, state)
    real(dp), intent(in) :: rho_value, velocity(3), pressure_value, gamma
    real(dp), intent(out) :: state(5)

    state(1) = rho_value
    state(2:4) = rho_value*velocity
    state(5) = pressure_value/(gamma-1.0_dp) + &
      0.5_dp*rho_value*sum(velocity*velocity)
  end subroutine make_conserved

  pure real(dp) function conserved_pressure(state, gamma) result(value)
    real(dp), intent(in) :: state(5), gamma
    real(dp) :: velocity(3)

    velocity = state(2:4)/state(1)
    value = (gamma-1.0_dp) * &
      (state(5)-0.5_dp*state(1)*sum(velocity*velocity))
  end function conserved_pressure

end program test_boundary_non_reflecting
