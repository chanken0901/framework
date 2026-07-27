program test_restart_slf
  use gp3d_types, only: dp, gp3d_grid_t, gp3d_state_t
  use gp3d_grid, only: gp3d_grid_init_bounds
  use gp3d_io, only: gp3d_output_config_t, gp3d_output_config_from_grid, &
    gp3d_write_field_complex3_slf
  use gp3d_restart, only: gp3d_restart_info_t, gp3d_restart_load
  implicit none

  type(gp3d_grid_t) :: global_grid, local_grid
  type(gp3d_state_t) :: source, restored, local_restored
  type(gp3d_output_config_t) :: cfg
  type(gp3d_restart_info_t) :: info
  character(len=*), parameter :: global_dir = "restart_roundtrip_output"
  character(len=*), parameter :: rank_dir = "restart_rank_roundtrip_output"
  character(len=256) :: global_file, rank_base
  integer :: i, j, k

  call gp3d_grid_init_bounds(global_grid, 4, 3, 4, -2.0_dp, 2.0_dp, &
    -1.5_dp, 1.5_dp, -1.0_dp, 3.0_dp)
  allocate(source%psi(4, 3, 4), restored%psi(4, 3, 4))
  do k = 1, 4
    do j = 1, 3
      do i = 1, 4
        source%psi(i,j,k) = cmplx(real(i + 10 * j + 100 * k, dp), &
          -real(2 * i + 3 * j + 5 * k, dp), kind=dp)
      end do
    end do
  end do

  call gp3d_output_config_from_grid(cfg, global_grid, global_dir, "restart_roundtrip")
  cfg%mpi_nprocs = 1
  call gp3d_write_field_complex3_slf(cfg, 7, 1.25_dp, source%psi)
  global_file = global_dir // "/field_000007.slf"
  restored%psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
  call gp3d_restart_load(global_file, restored, global_grid, info)
  call assert_true(info%loaded, "global restart was not marked loaded")
  call assert_true(info%step == 7, "global restart step was not preserved")
  call assert_true(abs(info%time - 1.25_dp) < 1.0e-14_dp, "global restart time was not preserved")
  call assert_close(restored%psi, source%psi, "global restart values differ")

  call gp3d_output_config_from_grid(cfg, global_grid, rank_dir, "restart_rank_roundtrip")
  cfg%use_mpi = .true.
  cfg%mpi_nprocs = 2
  cfg%local_nz = 2
  cfg%k_start = 1
  cfg%k_end = 2
  call gp3d_write_field_complex3_slf(cfg, 9, 2.5_dp, source%psi(:,:,1:2), rank=0)
  cfg%k_start = 3
  cfg%k_end = 4
  call gp3d_write_field_complex3_slf(cfg, 9, 2.5_dp, source%psi(:,:,3:4), rank=1)

  rank_base = rank_dir // "/field_000009.slf"
  call delete_if_exists(rank_base)
  restored%psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
  call gp3d_restart_load(rank_base, restored, global_grid, info)
  call assert_true(info%step == 9, "rank-family restart step was not preserved")
  call assert_true(info%source_nprocs == 2, "rank-family source process count is wrong")
  call assert_close(restored%psi, source%psi, "assembled rank-family restart values differ")

  call gp3d_grid_init_bounds(local_grid, 4, 3, 4, -2.0_dp, 2.0_dp, &
    -1.5_dp, 1.5_dp, -1.0_dp, 3.0_dp, rank=1, nprocs=2)
  allocate(local_restored%psi(4, 3, local_grid%local_nz))
  local_restored%psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
  call gp3d_restart_load(rank_base, local_restored, local_grid, info)
  call assert_close(local_restored%psi, source%psi(:,:,3:4), &
    "rank-family restart did not select the current local slab")

  write(*,'(a)') "SLF restart roundtrip passed"

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine assert_true

  subroutine assert_close(actual, expected, message)
    complex(dp), intent(in) :: actual(:,:,:), expected(:,:,:)
    character(len=*), intent(in) :: message

    if (maxval(abs(actual - expected)) > 1.0e-13_dp) error stop message
  end subroutine assert_close

  subroutine delete_if_exists(filename)
    character(len=*), intent(in) :: filename
    logical :: exists
    integer :: unit, ios

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=trim(filename), status="old", iostat=ios)
    if (ios == 0) close(unit, status="delete")
  end subroutine delete_if_exists

end program test_restart_slf
