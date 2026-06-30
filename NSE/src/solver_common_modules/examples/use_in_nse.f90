program use_in_nse
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config, print_simulation_config
  use mod_model_config, only : nse_config, init_nse_config, print_nse_config
  use mod_input_reader, only : read_all_inputs
  use mod_slf_output, only : write_meta_json, write_field_slf
  implicit none

  type(simulation_config) :: sim
  type(nse_config) :: nse
  real(dp), allocatable :: qvis(:,:,:,:)
  character(len=32) :: names(5)

  call init_simulation_config(sim)
  call init_nse_config(nse)
  call read_all_inputs('input_nse.dat', sim, nse=nse)

  call print_simulation_config(sim)
  call print_nse_config(nse)

  allocate(qvis(sim%nx, sim%ny, sim%nz, nse%nv))
  qvis = 0.0_dp
  names = [character(len=32) :: 'rho','u','v','w','p']

  call write_meta_json(sim, names)
  call write_field_slf(sim, step=0, time=0.0_dp, field=qvis, variable_names=names)
end program use_in_nse
