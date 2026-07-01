program use_in_gpe
  use mod_precision, only : dp
  use mod_common_config, only : simulation_config, init_simulation_config, print_simulation_config
  use mod_model_config, only : gpe_config, init_gpe_config, print_gpe_config
  use mod_input_reader, only : read_all_inputs
  use mod_slf_output, only : write_meta_json, write_field_slf
  implicit none

  type(simulation_config) :: sim
  type(gpe_config) :: gpe
  complex(dp), allocatable :: psi(:,:,:)

  call init_simulation_config(sim)
  call init_gpe_config(gpe)
  call read_all_inputs('input_gpe.dat', sim, gpe=gpe)

  call print_simulation_config(sim)
  call print_gpe_config(gpe)

  allocate(psi(sim%nx, sim%ny, sim%nz))
  psi = (1.0_dp, 0.0_dp)

  call write_meta_json(sim, [character(len=32) :: 'psi_real','psi_imag','rho','phase'])
  call write_field_slf(sim, step=0, time=0.0_dp, psi=psi)
end program use_in_gpe
