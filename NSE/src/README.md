# mod_slf_output.f90 update

This version follows the rule:

- NSE output: conservative variables only, e.g. `q(:,:,:,1:5) = [rho, rho_u, rho_v, rho_w, rho_E]`.
- GPE output: wave function only, stored as `psi_real` and `psi_imag`.

Derived quantities such as velocity, pressure, density, and phase should be computed in post-processing.

## NSE usage

```fortran
use mod_slf_output, only : write_nse_conserved_slf

call write_nse_conserved_slf(sim, step, time, q, rank=myrank)
```

## GPE usage

```fortran
use mod_slf_output, only : write_gpe_psi_slf

call write_gpe_psi_slf(sim, step, time, psi, rank=myrank)
```

The generic interface is also available:

```fortran
call write_field_slf(sim, step, time, q, variable_names, rank=myrank)
call write_field_slf(sim, step, time, psi, rank=myrank)
```
