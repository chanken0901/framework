from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
FRAMEWORK_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from case_input import (  # noqa: E402
    CaseInputError,
    _validate_solver_selection,
    derive_nse_hit_transport,
    render_nse,
)
from yaml_support import load_yaml  # noqa: E402


MANIFEST = (
    FRAMEWORK_ROOT / "SolverLibrary" / "GPE" / "gp3d" / "solver_manifest.yaml"
)
NSE_MANIFEST = FRAMEWORK_ROOT / "SolverLibrary" / "NSE" / "solver_manifest.yaml"


class CaseInputProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = load_yaml(MANIFEST)

    def test_rejects_profile_that_differs_from_case(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "use_mpi": True,
                "use_cuda": False,
                "mpi_processes": 8,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "does not match"):
            _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_rejects_multiple_mpi_processes_for_single_gpu(self) -> None:
        case = {
            "solver": {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_cuda": True,
                "mpi_processes": 8,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "must be 1"):
            _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_accepts_matching_cpu_mpi_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "use_mpi": True,
                "use_cuda": False,
                "mpi_processes": 8,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")

    def test_profile_derived_fields_may_be_omitted(self) -> None:
        case = {
            "solver": {
                "use_openmp": True,
                "mpi_processes": 8,
                "omp_threads": 2,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")

    def test_accepts_runtime_openmp_for_hybrid_capable_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
                "mpi_processes": 4,
                "omp_threads": 4,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")

    def test_rejects_runtime_openmp_for_non_openmp_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": True,
                "use_cuda": True,
                "mpi_processes": 1,
                "omp_threads": 2,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "without OpenMP"):
            _validate_solver_selection(case, self.manifest, "cuda_single")


class NseCaseInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = load_yaml(NSE_MANIFEST)

    @staticmethod
    def case(mpi_processes: int = 4) -> dict:
        return {
            "case_id": "case0001",
            "physics": {
                "model": "nse",
                "nse": {
                    "nv": 5,
                    "gamma": 1.4,
                    "rho0": 1.0,
                    "mach_number": 0.5,
                    "reynolds_number": 0.0,
                    "prandtl_number": 0.72,
                },
            },
            "flow": {"type": "tgv"},
            "grid": {
                "nx": 8,
                "ny": 8,
                "nz": 8,
                "nghost": 3,
                "x_min": 0.0,
                "x_max": 6.283185307179586,
                "y_min": 0.0,
                "y_max": 6.283185307179586,
                "z_min": 0.0,
                "z_max": 6.283185307179586,
            },
            "time": {
                "dt": 1.0e-4,
                "t_max": 1.0e-4,
                "nsteps": 1,
                "cfl": 0.5,
                "use_fixed_dt": False,
                "output_frequency": 1,
            },
            "solver": {
                "profile": "cpu_mpi",
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
                "mpi_processes": mpi_processes,
                "omp_threads": 2,
            },
            "numerics": {
                "convective_scheme": "KEEP6",
                "viscous_scheme": "none",
                "boundary_condition": "periodic",
                "time_integration": "RK3",
            },
            "output": {},
        }

    def test_renders_current_modular_nse_input(self) -> None:
        text = render_nse(self.case(), self.manifest, "cpu_mpi")

        self.assertIn('initial_condition = "taylor_green"', text)
        self.assertIn('backend = "cpu_mpi"', text)
        self.assertIn('convective_scheme = "keep6"', text)
        self.assertNotIn("hybrid_smooth_scheme", text)
        self.assertNotIn("hybrid_shock_scheme", text)
        self.assertNotIn("convective_order", text)
        self.assertIn('viscous_scheme = "none"', text)
        self.assertIn('boundary_condition = "periodic"', text)
        self.assertIn('time_integrator = "ssprk3"', text)

    def test_renders_second_order_keep_selection(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "KEEP2"

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "keep2"', text)

    def test_renders_imported_turbulence_embed(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "mode": "embed",
                "x_start": 2.0,
                "blend_cells": 4,
                "velocity_offset": [0.5, 0.0, 0.0],
                "background": {
                    "density": 1.0,
                    "velocity": [0.5, 0.0, 0.0],
                    "pressure": 0.75,
                },
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('initial_condition = "imported_turbulence"', text)
        self.assertIn(
            'imported_turbulence_file = "initial_data/turbulence.slf"', text
        )
        self.assertIn('imported_turbulence_mode = "embed"', text)
        self.assertIn("imported_turbulence_x_start = 2", text)
        self.assertIn("imported_turbulence_blend_cells = 4", text)
        self.assertIn("imported_turbulence_velocity_offset_x = 0.5", text)
        self.assertIn("imported_turbulence_background_u = 0.5", text)
        self.assertIn("imported_turbulence_background_p = 0.75", text)

    def test_resolves_imported_turbulence_file_from_case_directory(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "blend_cells": 0,
            },
        }
        runtime_root = Path.cwd() / "portable_runtime"
        case_dir = runtime_root / "cases" / "case0001"

        text = render_nse(
            case,
            self.manifest,
            "cpu_mpi",
            case_dir=case_dir,
            runtime_root=runtime_root,
        )

        self.assertIn(
            'imported_turbulence_file = '
            '"cases/case0001/initial_data/turbulence.slf"',
            text,
        )

    def test_rejects_tile_with_blending(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "turbulence_tile",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "blend_cells": 2,
            },
        }

        with self.assertRaisesRegex(CaseInputError, "requires blend_cells: 0"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_invalid_imported_turbulence_velocity(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "velocity_offset": [0.5, 0.0],
            },
        }

        with self.assertRaisesRegex(CaseInputError, "exactly three numbers"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_weno5z_roe_selection_for_cpu(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "WENO5Z_ROE"

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "weno5z_roe"', text)

    def test_renders_hybrid_keep_weno_selection(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "HYBRID"
        case["numerics"]["hybrid"] = {
            "smooth_scheme": "KEEP6",
            "shock_scheme": "WENO5Z_ROE",
            "sensor": "DUCROS_PRESSURE",
            "sensor_onset": 0.02,
            "sensor_full": 0.15,
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "hybrid"', text)
        self.assertIn('hybrid_smooth_scheme = "keep6"', text)
        self.assertIn('hybrid_shock_scheme = "weno5z_roe"', text)
        self.assertIn('hybrid_sensor = "ducros_pressure"', text)
        self.assertIn("hybrid_sensor_onset = 0.02", text)
        self.assertIn("hybrid_sensor_full = 0.14999999999999999", text)

    def test_rejects_invalid_hybrid_sensor_thresholds(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "hybrid"
        case["numerics"]["hybrid"] = {
            "sensor_onset": 0.10,
            "sensor_full": 0.05,
        }

        with self.assertRaisesRegex(CaseInputError, "sensor_onset < sensor_full"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_keep_without_order_suffix(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "KEEP"

        with self.assertRaisesRegex(CaseInputError, "KEEP2 or KEEP6"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_removed_convective_order_selector(self) -> None:
        case = self.case()
        case["numerics"]["convective_order"] = 2

        with self.assertRaisesRegex(CaseInputError, "no longer supported"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_legacy_flux_selector(self) -> None:
        case = self.case()
        case["numerics"]["flux"] = "KEEP6"

        with self.assertRaisesRegex(
            CaseInputError, r"numerics\.flux.*numerics\.convective_scheme"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_unused_reconstruction_selector(self) -> None:
        case = self.case()
        case["numerics"]["reconstruction"] = "hybrid"

        with self.assertRaisesRegex(
            CaseInputError,
            r"numerics\.reconstruction.*numerics\.convective_scheme",
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_accepts_nse_mpi_openmp_profile(self) -> None:
        _validate_solver_selection(self.case(), self.manifest, "cpu_mpi")

    def test_renders_distributed_spectral_hit_input(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "spectrum": "johnsen",
                "random_seed": 24680,
                "rms_velocity": 0.1,
                "peak_wavenumber": 4.0,
                "integral_length": 1.0,
                "kolmogorov_length": 0.02,
                "dealias_fraction": 2.0 / 3.0,
                "isotropy_mode": "projected_shell",
                "isotropy_k_cutoff": 2.5,
                "isotropy_tolerance": 1.0e-8,
                "isotropy_max_iterations": 80,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('initial_condition = "hit_spectral"', text)
        self.assertIn('hit_spectrum = "johnsen"', text)
        self.assertIn("hit_seed = 24680", text)
        self.assertIn("hit_rms_velocity = 0.1", text)
        self.assertIn("hit_peak_wavenumber = 4", text)
        self.assertIn('hit_isotropy_mode = "projected_shell"', text)
        self.assertIn("hit_isotropy_k_cutoff = 2.5", text)
        self.assertIn("hit_isotropy_tolerance = 1e-08", text)
        self.assertIn("hit_isotropy_max_iterations = 80", text)
        _validate_solver_selection(
            case, self.manifest, "cpu_mpi_2decomp_fftw"
        )

    def test_rejects_hit_with_profile_without_initial_fft_backend(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.1,
                "turbulent_reynolds_number": 30.0,
                "random_seed": 13579,
                "spectrum": {
                    "type": "johnsen",
                    "johnsen": {
                        "characteristic_length": 1.0,
                        "length_scale_ratio": 2.0,
                    },
                },
            },
        }

        with self.assertRaisesRegex(
            CaseInputError, "initial-condition FFT backend"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_target_driven_johnsen_hit_input(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.3,
                "turbulent_reynolds_number": 40.0,
                "random_seed": 24680,
                "spectrum": {
                    "type": "johnsen",
                    "johnsen": {
                        "characteristic_length": 2.0,
                        "length_scale_ratio": 2.5,
                    },
                    "pope": {"integral_length": 1.0},
                },
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")
        derived = derive_nse_hit_transport(0.3, 40.0, 2.0)

        self.assertIn('hit_spectrum = "johnsen"', text)
        self.assertIn("hit_turbulent_mach = 0.29999999999999999", text)
        self.assertIn("hit_turbulent_reynolds = 40", text)
        self.assertIn("hit_johnsen_length_scale_ratio = 2.5", text)
        self.assertIn("hit_peak_wavenumber = 2.5", text)
        self.assertIn(
            f"reynolds = {format(derived['solver_reynolds'], '.17g')}", text
        )

    def test_renders_target_driven_pope_hit_input(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.5,
                "taylor_reynolds_number": 30.0,
                "spectrum": {
                    "type": "pope",
                    "johnsen": {
                        "characteristic_length": 1.0,
                        "length_scale_ratio": 2.0,
                    },
                    "pope": {
                        "integral_length": 1.0,
                        "energy_constant": 1.5,
                        "large_scale_constant": 6.78,
                        "dissipation_constant": 0.4,
                        "large_scale_exponent": 2.0,
                        "dissipation_exponent": 5.2,
                    },
                },
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")
        derived = derive_nse_hit_transport(0.5, 30.0, 1.0)

        self.assertIn('hit_spectrum = "pope"', text)
        self.assertIn("hit_pope_energy_constant = 1.5", text)
        self.assertIn("hit_pope_large_scale_constant = 6.7800000000000002", text)
        self.assertIn("hit_pope_dissipation_exponent = 5.2000000000000002", text)
        self.assertIn(
            "hit_kolmogorov_length = "
            f"{format(derived['kolmogorov_length'], '.17g')}",
            text,
        )

    def test_rejects_unknown_hit_spectrum_type(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.5,
                "turbulent_reynolds_number": 30.0,
                "spectrum": {"type": "poppe"},
            },
        }

        with self.assertRaisesRegex(CaseInputError, "spectrum.type"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_rejects_hit_target_missing_reynolds_number(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.5,
                "spectrum": {
                    "type": "johnsen",
                    "johnsen": {
                        "characteristic_length": 1.0,
                        "length_scale_ratio": 2.0,
                    },
                },
            },
        }

        with self.assertRaisesRegex(CaseInputError, "turbulent_reynolds_number"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_renders_petersen_livescu_forcing(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "spectrum": "low_wavenumber",
                "fft_backend": "2decomp_fftw",
                "k_cutoff": 2.5,
                "target_dissipation": 0.1,
                "dilatational_ratio": 0.25,
                "denominator_floor": 1.0e-14,
                "max_coefficient": 20.0,
                "report_interval": 50,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)
        self.assertIn('forcing_spectrum = "low_wavenumber"', text)
        self.assertIn('forcing_fft_backend = "2decomp_fftw"', text)
        self.assertIn("forcing_k_cutoff = 2.5", text)
        self.assertIn("forcing_target_dissipation = 0.1", text)
        self.assertIn("forcing_dilatational_ratio = 0.25", text)
        self.assertIn("forcing_report_interval = 50", text)

    def test_normalizes_human_readable_forcing_selectors(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "Petersen-Livescu",
            "petersen_livescu": {
                "spectrum": "Low Wavenumber",
                "fft_backend": "2decomp-fftw",
                "target_dissipation": 0.1,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)
        self.assertIn('forcing_spectrum = "low_wavenumber"', text)
        self.assertIn('forcing_fft_backend = "2decomp_fftw"', text)

    def test_normalizes_full_wavenumber_forcing_alias(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "spectrum": "full_wavenumber",
                "target_dissipation": 0.1,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_spectrum = "full_spectrum"', text)

    def test_rejects_unknown_forcing_spectrum(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "spectrum": "full_wavenumbar",
                "target_dissipation": 0.1,
            },
        }

        with self.assertRaisesRegex(
            CaseInputError, "forcing.petersen_livescu.spectrum"
        ):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_rejects_forcing_with_profile_without_fft_backend(self) -> None:
        case = self.case()
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {"target_dissipation": 0.1},
        }

        with self.assertRaisesRegex(
            CaseInputError, "compatible staged profile"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_accepts_legacy_forcing_scheme_alias(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "scheme": "petersen_livescu",
            "target_dissipation": 0.1,
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)

    def test_ignores_inactive_petersen_livescu_settings(self) -> None:
        case = self.case()
        case["forcing"] = {
            "type": "none",
            "petersen_livescu": {
                "target_dissipation": 0.1,
                "fft_backend": "2decomp_fftw",
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('forcing_scheme = "none"', text)
        self.assertNotIn("forcing_target_dissipation", text)

    def test_rejects_unknown_forcing_type(self) -> None:
        case = self.case()
        case["forcing"] = {"type": "petersen_livecu"}

        with self.assertRaisesRegex(CaseInputError, "unknown forcing.type"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_misspelled_type_specific_forcing_key(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "target_disipation": 0.1,
            },
        }

        with self.assertRaisesRegex(
            CaseInputError, "unknown forcing.petersen_livescu key"
        ):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_rejects_mixed_nested_and_legacy_forcing_settings(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "target_dissipation": 0.1,
            "petersen_livescu": {"target_dissipation": 0.2},
        }

        with self.assertRaisesRegex(CaseInputError, "mixes the new"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_renders_nse_single_gpu_input(self) -> None:
        case = self.case(mpi_processes=1)
        case["solver"].update(
            {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
                "cuda_device": 2,
            }
        )
        text = render_nse(case, self.manifest, "cuda_single")

        self.assertIn('backend = "cuda"', text)
        self.assertIn("use_mpi = .false.", text)
        self.assertIn("use_openmp = .false.", text)
        self.assertIn("cuda_device = 2", text)
        _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_renders_weno5z_roe_for_cuda_profile(self) -> None:
        case = self.case(mpi_processes=1)
        case["solver"].update(
            {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
            }
        )
        case["numerics"]["convective_scheme"] = "WENO5Z_ROE"

        text = render_nse(case, self.manifest, "cuda_single")

        self.assertIn('convective_scheme = "weno5z_roe"', text)

    def test_renders_nse_multi_gpu_input(self) -> None:
        case = self.case(mpi_processes=2)
        case["solver"].update(
            {
                "profile": "cuda_mpi",
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
            }
        )

        text = render_nse(case, self.manifest, "cuda_mpi")

        self.assertIn('backend = "cuda_mpi"', text)
        self.assertIn("use_mpi = .true.", text)
        self.assertIn("use_openmp = .false.", text)
        _validate_solver_selection(case, self.manifest, "cuda_mpi")

    def test_rejects_too_few_nse_mpi_processes(self) -> None:
        with self.assertRaisesRegex(CaseInputError, "at least 2"):
            _validate_solver_selection(
                self.case(mpi_processes=1), self.manifest, "cpu_mpi"
            )


if __name__ == "__main__":
    unittest.main()
