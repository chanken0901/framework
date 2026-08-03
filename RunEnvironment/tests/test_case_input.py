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
                "flux": "KEEP6",
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
        self.assertNotIn("convective_order", text)
        self.assertIn('viscous_scheme = "none"', text)
        self.assertIn('boundary_condition = "periodic"', text)
        self.assertIn('time_integrator = "ssprk3"', text)

    def test_renders_second_order_keep_selection(self) -> None:
        case = self.case()
        case["numerics"]["flux"] = "KEEP2"

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "keep2"', text)

    def test_rejects_keep_without_order_suffix(self) -> None:
        case = self.case()
        case["numerics"]["flux"] = "KEEP"

        with self.assertRaisesRegex(CaseInputError, "KEEP2 or KEEP6"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_removed_convective_order_selector(self) -> None:
        case = self.case()
        case["numerics"]["convective_order"] = 2

        with self.assertRaisesRegex(CaseInputError, "no longer supported"):
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
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('initial_condition = "hit_spectral"', text)
        self.assertIn('hit_spectrum = "johnsen"', text)
        self.assertIn("hit_seed = 24680", text)
        self.assertIn("hit_rms_velocity = 0.1", text)
        self.assertIn("hit_peak_wavenumber = 4", text)
        _validate_solver_selection(
            case, self.manifest, "cpu_mpi_2decomp_fftw"
        )

    def test_renders_petersen_livescu_forcing(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "spectrum": "low_wavenumber",
            "fft_backend": "2decomp_fftw",
            "k_cutoff": 2.5,
            "target_dissipation": 0.1,
            "dilatational_ratio": 0.25,
            "denominator_floor": 1.0e-14,
            "max_coefficient": 20.0,
            "report_interval": 50,
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
            "spectrum": "Low Wavenumber",
            "fft_backend": "2decomp-fftw",
            "target_dissipation": 0.1,
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)
        self.assertIn('forcing_spectrum = "low_wavenumber"', text)
        self.assertIn('forcing_fft_backend = "2decomp_fftw"', text)

    def test_rejects_forcing_with_profile_without_fft_backend(self) -> None:
        case = self.case()
        case["forcing"] = {
            "type": "petersen_livescu",
            "target_dissipation": 0.1,
        }

        with self.assertRaisesRegex(
            CaseInputError, "cpu_mpi_2decomp_fftw"
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

    def test_rejects_too_few_nse_mpi_processes(self) -> None:
        with self.assertRaisesRegex(CaseInputError, "at least 4"):
            _validate_solver_selection(
                self.case(mpi_processes=2), self.manifest, "cpu_mpi"
            )


if __name__ == "__main__":
    unittest.main()
