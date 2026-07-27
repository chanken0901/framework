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
)
from yaml_support import load_yaml  # noqa: E402


MANIFEST = (
    FRAMEWORK_ROOT / "SolverLibrary" / "GPE" / "gp3d" / "solver_manifest.yaml"
)


class CaseInputProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = load_yaml(MANIFEST)

    def test_rejects_profile_that_differs_from_case(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "processes": 8,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "does not match"):
            _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_rejects_multiple_processes_for_single_gpu(self) -> None:
        case = {
            "solver": {
                "profile": "cuda_single",
                "processes": 8,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "must be 1"):
            _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_accepts_matching_cpu_mpi_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "processes": 8,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")


if __name__ == "__main__":
    unittest.main()
