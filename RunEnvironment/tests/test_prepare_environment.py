from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from prepare_environment import (  # noqa: E402
    EnvironmentError,
    _global_case_index_path,
    _numbered_case_destination,
    _parallel_features,
    _requested_case_destination,
    _slurm_script,
)
from yaml_support import load_yaml  # noqa: E402


class NumberedCaseDestinationTests(unittest.TestCase):
    def test_uses_first_available_model_case_number(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "gpe_case0001").mkdir()
            (root / "gpe_case0002.tar.gz").touch()

            output, case_id = _numbered_case_destination(root, "gpe")

            self.assertEqual(output, root / "gpe_case0003")
            self.assertEqual(case_id, "case0003")

    def test_supports_an_independent_model_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "gpe_case0001").mkdir()

            output, case_id = _numbered_case_destination(root, "nse")

            self.assertEqual(output, root / "nse_case0001")
            self.assertEqual(case_id, "case0001")

    def test_selects_an_existing_case_number_for_overwrite(self) -> None:
        root = Path("C:/ResearchRuns")

        output, case_id = _requested_case_destination(root, "nse", "case0015")

        self.assertEqual(output, root / "nse_case0015")
        self.assertEqual(case_id, "case0015")

    def test_normalizes_a_numeric_overwrite_case_id(self) -> None:
        root = Path("C:/ResearchRuns")

        output, case_id = _requested_case_destination(root, "gpe", "15")

        self.assertEqual(output, root / "gpe_case0015")
        self.assertEqual(case_id, "case0015")

    def test_rejects_an_invalid_overwrite_case_id(self) -> None:
        with self.assertRaisesRegex(EnvironmentError, "--case-id"):
            _requested_case_destination(Path("C:/ResearchRuns"), "nse", "latest")

    def test_global_index_must_be_outside_case_environment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "gpe_case0001"

            with self.assertRaises(EnvironmentError):
                _global_case_index_path(
                    {"case_index": str(output / "cases.csv")},
                    root / "environment.yaml",
                    output,
                )


class SlurmScriptTests(unittest.TestCase):
    def test_generated_job_runs_existing_executable_only(self) -> None:
        script = _slurm_script(
            {
                "job_name": "gpe-case",
                "nodes": 1,
                "tasks_per_node": 4,
                "modules": ["gcc", "openmpi"],
            },
            include_tests=False,
            use_mpi=True,
            use_openmp=True,
        )

        self.assertIn("python3 tools/run_case.py --run", script)
        self.assertIn('--processes "${SLURM_NTASKS}"', script)
        self.assertIn('--omp-threads "${SLURM_CPUS_PER_TASK}"', script)
        self.assertNotIn("#SBATCH --ntasks-per-node", script)
        self.assertNotIn("#SBATCH --cpus-per-task", script)
        self.assertNotIn("--build", script)
        self.assertNotIn("--all", script)


class ParallelFeatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        framework_root = SCRIPT_DIR.parents[1]
        cls.gpe_manifest = load_yaml(
            framework_root
            / "SolverLibrary"
            / "GPE"
            / "gp3d"
            / "solver_manifest.yaml"
        )

    def test_accepts_mpi_openmp_capability_without_counts(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
            }
        }

        self.assertEqual(
            _parallel_features(design, self.gpe_manifest, "cpu_mpi_fftw"),
            (True, True, False),
        )

    def test_rejects_cuda_flag_for_cpu_profile(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": True,
            }
        }

        with self.assertRaisesRegex(EnvironmentError, "use_cuda"):
            _parallel_features(design, self.gpe_manifest, "cpu_mpi_fftw")


class NseCaseTemplateTests(unittest.TestCase):
    def test_parallel_features_and_runtime_counts_are_template_fields(self) -> None:
        template = (
            SCRIPT_DIR / "case_templates" / "nse.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("use_openmp: {{use_openmp}}", template)
        self.assertIn("mpi_processes: {{mpi_processes}}", template)
        self.assertIn("omp_threads: {{omp_threads}}", template)
        self.assertNotIn("use_openmp: true", template)

    def test_flow_conditions_share_one_nse_template(self) -> None:
        template = (
            SCRIPT_DIR / "case_templates" / "nse.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("type: taylor_green", template)
        self.assertIn("taylor_green:", template)
        self.assertIn("hit:", template)

    def test_forcing_conditions_use_type_specific_mapping(self) -> None:
        template = (
            SCRIPT_DIR / "case_templates" / "nse.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("forcing:\n  # choices: none, petersen_livescu", template)
        self.assertIn("  type: none\n  petersen_livescu:\n", template)
        self.assertIn(
            "# choices: full_spectrum, low_wavenumber", template
        )
        self.assertIn("# choices: auto, 2decomp_fftw, cufft", template)
        self.assertIn("    target_dissipation: 0.1", template)
        self.assertNotIn("\n  target_dissipation: 0.1", template)


if __name__ == "__main__":
    unittest.main()
