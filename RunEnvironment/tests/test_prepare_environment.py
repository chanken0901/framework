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
    _slurm_script,
)


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
        )

        self.assertIn("python3 tools/run_case.py --run", script)
        self.assertNotIn("--build", script)
        self.assertNotIn("--all", script)


if __name__ == "__main__":
    unittest.main()
