from __future__ import annotations

import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import postprocess_case  # noqa: E402


class PostprocessStatisticsTests(unittest.TestCase):
    def test_hit_statistics_uses_derived_solver_reynolds(self) -> None:
        case = {
            "physics": {"nse": {"reynolds_number": 100.0}},
            "flow": {
                "type": "hit",
                "hit": {
                    "turbulent_mach_number": 0.5,
                    "turbulent_reynolds_number": 30.0,
                    "spectrum": {
                        "type": "pope",
                        "pope": {"integral_length": 1.0},
                    },
                },
            },
        }

        actual = postprocess_case._resolved_nse_reynolds(case)
        expected = 135.0 / ((1.5**0.5) * (0.5 / (3.0**0.5)))
        self.assertAlmostEqual(actual, expected, places=12)

    def test_statistics_command_uses_case_physics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tools = root / "tools"
            case_root = root / "cases" / "case0012"
            solver_tools = root / "SolverLibrary" / "NSE" / "tools"
            for path in (tools, case_root / "output", solver_tools):
                path.mkdir(parents=True, exist_ok=True)
            (root / "environment.lock.json").write_text(
                json.dumps(
                    {
                        "model": "nse",
                        "case_directory": "cases/case0012",
                    }
                ),
                encoding="utf-8",
            )
            (case_root / "case.yaml").write_text(
                "\n".join(
                    [
                        "physics:",
                        "  nse:",
                        "    gamma: 1.67",
                        "    reynolds_number: 250.0",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            statistics_tool = solver_tools / "nse_turbulence_statistics.py"
            statistics_tool.write_text("", encoding="ascii")
            original_file = postprocess_case.__file__
            postprocess_case.__file__ = str(tools / "postprocess_case.py")
            try:
                args = argparse.Namespace(
                    case_directory=None,
                    input_dir=None,
                    meta=None,
                    statistics_output=None,
                    gamma=None,
                    reynolds=None,
                    layout="auto",
                    steps=None,
                    density_floor=1.0e-12,
                    pressure_floor=1.0e-12,
                )
                output, command = postprocess_case.build_statistics_command(args)
            finally:
                postprocess_case.__file__ = original_file

            self.assertEqual(
                output,
                (case_root / "statistics" / "turbulence_statistics.csv").resolve(),
            )
            self.assertIn(str(statistics_tool), command)
            self.assertEqual(command[command.index("--steps") + 1], "all")
            self.assertEqual(command[command.index("--gamma") + 1], "1.67")
            self.assertEqual(command[command.index("--reynolds") + 1], "250.0")

    def test_paraview_command_passes_inspection_mode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tools = root / "tools"
            case_root = root / "cases" / "case0001"
            solver_tools = root / "SolverLibrary" / "GPE" / "gp3d" / "tools"
            for path in (tools, case_root / "output", solver_tools):
                path.mkdir(parents=True, exist_ok=True)
            (root / "environment.lock.json").write_text(
                json.dumps({"model": "gpe", "case_directory": "cases/case0001"}),
                encoding="utf-8",
            )
            (case_root / "case.yaml").write_text("physics: {}\n", encoding="ascii")
            converter = solver_tools / "slf_to_paraview_merged_cropghost.py"
            converter.write_text("", encoding="ascii")
            original_file = postprocess_case.__file__
            postprocess_case.__file__ = str(tools / "postprocess_case.py")
            try:
                args = argparse.Namespace(
                    case_directory=None,
                    input_dir=None,
                    output_dir=None,
                    meta=None,
                    gamma=None,
                    reynolds=None,
                    fields=None,
                    derive="auto",
                    layout="auto",
                    stride=2,
                    steps=None,
                    pvd_name="collection.pvd",
                    inspect_only=True,
                )
                _output, command = postprocess_case.build_command(args)
            finally:
                postprocess_case.__file__ = original_file

            self.assertIn("--inspect-only", command)
            self.assertEqual(command[command.index("--steps") + 1], "latest")
            self.assertEqual(command[command.index("--fields") + 1], "density,phase")


if __name__ == "__main__":
    unittest.main()
