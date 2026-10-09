from __future__ import annotations

import argparse
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import shutil
import struct
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import postprocess_case  # noqa: E402


class PostprocessStatisticsTests(unittest.TestCase):
    def test_fft_entrypoint_executes_and_dry_run_is_read_only(self):
        import numpy as np
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'cases/case0001/output'
            output.mkdir(parents=True)
            tool_dir = root / 'SolverLibrary/NSE/tools'
            tool_dir.mkdir(parents=True)
            source = SCRIPT_DIR.parents[1] / 'SolverLibrary/NSE/tools'
            for name in ('slf_fft.py', 'slf_to_paraview_merged_cropghost.py'):
                shutil.copy2(source / name, tool_dir / name)
            context = dict(root=root, model='nse', case_root=output.parent,
                           input_dir=output, meta=output/'meta.json', gamma=1.67, reynolds=None)
            with patch.object(postprocess_case, '_case_context', return_value=context):
                self.assertEqual(postprocess_case.main(['--task', 'fft', '--dry-run']), 0)
                self.assertFalse((output.parent/'fft').exists())
                data = np.ones((4,4,4,1), dtype='<f8')
                with (output/'field_000001.slf').open('wb') as f:
                    f.write(b'SLF1\0\0\0\0')
                    f.write(struct.pack('<iii', 1,2,3))
                    f.write(np.array(data.shape, dtype='<i4').tobytes())
                    f.write(np.array([1,0,4,4,4,0,1,0], dtype='<i4').tobytes())
                    f.write(struct.pack('<d6di', .1,0,4,0,4,0,4,1))
                    f.write(b'rho'.ljust(32,b'\0'))
                    f.write(data.tobytes(order='F'))
                context['meta'].write_text(json.dumps(dict(grid=[4,4,4])))
                self.assertEqual(postprocess_case.main(['--task','fft','--save-fft']), 0)
                self.assertTrue((output.parent/'fft/fft_000001_rho.npz').is_file())
                context['model'] = 'gpe'
                with self.assertRaises(SystemExit):
                    postprocess_case.main(['--task','fft','--dry-run'])

    def test_all_dispatches_three_tasks_and_fft_options(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tool = root/'SolverLibrary/NSE/tools/slf_fft.py'
            tool.parent.mkdir(parents=True)
            tool.touch()
            context = dict(root=root, model='nse', case_root=root/'cases/case1',
                           input_dir=root/'input', meta=root/'input/meta.json', gamma=1.67)
            args = argparse.Namespace(task='all', fft_output='custom', fft_field='p',
                fft_window='hann', fft_keep_mean=True, save_fft=True, fft_overwrite=True,
                steps='10:20:10', layout='rank')
            with patch.object(postprocess_case, '_case_context', return_value=context), \
                 patch.object(postprocess_case, 'build_command', return_value=(root, ['pv'])), \
                 patch.object(postprocess_case, 'build_statistics_command', return_value=(root, ['stats'])):
                commands = postprocess_case.build_commands(args)
            self.assertEqual([x[0] for x in commands], ['ParaView','Turbulence statistics','FFT'])
            command = commands[-1][2]
            for option, value in [('--field','p'),('--gamma','1.67'),('--window','hann'),('--steps','10:20:10')]:
                self.assertEqual(command[command.index(option)+1], value)
            self.assertIn('--save-fft', command)
            self.assertIn('--keep-mean', command)
            self.assertIn('--overwrite', command)
            self.assertEqual(commands[-1][1], root/'custom')

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
