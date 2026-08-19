from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from run_case import RunCaseError, _case_parallel_settings, _prepare_input  # noqa: E402


class RunCasePrepareTests(unittest.TestCase):
    def test_prepare_uses_lock_model_and_profile(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            tools_dir = root / "tools"
            build_dir = root / "ScriptLibrary" / "BuildSolver"
            solver_dir = root / "SolverLibrary" / "GPE" / "gp3d"
            for path in (case_dir, tools_dir, build_dir, solver_dir):
                path.mkdir(parents=True, exist_ok=True)

            (case_dir / "case.yaml").write_text(
                "physics:\n  model: gpe\n",
                encoding="utf-8",
            )
            (solver_dir / "solver_manifest.yaml").write_text(
                "schema_version: 1\n",
                encoding="utf-8",
            )
            (build_dir / "model_catalog.yaml").write_text(
                "\n".join(
                    [
                        "schema_version: 1",
                        "models:",
                        "  gpe:",
                        "    library_subpath: GPE/gp3d",
                        "    manifest: solver_manifest.yaml",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (tools_dir / "case_input.py").write_text(
                "\n".join(
                    [
                        "import argparse",
                        "from pathlib import Path",
                        "p = argparse.ArgumentParser()",
                        "p.add_argument('--case', required=True)",
                        "p.add_argument('--manifest', required=True)",
                        "p.add_argument('--model', required=True)",
                        "p.add_argument('--profile', required=True)",
                        "p.add_argument('--overwrite', action='store_true')",
                        "a = p.parse_args()",
                        "Path(a.case).with_name('input.nml').write_text(",
                        "    f'{a.model}:{a.profile}', encoding='ascii')",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            lock = {
                "case_directory": "cases/case0001",
                "input_name": "input.nml",
                "model": "gpe",
                "profile": "cuda_single",
            }

            output = _prepare_input(root, lock, force=True, dry_run=False)

            self.assertEqual(output.read_text(encoding="ascii"), "gpe:cuda_single")

    def test_current_input_is_not_regenerated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            case_path = case_dir / "case.yaml"
            input_path = case_dir / "input.nml"
            case_path.write_text("case_id: case0001\n", encoding="utf-8")
            input_path.write_text("current", encoding="ascii")
            input_path.touch()
            lock = {
                "case_directory": "cases/case0001",
                "input_name": "input.nml",
                "model": "gpe",
                "profile": "cuda_single",
            }

            output = _prepare_input(root, lock, force=False, dry_run=False)

            self.assertEqual(output.read_text(encoding="ascii"), "current")

    def test_parallel_settings_follow_case_yaml(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            (case_dir / "case.yaml").write_text(
                "\n".join(
                    [
                        "solver:",
                        "  use_mpi: true",
                        "  mpi_processes: 4",
                        "  use_openmp: true",
                        "  use_cuda: false",
                        "  omp_threads: 6",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            lock = {
                "case_directory": "cases/case0001",
                "use_mpi": True,
                "use_openmp": True,
                "openmp_capable": True,
            }

            self.assertEqual(_case_parallel_settings(root, lock), (4, 6, True))

    def test_profile_derived_case_fields_may_be_omitted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            (case_dir / "case.yaml").write_text(
                "solver:\n  use_openmp: false\n  mpi_processes: 4\n  omp_threads: 1\n",
                encoding="utf-8",
            )
            lock = {
                "case_directory": "cases/case0001",
                "use_mpi": True,
                "use_cuda": False,
                "use_openmp": True,
                "openmp_capable": True,
            }

            self.assertEqual(_case_parallel_settings(root, lock), (4, 1, False))

    def test_openmp_can_be_enabled_per_case_when_profile_is_capable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            (case_dir / "case.yaml").write_text(
                "solver:\n  use_openmp: true\n  mpi_processes: 4\n  omp_threads: 6\n",
                encoding="utf-8",
            )
            lock = {
                "case_directory": "cases/case0001",
                "use_mpi": True,
                "use_cuda": False,
                "use_openmp": False,
                "openmp_capable": True,
            }

            self.assertEqual(_case_parallel_settings(root, lock), (4, 6, True))

    def test_parallel_settings_reject_zero_threads(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            (case_dir / "case.yaml").write_text(
                "solver:\n  mpi_processes: 2\n  omp_threads: 0\n",
                encoding="utf-8",
            )
            lock = {
                "case_directory": "cases/case0001",
                "use_mpi": True,
                "use_openmp": True,
                "openmp_capable": True,
            }

            with self.assertRaisesRegex(RunCaseError, "must be positive"):
                _case_parallel_settings(root, lock)

    def test_rejects_removed_processes_key(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            (case_dir / "case.yaml").write_text(
                "solver:\n  processes: 4\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(RunCaseError, "no longer supported"):
                _case_parallel_settings(
                    root,
                    {
                        "case_directory": "cases/case0001",
                        "use_mpi": True,
                        "use_openmp": False,
                        "openmp_capable": True,
                    },
                )

    def test_parallel_settings_reject_environment_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            case_dir.mkdir(parents=True)
            (case_dir / "case.yaml").write_text(
                "\n".join(
                    [
                        "solver:",
                        "  use_mpi: false",
                        "  use_openmp: true",
                        "  use_cuda: false",
                        "  mpi_processes: 1",
                        "  omp_threads: 4",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(RunCaseError, "regenerate"):
                _case_parallel_settings(
                    root,
                    {
                        "case_directory": "cases/case0001",
                        "use_mpi": True,
                        "use_openmp": True,
                        "use_cuda": False,
                    },
                )


if __name__ == "__main__":
    unittest.main()
