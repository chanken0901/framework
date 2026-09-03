from __future__ import annotations

import copy
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from run_case import (  # noqa: E402
    RunCaseError,
    _case_parallel_settings,
    _prepare_input,
    _select_case_profile,
    _sync_global_case_index,
)
from global_case_index import GlobalCaseIndexError  # noqa: E402
from yaml_support import load_yaml  # noqa: E402


class RunCaseProfileSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = load_yaml(
            SCRIPT_DIR.parents[1]
            / "SolverLibrary"
            / "NSE"
            / "solver_manifest.yaml"
        )
        cls.lock = {
            "model": "nse",
            "profile": "cpu_mpi",
            "available_profiles": ["cpu_mpi", "cpu_mpi_2decomp_fftw"],
            "profile_explicit": False,
        }

    def test_basic_case_keeps_lightweight_mpi_profile(self) -> None:
        case = {"flow": {"type": "taylor_green"}, "forcing": {"type": "none"}}

        profile, requirements = _select_case_profile(
            case, self.manifest, self.lock
        )

        self.assertEqual(profile, "cpu_mpi")
        self.assertEqual(requirements, set())

    def test_hit_case_selects_distributed_fft_profile(self) -> None:
        case = {"flow": {"type": "hit"}, "forcing": {"type": "none"}}

        profile, requirements = _select_case_profile(
            case, self.manifest, self.lock
        )

        self.assertEqual(profile, "cpu_mpi_2decomp_fftw")
        self.assertEqual(requirements, {"hit_spectral"})

    def test_forcing_case_selects_distributed_fft_profile(self) -> None:
        case = {
            "flow": {"type": "taylor_green"},
            "forcing": {"type": "petersen_livescu"},
        }

        profile, requirements = _select_case_profile(
            case, self.manifest, self.lock
        )

        self.assertEqual(profile, "cpu_mpi_2decomp_fftw")
        self.assertEqual(requirements, {"forcing_fft"})

    def test_hit_and_forcing_require_both_fft_capabilities(self) -> None:
        case = {
            "flow": {"type": "homogeneous-isotropic-turbulence"},
            "forcing": {"type": "petersen_livescu"},
        }

        profile, requirements = _select_case_profile(
            case, self.manifest, self.lock
        )

        self.assertEqual(profile, "cpu_mpi_2decomp_fftw")
        self.assertEqual(requirements, {"hit_spectral", "forcing_fft"})

    def test_selects_least_specialized_satisfying_profile(self) -> None:
        manifest = copy.deepcopy(self.manifest)
        hit_only = copy.deepcopy(manifest["profiles"]["cpu_mpi"])
        hit_only["capabilities"] = ["nse_basic", "hit_spectral"]
        hit_only["cmake"]["NSE_INIT_FFT_BACKEND"] = "2decomp_fftw"
        manifest["profiles"]["cpu_mpi_hit_only"] = hit_only
        lock = {
            **self.lock,
            "available_profiles": [
                "cpu_mpi",
                "cpu_mpi_2decomp_fftw",
                "cpu_mpi_hit_only",
            ],
        }

        profile, _ = _select_case_profile(
            {"flow": {"type": "hit"}, "forcing": {"type": "none"}},
            manifest,
            lock,
        )

        self.assertEqual(profile, "cpu_mpi_hit_only")

    def test_old_environment_requests_one_time_regeneration(self) -> None:
        case = {"flow": {"type": "hit"}, "forcing": {"type": "none"}}
        old_lock = {"model": "nse", "profile": "cpu_mpi"}

        with self.assertRaisesRegex(RunCaseError, "Regenerate"):
            _select_case_profile(case, self.manifest, old_lock)

    def test_explicit_minimal_profile_is_not_silently_overridden(self) -> None:
        case = {"flow": {"type": "hit"}, "forcing": {"type": "none"}}
        explicit_lock = {**self.lock, "profile_explicit": True}

        with self.assertRaisesRegex(RunCaseError, "case capabilities"):
            _select_case_profile(case, self.manifest, explicit_lock)


class RunCasePrepareTests(unittest.TestCase):
    def test_locked_global_index_does_not_block_case_preparation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            errors = io.StringIO()
            with patch(
                "run_case.sync_environment_case",
                side_effect=GlobalCaseIndexError("case_index.csv is locked"),
            ), redirect_stderr(errors):
                index = _sync_global_case_index(
                    root,
                    {"case_index_path": "../case_index.csv"},
                    dry_run=False,
                )

            self.assertEqual(index, (root / "../case_index.csv").resolve())
            self.assertIn("[WARNING]", errors.getvalue())
            self.assertIn("Continuing", errors.getvalue())

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

    def test_newer_extension_regenerates_solver_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "cases" / "case0001"
            tools_dir = root / "tools"
            build_dir = root / "ScriptLibrary" / "BuildSolver"
            solver_dir = root / "SolverLibrary" / "GPE" / "gp3d"
            for path in (case_dir, tools_dir, build_dir, solver_dir):
                path.mkdir(parents=True, exist_ok=True)
            case_path = case_dir / "case.yaml"
            extension_path = case_dir / "chemistry.yaml"
            input_path = case_dir / "input.nml"
            case_path.write_text(
                "schema_version: 2\nphysics:\n  model: gpe\n"
                "extensions:\n  chemistry: chemistry.yaml\n",
                encoding="utf-8",
            )
            extension_path.write_text(
                "schema_version: 1\nextension: chemistry\n"
                "config:\n  model: none\n",
                encoding="utf-8",
            )
            input_path.write_text("old", encoding="ascii")
            current = input_path.stat().st_mtime_ns
            os.utime(extension_path, ns=(current + 1_000_000_000,) * 2)
            (tools_dir / "case_input.py").write_text("", encoding="utf-8")
            (solver_dir / "solver_manifest.yaml").write_text(
                "schema_version: 1\n", encoding="utf-8"
            )
            (build_dir / "model_catalog.yaml").write_text(
                "schema_version: 1\nmodels:\n  gpe:\n"
                "    library_subpath: GPE/gp3d\n"
                "    manifest: solver_manifest.yaml\n",
                encoding="utf-8",
            )
            lock = {
                "case_directory": "cases/case0001",
                "input_name": "input.nml",
                "model": "gpe",
                "profile": "cuda_single",
            }

            with patch(
                "run_case.subprocess.run", return_value=SimpleNamespace(returncode=0)
            ) as run:
                _prepare_input(root, lock, force=False, dry_run=False)

            run.assert_called_once()

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
