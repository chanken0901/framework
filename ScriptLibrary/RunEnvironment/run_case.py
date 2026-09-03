#!/usr/bin/env python3
"""Build, test, and run one generated execution environment."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from case_configuration import CaseConfigurationError, resolve_case_configuration
from global_case_index import (
    GlobalCaseIndexError,
    case_index_path,
    sync_environment_case,
)
from profile_selection import ProfileSelectionError, select_case_profile
from yaml_support import YamlFormatError, load_yaml


class RunCaseError(RuntimeError):
    """Raised when the generated environment is incomplete."""


def _load_lock(root: Path) -> dict[str, Any]:
    path = root / "environment.lock.json"
    if not path.is_file():
        raise RunCaseError(f"environment lock file not found: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise RunCaseError(f"invalid environment lock file: {path}")
    return value


def _model_manifest(root: Path, model: str) -> Path:
    catalog_path = root / "ScriptLibrary" / "BuildSolver" / "model_catalog.yaml"
    if not catalog_path.is_file():
        raise RunCaseError(f"model catalog not found: {catalog_path}")
    catalog = load_yaml(catalog_path)
    if not isinstance(catalog, dict):
        raise RunCaseError(f"invalid model catalog: {catalog_path}")
    models = catalog.get("models")
    if not isinstance(models, dict) or not isinstance(models.get(model), dict):
        raise RunCaseError(f"model {model!r} is not registered in {catalog_path}")
    record = models[model]
    library_subpath = record.get("library_subpath")
    manifest_name = record.get("manifest")
    if not library_subpath or not manifest_name:
        raise RunCaseError(f"model {model!r} has no solver manifest mapping")
    manifest = root / "SolverLibrary" / str(library_subpath) / str(manifest_name)
    if not manifest.is_file():
        raise RunCaseError(f"solver manifest not found: {manifest}")
    return manifest


def _select_case_profile(
    case: dict[str, Any], manifest: dict[str, Any], lock: dict[str, Any]
) -> tuple[str, set[str]]:
    try:
        return select_case_profile(case, manifest, lock)
    except ProfileSelectionError as exc:
        raise RunCaseError(str(exc)) from exc


def _resolve_case_profile(
    root: Path, lock: dict[str, Any]
) -> tuple[str, set[str]]:
    case_path = root / str(lock["case_directory"]) / "case.yaml"
    if not case_path.is_file():
        raise RunCaseError(f"case design not found: {case_path}")
    case = resolve_case_configuration(case_path).document
    manifest_path = _model_manifest(root, str(lock["model"]))
    manifest = load_yaml(manifest_path)
    if not isinstance(manifest, dict):
        raise RunCaseError(f"invalid solver manifest: {manifest_path}")
    return _select_case_profile(case, manifest, lock)


def _sync_global_case_index(
    root: Path, lock: dict[str, Any], *, dry_run: bool
) -> Path:
    index_path = case_index_path(root, lock)
    if dry_run:
        print(f"[DRY-RUN] Would update global case index: {index_path}")
        return index_path
    try:
        sync_environment_case(root, lock, index_path)
    except (GlobalCaseIndexError, OSError) as exc:
        # The shared CSV is a rebuildable catalogue, not a solver input.  In
        # particular, Windows applications such as Excel may temporarily lock
        # the destination and reject the atomic replace.  Do not prevent input
        # generation or a simulation for an auxiliary-index failure.
        print(f"[WARNING] Could not update global case index: {exc}", file=sys.stderr)
        print(
            "[WARNING] Continuing without updating case_index.csv; "
            "it will be synchronized from case.yaml on the next run.",
            file=sys.stderr,
        )
        return index_path
    print(f"[OK] Updated global case index: {index_path}")
    return index_path


def _prepare_input(
    root: Path,
    lock: dict[str, Any],
    *,
    force: bool,
    dry_run: bool,
) -> Path:
    case_dir = root / str(lock["case_directory"])
    case_path = case_dir / "case.yaml"
    input_path = case_dir / str(lock["input_name"])
    resolved_path = case_dir / "resolved_case.yaml"
    if not case_path.is_file():
        raise RunCaseError(f"case design not found: {case_path}")

    configuration = resolve_case_configuration(case_path)
    stale = (
        force
        or not input_path.is_file()
        or (bool(configuration.extension_paths) and not resolved_path.is_file())
        or any(
            source.stat().st_mtime_ns > input_path.stat().st_mtime_ns
            for source in configuration.source_paths
        )
    )
    if not stale:
        return input_path

    generator = root / "tools" / "case_input.py"
    if not generator.is_file():
        raise RunCaseError(f"case input generator not found: {generator}")
    model = str(lock["model"])
    command = [
        sys.executable,
        str(generator),
        "--case",
        str(case_path),
        "--manifest",
        str(_model_manifest(root, model)),
        "--model",
        model,
        "--profile",
        str(lock["profile"]),
    ]
    if input_path.exists():
        command.append("--overwrite")
    print("[CMD] " + subprocess.list2cmdline(command), flush=True)
    if dry_run:
        return input_path
    result = subprocess.run(command, cwd=str(root), check=False)
    if result.returncode != 0:
        raise RunCaseError(
            f"failed to generate solver input with exit code {result.returncode}"
        )
    return input_path


def _case_parallel_settings(
    root: Path, lock: dict[str, Any]
) -> tuple[int, int, bool]:
    case_path = root / str(lock["case_directory"]) / "case.yaml"
    if not case_path.is_file():
        raise RunCaseError(f"case design not found: {case_path}")
    case = resolve_case_configuration(case_path).document
    solver = case.get("solver", {})
    if not isinstance(solver, dict):
        raise RunCaseError("case solver section must be a YAML mapping")
    if "processes" in solver:
        raise RunCaseError(
            "solver.processes is no longer supported; use solver.mpi_processes"
        )

    # Older generated cases repeated these profile-derived values. Accept them
    # as compatibility assertions, but new case templates omit both fields.
    for key in ("use_mpi", "use_cuda"):
        expected = bool(lock.get(key, False))
        if key not in solver:
            continue
        actual = solver[key]
        if not isinstance(actual, bool):
            raise RunCaseError(f"solver.{key} must be true or false")
        if actual != expected:
            raise RunCaseError(
                f"solver.{key} does not match the generated environment; "
                "update the environment design and regenerate the environment"
            )

    use_openmp = solver.get("use_openmp", bool(lock.get("use_openmp", False)))
    if not isinstance(use_openmp, bool):
        raise RunCaseError("solver.use_openmp must be true or false")
    openmp_capable = bool(
        lock.get("openmp_capable", lock.get("use_openmp", False))
    )
    if use_openmp and not openmp_capable:
        raise RunCaseError(
            "solver.use_openmp=true requires an OpenMP-capable solver profile"
        )

    try:
        processes = int(solver.get("mpi_processes", 1))
        omp_threads = int(solver.get("omp_threads", 1))
    except (TypeError, ValueError) as exc:
        raise RunCaseError(
            "solver.mpi_processes and solver.omp_threads must be integers"
        ) from exc
    if processes < 1 or omp_threads < 1:
        raise RunCaseError(
            "solver.mpi_processes and solver.omp_threads must be positive"
        )
    use_mpi = bool(lock.get("use_mpi", False))
    if not use_mpi and processes != 1:
        raise RunCaseError("solver.mpi_processes must be 1 when MPI is disabled")
    if not use_openmp and omp_threads != 1:
        raise RunCaseError("solver.omp_threads must be 1 when OpenMP is disabled")
    return processes, omp_threads, use_openmp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the local BuildSolver copy for this case."
    )
    parser.add_argument(
        "--prepare",
        action="store_true",
        help="Regenerate the solver input from case.yaml and environment.lock.json.",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--configure", action="store_true", help="Configure CMake only.")
    parser.add_argument("--build", action="store_true", help="Configure and build only.")
    parser.add_argument("--test", action="store_true", help="Build and run tests.")
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run the existing executable without configuring or building.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Build and run; include tests when they were staged.",
    )
    parser.add_argument("--clean-first", action="store_true")
    parser.add_argument("--configuration", choices=["Debug", "Release", "RelWithDebInfo"])
    parser.add_argument("--processes", type=int)
    parser.add_argument("--omp-threads", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    try:
        lock = _load_lock(root)
        profile, requirements = _resolve_case_profile(root, lock)
        effective_lock = dict(lock)
        effective_lock["profile"] = profile
        if profile != str(lock.get("profile", "")):
            reason = ", ".join(sorted(requirements))
            print(
                f"[OK] Selected compatible solver profile: {profile} "
                f"(case requires: {reason})"
            )
        _sync_global_case_index(root, effective_lock, dry_run=args.dry_run)
        include_tests = bool(lock.get("include_tests", False))
        if args.test and not include_tests:
            raise RunCaseError(
                "tests were not staged; regenerate with solver.include_tests: true"
            )
        runner = root / "ScriptLibrary" / "BuildSolver" / "build_model.py"
        design = root / "ScriptLibrary" / "BuildSolver" / "build.local.yaml"
        case_dir = root / str(lock["case_directory"])
        input_path = _prepare_input(
            root,
            effective_lock,
            force=args.prepare,
            dry_run=args.dry_run,
        )
        prepare_only = args.prepare and not any(
            [
                args.validate_only,
                args.configure,
                args.build,
                args.test,
                args.run,
                args.all,
            ]
        )
        if prepare_only:
            return 0
        if (
            not runner.is_file()
            or not design.is_file()
            or (not input_path.is_file() and not args.dry_run)
        ):
            raise RunCaseError("generated environment is missing runner, design, or input")

        case_processes, case_omp_threads, case_use_openmp = (
            _case_parallel_settings(root, effective_lock)
        )
        processes = args.processes if args.processes is not None else case_processes
        omp_threads = (
            args.omp_threads
            if args.omp_threads is not None
            else case_omp_threads
        )
        if processes < 1 or omp_threads < 1:
            raise RunCaseError("runtime parallel counts must be positive")
        if not bool(lock.get("use_mpi", False)) and processes != 1:
            raise RunCaseError("--processes must be 1 when MPI is disabled")
        if not case_use_openmp and omp_threads != 1:
            raise RunCaseError(
                "--omp-threads must be 1 when solver.use_openmp is false"
            )
        if (
            str(lock.get("model", "")).lower() == "nse"
            and bool(lock.get("use_mpi", False))
            and processes < 4
        ):
            raise RunCaseError(
                "the current NSE y-z decomposition requires at least 4 MPI processes"
            )
        command = [
            sys.executable,
            str(runner),
            str(design),
            "--model",
            str(lock["model"]),
            "--profile",
            profile,
            "--input-file",
            str(input_path),
            "--run-dir",
            str(case_dir),
            "--processes",
            str(processes),
            "--omp-threads",
            str(omp_threads),
        ]
        if args.configuration:
            command.extend(["--configuration", args.configuration])
        if args.validate_only:
            command.append("--validate-only")
        if args.all:
            if include_tests:
                command.append("--all")
            else:
                command.extend(["--build", "--run"])
        else:
            for enabled, option in (
                (args.configure, "--configure"),
                (args.build, "--build"),
                (args.test, "--test"),
                (args.run, "--run"),
            ):
                if enabled:
                    command.append(option)
        if not any(
            [
                args.validate_only,
                args.configure,
                args.build,
                args.test,
                args.run,
                args.all,
            ]
        ):
            command.append("--build")
        if args.clean_first:
            command.append("--clean-first")
        if args.dry_run:
            command.append("--dry-run")

        print("[CMD] " + subprocess.list2cmdline(command), flush=True)
        return subprocess.run(command, cwd=str(root), check=False).returncode
    except (
        KeyError,
        OSError,
        RunCaseError,
        GlobalCaseIndexError,
        CaseConfigurationError,
        YamlFormatError,
        json.JSONDecodeError,
    ) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
