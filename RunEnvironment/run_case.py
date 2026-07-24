#!/usr/bin/env python3
"""Build, test, and run one generated execution environment."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the local BuildSolver copy for this case."
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--configure", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--all", action="store_true")
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
        include_tests = bool(lock.get("include_tests", False))
        if args.test and not include_tests:
            raise RunCaseError(
                "tests were not staged; regenerate with model.include_tests: true"
            )
        runner = root / "ScriptLibrary" / "BuildSolver" / "build_model.py"
        design = root / "ScriptLibrary" / "BuildSolver" / "build.local.yaml"
        case_dir = root / str(lock["case_directory"])
        input_path = case_dir / str(lock["input_name"])
        if not runner.is_file() or not design.is_file() or not input_path.is_file():
            raise RunCaseError("generated environment is missing runner, design, or input")

        command = [
            sys.executable,
            str(runner),
            str(design),
            "--model",
            str(lock["model"]),
            "--profile",
            str(lock["profile"]),
            "--input-file",
            str(input_path),
            "--run-dir",
            str(case_dir),
            "--processes",
            str(args.processes or lock.get("processes", 1)),
            "--omp-threads",
            str(args.omp_threads or lock.get("omp_threads", 1)),
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
    except (KeyError, OSError, RunCaseError, json.JSONDecodeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
