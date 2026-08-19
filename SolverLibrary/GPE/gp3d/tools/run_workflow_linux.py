#!/usr/bin/env python3
"""Run the JSON-described CMake workflow on Linux/HPC systems."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


def full_path(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def cmake_value(value: object) -> str:
    if isinstance(value, bool):
        return "ON" if value else "OFF"
    return str(value)


def enabled(value: object) -> bool:
    return str(value).strip().upper() in {"1", "ON", "TRUE", "YES"}


def run(command: list[str], cwd: Path | None, env: dict[str, str], dry_run: bool) -> None:
    print("+ " + shlex.join(command), flush=True)
    if not dry_run:
        subprocess.run(command, cwd=cwd, env=env, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs="?", default="workflow.json")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as stream:
        config = json.load(stream)
    if config.get("schema_version") != 1:
        raise ValueError(f"unsupported workflow schema_version: {config.get('schema_version')}")

    config_dir = config_path.parent
    project_root = full_path(config_dir, config["project_root"])
    cmake = config["cmake"]
    source_dir = full_path(project_root, cmake["source_directory"])
    build_dir = full_path(project_root, cmake["build_directory"])
    run_config = config["run"]
    working_dir = full_path(project_root, run_config["working_directory"])
    definitions = cmake["definitions"]
    use_mpi = enabled(definitions.get("USE_MPI", False))
    gpu_backend = str(definitions.get("GPU_BACKEND", "none")).strip().lower()

    if gpu_backend == "cuda" and use_mpi:
        raise ValueError("GPU_BACKEND=cuda is single-GPU only; set USE_MPI=OFF")
    if gpu_backend == "cufftmp" and not use_mpi:
        raise ValueError("GPU_BACKEND=cufftmp requires USE_MPI=ON")
    if run_config["use_mpi_launcher"] and not use_mpi:
        raise ValueError("MPI launcher requested while USE_MPI is disabled")

    print("GP3D Linux workflow")
    print(f"  config : {config_path}")
    print(f"  source : {source_dir}")
    print(f"  build  : {build_dir}")
    print(f"  MPI    : {use_mpi}")
    print(f"  GPU    : {gpu_backend}")

    env = os.environ.copy()
    for name, value in run_config["environment"].items():
        env[name] = str(value)
        print(f"  env    : {name}={value}")

    stages = config["stages"]
    cmake_command = str(cmake["command"])
    if stages["configure"]:
        command = [cmake_command, "-S", str(source_dir), "-B", str(build_dir)]
        if cmake["generator"]:
            command.extend(["-G", str(cmake["generator"])])
        command.extend(f"-D{name}={cmake_value(value)}" for name, value in definitions.items())
        command.extend(str(value) for value in cmake["configure_arguments"])
        run(command, None, env, args.dry_run)

    if stages["build"]:
        build = config["build"]
        command = [cmake_command, "--build", str(build_dir)]
        if build["configuration"]:
            command.extend(["--config", str(build["configuration"])])
        if build["target"]:
            command.extend(["--target", str(build["target"])])
        if int(build["parallel_jobs"]) > 0:
            command.extend(["--parallel", str(build["parallel_jobs"])])
        if build["clean_first"]:
            command.append("--clean-first")
        command.extend(str(value) for value in build["arguments"])
        run(command, None, env, args.dry_run)

    if stages["test"]:
        test = config["test"]
        command = [str(test["command"]), "--test-dir", str(build_dir)]
        if test["configuration"]:
            command.extend(["--build-config", str(test["configuration"])])
        command.extend(str(value) for value in test["arguments"])
        run(command, None, env, args.dry_run)

    if stages["run"]:
        executable = str(run_config["executable"])
        if executable == "auto":
            if gpu_backend == "cufftmp":
                executable = "gp3d_cufftmp"
            elif gpu_backend == "cuda":
                executable = "gp3d_cuda"
            elif use_mpi:
                executable = "gp3d_mpi"
            else:
                executable = "gp3d_sequential"
        executable_path = full_path(build_dir, executable)
        program_args: list[str] = []
        if run_config["input_file"]:
            program_args.append(str(full_path(project_root, run_config["input_file"])))
        program_args.extend(str(value) for value in run_config["program_arguments"])

        if run_config["use_mpi_launcher"]:
            command = [
                str(run_config["launcher"]),
                str(run_config["process_option"]),
                str(run_config["processes"]),
            ]
            command.extend(str(value) for value in run_config["launcher_arguments"])
            command.append(str(executable_path))
            command.extend(program_args)
        else:
            command = [str(executable_path), *program_args]
        run(command, working_dir, env, args.dry_run)

    print("Workflow completed successfully.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, subprocess.CalledProcessError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
