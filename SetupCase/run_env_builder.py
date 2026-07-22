#!/usr/bin/env python3
"""ケースとソルバーマニフェストから最小構成の実行環境を生成する。

``solver.profile`` が要求するコンポーネントだけをSolverLibraryからコピーし、
namelist、CMakeワークフロー、チェックサム、来歴を一つの生成ディレクトリへ
まとめる。Fortranモジュールの重複と依存欠落もコピー前に検査する。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from gp3d_case import CaseValidationError, case_digest, generate_input_namelist, nested
from yaml_support import YamlFormatError, load_yaml


FRAMEWORK_ROOT = Path(__file__).resolve().parents[2]
GENERATED_MARKER = ".generated_run_env.json"
# 選択したFortranソースだけで内部依存が閉じているかを確認するための字句パターン。
MODULE_PATTERN = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)
USE_PATTERN = re.compile(
    r"^\s*use(?:\s*,\s*[^:]*)?\s*(?:::\s*)?([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)


class BuildEnvironmentError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _native_copy_path(path: Path) -> str:
    """Use Win32 extended paths so generated environments may live under long NAS paths."""
    resolved = str(path.resolve())
    if os.name != "nt" or resolved.startswith("\\\\?\\"):
        return resolved
    if resolved.startswith("\\\\"):
        return "\\\\?\\UNC\\" + resolved[2:]
    return "\\\\?\\" + resolved


def copy_file(source: Path, destination: Path) -> None:
    try:
        shutil.copy2(_native_copy_path(source), _native_copy_path(destination))
    except OSError as exc:
        raise BuildEnvironmentError(
            f"failed to copy {source} to {destination}: {exc}"
        ) from exc


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise BuildEnvironmentError(f"{label} must be a YAML mapping")
    return value


def _path_within(root: Path, relative: str) -> Path:
    root_resolved = root.resolve()
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root_resolved)
    except ValueError as exc:
        raise BuildEnvironmentError(f"manifest path escapes solver root: {relative}") from exc
    if not candidate.is_file():
        raise BuildEnvironmentError(f"manifest file does not exist: {relative}")
    return candidate


def resolve_component_files(
    solver_root: Path,
    manifest: dict[str, Any],
    profile_name: str,
    include_tests: bool,
) -> tuple[list[str], list[str]]:
    """プロファイルをコンポーネントへ展開し、安全な相対ファイル一覧を返す。"""
    profiles = _mapping(manifest.get("profiles"), "manifest.profiles")
    if profile_name not in profiles:
        raise BuildEnvironmentError(
            f"unknown profile {profile_name!r}; available profiles: {sorted(profiles)}"
        )
    profile = _mapping(profiles[profile_name], f"profile {profile_name}")
    component_names = list(profile.get("components") or [])
    if include_tests:
        component_names.extend(profile.get("test_components") or [])

    components = _mapping(manifest.get("components"), "manifest.components")
    selected: list[str] = []
    seen: set[str] = set()
    for component_name in component_names:
        if component_name not in components:
            raise BuildEnvironmentError(
                f"profile {profile_name} references unknown component {component_name!r}"
            )
        component = _mapping(components[component_name], f"component {component_name}")
        files = component.get("files") or []
        if not isinstance(files, list):
            raise BuildEnvironmentError(f"component {component_name}.files must be a list")
        for relative in files:
            relative_text = str(relative).replace("\\", "/")
            _path_within(solver_root, relative_text)
            if relative_text not in seen:
                seen.add(relative_text)
                selected.append(relative_text)

    cmake_file = str(manifest.get("cmake_file", "CMakeLists.txt"))
    _path_within(solver_root, cmake_file)
    if cmake_file not in seen:
        selected.insert(0, cmake_file)
    return selected, component_names


def inspect_fortran_dependencies(
    solver_root: Path, selected_files: list[str]
) -> dict[str, dict[str, list[str]]]:
    """選択ソースのmodule/useを調べ、重複定義と未解決gp3d依存を拒否する。"""
    records: dict[str, dict[str, list[str]]] = {}
    providers: dict[str, str] = {}
    for relative in selected_files:
        if Path(relative).suffix.lower() not in {".f90", ".f95", ".f03", ".f08"}:
            continue
        text = (solver_root / relative).read_text(encoding="utf-8")
        modules = sorted({name.lower() for name in MODULE_PATTERN.findall(text)})
        uses = sorted({name.lower() for name in USE_PATTERN.findall(text)})
        records[relative] = {"modules": modules, "uses": uses}
        for module in modules:
            if module in providers:
                raise BuildEnvironmentError(
                    f"selected sources provide duplicate Fortran module {module}: "
                    f"{providers[module]} and {relative}"
                )
            providers[module] = relative

    missing: list[str] = []
    for relative, record in records.items():
        for module in record["uses"]:
            if module.startswith("gp3d_") and module not in providers:
                missing.append(f"{relative} uses {module}")
    if missing:
        raise BuildEnvironmentError(
            "selected profile has unresolved GP3D module dependencies:\n  - "
            + "\n  - ".join(missing)
        )
    return records


def git_state(path: Path) -> dict[str, Any]:
    """利用可能な場合だけGitコミットとdirty状態を来歴用に取得する。"""
    def run(*args: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(path), *args],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        return result.stdout.strip() if result.returncode == 0 else ""

    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain") if commit else ""
    return {"commit": commit or None, "dirty": bool(status) if commit else None}


def resolve_solver(
    project_path: Path,
    project: dict[str, Any],
    machine: dict[str, Any],
    solver_id: str,
    solver_override: Path | None,
    library_override: Path | None,
) -> tuple[Path, Path]:
    if solver_override is not None:
        solver_root = solver_override.resolve()
        manifest_path = solver_root / "solver_manifest.yaml"
        return solver_root, manifest_path

    catalog = _mapping(project.get("solver_catalog"), "project.solver_catalog")
    if solver_id not in catalog:
        raise BuildEnvironmentError(
            f"solver {solver_id!r} is not in project_schema.yaml; available: {sorted(catalog)}"
        )
    entry = _mapping(catalog[solver_id], f"solver_catalog.{solver_id}")
    if entry.get("root"):
        solver_root = (project_path.parent / str(entry["root"])).resolve()
    else:
        machine_paths = _mapping(machine.get("paths", {}), "machine.paths")
        library_text = (
            str(library_override) if library_override is not None else
            str(machine_paths.get("solver_library_root") or os.environ.get("SOLVER_LIBRARY_ROOT", ""))
        )
        if not library_text:
            raise BuildEnvironmentError(
                "solver library root is not configured; set machine.paths.solver_library_root, "
                "SOLVER_LIBRARY_ROOT, or --solver-library"
            )
        library_root = Path(os.path.expandvars(library_text))
        if not library_root.is_absolute():
            library_root = project_path.parent / library_root
        subpath = str(entry.get("library_subpath", solver_id))
        solver_root = (library_root / subpath).resolve()
    manifest_path = solver_root / str(entry.get("manifest", "solver_manifest.yaml"))
    if not manifest_path.is_file():
        raise BuildEnvironmentError(f"solver manifest not found: {manifest_path}")
    return solver_root, manifest_path


def workflow_config(
    case: dict[str, Any],
    profile: dict[str, Any],
    machine: dict[str, Any],
    include_tests: bool,
) -> dict[str, Any]:
    """ケース、プロファイル、マシン設定を実行用workflow.jsonへ統合する。"""
    machine_cmake = _mapping(machine.get("cmake", {}), "machine.cmake")
    machine_tools = _mapping(machine.get("tools", {}), "machine.tools")
    libraries = _mapping(machine.get("libraries", {}), "machine.libraries")
    definitions = dict(_mapping(profile.get("cmake", {}), "profile.cmake"))
    definitions["CMAKE_BUILD_TYPE"] = str(machine_cmake.get("build_type", "Release"))
    compiler = _mapping(machine.get("compiler", {}), "machine.compiler")
    definitions["CMAKE_Fortran_COMPILER"] = str(
        machine_cmake.get("fortran_compiler", compiler.get("fortran", "gfortran"))
    )
    definitions["BUILD_TESTING"] = include_tests

    if str(definitions.get("FFT_BACKEND", "")).lower() == "fftw" and libraries.get("fftw_root"):
        definitions["FFTW_ROOT"] = libraries["fftw_root"]
    if bool(definitions.get("USE_MPI", False)) and libraries.get("mpi_root"):
        definitions["MPI_ROOT"] = libraries["mpi_root"]
    gpu_backend = str(definitions.get("GPU_BACKEND", "none")).lower()
    use_cuda = gpu_backend in {"cuda", "cufftmp"}
    if use_cuda:
        if libraries.get("cuda_compiler"):
            definitions["CMAKE_CUDA_COMPILER"] = libraries["cuda_compiler"]
        if libraries.get("cuda_toolkit_root"):
            definitions["CUDAToolkit_ROOT"] = libraries["cuda_toolkit_root"]
        if libraries.get("cuda_architectures") is not None:
            definitions["GP3D_CUDA_ARCHITECTURES"] = libraries["cuda_architectures"]
    if gpu_backend == "cufftmp":
        if libraries.get("cufftmp_root"):
            definitions["CUFFTMP_ROOT"] = libraries["cufftmp_root"]
        if libraries.get("nvshmem_root"):
            definitions["NVSHMEM_ROOT"] = libraries["nvshmem_root"]
        if libraries.get("cufftmp_api"):
            definitions["CUFFTMP_API"] = libraries["cufftmp_api"]

    use_mpi = bool(definitions.get("USE_MPI", False))
    processes = int(nested(case, "solver.processes", 1))
    environment: dict[str, Any] = {}
    if bool(nested(case, "solver.use_openmp", False)):
        environment["OMP_NUM_THREADS"] = int(nested(case, "solver.omp_threads", 1))
    if gpu_backend == "cufftmp" and libraries.get("nvshmem_symmetric_size"):
        environment["NVSHMEM_SYMMETRIC_SIZE"] = str(libraries["nvshmem_symmetric_size"])

    fresh = bool(machine_cmake.get("configure_fresh", False))
    return {
        "$schema": "./workflow.schema.json",
        "schema_version": 1,
        "project_root": ".",
        "stages": {
            "configure": True,
            "build": True,
            "test": include_tests,
            "run": bool(nested(case, "execution.run", True)),
        },
        "cmake": {
            "command": str(machine_cmake.get("command", "cmake")),
            "source_directory": ".",
            "build_directory": "build",
            "generator": str(machine_cmake.get("generator", "Ninja")),
            "definitions": definitions,
            "configure_arguments": ["--fresh"] if fresh else [],
        },
        "build": {
            "configuration": str(machine_cmake.get("build_type", "Release")),
            "target": "",
            "parallel_jobs": 0,
            "clean_first": False,
            "arguments": [],
        },
        "test": {
            "command": str(machine_tools.get("ctest", "ctest")),
            "configuration": str(machine_cmake.get("build_type", "Release")),
            "arguments": ["--output-on-failure"],
        },
        "run": {
            "working_directory": ".",
            "executable": "auto",
            "input_file": "input.nml",
            "use_mpi_launcher": use_mpi,
            "launcher": str(machine_tools.get("mpi_launcher", "mpiexec")),
            "process_option": "-n",
            "processes": processes,
            "launcher_arguments": [],
            "program_arguments": [],
            "environment": environment,
        },
    }


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _generated_readme(case_id: str, profile_name: str, include_tests: bool) -> str:
    tests = "有効" if include_tests else "無効"
    return f"""# 自動生成されたGP3D実行環境

- ケース: `{case_id}`
- ソルバープロファイル: `{profile_name}`
- CTest: `{tests}`

このディレクトリは自動生成されています。コピーされたソースコードを直接変更せず、
原本の`case.yaml`またはSolverLibraryのマニフェストを編集してから再生成してください。

## Windowsでの実行

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\\tools\\run_workflow.ps1 .\\workflow.json
```

## CMakeを手動実行する場合

確定したCMake定義は`workflow.json`で確認できます。計算結果は相対パス`output/`へ
出力され、生成ソースのマニフェストには含まれません。
"""


def build_environment(args: argparse.Namespace) -> Path:
    """依存解決、検証、コピー、成果物生成を一時ディレクトリ上で完結させる。"""
    project_path = Path(args.project).resolve()
    case_path = Path(args.case).resolve()
    project = _mapping(load_yaml(project_path), "project schema")
    case = _mapping(load_yaml(case_path), "case YAML")
    machine_relative = nested(project, "paths.default_machine", "config/machine.local.yaml")
    machine_path = Path(args.machine).resolve() if args.machine else (
        project_path.parent / str(machine_relative)
    ).resolve()
    machine = _mapping(load_yaml(machine_path), "machine YAML")
    solver_id = str(nested(case, "solver.implementation", "gp3d"))
    solver_root, manifest_path = resolve_solver(
        project_path,
        project,
        machine,
        solver_id,
        Path(args.solver_root) if args.solver_root else None,
        Path(args.solver_library) if args.solver_library else None,
    )
    manifest = _mapping(load_yaml(manifest_path), "solver manifest")
    profile_name = args.profile or str(nested(case, "solver.profile", ""))
    if not profile_name:
        raise BuildEnvironmentError("solver.profile is required in case.yaml")

    selected_files, selected_components = resolve_component_files(
        solver_root, manifest, profile_name, args.include_tests
    )
    dependency_records = inspect_fortran_dependencies(solver_root, selected_files)
    generate_input_namelist(case, manifest, profile_name)

    output_root = project_path.parent / str(
        nested(project, "paths.generated_run_environments", "generated")
    )
    default_name = f"{case['case_id']}_{profile_name}"
    output = Path(args.output).resolve() if args.output else (output_root / default_name).resolve()
    if args.dry_run:
        print(f"case={case['case_id']}")
        print(f"profile={profile_name}")
        print(f"solver_root={solver_root}")
        print(f"output={output}")
        for relative in selected_files:
            print(f"  {relative}")
        return output

    if output.exists():
        if not args.overwrite:
            raise BuildEnvironmentError(f"output already exists; use --overwrite: {output}")
        if not (output / GENERATED_MARKER).is_file():
            raise BuildEnvironmentError(
                f"refusing to replace a directory without {GENERATED_MARKER}: {output}"
            )
        shutil.rmtree(output)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.tmp-{uuid.uuid4().hex}"
    temporary.mkdir(parents=True)
    try:
        source_entries: list[dict[str, Any]] = []
        for relative in selected_files:
            source = solver_root / relative
            destination = temporary / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            copy_file(source, destination)
            entry: dict[str, Any] = {
                "path": relative,
                "size": source.stat().st_size,
                "sha256": sha256_file(source),
            }
            if relative in dependency_records:
                entry.update(dependency_records[relative])
            source_entries.append(entry)

        input_text = generate_input_namelist(case, manifest, profile_name)
        (temporary / "input.nml").write_text(input_text, encoding="ascii")
        copy_file(case_path, temporary / "case.yaml")
        copy_file(machine_path, temporary / "machine.yaml")
        copy_file(project_path, temporary / "project_schema.yaml")
        copy_file(manifest_path, temporary / "solver_manifest.yaml")

        profile = _mapping(manifest["profiles"][profile_name], f"profile {profile_name}")
        workflow = workflow_config(case, profile, machine, args.include_tests)
        _write_json(temporary / "workflow.json", workflow)

        generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        provenance = {
            "schema_version": 1,
            "generated_at_utc": generated_at,
            "case_id": case["case_id"],
            "case_sha256": case_digest(case),
            "solver_id": manifest.get("solver_id"),
            "solver_api_version": manifest.get("solver_api_version"),
            "profile": profile_name,
            "machine_id": machine.get("machine_id", machine.get("machine_name")),
            "library_git": git_state(solver_root),
            "solver_manifest_sha256": sha256_file(manifest_path),
        }
        source_manifest = {
            "schema_version": 1,
            "solver_id": manifest.get("solver_id"),
            "profile": profile_name,
            "components": selected_components,
            "files": source_entries,
        }
        _write_json(temporary / "provenance.json", provenance)
        _write_json(temporary / "source_manifest.json", source_manifest)
        _write_json(
            temporary / GENERATED_MARKER,
            {"schema_version": 1, "generated_at_utc": generated_at, "case_id": case["case_id"]},
        )
        (temporary / "README.md").write_text(
            _generated_readme(str(case["case_id"]), profile_name, args.include_tests),
            encoding="utf-8",
        )
        temporary.replace(output)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise

    print(f"[OK] Generated run environment: {output}")
    print(f"     profile: {profile_name}")
    print(f"     files:   {len(selected_files)}")
    return output


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a minimal GP3D run environment from case and solver YAML manifests."
    )
    parser.add_argument("--case", required=True, help="Path to case.yaml")
    parser.add_argument(
        "--project",
        default="project_schema.yaml",
        help="Path to project_schema.yaml",
    )
    parser.add_argument("--machine", help="Override machine YAML")
    parser.add_argument("--solver-root", help="Override SolverLibrary GP3D root")
    parser.add_argument("--solver-library", help="Override the SolverLibrary root")
    parser.add_argument("--profile", help="Override solver.profile from case.yaml")
    parser.add_argument("--output", help="Generated run environment directory")
    parser.add_argument("--include-tests", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    try:
        build_environment(parse_args(argv))
    except (BuildEnvironmentError, CaseValidationError, YamlFormatError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
