#!/usr/bin/env python3
"""Create a self-contained run environment from a NAS framework library."""

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

from case_input import CaseInputError, _profile_settings, render_case_input
from environment_options import (
    OPTION_CATEGORIES,
    OptionCatalogError,
    catalog_as_json,
    catalog_as_text,
    load_option_catalogs,
    resolve_design_options,
    validate_design_selections,
)
from global_case_index import GlobalCaseIndexError, sync_environment_case
from yaml_support import YamlFormatError, load_yaml


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DESIGN = SCRIPT_DIR / "environment.gpe.yaml"
DEFAULT_OPTION_CATALOG = SCRIPT_DIR / "environment_options.yaml"
GENERATED_MARKER = ".generated_run_environment.json"
FORTRAN_SUFFIXES = {".f90", ".f95", ".f03", ".f08"}
MODULE_PATTERN = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)
USE_PATTERN = re.compile(
    r"^\s*use(?:\s*,\s*[^:]*)?\s*(?:::\s*)?([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)


class EnvironmentError(RuntimeError):
    """Raised when an execution environment cannot be generated safely."""


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise EnvironmentError(f"{label} must be a YAML mapping")
    return value


def _sequence(value: Any, label: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise EnvironmentError(f"{label} must be a YAML list")
    return value


def _expand_path(value: str, base: Path, label: str) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(value))
    unresolved = re.search(r"\$\{[^}]+\}|%[^%]+%", expanded)
    if unresolved:
        raise EnvironmentError(
            f"undefined environment variable in {label}: {unresolved.group(0)}"
        )
    path = Path(expanded)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _numbered_case_destination(
    root: Path,
    model: str,
    width: int = 4,
    start: int = 1,
) -> tuple[Path, str]:
    if width < 1:
        raise EnvironmentError("destination.case_number_width must be positive")
    if start < 1:
        raise EnvironmentError("destination.case_number_start must be positive")

    model_name = re.sub(r"[^0-9A-Za-z_-]+", "_", model.strip().lower()).strip("_")
    if not model_name:
        raise EnvironmentError("model.name cannot form a case directory name")

    number = start
    while True:
        case_id = f"case{number:0{width}d}"
        output = root / f"{model_name}_{case_id}"
        archive_candidates = (
            Path(f"{output}.zip"),
            Path(f"{output}.tar"),
            Path(f"{output}.tar.gz"),
            Path(f"{output}.tgz"),
        )
        if not output.exists() and not any(path.exists() for path in archive_candidates):
            return output, case_id
        number += 1


def _requested_case_destination(
    root: Path,
    model: str,
    requested_case_id: str,
    width: int = 4,
    start: int = 1,
) -> tuple[Path, str]:
    if width < 1:
        raise EnvironmentError("destination.case_number_width must be positive")
    if start < 1:
        raise EnvironmentError("destination.case_number_start must be positive")

    model_name = re.sub(r"[^0-9A-Za-z_-]+", "_", model.strip().lower()).strip("_")
    if not model_name:
        raise EnvironmentError("model.name cannot form a case directory name")

    case_text = requested_case_id.strip().lower()
    if case_text.startswith("case"):
        case_text = case_text[4:]
    if not case_text.isdigit():
        raise EnvironmentError(
            "--case-id must be a case number such as case0015 or 0015"
        )
    number = int(case_text)
    if number < start:
        raise EnvironmentError(
            f"--case-id must be case{start:0{width}d} or greater"
        )

    case_id = f"case{number:0{width}d}"
    return root / f"{model_name}_{case_id}", case_id


def _global_case_index_path(
    destination_cfg: dict[str, Any],
    design_path: Path,
    output: Path,
) -> Path:
    configured = destination_cfg.get("case_index")
    if configured:
        index_path = _expand_path(
            str(configured), design_path.parent, "destination.case_index"
        )
    else:
        index_path = output.parent / "case_index.csv"
    try:
        index_path.relative_to(output)
    except ValueError:
        return index_path
    raise EnvironmentError(
        "destination.case_index must be outside the generated execution environment"
    )


def _portable_path_reference(path: Path, root: Path) -> str:
    try:
        return os.path.relpath(path, root).replace("\\", "/")
    except ValueError:
        return str(path)


def _safe_source(root: Path, relative: str, label: str) -> Path:
    clean = str(relative).replace("\\", "/")
    candidate = (root / clean).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as exc:
        raise EnvironmentError(f"{label} escapes framework root: {relative}") from exc
    if not candidate.is_file():
        raise EnvironmentError(f"{label} not found: {candidate}")
    return candidate


def _safe_solver_file(solver_root: Path, relative: str) -> Path:
    clean = str(relative).replace("\\", "/")
    candidate = (solver_root / clean).resolve()
    try:
        candidate.relative_to(solver_root.resolve())
    except ValueError as exc:
        raise EnvironmentError(f"solver path escapes solver root: {relative}") from exc
    if not candidate.is_file():
        raise EnvironmentError(f"solver file not found: {candidate}")
    return candidate


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _copy_file(
    source: Path,
    destination: Path,
    records: list[dict[str, Any]],
    framework_root: Path,
    category: str,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    try:
        source_name = source.resolve().relative_to(framework_root.resolve()).as_posix()
    except ValueError:
        source_name = str(source.resolve())
    destination_name = destination.as_posix()
    for parent in destination.parents:
        if parent.name.startswith(".") and ".tmp-" in parent.name:
            destination_name = destination.relative_to(parent).as_posix()
            break
    records.append(
        {
            "category": category,
            "source": source_name,
            "destination": destination_name,
            "size": source.stat().st_size,
            "sha256": _sha256(source),
        }
    )


def _tree_files(root: Path) -> list[Path]:
    if not root.is_dir():
        raise EnvironmentError(f"required tool directory not found: {root}")
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and path.suffix.lower() not in {".pyc", ".pyo"}
        and path.name != "build.local.yaml"
    )


def _selected_solver_files(
    solver_root: Path,
    manifest: dict[str, Any],
    profile_name: str,
    include_tests: bool,
) -> tuple[list[str], list[str]]:
    profiles = _mapping(manifest.get("profiles"), "solver manifest.profiles")
    if profile_name not in profiles:
        raise EnvironmentError(
            f"unknown profile {profile_name!r}; available: {sorted(profiles)}"
        )
    profile = _mapping(profiles[profile_name], f"solver profile {profile_name}")
    component_names = [
        str(value)
        for value in _sequence(profile.get("components"), "profile.components")
    ]
    if include_tests:
        component_names.extend(
            str(value)
            for value in _sequence(
                profile.get("test_components"), "profile.test_components"
            )
        )

    components = _mapping(manifest.get("components"), "solver manifest.components")
    selected: list[str] = []
    seen: set[str] = set()
    for name in component_names:
        if name not in components:
            raise EnvironmentError(f"profile references unknown component: {name}")
        component = _mapping(components[name], f"component {name}")
        for value in _sequence(component.get("files"), f"component {name}.files"):
            relative = str(value).replace("\\", "/")
            _safe_solver_file(solver_root, relative)
            if relative not in seen:
                seen.add(relative)
                selected.append(relative)

    distribution = _mapping(
        manifest.get("distribution", {}), "solver manifest.distribution"
    )
    support_files = [
        str(value).replace("\\", "/")
        for value in _sequence(
            distribution.get("required_files"), "distribution.required_files"
        )
    ]
    cmake_file = str(manifest.get("cmake_file") or "CMakeLists.txt").replace("\\", "/")
    support_files.insert(0, cmake_file)
    for relative in support_files:
        _safe_solver_file(solver_root, relative)
        if relative not in seen:
            seen.add(relative)
            selected.insert(0, relative)
    return selected, component_names


def _inspect_dependencies(
    solver_root: Path,
    manifest: dict[str, Any],
    selected_files: list[str],
) -> dict[str, dict[str, list[str]]]:
    components = _mapping(manifest.get("components"), "solver manifest.components")
    all_providers: set[str] = set()
    for component in components.values():
        record = _mapping(component, "solver component")
        for value in _sequence(record.get("files"), "component.files"):
            path = _safe_solver_file(solver_root, str(value))
            if path.suffix.lower() in FORTRAN_SUFFIXES:
                text = path.read_text(encoding="utf-8", errors="replace")
                all_providers.update(name.lower() for name in MODULE_PATTERN.findall(text))

    providers: dict[str, str] = {}
    records: dict[str, dict[str, list[str]]] = {}
    for relative in selected_files:
        path = _safe_solver_file(solver_root, relative)
        if path.suffix.lower() not in FORTRAN_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        modules = sorted({name.lower() for name in MODULE_PATTERN.findall(text)})
        uses = sorted({name.lower() for name in USE_PATTERN.findall(text)})
        records[relative] = {"modules": modules, "uses": uses}
        for module in modules:
            if module in providers and providers[module] != relative:
                raise EnvironmentError(
                    f"Fortran module {module} is defined by both "
                    f"{providers[module]} and {relative}"
                )
            providers[module] = relative

    missing = [
        f"{relative} uses {module}"
        for relative, record in records.items()
        for module in record["uses"]
        if module in all_providers and module not in providers
    ]
    if missing:
        raise EnvironmentError(
            "selected profile has missing internal modules:\n  - "
            + "\n  - ".join(missing)
        )
    return records


def _git_state(path: Path) -> dict[str, Any]:
    def run(*arguments: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(path), *arguments],
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


def _local_build_design(
    model: str,
    profile: str,
    configuration: str,
    include_tests: bool,
    parallel_jobs: int,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "design_id": "generated_local_execution_environment",
        "selected_model": model,
        "model_catalog": "model_catalog.yaml",
        "machine_profile": "../../config/machine.yaml",
        "solver_library_root": "../../SolverLibrary",
        "models": {
            model: {
                "profile": profile,
                "input_file": None,
                "cmake_overrides": {},
            }
        },
        "build": {
            "configuration": configuration,
            "output_root": "../../build",
            "directory_name": "{model}-{profile}-{configuration}",
            "warnings": True,
            "tests": include_tests,
            "parallel_jobs": parallel_jobs,
            "configure_fresh": True,
        },
        "run": {
            "enabled": False,
            "working_directory": None,
            "mpi_processes": 1,
            "omp_threads": 1,
            "launcher_arguments": [],
            "program_arguments": [],
            "environment": {},
        },
    }


def _parallel_features(
    design: dict[str, Any],
    manifest: dict[str, Any],
    profile: str,
) -> tuple[bool, bool, bool]:
    parallel = _mapping(design.get("parallel"), "parallel")
    values: dict[str, bool] = {}
    for name in ("use_mpi", "use_openmp", "use_cuda"):
        value = parallel.get(name)
        if not isinstance(value, bool):
            raise EnvironmentError(f"parallel.{name} must be true or false")
        values[name] = value

    _, profile_mpi, profile_openmp, backend = _profile_settings(manifest, profile)
    profile_cuda = backend in {"cuda", "cufft", "cufftmp"}
    if values["use_mpi"] != profile_mpi:
        raise EnvironmentError(
            f"parallel.use_mpi={values['use_mpi']} does not match solver "
            f"profile {profile!r} (use_mpi={profile_mpi})"
        )
    if values["use_cuda"] != profile_cuda:
        raise EnvironmentError(
            f"parallel.use_cuda={values['use_cuda']} does not match solver "
            f"profile {profile!r} (use_cuda={profile_cuda})"
        )
    if values["use_openmp"] and not profile_openmp:
        raise EnvironmentError(
            f"parallel.use_openmp=true requires an OpenMP-capable profile; "
            f"{profile!r} is not OpenMP-capable"
        )
    return values["use_mpi"], values["use_openmp"], values["use_cuda"]


def _slurm_script(
    scheduler: dict[str, Any],
    include_tests: bool,
    use_mpi: bool,
    use_openmp: bool,
) -> str:
    directives = [
        ("job-name", scheduler.get("job_name", "solver-case")),
        ("account", scheduler.get("account")),
        ("partition", scheduler.get("partition")),
        ("time", scheduler.get("time_limit", "01:00:00")),
        ("output", scheduler.get("output", "slurm-%j.out")),
    ]
    lines = ["#!/bin/bash"]
    for name, value in directives:
        if value is not None and value != "":
            lines.append(f"#SBATCH --{name}={value}")
    for value in _sequence(
        scheduler.get("extra_directives"), "scheduler.extra_directives"
    ):
        lines.append(f"#SBATCH {value}")
    lines.extend(["", "set -euo pipefail", 'cd "$(dirname "$0")"'])
    modules = [
        str(value)
        for value in _sequence(scheduler.get("modules"), "scheduler.modules")
    ]
    if modules:
        lines.append("module purge")
        lines.extend(f"module load {name}" for name in modules)
    lines.append("")
    run_command = "python3 tools/run_case.py --run"
    if use_mpi:
        lines.append(': "${SLURM_NTASKS:?submit with --ntasks or --ntasks-per-node}"')
        run_command += ' --processes "${SLURM_NTASKS}"'
    if use_openmp:
        lines.append(': "${SLURM_CPUS_PER_TASK:?submit with --cpus-per-task}"')
        run_command += ' --omp-threads "${SLURM_CPUS_PER_TASK}"'
    lines.extend([run_command, ""])
    return "\n".join(lines)


def _generated_readme(
    model: str,
    profile: str,
    case_id: str,
    use_mpi: bool,
    use_openmp: bool,
) -> str:
    enabled_features = ", ".join(
        name
        for name, enabled in (("MPI", use_mpi), ("OpenMP", use_openmp))
        if enabled
    ) or "serial"
    return f"""# Generated execution environment

This directory is a disposable execution copy. The canonical source remains in
the framework NAS.

- Model: `{model}`
- Solver profile: `{profile}`
- Case: `{case_id}`
- Parallel features: `{enabled_features}`
- Shared case index: `../case_index.csv`

## Workstation

```powershell
python .\\tools\\run_case.py --prepare
python .\\tools\\run_case.py --validate-only
python .\\tools\\run_case.py --build
python .\\tools\\run_case.py --run
```

Set `solver.mpi_processes` and `solver.omp_threads` in
`cases/{case_id}/case.yaml` before running. Command-line options
`--processes` and `--omp-threads` temporarily override those values.

## ParaView preview

Run these commands from this generated environment root. First verify the
postprocessing tool, display the child command, and inspect the selected SLFs:

```powershell
Test-Path .\\environment.lock.json
python .\\tools\\postprocess_case.py --version
python .\\tools\\postprocess_case.py --dry-run
python .\\tools\\postprocess_case.py --inspect-only --steps latest
```

Version 2.0.0 or later supports the same step selection and inspection flow for
both GPE and NSE. Convert the latest complete SLF step with a spatial stride of
two:

```powershell
python .\\tools\\postprocess_case.py
```

Convert selected full-resolution steps:

```powershell
python .\\tools\\postprocess_case.py --steps 0,1000 --stride 1
```

The step syntax is `all`, `latest`, a comma-separated list such as
`0,500,1000`, or an inclusive range such as `0:1000:100`. ParaView output is
written below `cases/{case_id}/paraview`; open `collection.pvd`. Incomplete MPI
steps are skipped. A generated environment is a copy and does not automatically
follow later FrameWork updates; regenerate it when the reported version is older.

## NSE turbulence statistics

Compute a CSV time series from every complete NSE SLF step. Gamma and the
reference Reynolds number are read from `cases/{case_id}/case.yaml`:

```powershell
python .\\tools\\postprocess_case.py --task statistics
```

Run ParaView conversion and statistics together:

```powershell
python .\\tools\\postprocess_case.py --task all
```

## Linux / HPC

```bash
python3 tools/run_case.py --prepare
python3 tools/run_case.py --validate-only
python3 tools/run_case.py --build
python3 tools/run_case.py --run
python3 tools/postprocess_case.py
python3 tools/postprocess_case.py --task statistics
```

Build once before submitting production jobs. When `submit.slurm` exists, it
runs the existing executable with `python3 tools/run_case.py --run`; submit it
with `sbatch submit.slurm`. Calculation output is written below
`cases/{case_id}/output` because the case directory is used as the run working
directory. The shared index is outside this environment and is refreshed from
`case.yaml` before prepare, validation, build, test, and run operations.
"""


def _resolve_case_id(before: set[str], cases_root: Path) -> str:
    after = {
        path.name
        for path in cases_root.iterdir()
        if path.is_dir() and (path / "case.yaml").is_file()
    }
    created = sorted(after - before)
    if len(created) != 1:
        raise EnvironmentError(
            f"case generator created {len(created)} cases; expected exactly one"
        )
    return created[0]


def _create_case(
    temporary: Path,
    design_path: Path,
    case_cfg: dict[str, Any],
    model: str,
    profile: str,
    manifest_path: Path,
    use_mpi: bool,
    use_openmp: bool,
    use_cuda: bool,
) -> tuple[str, str]:
    cases_root = temporary / "cases"
    cases_root.mkdir(parents=True, exist_ok=True)
    create = bool(case_cfg.get("create", True))
    if create:
        generator = (
            temporary
            / "ScriptLibrary"
            / "SetupCase"
            / "create_case_from_template.py"
        )
        template = temporary / "templates" / "case_template.yaml"
        before = {path.name for path in cases_root.iterdir() if path.is_dir()}
        command = [
            sys.executable,
            str(generator),
            "--cases-root",
            "cases",
            "--case-index",
            ".case_index.build.csv",
            "--template",
            "templates/case_template.yaml",
            "--label",
            str(case_cfg.get("label", "baseline")),
            "--description",
            str(case_cfg.get("description", "")),
            "--project-name",
            str(case_cfg.get("project_name", temporary.name)),
            "--data-root",
            str(case_cfg.get("data_root", "results")),
            "--model",
            model,
            "--solver-profile",
            profile,
        ]
        command.append("--use-mpi" if use_mpi else "--no-use-mpi")
        command.append("--use-openmp" if use_openmp else "--no-use-openmp")
        command.append("--use-cuda" if use_cuda else "--no-use-cuda")
        command.extend(
            [
                "--mpi-processes",
                "4" if use_mpi else "1",
                "--omp-threads",
                "1",
            ]
        )
        if case_cfg.get("id"):
            command.extend(["--case-id", str(case_cfg["id"])])
        result = subprocess.run(command, cwd=str(temporary), check=False)
        if result.returncode != 0:
            raise EnvironmentError("SetupCase/create_case_from_template.py failed")
        (temporary / ".case_index.build.csv").unlink(missing_ok=True)
        case_id = _resolve_case_id(before, cases_root)
    else:
        source_text = str(case_cfg.get("source") or "")
        if not source_text:
            raise EnvironmentError("case.source is required when case.create is false")
        source = _expand_path(source_text, design_path.parent, "case.source")
        if not source.is_file():
            raise EnvironmentError(f"case source not found: {source}")
        case_document = _mapping(load_yaml(source), "source case YAML")
        case_id = str(case_document.get("case_id") or case_cfg.get("id") or "")
        if not case_id:
            raise EnvironmentError("source case YAML has no case_id")
        case_dir = cases_root / case_id
        case_dir.mkdir(parents=True)
        shutil.copy2(source, case_dir / "case.yaml")
        (case_dir / "notes.md").write_text(f"# Notes for {case_id}\n", encoding="utf-8")

    case_path = cases_root / case_id / "case.yaml"
    input_name, input_text = render_case_input(
        case_path, manifest_path, model, profile
    )
    (case_path.parent / input_name).write_text(input_text, encoding="ascii")
    return case_id, input_name


def _create_archive(
    output: Path,
    design_path: Path,
    archive_cfg: dict[str, Any],
    cli_enabled: bool,
    cli_format: str | None,
) -> Path | None:
    enabled = cli_enabled or bool(archive_cfg.get("enabled", False))
    if not enabled:
        return None
    archive_format = str(cli_format or archive_cfg.get("format", "gztar"))
    if archive_format not in {"zip", "gztar"}:
        raise EnvironmentError("archive format must be zip or gztar")
    destination_text = archive_cfg.get("output")
    if destination_text:
        base = _expand_path(str(destination_text), design_path.parent, "archive.output")
    else:
        base = output.parent / output.name
    suffixes = {".zip", ".gz", ".tgz", ".tar"}
    while base.suffix.lower() in suffixes:
        base = base.with_suffix("")
    archive = shutil.make_archive(
        str(base), archive_format, root_dir=str(output.parent), base_dir=output.name
    )
    return Path(archive)


def prepare(args: argparse.Namespace) -> Path:
    design_path = (
        Path(args.design).resolve() if args.design else DEFAULT_DESIGN.resolve()
    )
    source_design = _mapping(load_yaml(design_path), "environment design")
    if source_design.get("schema_version") != 1:
        raise EnvironmentError("environment design.schema_version must be 1")
    design, option_state = resolve_design_options(
        source_design, design_path, DEFAULT_OPTION_CATALOG
    )

    if args.framework_root:
        _mapping(design.setdefault("source", {}), "source")[
            "framework_root"
        ] = args.framework_root
    if args.output:
        _mapping(design.setdefault("destination", {}), "destination")[
            "root"
        ] = args.output
    if args.model:
        _mapping(design.setdefault("model", {}), "model")[
            "name"
        ] = args.model
    if args.profile:
        _mapping(design.setdefault("model", {}), "model")[
            "profile"
        ] = args.profile
    validate_design_selections(design, option_state)

    source_cfg = _mapping(design.get("source"), "source")
    source_text = source_cfg.get("framework_root")
    if not source_text:
        raise EnvironmentError("source.framework_root is required")
    framework_root = _expand_path(
        str(source_text), design_path.parent, "source.framework_root"
    )
    if not framework_root.is_dir():
        raise EnvironmentError(f"framework root not found: {framework_root}")

    model_cfg = _mapping(design.get("model"), "model")
    model = str(model_cfg.get("name") or "").lower()
    profile = str(model_cfg.get("profile") or "")
    include_tests = bool(model_cfg.get("include_tests", False))
    if not model or not profile:
        raise EnvironmentError("model.name and model.profile are required")

    catalog_relative = str(
        source_cfg.get(
            "model_catalog", "ScriptLibrary/BuildSolver/model_catalog.yaml"
        )
    )
    catalog_path = _safe_source(
        framework_root, catalog_relative, "source.model_catalog"
    )
    catalog = _mapping(load_yaml(catalog_path), "model catalog")
    models = _mapping(catalog.get("models"), "model catalog.models")
    if model not in models:
        raise EnvironmentError(f"unknown model {model!r}; available: {sorted(models)}")
    catalog_entry = _mapping(models[model], f"model catalog.{model}")
    solver_subpath = str(catalog_entry.get("library_subpath") or "")
    solver_root = (framework_root / "SolverLibrary" / solver_subpath).resolve()
    if not solver_root.is_dir():
        raise EnvironmentError(f"solver root not found: {solver_root}")
    manifest_name = str(catalog_entry.get("manifest") or "solver_manifest.yaml")
    manifest_path = _safe_solver_file(solver_root, manifest_name)
    manifest = _mapping(load_yaml(manifest_path), "solver manifest")
    selected_files, components = _selected_solver_files(
        solver_root, manifest, profile, include_tests
    )
    dependencies = _inspect_dependencies(solver_root, manifest, selected_files)
    use_mpi, use_openmp, use_cuda = _parallel_features(
        design, manifest, profile
    )

    destination_cfg = _mapping(design.get("destination"), "destination")
    destination_text = destination_cfg.get("root")
    if not destination_text:
        raise EnvironmentError("destination.root is required")
    destination_root = _expand_path(
        str(destination_text), design_path.parent, "destination.root"
    )

    case_cfg = _mapping(design.get("case"), "case")
    auto_case_number = bool(destination_cfg.get("auto_case_number", False))
    if args.output:
        auto_case_number = False
    if auto_case_number:
        width = int(destination_cfg.get("case_number_width", 4))
        start = int(destination_cfg.get("case_number_start", 1))
        if args.overwrite:
            if not args.case_id:
                raise EnvironmentError(
                    "--overwrite with automatic case numbering requires "
                    "--case-id (for example --case-id case0015), or use "
                    "--output with the complete destination path"
                )
            output, allocated_case_id = _requested_case_destination(
                destination_root, model, args.case_id, width, start
            )
            if not output.is_dir():
                raise EnvironmentError(
                    f"overwrite target does not exist: {output}"
                )
        else:
            if args.case_id:
                raise EnvironmentError(
                    "--case-id selects an existing automatically numbered "
                    "environment and must be used with --overwrite"
                )
            output, allocated_case_id = _numbered_case_destination(
                destination_root, model, width, start
            )
        case_cfg["id"] = allocated_case_id
    else:
        if args.case_id:
            raise EnvironmentError(
                "--case-id is only valid with automatic case numbering; "
                "use --output to select an explicit destination path"
            )
        output = destination_root
    global_case_index = _global_case_index_path(
        destination_cfg, design_path, output
    )

    target_cfg = _mapping(design.get("target"), "target")
    machine_relative = str(target_cfg.get("machine_profile") or "")
    if not machine_relative:
        raise EnvironmentError("target.machine_profile is required")
    machine_path = _safe_source(
        framework_root, machine_relative, "target.machine_profile"
    )
    template_relative = str(case_cfg.get("template") or "")
    if bool(case_cfg.get("create", True)) and not template_relative:
        raise EnvironmentError("case.template is required when case.create is true")

    execution = _mapping(design.get("execution", {}), "execution")
    parallel_jobs = int(execution.get("parallel_jobs", 8))
    if parallel_jobs < 1:
        raise EnvironmentError("execution.parallel_jobs must be positive")

    if args.dry_run:
        print("[DRY-RUN] Execution environment plan")
        if option_state["selections"]:
            selected_text = ", ".join(
                f"{category}={record['id']}"
                for category, record in option_state["selections"].items()
            )
            print(f"  selections:  {selected_text}")
        print(f"  source:      {framework_root}")
        print(f"  destination: {output}")
        print(f"  case index:  {global_case_index}")
        print(f"  model:       {model}")
        print(f"  profile:     {profile}")
        print(
            "  parallel:    "
            f"MPI={use_mpi}, OpenMP={use_openmp}, CUDA={use_cuda}"
        )
        print(f"  components:  {', '.join(components)}")
        print(f"  solver files:{len(selected_files)}")
        return output

    if output.exists():
        if not args.overwrite:
            raise EnvironmentError(f"destination already exists; use --overwrite: {output}")
        if not (output / GENERATED_MARKER).is_file():
            raise EnvironmentError(
                f"refusing to replace directory without {GENERATED_MARKER}: {output}"
            )
        shutil.rmtree(output)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.tmp-{uuid.uuid4().hex}"
    temporary.mkdir(parents=True)
    source_records: list[dict[str, Any]] = []
    try:
        solver_destination = temporary / "SolverLibrary" / solver_subpath
        for relative in selected_files:
            source = _safe_solver_file(solver_root, relative)
            destination = solver_destination / relative
            _copy_file(
                source, destination, source_records, framework_root, "solver"
            )
            if relative in dependencies:
                source_records[-1].update(dependencies[relative])
        _copy_file(
            manifest_path,
            solver_destination / manifest_name,
            source_records,
            framework_root,
            "solver_manifest",
        )

        build_tools = framework_root / "ScriptLibrary" / "BuildSolver"
        for source in _tree_files(build_tools):
            relative = source.relative_to(build_tools)
            _copy_file(
                source,
                temporary / "ScriptLibrary" / "BuildSolver" / relative,
                source_records,
                framework_root,
                "build_tool",
            )
        setup_generator = _safe_source(
            framework_root,
            "ScriptLibrary/SetupCase/create_case_from_template.py",
            "SetupCase generator",
        )
        _copy_file(
            setup_generator,
            temporary
            / "ScriptLibrary"
            / "SetupCase"
            / "create_case_from_template.py",
            source_records,
            framework_root,
            "case_tool",
        )
        for name in (
            "run_case.py",
            "postprocess_case.py",
            "case_input.py",
            "global_case_index.py",
            "yaml_support.py",
        ):
            source = SCRIPT_DIR / name
            _copy_file(
                source,
                temporary / "tools" / name,
                source_records,
                framework_root,
                "runtime_tool",
            )
        _copy_file(
            machine_path,
            temporary / "config" / "machine.yaml",
            source_records,
            framework_root,
            "machine_profile",
        )
        if bool(case_cfg.get("create", True)):
            template_path = _safe_source(
                framework_root, template_relative, "case.template"
            )
            _copy_file(
                template_path,
                temporary / "templates" / "case_template.yaml",
                source_records,
                framework_root,
                "case_template",
            )

        local_design = _local_build_design(
            model,
            profile,
            str(execution.get("configuration", "Release")),
            include_tests,
            parallel_jobs,
        )
        _write_json(
            temporary
            / "ScriptLibrary"
            / "BuildSolver"
            / "build.local.yaml",
            local_design,
        )
        shutil.copy2(design_path, temporary / "environment.source.yaml")
        _write_json(temporary / "environment.resolved.yaml", design)
        local_manifest = solver_destination / manifest_name
        effective_case_cfg = dict(case_cfg)
        case_id, input_name = _create_case(
            temporary,
            design_path,
            effective_case_cfg,
            model,
            profile,
            local_manifest,
            use_mpi,
            use_openmp,
            use_cuda,
        )

        scheduler = _mapping(design.get("scheduler", {}), "scheduler")
        if bool(scheduler.get("enabled", False)):
            (temporary / "submit.slurm").write_text(
                _slurm_script(
                    scheduler, include_tests, use_mpi, use_openmp
                ),
                encoding="utf-8",
                newline="\n",
            )

        generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        lock = {
            "schema_version": 1,
            "generated_at_utc": generated_at,
            "selections": {
                category: {
                    "id": record["id"],
                    "catalog_id": record["catalog_id"],
                }
                for category, record in option_state["selections"].items()
            },
            "model": model,
            "profile": profile,
            "include_tests": include_tests,
            "case_id": case_id,
            "case_directory": f"cases/{case_id}",
            "input_name": input_name,
            "case_index_path": _portable_path_reference(
                global_case_index, output
            ),
            "use_mpi": use_mpi,
            "use_openmp": use_openmp,
            "use_cuda": use_cuda,
        }
        provenance = {
            "schema_version": 1,
            "generated_at_utc": generated_at,
            "framework_root_at_generation": str(framework_root),
            "environment_design_sha256": _sha256(design_path),
            "option_catalogs": [
                {
                    **record,
                    "sha256": _sha256(Path(record["path"])),
                }
                for record in option_state["catalogs"]
            ],
            "selections": {
                category: {
                    "id": record["id"],
                    "catalog_id": record["catalog_id"],
                }
                for category, record in option_state["selections"].items()
            },
            "model_catalog_sha256": _sha256(catalog_path),
            "solver_manifest_sha256": _sha256(manifest_path),
            "solver_git": _git_state(solver_root),
            "model": model,
            "profile": profile,
            "components": components,
            "files": source_records,
        }
        _write_json(temporary / "environment.lock.json", lock)
        _write_json(temporary / "provenance.json", provenance)
        _write_json(
            temporary / GENERATED_MARKER,
            {"schema_version": 1, "generated_at_utc": generated_at},
        )
        (temporary / "README.md").write_text(
            _generated_readme(
                model, profile, case_id, use_mpi, use_openmp
            ),
            encoding="utf-8",
        )
        temporary.replace(output)
        sync_environment_case(output, lock, global_case_index)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise

    archive_cfg = _mapping(design.get("archive", {}), "archive")
    archive = _create_archive(
        output, design_path, archive_cfg, args.archive, args.archive_format
    )
    print(f"[OK] Generated execution environment: {output}")
    print(f"     model/profile: {model}/{profile}")
    print(f"     case:          {case_id}")
    print(f"     case index:    {global_case_index}")
    print(f"     solver files:  {len(selected_files)}")
    if archive:
        print(f"[OK] Portable archive: {archive}")
    return output


def list_options(args: argparse.Namespace) -> None:
    design_path = (
        Path(args.design).resolve() if args.design else DEFAULT_DESIGN.resolve()
    )
    design = _mapping(load_yaml(design_path), "environment design")
    if design.get("schema_version") != 1:
        raise EnvironmentError("environment design.schema_version must be 1")
    catalog = load_option_catalogs(
        design, design_path, DEFAULT_OPTION_CATALOG
    )
    category = None if args.list_options == "all" else args.list_options
    if category is not None and category not in OPTION_CATEGORIES:
        raise EnvironmentError(
            f"unknown option category {category!r}; "
            f"available: {list(OPTION_CATEGORIES)}"
        )
    if args.options_format == "json":
        print(catalog_as_json(catalog, category))
    else:
        print(catalog_as_text(catalog, category))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy selected framework sources into an external run environment."
    )
    parser.add_argument("design", nargs="?", help="Environment YAML")
    parser.add_argument("--framework-root", help="Override source.framework_root")
    parser.add_argument("--output", help="Override destination.root")
    parser.add_argument("--model", help="Override model.name")
    parser.add_argument("--profile")
    parser.add_argument(
        "--case-id",
        help="Existing automatic case number selected with --overwrite, for example case0015",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--archive", action="store_true")
    parser.add_argument("--archive-format", choices=["zip", "gztar"])
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--list-options",
        nargs="?",
        const="all",
        metavar="CATEGORY",
        help="List all choices, or one category, without generating an environment",
    )
    parser.add_argument(
        "--options-format",
        choices=["text", "json"],
        default="text",
        help="Output format used with --list-options",
    )
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        if args.list_options is not None:
            list_options(args)
        else:
            prepare(args)
        return 0
    except (
        CaseInputError,
        EnvironmentError,
        GlobalCaseIndexError,
        OptionCatalogError,
        OSError,
        ValueError,
        YamlFormatError,
    ) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
