#!/usr/bin/env python3
"""Validate an NSE YAML build design, generate CMake cache, build, and run.

Generated files and compiler working directories are placed outside SolverLibrary.
This keeps the NAS copy read-mostly and avoids cmd.exe using a UNC working directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from yaml_support import YamlFormatError, load_yaml


SCRIPT_PATH = Path(__file__).resolve()
SOLVER_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_DESIGN = SOLVER_ROOT / "config" / "build.yaml"
FORTRAN_SUFFIXES = {".f90", ".f95", ".f03", ".f08"}
MODULE_PATTERN = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)
USE_PATTERN = re.compile(
    r"^\s*use(?:\s*,\s*[^:]*)?\s*(?:::\s*)?([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)


class BuildDesignError(RuntimeError):
    """Raised when the YAML design is incomplete or internally inconsistent."""


@dataclass(frozen=True)
class ResolvedBuild:
    design_path: Path
    profile_path: Path
    catalog_path: Path
    design: dict[str, Any]
    profile: dict[str, Any]
    catalog: dict[str, Any]
    configuration: str
    build_dir: Path
    executable_name: str
    library_sources: tuple[str, ...]
    main_source: str
    modules: tuple[str, ...]
    features: tuple[str, ...]
    cmake_variables: dict[str, Any]


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise BuildDesignError(f"{label} must be a YAML mapping")
    return value


def _list(value: Any, label: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise BuildDesignError(f"{label} must be a YAML list")
    return value


def _schema_one(document: dict[str, Any], label: str) -> None:
    if document.get("schema_version") != 1:
        raise BuildDesignError(f"{label}.schema_version must be 1")


def _relative_file(base: Path, text: str, label: str) -> Path:
    candidate = Path(os.path.expandvars(os.path.expanduser(text)))
    if not candidate.is_absolute():
        candidate = base / candidate
    candidate = candidate.resolve()
    if not candidate.is_file():
        raise BuildDesignError(f"{label} not found: {candidate}")
    return candidate


def _source_file(relative: str) -> Path:
    root = SOLVER_ROOT.resolve()
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise BuildDesignError(f"module path escapes the solver root: {relative}") from exc
    if not candidate.is_file():
        raise BuildDesignError(f"module source not found: {relative}")
    return candidate


def _expand_output_path(text: str, base: Path) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(text))
    unresolved = re.search(r"\$\{[^}]+\}|%[^%]+%", expanded)
    if unresolved:
        raise BuildDesignError(f"undefined environment variable in output path: {unresolved.group(0)}")
    path = Path(expanded)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _safe_name(value: Any, label: str) -> str:
    text = str(value or "").strip()
    if not text or not re.fullmatch(r"[A-Za-z0-9_.-]+", text):
        raise BuildDesignError(f"{label} must contain only letters, numbers, '.', '_' or '-'")
    return text


def _resolve_modules(
    catalog: dict[str, Any], module_set_name: str, includes: list[str], excludes: list[str]
) -> tuple[list[str], dict[str, dict[str, Any]]]:
    modules_raw = _mapping(catalog.get("modules"), "module catalog.modules")
    modules = {
        str(name): _mapping(record, f"module {name}") for name, record in modules_raw.items()
    }
    sets = _mapping(catalog.get("module_sets"), "module catalog.module_sets")
    if module_set_name not in sets:
        raise BuildDesignError(
            f"unknown module set {module_set_name!r}; available: {sorted(sets)}"
        )
    module_set = _mapping(sets[module_set_name], f"module set {module_set_name}")
    roots = [str(value) for value in _list(module_set.get("entrypoints"), "entrypoints")]
    roots.extend(includes)
    excluded = set(excludes)
    state: dict[str, int] = {}
    ordered: list[str] = []

    def visit(name: str, chain: tuple[str, ...]) -> None:
        if name in excluded:
            parent = chain[-1] if chain else "module set"
            raise BuildDesignError(f"excluded module {name!r} is required by {parent!r}")
        if name not in modules:
            raise BuildDesignError(f"unknown module {name!r}")
        if state.get(name) == 2:
            return
        if state.get(name) == 1:
            raise BuildDesignError("module dependency cycle: " + " -> ".join((*chain, name)))
        state[name] = 1
        requires = _list(modules[name].get("requires"), f"module {name}.requires")
        for dependency in requires:
            visit(str(dependency), (*chain, name))
        state[name] = 2
        ordered.append(name)

    for root in roots:
        visit(root, ())
    return ordered, modules


def _inspect_fortran_sources(
    selected: list[str], modules: dict[str, dict[str, Any]]
) -> None:
    all_providers: dict[str, str] = {}
    for module_name, record in modules.items():
        relative = str(record.get("path", ""))
        path = _source_file(relative)
        if path.suffix.lower() not in FORTRAN_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for provided in MODULE_PATTERN.findall(text):
            all_providers[provided.lower()] = module_name

    selected_providers: dict[str, str] = {}
    selected_uses: list[tuple[str, str]] = []
    for module_name in selected:
        relative = str(modules[module_name].get("path", ""))
        path = _source_file(relative)
        if path.suffix.lower() not in FORTRAN_SUFFIXES:
            raise BuildDesignError(f"unsupported NSE source extension: {relative}")
        text = path.read_text(encoding="utf-8", errors="replace")
        for provided in MODULE_PATTERN.findall(text):
            key = provided.lower()
            if key in selected_providers and selected_providers[key] != relative:
                raise BuildDesignError(
                    f"Fortran module {key} is defined by both "
                    f"{selected_providers[key]} and {relative}"
                )
            selected_providers[key] = relative
        selected_uses.extend((relative, used.lower()) for used in USE_PATTERN.findall(text))

    missing = [
        f"{relative} uses {used}"
        for relative, used in selected_uses
        if used in all_providers and used not in selected_providers
    ]
    if missing:
        raise BuildDesignError("selected modules have missing internal dependencies:\n  - " + "\n  - ".join(missing))


def _platform_name() -> str:
    if os.name == "nt":
        return "windows"
    if sys.platform.startswith("linux"):
        return "linux"
    return sys.platform.lower()


def _cmake_path(path: Path | str) -> str:
    return str(path).replace("\\", "/")


def _resolve_build(args: argparse.Namespace) -> ResolvedBuild:
    design_path = Path(args.design or DEFAULT_DESIGN).resolve()
    if not design_path.is_file():
        raise BuildDesignError(f"build design not found: {design_path}")
    design = _mapping(load_yaml(design_path), "build design")
    _schema_one(design, "build design")

    profile_text = args.profile or design.get("machine_profile")
    if not profile_text:
        raise BuildDesignError("machine_profile is required")
    profile_path = _relative_file(design_path.parent, str(profile_text), "machine profile")
    profile = _mapping(load_yaml(profile_path), "machine profile")
    _schema_one(profile, "machine profile")
    expected_platform = str(profile.get("platform", "")).lower()
    if expected_platform and expected_platform != _platform_name():
        raise BuildDesignError(
            f"profile platform is {expected_platform}, but this machine is {_platform_name()}"
        )

    catalog_text = str(design.get("module_catalog") or "module_catalog.yaml")
    catalog_path = _relative_file(design_path.parent, catalog_text, "module catalog")
    catalog = _mapping(load_yaml(catalog_path), "module catalog")
    _schema_one(catalog, "module catalog")

    solver = _mapping(design.get("solver"), "build design.solver")
    module_set = str(solver.get("module_set") or "default")
    include_modules = [
        str(value) for value in _list(solver.get("include_modules"), "solver.include_modules")
    ]
    exclude_modules = [
        str(value) for value in _list(solver.get("exclude_modules"), "solver.exclude_modules")
    ]
    selected, modules = _resolve_modules(
        catalog, module_set, include_modules, exclude_modules
    )
    _inspect_fortran_sources(selected, modules)

    executable_modules = [
        name for name in selected if str(modules[name].get("role", "library")) == "executable"
    ]
    if len(executable_modules) != 1:
        raise BuildDesignError(
            f"exactly one executable module is required; selected: {executable_modules}"
        )
    main_module = executable_modules[0]
    library_sources = tuple(
        str(modules[name]["path"])
        for name in selected
        if str(modules[name].get("role", "library")) == "library"
    )
    main_source = str(modules[main_module]["path"])
    features = sorted(
        {
            str(feature).lower()
            for name in selected
            for feature in _list(modules[name].get("features"), f"module {name}.features")
        }
    )

    mpi = _mapping(profile.get("mpi"), "machine profile.mpi")
    openmp = _mapping(profile.get("openmp"), "machine profile.openmp")
    if "mpi" in features and not bool(mpi.get("enabled", False)):
        raise BuildDesignError("selected NSE modules require MPI, but mpi.enabled is false")
    if "openmp" in features and not bool(openmp.get("enabled", False)):
        raise BuildDesignError("selected NSE modules require OpenMP, but openmp.enabled is false")

    build = _mapping(design.get("build"), "build design.build")
    configuration = str(args.configuration or build.get("configuration") or "Release")
    if configuration not in {"Debug", "Release", "RelWithDebInfo"}:
        raise BuildDesignError("build.configuration must be Debug, Release, or RelWithDebInfo")
    profile_id = _safe_name(profile.get("profile_id"), "profile_id")
    design_id = _safe_name(design.get("design_id"), "design_id")
    executable_name = _safe_name(solver.get("executable_name") or "solver", "executable_name")

    if args.build_dir:
        build_dir = Path(args.build_dir).resolve()
    else:
        output_root = _expand_output_path(
            str(build.get("output_root") or "build/generated"), design_path.parent
        )
        directory_template = str(build.get("directory_name") or "{profile_id}-{configuration}")
        try:
            directory_name = directory_template.format(
                profile_id=profile_id,
                configuration=configuration.lower(),
                design_id=design_id,
            )
        except KeyError as exc:
            raise BuildDesignError(f"unknown placeholder in build.directory_name: {exc}") from exc
        directory_name = _safe_name(directory_name, "build.directory_name")
        build_dir = (output_root / directory_name).resolve()
    if os.name == "nt" and str(build_dir).startswith("\\\\"):
        raise BuildDesignError(
            "build directory must be local or a mapped drive on Windows; UNC build directories "
            "cause cmd.exe failures"
        )

    compiler = _mapping(profile.get("compiler"), "machine profile.compiler")
    compiler_command = str(compiler.get("fortran") or "gfortran")
    cmake_variables: dict[str, Any] = {
        "CMAKE_Fortran_COMPILER": compiler_command,
        "CMAKE_BUILD_TYPE": configuration,
        "NSE_CORE_SOURCES": list(library_sources),
        "NSE_MAIN_SOURCE": main_source,
        "NSE_OUTPUT_NAME": executable_name,
        "NSE_MPI_PROVIDER": str(mpi.get("provider") or "AUTO").upper(),
        "NSE_ENABLE_OPENMP": bool(openmp.get("enabled", False)),
        "NSE_ENABLE_WARNINGS": bool(build.get("warnings", True)),
        "NSE_COMMON_FLAGS": [str(value) for value in _list(compiler.get("common_flags"), "compiler.common_flags")],
        "NSE_WARNING_FLAGS": [str(value) for value in _list(compiler.get("warning_flags"), "compiler.warning_flags")],
        "NSE_DEBUG_FLAGS": [str(value) for value in _list(compiler.get("debug_flags"), "compiler.debug_flags")],
        "NSE_RELEASE_FLAGS": [str(value) for value in _list(compiler.get("release_flags"), "compiler.release_flags")],
        "NSE_RELWITHDEBINFO_FLAGS": [
            str(value)
            for value in _list(compiler.get("relwithdebinfo_flags"), "compiler.relwithdebinfo_flags")
        ],
        "NSE_INIT_FFT_BACKEND": (
            "2decomp_fftw" if "distributed_fft" in features else "none"
        ),
    }
    if mpi.get("root"):
        cmake_variables["MSMPI_ROOT"] = str(mpi["root"])
    libraries = _mapping(profile.get("libraries") or {}, "machine profile.libraries")
    if libraries.get("decomp2d_root"):
        cmake_variables["NSE_2DECOMP_ROOT"] = str(libraries["decomp2d_root"])

    return ResolvedBuild(
        design_path=design_path,
        profile_path=profile_path,
        catalog_path=catalog_path,
        design=design,
        profile=profile,
        catalog=catalog,
        configuration=configuration,
        build_dir=build_dir,
        executable_name=executable_name,
        library_sources=library_sources,
        main_source=main_source,
        modules=tuple(selected),
        features=tuple(features),
        cmake_variables=cmake_variables,
    )


def _cmake_bracket(value: str) -> str:
    marker = "="
    while f"]{marker}]" in value:
        marker += "="
    return f"[{marker}[{value}]{marker}]"


def _cache_value(value: Any) -> tuple[str, str]:
    if isinstance(value, bool):
        return ("BOOL", "ON" if value else "OFF")
    if isinstance(value, list):
        return ("STRING", ";".join(str(item) for item in value))
    return ("STRING", str(value))


def _initial_cache_text(resolved: ResolvedBuild) -> str:
    lines = [
        "# Generated by tools/build_from_yaml.py. Do not edit.",
        f"# Design: {_cmake_path(resolved.design_path)}",
        "",
    ]
    path_variables = {"CMAKE_Fortran_COMPILER", "MSMPI_ROOT"}
    for name, raw_value in resolved.cmake_variables.items():
        cache_type, value = _cache_value(raw_value)
        if name in path_variables:
            cache_type = "FILEPATH" if name == "CMAKE_Fortran_COMPILER" else "PATH"
            value = value.replace("\\", "/")
        lines.append(
            f"set({name} {_cmake_bracket(value)} CACHE {cache_type} \"Generated from YAML\" FORCE)"
        )
    return "\n".join(lines) + "\n"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolved_json(resolved: ResolvedBuild) -> dict[str, Any]:
    build = _mapping(resolved.design["build"], "build")
    run = _mapping(resolved.design.get("run", {}), "run")
    return {
        "schema_version": 1,
        "design_id": resolved.design.get("design_id"),
        "design": _cmake_path(resolved.design_path),
        "design_sha256": _sha256(resolved.design_path),
        "profile": _cmake_path(resolved.profile_path),
        "profile_sha256": _sha256(resolved.profile_path),
        "module_catalog": _cmake_path(resolved.catalog_path),
        "module_catalog_sha256": _sha256(resolved.catalog_path),
        "solver_root": _cmake_path(SOLVER_ROOT),
        "build_directory": _cmake_path(resolved.build_dir),
        "configuration": resolved.configuration,
        "executable_name": resolved.executable_name,
        "selected_modules": list(resolved.modules),
        "required_features": list(resolved.features),
        "library_sources": list(resolved.library_sources),
        "main_source": resolved.main_source,
        "cmake_variables": resolved.cmake_variables,
        "parallel_jobs": int(build.get("parallel_jobs", 1)),
        "run": run,
    }


def _write_generated_files(resolved: ResolvedBuild) -> tuple[Path, Path]:
    generated = resolved.build_dir / "generated"
    generated.mkdir(parents=True, exist_ok=True)
    cache_path = generated / "NSEInitialCache.cmake"
    resolved_path = generated / "resolved_build.json"
    cache_path.write_text(_initial_cache_text(resolved), encoding="utf-8")
    resolved_path.write_text(
        json.dumps(_resolved_json(resolved), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return cache_path, resolved_path


def _display_command(command: list[str]) -> str:
    return subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)


def _run_command(command: list[str], cwd: Path, env: dict[str, str] | None, dry_run: bool) -> None:
    print(f"[CMD] {_display_command(command)}", flush=True)
    print(f"      cwd={cwd}", flush=True)
    if dry_run:
        return
    result = subprocess.run(command, cwd=str(cwd), env=env, check=False)
    if result.returncode != 0:
        raise BuildDesignError(
            f"command failed with exit code {result.returncode}: {_display_command(command)}"
        )


def _tool(profile: dict[str, Any], name: str, default: str) -> str:
    tools = _mapping(profile.get("tools", {}), "machine profile.tools")
    return str(tools.get(name) or default)


def _configure_command(resolved: ResolvedBuild, cache_path: Path, fresh: bool) -> list[str]:
    command = [_tool(resolved.profile, "cmake", "cmake")]
    if fresh:
        command.append("--fresh")
    command.extend(
        [
            "-S",
            str(SOLVER_ROOT),
            "-B",
            str(resolved.build_dir),
            "-G",
            str(resolved.profile.get("generator") or "Ninja"),
            "-C",
            str(cache_path),
        ]
    )
    return command


def _build_command(resolved: ResolvedBuild, clean_first: bool) -> list[str]:
    build = _mapping(resolved.design["build"], "build")
    command = [
        _tool(resolved.profile, "cmake", "cmake"),
        "--build",
        str(resolved.build_dir),
        "--config",
        resolved.configuration,
        "--parallel",
        str(int(build.get("parallel_jobs", 1))),
    ]
    if clean_first:
        command.append("--clean-first")
    return command


def _resolve_input(args: argparse.Namespace, resolved: ResolvedBuild) -> Path:
    run = _mapping(resolved.design.get("run", {}), "run")
    if args.input_file:
        path = Path(args.input_file)
        if not path.is_absolute():
            path = Path.cwd() / path
    elif run.get("input_file"):
        path = resolved.design_path.parent / str(run["input_file"])
    else:
        raise BuildDesignError(
            "run requires input.dat; set run.input_file in build.yaml or pass --input-file"
        )
    path = path.resolve()
    if not path.is_file():
        raise BuildDesignError(f"run input file not found: {path}")
    return path


def _run_solver(args: argparse.Namespace, resolved: ResolvedBuild) -> None:
    suffix = ".exe" if os.name == "nt" else ""
    executable = resolved.build_dir / "bin" / f"{resolved.executable_name}{suffix}"
    if not executable.is_file():
        raise BuildDesignError(f"solver executable not found; build first: {executable}")
    input_path = _resolve_input(args, resolved)
    run_dir = resolved.build_dir / "run"
    if not args.dry_run:
        run_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, run_dir / "input.dat")

    run = _mapping(resolved.design.get("run", {}), "run")
    mpi = _mapping(resolved.profile.get("mpi", {}), "machine profile.mpi")
    processes = int(args.processes or run.get("mpi_processes", 1))
    if processes < 1:
        raise BuildDesignError("MPI process count must be at least 1")
    arguments = [str(value) for value in _list(run.get("program_arguments"), "run.program_arguments")]
    command: list[str] = []
    if bool(mpi.get("enabled", False)):
        command.extend(
            [
                _tool(resolved.profile, "mpi_launcher", "mpiexec"),
                _tool(resolved.profile, "mpi_process_option", "-n"),
                str(processes),
            ]
        )
    command.extend([str(executable), *arguments])
    environment = os.environ.copy()
    omp_threads = int(args.omp_threads or run.get("omp_threads", 1))
    if omp_threads < 1:
        raise BuildDesignError("OpenMP thread count must be at least 1")
    environment["OMP_NUM_THREADS"] = str(omp_threads)
    _run_command(command, run_dir, environment, args.dry_run)


def _summary(resolved: ResolvedBuild) -> None:
    print("[OK] NSE build design is valid")
    print(f"     design:        {resolved.design_path}")
    print(f"     profile:       {resolved.profile.get('profile_id')}")
    print(f"     configuration: {resolved.configuration}")
    print(f"     build dir:     {resolved.build_dir}")
    print(f"     modules:       {len(resolved.modules)}")
    print(f"     features:      {', '.join(resolved.features)}", flush=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate and execute an NSE CMake build from a YAML design."
    )
    parser.add_argument("design", nargs="?", help="Build YAML; default: config/build.yaml")
    parser.add_argument("--profile", help="Override machine_profile with another YAML file")
    parser.add_argument("--configuration", choices=["Debug", "Release", "RelWithDebInfo"])
    parser.add_argument("--build-dir", help="Override the generated local build directory")
    parser.add_argument("--input-file", help="input.dat source used by --run")
    parser.add_argument("--processes", type=int, help="MPI process count used by --run")
    parser.add_argument("--omp-threads", type=int, help="OpenMP thread count used by --run")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--configure", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--all", action="store_true", help="Configure, build, and run")
    parser.add_argument("--clean-first", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    fresh_group = parser.add_mutually_exclusive_group()
    fresh_group.add_argument("--fresh", dest="fresh", action="store_true")
    fresh_group.add_argument("--no-fresh", dest="fresh", action="store_false")
    parser.set_defaults(fresh=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.validate_only and (
            args.generate_only or args.configure or args.build or args.run or args.all
        ):
            raise BuildDesignError("--validate-only cannot be combined with build stages")
        if args.generate_only and (args.configure or args.build or args.run or args.all):
            raise BuildDesignError("--generate-only cannot be combined with build stages")

        resolved = _resolve_build(args)
        _summary(resolved)
        if args.validate_only:
            return 0

        explicit_stages = args.configure or args.build or args.run or args.all
        do_build = bool(args.build or args.all)
        do_run = bool(args.run or args.all)
        do_configure = bool(args.configure or do_build)
        if not explicit_stages and not args.generate_only:
            do_configure = True
            do_build = True
            run = _mapping(resolved.design.get("run", {}), "run")
            do_run = bool(run.get("enabled", False))

        build = _mapping(resolved.design["build"], "build")
        fresh = bool(build.get("configure_fresh", True)) if args.fresh is None else args.fresh
        cache_path = resolved.build_dir / "generated" / "NSEInitialCache.cmake"
        if args.dry_run:
            print(f"[DRY-RUN] generate {cache_path}")
        else:
            cache_path, resolved_path = _write_generated_files(resolved)
            print(f"[OK] Generated CMake cache: {cache_path}")
            print(f"[OK] Resolved build plan:  {resolved_path}", flush=True)

        if args.generate_only:
            return 0
        if do_configure:
            if not args.dry_run:
                resolved.build_dir.mkdir(parents=True, exist_ok=True)
            _run_command(
                _configure_command(resolved, cache_path, fresh),
                resolved.build_dir,
                None,
                args.dry_run,
            )
        if do_build:
            _run_command(
                _build_command(resolved, args.clean_first),
                resolved.build_dir,
                None,
                args.dry_run,
            )
        if do_run:
            _run_solver(args, resolved)
        return 0
    except (BuildDesignError, YamlFormatError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
