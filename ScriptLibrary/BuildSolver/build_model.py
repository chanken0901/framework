#!/usr/bin/env python3
"""Build and run an NSE or GPE solver selected by one YAML design.

The controller reads a model catalog and each solver's public manifest. CMake
cache files and all compiler outputs are generated in a local build directory,
so SolverLibrary may remain on a read-mostly NAS share.
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
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from yaml_support import YamlFormatError, load_yaml


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DESIGN = SCRIPT_DIR / "build.yaml"
FORTRAN_SUFFIXES = {".f90", ".f95", ".f03", ".f08"}
MODULE_PATTERN = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)
USE_PATTERN = re.compile(
    r"^\s*use(?:\s*,\s*[^:]*)?\s*(?:::\s*)?([a-z][a-z0-9_]*)",
    re.IGNORECASE | re.MULTILINE,
)
CPP_DIRECTIVE_PATTERN = re.compile(
    r"^\s*#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)$",
    re.IGNORECASE,
)
CPP_DEFINED_PATTERN = re.compile(
    r"\bdefined\s*(?:\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)|"
    r"([A-Za-z_][A-Za-z0-9_]*))"
)
CPP_IDENTIFIER_PATTERN = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")


class ModelBuildError(RuntimeError):
    """Raised when the unified build design cannot be resolved safely."""


def _cpp_condition(expression: str, defines: set[str], label: str) -> bool:
    def replace_defined(match: re.Match[str]) -> str:
        name = match.group(1) or match.group(2)
        return "1" if name in defines else "0"

    resolved = CPP_DEFINED_PATTERN.sub(replace_defined, expression.strip())
    resolved = CPP_IDENTIFIER_PATTERN.sub(
        lambda match: "1" if match.group(0) in defines else "0", resolved
    )
    resolved = resolved.replace("&&", " and ").replace("||", " or ")
    resolved = re.sub(r"!(?!=)", " not ", resolved)
    tokens = re.findall(r"\d+|and|or|not|\(|\)|\S+", resolved)
    if not tokens or any(
        token not in {"0", "1", "and", "or", "not", "(", ")"}
        for token in tokens
    ):
        raise ModelBuildError(
            f"unsupported Fortran preprocessor condition in {label}: {expression.strip()}"
        )
    return bool(eval(" ".join(tokens), {"__builtins__": {}}, {}))


def _active_fortran_source(
    text: str, defines: set[str], label: str
) -> str:
    active = True
    stack: list[tuple[bool, bool, bool]] = []
    output: list[str] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        match = CPP_DIRECTIVE_PATTERN.match(line)
        if match is None:
            if active:
                output.append(line)
            continue

        directive = match.group(1).lower()
        argument = match.group(2).strip()
        if directive in {"ifdef", "ifndef", "if"}:
            if directive == "ifdef":
                condition = argument in defines
            elif directive == "ifndef":
                condition = argument not in defines
            else:
                condition = _cpp_condition(argument, defines, label)
            stack.append((active, condition, False))
            active = active and condition
        elif directive == "elif":
            if not stack:
                raise ModelBuildError(
                    f"unmatched #elif in {label}:{line_number}"
                )
            parent_active, branch_taken, else_seen = stack[-1]
            if else_seen:
                raise ModelBuildError(f"#elif after #else in {label}:{line_number}")
            condition = False if branch_taken else _cpp_condition(
                argument, defines, label
            )
            stack[-1] = (parent_active, branch_taken or condition, False)
            active = parent_active and condition
        elif directive == "else":
            if not stack:
                raise ModelBuildError(
                    f"unmatched #else in {label}:{line_number}"
                )
            parent_active, branch_taken, else_seen = stack[-1]
            if else_seen:
                raise ModelBuildError(f"duplicate #else in {label}:{line_number}")
            condition = not branch_taken
            stack[-1] = (parent_active, True, True)
            active = parent_active and condition
        else:
            if not stack:
                raise ModelBuildError(
                    f"unmatched #endif in {label}:{line_number}"
                )
            parent_active, _, _ = stack.pop()
            active = parent_active

    if stack:
        raise ModelBuildError(f"unterminated preprocessor condition in {label}")
    return "\n".join(output)


@dataclass(frozen=True)
class ResolvedBuild:
    design_path: Path
    catalog_path: Path
    machine_path: Path
    manifest_path: Path
    solver_library_root: Path
    solver_root: Path
    build_dir: Path
    design: dict[str, Any]
    machine: dict[str, Any]
    manifest: dict[str, Any]
    model: str
    solver_profile: str
    adapter: str
    configuration: str
    executable: str
    selected_components: tuple[str, ...]
    selected_files: tuple[str, ...]
    cmake_variables: dict[str, Any]
    use_mpi: bool
    use_openmp: bool
    gpu_backend: str
    input_name: str
    pass_input_as_argument: bool


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ModelBuildError(f"{label} must be a YAML mapping")
    return value


def _list(value: Any, label: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ModelBuildError(f"{label} must be a YAML list")
    return value


def _schema_one(document: dict[str, Any], label: str) -> None:
    if document.get("schema_version") != 1:
        raise ModelBuildError(f"{label}.schema_version must be 1")


def _platform_name() -> str:
    if os.name == "nt":
        return "windows"
    if sys.platform.startswith("linux"):
        return "linux"
    return sys.platform.lower()


def _safe_name(value: Any, label: str) -> str:
    text = str(value or "").strip()
    if not text or not re.fullmatch(r"[A-Za-z0-9_.-]+", text):
        raise ModelBuildError(
            f"{label} must contain only letters, numbers, '.', '_' or '-'"
        )
    return text


def _expand_path(text: str, base: Path, label: str) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(text))
    unresolved = re.search(r"\$\{[^}]+\}|%[^%]+%", expanded)
    if unresolved:
        raise ModelBuildError(
            f"undefined environment variable in {label}: {unresolved.group(0)}"
        )
    path = Path(expanded)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _existing_file(base: Path, text: str, label: str) -> Path:
    path = _expand_path(text, base, label)
    if not path.is_file():
        raise ModelBuildError(f"{label} not found: {path}")
    return path


def _path_within(root: Path, relative: str, require_file: bool = True) -> Path:
    root_resolved = root.resolve()
    candidate = (root_resolved / relative).resolve()
    try:
        candidate.relative_to(root_resolved)
    except ValueError as exc:
        raise ModelBuildError(f"manifest path escapes solver root: {relative}") from exc
    if require_file and not candidate.is_file():
        raise ModelBuildError(f"manifest file does not exist: {relative}")
    return candidate


def _enabled(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().upper() in {"1", "ON", "TRUE", "YES"}


def _resolve_components(
    solver_root: Path,
    manifest: dict[str, Any],
    profile_name: str,
    include_tests: bool,
) -> tuple[list[str], list[str], list[str], list[str]]:
    profiles = _mapping(manifest.get("profiles"), "solver manifest.profiles")
    if profile_name not in profiles:
        raise ModelBuildError(
            f"unknown solver profile {profile_name!r}; available: {sorted(profiles)}"
        )
    profile = _mapping(profiles[profile_name], f"solver profile {profile_name}")
    component_names = [
        str(value)
        for value in _list(profile.get("components"), f"profile {profile_name}.components")
    ]
    if include_tests:
        component_names.extend(
            str(value)
            for value in _list(
                profile.get("test_components"), f"profile {profile_name}.test_components"
            )
        )

    components = _mapping(manifest.get("components"), "solver manifest.components")
    selected_files: list[str] = []
    library_sources: list[str] = []
    executable_sources: list[str] = []
    seen: set[str] = set()
    for component_name in component_names:
        if component_name not in components:
            raise ModelBuildError(
                f"profile {profile_name} references unknown component {component_name!r}"
            )
        component = _mapping(components[component_name], f"component {component_name}")
        files = _list(component.get("files"), f"component {component_name}.files")
        role = str(component.get("cmake_role") or "")
        for value in files:
            relative = str(value).replace("\\", "/")
            path = _path_within(solver_root, relative)
            if relative not in seen:
                seen.add(relative)
                selected_files.append(relative)
            if path.suffix.lower() in FORTRAN_SUFFIXES:
                if role == "library" and relative not in library_sources:
                    library_sources.append(relative)
                elif role == "executable" and relative not in executable_sources:
                    executable_sources.append(relative)

    cmake_file = str(manifest.get("cmake_file") or "CMakeLists.txt")
    _path_within(solver_root, cmake_file)
    return component_names, selected_files, library_sources, executable_sources


def _inspect_fortran_dependencies(
    solver_root: Path,
    manifest: dict[str, Any],
    selected_files: list[str],
    preprocessor_defines: set[str] | None = None,
) -> None:
    defines = preprocessor_defines or set()
    all_provider_names: set[str] = set()
    components = _mapping(manifest.get("components"), "solver manifest.components")
    for component in components.values():
        record = _mapping(component, "component")
        for value in _list(record.get("files"), "component.files"):
            relative = str(value).replace("\\", "/")
            path = _path_within(solver_root, relative, require_file=False)
            if not path.is_file():
                # Generated run environments intentionally contain only the
                # source files selected by one solver profile.
                continue
            if path.suffix.lower() not in FORTRAN_SUFFIXES:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            all_provider_names.update(name.lower() for name in MODULE_PATTERN.findall(text))

    providers: dict[str, str] = {}
    uses: list[tuple[str, str]] = []
    for relative in selected_files:
        path = _path_within(solver_root, relative)
        if path.suffix.lower() not in FORTRAN_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        active_text = _active_fortran_source(text, defines, relative)
        for provided in MODULE_PATTERN.findall(active_text):
            key = provided.lower()
            if key in providers and providers[key] != relative:
                raise ModelBuildError(
                    f"selected sources define Fortran module {key} twice: "
                    f"{providers[key]} and {relative}"
                )
            providers[key] = relative
        uses.extend(
            (relative, used.lower()) for used in USE_PATTERN.findall(active_text)
        )

    missing = [
        f"{relative} uses {used}"
        for relative, used in uses
        if used in all_provider_names and used not in providers
    ]
    if missing:
        raise ModelBuildError(
            "selected solver profile has missing internal modules:\n  - "
            + "\n  - ".join(missing)
        )


def _model_machine_variables(
    adapter: str,
    machine: dict[str, Any],
    build: dict[str, Any],
    profile: dict[str, Any],
    library_sources: list[str],
    executable_sources: list[str],
    executable: str,
    build_dir: Path,
    cmake_variables: dict[str, Any],
) -> tuple[bool, bool, str]:
    compiler = _mapping(machine.get("compiler"), "machine profile.compiler")
    mpi = _mapping(machine.get("mpi", {}), "machine profile.mpi")
    openmp = _mapping(machine.get("openmp", {}), "machine profile.openmp")
    libraries = _mapping(machine.get("libraries", {}), "machine profile.libraries")
    execution = _mapping(profile.get("execution", {}), "solver profile.execution")

    common_flags = [
        str(value) for value in _list(compiler.get("common_flags"), "compiler.common_flags")
    ]
    warning_flags = [
        str(value) for value in _list(compiler.get("warning_flags"), "compiler.warning_flags")
    ]
    debug_flags = [
        str(value) for value in _list(compiler.get("debug_flags"), "compiler.debug_flags")
    ]
    release_flags = [
        str(value) for value in _list(compiler.get("release_flags"), "compiler.release_flags")
    ]
    relwithdebinfo_flags = [
        str(value)
        for value in _list(
            compiler.get("relwithdebinfo_flags"), "compiler.relwithdebinfo_flags"
        )
    ]
    warnings = bool(build.get("warnings", True))
    cmake_variables["CMAKE_RUNTIME_OUTPUT_DIRECTORY"] = str(build_dir / "bin")

    if adapter == "nse_cmake_v1":
        if not library_sources or len(executable_sources) != 1:
            raise ModelBuildError(
                "nse_cmake_v1 requires library components and one executable source"
            )
        use_mpi = bool(execution.get("use_mpi", True))
        use_openmp = bool(execution.get("use_openmp", True))
        gpu_backend = str(cmake_variables.get("NSE_GPU_BACKEND", "none")).lower()
        if use_openmp and not bool(openmp.get("enabled", False)):
            raise ModelBuildError("NSE profile requires OpenMP, but machine openmp.enabled is false")
        if gpu_backend not in {"none", "cuda"}:
            raise ModelBuildError(
                f"unsupported NSE GPU backend {gpu_backend!r}; use none or cuda"
            )
        # NSE supports both cuda_single and MPI+CUDA profiles.  The profile's
        # execution metadata and selected main program determine whether MPI
        # is enabled; GPU_BACKEND=cuda by itself does not imply single GPU.
        cmake_variables.update(
            {
                "NSE_CORE_SOURCES": library_sources,
                "NSE_MAIN_SOURCE": executable_sources[0],
                "NSE_OUTPUT_NAME": executable,
                "NSE_USE_MPI": use_mpi,
                "NSE_MPI_PROVIDER": str(mpi.get("provider") or "AUTO").upper(),
                "NSE_ENABLE_OPENMP": use_openmp,
                "NSE_ENABLE_WARNINGS": warnings,
                "NSE_COMMON_FLAGS": common_flags,
                "NSE_WARNING_FLAGS": warning_flags,
                "NSE_DEBUG_FLAGS": debug_flags,
                "NSE_RELEASE_FLAGS": release_flags,
                "NSE_RELWITHDEBINFO_FLAGS": relwithdebinfo_flags,
            }
        )
        if use_mpi and mpi.get("root"):
            cmake_variables["MSMPI_ROOT"] = str(mpi["root"])
        if libraries.get("decomp2d_root"):
            cmake_variables["NSE_2DECOMP_ROOT"] = str(
                libraries["decomp2d_root"]
            )
        if gpu_backend == "cuda":
            if libraries.get("cuda_compiler"):
                cmake_variables["CMAKE_CUDA_COMPILER"] = str(
                    libraries["cuda_compiler"]
                )
            if libraries.get("cuda_toolkit_root"):
                cmake_variables["CUDAToolkit_ROOT"] = str(
                    libraries["cuda_toolkit_root"]
                )
            if libraries.get("cuda_architectures") is not None:
                cmake_variables["NSE_CUDA_ARCHITECTURES"] = libraries[
                    "cuda_architectures"
                ]
        return use_mpi, use_openmp, gpu_backend

    if adapter != "gp3d_cmake_v1":
        raise ModelBuildError(f"unsupported build adapter: {adapter}")

    use_mpi = _enabled(cmake_variables.get("USE_MPI", False))
    use_openmp = _enabled(cmake_variables.get("USE_OPENMP", False))
    gpu_backend = str(cmake_variables.get("GPU_BACKEND", "none")).lower()
    if use_openmp and not bool(openmp.get("enabled", False)):
        raise ModelBuildError(
            "GPE profile requires OpenMP, but machine openmp.enabled is false"
        )
    if gpu_backend == "cuda" and use_mpi:
        raise ModelBuildError("GPE cuda profile is single-GPU and cannot enable MPI")
    if gpu_backend == "cufftmp" and not use_mpi:
        raise ModelBuildError("GPE cuFFTMp profile requires MPI")
    if gpu_backend == "cufftmp" and os.name == "nt":
        raise ModelBuildError("GPE cuFFTMp profile can only be built on Linux")
    fft_decomposition = str(
        cmake_variables.get("FFT_DECOMPOSITION", "slab")
    ).lower()
    if fft_decomposition not in {"slab", "pencil"}:
        raise ModelBuildError("GPE FFT_DECOMPOSITION must be slab or pencil")
    if fft_decomposition == "pencil":
        if not use_mpi:
            raise ModelBuildError("GPE pencil FFT requires MPI")
        if gpu_backend == "cuda":
            raise ModelBuildError(
                "GPE pencil FFT is unavailable for the single-GPU CUDA profile"
            )
        if (
            gpu_backend == "none"
            and str(cmake_variables.get("FFT_BACKEND", "dft")).lower()
            != "fftw"
        ):
            raise ModelBuildError("GPE CPU pencil FFT requires FFT_BACKEND=fftw")
        if (
            gpu_backend == "cufftmp"
            and str(cmake_variables.get("CUFFTMP_API", "auto")).lower()
            == "legacy"
        ):
            raise ModelBuildError(
                "GPE cuFFTMp pencil FFT requires the modern custom-decomposition API"
            )

    cmake_variables.update(
        {
            "GP3D_ENABLE_WARNINGS": warnings,
            "GP3D_COMMON_FLAGS": common_flags,
            "GP3D_WARNING_FLAGS": warning_flags,
            "GP3D_DEBUG_FLAGS": debug_flags,
            "GP3D_RELEASE_FLAGS": release_flags,
            "GP3D_RELWITHDEBINFO_FLAGS": relwithdebinfo_flags,
        }
    )
    if use_mpi and mpi.get("root"):
        cmake_variables["MPI_ROOT"] = str(mpi["root"])
    if str(cmake_variables.get("FFT_BACKEND", "dft")).lower() == "fftw":
        if libraries.get("fftw_root"):
            cmake_variables["FFTW_ROOT"] = str(libraries["fftw_root"])
    if (
        fft_decomposition == "pencil"
        and gpu_backend == "none"
        and libraries.get("decomp2d_root")
    ):
        cmake_variables["GP3D_2DECOMP_ROOT"] = str(libraries["decomp2d_root"])
    if gpu_backend in {"cuda", "cufftmp"}:
        if libraries.get("cuda_compiler"):
            cmake_variables["CMAKE_CUDA_COMPILER"] = str(libraries["cuda_compiler"])
        if libraries.get("cuda_toolkit_root"):
            cmake_variables["CUDAToolkit_ROOT"] = str(libraries["cuda_toolkit_root"])
        if libraries.get("cuda_architectures") is not None:
            cmake_variables["GP3D_CUDA_ARCHITECTURES"] = libraries["cuda_architectures"]
    if gpu_backend == "cufftmp":
        if libraries.get("cufftmp_root"):
            cmake_variables["CUFFTMP_ROOT"] = str(libraries["cufftmp_root"])
        if libraries.get("nvshmem_root"):
            cmake_variables["NVSHMEM_ROOT"] = str(libraries["nvshmem_root"])
        if libraries.get("cufftmp_api"):
            cmake_variables["CUFFTMP_API"] = str(libraries["cufftmp_api"])
    return use_mpi, use_openmp, gpu_backend


def _resolve(args: argparse.Namespace) -> ResolvedBuild:
    design_path = Path(args.design or DEFAULT_DESIGN).resolve()
    if not design_path.is_file():
        raise ModelBuildError(f"build design not found: {design_path}")
    design = _mapping(load_yaml(design_path), "build design")
    _schema_one(design, "build design")

    catalog_path = _existing_file(
        design_path.parent,
        str(design.get("model_catalog") or "model_catalog.yaml"),
        "model catalog",
    )
    catalog = _mapping(load_yaml(catalog_path), "model catalog")
    _schema_one(catalog, "model catalog")
    models_catalog = _mapping(catalog.get("models"), "model catalog.models")

    model = _safe_name(args.model or design.get("selected_model"), "selected_model").lower()
    if model not in models_catalog:
        raise ModelBuildError(f"unknown model {model!r}; available: {sorted(models_catalog)}")
    catalog_entry = _mapping(models_catalog[model], f"model catalog.{model}")
    design_models = _mapping(design.get("models"), "build design.models")
    model_design = _mapping(design_models.get(model, {}), f"build design.models.{model}")

    machine_text = args.machine_profile or design.get("machine_profile")
    if not machine_text:
        raise ModelBuildError("machine_profile is required")
    machine_path = _existing_file(
        design_path.parent, str(machine_text), "machine profile"
    )
    machine = _mapping(load_yaml(machine_path), "machine profile")
    _schema_one(machine, "machine profile")
    expected_platform = str(machine.get("platform") or "").lower()
    if expected_platform and expected_platform != _platform_name():
        raise ModelBuildError(
            f"machine profile is for {expected_platform}, but this host is {_platform_name()}"
        )

    library_text = (
        args.solver_library_root
        or os.environ.get("SOLVER_LIBRARY_ROOT")
        or design.get("solver_library_root")
    )
    if not library_text:
        raise ModelBuildError(
            "solver_library_root is required in YAML, SOLVER_LIBRARY_ROOT, or command line"
        )
    solver_library_root = _expand_path(
        str(library_text), design_path.parent, "solver_library_root"
    )
    if not solver_library_root.is_dir():
        raise ModelBuildError(f"SolverLibrary root not found: {solver_library_root}")
    solver_root = (solver_library_root / str(catalog_entry["library_subpath"])).resolve()
    if not solver_root.is_dir():
        raise ModelBuildError(f"solver root not found for model {model}: {solver_root}")
    manifest_path = _existing_file(
        solver_root,
        str(catalog_entry.get("manifest") or "solver_manifest.yaml"),
        "solver manifest",
    )
    manifest = _mapping(load_yaml(manifest_path), "solver manifest")
    _schema_one(manifest, "solver manifest")
    if str(manifest.get("model") or "").lower() != model:
        raise ModelBuildError(
            f"solver manifest model {manifest.get('model')!r} does not match {model!r}"
        )
    adapter = str(manifest.get("build_adapter") or "")
    required_adapter = str(catalog_entry.get("required_adapter") or "")
    if required_adapter and adapter != required_adapter:
        raise ModelBuildError(
            f"model {model} requires adapter {required_adapter}, manifest provides {adapter}"
        )

    solver_profile = _safe_name(
        args.solver_profile
        or model_design.get("profile")
        or catalog_entry.get("default_profile"),
        f"models.{model}.profile",
    )
    build = _mapping(design.get("build"), "build design.build")
    include_tests = bool(build.get("tests", False) or args.test or args.all)
    profiles = _mapping(manifest.get("profiles"), "manifest.profiles")
    if solver_profile not in profiles:
        raise ModelBuildError(
            f"unknown solver profile {solver_profile!r}; available: {sorted(profiles)}"
        )
    profile = _mapping(
        profiles[solver_profile], f"solver profile {solver_profile}"
    )
    components, selected_files, library_sources, executable_sources = _resolve_components(
        solver_root, manifest, solver_profile, include_tests
    )
    preprocessor_defines = {
        str(value)
        for value in _list(
            profile.get("fortran_preprocessor_defines"),
            f"solver profile {solver_profile}.fortran_preprocessor_defines",
        )
    }
    _inspect_fortran_dependencies(
        solver_root, manifest, selected_files, preprocessor_defines
    )
    executable = _safe_name(profile.get("executable"), "solver profile.executable")

    configuration = str(args.configuration or build.get("configuration") or "Release")
    if configuration not in {"Debug", "Release", "RelWithDebInfo"}:
        raise ModelBuildError("build.configuration must be Debug, Release, or RelWithDebInfo")
    if args.build_dir:
        build_dir = Path(args.build_dir).resolve()
    else:
        output_root = _expand_path(
            str(build.get("output_root") or "build"), design_path.parent, "build.output_root"
        )
        template = str(build.get("directory_name") or "{model}-{profile}-{configuration}")
        try:
            directory_name = template.format(
                model=model,
                profile=solver_profile,
                configuration=configuration.lower(),
                design_id=_safe_name(design.get("design_id"), "design_id"),
            )
        except KeyError as exc:
            raise ModelBuildError(f"unknown build.directory_name placeholder: {exc}") from exc
        build_dir = (output_root / _safe_name(directory_name, "build.directory_name")).resolve()
    if os.name == "nt" and str(build_dir).startswith("\\\\"):
        raise ModelBuildError(
            "Windows build directory must be local or a mapped drive, not a UNC path"
        )

    compiler = _mapping(machine.get("compiler"), "machine profile.compiler")
    cmake_variables: dict[str, Any] = {
        "CMAKE_Fortran_COMPILER": str(compiler.get("fortran") or "gfortran"),
        "CMAKE_BUILD_TYPE": configuration,
        "BUILD_TESTING": include_tests,
    }
    cmake_variables.update(_mapping(profile.get("cmake", {}), "solver profile.cmake"))
    cmake_variables.update(
        _mapping(model_design.get("cmake_overrides", {}), f"models.{model}.cmake_overrides")
    )
    use_mpi, use_openmp, gpu_backend = _model_machine_variables(
        adapter,
        machine,
        build,
        profile,
        library_sources,
        executable_sources,
        executable,
        build_dir,
        cmake_variables,
    )

    input_config = _mapping(manifest.get("input", {}), "solver manifest.input")
    input_name = str(input_config.get("default_name") or "input.dat")
    pass_input = bool(input_config.get("pass_as_argument", False))
    return ResolvedBuild(
        design_path=design_path,
        catalog_path=catalog_path,
        machine_path=machine_path,
        manifest_path=manifest_path,
        solver_library_root=solver_library_root,
        solver_root=solver_root,
        build_dir=build_dir,
        design=design,
        machine=machine,
        manifest=manifest,
        model=model,
        solver_profile=solver_profile,
        adapter=adapter,
        configuration=configuration,
        executable=executable,
        selected_components=tuple(components),
        selected_files=tuple(selected_files),
        cmake_variables=cmake_variables,
        use_mpi=use_mpi,
        use_openmp=use_openmp,
        gpu_backend=gpu_backend,
        input_name=input_name,
        pass_input_as_argument=pass_input,
    )


def _cmake_bracket(value: str) -> str:
    marker = "="
    while f"]{marker}]" in value:
        marker += "="
    return f"[{marker}[{value}]{marker}]"


def _cache_type(name: str, value: Any) -> tuple[str, str]:
    if isinstance(value, bool):
        return "BOOL", "ON" if value else "OFF"
    if isinstance(value, list):
        return "STRING", ";".join(str(item) for item in value)
    text = str(value)
    file_paths = {"CMAKE_Fortran_COMPILER", "CMAKE_CUDA_COMPILER", "FFTW3_LIBRARY"}
    path_names = {
        "CMAKE_RUNTIME_OUTPUT_DIRECTORY",
        "MSMPI_ROOT",
        "MPI_ROOT",
        "FFTW_ROOT",
        "NSE_2DECOMP_ROOT",
        "CUDAToolkit_ROOT",
        "CUFFTMP_ROOT",
        "NVSHMEM_ROOT",
    }
    if name in file_paths:
        return "FILEPATH", text.replace("\\", "/")
    if name in path_names:
        return "PATH", text.replace("\\", "/")
    return "STRING", text


def _initial_cache(resolved: ResolvedBuild) -> str:
    lines = [
        "# Generated by ScriptLibrary/BuildSolver/build_model.py. Do not edit.",
        f"# Model: {resolved.model}",
        f"# Solver profile: {resolved.solver_profile}",
        "",
    ]
    for name, raw_value in resolved.cmake_variables.items():
        if raw_value is None:
            continue
        cache_type, value = _cache_type(name, raw_value)
        lines.append(
            f"set({name} {_cmake_bracket(value)} CACHE {cache_type} "
            '"Generated from unified YAML" FORCE)'
        )
    return "\n".join(lines) + "\n"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_path(path: Path) -> str:
    return str(path).replace("\\", "/")


def _resolved_document(resolved: ResolvedBuild) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "model": resolved.model,
        "solver_profile": resolved.solver_profile,
        "adapter": resolved.adapter,
        "configuration": resolved.configuration,
        "executable": resolved.executable,
        "use_mpi": resolved.use_mpi,
        "use_openmp": resolved.use_openmp,
        "gpu_backend": resolved.gpu_backend,
        "design": _json_path(resolved.design_path),
        "design_sha256": _sha256(resolved.design_path),
        "model_catalog": _json_path(resolved.catalog_path),
        "model_catalog_sha256": _sha256(resolved.catalog_path),
        "machine_profile": _json_path(resolved.machine_path),
        "machine_profile_sha256": _sha256(resolved.machine_path),
        "solver_manifest": _json_path(resolved.manifest_path),
        "solver_manifest_sha256": _sha256(resolved.manifest_path),
        "solver_root": _json_path(resolved.solver_root),
        "build_directory": _json_path(resolved.build_dir),
        "selected_components": list(resolved.selected_components),
        "selected_files": list(resolved.selected_files),
        "cmake_variables": resolved.cmake_variables,
    }


def _write_generated(resolved: ResolvedBuild) -> tuple[Path, Path]:
    directory = resolved.build_dir / "generated"
    directory.mkdir(parents=True, exist_ok=True)
    cache_path = directory / "InitialCache.cmake"
    resolved_path = directory / "resolved_build.json"
    cache_path.write_text(_initial_cache(resolved), encoding="utf-8")
    resolved_path.write_text(
        json.dumps(_resolved_document(resolved), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return cache_path, resolved_path


def _tool(machine: dict[str, Any], name: str, default: str) -> str:
    tools = _mapping(machine.get("tools", {}), "machine profile.tools")
    return str(tools.get(name) or default)


def _display(command: list[str]) -> str:
    return subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)


def _run_command(
    command: list[str], cwd: Path, environment: dict[str, str], dry_run: bool
) -> None:
    print(f"[CMD] {_display(command)}", flush=True)
    print(f"      cwd={cwd}", flush=True)
    if dry_run:
        return
    result = subprocess.run(command, cwd=str(cwd), env=environment, check=False)
    if result.returncode != 0:
        raise ModelBuildError(
            f"command failed with exit code {result.returncode}: {_display(command)}"
        )


def _base_environment(resolved: ResolvedBuild, args: argparse.Namespace) -> dict[str, str]:
    environment = os.environ.copy()
    machine_environment = _mapping(
        resolved.machine.get("environment", {}), "machine profile.environment"
    )
    run = _mapping(resolved.design.get("run", {}), "build design.run")
    run_environment = _mapping(run.get("environment", {}), "run.environment")
    for source in (machine_environment, run_environment):
        for name, value in source.items():
            environment[str(name)] = str(value)
    if resolved.use_openmp:
        threads = int(args.omp_threads or run.get("omp_threads", 1))
        if threads < 1:
            raise ModelBuildError("OpenMP thread count must be at least 1")
        environment["OMP_NUM_THREADS"] = str(threads)
    libraries = _mapping(resolved.machine.get("libraries", {}), "machine libraries")
    if resolved.gpu_backend == "cufftmp" and libraries.get("nvshmem_symmetric_size"):
        environment["NVSHMEM_SYMMETRIC_SIZE"] = str(
            libraries["nvshmem_symmetric_size"]
        )
    return environment


def _with_visual_studio(environment: dict[str, str], cwd: Path) -> dict[str, str]:
    if os.name != "nt":
        return environment
    if shutil.which("cl.exe", path=environment.get("PATH")) and shutil.which(
        "lib.exe", path=environment.get("PATH")
    ):
        return environment
    program_files = environment.get("ProgramFiles(x86)") or os.environ.get(
        "ProgramFiles(x86)", ""
    )
    vswhere = Path(program_files) / "Microsoft Visual Studio" / "Installer" / "vswhere.exe"
    if not vswhere.is_file():
        raise ModelBuildError("CUDA build requires Visual Studio C++ tools; vswhere.exe was not found")
    query = subprocess.run(
        [
            str(vswhere),
            "-latest",
            "-products",
            "*",
            "-requires",
            "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
            "-property",
            "installationPath",
        ],
        cwd=str(cwd),
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    installation = query.stdout.strip().splitlines()
    if query.returncode != 0 or not installation:
        raise ModelBuildError("Visual Studio x64 C++ build tools were not found")
    vcvars = Path(installation[0]) / "VC" / "Auxiliary" / "Build" / "vcvars64.bat"
    if not vcvars.is_file():
        raise ModelBuildError(f"vcvars64.bat was not found: {vcvars}")
    comspec = environment.get("ComSpec") or os.environ.get("ComSpec") or "cmd.exe"
    marker = "__BUILD_MODEL_VISUAL_STUDIO_ENV__"
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="mbcs",
        suffix=".cmd",
        prefix="load-vs-env-",
        dir=cwd,
        delete=False,
    ) as command_file:
        command_file.write("@echo off\n")
        command_file.write(f'call "{vcvars}"\n')
        command_file.write("if errorlevel 1 exit /b %errorlevel%\n")
        command_file.write(f"echo {marker}\n")
        command_file.write("set\n")
        command_path = Path(command_file.name)
    try:
        loaded = subprocess.run(
            [comspec, "/d", "/s", "/c", command_path.name],
            cwd=str(cwd),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="mbcs",
            errors="replace",
            check=False,
        )
    finally:
        command_path.unlink(missing_ok=True)
    if loaded.returncode != 0:
        details = (loaded.stderr.strip() or loaded.stdout.strip())[-2000:]
        suffix = f":\n{details}" if details else ""
        raise ModelBuildError(
            f"failed to initialize the Visual Studio x64 environment{suffix}"
        )
    if marker not in loaded.stdout:
        raise ModelBuildError(
            "Visual Studio environment initialized without returning environment variables"
        )
    result = environment.copy()
    environment_text = loaded.stdout.split(marker, 1)[1]
    for line in environment_text.splitlines():
        if "=" in line:
            name, value = line.split("=", 1)
            result[name] = value
    if not shutil.which("cl.exe", path=result.get("PATH")):
        raise ModelBuildError("Visual Studio environment loaded, but cl.exe is unavailable")
    print(f"[OK] Visual Studio environment: {vcvars}", flush=True)
    return result


def _configure_command(
    resolved: ResolvedBuild, cache_path: Path, fresh: bool
) -> list[str]:
    command = [_tool(resolved.machine, "cmake", "cmake")]
    if fresh:
        command.append("--fresh")
    command.extend(
        [
            "-S",
            str(resolved.solver_root),
            "-B",
            str(resolved.build_dir),
            "-G",
            str(resolved.machine.get("generator") or "Ninja"),
            "-C",
            str(cache_path),
        ]
    )
    return command


def _build_command(resolved: ResolvedBuild, clean_first: bool) -> list[str]:
    build = _mapping(resolved.design["build"], "build design.build")
    command = [
        _tool(resolved.machine, "cmake", "cmake"),
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


def _test_command(resolved: ResolvedBuild) -> list[str]:
    return [
        _tool(resolved.machine, "ctest", "ctest"),
        "--test-dir",
        str(resolved.build_dir),
        "--build-config",
        resolved.configuration,
        "--output-on-failure",
    ]


def _model_design(resolved: ResolvedBuild) -> dict[str, Any]:
    models = _mapping(resolved.design.get("models"), "build design.models")
    return _mapping(models.get(resolved.model, {}), f"models.{resolved.model}")


def _input_path(args: argparse.Namespace, resolved: ResolvedBuild) -> Path:
    model_design = _model_design(resolved)
    value = args.input_file or model_design.get("input_file")
    if not value:
        raise ModelBuildError(
            f"{resolved.model} execution requires an input file; use --input-file or "
            f"models.{resolved.model}.input_file"
        )
    path = Path(str(value))
    if not path.is_absolute():
        path = (Path.cwd() if args.input_file else resolved.design_path.parent) / path
    path = path.resolve()
    if not path.is_file():
        raise ModelBuildError(f"input file not found: {path}")
    return path


def _executable_path(resolved: ResolvedBuild, dry_run: bool) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    name = resolved.executable + suffix
    candidates = [
        resolved.build_dir / "bin" / name,
        resolved.build_dir / "bin" / resolved.configuration / name,
        resolved.build_dir / name,
        resolved.build_dir / resolved.configuration / name,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    if dry_run:
        return candidates[0]
    raise ModelBuildError(
        "solver executable not found; expected one of: "
        + ", ".join(str(value) for value in candidates)
    )


def _run_solver(
    args: argparse.Namespace,
    resolved: ResolvedBuild,
    environment: dict[str, str],
) -> None:
    source_input = _input_path(args, resolved)
    run = _mapping(resolved.design.get("run", {}), "build design.run")
    run_dir_value = args.run_dir or run.get("working_directory")
    if run_dir_value:
        run_dir = Path(str(run_dir_value))
        if not run_dir.is_absolute():
            base = Path.cwd() if args.run_dir else resolved.design_path.parent
            run_dir = base / run_dir
        run_dir = run_dir.resolve()
    else:
        run_dir = resolved.build_dir / "run"
    staged_input = run_dir / resolved.input_name
    if not args.dry_run:
        run_dir.mkdir(parents=True, exist_ok=True)
        if source_input != staged_input.resolve():
            shutil.copy2(source_input, staged_input)
    executable = _executable_path(resolved, args.dry_run)
    program_arguments: list[str] = []
    if resolved.pass_input_as_argument:
        program_arguments.append(str(staged_input))
    program_arguments.extend(
        str(value) for value in _list(run.get("program_arguments"), "run.program_arguments")
    )
    command: list[str] = []
    if resolved.use_mpi:
        processes = int(args.processes or run.get("mpi_processes", 1))
        if processes < 1:
            raise ModelBuildError("MPI process count must be at least 1")
        command.extend(
            [
                _tool(resolved.machine, "mpi_launcher", "mpiexec"),
                _tool(resolved.machine, "mpi_process_option", "-n"),
                str(processes),
            ]
        )
        command.extend(
            str(value)
            for value in _list(run.get("launcher_arguments"), "run.launcher_arguments")
        )
    command.extend([str(executable), *program_arguments])
    _run_command(command, run_dir, environment, args.dry_run)


def _summary(resolved: ResolvedBuild) -> None:
    print("[OK] Unified build design is valid")
    print(f"     model:          {resolved.model}")
    print(f"     solver profile: {resolved.solver_profile}")
    print(f"     adapter:        {resolved.adapter}")
    print(f"     configuration:  {resolved.configuration}")
    print(f"     solver root:    {resolved.solver_root}")
    print(f"     build dir:      {resolved.build_dir}")
    print(
        "     MPI/OpenMP/GPU: "
        f"{resolved.use_mpi}/{resolved.use_openmp}/{resolved.gpu_backend}",
        flush=True,
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and run NSE or GPE from one model-selecting YAML design."
    )
    parser.add_argument("design", nargs="?", help="Build YAML; default: build.yaml")
    parser.add_argument(
        "--model",
        help="Override selected_model; the value must exist in model_catalog.yaml",
    )
    parser.add_argument(
        "--solver-profile", "--profile", dest="solver_profile", help="Override solver profile"
    )
    parser.add_argument("--machine-profile", help="Override machine profile YAML")
    parser.add_argument("--solver-library-root", help="Override SolverLibrary root")
    parser.add_argument("--configuration", choices=["Debug", "Release", "RelWithDebInfo"])
    parser.add_argument("--build-dir", help="Override local build directory")
    parser.add_argument("--input-file", help="Model input file used by --run")
    parser.add_argument(
        "--run-dir",
        help="Run working directory; defaults to <build-dir>/run",
    )
    parser.add_argument("--processes", type=int, help="MPI process count")
    parser.add_argument("--omp-threads", type=int, help="OpenMP threads per rank")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--list-models", action="store_true")
    parser.add_argument("--list-profiles", action="store_true")
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--configure", action="store_true", help="Configure CMake only")
    parser.add_argument("--build", action="store_true", help="Configure and build only")
    parser.add_argument("--test", action="store_true", help="Build and run tests")
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run an existing executable without configuring or building",
    )
    parser.add_argument("--all", action="store_true", help="Configure, build, test, and run")
    parser.add_argument("--clean-first", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    fresh = parser.add_mutually_exclusive_group()
    fresh.add_argument("--fresh", dest="fresh", action="store_true")
    fresh.add_argument("--no-fresh", dest="fresh", action="store_false")
    parser.set_defaults(fresh=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.validate_only and any(
            [args.generate_only, args.configure, args.build, args.test, args.run, args.all]
        ):
            raise ModelBuildError("--validate-only cannot be combined with build stages")
        if args.generate_only and any(
            [args.configure, args.build, args.test, args.run, args.all]
        ):
            raise ModelBuildError("--generate-only cannot be combined with build stages")

        resolved = _resolve(args)
        _summary(resolved)
        if args.list_models:
            catalog = _mapping(load_yaml(resolved.catalog_path), "model catalog")
            models = _mapping(catalog.get("models"), "model catalog.models")
            print("Available models:")
            for name, value in models.items():
                record = _mapping(value, f"model {name}")
                print(f"  {name}: {record.get('description', '')}")
        if args.list_profiles:
            profiles = _mapping(resolved.manifest.get("profiles"), "solver manifest.profiles")
            print(f"Available profiles for {resolved.model}:")
            for name, value in profiles.items():
                record = _mapping(value, f"profile {name}")
                print(f"  {name}: {record.get('description', '')}")
        if args.list_models or args.list_profiles:
            return 0
        if args.validate_only:
            return 0

        explicit = any([args.configure, args.build, args.test, args.run, args.all])
        do_build = bool(args.build or args.test or args.all)
        do_test = bool(args.test or args.all)
        do_run = bool(args.run or args.all)
        do_configure = bool(args.configure or do_build or do_test)
        if not explicit and not args.generate_only:
            do_configure = True
            do_build = True
            run = _mapping(resolved.design.get("run", {}), "build design.run")
            do_run = bool(run.get("enabled", False))

        run_only = do_run and not any([do_configure, do_build, do_test])
        if run_only:
            environment = _base_environment(resolved, args)
            _run_solver(args, resolved, environment)
            return 0

        build = _mapping(resolved.design["build"], "build design.build")
        fresh = bool(build.get("configure_fresh", True)) if args.fresh is None else args.fresh
        cache_path = resolved.build_dir / "generated" / "InitialCache.cmake"
        if args.dry_run:
            print(f"[DRY-RUN] generate {cache_path}")
        else:
            cache_path, resolved_path = _write_generated(resolved)
            print(f"[OK] Generated CMake cache: {cache_path}")
            print(f"[OK] Resolved build plan:  {resolved_path}", flush=True)
        if args.generate_only:
            return 0

        environment = _base_environment(resolved, args)
        if resolved.gpu_backend == "cuda" and os.name == "nt" and not args.dry_run:
            resolved.build_dir.mkdir(parents=True, exist_ok=True)
            environment = _with_visual_studio(environment, resolved.build_dir)
        if do_configure:
            if not args.dry_run:
                resolved.build_dir.mkdir(parents=True, exist_ok=True)
            _run_command(
                _configure_command(resolved, cache_path, fresh),
                resolved.build_dir,
                environment,
                args.dry_run,
            )
        if do_build:
            _run_command(
                _build_command(resolved, args.clean_first),
                resolved.build_dir,
                environment,
                args.dry_run,
            )
        if do_test:
            _run_command(
                _test_command(resolved),
                resolved.build_dir,
                environment,
                args.dry_run,
            )
        if do_run:
            _run_solver(args, resolved, environment)
        return 0
    except (ModelBuildError, YamlFormatError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
