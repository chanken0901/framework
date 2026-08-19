#!/usr/bin/env python3
"""フレームワークのケースYAMLを検証し、GP3D用namelistへ変換する。

このファイルはYAML側の名前とFortran側の入力契約をつなぐアダプターである。
物理計算は行わず、未対応の初期条件やパラメータを実行前に検出する。
"""

from __future__ import annotations

import hashlib
import json
from typing import Any


# flow.type に記述できる初期条件名と互換エイリアス。
SUPPORTED_INITIAL_CONDITIONS = {
    "gaussian",
    "uniform_vortex",
    "vortex_uniform",
    "tf_vortex",
    "thomas_fermi_vortex",
    "vortex",
    "vortex_ring",
    "ring",
    "vortex_tangle",
    "random_vortices",
    "quantum_turbulence",
    "ring_tangle",
    "vortex_ring_tangle",
    "random_rings",
    "quantum_taylor_green",
    "taylor_green",
    "tg",
    "restart_slf",
    "restart",
    "slf",
}

# Fortranの &gpe namelistへそのまま渡せる正式なキー。
GPE_KEYS = [
    "use_dimensionless_parameters",
    "alpha",
    "beta",
    "g",
    "sigma0",
    "wx",
    "wy",
    "wz",
    "hbar",
    "mass",
    "norm",
    "imaginary_time",
    "mu",
    "healing_length",
    "vortex_charge",
    "vortex_x0",
    "vortex_y0",
    "ring_radius",
    "ring_z0",
    "phase_noise",
    "tangle_nlines",
    "tangle_nrings",
    "ring_radius_min",
    "ring_radius_max",
    "random_seed",
    "tg_velocity_amplitude",
    "tg_winding",
    "tg_auto_winding",
    "argle_enabled",
    "argle_write_seed",
    "argle_steps",
    "argle_output_every",
    "argle_dtau",
    "argle_tolerance",
]

LEGACY_GPE_ALIASES = {
    "interaction_strength": "g",
    "chemical_potential": "mu",
}

UNSUPPORTED_LEGACY_GPE_KEYS = {"density0", "damping"}


class CaseValidationError(ValueError):
    pass


def nested(data: dict[str, Any], path: str, default: Any = None) -> Any:
    current: Any = data
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def _required(data: dict[str, Any], path: str, errors: list[str]) -> Any:
    value = nested(data, path)
    if value is None or value == "":
        errors.append(f"missing required value: {path}")
    return value


def validate_case(case: dict[str, Any], manifest: dict[str, Any], profile_name: str) -> None:
    """必須項目、数値範囲、初期条件、実行プロファイルの整合性を検査する。"""
    errors: list[str] = []
    if not isinstance(case, dict):
        raise CaseValidationError("case YAML must contain a top-level mapping")

    _required(case, "case_id", errors)
    model = str(_required(case, "physics.model", errors) or "").strip().lower()
    if model not in {"gpe", "gross_pitaevskii", "gross-pitaevskii"}:
        errors.append(f"physics.model must be gpe for solver gp3d; got {model!r}")

    profiles = manifest.get("profiles", {})
    if profile_name not in profiles:
        errors.append(
            f"unknown solver.profile {profile_name!r}; available profiles: {sorted(profiles)}"
        )

    initial_condition = str(_required(case, "flow.type", errors) or "").strip().lower()
    if initial_condition and initial_condition not in SUPPORTED_INITIAL_CONDITIONS:
        errors.append(f"unsupported flow.type for gp3d: {initial_condition!r}")

    for axis in "xyz":
        n = _required(case, f"grid.n{axis}", errors)
        lower = _required(case, f"grid.{axis}_min", errors)
        upper = _required(case, f"grid.{axis}_max", errors)
        if n is not None and (not isinstance(n, int) or isinstance(n, bool) or n < 2):
            errors.append(f"grid.n{axis} must be an integer >= 2")
        if isinstance(lower, (int, float)) and isinstance(upper, (int, float)) and lower >= upper:
            errors.append(f"grid.{axis}_min must be smaller than grid.{axis}_max")

    dt = _required(case, "time.dt", errors)
    nsteps = _required(case, "time.nsteps", errors)
    output_frequency = _required(case, "time.output_frequency", errors)
    if isinstance(dt, (int, float)) and dt <= 0:
        errors.append("time.dt must be positive")
    if nsteps is not None and (not isinstance(nsteps, int) or isinstance(nsteps, bool) or nsteps < 0):
        errors.append("time.nsteps must be a non-negative integer")
    if output_frequency is not None and (
        not isinstance(output_frequency, int) or isinstance(output_frequency, bool) or output_frequency < 0
    ):
        errors.append("time.output_frequency must be a non-negative integer")

    gpe = nested(case, "physics.gpe", {})
    if not isinstance(gpe, dict):
        errors.append("physics.gpe must be a mapping")
    else:
        unsupported = sorted(UNSUPPORTED_LEGACY_GPE_KEYS.intersection(gpe))
        if unsupported:
            errors.append(
                "unsupported legacy GPE keys have no GP3D namelist equivalent: "
                + ", ".join(unsupported)
            )
        allowed = set(GPE_KEYS) | set(LEGACY_GPE_ALIASES)
        unknown = sorted(set(gpe) - allowed)
        if unknown:
            errors.append("unknown physics.gpe keys: " + ", ".join(unknown))

    restart_file = nested(case, "restart.file")
    if initial_condition in {"restart_slf", "restart", "slf"} and not restart_file:
        errors.append("restart.file is required when flow.type selects restart")

    profile = profiles.get(profile_name, {})
    definitions = profile.get("cmake", {}) if isinstance(profile, dict) else {}
    use_mpi = bool(definitions.get("USE_MPI", False))
    openmp_capable = bool(definitions.get("USE_OPENMP", False))
    gpu_backend = str(definitions.get("GPU_BACKEND", "none")).lower()
    use_cuda = gpu_backend in {"cuda", "cufftmp"}
    processes = nested(case, "solver.mpi_processes", 1)
    if not isinstance(processes, int) or isinstance(processes, bool) or processes < 1:
        errors.append("solver.mpi_processes must be a positive integer")
    elif not use_mpi and processes != 1:
        errors.append(
            f"solver.mpi_processes must be 1 for non-MPI profile {profile_name}"
        )
    if use_mpi and gpu_backend == "cuda":
        errors.append("the single-GPU CUDA profile cannot be combined with MPI")
    use_openmp = nested(case, "solver.use_openmp", False)
    omp_threads = nested(case, "solver.omp_threads", 1)
    if not isinstance(use_openmp, bool):
        errors.append("solver.use_openmp must be true or false")
    elif use_openmp and not openmp_capable:
        errors.append(
            f"solver.use_openmp requires an OpenMP-capable profile; {profile_name} is not one"
        )
    if not isinstance(omp_threads, int) or isinstance(omp_threads, bool) or omp_threads < 1:
        errors.append("solver.omp_threads must be a positive integer")

    if errors:
        raise CaseValidationError("invalid case YAML:\n  - " + "\n  - ".join(errors))


def _fortran_value(value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return '"' + value.replace('"', '""') + '"'
    if isinstance(value, float):
        return format(value, ".17g")
    if isinstance(value, int):
        return str(value)
    raise TypeError(f"unsupported namelist value: {value!r}")


def _append_values(lines: list[str], values: list[tuple[str, Any]]) -> None:
    for key, value in values:
        if value is not None and value != "":
            lines.append(f"  {key} = {_fortran_value(value)}")


def generate_input_namelist(
    case: dict[str, Any], manifest: dict[str, Any], profile_name: str
) -> str:
    """検証済みケースをFortranの ``&simulation`` と ``&gpe`` へ変換する。"""
    validate_case(case, manifest, profile_name)
    profile = manifest["profiles"][profile_name]
    definitions = profile["cmake"]
    use_mpi = bool(definitions.get("USE_MPI", False))
    gpu_backend = str(definitions.get("GPU_BACKEND", "none")).lower()
    use_cuda = gpu_backend in {"cuda", "cufftmp"}
    fft_backend = str(definitions.get("FFT_BACKEND", "dft"))
    if gpu_backend == "cufftmp":
        backend = "cufftmp"
    elif gpu_backend == "cuda":
        backend = "cufft"
    else:
        backend = fft_backend

    simulation = [
        ("equation", "GPE"),
        ("case_name", str(case["case_id"])),
        ("initial_condition", str(nested(case, "flow.type"))),
        ("restart_file", nested(case, "restart.file")),
        ("nx", nested(case, "grid.nx")),
        ("ny", nested(case, "grid.ny")),
        ("nz", nested(case, "grid.nz")),
        ("x_min", nested(case, "grid.x_min")),
        ("x_max", nested(case, "grid.x_max")),
        ("y_min", nested(case, "grid.y_min")),
        ("y_max", nested(case, "grid.y_max")),
        ("z_min", nested(case, "grid.z_min")),
        ("z_max", nested(case, "grid.z_max")),
        ("dt", nested(case, "time.dt")),
        ("t_max", nested(case, "time.t_max")),
        ("cfl", nested(case, "time.cfl")),
        ("nsteps", nested(case, "time.nsteps")),
        ("output_frequency", nested(case, "time.output_frequency")),
        ("output_dir", nested(case, "output.directory", "output")),
        ("output_format", nested(case, "output.format", "slf")),
        ("backend", backend),
        ("write_initial", nested(case, "output.write_initial", True)),
        ("write_meta", nested(case, "output.write_meta", True)),
        ("timing_enabled", nested(case, "output.timing_enabled", False)),
        ("use_mpi", use_mpi),
        ("use_openmp", nested(case, "solver.use_openmp", False)),
        ("use_cuda", use_cuda),
        ("cuda_device", nested(case, "solver.cuda_device", 0)),
        ("rank", 0),
        ("nprocs", nested(case, "solver.mpi_processes", 1)),
    ]

    gpe = dict(nested(case, "physics.gpe", {}))
    for legacy, canonical in LEGACY_GPE_ALIASES.items():
        if canonical not in gpe and legacy in gpe:
            gpe[canonical] = gpe[legacy]

    lines = [
        "! Automatically generated from case.yaml.",
        "! Do not edit this file directly.",
        "",
        "&simulation",
    ]
    _append_values(lines, simulation)
    lines.extend(["/", "", "&gpe"])
    _append_values(lines, [(key, gpe.get(key)) for key in GPE_KEYS])
    lines.extend(["/", ""])
    return "\n".join(lines)


def case_digest(case: dict[str, Any]) -> str:
    """再現性記録に使う、キー順序に依存しないケースSHA-256を返す。"""
    canonical = json.dumps(case, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
