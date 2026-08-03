#!/usr/bin/env python3
"""Render model-specific Fortran namelists from a framework case.yaml."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from yaml_support import YamlFormatError, load_yaml


class CaseInputError(ValueError):
    """Raised when a case cannot be converted to a solver input file."""


GPE_KEYS = (
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
)

NSE_KEYS = (
    "nv",
    "nghost",
    "gamma",
    "cfl",
    "small_rho",
    "small_p",
    "rho0",
    "mach",
    "reynolds",
    "prandtl",
    "convective_scheme",
    "viscous_scheme",
    "boundary_condition",
    "time_integrator",
    "hit_spectrum",
    "hit_seed",
    "hit_rms_velocity",
    "hit_peak_wavenumber",
    "hit_integral_length",
    "hit_kolmogorov_length",
    "hit_dealias_fraction",
    "forcing_scheme",
    "forcing_spectrum",
    "forcing_fft_backend",
    "forcing_k_cutoff",
    "forcing_target_dissipation",
    "forcing_dilatational_ratio",
    "forcing_denominator_floor",
    "forcing_max_coefficient",
    "forcing_report_interval",
)

NSE_ALIASES = {
    "mach_number": "mach",
    "reynolds_number": "reynolds",
    "prandtl_number": "prandtl",
}

GPE_ALIASES = {
    "interaction_strength": "g",
    "chemical_potential": "mu",
}


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise CaseInputError(f"{label} must be a YAML mapping")
    return value


def nested(data: dict[str, Any], path: str, default: Any = None) -> Any:
    current: Any = data
    for name in path.split("."):
        if not isinstance(current, dict) or name not in current:
            return default
        current = current[name]
    return current


def _required(case: dict[str, Any], path: str) -> Any:
    value = nested(case, path)
    if value is None or value == "":
        raise CaseInputError(f"case is missing required value: {path}")
    return value


def _fortran(value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return '"' + value.replace('"', '""') + '"'
    if isinstance(value, float):
        return format(value, ".17g")
    if isinstance(value, int):
        return str(value)
    raise CaseInputError(f"unsupported namelist value: {value!r}")


def _append(lines: list[str], values: list[tuple[str, Any]]) -> None:
    for key, value in values:
        if value is not None and value != "":
            lines.append(f"  {key} = {_fortran(value)}")


def _canonical_selector(value: str) -> str:
    """Normalize human-readable selector spelling to the solver convention."""
    return "_".join(value.strip().lower().replace("-", " ").split())


def _profile_settings(
    manifest: dict[str, Any], profile_name: str
) -> tuple[dict[str, Any], bool, bool, str]:
    profiles = _mapping(manifest.get("profiles"), "solver manifest.profiles")
    if profile_name not in profiles:
        raise CaseInputError(
            f"unknown solver profile {profile_name!r}; available: {sorted(profiles)}"
        )
    profile = _mapping(profiles[profile_name], f"solver profile {profile_name}")
    cmake = _mapping(profile.get("cmake", {}), f"profile {profile_name}.cmake")
    execution = _mapping(
        profile.get("execution", {}), f"profile {profile_name}.execution"
    )
    use_mpi = bool(cmake.get("USE_MPI", execution.get("use_mpi", False)))
    use_openmp = bool(
        cmake.get(
            "USE_OPENMP",
            cmake.get("NSE_ENABLE_OPENMP", execution.get("use_openmp", False)),
        )
    )
    gpu_backend = str(
        cmake.get("GPU_BACKEND", cmake.get("NSE_GPU_BACKEND", "none"))
    ).lower()
    model = str(manifest.get("model", "")).lower()
    if gpu_backend == "cufftmp":
        backend = "cufftmp"
    elif gpu_backend == "cuda":
        backend = "cuda" if model == "nse" else "cufft"
    elif model == "nse":
        backend = "cpu_mpi" if use_mpi else "serial"
    else:
        backend = str(cmake.get("FFT_BACKEND", nested(profile, "backend", "serial")))
    return profile, use_mpi, use_openmp, backend


def _validate_solver_selection(
    case: dict[str, Any], manifest: dict[str, Any], profile_name: str
) -> None:
    configured_profile = nested(case, "solver.profile")
    if configured_profile not in {None, ""} and str(configured_profile) != profile_name:
        raise CaseInputError(
            f"case solver.profile {configured_profile!r} does not match requested "
            f"profile {profile_name!r}; change the environment design and regenerate "
            "the execution environment"
        )

    profile, use_mpi, profile_openmp, backend = _profile_settings(
        manifest, profile_name
    )
    use_cuda = backend in {"cuda", "cufft", "cufftmp"}
    solver = _mapping(nested(case, "solver", {}), "case solver")
    if "processes" in solver:
        raise CaseInputError(
            "solver.processes is no longer supported; use solver.mpi_processes"
        )
    try:
        processes = int(solver.get("mpi_processes", 1))
    except (TypeError, ValueError) as exc:
        raise CaseInputError("case solver.mpi_processes must be an integer") from exc
    if processes < 1:
        raise CaseInputError("case solver.mpi_processes must be positive")
    if not use_mpi and processes != 1:
        raise CaseInputError(
            f"profile {profile_name!r} does not use MPI, so solver.mpi_processes "
            "must be 1"
        )
    if (
        str(manifest.get("model", "")).lower() == "nse"
        and use_mpi
        and processes < 4
    ):
        raise CaseInputError(
            "the current NSE y-z decomposition requires at least 4 MPI processes"
        )
    requested_mpi = solver.get("use_mpi", use_mpi)
    requested_openmp = solver.get("use_openmp", False)
    requested_cuda = solver.get("use_cuda", use_cuda)
    for name, value in (
        ("use_mpi", requested_mpi),
        ("use_openmp", requested_openmp),
        ("use_cuda", requested_cuda),
    ):
        if not isinstance(value, bool):
            raise CaseInputError(f"case solver.{name} must be true or false")
    if requested_mpi != use_mpi:
        raise CaseInputError(
            f"case solver.use_mpi does not match profile {profile_name!r}; "
            "regenerate the execution environment"
        )
    if requested_cuda != use_cuda:
        raise CaseInputError(
            f"case solver.use_cuda does not match profile {profile_name!r}; "
            "regenerate the execution environment"
        )
    if requested_openmp and not profile_openmp:
        raise CaseInputError(
            f"profile {profile_name!r} was built without OpenMP support"
        )
    omp_threads = nested(case, "solver.omp_threads", 1)
    if not isinstance(omp_threads, int) or isinstance(omp_threads, bool) or omp_threads < 1:
        raise CaseInputError("case solver.omp_threads must be a positive integer")


def _common_values(
    case: dict[str, Any],
    equation: str,
    use_mpi: bool,
    use_openmp: bool,
    backend: str,
) -> list[tuple[str, Any]]:
    processes = int(nested(case, "solver.mpi_processes", 1))
    initial_condition = str(_required(case, "flow.type")).strip()
    if equation.upper() == "NSE":
        initial_condition = {
            "tgv": "taylor_green",
            "taylor-green": "taylor_green",
            "taylor_green_vortex": "taylor_green",
            "hit": "hit_spectral",
            "homogeneous_isotropic_turbulence": "hit_spectral",
            "homogeneous-isotropic-turbulence": "hit_spectral",
        }.get(initial_condition.lower(), initial_condition.lower())

    values = [
        ("equation", equation),
        ("case_name", str(_required(case, "case_id"))),
        ("initial_condition", initial_condition),
        ("nx", _required(case, "grid.nx")),
        ("ny", _required(case, "grid.ny")),
        ("nz", _required(case, "grid.nz")),
        ("nghost", nested(case, "grid.nghost")),
        ("x_min", _required(case, "grid.x_min")),
        ("x_max", _required(case, "grid.x_max")),
        ("y_min", _required(case, "grid.y_min")),
        ("y_max", _required(case, "grid.y_max")),
        ("z_min", _required(case, "grid.z_min")),
        ("z_max", _required(case, "grid.z_max")),
        ("dt", _required(case, "time.dt")),
        ("t_max", nested(case, "time.t_max", 0.0)),
        ("nsteps", _required(case, "time.nsteps")),
        ("cfl", nested(case, "time.cfl", 0.0)),
        ("use_fixed_dt", nested(case, "time.use_fixed_dt", True)),
        ("output_frequency", _required(case, "time.output_frequency")),
        ("output_dir", nested(case, "output.directory", "output")),
        ("output_format", nested(case, "output.format", "slf")),
        ("precision_name", nested(case, "output.precision", "float64")),
        ("write_initial", nested(case, "output.write_initial", True)),
        ("write_meta", nested(case, "output.write_meta", True)),
        ("backend", backend),
        ("use_mpi", use_mpi),
        ("use_openmp", nested(case, "solver.use_openmp", use_openmp)),
        ("rank", 0),
        ("nprocs", processes),
    ]
    return values


def render_nse(
    case: dict[str, Any], manifest: dict[str, Any], profile_name: str
) -> str:
    profile, use_mpi, use_openmp, backend = _profile_settings(
        manifest, profile_name
    )
    nse = dict(_mapping(nested(case, "physics.nse", {}), "physics.nse"))
    numerics = _mapping(nested(case, "numerics", {}), "numerics")
    if "convective_order" in nse or "convective_order" in numerics:
        raise CaseInputError(
            "convective_order is no longer supported; set "
            "numerics.convective_scheme to KEEP2 or KEEP6"
        )
    for source, target in NSE_ALIASES.items():
        if target not in nse and source in nse:
            nse[target] = nse[source]
    numerical_aliases = {
        "convective_scheme": (
            numerics.get("convective_scheme"),
            numerics.get("flux"),
        ),
        "viscous_scheme": (numerics.get("viscous_scheme"),),
        "boundary_condition": (numerics.get("boundary_condition"),),
        "time_integrator": (
            numerics.get("time_integrator"),
            numerics.get("time_integration"),
        ),
    }
    for target, candidates in numerical_aliases.items():
        if target not in nse or nse[target] in {None, ""}:
            for value in candidates:
                if value not in {None, ""}:
                    nse[target] = value
                    break
    if "cfl" not in nse:
        nse["cfl"] = nested(case, "time.cfl")
    hit = _mapping(nested(case, "flow.hit", {}), "flow.hit")
    hit_aliases = {
        "hit_spectrum": ("spectrum",),
        "hit_seed": ("random_seed",),
        "hit_rms_velocity": ("rms_velocity", "turbulent_mach_number"),
        "hit_peak_wavenumber": ("peak_wavenumber",),
        "hit_integral_length": ("integral_length",),
        "hit_kolmogorov_length": ("kolmogorov_length",),
        "hit_dealias_fraction": ("dealias_fraction",),
    }
    for target, candidates in hit_aliases.items():
        if target not in nse or nse[target] in {None, ""}:
            for source in candidates:
                value = hit.get(source)
                if value not in {None, ""}:
                    nse[target] = value
                    break
    forcing = _mapping(nested(case, "forcing", {}), "forcing")
    forcing_aliases = {
        "forcing_scheme": ("type", "scheme"),
        "forcing_spectrum": ("spectrum",),
        "forcing_fft_backend": ("fft_backend",),
        "forcing_k_cutoff": ("k_cutoff",),
        "forcing_target_dissipation": ("target_dissipation",),
        "forcing_dilatational_ratio": ("dilatational_ratio",),
        "forcing_denominator_floor": ("denominator_floor",),
        "forcing_max_coefficient": ("max_coefficient",),
        "forcing_report_interval": ("report_interval",),
    }
    for target, candidates in forcing_aliases.items():
        if target not in nse or nse[target] in {None, ""}:
            for source in candidates:
                value = forcing.get(source)
                if value not in {None, ""}:
                    nse[target] = value
                    break
    for key in (
        "convective_scheme",
        "viscous_scheme",
        "boundary_condition",
        "time_integrator",
        "hit_spectrum",
        "forcing_scheme",
        "forcing_spectrum",
        "forcing_fft_backend",
    ):
        if isinstance(nse.get(key), str):
            nse[key] = _canonical_selector(nse[key])
    forcing_backend = str(
        _mapping(profile.get("cmake", {}), f"profile {profile_name}.cmake").get(
            "NSE_FORCING_FFT_BACKEND", "none"
        )
    ).strip().lower()
    forcing_scheme = str(nse.get("forcing_scheme", "none"))
    requested_forcing_backend = str(nse.get("forcing_fft_backend", "auto"))
    if forcing_scheme != "none":
        if forcing_backend == "none":
            raise CaseInputError(
                f"forcing.type={forcing_scheme!r} requires a forcing FFT backend; "
                f"profile {profile_name!r} has none. Use "
                "solver.profile=cpu_mpi_2decomp_fftw for CPU/MPI or "
                "solver.profile=cuda_single for a single GPU, then regenerate "
                "the execution environment"
            )
        if requested_forcing_backend not in {"auto", forcing_backend}:
            raise CaseInputError(
                f"forcing.fft_backend={requested_forcing_backend!r} does not match "
                f"profile {profile_name!r} backend {forcing_backend!r}"
            )
    if nse.get("convective_scheme") == "keep":
        raise CaseInputError(
            "convective_scheme=KEEP is no longer supported; use KEEP2 or KEEP6"
        )
    if nse.get("time_integrator") in {"rk3", "ssp_rk3", "ssp-rk3"}:
        nse["time_integrator"] = "ssprk3"
    lines = [
        "! Automatically generated from case.yaml.",
        "! Edit case.yaml and regenerate this file.",
        "",
        "&simulation",
    ]
    common = _common_values(case, "NSE", use_mpi, use_openmp, backend)
    common.append(("cuda_device", nested(case, "solver.cuda_device", 0)))
    _append(lines, common)
    lines.extend(["/", "", "&nse"])
    _append(lines, [(key, nse.get(key)) for key in NSE_KEYS])
    lines.extend(["/", ""])
    return "\n".join(lines)


def render_gpe(
    case: dict[str, Any], manifest: dict[str, Any], profile_name: str
) -> str:
    _, use_mpi, use_openmp, backend = _profile_settings(manifest, profile_name)
    gpe = dict(_mapping(nested(case, "physics.gpe", {}), "physics.gpe"))
    for source, target in GPE_ALIASES.items():
        if target not in gpe and source in gpe:
            gpe[target] = gpe[source]
    common = _common_values(case, "GPE", use_mpi, use_openmp, backend)
    common.extend(
        [
            ("restart_file", nested(case, "restart.file")),
            ("timing_enabled", nested(case, "output.timing_enabled", False)),
            ("use_cuda", backend in {"cufft", "cufftmp"}),
            ("cuda_device", nested(case, "solver.cuda_device", 0)),
        ]
    )
    lines = [
        "! Automatically generated from case.yaml.",
        "! Edit case.yaml and regenerate this file.",
        "",
        "&simulation",
    ]
    _append(lines, common)
    lines.extend(["/", "", "&gpe"])
    _append(lines, [(key, gpe.get(key)) for key in GPE_KEYS])
    lines.extend(["/", ""])
    return "\n".join(lines)


def render_case_input(
    case_path: Path,
    manifest_path: Path,
    model: str,
    profile_name: str,
) -> tuple[str, str]:
    case = _mapping(load_yaml(case_path), "case YAML")
    manifest = _mapping(load_yaml(manifest_path), "solver manifest")
    actual_model = str(nested(case, "physics.model", model)).strip().lower()
    allowed_models = {model.lower()}
    if model.lower() == "gpe":
        allowed_models.update({"gross_pitaevskii", "gross-pitaevskii"})
    if actual_model not in allowed_models:
        raise CaseInputError(
            f"case physics.model {actual_model!r} does not match selected model {model!r}"
        )
    _validate_solver_selection(case, manifest, profile_name)
    input_cfg = _mapping(manifest.get("input"), "solver manifest.input")
    input_name = str(input_cfg.get("default_name") or "input.dat")
    if model.lower() == "nse":
        return input_name, render_nse(case, manifest, profile_name)
    if model.lower() == "gpe":
        return input_name, render_gpe(case, manifest, profile_name)
    raise CaseInputError(f"unsupported model: {model}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a solver input namelist from case.yaml."
    )
    parser.add_argument("--case", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--model", required=True, choices=["nse", "gpe"])
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        case_path = Path(args.case).resolve()
        input_name, text = render_case_input(
            case_path,
            Path(args.manifest).resolve(),
            args.model,
            args.profile,
        )
        output = Path(args.output).resolve() if args.output else case_path.parent / input_name
        if output.exists() and not args.overwrite:
            raise CaseInputError(f"output already exists; use --overwrite: {output}")
        output.write_text(text, encoding="ascii")
        print(f"[OK] Generated solver input: {output}")
        return 0
    except (CaseInputError, YamlFormatError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
