#!/usr/bin/env python3
"""Render model-specific Fortran namelists from a framework case.yaml."""

from __future__ import annotations

import argparse
import math
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
    "hybrid_smooth_scheme",
    "hybrid_shock_scheme",
    "hybrid_sensor",
    "hybrid_sensor_onset",
    "hybrid_sensor_full",
    "viscous_scheme",
    "boundary_condition",
    "time_integrator",
    "hit_spectrum",
    "hit_seed",
    "hit_turbulent_mach",
    "hit_turbulent_reynolds",
    "hit_rms_velocity",
    "hit_peak_wavenumber",
    "hit_integral_length",
    "hit_kolmogorov_length",
    "hit_johnsen_length_scale_ratio",
    "hit_pope_energy_constant",
    "hit_pope_large_scale_constant",
    "hit_pope_dissipation_constant",
    "hit_pope_large_scale_exponent",
    "hit_pope_dissipation_exponent",
    "hit_dealias_fraction",
    "hit_isotropy_mode",
    "hit_isotropy_k_cutoff",
    "hit_isotropy_tolerance",
    "hit_isotropy_max_iterations",
    "imported_turbulence_file",
    "imported_turbulence_mode",
    "imported_turbulence_x_start",
    "imported_turbulence_blend_cells",
    "imported_turbulence_velocity_offset_x",
    "imported_turbulence_velocity_offset_y",
    "imported_turbulence_velocity_offset_z",
    "imported_turbulence_background_rho",
    "imported_turbulence_background_u",
    "imported_turbulence_background_v",
    "imported_turbulence_background_w",
    "imported_turbulence_background_p",
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

NSE_HIT_SPECTRUM_SCHEMAS: dict[str, dict[str, Any]] = {
    "johnsen": {
        "parameters": ("characteristic_length", "length_scale_ratio"),
        "required": ("characteristic_length", "length_scale_ratio"),
        "namelist": {
            "characteristic_length": "hit_integral_length",
            "length_scale_ratio": "hit_johnsen_length_scale_ratio",
        },
    },
    "pope": {
        "parameters": (
            "integral_length",
            "energy_constant",
            "large_scale_constant",
            "dissipation_constant",
            "large_scale_exponent",
            "dissipation_exponent",
        ),
        "required": ("integral_length",),
        "namelist": {
            "integral_length": "hit_integral_length",
            "energy_constant": "hit_pope_energy_constant",
            "large_scale_constant": "hit_pope_large_scale_constant",
            "dissipation_constant": "hit_pope_dissipation_constant",
            "large_scale_exponent": "hit_pope_large_scale_exponent",
            "dissipation_exponent": "hit_pope_dissipation_exponent",
        },
    },
}
NSE_HIT_SPECTRUM_TYPES = set(NSE_HIT_SPECTRUM_SCHEMAS)

NSE_FORCING_SCHEMAS: dict[str, dict[str, Any]] = {
    "petersen_livescu": {
        "parameters": (
            "spectrum",
            "fft_backend",
            "k_cutoff",
            "target_dissipation",
            "dilatational_ratio",
            "denominator_floor",
            "max_coefficient",
            "report_interval",
        ),
        "required": ("target_dissipation",),
        "selectors": {
            "spectrum": {
                "allowed": ("full_spectrum", "low_wavenumber"),
                "aliases": {"full_wavenumber": "full_spectrum"},
            },
            "fft_backend": {
                "allowed": ("auto", "2decomp_fftw", "cufft"),
                "aliases": {},
            },
        },
        "namelist": {
            "spectrum": "forcing_spectrum",
            "fft_backend": "forcing_fft_backend",
            "k_cutoff": "forcing_k_cutoff",
            "target_dissipation": "forcing_target_dissipation",
            "dilatational_ratio": "forcing_dilatational_ratio",
            "denominator_floor": "forcing_denominator_floor",
            "max_coefficient": "forcing_max_coefficient",
            "report_interval": "forcing_report_interval",
        },
    }
}
NSE_FORCING_TYPES = {"none", *NSE_FORCING_SCHEMAS}

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


def _positive_float(value: Any, label: str, *, allow_zero: bool = False) -> float:
    if isinstance(value, bool):
        raise CaseInputError(f"{label} must be a number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise CaseInputError(f"{label} must be a number") from exc
    if not math.isfinite(number):
        raise CaseInputError(f"{label} must be finite")
    outside_range = number < 0.0 if allow_zero else number <= 0.0
    if outside_range:
        comparison = "non-negative" if allow_zero else "positive"
        raise CaseInputError(f"{label} must be {comparison}")
    return number


def _finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise CaseInputError(f"{label} must be a number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise CaseInputError(f"{label} must be a number") from exc
    if not math.isfinite(number):
        raise CaseInputError(f"{label} must be finite")
    return number


def _vector3(value: Any, label: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise CaseInputError(f"{label} must contain exactly three numbers")
    return tuple(
        _finite_float(component, f"{label}[{index}]")
        for index, component in enumerate(value)
    )


def derive_nse_hit_transport(
    turbulent_mach: Any,
    turbulent_reynolds: Any,
    characteristic_length: Any,
) -> dict[str, float]:
    """Resolve the original HIT Re_lambda scaling in solver units."""

    mach = _positive_float(turbulent_mach, "HIT turbulent Mach number")
    re_lambda = _positive_float(
        turbulent_reynolds, "HIT turbulent Reynolds number (Re_lambda)"
    )
    length = _positive_float(characteristic_length, "HIT characteristic length")
    component_rms = mach / (3.0**0.5)
    integral_reynolds = 3.0 * re_lambda * re_lambda / 20.0
    velocity_scale = (1.5**0.5) * component_rms
    kinematic_viscosity = velocity_scale * length / integral_reynolds
    return {
        "turbulent_mach": mach,
        "turbulent_reynolds": re_lambda,
        "component_rms": component_rms,
        "integral_reynolds": integral_reynolds,
        "kinematic_viscosity": kinematic_viscosity,
        "solver_reynolds": 1.0 / kinematic_viscosity,
        "taylor_microscale": length * (10.0 / integral_reynolds) ** 0.5,
        "kolmogorov_length": length * integral_reynolds ** (-0.75),
    }


def _resolve_legacy_nse_hit(hit: dict[str, Any]) -> dict[str, Any]:
    values: dict[str, Any] = {}
    aliases = {
        "hit_spectrum": ("spectrum",),
        "hit_seed": ("random_seed",),
        "hit_rms_velocity": ("rms_velocity",),
        "hit_peak_wavenumber": ("peak_wavenumber",),
        "hit_integral_length": ("integral_length",),
        "hit_kolmogorov_length": ("kolmogorov_length",),
        "hit_dealias_fraction": ("dealias_fraction",),
        "hit_isotropy_mode": ("isotropy_mode",),
        "hit_isotropy_k_cutoff": ("isotropy_k_cutoff",),
        "hit_isotropy_tolerance": ("isotropy_tolerance",),
        "hit_isotropy_max_iterations": ("isotropy_max_iterations",),
    }
    for target, candidates in aliases.items():
        for source in candidates:
            value = hit.get(source)
            if value not in {None, ""}:
                values[target] = value
                break

    turbulent_mach = hit.get("turbulent_mach_number")
    turbulent_reynolds = hit.get(
        "turbulent_reynolds_number", hit.get("taylor_reynolds_number")
    )
    if turbulent_mach not in {None, ""} and turbulent_reynolds not in {None, ""}:
        length = hit.get("integral_length")
        if length in {None, ""}:
            raise CaseInputError(
                "legacy HIT targets require flow.hit.integral_length"
            )
        derived = derive_nse_hit_transport(
            turbulent_mach, turbulent_reynolds, length
        )
        values.update(
            {
                "mach": derived["turbulent_mach"],
                "reynolds": derived["solver_reynolds"],
                "hit_turbulent_mach": derived["turbulent_mach"],
                "hit_turbulent_reynolds": derived["turbulent_reynolds"],
                "hit_rms_velocity": derived["component_rms"],
                "hit_kolmogorov_length": derived["kolmogorov_length"],
            }
        )
    elif turbulent_mach not in {None, ""}:
        values["mach"] = _positive_float(
            turbulent_mach, "flow.hit.turbulent_mach_number"
        )
    elif turbulent_reynolds not in {None, ""}:
        raise CaseInputError(
            "flow.hit.turbulent_reynolds_number also requires "
            "flow.hit.turbulent_mach_number"
        )
    return values


def _resolve_nse_hit(case: dict[str, Any]) -> dict[str, Any]:
    flow_type = _canonical_selector(str(nested(case, "flow.type", "")))
    if flow_type not in {
        "hit",
        "hit_spectral",
        "homogeneous_isotropic_turbulence",
    }:
        return {}

    hit = _mapping(nested(case, "flow.hit", {}), "flow.hit")
    spectrum_value = hit.get("spectrum")
    if not isinstance(spectrum_value, dict):
        return _resolve_legacy_nse_hit(hit)

    common_keys = {
        "turbulent_mach_number",
        "turbulent_reynolds_number",
        "taylor_reynolds_number",
        "random_seed",
        "dealias_fraction",
        "isotropy_mode",
        "isotropy_k_cutoff",
        "isotropy_tolerance",
        "isotropy_max_iterations",
    }
    unknown_hit = sorted(set(hit) - common_keys - {"spectrum"})
    if unknown_hit:
        raise CaseInputError(
            f"unknown flow.hit key(s): {', '.join(unknown_hit)}"
        )

    turbulent_mach = hit.get("turbulent_mach_number")
    turbulent_reynolds = hit.get("turbulent_reynolds_number")
    taylor_reynolds = hit.get("taylor_reynolds_number")
    if turbulent_reynolds not in {None, ""} and taylor_reynolds not in {None, ""}:
        turbulent_value = _positive_float(
            turbulent_reynolds, "flow.hit.turbulent_reynolds_number"
        )
        taylor_value = _positive_float(
            taylor_reynolds, "flow.hit.taylor_reynolds_number"
        )
        if turbulent_value != taylor_value:
            raise CaseInputError(
                "flow.hit turbulent_reynolds_number and "
                "taylor_reynolds_number disagree"
            )
        turbulent_reynolds = turbulent_value
    if turbulent_reynolds in {None, ""}:
        turbulent_reynolds = taylor_reynolds
    if turbulent_mach in {None, ""} or turbulent_reynolds in {None, ""}:
        raise CaseInputError(
            "flow.hit requires turbulent_mach_number and "
            "turbulent_reynolds_number"
        )

    spectrum = _mapping(spectrum_value, "flow.hit.spectrum")
    allowed_spectrum_keys = {"type", *NSE_HIT_SPECTRUM_TYPES}
    unknown_spectrum = sorted(set(spectrum) - allowed_spectrum_keys)
    if unknown_spectrum:
        raise CaseInputError(
            "unknown flow.hit.spectrum key(s): "
            f"{', '.join(unknown_spectrum)}"
        )
    raw_type = spectrum.get("type")
    if raw_type in {None, ""}:
        raise CaseInputError("flow.hit.spectrum.type is required")
    spectrum_type = _canonical_selector(str(raw_type))
    spectrum_type = {"k4_gaussian": "johnsen"}.get(
        spectrum_type, spectrum_type
    )
    if spectrum_type not in NSE_HIT_SPECTRUM_TYPES:
        raise CaseInputError(
            f"unknown flow.hit.spectrum.type={raw_type!r}; available: "
            f"{sorted(NSE_HIT_SPECTRUM_TYPES)}"
        )
    if spectrum_type not in spectrum:
        raise CaseInputError(
            f"flow.hit.spectrum.type={spectrum_type!r} requires the "
            f"flow.hit.spectrum.{spectrum_type} mapping"
        )

    schema = NSE_HIT_SPECTRUM_SCHEMAS[spectrum_type]
    selected = _mapping(
        spectrum[spectrum_type], f"flow.hit.spectrum.{spectrum_type}"
    )
    parameter_keys = set(schema["parameters"])
    unknown_selected = sorted(set(selected) - parameter_keys)
    if unknown_selected:
        raise CaseInputError(
            f"unknown flow.hit.spectrum.{spectrum_type} key(s): "
            f"{', '.join(unknown_selected)}; available: "
            f"{sorted(parameter_keys)}"
        )
    for required in schema["required"]:
        if selected.get(required) in {None, ""}:
            raise CaseInputError(
                f"flow.hit.spectrum.{spectrum_type}.{required} is required"
            )

    length_key = "characteristic_length"
    if spectrum_type == "pope":
        length_key = "integral_length"
    length = _positive_float(
        selected[length_key],
        f"flow.hit.spectrum.{spectrum_type}.{length_key}",
    )
    derived = derive_nse_hit_transport(
        turbulent_mach, turbulent_reynolds, length
    )
    values: dict[str, Any] = {
        "hit_spectrum": spectrum_type,
        "hit_turbulent_mach": derived["turbulent_mach"],
        "hit_turbulent_reynolds": derived["turbulent_reynolds"],
        "hit_rms_velocity": derived["component_rms"],
        "hit_integral_length": length,
        "hit_kolmogorov_length": derived["kolmogorov_length"],
        "mach": derived["turbulent_mach"],
        "reynolds": derived["solver_reynolds"],
    }
    common_namelist = {
        "random_seed": "hit_seed",
        "dealias_fraction": "hit_dealias_fraction",
        "isotropy_mode": "hit_isotropy_mode",
        "isotropy_k_cutoff": "hit_isotropy_k_cutoff",
        "isotropy_tolerance": "hit_isotropy_tolerance",
        "isotropy_max_iterations": "hit_isotropy_max_iterations",
    }
    for source, target in common_namelist.items():
        if hit.get(source) not in {None, ""}:
            values[target] = hit[source]
    for source, target in schema["namelist"].items():
        if selected.get(source) not in {None, ""}:
            values[target] = selected[source]
    values["hit_integral_length"] = length

    if spectrum_type == "johnsen":
        ratio = _positive_float(
            selected["length_scale_ratio"],
            "flow.hit.spectrum.johnsen.length_scale_ratio",
        )
        values["hit_johnsen_length_scale_ratio"] = ratio
        values["hit_peak_wavenumber"] = 2.0 * ratio / length
    else:
        positive_parameters = {
            "energy_constant",
            "large_scale_constant",
            "dissipation_constant",
            "dissipation_exponent",
        }
        for name in positive_parameters.intersection(selected):
            normalized = _positive_float(
                selected[name], f"flow.hit.spectrum.pope.{name}"
            )
            values[schema["namelist"][name]] = normalized
        if "large_scale_exponent" in selected:
            normalized = _positive_float(
                selected["large_scale_exponent"],
                "flow.hit.spectrum.pope.large_scale_exponent",
                allow_zero=True,
            )
            values["hit_pope_large_scale_exponent"] = normalized
    return values


def _resolve_nse_imported_turbulence(
    case: dict[str, Any],
    *,
    case_dir: Path | None = None,
    runtime_root: Path | None = None,
) -> dict[str, Any]:
    flow_type = _canonical_selector(str(nested(case, "flow.type", "")))
    aliases = {
        "imported_turbulence",
        "turbulence_import",
        "turbulence_embed",
        "turbulence_tile",
    }
    if flow_type not in aliases:
        return {}

    imported = _mapping(
        nested(case, "flow.imported_turbulence", {}),
        "flow.imported_turbulence",
    )
    allowed = {
        "file",
        "mode",
        "x_start",
        "blend_cells",
        "velocity_offset",
        "background",
    }
    unknown = sorted(set(imported) - allowed)
    if unknown:
        raise CaseInputError(
            "unknown flow.imported_turbulence key(s): " + ", ".join(unknown)
        )

    source_file = imported.get("file")
    if not isinstance(source_file, str) or not source_file.strip():
        raise CaseInputError("flow.imported_turbulence.file is required")

    implied_mode = "tile" if flow_type == "turbulence_tile" else "embed"
    mode = _canonical_selector(str(imported.get("mode", implied_mode)))
    if mode not in {"embed", "tile"}:
        raise CaseInputError(
            "flow.imported_turbulence.mode must be EMBED or TILE"
        )

    blend_cells = imported.get("blend_cells", 0)
    if (
        not isinstance(blend_cells, int)
        or isinstance(blend_cells, bool)
        or blend_cells < 0
    ):
        raise CaseInputError(
            "flow.imported_turbulence.blend_cells must be a non-negative integer"
        )
    if mode == "tile" and blend_cells != 0:
        raise CaseInputError(
            "flow.imported_turbulence tile mode requires blend_cells: 0"
        )

    velocity_offset = _vector3(
        imported.get("velocity_offset", [0.0, 0.0, 0.0]),
        "flow.imported_turbulence.velocity_offset",
    )
    background = _mapping(
        imported.get("background", {}),
        "flow.imported_turbulence.background",
    )
    allowed_background = {"density", "velocity", "pressure"}
    unknown_background = sorted(set(background) - allowed_background)
    if unknown_background:
        raise CaseInputError(
            "unknown flow.imported_turbulence.background key(s): "
            + ", ".join(unknown_background)
        )

    gamma = _positive_float(
        nested(case, "physics.nse.gamma", 1.4), "physics.nse.gamma"
    )
    if gamma <= 1.0:
        raise CaseInputError("physics.nse.gamma must be greater than 1")
    rho0 = _positive_float(
        nested(case, "physics.nse.rho0", 1.0), "physics.nse.rho0"
    )
    background_rho = _positive_float(
        background.get("density", rho0),
        "flow.imported_turbulence.background.density",
    )
    background_pressure = _positive_float(
        background.get("pressure", 1.0 / gamma),
        "flow.imported_turbulence.background.pressure",
    )
    background_velocity = _vector3(
        background.get("velocity", [0.0, 0.0, 0.0]),
        "flow.imported_turbulence.background.velocity",
    )
    x_start = _finite_float(
        imported.get("x_start", nested(case, "grid.x_min", 0.0)),
        "flow.imported_turbulence.x_start",
    )

    source_file = source_file.strip()
    if case_dir is not None and runtime_root is not None:
        source_path = Path(source_file)
        if not source_path.is_absolute():
            source_path = (case_dir / source_path).resolve()
            try:
                source_file = source_path.relative_to(
                    runtime_root.resolve()
                ).as_posix()
            except ValueError:
                source_file = source_path.as_posix()
        else:
            source_file = source_path.as_posix()

    return {
        "imported_turbulence_file": source_file,
        "imported_turbulence_mode": mode,
        "imported_turbulence_x_start": x_start,
        "imported_turbulence_blend_cells": blend_cells,
        "imported_turbulence_velocity_offset_x": velocity_offset[0],
        "imported_turbulence_velocity_offset_y": velocity_offset[1],
        "imported_turbulence_velocity_offset_z": velocity_offset[2],
        "imported_turbulence_background_rho": background_rho,
        "imported_turbulence_background_u": background_velocity[0],
        "imported_turbulence_background_v": background_velocity[1],
        "imported_turbulence_background_w": background_velocity[2],
        "imported_turbulence_background_p": background_pressure,
    }


def _resolve_nse_forcing(case: dict[str, Any]) -> dict[str, Any]:
    forcing = _mapping(nested(case, "forcing", {}), "forcing")
    legacy_parameters = set(
        NSE_FORCING_SCHEMAS["petersen_livescu"]["parameters"]
    )
    allowed_root = {
        "type",
        "scheme",
        *NSE_FORCING_SCHEMAS,
        *legacy_parameters,
    }
    unknown_root = sorted(set(forcing) - allowed_root)
    if unknown_root:
        raise CaseInputError(
            "unknown forcing key(s): "
            f"{', '.join(unknown_root)}; available root keys are "
            "type and petersen_livescu"
        )

    if "type" in forcing and "scheme" in forcing:
        type_value = _canonical_selector(str(forcing["type"]))
        scheme_value = _canonical_selector(str(forcing["scheme"]))
        if type_value != scheme_value:
            raise CaseInputError(
                "forcing.type and legacy forcing.scheme select different types"
            )
    raw_type = forcing.get("type", forcing.get("scheme", "none"))
    forcing_type = _canonical_selector(str(raw_type))
    if forcing_type not in NSE_FORCING_TYPES:
        raise CaseInputError(
            f"unknown forcing.type={raw_type!r}; available: "
            f"{sorted(NSE_FORCING_TYPES)}"
        )

    values: dict[str, Any] = {"forcing_scheme": forcing_type}
    if forcing_type == "none":
        return values

    schema = NSE_FORCING_SCHEMAS[forcing_type]
    parameter_keys = set(schema["parameters"])
    legacy_keys = sorted(legacy_parameters.intersection(forcing))
    if forcing_type in forcing:
        if legacy_keys:
            raise CaseInputError(
                f"forcing mixes the new forcing.{forcing_type} mapping with "
                f"legacy flat key(s): {', '.join(legacy_keys)}"
            )
        selected = dict(
            _mapping(forcing[forcing_type], f"forcing.{forcing_type}")
        )
        unknown_selected = sorted(set(selected) - parameter_keys)
        if unknown_selected:
            raise CaseInputError(
                f"unknown forcing.{forcing_type} key(s): "
                f"{', '.join(unknown_selected)}; available: "
                f"{sorted(parameter_keys)}"
            )
    elif forcing_type == "petersen_livescu" and legacy_keys:
        selected = dict(forcing)
    else:
        raise CaseInputError(
            f"forcing.type={forcing_type!r} requires the "
            f"forcing.{forcing_type} mapping"
        )

    for required in schema["required"]:
        if selected.get(required) in {None, ""}:
            raise CaseInputError(
                f"forcing.{forcing_type}.{required} is required"
            )

    for parameter, selector in schema.get("selectors", {}).items():
        value = selected.get(parameter)
        if value in {None, ""}:
            continue
        normalized = _canonical_selector(str(value))
        normalized = selector.get("aliases", {}).get(normalized, normalized)
        if normalized not in selector["allowed"]:
            raise CaseInputError(
                f"forcing.{forcing_type}.{parameter}={value!r} is invalid; "
                f"available: {list(selector['allowed'])}"
            )
        selected[parameter] = normalized

    for source, target in schema["namelist"].items():
        value = selected.get(source)
        if value not in {None, ""}:
            values[target] = value
    return values


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
        if model == "nse":
            backend = "cuda_mpi" if use_mpi else "cuda"
        else:
            backend = "cufft"
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
    use_cuda = backend in {"cuda", "cuda_mpi", "cufft", "cufftmp"}
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
        and processes < 2
    ):
        raise CaseInputError(
            "the NSE y-z decomposition requires at least 2 MPI processes"
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
        initial_condition = _canonical_selector(initial_condition)
        initial_condition = {
            "tgv": "taylor_green",
            "taylor_green_vortex": "taylor_green",
            "hit": "hit_spectral",
            "homogeneous_isotropic_turbulence": "hit_spectral",
            "turbulence_import": "imported_turbulence",
            "turbulence_embed": "imported_turbulence",
            "turbulence_tile": "imported_turbulence",
        }.get(initial_condition, initial_condition)

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
    case: dict[str, Any],
    manifest: dict[str, Any],
    profile_name: str,
    *,
    case_dir: Path | None = None,
    runtime_root: Path | None = None,
) -> str:
    profile, use_mpi, use_openmp, backend = _profile_settings(
        manifest, profile_name
    )
    nse = dict(_mapping(nested(case, "physics.nse", {}), "physics.nse"))
    numerics = _mapping(nested(case, "numerics", {}), "numerics")
    hybrid = _mapping(numerics.get("hybrid", {}), "numerics.hybrid")
    for legacy_key in ("flux", "reconstruction"):
        if legacy_key in numerics:
            raise CaseInputError(
                f"numerics.{legacy_key} is no longer supported; use "
                "numerics.convective_scheme"
            )
    if "convective_order" in nse or "convective_order" in numerics:
        raise CaseInputError(
            "convective_order is no longer supported; set "
            "numerics.convective_scheme to KEEP2, KEEP6, WENO5Z_ROE, "
            "or HYBRID"
        )
    for source, target in NSE_ALIASES.items():
        if target not in nse and source in nse:
            nse[target] = nse[source]
    numerical_aliases = {
        "convective_scheme": (
            numerics.get("convective_scheme"),
        ),
        "hybrid_smooth_scheme": (
            numerics.get("hybrid_smooth_scheme"),
            hybrid.get("smooth_scheme"),
        ),
        "hybrid_shock_scheme": (
            numerics.get("hybrid_shock_scheme"),
            hybrid.get("shock_scheme"),
        ),
        "hybrid_sensor": (
            numerics.get("hybrid_sensor"),
            hybrid.get("sensor"),
        ),
        "hybrid_sensor_onset": (
            numerics.get("hybrid_sensor_onset"),
            hybrid.get("sensor_onset"),
        ),
        "hybrid_sensor_full": (
            numerics.get("hybrid_sensor_full"),
            hybrid.get("sensor_full"),
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
    hit_values = _resolve_nse_hit(case)
    nse.update(hit_values)
    imported_turbulence_values = _resolve_nse_imported_turbulence(
        case, case_dir=case_dir, runtime_root=runtime_root
    )
    nse.update(imported_turbulence_values)
    forcing_values = _resolve_nse_forcing(case)
    for target, value in forcing_values.items():
        if target not in nse or nse[target] in {None, ""}:
            nse[target] = value
    for key in (
        "convective_scheme",
        "hybrid_smooth_scheme",
        "hybrid_shock_scheme",
        "hybrid_sensor",
        "viscous_scheme",
        "boundary_condition",
        "time_integrator",
        "hit_spectrum",
        "hit_isotropy_mode",
        "imported_turbulence_mode",
        "forcing_scheme",
        "forcing_spectrum",
        "forcing_fft_backend",
    ):
        if isinstance(nse.get(key), str):
            nse[key] = _canonical_selector(nse[key])
    profile_cmake = _mapping(
        profile.get("cmake", {}), f"profile {profile_name}.cmake"
    )
    initial_backend = str(
        profile_cmake.get("NSE_INIT_FFT_BACKEND", "none")
    ).strip().lower()
    flow_type = _canonical_selector(str(_required(case, "flow.type")))
    if flow_type in {
        "hit",
        "hit_spectral",
        "homogeneous_isotropic_turbulence",
    } and initial_backend == "none":
        raise CaseInputError(
            f"flow.type={flow_type!r} requires an initial-condition FFT backend; "
            f"profile {profile_name!r} has none. Remove an incompatible explicit "
            "solver.profile and regenerate the execution environment from the "
            "current FrameWork so run_case.py can select a compatible staged "
            "profile"
        )
    forcing_backend = str(
        profile_cmake.get("NSE_FORCING_FFT_BACKEND", "none")
    ).strip().lower()
    forcing_scheme = str(nse.get("forcing_scheme", "none"))
    requested_forcing_backend = str(nse.get("forcing_fft_backend", "auto"))
    if forcing_scheme != "none":
        if forcing_backend == "none":
            raise CaseInputError(
                f"forcing.type={forcing_scheme!r} requires a forcing FFT backend; "
                f"profile {profile_name!r} has none. Remove an incompatible "
                "explicit solver.profile and regenerate the execution environment "
                "from the current FrameWork so run_case.py can select a compatible "
                "staged profile"
            )
        if requested_forcing_backend not in {"auto", forcing_backend}:
            raise CaseInputError(
                "forcing.petersen_livescu.fft_backend="
                f"{requested_forcing_backend!r} does not match "
                f"profile {profile_name!r} backend {forcing_backend!r}"
            )
    if nse.get("convective_scheme") == "keep":
        raise CaseInputError(
            "convective_scheme=KEEP is no longer supported; use KEEP2 or KEEP6"
        )
    convective_scheme = str(nse.get("convective_scheme", "keep6"))
    allowed_convective_schemes = {"keep2", "keep6", "weno5z_roe", "hybrid"}
    if convective_scheme not in allowed_convective_schemes:
        raise CaseInputError(
            f"unsupported numerics.convective_scheme={convective_scheme!r}; "
            "use KEEP2, KEEP6, WENO5Z_ROE, or HYBRID"
        )
    if convective_scheme == "hybrid":
        hybrid_defaults = {
            "hybrid_smooth_scheme": "keep6",
            "hybrid_shock_scheme": "weno5z_roe",
            "hybrid_sensor": "ducros_pressure",
            "hybrid_sensor_onset": 0.01,
            "hybrid_sensor_full": 0.10,
        }
        for key, value in hybrid_defaults.items():
            if nse.get(key) in {None, ""}:
                nse[key] = value
        allowed_leaf_schemes = {"keep2", "keep6", "weno5z_roe"}
        for key in ("hybrid_smooth_scheme", "hybrid_shock_scheme"):
            if nse[key] not in allowed_leaf_schemes:
                raise CaseInputError(
                    f"unsupported numerics.hybrid.{key[len('hybrid_'):]}="
                    f"{nse[key]!r}; use KEEP2, KEEP6, or WENO5Z_ROE"
                )
        if nse["hybrid_sensor"] != "ducros_pressure":
            raise CaseInputError(
                f"unsupported numerics.hybrid.sensor={nse['hybrid_sensor']!r}; "
                "use DUCROS_PRESSURE"
            )
        onset = _positive_float(
            nse["hybrid_sensor_onset"],
            "numerics.hybrid.sensor_onset",
            allow_zero=True,
        )
        full = _positive_float(
            nse["hybrid_sensor_full"], "numerics.hybrid.sensor_full"
        )
        if full <= onset:
            raise CaseInputError(
                "numerics.hybrid requires 0 <= sensor_onset < sensor_full"
            )
        nse["hybrid_sensor_onset"] = onset
        nse["hybrid_sensor_full"] = full
    else:
        for key in (
            "hybrid_smooth_scheme",
            "hybrid_shock_scheme",
            "hybrid_sensor",
            "hybrid_sensor_onset",
            "hybrid_sensor_full",
        ):
            nse.pop(key, None)
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
        runtime_root = manifest_path.parent.parent.parent
        return input_name, render_nse(
            case,
            manifest,
            profile_name,
            case_dir=case_path.parent,
            runtime_root=runtime_root,
        )
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
