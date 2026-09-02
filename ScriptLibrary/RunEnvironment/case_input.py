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
    "boundary_x_min",
    "boundary_x_max",
    "boundary_y_min",
    "boundary_y_max",
    "boundary_z_min",
    "boundary_z_max",
    "boundary_x_min_reference_rho",
    "boundary_x_min_reference_u",
    "boundary_x_min_reference_v",
    "boundary_x_min_reference_w",
    "boundary_x_min_reference_p",
    "boundary_x_max_reference_rho",
    "boundary_x_max_reference_u",
    "boundary_x_max_reference_v",
    "boundary_x_max_reference_w",
    "boundary_x_max_reference_p",
    "boundary_y_min_reference_rho",
    "boundary_y_min_reference_u",
    "boundary_y_min_reference_v",
    "boundary_y_min_reference_w",
    "boundary_y_min_reference_p",
    "boundary_y_max_reference_rho",
    "boundary_y_max_reference_u",
    "boundary_y_max_reference_v",
    "boundary_y_max_reference_w",
    "boundary_y_max_reference_p",
    "boundary_z_min_reference_rho",
    "boundary_z_min_reference_u",
    "boundary_z_min_reference_v",
    "boundary_z_min_reference_w",
    "boundary_z_min_reference_p",
    "boundary_z_max_reference_rho",
    "boundary_z_max_reference_u",
    "boundary_z_max_reference_v",
    "boundary_z_max_reference_w",
    "boundary_z_max_reference_p",
    "boundary_relaxation_strength",
    "boundary_length_scale",
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
    "planar_shock_position",
    "planar_shock_direction",
    "planar_shock_mach",
    "planar_shock_upstream_rho",
    "planar_shock_upstream_u",
    "planar_shock_upstream_v",
    "planar_shock_upstream_w",
    "planar_shock_upstream_p",
    "planar_shock_downstream_rho",
    "planar_shock_downstream_u",
    "planar_shock_downstream_v",
    "planar_shock_downstream_w",
    "planar_shock_downstream_p",
    "shock_tube_diaphragm_position",
    "shock_tube_driver_rho",
    "shock_tube_driver_u",
    "shock_tube_driver_v",
    "shock_tube_driver_w",
    "shock_tube_driver_p",
    "shock_tube_driven_rho",
    "shock_tube_driven_u",
    "shock_tube_driven_v",
    "shock_tube_driven_w",
    "shock_tube_driven_p",
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

NSE_BOUNDARY_FACES = (
    "x_min",
    "x_max",
    "y_min",
    "y_max",
    "z_min",
    "z_max",
)

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
    if isinstance(value, (list, tuple)):
        if not value:
            raise CaseInputError("empty namelist sequences are not supported")
        return ", ".join(_fortran(item) for item in value)
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


def _nse_primitive_state(
    raw: Any, label: str
) -> tuple[float, tuple[float, float, float], float]:
    state = _mapping(raw, label)
    allowed = {"density", "velocity", "pressure"}
    unknown = sorted(set(state) - allowed)
    if unknown:
        raise CaseInputError(f"unknown {label} key(s): " + ", ".join(unknown))
    missing = sorted(allowed - set(state))
    if missing:
        raise CaseInputError(f"{label} is missing: " + ", ".join(missing))
    return (
        _positive_float(state["density"], f"{label}.density"),
        _vector3(state["velocity"], f"{label}.velocity"),
        _positive_float(state["pressure"], f"{label}.pressure"),
    )


def _resolve_nse_planar_shock(
    case: dict[str, Any],
) -> tuple[
    dict[str, Any],
    dict[str, tuple[float, tuple[float, float, float], float]],
]:
    flow_type = _canonical_selector(str(nested(case, "flow.type", "")))
    aliases = {
        "shock_turbulence_interaction",
        "planar_shock_turbulence",
        "shock_turbulence",
    }
    if flow_type not in aliases:
        return {}, {}

    shock = _mapping(
        nested(case, "flow.planar_shock", {}), "flow.planar_shock"
    )
    allowed = {
        "position",
        "propagation_direction",
        "upstream",
        "mach_number",
        "downstream",
    }
    unknown = sorted(set(shock) - allowed)
    if unknown:
        raise CaseInputError(
            "unknown flow.planar_shock key(s): " + ", ".join(unknown)
        )

    position = _finite_float(
        shock.get("position"), "flow.planar_shock.position"
    )
    x_min = _finite_float(_required(case, "grid.x_min"), "grid.x_min")
    x_max = _finite_float(_required(case, "grid.x_max"), "grid.x_max")
    nx = _required(case, "grid.nx")
    if not isinstance(nx, int) or isinstance(nx, bool) or nx <= 0:
        raise CaseInputError("grid.nx must be a positive integer")
    if not x_min < position < x_max:
        raise CaseInputError(
            "flow.planar_shock.position must lie strictly inside the x domain"
        )
    dx = (x_max - x_min) / nx
    face_index = round((position - x_min) / dx)
    aligned_position = x_min + face_index * dx
    if not math.isclose(position, aligned_position, rel_tol=1.0e-10, abs_tol=1.0e-12):
        raise CaseInputError(
            "flow.planar_shock.position must lie on a target x-cell boundary"
        )

    direction = _canonical_selector(
        str(shock.get("propagation_direction", "positive_x"))
    )
    direction = {
        "positive": "positive_x",
        "plus_x": "positive_x",
        "x_plus": "positive_x",
        "negative": "negative_x",
        "minus_x": "negative_x",
        "x_minus": "negative_x",
    }.get(direction, direction)
    if direction not in {"positive_x", "negative_x"}:
        raise CaseInputError(
            "flow.planar_shock.propagation_direction must be POSITIVE_X "
            "or NEGATIVE_X"
        )

    upstream = _nse_primitive_state(
        shock.get("upstream", {}), "flow.planar_shock.upstream"
    )
    mach_value = shock.get("mach_number")
    downstream_value = shock.get("downstream")
    has_mach = mach_value is not None and mach_value != ""
    has_downstream = downstream_value is not None and downstream_value != ""
    if has_mach == has_downstream:
        raise CaseInputError(
            "flow.planar_shock must specify exactly one of mach_number or downstream"
        )

    gamma = _positive_float(
        nested(case, "physics.nse.gamma", 1.4), "physics.nse.gamma"
    )
    if gamma <= 1.0:
        raise CaseInputError("physics.nse.gamma must be greater than 1")
    shock_mach = -1.0
    if has_downstream:
        downstream = _nse_primitive_state(
            shock["downstream"], "flow.planar_shock.downstream"
        )
    else:
        shock_mach = _positive_float(
            shock["mach_number"], "flow.planar_shock.mach_number"
        )
        if shock_mach <= 1.0:
            raise CaseInputError(
                "flow.planar_shock.mach_number must be greater than 1"
            )
        rho1, velocity1, pressure1 = upstream
        mach2 = shock_mach * shock_mach
        density_ratio = (
            (gamma + 1.0) * mach2
            / ((gamma - 1.0) * mach2 + 2.0)
        )
        pressure_ratio = 1.0 + (
            2.0 * gamma / (gamma + 1.0) * (mach2 - 1.0)
        )
        sound1 = math.sqrt(gamma * pressure1 / rho1)
        sign = 1.0 if direction == "positive_x" else -1.0
        velocity2 = list(velocity1)
        velocity2[0] += (
            sign * shock_mach * sound1 * (1.0 - 1.0 / density_ratio)
        )
        downstream = (
            rho1 * density_ratio,
            (velocity2[0], velocity2[1], velocity2[2]),
            pressure1 * pressure_ratio,
        )

    rho1, velocity1, pressure1 = upstream
    rho2, velocity2, pressure2 = downstream
    values = {
        "planar_shock_position": aligned_position,
        "planar_shock_direction": direction,
        "planar_shock_mach": shock_mach,
        "planar_shock_upstream_rho": rho1,
        "planar_shock_upstream_u": velocity1[0],
        "planar_shock_upstream_v": velocity1[1],
        "planar_shock_upstream_w": velocity1[2],
        "planar_shock_upstream_p": pressure1,
        "planar_shock_downstream_rho": rho2,
        "planar_shock_downstream_u": velocity2[0],
        "planar_shock_downstream_v": velocity2[1],
        "planar_shock_downstream_w": velocity2[2],
        "planar_shock_downstream_p": pressure2,
    }
    states = {
        "planar_shock.upstream": upstream,
        "planar_shock.downstream": downstream,
    }
    return values, states


def _resolve_nse_shock_tube(
    case: dict[str, Any],
) -> tuple[
    dict[str, Any],
    dict[str, tuple[float, tuple[float, float, float], float]],
]:
    flow_type = _canonical_selector(str(nested(case, "flow.type", "")))
    aliases = {
        "shock_tube_turbulence_interaction",
        "shock_tube_turbulence",
        "finite_driver_shock_turbulence",
    }
    if flow_type not in aliases:
        return {}, {}

    tube = _mapping(nested(case, "flow.shock_tube", {}), "flow.shock_tube")
    allowed = {"diaphragm_position", "driver", "driven"}
    unknown = sorted(set(tube) - allowed)
    if unknown:
        raise CaseInputError(
            "unknown flow.shock_tube key(s): " + ", ".join(unknown)
        )

    position = _finite_float(
        tube.get("diaphragm_position"), "flow.shock_tube.diaphragm_position"
    )
    x_min = _finite_float(_required(case, "grid.x_min"), "grid.x_min")
    x_max = _finite_float(_required(case, "grid.x_max"), "grid.x_max")
    nx = _required(case, "grid.nx")
    if not isinstance(nx, int) or isinstance(nx, bool) or nx <= 0:
        raise CaseInputError("grid.nx must be a positive integer")
    if not x_min < position < x_max:
        raise CaseInputError(
            "flow.shock_tube.diaphragm_position must lie strictly inside "
            "the x domain"
        )
    dx = (x_max - x_min) / nx
    face_index = round((position - x_min) / dx)
    aligned_position = x_min + face_index * dx
    if not math.isclose(position, aligned_position, rel_tol=1.0e-10, abs_tol=1.0e-12):
        raise CaseInputError(
            "flow.shock_tube.diaphragm_position must lie on a target "
            "x-cell boundary"
        )

    driver = _nse_primitive_state(
        tube.get("driver", {}), "flow.shock_tube.driver"
    )
    driven = _nse_primitive_state(
        tube.get("driven", {}), "flow.shock_tube.driven"
    )
    if driver[2] <= driven[2]:
        raise CaseInputError(
            "flow.shock_tube.driver.pressure must be greater than "
            "flow.shock_tube.driven.pressure"
        )
    if not math.isclose(driver[1][0], 0.0, rel_tol=0.0, abs_tol=1.0e-14):
        raise CaseInputError(
            "flow.shock_tube.driver.velocity x component must be zero at "
            "the reflective closed end"
        )

    rho4, velocity4, pressure4 = driver
    rho1, velocity1, pressure1 = driven
    values = {
        "shock_tube_diaphragm_position": aligned_position,
        "shock_tube_driver_rho": rho4,
        "shock_tube_driver_u": velocity4[0],
        "shock_tube_driver_v": velocity4[1],
        "shock_tube_driver_w": velocity4[2],
        "shock_tube_driver_p": pressure4,
        "shock_tube_driven_rho": rho1,
        "shock_tube_driven_u": velocity1[0],
        "shock_tube_driven_v": velocity1[1],
        "shock_tube_driven_w": velocity1[2],
        "shock_tube_driven_p": pressure1,
    }
    states = {
        "shock_tube.driver": driver,
        "shock_tube.driven": driven,
    }
    return values, states


def _resolve_nse_imported_turbulence(
    case: dict[str, Any],
    *,
    case_dir: Path | None = None,
    runtime_root: Path | None = None,
    shock_states: dict[
        str, tuple[float, tuple[float, float, float], float]
    ] | None = None,
) -> dict[str, Any]:
    flow_type = _canonical_selector(str(nested(case, "flow.type", "")))
    aliases = {
        "imported_turbulence",
        "turbulence_import",
        "turbulence_embed",
        "turbulence_tile",
        "shock_turbulence_interaction",
        "planar_shock_turbulence",
        "shock_turbulence",
        "shock_tube_turbulence_interaction",
        "shock_tube_turbulence",
        "finite_driver_shock_turbulence",
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

    planar_shock_interaction = flow_type in {
        "shock_turbulence_interaction",
        "planar_shock_turbulence",
        "shock_turbulence",
    }
    shock_tube_interaction = flow_type in {
        "shock_tube_turbulence_interaction",
        "shock_tube_turbulence",
        "finite_driver_shock_turbulence",
    }
    shock_interaction = planar_shock_interaction or shock_tube_interaction
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
    if shock_interaction and mode != "embed":
        raise CaseInputError(
            "shock-turbulence interaction requires imported_turbulence.mode=EMBED"
        )

    background_reference = None
    background_label = ""
    if planar_shock_interaction:
        background_reference = (shock_states or {}).get("planar_shock.upstream")
        background_label = "flow.planar_shock.upstream"
        if background_reference is None:
            raise CaseInputError(
                "shock-turbulence interaction requires flow.planar_shock"
            )
    elif shock_tube_interaction:
        background_reference = (shock_states or {}).get("shock_tube.driven")
        background_label = "flow.shock_tube.driven"
        if background_reference is None:
            raise CaseInputError(
                "shock-tube turbulence interaction requires flow.shock_tube"
            )

    velocity_offset = _vector3(
        imported.get(
            "velocity_offset",
            background_reference[1]
            if background_reference is not None
            else [0.0, 0.0, 0.0],
        ),
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
        background.get(
            "density",
            background_reference[0] if background_reference is not None else rho0,
        ),
        "flow.imported_turbulence.background.density",
    )
    background_pressure = _positive_float(
        background.get(
            "pressure",
            background_reference[2]
            if background_reference is not None
            else 1.0 / gamma,
        ),
        "flow.imported_turbulence.background.pressure",
    )
    background_velocity = _vector3(
        background.get(
            "velocity",
            background_reference[1]
            if background_reference is not None
            else [0.0, 0.0, 0.0],
        ),
        "flow.imported_turbulence.background.velocity",
    )
    x_start = _finite_float(
        imported.get("x_start", nested(case, "grid.x_min", 0.0)),
        "flow.imported_turbulence.x_start",
    )
    if background_reference is not None:
        actual = (background_rho, background_velocity, background_pressure)
        flattened_actual = (actual[0], *actual[1], actual[2])
        flattened_reference = (
            background_reference[0],
            *background_reference[1],
            background_reference[2],
        )
        if any(
            not math.isclose(a, b, rel_tol=1.0e-12, abs_tol=1.0e-14)
            for a, b in zip(flattened_actual, flattened_reference)
        ):
            raise CaseInputError(
                "flow.imported_turbulence.background must equal "
                f"{background_label}"
            )

    source_file = source_file.strip()
    if case_dir is not None:
        # BuildSolver launches the executable with case_dir as its working
        # directory.  Keep files contained in the case portable by writing a
        # case-relative path; using runtime_root here would duplicate
        # ``cases/<case_id>`` when the Fortran runtime opens the file.
        resolved_case_dir = case_dir.resolve()
        source_path = Path(source_file)
        if not source_path.is_absolute():
            source_path = (resolved_case_dir / source_path).resolve()
        else:
            source_path = source_path.resolve()
        try:
            source_file = source_path.relative_to(resolved_case_dir).as_posix()
        except ValueError:
            # An explicitly external data file cannot be represented by a
            # portable path inside the case, so retain its absolute location.
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


def _resolve_nse_boundary(
    case: dict[str, Any],
    *,
    shock_states: dict[
        str, tuple[float, tuple[float, float, float], float]
    ] | None = None,
) -> dict[str, Any]:
    numerics = _mapping(nested(case, "numerics", {}), "numerics")
    physics_nse = _mapping(nested(case, "physics.nse", {}), "physics.nse")
    legacy_numerics = numerics.get("boundary_condition")
    legacy_physics = physics_nse.get("boundary_condition")
    boundary = case.get("boundary")

    if boundary is None:
        legacy = legacy_numerics
        if legacy in {None, ""}:
            legacy = legacy_physics
        if legacy in {None, ""}:
            legacy = "periodic"
        normalized = _canonical_selector(str(legacy))
        if normalized != "periodic":
            raise CaseInputError(
                "legacy numerics.boundary_condition supports only PERIODIC; "
                "use the top-level boundary section for face-specific boundaries"
            )
        return {"boundary_condition": "periodic"}

    if legacy_numerics not in {None, ""} or legacy_physics not in {None, ""}:
        raise CaseInputError(
            "boundary and legacy numerics.boundary_condition/"
            "physics.nse.boundary_condition cannot be specified together"
        )

    boundary_map = _mapping(boundary, "boundary")
    allowed_boundary = {"faces", "reference_states", "non_reflecting"}
    unknown_boundary = sorted(set(boundary_map) - allowed_boundary)
    if unknown_boundary:
        raise CaseInputError(
            "unknown boundary key(s): " + ", ".join(unknown_boundary)
        )

    faces = _mapping(boundary_map.get("faces", {}), "boundary.faces")
    missing_faces = [face for face in NSE_BOUNDARY_FACES if face not in faces]
    unknown_faces = sorted(set(faces) - set(NSE_BOUNDARY_FACES))
    if missing_faces:
        raise CaseInputError(
            "boundary.faces must explicitly define all six physical faces; missing: "
            + ", ".join(missing_faces)
        )
    if unknown_faces:
        raise CaseInputError(
            "unknown boundary.faces key(s): " + ", ".join(unknown_faces)
        )

    reference_states = _mapping(
        boundary_map.get("reference_states", {}),
        "boundary.reference_states",
    )
    non_reflecting = _mapping(
        boundary_map.get("non_reflecting", {}),
        "boundary.non_reflecting",
    )
    allowed_non_reflecting = {
        "formulation",
        "relaxation_strength",
        "length_scale",
    }
    unknown_non_reflecting = sorted(
        set(non_reflecting) - allowed_non_reflecting
    )
    if unknown_non_reflecting:
        raise CaseInputError(
            "unknown boundary.non_reflecting key(s): "
            + ", ".join(unknown_non_reflecting)
        )
    formulation = _canonical_selector(
        str(non_reflecting.get("formulation", "characteristic_relaxation"))
    )
    if formulation != "characteristic_relaxation":
        raise CaseInputError(
            "boundary.non_reflecting.formulation must be "
            "CHARACTERISTIC_RELAXATION"
        )
    relaxation = _positive_float(
        non_reflecting.get("relaxation_strength", 0.1),
        "boundary.non_reflecting.relaxation_strength",
        allow_zero=True,
    )
    length_value = non_reflecting.get("length_scale", "auto")
    if isinstance(length_value, str):
        if _canonical_selector(length_value) != "auto":
            raise CaseInputError(
                "boundary.non_reflecting.length_scale must be AUTO or positive"
            )
        length_scale = -1.0
    else:
        length_scale = _positive_float(
            length_value, "boundary.non_reflecting.length_scale"
        )

    resolved_reference_states: dict[
        str, tuple[float, tuple[float, float, float], float]
    ] = {}
    for reference_name, raw_reference in reference_states.items():
        if not isinstance(reference_name, str) or not reference_name.strip():
            raise CaseInputError(
                "boundary.reference_states names must be non-empty strings"
            )
        reference_name = reference_name.strip()
        reference = _mapping(
            raw_reference,
            f"boundary.reference_states.{reference_name}",
        )
        if set(reference) == {"source"}:
            source = str(reference["source"]).strip().lower()
            if source not in (shock_states or {}):
                raise CaseInputError(
                    f"boundary.reference_states.{reference_name}.source="
                    f"{source!r} is not available"
                )
            resolved_reference_states[reference_name] = (shock_states or {})[
                source
            ]
            continue
        allowed_reference = {"density", "velocity", "pressure"}
        unknown_reference = sorted(set(reference) - allowed_reference)
        if unknown_reference:
            raise CaseInputError(
                f"unknown boundary.reference_states.{reference_name} key(s): "
                + ", ".join(unknown_reference)
            )
        missing_reference = sorted(allowed_reference - set(reference))
        if missing_reference:
            raise CaseInputError(
                f"boundary.reference_states.{reference_name} is missing: "
                + ", ".join(missing_reference)
            )
        density = _positive_float(
            reference["density"],
            f"boundary.reference_states.{reference_name}.density",
        )
        velocity = _vector3(
            reference["velocity"],
            f"boundary.reference_states.{reference_name}.velocity",
        )
        pressure = _positive_float(
            reference["pressure"],
            f"boundary.reference_states.{reference_name}.pressure",
        )
        resolved_reference_states[reference_name] = (
            density,
            velocity,
            pressure,
        )

    values: dict[str, Any] = {
        "boundary_relaxation_strength": relaxation,
        "boundary_length_scale": length_scale,
    }
    face_types: dict[str, str] = {}
    for face in NSE_BOUNDARY_FACES:
        face_config = _mapping(faces[face], f"boundary.faces.{face}")
        allowed_face = {"type", "reference_state"}
        unknown_face = sorted(set(face_config) - allowed_face)
        if unknown_face:
            raise CaseInputError(
                f"unknown boundary.faces.{face} key(s): "
                + ", ".join(unknown_face)
            )
        if "type" not in face_config:
            raise CaseInputError(f"boundary.faces.{face}.type is required")
        face_type = _canonical_selector(str(face_config["type"]))
        if face_type not in {
            "periodic",
            "non_reflecting",
            "reflective",
            "dirichlet",
        }:
            raise CaseInputError(
                f"boundary.faces.{face}.type must be PERIODIC, "
                "NON_REFLECTING, REFLECTIVE, or DIRICHLET"
            )
        face_types[face] = face_type
        values[f"boundary_{face}"] = face_type

        reference_name = face_config.get("reference_state")
        if face_type not in {"non_reflecting", "dirichlet"}:
            if reference_name not in {None, ""}:
                raise CaseInputError(
                    f"boundary.faces.{face}.reference_state is only valid for "
                    "NON_REFLECTING or DIRICHLET"
                )
            continue
        if not isinstance(reference_name, str) or not reference_name.strip():
            raise CaseInputError(
                f"boundary.faces.{face}.reference_state is required for "
                "NON_REFLECTING or DIRICHLET"
            )
        reference_name = reference_name.strip()
        if reference_name not in resolved_reference_states:
            raise CaseInputError(
                f"boundary.faces.{face}.reference_state={reference_name!r} "
                "is not defined in boundary.reference_states"
            )
        density, velocity, pressure = resolved_reference_states[reference_name]
        values[f"boundary_{face}_reference_rho"] = density
        values[f"boundary_{face}_reference_u"] = velocity[0]
        values[f"boundary_{face}_reference_v"] = velocity[1]
        values[f"boundary_{face}_reference_w"] = velocity[2]
        values[f"boundary_{face}_reference_p"] = pressure

    for lower, upper, direction in (
        ("x_min", "x_max", "x"),
        ("y_min", "y_max", "y"),
        ("z_min", "z_max", "z"),
    ):
        if (face_types[lower] == "periodic") != (
            face_types[upper] == "periodic"
        ):
            raise CaseInputError(
                f"periodic {direction} boundaries must be specified on both "
                f"{lower} and {upper}"
            )

    return values


def _validate_nse_shock_driver(
    shock_values: dict[str, Any], boundary_values: dict[str, Any]
) -> None:
    """Require the inflow reservoir to reproduce the post-shock state exactly."""
    if not shock_values:
        return

    direction = shock_values["planar_shock_direction"]
    driver_face = "x_min" if direction == "positive_x" else "x_max"
    if boundary_values.get(f"boundary_{driver_face}") != "dirichlet":
        raise CaseInputError(
            "shock-turbulence interaction requires boundary.faces."
            f"{driver_face}.type=DIRICHLET on the post-shock driver face"
        )

    component_pairs = (
        ("rho", "rho"),
        ("u", "u"),
        ("v", "v"),
        ("w", "w"),
        ("p", "p"),
    )
    for boundary_component, shock_component in component_pairs:
        actual = boundary_values.get(
            f"boundary_{driver_face}_reference_{boundary_component}"
        )
        expected = shock_values[f"planar_shock_downstream_{shock_component}"]
        if actual is None or not math.isclose(
            float(actual), float(expected), rel_tol=1.0e-12, abs_tol=1.0e-14
        ):
            raise CaseInputError(
                "shock-turbulence interaction requires the Dirichlet driver "
                "reference state to equal flow.planar_shock.downstream; use "
                "boundary.reference_states.<name>.source: "
                "planar_shock.downstream"
            )


def _validate_nse_shock_tube_configuration(
    tube_values: dict[str, Any],
    imported_values: dict[str, Any],
    boundary_values: dict[str, Any],
) -> None:
    if not tube_values:
        return

    diaphragm = float(tube_values["shock_tube_diaphragm_position"])
    turbulence_start = float(imported_values["imported_turbulence_x_start"])
    if diaphragm > turbulence_start and not math.isclose(
        diaphragm, turbulence_start, rel_tol=1.0e-12, abs_tol=1.0e-14
    ):
        raise CaseInputError(
            "flow.shock_tube.diaphragm_position must not lie inside or "
            "downstream of the imported turbulence block"
        )

    if boundary_values.get("boundary_x_min") != "reflective":
        raise CaseInputError(
            "shock-tube turbulence interaction requires "
            "boundary.faces.x_min.type=REFLECTIVE for the closed driver end"
        )
    if boundary_values.get("boundary_x_max") != "non_reflecting":
        raise CaseInputError(
            "shock-tube turbulence interaction requires "
            "boundary.faces.x_max.type=NON_REFLECTING"
        )

    for component in ("rho", "u", "v", "w", "p"):
        actual = boundary_values.get(f"boundary_x_max_reference_{component}")
        expected = tube_values[f"shock_tube_driven_{component}"]
        if actual is None or not math.isclose(
            float(actual), float(expected), rel_tol=1.0e-12, abs_tol=1.0e-14
        ):
            raise CaseInputError(
                "shock-tube turbulence interaction requires the x_max "
                "non-reflecting reference state to equal flow.shock_tube.driven; "
                "use boundary.reference_states.<name>.source: shock_tube.driven"
            )


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
            "planar_shock_turbulence": "shock_turbulence_interaction",
            "shock_turbulence": "shock_turbulence_interaction",
            "shock_tube_turbulence": "shock_tube_turbulence_interaction",
            "finite_driver_shock_turbulence": (
                "shock_tube_turbulence_interaction"
            ),
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
    shock_values, shock_states = _resolve_nse_planar_shock(case)
    nse.update(shock_values)
    tube_values, tube_states = _resolve_nse_shock_tube(case)
    nse.update(tube_values)
    flow_states = dict(shock_states)
    flow_states.update(tube_states)
    imported_turbulence_values = _resolve_nse_imported_turbulence(
        case,
        case_dir=case_dir,
        runtime_root=runtime_root,
        shock_states=flow_states,
    )
    nse.update(imported_turbulence_values)
    boundary_values = _resolve_nse_boundary(case, shock_states=flow_states)
    _validate_nse_shock_driver(shock_values, boundary_values)
    _validate_nse_shock_tube_configuration(
        tube_values, imported_turbulence_values, boundary_values
    )
    if case.get("boundary") is not None:
        nse.pop("boundary_condition", None)
    nse.update(boundary_values)
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
        "planar_shock_direction",
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
        face_types = [
            str(nse.get(f"boundary_{face}", nse.get("boundary_condition", "periodic")))
            for face in NSE_BOUNDARY_FACES
        ]
        if any(face_type != "periodic" for face_type in face_types):
            raise CaseInputError(
                "Petersen-Livescu forcing requires periodic boundaries on all "
                "six physical faces"
            )
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


def _multicomponent_grid_settings(
    case: dict[str, Any], stage_name: str
) -> tuple[list[int], list[float]]:
    grid = _mapping(nested(case, "grid", {}), "grid")
    dimensions: list[int] = []
    for name in ("nx", "ny", "nz"):
        value = grid.get(name)
        if isinstance(value, bool) or not isinstance(value, int) or value < 2:
            raise CaseInputError(f"grid.{name} must be an integer of at least 2")
        dimensions.append(value)
    extents = [
        _finite_float(_required(case, f"grid.{name}"), f"grid.{name}")
        for name in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
    ]
    if any(extents[index+1] <= extents[index] for index in (0, 2, 4)):
        raise CaseInputError(f"{stage_name} grid extents must be positive")
    return dimensions, extents


def _mass_fraction_vector(value: Any, label: str, nspecies: int) -> list[float]:
    if not isinstance(value, (list, tuple)) or len(value) != nspecies:
        raise CaseInputError(f"{label} must contain exactly {nspecies} values")
    fractions = [
        _positive_float(item, f"{label}[{index}]", allow_zero=True)
        for index, item in enumerate(value)
    ]
    if abs(sum(fractions)-1.0) > 1.0e-12:
        raise CaseInputError(f"{label} must sum to one")
    return fractions


def _nasa7_coefficients(value: Any, label: str) -> list[float]:
    if not isinstance(value, (list, tuple)) or len(value) != 7:
        raise CaseInputError(f"{label} must contain exactly seven coefficients")
    return [
        _finite_float(coefficient, f"{label}[{index}]")
        for index, coefficient in enumerate(value)
    ]


def _thermally_perfect_settings(
    thermodynamics: dict[str, Any], species: list[str]
) -> dict[str, Any]:
    universal_gas_constant = _positive_float(
        thermodynamics.get("universal_gas_constant", 8314.46261815324),
        "thermodynamics.universal_gas_constant",
    )
    temperature_min = _positive_float(
        thermodynamics.get("temperature_min", 200.0),
        "thermodynamics.temperature_min",
    )
    temperature_max = _positive_float(
        thermodynamics.get("temperature_max", 6000.0),
        "thermodynamics.temperature_max",
    )
    if temperature_max <= temperature_min:
        raise CaseInputError(
            "thermodynamics.temperature_max must exceed temperature_min"
        )
    temperature_tolerance = _positive_float(
        thermodynamics.get("temperature_tolerance", 1.0e-10),
        "thermodynamics.temperature_tolerance",
    )
    if temperature_tolerance >= 1.0e-3:
        raise CaseInputError(
            "thermodynamics.temperature_tolerance must be less than 1e-3"
        )

    species_data = _mapping(
        thermodynamics.get("species_data"),
        "thermodynamics.species_data",
    )
    missing = [name for name in species if name not in species_data]
    extra = [str(name) for name in species_data if name not in species]
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing={missing}")
        if extra:
            details.append(f"extra={extra}")
        raise CaseInputError(
            "thermodynamics.species_data keys must exactly match "
            f"physics.multicomponent.species ({', '.join(details)})"
        )

    molecular_weights: list[float] = []
    temperature_midpoints: list[float] = []
    nasa_low: list[float] = []
    nasa_high: list[float] = []
    for name in species:
        label = f"thermodynamics.species_data.{name}"
        properties = _mapping(species_data[name], label)
        molecular_weight = _positive_float(
            properties.get("molecular_weight"), f"{label}.molecular_weight"
        )
        midpoint = _positive_float(
            properties.get("temperature_midpoint", 1000.0),
            f"{label}.temperature_midpoint",
        )
        if not temperature_min < midpoint < temperature_max:
            raise CaseInputError(
                f"{label}.temperature_midpoint must lie inside the "
                "thermodynamics temperature range"
            )
        low = _nasa7_coefficients(properties.get("nasa7_low"), f"{label}.nasa7_low")
        high = _nasa7_coefficients(
            properties.get("nasa7_high"), f"{label}.nasa7_high"
        )

        species_gas_constant = universal_gas_constant / molecular_weight
        for range_name, coefficients, samples in (
            (
                "nasa7_low",
                low,
                (temperature_min, 0.5 * (temperature_min + midpoint), midpoint),
            ),
            (
                "nasa7_high",
                high,
                (midpoint, 0.5 * (midpoint + temperature_max), temperature_max),
            ),
        ):
            for temperature in samples:
                cp_over_r = sum(
                    coefficients[index] * temperature**index
                    for index in range(5)
                )
                if not math.isfinite(cp_over_r) or cp_over_r <= 1.0:
                    raise CaseInputError(
                        f"{label}.{range_name} produces non-positive cv "
                        f"at T={temperature:g}"
                    )
                cp_value = species_gas_constant * cp_over_r
                if not math.isfinite(cp_value):
                    raise CaseInputError(
                        f"{label}.{range_name} produces non-finite cp"
                    )

        molecular_weights.append(molecular_weight)
        temperature_midpoints.append(midpoint)
        nasa_low.extend(low)
        nasa_high.extend(high)

    return {
        "thermo_species_names": species,
        "universal_gas_constant": universal_gas_constant,
        "temperature_min": temperature_min,
        "temperature_max": temperature_max,
        "temperature_tolerance": temperature_tolerance,
        "molecular_weights": molecular_weights,
        "temperature_midpoints": temperature_midpoints,
        "nasa_low_coefficients": nasa_low,
        "nasa_high_coefficients": nasa_high,
    }


def render_nse_multicomponent(
    case: dict[str, Any], profile_name: str | None = None
) -> str:
    physics = _mapping(
        nested(case, "physics.multicomponent", {}),
        "physics.multicomponent",
    )
    raw_species = physics.get("species")
    if not isinstance(raw_species, list) or not raw_species:
        raise CaseInputError(
            "physics.multicomponent.species must be a non-empty YAML list"
        )
    species = [str(value).strip() for value in raw_species]
    if any(not value for value in species):
        raise CaseInputError("multicomponent species names must not be empty")
    if any(len(value) > 32 for value in species):
        raise CaseInputError("multicomponent species names must not exceed 32 characters")
    if len(species) > 64:
        raise CaseInputError("multicomponent NSE supports at most 64 species")
    if len(set(species)) != len(species):
        raise CaseInputError("multicomponent species names must be unique")

    simulation_mode = _canonical_selector(str(physics.get("mode", "foundation")))
    if simulation_mode not in {
        "foundation",
        "passive_scalar",
        "inviscid_euler",
        "thermally_perfect_euler",
    }:
        raise CaseInputError(
            "physics.multicomponent.mode must be foundation, passive_scalar, "
            "inviscid_euler, or thermally_perfect_euler"
        )
    if simulation_mode == "passive_scalar" and len(species) < 2:
        raise CaseInputError(
            "passive_scalar mode requires at least tracer and carrier species"
        )
    expected_modes = {
        "cpu_serial_foundation": "foundation",
        "cpu_serial_passive_scalar": "passive_scalar",
        "cpu_serial_inviscid": "inviscid_euler",
        "cpu_serial_thermally_perfect": "thermally_perfect_euler",
    }
    if profile_name in expected_modes and simulation_mode != expected_modes[profile_name]:
        raise CaseInputError(
            f"solver profile {profile_name!r} requires "
            f"physics.multicomponent.mode={expected_modes[profile_name]!r}"
        )

    thermodynamics = _mapping(
        nested(case, "thermodynamics", {}), "thermodynamics"
    )
    transport = _mapping(nested(case, "transport", {}), "transport")
    chemistry = _mapping(nested(case, "chemistry", {}), "chemistry")
    models = {
        "thermodynamics": _canonical_selector(
            str(thermodynamics.get("model", "calorically_perfect"))
        ),
        "transport": _canonical_selector(str(transport.get("model", "none"))),
        "chemistry": _canonical_selector(str(chemistry.get("model", "none"))),
    }
    expected_thermodynamics = (
        "thermally_perfect"
        if simulation_mode == "thermally_perfect_euler"
        else "calorically_perfect"
    )
    supported = {
        "thermodynamics": expected_thermodynamics,
        "transport": "none",
        "chemistry": "none",
    }
    for category, expected in supported.items():
        if models[category] != expected:
            raise CaseInputError(
                f"current multicomponent stages require {category}.model="
                f"{expected!r}, got {models[category]!r}"
            )

    stage_names = {
        "foundation": "Stage-0 foundation",
        "passive_scalar": "Stage-1 passive-scalar advection",
        "inviscid_euler": "Stage-2 non-reacting inviscid multicomponent Euler",
        "thermally_perfect_euler": (
            "Stage-3 thermally-perfect non-reacting multicomponent Euler"
        ),
    }
    lines = [
        "! Automatically generated from case.yaml.",
        f"! {stage_names[simulation_mode]}.",
        "",
        "&multicomponent",
    ]
    _append(
        lines,
        [
            ("nspecies", len(species)),
            ("species_names", species),
            ("simulation_mode", simulation_mode),
            ("thermodynamics_model", models["thermodynamics"]),
            ("transport_model", models["transport"]),
            ("chemistry_model", models["chemistry"]),
        ],
    )
    lines.extend(["/", ""])

    thermally_perfect: dict[str, Any] | None = None
    if simulation_mode == "thermally_perfect_euler":
        thermally_perfect = _thermally_perfect_settings(
            thermodynamics, species
        )
        lines.append("&thermally_perfect")
        _append(lines, list(thermally_perfect.items()))
        lines.extend(["/", ""])

    if simulation_mode == "passive_scalar":
        flow = _mapping(nested(case, "flow", {}), "flow")
        flow_type = _canonical_selector(str(flow.get("type", "")))
        if flow_type != "passive_scalar_advection":
            raise CaseInputError(
                "stage-1 requires flow.type='passive_scalar_advection'"
            )
        initial = _mapping(
            nested(case, "flow.passive_scalar", {}),
            "flow.passive_scalar",
        )
        time = _mapping(nested(case, "time", {}), "time")
        numerics = _mapping(nested(case, "numerics", {}), "numerics")
        output = _mapping(nested(case, "output", {}), "output")

        dimensions, extents = _multicomponent_grid_settings(
            case, "passive-scalar"
        )
        velocity = _vector3(_required(case, "flow.velocity"), "flow.velocity")
        center = _vector3(
            initial.get("tracer_center", [0.25, 0.5, 0.5]),
            "flow.passive_scalar.tracer_center",
        )
        if any(
            center[axis] < extents[2*axis]
            or center[axis] > extents[2*axis+1]
            for axis in range(3)
        ):
            raise CaseInputError(
                "flow.passive_scalar.tracer_center must lie inside the domain"
            )
        tracer_background = _positive_float(
            initial.get("tracer_background", 0.05),
            "flow.passive_scalar.tracer_background",
            allow_zero=True,
        )
        tracer_amplitude = _positive_float(
            initial.get("tracer_amplitude", 0.90),
            "flow.passive_scalar.tracer_amplitude",
            allow_zero=True,
        )
        if tracer_background + tracer_amplitude > 1.0:
            raise CaseInputError(
                "passive-scalar tracer_background + tracer_amplitude must not exceed 1"
            )
        tracer_width = _positive_float(
            initial.get("tracer_width", 0.08),
            "flow.passive_scalar.tracer_width",
        )
        nsteps = time.get("nsteps")
        if isinstance(nsteps, bool) or not isinstance(nsteps, int) or nsteps < 0:
            raise CaseInputError("time.nsteps must be a non-negative integer")
        cfl = _positive_float(time.get("cfl", 0.45), "time.cfl")
        if cfl > 1.0:
            raise CaseInputError("stage-1 passive-scalar CFL must not exceed 1")
        fixed_dt = _positive_float(time.get("dt", 0.0), "time.dt", allow_zero=True)
        cell_widths = [
            (extents[2*axis+1]-extents[2*axis]) / dimensions[axis]
            for axis in range(3)
        ]
        advection_rate = sum(
            abs(velocity[axis]) / cell_widths[axis] for axis in range(3)
        )
        if fixed_dt == 0.0 and advection_rate == 0.0:
            raise CaseInputError(
                "zero passive-scalar velocity requires a positive time.dt"
            )
        if fixed_dt * advection_rate > 1.0 + 1.0e-14:
            raise CaseInputError("time.dt violates the stage-1 upwind CFL limit")

        initial_condition = _canonical_selector(
            str(initial.get("initial_condition", "gaussian"))
        )
        advection_scheme = _canonical_selector(
            str(numerics.get("convective_scheme", "upwind1"))
        )
        boundary_condition = _canonical_selector(
            str(numerics.get("boundary_condition", "periodic"))
        )
        time_integrator = _canonical_selector(
            str(numerics.get("time_integration", "ssprk3"))
        )
        supported_values = {
            "flow.passive_scalar.initial_condition": (initial_condition, "gaussian"),
            "numerics.convective_scheme": (advection_scheme, "upwind1"),
            "numerics.boundary_condition": (boundary_condition, "periodic"),
            "numerics.time_integration": (time_integrator, "ssprk3"),
        }
        for label, (actual, expected) in supported_values.items():
            if actual != expected:
                raise CaseInputError(
                    f"stage-1 requires {label}={expected!r}, got {actual!r}"
                )
        write_final = output.get("write_final", True)
        if not isinstance(write_final, bool):
            raise CaseInputError("output.write_final must be true or false")
        output_file = str(output.get("filename", "passive_scalar_final.csv")).strip()
        if write_final and not output_file:
            raise CaseInputError("output.filename must not be empty")

        lines.append("&passive_scalar")
        _append(
            lines,
            [
                ("nx", dimensions[0]),
                ("ny", dimensions[1]),
                ("nz", dimensions[2]),
                ("x_min", extents[0]),
                ("x_max", extents[1]),
                ("y_min", extents[2]),
                ("y_max", extents[3]),
                ("z_min", extents[4]),
                ("z_max", extents[5]),
                ("velocity", velocity),
                ("cfl", cfl),
                ("dt", fixed_dt),
                ("nsteps", nsteps),
                ("initial_condition", initial_condition),
                ("tracer_background", tracer_background),
                ("tracer_amplitude", tracer_amplitude),
                ("tracer_center", center),
                ("tracer_width", tracer_width),
                ("advection_scheme", advection_scheme),
                ("boundary_condition", boundary_condition),
                ("time_integrator", time_integrator),
                ("write_final", write_final),
                ("output_file", output_file),
            ],
        )
        lines.extend(["/", ""])
    elif simulation_mode in {"inviscid_euler", "thermally_perfect_euler"}:
        stage_number = 3 if simulation_mode == "thermally_perfect_euler" else 2
        flow = _mapping(nested(case, "flow", {}), "flow")
        flow_type = _canonical_selector(str(flow.get("type", "")))
        if flow_type != "multispecies_sod":
            raise CaseInputError(
                f"stage-{stage_number} requires flow.type='multispecies_sod'"
            )
        initial = _mapping(
            nested(case, "flow.multispecies_sod", {}),
            "flow.multispecies_sod",
        )
        left = _mapping(initial.get("left"), "flow.multispecies_sod.left")
        right = _mapping(initial.get("right"), "flow.multispecies_sod.right")
        time = _mapping(nested(case, "time", {}), "time")
        numerics = _mapping(nested(case, "numerics", {}), "numerics")
        output = _mapping(nested(case, "output", {}), "output")
        dimensions, extents = _multicomponent_grid_settings(
            case, "multicomponent Euler"
        )

        gamma = 1.4
        if simulation_mode == "inviscid_euler":
            gamma = _positive_float(
                thermodynamics.get("gamma", 1.4), "thermodynamics.gamma"
            )
            if gamma <= 1.0:
                raise CaseInputError("thermodynamics.gamma must exceed one")
        interface_location = _finite_float(
            initial.get("interface_location", 0.5),
            "flow.multispecies_sod.interface_location",
        )
        if not extents[0] < interface_location < extents[1]:
            raise CaseInputError(
                "flow.multispecies_sod.interface_location must lie inside x"
            )
        left_density = _positive_float(
            left.get("density"), "flow.multispecies_sod.left.density"
        )
        left_velocity = _vector3(
            left.get("velocity"), "flow.multispecies_sod.left.velocity"
        )
        left_pressure = _positive_float(
            left.get("pressure"), "flow.multispecies_sod.left.pressure"
        )
        left_fractions = _mass_fraction_vector(
            left.get("mass_fractions"),
            "flow.multispecies_sod.left.mass_fractions",
            len(species),
        )
        right_density = _positive_float(
            right.get("density"), "flow.multispecies_sod.right.density"
        )
        right_velocity = _vector3(
            right.get("velocity"), "flow.multispecies_sod.right.velocity"
        )
        right_pressure = _positive_float(
            right.get("pressure"), "flow.multispecies_sod.right.pressure"
        )
        right_fractions = _mass_fraction_vector(
            right.get("mass_fractions"),
            "flow.multispecies_sod.right.mass_fractions",
            len(species),
        )
        if thermally_perfect is not None:
            for side_name, density, pressure_value, fractions in (
                ("left", left_density, left_pressure, left_fractions),
                ("right", right_density, right_pressure, right_fractions),
            ):
                mixture_gas_constant = thermally_perfect[
                    "universal_gas_constant"
                ] * sum(
                    fractions[index]
                    / thermally_perfect["molecular_weights"][index]
                    for index in range(len(species))
                )
                initial_temperature = pressure_value / (
                    density * mixture_gas_constant
                )
                if not (
                    thermally_perfect["temperature_min"]
                    <= initial_temperature
                    <= thermally_perfect["temperature_max"]
                ):
                    raise CaseInputError(
                        f"flow.multispecies_sod.{side_name} gives "
                        f"T={initial_temperature:g}, outside the configured "
                        "thermodynamics temperature range"
                    )

        nsteps = time.get("nsteps")
        if isinstance(nsteps, bool) or not isinstance(nsteps, int) or nsteps < 0:
            raise CaseInputError("time.nsteps must be a non-negative integer")
        cfl = _positive_float(time.get("cfl",0.35), "time.cfl")
        if cfl > 1.0:
            raise CaseInputError(
                f"stage-{stage_number} multicomponent Euler CFL must not exceed 1"
            )
        fixed_dt = _positive_float(
            time.get("dt",0.0), "time.dt", allow_zero=True
        )
        initial_condition = _canonical_selector(
            str(initial.get("initial_condition","multispecies_sod_x"))
        )
        riemann_solver = _canonical_selector(
            str(numerics.get("convective_scheme","rusanov1"))
        )
        boundary_condition = _canonical_selector(
            str(numerics.get("boundary_condition","periodic"))
        )
        time_integrator = _canonical_selector(
            str(numerics.get("time_integration","ssprk3"))
        )
        supported_values = {
            "flow.multispecies_sod.initial_condition": (
                initial_condition, "multispecies_sod_x"
            ),
            "numerics.convective_scheme": (riemann_solver,"rusanov1"),
            "numerics.boundary_condition": (boundary_condition,"periodic"),
            "numerics.time_integration": (time_integrator,"ssprk3"),
        }
        for label, (actual, expected) in supported_values.items():
            if actual != expected:
                raise CaseInputError(
                    f"stage-{stage_number} requires {label}={expected!r}, "
                    f"got {actual!r}"
                )
        write_final = output.get("write_final",True)
        if not isinstance(write_final,bool):
            raise CaseInputError("output.write_final must be true or false")
        output_file = str(
            output.get("filename","multicomponent_euler_final.csv")
        ).strip()
        if write_final and not output_file:
            raise CaseInputError("output.filename must not be empty")

        lines.append("&multicomponent_euler")
        _append(
            lines,
            [
                ("nx",dimensions[0]),
                ("ny",dimensions[1]),
                ("nz",dimensions[2]),
                ("x_min",extents[0]),
                ("x_max",extents[1]),
                ("y_min",extents[2]),
                ("y_max",extents[3]),
                ("z_min",extents[4]),
                ("z_max",extents[5]),
                ("gamma",gamma if simulation_mode == "inviscid_euler" else None),
                ("cfl",cfl),
                ("dt",fixed_dt),
                ("nsteps",nsteps),
                ("initial_condition",initial_condition),
                ("interface_location",interface_location),
                ("left_density",left_density),
                ("left_velocity",left_velocity),
                ("left_pressure",left_pressure),
                ("left_mass_fractions",left_fractions),
                ("right_density",right_density),
                ("right_velocity",right_velocity),
                ("right_pressure",right_pressure),
                ("right_mass_fractions",right_fractions),
                ("riemann_solver",riemann_solver),
                ("boundary_condition",boundary_condition),
                ("time_integrator",time_integrator),
                ("write_final",write_final),
                ("output_file",output_file),
            ],
        )
        lines.extend(["/",""])
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
    if model.lower() == "nse_multicomponent":
        return input_name, render_nse_multicomponent(case, profile_name)
    raise CaseInputError(f"unsupported model: {model}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a solver input namelist from case.yaml."
    )
    parser.add_argument("--case", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument(
        "--model",
        required=True,
        help="Model identifier declared by the selected solver manifest",
    )
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
