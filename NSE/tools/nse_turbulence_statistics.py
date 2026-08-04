#!/usr/bin/env python3
"""Compute time-resolved turbulence statistics from NSE SLF output."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

from slf_to_paraview_merged_cropghost import (
    build_rank_range_lookup,
    complete_step_groups,
    crop_center_to_shape,
    discover_files,
    get_global_grid,
    get_origin_spacing,
    group_by_step,
    load_meta,
    parse_step_rank,
    rank_slices,
    read_slf,
    select_step_groups,
)


CSV_COLUMNS = (
    "step",
    "time",
    "mean_density",
    "mean_pressure",
    "mean_sound_speed",
    "mean_u",
    "mean_v",
    "mean_w",
    "rms_u_fluctuation",
    "rms_v_fluctuation",
    "rms_w_fluctuation",
    "rms_velocity_fluctuation",
    "turbulent_kinetic_energy",
    "reynolds_stress_uv",
    "reynolds_stress_uw",
    "reynolds_stress_vw",
    "isotropy_error",
    "kinematic_viscosity",
    "dissipation_rate",
    "integral_length_scale",
    "taylor_microscale",
    "kolmogorov_length_scale",
    "reynolds_integral",
    "reynolds_taylor",
    "turbulent_mach_number",
    "spectral_energy_relative_error",
)


def _nse_indices(names: list[str]) -> tuple[int, int, int, int, int]:
    lookup = {name.lower(): index for index, name in enumerate(names)}
    required = ("rho", "rho_u", "rho_v", "rho_w", "rho_e")
    missing = [name for name in required if name not in lookup]
    if missing:
        raise ValueError(
            f"NSE statistics require {required}; missing={missing}, found={names}"
        )
    return tuple(lookup[name] for name in required)  # type: ignore[return-value]


def _physical_local_data(
    slf,
    rank: int | None,
    rank_lookup: dict[int, dict],
    global_shape: tuple[int, int, int],
) -> tuple[np.ndarray, tuple[slice, slice, slice]]:
    if rank is not None and rank_lookup and rank not in rank_lookup:
        raise ValueError(f"{slf.path}: rank {rank} is missing from meta.json")

    if rank is None:
        slices = tuple(slice(0, value) for value in global_shape)
    else:
        slices = rank_slices(rank, slf, rank_lookup)

    expected = tuple(axis.stop - axis.start for axis in slices)
    local_data = slf.data
    if slf.shape[:3] != expected:
        print(
            f"INFO: {slf.path.name}: cropping ghost cells "
            f"{slf.shape[:3]} -> {expected}",
            file=sys.stderr,
        )
        local_data = crop_center_to_shape(local_data, expected, slf.path)
    return local_data, slices


def merge_velocity_state(
    files: list[Path],
    meta: dict,
    rank_lookup: dict[int, dict],
    gamma: float,
    density_floor: float,
    pressure_floor: float,
) -> tuple[
    int,
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, float],
    tuple[float, float, float],
]:
    """Merge one output step and retain only the primitive velocity arrays."""
    first = read_slf(files[0])
    step, _ = parse_step_rank(files[0], first)
    global_shape = get_global_grid(meta, first)
    u = np.empty(global_shape, dtype=np.float64)
    v = np.empty(global_shape, dtype=np.float64)
    w = np.empty(global_shape, dtype=np.float64)
    filled = np.zeros(global_shape, dtype=bool)

    density_sum = 0.0
    pressure_sum = 0.0
    sound_speed_sum = 0.0
    point_count = 0
    time = float(first.time)

    for file_index, path in enumerate(files):
        slf = first if file_index == 0 else read_slf(path)
        step_read, rank = parse_step_rank(path, slf)
        if step_read != step:
            raise ValueError(f"Internal error: step {step_read} is mixed with {step}")
        if not math.isclose(float(slf.time), time, rel_tol=1.0e-12, abs_tol=1.0e-14):
            raise ValueError(f"{path}: rank files have inconsistent time values")

        data, slices = _physical_local_data(slf, rank, rank_lookup, global_shape)
        rho_i, rhou_i, rhov_i, rhow_i, rhoe_i = _nse_indices(slf.names)
        rho = data[..., rho_i]
        if not np.all(np.isfinite(data)):
            raise ValueError(f"{path}: SLF contains NaN or infinite values")
        minimum_density = float(np.min(rho))
        if minimum_density <= density_floor:
            raise ValueError(
                f"{path}: density minimum {minimum_density:.8e} is not above "
                f"the floor {density_floor:.8e}"
            )

        local_u = data[..., rhou_i] / rho
        local_v = data[..., rhov_i] / rho
        local_w = data[..., rhow_i] / rho
        pressure = (gamma - 1.0) * (
            data[..., rhoe_i]
            - 0.5 * rho * (local_u * local_u + local_v * local_v + local_w * local_w)
        )
        minimum_pressure = float(np.min(pressure))
        if minimum_pressure <= pressure_floor:
            raise ValueError(
                f"{path}: pressure minimum {minimum_pressure:.8e} is not above "
                f"the floor {pressure_floor:.8e}"
            )
        sound_speed = np.sqrt(gamma * pressure / rho)

        if np.any(filled[slices]):
            raise ValueError(f"{path}: rank ranges overlap previously merged data")
        u[slices] = local_u
        v[slices] = local_v
        w[slices] = local_w
        filled[slices] = True

        density_sum += float(np.sum(rho, dtype=np.float64))
        pressure_sum += float(np.sum(pressure, dtype=np.float64))
        sound_speed_sum += float(np.sum(sound_speed, dtype=np.float64))
        point_count += int(rho.size)

    missing = int((~filled).sum())
    if missing:
        raise ValueError(f"Step {step} has {missing} unfilled cells after rank merge")
    expected_count = int(np.prod(global_shape))
    if point_count != expected_count:
        raise ValueError(
            f"Step {step} merged {point_count} points; expected {expected_count}"
        )

    thermodynamics = {
        "mean_density": density_sum / point_count,
        "mean_pressure": pressure_sum / point_count,
        "mean_sound_speed": sound_speed_sum / point_count,
    }
    _, spacing = get_origin_spacing(meta, first)
    return step, time, u, v, w, thermodynamics, spacing


def _isotropy_error(reynolds: np.ndarray) -> float:
    trace = float(np.trace(reynolds))
    if trace <= np.finfo(np.float64).tiny:
        return math.nan
    target = trace / 3.0
    diagonal_error = max(
        abs(float(reynolds[index, index]) / target - 1.0) for index in range(3)
    )
    off_diagonal_error = max(
        abs(float(reynolds[i, j]) / target)
        for i in range(3)
        for j in range(i + 1, 3)
    )
    return max(diagonal_error, off_diagonal_error)


def spectral_length_and_dissipation(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    spacing: tuple[float, float, float],
    kinematic_viscosity: float,
) -> tuple[float, float, float]:
    """Return integral scale, dissipation, and spectral energy consistency."""
    shape = u.shape
    point_count = int(np.prod(shape))
    transform_u = np.fft.fftn(u)
    transform_v = np.fft.fftn(v)
    transform_w = np.fft.fftn(w)

    kx = 2.0 * np.pi * np.fft.fftfreq(shape[0], d=spacing[0])
    ky = 2.0 * np.pi * np.fft.fftfreq(shape[1], d=spacing[1])
    kz = 2.0 * np.pi * np.fft.fftfreq(shape[2], d=spacing[2])
    kx_2d = kx[:, None]
    ky_2d = ky[None, :]

    energy_sum = 0.0
    integral_numerator = 0.0
    dissipation_sum = 0.0
    for iz, kz_value in enumerate(kz):
        u_hat = transform_u[:, :, iz]
        v_hat = transform_v[:, :, iz]
        w_hat = transform_w[:, :, iz]
        velocity_norm = (
            np.abs(u_hat) ** 2 + np.abs(v_hat) ** 2 + np.abs(w_hat) ** 2
        )
        k_squared = kx_2d * kx_2d + ky_2d * ky_2d + kz_value * kz_value
        retained = k_squared > 0.0
        mode_energy = 0.5 * velocity_norm
        energy_sum += float(np.sum(mode_energy[retained], dtype=np.float64))
        integral_numerator += float(
            np.sum(
                mode_energy[retained] / np.sqrt(k_squared[retained]),
                dtype=np.float64,
            )
        )

        wave_dot_velocity = (
            kx_2d * u_hat + ky_2d * v_hat + kz_value * w_hat
        )
        dissipation_sum += float(
            np.sum(
                k_squared * velocity_norm
                + (1.0 / 3.0) * np.abs(wave_dot_velocity) ** 2,
                dtype=np.float64,
            )
        )

    normalization = float(point_count) ** 2
    spectral_energy = energy_sum / normalization
    physical_energy = 0.5 * float(
        np.mean(u * u + v * v + w * w, dtype=np.float64)
    )
    if energy_sum <= np.finfo(np.float64).tiny:
        integral_scale = math.nan
    else:
        integral_scale = 0.75 * np.pi * integral_numerator / energy_sum
    dissipation = kinematic_viscosity * dissipation_sum / normalization
    relative_error = abs(spectral_energy - physical_energy) / max(
        physical_energy, np.finfo(np.float64).tiny
    )
    return integral_scale, dissipation, relative_error


def calculate_statistics(
    step: int,
    time: float,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    thermodynamics: dict[str, float],
    spacing: tuple[float, float, float],
    reynolds_number: float,
) -> dict[str, float | int]:
    if reynolds_number <= 0.0:
        raise ValueError("reynolds_number must be positive")
    if u.shape != v.shape or u.shape != w.shape:
        raise ValueError("velocity components must have the same shape")

    mean_u = float(np.mean(u, dtype=np.float64))
    mean_v = float(np.mean(v, dtype=np.float64))
    mean_w = float(np.mean(w, dtype=np.float64))
    u -= mean_u
    v -= mean_v
    w -= mean_w

    reynolds = np.array(
        [
            [np.mean(u * u), np.mean(u * v), np.mean(u * w)],
            [np.mean(u * v), np.mean(v * v), np.mean(v * w)],
            [np.mean(u * w), np.mean(v * w), np.mean(w * w)],
        ],
        dtype=np.float64,
    )
    rms_u = math.sqrt(max(float(reynolds[0, 0]), 0.0))
    rms_v = math.sqrt(max(float(reynolds[1, 1]), 0.0))
    rms_w = math.sqrt(max(float(reynolds[2, 2]), 0.0))
    variance_sum = float(np.trace(reynolds))
    velocity_rms = math.sqrt(max(variance_sum, 0.0))
    component_rms = math.sqrt(max(variance_sum / 3.0, 0.0))
    kinetic_energy = 0.5 * variance_sum

    mean_density = thermodynamics["mean_density"]
    kinematic_viscosity = 1.0 / (reynolds_number * mean_density)
    integral_scale, dissipation, energy_error = spectral_length_and_dissipation(
        u, v, w, spacing, kinematic_viscosity
    )
    if dissipation <= np.finfo(np.float64).tiny:
        taylor_scale = math.nan
        kolmogorov_scale = math.nan
    else:
        taylor_scale = math.sqrt(
            15.0 * kinematic_viscosity * component_rms**2 / dissipation
        )
        kolmogorov_scale = (kinematic_viscosity**3 / dissipation) ** 0.25

    reynolds_integral = (
        component_rms * integral_scale / kinematic_viscosity
        if math.isfinite(integral_scale)
        else math.nan
    )
    reynolds_taylor = (
        component_rms * taylor_scale / kinematic_viscosity
        if math.isfinite(taylor_scale)
        else math.nan
    )
    sound_speed = thermodynamics["mean_sound_speed"]
    turbulent_mach = velocity_rms / sound_speed

    return {
        "step": step,
        "time": time,
        **thermodynamics,
        "mean_u": mean_u,
        "mean_v": mean_v,
        "mean_w": mean_w,
        "rms_u_fluctuation": rms_u,
        "rms_v_fluctuation": rms_v,
        "rms_w_fluctuation": rms_w,
        "rms_velocity_fluctuation": velocity_rms,
        "turbulent_kinetic_energy": kinetic_energy,
        "reynolds_stress_uv": float(reynolds[0, 1]),
        "reynolds_stress_uw": float(reynolds[0, 2]),
        "reynolds_stress_vw": float(reynolds[1, 2]),
        "isotropy_error": _isotropy_error(reynolds),
        "kinematic_viscosity": kinematic_viscosity,
        "dissipation_rate": dissipation,
        "integral_length_scale": integral_scale,
        "taylor_microscale": taylor_scale,
        "kolmogorov_length_scale": kolmogorov_scale,
        "reynolds_integral": reynolds_integral,
        "reynolds_taylor": reynolds_taylor,
        "turbulent_mach_number": turbulent_mach,
        "spectral_energy_relative_error": energy_error,
    }


def _write_csv(path: Path, records: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    key: record[key]
                    if key == "step"
                    else format(float(record[key]), ".17g")
                    for key in CSV_COLUMNS
                }
            )


def _write_metadata(
    csv_path: Path,
    input_path: Path,
    meta_path: Path | None,
    gamma: float,
    reynolds_number: float,
    records: list[dict[str, float | int]],
) -> Path:
    metadata_path = csv_path.with_name(f"{csv_path.stem}_metadata.json")
    metadata = {
        "schema_version": 1,
        "source": {
            "slf": str(input_path.resolve()),
            "meta": str(meta_path.resolve()) if meta_path is not None else None,
        },
        "parameters": {
            "gamma": gamma,
            "reference_reynolds_number": reynolds_number,
        },
        "definitions": {
            "velocity_fluctuations": "Volume RMS about each component's volume mean.",
            "integral_length_scale": "(3*pi/4) * sum(E(k)/k) / sum(E(k)), excluding k=0.",
            "dissipation_rate": "nu * <2*S_ij*S_ij - (2/3)*(div u)^2>, evaluated spectrally.",
            "taylor_microscale": "sqrt(15*nu*u_prime^2/epsilon), with u_prime^2=<u_i'u_i'>/3.",
            "kolmogorov_length_scale": "(nu^3/epsilon)^(1/4).",
            "reynolds_integral": "u_prime*L/nu.",
            "reynolds_taylor": "u_prime*lambda/nu.",
            "turbulent_mach_number": "sqrt(<u_i'u_i'>)/<sqrt(gamma*p/rho)>.",
            "kinematic_viscosity": "1/(reference_reynolds_number*mean_density).",
            "isotropy_error": "Maximum normalized diagonal or off-diagonal Reynolds-stress error.",
        },
        "columns": list(CSV_COLUMNS),
        "steps": [int(record["step"]) for record in records],
    }
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=True)
        handle.write("\n")
    return metadata_path


def _default_meta_path(input_path: Path) -> Path | None:
    candidates = (
        (input_path / "meta.json", input_path.parent / "meta.json")
        if input_path.is_dir()
        else (input_path.parent / "meta.json", input_path.parent.parent / "meta.json")
    )
    return next((path for path in candidates if path.is_file()), None)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute an NSE turbulence-statistics time series from SLF output."
    )
    parser.add_argument("input", help="SLF file or output directory")
    parser.add_argument(
        "-o",
        "--output",
        default="turbulence_statistics.csv",
        help="Output CSV path",
    )
    parser.add_argument("--meta", help="meta.json path; inferred from input by default")
    parser.add_argument(
        "--steps",
        default="all",
        help="all, latest, comma-separated steps, or inclusive start:stop:stride",
    )
    parser.add_argument(
        "--layout", choices=("auto", "global", "rank"), default="auto"
    )
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--reynolds", type=float, required=True)
    parser.add_argument("--density-floor", type=float, default=1.0e-12)
    parser.add_argument("--pressure-floor", type=float, default=1.0e-12)
    args = parser.parse_args(argv)

    if args.gamma <= 1.0:
        parser.error("--gamma must be greater than one")
    if args.reynolds <= 0.0:
        parser.error("--reynolds must be positive")
    if args.density_floor < 0.0 or args.pressure_floor < 0.0:
        parser.error("density and pressure floors must be non-negative")

    input_path = Path(args.input)
    meta_path = Path(args.meta) if args.meta else _default_meta_path(input_path)
    try:
        meta = load_meta(meta_path) if meta_path is not None else {}
        equation = str(meta.get("equation", "nse")).strip().lower()
        if equation != "nse":
            raise ValueError(f"NSE statistics cannot process equation={equation!r}")
        case_name = str(meta.get("case_name", "")).strip() or None
        files = discover_files(
            input_path, case_name=case_name, layout=args.layout, meta=meta
        )
        if not files:
            raise ValueError(f"No SLF files found: {input_path}")
        groups = complete_step_groups(group_by_step(files), meta)
        groups = select_step_groups(groups, args.steps)
        if not groups:
            raise ValueError("No complete SLF steps are available")
        rank_lookup = build_rank_range_lookup(meta)

        records: list[dict[str, float | int]] = []
        for step, step_files in groups.items():
            (
                step_read,
                time,
                u,
                v,
                w,
                thermodynamics,
                spacing,
            ) = merge_velocity_state(
                step_files,
                meta,
                rank_lookup,
                args.gamma,
                args.density_floor,
                args.pressure_floor,
            )
            record = calculate_statistics(
                step_read,
                time,
                u,
                v,
                w,
                thermodynamics,
                spacing,
                args.reynolds,
            )
            records.append(record)
            print(
                f"[OK] step={step_read} time={time:.8e} "
                f"u_rms={record['rms_u_fluctuation']:.8e} "
                f"L={record['integral_length_scale']:.8e} "
                f"Re_lambda={record['reynolds_taylor']:.8e} "
                f"Mt={record['turbulent_mach_number']:.8e}"
            )

        output_path = Path(args.output)
        _write_csv(output_path, records)
        metadata_path = _write_metadata(
            output_path,
            input_path,
            meta_path,
            args.gamma,
            args.reynolds,
            records,
        )
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        parser.error(str(exc))

    print(f"Wrote turbulence time series: {output_path}")
    print(f"Wrote statistics definitions: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
