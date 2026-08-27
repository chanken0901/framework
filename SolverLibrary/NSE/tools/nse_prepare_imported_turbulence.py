#!/usr/bin/env python3
"""Create a ghost-free, MPI-layout-independent NSE turbulence SLF."""

from __future__ import annotations

import argparse
import json
import struct
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
    select_variables,
)


REQUIRED_FIELDS = ("rho", "rho_u", "rho_v", "rho_w", "rho_E")


def _default_meta_path(input_path: Path) -> Path | None:
    candidate = input_path / "meta.json" if input_path.is_dir() else input_path.parent / "meta.json"
    return candidate if candidate.exists() else None


def _select_step(groups: dict[int, list[Path]], requested: str) -> int:
    if not groups:
        raise ValueError("no complete SLF steps are available")
    if requested.lower() == "latest":
        return max(groups)
    try:
        step = int(requested)
    except ValueError as exc:
        raise ValueError("--step must be an integer or 'latest'") from exc
    if step not in groups:
        raise ValueError(f"requested SLF step {step} is not complete or does not exist")
    return step


def _write_header(
    stream,
    shape4: tuple[int, int, int, int],
    names: list[str],
    step: int,
    time: float,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
) -> int:
    nx, ny, nz, nvar = shape4
    if nvar != len(names):
        raise ValueError("SLF variable-name count does not match the field")
    shape = np.asarray([nx, ny, nz, nvar], dtype="<i4")
    metadata = np.asarray([step, 0, nx, ny, nz, 0, 1, 0], dtype="<i4")
    bounds = (
        origin[0], origin[0] + spacing[0] * nx,
        origin[1], origin[1] + spacing[1] * ny,
        origin[2], origin[2] + spacing[2] * nz,
    )

    stream.write(b"SLF1\x00\x00\x00\x00")
    stream.write(struct.pack("<i", 1))
    stream.write(struct.pack("<i", 2))
    stream.write(struct.pack("<i", 4))
    stream.write(shape.tobytes())
    stream.write(metadata.tobytes())
    stream.write(struct.pack("<d", float(time)))
    stream.write(struct.pack("<6d", *bounds))
    stream.write(struct.pack("<i", nvar))
    for name in names:
        encoded = name.encode("ascii")
        if len(encoded) > 32:
            raise ValueError(f"SLF variable name is longer than 32 bytes: {name!r}")
        stream.write(encoded.ljust(32, b" "))
    return stream.tell()


def write_global_slf(
    output_path: Path,
    data: np.ndarray,
    names: list[str],
    step: int,
    time: float,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
) -> None:
    if data.ndim != 4:
        raise ValueError(f"expected a four-dimensional field; got shape {data.shape}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as stream:
        _write_header(
            stream,
            tuple(map(int, data.shape)),
            names,
            step,
            time,
            origin,
            spacing,
        )
        field = np.asarray(data, dtype="<f8")
        stream.write(np.ravel(field, order="F").tobytes(order="C"))


def validate_conserved_field(data: np.ndarray, gamma: float) -> None:
    if not np.isfinite(data).all():
        raise ValueError("merged turbulence field contains NaN or infinite values")
    rho = data[..., 0]
    if np.any(rho <= 0.0):
        raise ValueError("merged turbulence field contains non-positive density")
    kinetic = 0.5 * (
        data[..., 1] ** 2 + data[..., 2] ** 2 + data[..., 3] ** 2
    ) / rho
    pressure = (gamma - 1.0) * (data[..., 4] - kinetic)
    if not np.isfinite(pressure).all() or np.any(pressure <= 0.0):
        raise ValueError("merged turbulence field contains non-positive pressure")


def prepare_imported_turbulence(
    input_path: Path,
    output_path: Path,
    *,
    meta_path: Path | None = None,
    step: str = "latest",
    layout: str = "auto",
    gamma: float = 1.4,
) -> tuple[int, tuple[int, int, int]]:
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    if gamma <= 1.0 or not np.isfinite(gamma):
        raise ValueError("gamma must be finite and greater than 1")

    if meta_path is None:
        meta_path = _default_meta_path(input_path)
    meta = load_meta(meta_path) if meta_path is not None else {}
    files = discover_files(input_path, layout=layout, meta=meta)
    if not files:
        raise ValueError(f"no SLF files found: {input_path}")

    rank_files = [path for path in files if parse_step_rank(path)[1] is not None]
    rank_lookup = build_rank_range_lookup(meta)
    if rank_files and not rank_lookup:
        raise ValueError(
            "rank-wise SLF input requires meta.json with parallel.rank_ranges"
        )

    groups = complete_step_groups(group_by_step(files), meta)
    selected_step = _select_step(groups, step)
    selected_files = groups[selected_step]
    first = read_slf(selected_files[0])
    origin, spacing = get_origin_spacing(meta, first)
    global_shape = get_global_grid(meta, first)
    output_shape = (*global_shape, len(REQUIRED_FIELDS))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as stream:
        data_offset = _write_header(
            stream,
            output_shape,
            list(REQUIRED_FIELDS),
            selected_step,
            first.time,
            origin,
            spacing,
        )
        total_bytes = int(np.prod(output_shape, dtype=np.int64)) * 8
        stream.truncate(data_offset + total_bytes)

    output = np.memmap(
        output_path,
        mode="r+",
        dtype="<f8",
        offset=data_offset,
        shape=output_shape,
        order="F",
    )
    output.fill(0.0)
    try:
        for index, path in enumerate(selected_files):
            slf = first if index == 0 else read_slf(path)
            step_read, rank = parse_step_rank(path, slf)
            if step_read != selected_step:
                raise RuntimeError("internal SLF step-selection mismatch")
            if abs(slf.time - first.time) > 1.0e-12 * max(1.0, abs(first.time)):
                raise ValueError(f"{path}: time does not match the other rank files")

            if rank is None and len(selected_files) == 1:
                if slf.shape[:3] != global_shape:
                    raise ValueError(
                        f"{path}: a global SLF must be ghost-free with shape {global_shape}"
                    )
                slices = (slice(0, global_shape[0]), slice(0, global_shape[1]), slice(0, global_shape[2]))
                local = slf.data
            else:
                slices = rank_slices(rank, slf, rank_lookup)
                expected = tuple(item.stop - item.start for item in slices)
                local = crop_center_to_shape(slf.data, expected, path)

            _names, local = select_variables(
                slf.names,
                local,
                mode="none",
                requested=list(REQUIRED_FIELDS),
                gamma=gamma,
            )
            validate_conserved_field(local, gamma)
            output[slices[0], slices[1], slices[2], :] = local

        output.flush()
        for k_start in range(0, global_shape[2], 8):
            k_end = min(k_start + 8, global_shape[2])
            validate_conserved_field(output[:, :, k_start:k_end, :], gamma)
    finally:
        del output
    return selected_step, global_shape


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Merge NSE rank-wise SLFs, remove ghost cells, and write one "
            "portable SLF for flow.type=imported_turbulence"
        )
    )
    parser.add_argument("input", type=Path, help="SLF file or output directory")
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--meta", type=Path, help="meta.json; inferred beside input")
    parser.add_argument("--step", default="latest", help="integer step or latest")
    parser.add_argument(
        "--layout", choices=("auto", "rank", "global"), default="auto"
    )
    parser.add_argument("--gamma", type=float, default=1.4)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        step, shape = prepare_imported_turbulence(
            args.input,
            args.output,
            meta_path=args.meta,
            step=args.step,
            layout=args.layout,
            gamma=args.gamma,
        )
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[ERROR] {exc}")
        return 1
    print(
        f"[OK] Prepared imported turbulence: {args.output.resolve()} "
        f"(step={step}, grid={shape[0]}x{shape[1]}x{shape[2]})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
