#!/usr/bin/env python3
#python slf_to_paraview_merged_cropghost.py output -o paraview --meta output/meta.json --derive nse
"""
Convert SolverLibrary rank-wise SLF files into full-domain ParaView VTI files.

This version MERGES all rank files belonging to the same step into one global
3D field before writing .vti. ParaView can then visualize the whole domain as a
single ImageData file.

Required for MPI output:
  meta.json must contain parallel.rank_ranges, e.g.

  "parallel": {
    "rank_ranges": [
      {"rank":0,"i_start":1,"i_end":64,"j_start":1,"j_end":32,"k_start":1,"k_end":32},
      ...
    ]
  }

Usage:
  python slf_to_paraview_merged.py output -o paraview --meta output/meta.json --derive nse
  python slf_to_paraview_merged.py output -o paraview --meta output/meta.json --derive gpe
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import struct
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class SLFData:
    path: Path
    version: int
    dtype_code: int
    shape: tuple[int, int, int, int]
    meta: np.ndarray
    time: float
    bounds: tuple[float, float, float, float, float, float]
    names: list[str]
    data: np.ndarray


def _read_exact(f, nbytes: int) -> bytes:
    b = f.read(nbytes)
    if len(b) != nbytes:
        raise EOFError(f"Unexpected EOF: wanted {nbytes} bytes, got {len(b)} bytes")
    return b


def read_slf(path: str | Path) -> SLFData:
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 8)
        if magic[:4] != b"SLF1":
            raise ValueError(f"{path}: not an SLF1 file. magic={magic!r}")

        version = struct.unpack("<i", _read_exact(f, 4))[0]
        dtype_code = struct.unpack("<i", _read_exact(f, 4))[0]
        _ndim = struct.unpack("<i", _read_exact(f, 4))[0]
        shape_arr = np.frombuffer(_read_exact(f, 4 * 4), dtype="<i4").copy()
        meta = np.frombuffer(_read_exact(f, 8 * 4), dtype="<i4").copy()
        time = struct.unpack("<d", _read_exact(f, 8))[0]
        bounds = struct.unpack("<6d", _read_exact(f, 6 * 8))
        nvar = struct.unpack("<i", _read_exact(f, 4))[0]

        nx, ny, nz, nvar_shape = map(int, shape_arr)
        if nvar != nvar_shape:
            raise ValueError(f"{path}: nvar mismatch: header nvar={nvar}, shape nvar={nvar_shape}")
        if dtype_code != 2:
            raise ValueError(f"{path}: unsupported dtype_code={dtype_code}; expected float64")

        names = []
        for _ in range(nvar):
            raw = _read_exact(f, 32)
            names.append(raw.split(b"\x00", 1)[0].decode("ascii", errors="ignore").strip())

        count = nx * ny * nz * nvar
        raw_data = np.fromfile(f, dtype="<f8", count=count)
        if raw_data.size != count:
            raise EOFError(f"{path}: data size mismatch: expected {count}, got {raw_data.size}")

    data = raw_data.reshape((nx, ny, nz, nvar), order="F")
    return SLFData(path, version, dtype_code, (nx, ny, nz, nvar), meta, time,
                   tuple(float(x) for x in bounds), names, data)


def parse_step_rank(path: Path, slf: SLFData | None = None) -> tuple[int, int | None]:
    stem = path.stem
    # Examples:
    #   field_000100_rank00003.slf
    #   3d_result_step_00000_rank00003.fbn converted/named as .slf
    #   step_000100_rank00003.slf
    rank = None
    mr = re.search(r"rank(\d+)", stem)
    if mr:
        rank = int(mr.group(1))

    ms = re.search(r"(?:field_|step_|step)(\d+)", stem)
    if ms:
        return int(ms.group(1)), rank

    # Fallback: last integer not belonging to rank
    cleaned = re.sub(r"rank\d+", "", stem)
    nums = re.findall(r"\d+", cleaned)
    if nums:
        return int(nums[-1]), rank

    if slf is not None:
        step = int(slf.meta[0])
        if rank is None and int(slf.meta[1]) >= 0:
            rank = int(slf.meta[1])
        return step, rank
    return 0, rank


def load_meta(meta_path: Path | None) -> dict:
    if meta_path is None:
        return {}
    if not meta_path.exists():
        raise FileNotFoundError(f"meta.json not found: {meta_path}")
    with meta_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def build_rank_range_lookup(meta: dict) -> dict[int, dict]:
    parallel = meta.get("parallel", {}) if isinstance(meta, dict) else {}
    ranges = parallel.get("rank_ranges", [])
    lookup: dict[int, dict] = {}
    for item in ranges:
        try:
            lookup[int(item["rank"])] = item
        except Exception:
            continue
    return lookup


def get_global_grid(meta: dict, first: SLFData) -> tuple[int, int, int]:
    grid = meta.get("grid")
    if isinstance(grid, list) and len(grid) >= 3:
        return int(grid[0]), int(grid[1]), int(grid[2])
    if isinstance(grid, dict):
        return int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
    # Serial fallback
    return first.shape[:3]


def get_origin_spacing(meta: dict, first: SLFData) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    origin_raw = meta.get("origin", [first.bounds[0], first.bounds[2], first.bounds[4]])
    origin = (float(origin_raw[0]), float(origin_raw[1]), float(origin_raw[2]))

    spacing_raw = meta.get("spacing")
    if spacing_raw is not None:
        spacing = (float(spacing_raw[0]), float(spacing_raw[1]), float(spacing_raw[2]))
    else:
        nx, ny, nz = get_global_grid(meta, first)
        domain = meta.get("domain_length", [
            first.bounds[1] - first.bounds[0],
            first.bounds[3] - first.bounds[2],
            first.bounds[5] - first.bounds[4],
        ])
        spacing = (float(domain[0]) / nx, float(domain[1]) / ny, float(domain[2]) / nz)
    return origin, spacing


def rank_slices(rank: int | None, slf: SLFData, meta: dict, rank_lookup: dict[int, dict]) -> tuple[slice, slice, slice]:
    if rank is None or rank not in rank_lookup:
        # Serial or no metadata: place at origin of global array.
        nx, ny, nz, _ = slf.shape
        return slice(0, nx), slice(0, ny), slice(0, nz)

    rr = rank_lookup[rank]
    ist = int(rr.get("i_start", 1)); ien = int(rr.get("i_end", slf.shape[0]))
    jst = int(rr.get("j_start", 1)); jen = int(rr.get("j_end", slf.shape[1]))
    kst = int(rr.get("k_start", 1)); ken = int(rr.get("k_end", slf.shape[2]))
    return slice(ist - 1, ien), slice(jst - 1, jen), slice(kst - 1, ken)


def add_derived_variables(names: list[str], data: np.ndarray, mode: str, gamma: float = 1.4) -> tuple[list[str], np.ndarray]:
    names = list(names)
    if mode == "none":
        return names, data

    derived = []
    derived_names = []

    if mode == "gpe":
        lookup = {name: i for i, name in enumerate(names)}
        if "psi_real" not in lookup or "psi_imag" not in lookup:
            raise ValueError("GPE derivation requires psi_real and psi_imag")
        re_part = data[..., lookup["psi_real"]]
        im_part = data[..., lookup["psi_imag"]]
        derived_names += ["rho", "phase"]
        derived += [re_part * re_part + im_part * im_part, np.arctan2(im_part, re_part)]

    elif mode == "nse":
        lookup = {name: i for i, name in enumerate(names)}
        required = ["rho", "rho_u", "rho_v", "rho_w", "rho_E"]
        if not all(v in lookup for v in required):
            raise ValueError(f"NSE derivation requires {required}; found {names}")
        rho = data[..., lookup["rho"]]
        rho_safe = np.where(np.abs(rho) > 1.0e-300, rho, np.nan)
        u = data[..., lookup["rho_u"]] / rho_safe
        v = data[..., lookup["rho_v"]] / rho_safe
        w = data[..., lookup["rho_w"]] / rho_safe
        E = data[..., lookup["rho_E"]]
        p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v + w * w))
        derived_names += ["u", "v", "w", "p"]
        derived += [u, v, w, p]
    else:
        raise ValueError(f"Unknown derive mode: {mode}")

    extra = np.stack(derived, axis=-1)
    return names + derived_names, np.concatenate([data, extra], axis=-1)


def crop_center_to_shape(arr: np.ndarray, expected: tuple[int, int, int], path: Path | None = None) -> np.ndarray:
    """Crop symmetric ghost cells so local data matches the physical rank range.

    The Fortran solver often writes q including ghost cells, e.g.
    local shape (70, 38, 38) for physical shape (64, 32, 32) with nghost=3.
    For ParaView full-domain visualization, only physical cells should be pasted
    into the global array. This function removes equal padding on both sides.
    """
    local = arr.shape[:3]
    slices = []
    for axis, (nlocal, nexp) in enumerate(zip(local, expected)):
        if nlocal == nexp:
            slices.append(slice(None))
            continue
        diff = nlocal - nexp
        if diff < 0 or diff % 2 != 0:
            loc = f"{path}: " if path is not None else ""
            raise ValueError(
                f"{loc}local data shape {local} cannot be cropped to meta rank range {expected}.\n"
                "The difference must be a non-negative even number on each axis. "
                "Check meta.json rank_ranges or write physical-cell data only."
            )
        g = diff // 2
        slices.append(slice(g, g + nexp))

    cropped = arr[slices[0], slices[1], slices[2], :]
    if cropped.shape[:3] != expected:
        loc = f"{path}: " if path is not None else ""
        raise ValueError(f"{loc}internal crop error: got {cropped.shape[:3]}, expected {expected}")
    return cropped


def merge_step(files: list[Path], meta: dict, rank_lookup: dict[int, dict], derive: str, gamma: float) -> tuple[int, float, list[str], np.ndarray, tuple[float, float, float], tuple[float, float, float]]:
    records = []
    for path in files:
        slf = read_slf(path)
        step, rank = parse_step_rank(path, slf)
        records.append((path, slf, step, rank))

    first = records[0][1]
    step = records[0][2]
    time = first.time
    gx, gy, gz = get_global_grid(meta, first)
    nvar = first.shape[3]
    names = first.names

    global_data = np.empty((gx, gy, gz, nvar), dtype=np.float64)
    global_data.fill(np.nan)
    filled = np.zeros((gx, gy, gz), dtype=bool)

    for path, slf, _step, rank in records:
        if _step != step:
            raise ValueError("Internal error: mixed steps in merge_step")
        if slf.shape[3] != nvar:
            raise ValueError(f"{path}: nvar mismatch; expected {nvar}, got {slf.shape[3]}")
        if slf.names != names:
            print(f"WARNING: variable names differ in {path}; using names from first file", file=sys.stderr)

        sx, sy, sz = rank_slices(rank, slf, meta, rank_lookup)
        expected = (sx.stop - sx.start, sy.stop - sy.start, sz.stop - sz.start)
        local_shape = slf.shape[:3]

        local_data = slf.data
        if local_shape != expected:
            print(
                f"INFO: {path.name}: cropping ghost cells {local_shape} -> {expected}",
                file=sys.stderr,
            )
            local_data = crop_center_to_shape(local_data, expected, path)

        global_data[sx, sy, sz, :] = local_data
        filled[sx, sy, sz] = True

    missing = int((~filled).sum())
    if missing:
        print(f"WARNING: merged step {step} has {missing} unfilled cells. They are NaN in the VTI.", file=sys.stderr)

    names, global_data = add_derived_variables(names, global_data, derive, gamma=gamma)
    origin, spacing = get_origin_spacing(meta, first)
    return step, time, names, global_data, origin, spacing


def _sanitize_name(name: str, fallback: str) -> str:
    name = name.strip() or fallback
    name = re.sub(r"[^0-9A-Za-z_]+", "_", name)
    if name[0].isdigit():
        name = "var_" + name
    return name


def write_vti(out_path: Path, names: list[str], data: np.ndarray, origin: tuple[float, float, float], spacing: tuple[float, float, float], time: float, step: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nx, ny, nz, nvar = data.shape
    x0, y0, z0 = origin
    dx, dy, dz = spacing
    whole_extent = f"0 {nx} 0 {ny} 0 {nz}"

    vtk = ET.Element("VTKFile", {"type": "ImageData", "version": "1.0", "byte_order": "LittleEndian", "header_type": "UInt64"})
    image = ET.SubElement(vtk, "ImageData", {"WholeExtent": whole_extent, "Origin": f"{x0:.17g} {y0:.17g} {z0:.17g}", "Spacing": f"{dx:.17g} {dy:.17g} {dz:.17g}"})
    piece = ET.SubElement(image, "Piece", {"Extent": whole_extent})

    field_data = ET.SubElement(piece, "FieldData")
    ET.SubElement(field_data, "DataArray", {"type": "Float64", "Name": "TimeValue", "NumberOfTuples": "1", "format": "ascii"}).text = f"{time:.17e}"
    ET.SubElement(field_data, "DataArray", {"type": "Int32", "Name": "Step", "NumberOfTuples": "1", "format": "ascii"}).text = str(int(step))

    container = ET.SubElement(piece, "CellData")
    ET.SubElement(piece, "PointData")

    arrays = []
    offset = 0
    for ivar in range(nvar):
        arr = np.asarray(data[..., ivar], dtype="<f8")
        flat = np.ravel(arr, order="F")
        arrays.append(flat)
        ET.SubElement(container, "DataArray", {
            "type": "Float64",
            "Name": _sanitize_name(names[ivar] if ivar < len(names) else "", f"var{ivar+1}"),
            "NumberOfComponents": "1",
            "format": "appended",
            "offset": str(offset),
        })
        offset += 8 + flat.nbytes

    appended = ET.SubElement(vtk, "AppendedData", {"encoding": "raw"})
    appended.text = "_"
    xml = ET.tostring(vtk, encoding="utf-8", xml_declaration=True, short_empty_elements=False)
    marker = b"_</AppendedData>"
    if marker not in xml:
        raise RuntimeError("Internal XML marker not found")
    prefix, suffix = xml.split(marker, 1)

    with out_path.open("wb") as f:
        f.write(prefix)
        f.write(b"_")
        for flat in arrays:
            f.write(struct.pack("<Q", flat.nbytes))
            f.write(flat.tobytes(order="C"))
        f.write(b"</AppendedData>")
        f.write(suffix)


def write_pvd(out_path: Path, datasets: list[tuple[float, Path]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    vtk = ET.Element("VTKFile", {"type": "Collection", "version": "0.1", "byte_order": "LittleEndian"})
    collection = ET.SubElement(vtk, "Collection")
    for time, path in datasets:
        ET.SubElement(collection, "DataSet", {"timestep": f"{time:.17g}", "group": "", "part": "0", "file": path.as_posix()})
    ET.ElementTree(vtk).write(out_path, encoding="utf-8", xml_declaration=True)


def discover_files(input_path: Path) -> list[Path]:
    if input_path.is_dir():
        files = sorted(Path(p) for p in glob.glob(str(input_path / "*.slf")))
    else:
        files = [input_path]
    return files


def group_by_step(files: list[Path]) -> dict[int, list[Path]]:
    groups: dict[int, list[Path]] = {}
    for p in files:
        step, _rank = parse_step_rank(p)
        groups.setdefault(step, []).append(p)
    return dict(sorted(groups.items()))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Merge rank-wise SolverLibrary .slf files into full-domain ParaView .vti files")
    parser.add_argument("input", help="Input .slf file or directory containing rank-wise .slf files")
    parser.add_argument("-o", "--output-dir", default="paraview", help="Output directory")
    parser.add_argument("--meta", default=None, help="Path to meta.json. Default: input_dir/meta.json or parent/meta.json if found")
    parser.add_argument("--derive", choices=["none", "nse", "gpe"], default="none")
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--pvd-name", default="collection.pvd")
    args = parser.parse_args(argv)

    in_path = Path(args.input)
    files = discover_files(in_path)
    if not files:
        print(f"No .slf files found: {in_path}", file=sys.stderr)
        return 1

    if args.meta:
        meta_path = Path(args.meta)
    else:
        candidates = []
        if in_path.is_dir():
            candidates = [in_path / "meta.json", in_path.parent / "meta.json"]
        else:
            candidates = [in_path.parent / "meta.json", in_path.parent.parent / "meta.json"]
        meta_path = next((p for p in candidates if p.exists()), None)

    meta = load_meta(meta_path) if meta_path is not None else {}
    rank_lookup = build_rank_range_lookup(meta)
    if len(files) > 1 and not rank_lookup:
        raise RuntimeError("Rank-wise merge requires meta.json with parallel.rank_ranges")

    groups = group_by_step(files)
    out_dir = Path(args.output_dir)
    datasets: list[tuple[float, Path]] = []

    for step, step_files in groups.items():
        step_i, time, names, data, origin, spacing = merge_step(step_files, meta, rank_lookup, args.derive, args.gamma)
        out = out_dir / f"field_{step_i:06d}.vti"
        write_vti(out, names, data, origin, spacing, time, step_i)
        datasets.append((time, out.relative_to(out_dir)))
        print(f"Wrote merged full-domain VTI: {out}  shape={data.shape}")

    pvd = out_dir / args.pvd_name
    write_pvd(pvd, datasets)
    print(f"Wrote: {pvd}")
    print("Open collection.pvd in ParaView to visualize the full computational domain.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
