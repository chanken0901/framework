#!/usr/bin/env python3
# python slf_to_paraview_merged_cropghost.py output -o paraview --meta output/meta.json --derive auto
"""rankごとのSLFを結合し、ParaView用の全領域VTIへ変換する。

同じ時刻ステップに属するrankファイルを ``meta.json`` の領域情報に従って
全体配列へ配置する。GPEではpsi_real/psi_imagからdensity、phase、abs_psiを
導出でき、NSEでは保存された保存変数からu、v、w、pを導出できる。
過去データにghost cellが含まれる場合は中央領域を切り出す。

Usage:
  python tools/slf_to_paraview_merged_cropghost.py output -o paraview --meta output/meta.json --derive auto
  python tools/slf_to_paraview_merged_cropghost.py output -o preview --derive nse --steps latest --fields rho,u,v,w,p --stride 2
  python tools/slf_to_paraview_merged_cropghost.py output -o paraview --derive nse --steps 0:1000:100
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


TOOL_VERSION = "2.1.0"


@dataclass
class SLFData:
    """SLFヘッダー情報とFortran順の変数配列をまとめる読み込み結果。"""
    path: Path
    version: int
    dtype_code: int
    shape: tuple[int, int, int, int]
    meta: np.ndarray
    time: float
    bounds: tuple[float, float, float, float, float, float]
    names: list[str]
    data: np.ndarray
    grid_header: tuple[int, ...] | None = None
    grid_values: tuple[float, ...] | None = None


def _read_exact(f, nbytes: int) -> bytes:
    b = f.read(nbytes)
    if len(b) != nbytes:
        raise EOFError(f"Unexpected EOF: wanted {nbytes} bytes, got {len(b)} bytes")
    return b


def read_slf(path: str | Path, *, allow_nonuniform: bool = False) -> SLFData:
    """単一SLFファイルの固定ヘッダー、変数名、float64配列を読み込む。"""
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
        if int(meta[7]) == 1 and not allow_nonuniform:
            raise ValueError("Nonuniform SLF: this postprocessor/FFT requires uniform spacing. "
                             "Use the ParaView converter for mapped visualization; standard FFT remains unsupported.")
        time = struct.unpack("<d", _read_exact(f, 8))[0]
        bounds = struct.unpack("<6d", _read_exact(f, 6 * 8))
        nvar = struct.unpack("<i", _read_exact(f, 4))[0]

        nx, ny, nz, nvar_shape = map(int, shape_arr)
        if min(nx, ny, nz, nvar_shape) <= 0:
            raise ValueError(f'{path}: invalid SLF shape')
        if nvar != nvar_shape:
            raise ValueError(f"{path}: nvar mismatch: header nvar={nvar}, shape nvar={nvar_shape}")
        if dtype_code != 2:
            raise ValueError(f"{path}: unsupported dtype_code={dtype_code}; expected float64")

        names = []
        for _ in range(nvar):
            raw = _read_exact(f, 32)
            names.append(raw.split(b"\x00", 1)[0].decode("ascii", errors="ignore").strip())

        count = nx * ny * nz * nvar
        if path.stat().st_size - f.tell() < count * 8:
            raise EOFError(f'{path}: truncated SLF field payload')
        raw_data = np.fromfile(f, dtype="<f8", count=count)
        if raw_data.size != count:
            raise EOFError(f"{path}: data size mismatch: expected {count}, got {raw_data.size}")
        physics = f.read(80)
        marker = f.read(8)
        grid_header = grid_values = None
        if marker == b"NSEGRID1" or int(meta[7]) == 1:
            if not allow_nonuniform:
                raise ValueError("Nonuniform SLF: uniform-grid postprocessing/FFT is not applicable")
            if marker != b"NSEGRID1" or not physics.startswith(b"NSEPAR1\0"):
                raise ValueError("Nonuniform SLF lacks grid/physical metadata")
            grid_header = struct.unpack("<15i", _read_exact(f, 60))
            grid_values = struct.unpack("<15d", _read_exact(f, 120))
            if _read_exact(f, 8) != b"NSREND1\0" or f.read(1):
                raise ValueError("Invalid mapped SLF end marker/length")

    data = raw_data.reshape((nx, ny, nz, nvar), order="F")
    result = SLFData(path, version, dtype_code, (nx, ny, nz, nvar), meta, time, tuple(float(x) for x in bounds), names, data, grid_header, grid_values)
    if grid_header is not None:
        validate_mapped_header(result, physics)
    return result


def validate_mapped_header(slf: SLFData, physics: bytes) -> None:
    h, r = slf.grid_header, slf.grid_values
    ng = int(slf.meta[5])
    valid = (h[0] == 1 and h[1] >= 0 and min(h[2:5]) >= 3 and
             1 <= h[5] <= h[6] <= h[3] and 1 <= h[7] <= h[8] <= h[4] and
             h[9] in (0, 1) and h[10] in (0, 1, 2) and h[11] == 1 and
             all(v in (0, 1) for v in h[12:15]) and 0 <= ng <= min(h[2:5]) and
             0 <= slf.meta[1] < slf.meta[6] <= 100000 and
             slf.version == 1 and slf.meta[7] == 1 and slf.shape[3] == 5 and
             [n.lower() for n in slf.names] == ['rho', 'rho_u', 'rho_v', 'rho_w', 'rho_e'])
    valid = valid and (slf.shape[:3] == (h[2]+2*ng, h[6]-h[5]+1+2*ng, h[8]-h[7]+1+2*ng))
    valid = valid and (tuple(slf.meta[2:5]) == h[2:5] and slf.meta[0] == h[1] and
                      slf.time == r[0] and slf.bounds == r[6:12] and
                      struct.unpack('<5d', physics[8:48]) == r[1:6])
    valid = valid and (np.all(np.isfinite(r)) and r[0] >= 0 and r[1] > 1 and
                      r[2] >= 0 and min(r[3:5]) > 0 and r[5] >= 0 and
                      all(r[i+1] > r[i] for i in (6, 8, 10)) and
                      all(0 <= v <= 20 for v in r[12:15]))
    if valid:
        viscosity = physics[48:80].split(b'\0', 1)[0].decode('ascii').strip()
        valid = viscosity == ('none', 'fv2', 'central6')[h[10]] and (h[10] == 0 or r[2] > 0)
    if not valid:
        raise ValueError(f'{slf.path}: invalid/inconsistent mapped SLF metadata')


def mapped_axes(slf: SLFData):
    """Physical edges, volumes and solver quadrature, matching mod_grid_axis."""
    h, r = slf.grid_header, slf.grid_values
    edges, widths, metrics = [], [], []
    for a, n in enumerate(h[2:5]):
        lower, upper = r[6+2*a:8+2*a]
        beta = r[12+a]
        s = np.arange(n+1, dtype=float)/n
        fraction = s if beta < np.sqrt(np.finfo(float).eps) else .5*(1+np.sinh(beta*(2*s-1))/np.sinh(beta))
        e = lower+(upper-lower)*fraction
        e[0], e[-1] = lower, upper
        d = np.diff(e)
        extended = np.empty(n+7)
        extended[3:n+4] = e
        for j in range(1, 4):
            extended[3-j] = extended[4-j] - (d[n-j] if h[12+a] else d[j-1])
            extended[n+3+j] = extended[n+2+j] + (d[j-1] if h[12+a] else d[n-j])
        c = extended[:-1]+np.diff(extended)/2
        i = np.arange(n)+3
        metric = .75*(c[i+1]-c[i-1])-.15*(c[i+2]-c[i-2])+(c[i+3]-c[i-3])/60
        if np.any(d <= 0) or (h[9] and np.any(metric <= 0)):
            raise ValueError('Nonpositive mapped grid width/metric')
        edges.append(e); widths.append(d); metrics.append(metric if h[9] else d)
    product = lambda v: v[0][:, None, None]*v[1][None, :, None]*v[2][None, None, :]
    return edges, product(widths), product(metrics)


def merge_mapped(files: list[Path]):
    """Strict rank merge using SLF metadata, without relying on meta.json."""
    first = read_slf(files[0], allow_nonuniform=True)
    if first.grid_header is None:
        raise ValueError('Expected mapped SLF')
    h = first.grid_header
    if len(files) != first.meta[6]:
        raise ValueError('Missing mapped SLF rank files')
    result = np.empty((*h[2:5], 5))
    filled = np.zeros(h[3:5], dtype=bool)
    ranks = set()
    common = (0, 1, 2, 3, 4, 9, 10, 11, 12, 13, 14)
    for index, path in enumerate(files):
        slf = first if index == 0 else read_slf(path, allow_nonuniform=True)
        p = slf.grid_header
        if (p is None or tuple(p[i] for i in common) != tuple(h[i] for i in common) or
                slf.grid_values != first.grid_values or slf.meta[6] != first.meta[6]):
            raise ValueError('Mixed grid, step or physical parameters in mapped SLF ranks')
        rank = int(slf.meta[1])
        ys, zs = slice(p[5]-1, p[6]), slice(p[7]-1, p[8])
        if rank in ranks or np.any(filled[ys, zs]):
            raise ValueError('Duplicate/overlapping mapped SLF ranks')
        ranks.add(rank)
        local = crop_center_to_shape(slf.data, (h[2], p[6]-p[5]+1, p[8]-p[7]+1), path)
        if not np.all(np.isfinite(local)):
            raise ValueError('Nonfinite mapped SLF state')
        result[:, ys, zs, :] = local
        filled[ys, zs] = True
    if not np.all(filled):
        raise ValueError('Missing mapped SLF domain coverage')
    return first, result


def write_mapped_vtr(path: Path, first: SLFData, data: np.ndarray, mode: str,
                     requested: list[str] | None, stride: int) -> None:
    if stride != 1:
        raise ValueError('Nonuniform SLF currently requires --stride 1 (no implicit resampling)')
    edges, volume, weight = mapped_axes(first)
    names, data = select_variables(first.names, data, mode, requested, gamma=first.grid_values[1])
    extent = ' '.join(f'0 {n}' for n in first.grid_header[2:5])
    vtk = ET.Element('VTKFile', type='RectilinearGrid', version='1.0', byte_order='LittleEndian', header_type='UInt64')
    grid = ET.SubElement(vtk, 'RectilinearGrid', WholeExtent=extent)
    field = ET.SubElement(grid, 'FieldData')
    ET.SubElement(field, 'DataArray', type='Float64', Name='TimeValue', NumberOfTuples='1', format='ascii').text = format(first.time, '.17g')
    ET.SubElement(field, 'DataArray', type='Int32', Name='Step', NumberOfTuples='1', format='ascii').text = str(first.grid_header[1])
    piece = ET.SubElement(grid, 'Piece', Extent=extent)
    cell = ET.SubElement(piece, 'CellData')
    ET.SubElement(piece, 'PointData')
    arrays = [data[..., i] for i in range(len(names))]+[volume, weight]
    labels = names+['cell_volume', 'integration_weight']
    offset = 0
    for label, array in zip(labels, arrays):
        ET.SubElement(cell, 'DataArray', type='Float64', Name=label, format='appended', offset=str(offset))
        offset += 8+array.size*8
    coordinates = ET.SubElement(piece, 'Coordinates')
    for a, e in enumerate(edges):
        ET.SubElement(coordinates, 'DataArray', type='Float64', Name='xyz'[a], format='appended', offset=str(offset))
        offset += 8+e.size*8
    ET.SubElement(vtk, 'AppendedData', encoding='raw').text = '_'
    xml = ET.tostring(vtk, encoding='utf-8', xml_declaration=True, short_empty_elements=False)
    prefix, suffix = xml.split(b'_</AppendedData>', 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('wb') as f:
        f.write(prefix+b'_')
        for array in arrays+edges:
            raw = np.asarray(array, dtype='<f8').ravel(order='F').tobytes()
            f.write(struct.pack('<Q', len(raw))); f.write(raw)
        f.write(b'</AppendedData>'+suffix)


def parse_step_rank(path: Path, slf: SLFData | None = None) -> tuple[int, int | None]:
    stem = path.stem
    rank = None
    mr = re.search(r"rank(\d+)", stem)
    if mr:
        rank = int(mr.group(1))

    ms = re.search(r"(?:field_|step_|step)(\d+)", stem)
    if ms:
        return int(ms.group(1)), rank

    cleaned = re.sub(r"rank\d+", "", stem)
    nums = re.findall(r"\d+", cleaned)
    if nums:
        return int(nums[-1]), rank

    if slf is not None:
        step = int(slf.meta[0])
        if rank is None and int(slf.meta[1]) > 0:
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


def infer_derive_mode(requested: str, meta: dict, first: SLFData | None = None) -> str:
    if requested != "auto":
        return requested

    equation = str(meta.get("equation", "")).lower()
    if equation == "gpe":
        return "gpe"
    if equation == "nse":
        return "nse"

    if first is not None:
        names = {name.lower() for name in first.names}
        if {"psi_real", "psi_imag"} <= names:
            return "gpe"
        if {"rho", "rho_u", "rho_v", "rho_w", "rho_e"} <= names:
            return "nse"
    return "none"


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
    if first.grid_header is not None:
        return first.grid_header[2:5]
    grid = meta.get("grid")
    if isinstance(grid, list) and len(grid) >= 3:
        return int(grid[0]), int(grid[1]), int(grid[2])
    if isinstance(grid, dict):
        return int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
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


def rank_slices(rank: int | None, slf: SLFData, rank_lookup: dict[int, dict]) -> tuple[slice, slice, slice]:
    if rank is None or rank not in rank_lookup:
        nx, ny, nz, _ = slf.shape
        return slice(0, nx), slice(0, ny), slice(0, nz)

    rr = rank_lookup[rank]
    ist = int(rr.get("i_start", 1)); ien = int(rr.get("i_end", slf.shape[0]))
    jst = int(rr.get("j_start", 1)); jen = int(rr.get("j_end", slf.shape[1]))
    kst = int(rr.get("k_start", 1)); ken = int(rr.get("k_end", slf.shape[2]))
    return slice(ist - 1, ien), slice(jst - 1, jen), slice(kst - 1, ken)


def add_derived_variables(names: list[str], data: np.ndarray, mode: str, gamma: float = 1.4) -> tuple[list[str], np.ndarray]:
    """保存済み一次変数からGPEまたはNSEの可視化用変数を導出する。"""
    names = list(names)
    if mode == "none":
        return names, data

    derived = []
    derived_names = []
    lookup = {name.lower(): i for i, name in enumerate(names)}

    if mode == "gpe":
        if "psi_real" not in lookup or "psi_imag" not in lookup:
            raise ValueError(f"GPE derivation requires psi_real and psi_imag; found {names}")
        re_part = data[..., lookup["psi_real"]]
        im_part = data[..., lookup["psi_imag"]]
        density = re_part * re_part + im_part * im_part
        derived_names += ["density", "phase", "abs_psi"]
        derived += [density, np.arctan2(im_part, re_part), np.sqrt(density)]

    elif mode == "nse":
        required = ["rho", "rho_u", "rho_v", "rho_w", "rho_e"]
        if not all(v in lookup for v in required):
            raise ValueError(f"NSE derivation requires {required}; found {names}")
        rho = data[..., lookup["rho"]]
        rho_safe = np.where(np.abs(rho) > 1.0e-300, rho, np.nan)
        u = data[..., lookup["rho_u"]] / rho_safe
        v = data[..., lookup["rho_v"]] / rho_safe
        w = data[..., lookup["rho_w"]] / rho_safe
        E = data[..., lookup["rho_e"]]
        p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v + w * w))
        derived_names += ["u", "v", "w", "p"]
        derived += [u, v, w, p]
    else:
        raise ValueError(f"Unknown derive mode: {mode}")

    extra = np.stack(derived, axis=-1)
    return names + derived_names, np.concatenate([data, extra], axis=-1)


def select_variables(
    names: list[str],
    data: np.ndarray,
    mode: str,
    requested: list[str] | None,
    gamma: float = 1.4,
) -> tuple[list[str], np.ndarray]:
    if requested is None:
        return add_derived_variables(names, data, mode, gamma=gamma)

    lookup = {name.lower(): i for i, name in enumerate(names)}
    derived: dict[str, np.ndarray] = {}

    if mode == "gpe":
        if "psi_real" not in lookup or "psi_imag" not in lookup:
            raise ValueError(f"GPE derivation requires psi_real and psi_imag; found {names}")
        re_part = data[..., lookup["psi_real"]]
        im_part = data[..., lookup["psi_imag"]]
        density = re_part * re_part + im_part * im_part
        derived = {
            "density": density,
            "phase": np.arctan2(im_part, re_part),
            "abs_psi": np.sqrt(density),
        }
    elif mode == "nse":
        required = ["rho", "rho_u", "rho_v", "rho_w", "rho_e"]
        if not all(v in lookup for v in required):
            raise ValueError(f"NSE derivation requires {required}; found {names}")
        rho = data[..., lookup["rho"]]
        rho_safe = np.where(np.abs(rho) > 1.0e-300, rho, np.nan)
        u = data[..., lookup["rho_u"]] / rho_safe
        v = data[..., lookup["rho_v"]] / rho_safe
        w = data[..., lookup["rho_w"]] / rho_safe
        energy = data[..., lookup["rho_e"]]
        derived = {
            "u": u,
            "v": v,
            "w": w,
            "p": (gamma - 1.0) * (energy - 0.5 * rho * (u * u + v * v + w * w)),
        }

    selected_names: list[str] = []
    selected_arrays: list[np.ndarray] = []
    for requested_name in requested:
        key = requested_name.lower()
        if key in lookup:
            selected_names.append(names[lookup[key]])
            selected_arrays.append(data[..., lookup[key]])
        elif key in derived:
            selected_names.append(key)
            selected_arrays.append(derived[key])
        else:
            available = list(names) + list(derived)
            raise ValueError(f"Unknown output field {requested_name!r}; available fields: {available}")

    return selected_names, np.stack(selected_arrays, axis=-1)


def crop_center_to_shape(arr: np.ndarray, expected: tuple[int, int, int], path: Path | None = None) -> np.ndarray:
    local = arr.shape[:3]
    slices = []
    for nlocal, nexp in zip(local, expected):
        if nlocal == nexp:
            slices.append(slice(None))
            continue
        diff = nlocal - nexp
        if diff < 0 or diff % 2 != 0:
            loc = f"{path}: " if path is not None else ""
            raise ValueError(
                f"{loc}local data shape {local} cannot be cropped to meta rank range {expected}. "
                "The difference must be a non-negative even number on each axis."
            )
        g = diff // 2
        slices.append(slice(g, g + nexp))
    return arr[slices[0], slices[1], slices[2], :]


def merge_step(
    files: list[Path],
    meta: dict,
    rank_lookup: dict[int, dict],
    derive: str,
    gamma: float,
    requested_fields: list[str] | None,
    stride: int,
) -> tuple[int, float, list[str], np.ndarray, tuple[float, float, float], tuple[float, float, float]]:
    """一つのステップに属する全rankを逐次読み込み、全体3次元場へ結合する。"""
    first_path = files[0]
    first = read_slf(first_path)
    step, first_rank = parse_step_rank(first_path, first)
    time = first.time
    gx, gy, gz = get_global_grid(meta, first)
    nvar = first.shape[3]
    names = first.names

    is_single_global = len(files) == 1 and first_rank is None and first.shape[:3] == (gx, gy, gz)
    if is_single_global:
        global_data = first.data
    else:
        global_data = np.empty((gx, gy, gz, nvar), dtype=np.float64)
        global_data.fill(np.nan)
        filled = np.zeros((gx, gy, gz), dtype=bool)

        for index, path in enumerate(files):
            slf = first if index == 0 else read_slf(path)
            step_read, rank = parse_step_rank(path, slf)
            if step_read != step:
                raise ValueError("Internal error: mixed steps in merge_step")
            if slf.shape[3] != nvar:
                raise ValueError(f"{path}: nvar mismatch; expected {nvar}, got {slf.shape[3]}")
            if [x.lower() for x in slf.names] != [x.lower() for x in names]:
                print(f"WARNING: variable names differ in {path}; using names from first file", file=sys.stderr)

            sx, sy, sz = rank_slices(rank, slf, rank_lookup)
            expected = (sx.stop - sx.start, sy.stop - sy.start, sz.stop - sz.start)
            local_data = slf.data
            if slf.shape[:3] != expected:
                print(f"INFO: {path.name}: cropping ghost cells {slf.shape[:3]} -> {expected}", file=sys.stderr)
                local_data = crop_center_to_shape(local_data, expected, path)

            global_data[sx, sy, sz, :] = local_data
            filled[sx, sy, sz] = True

        missing = int((~filled).sum())
        if missing:
            print(f"WARNING: merged step {step} has {missing} unfilled cells. They are NaN in the VTI.", file=sys.stderr)

    sampled_data = global_data[::stride, ::stride, ::stride, :]
    names, sampled_data = select_variables(names, sampled_data, derive, requested_fields, gamma=gamma)
    origin, spacing = get_origin_spacing(meta, first)
    # SLF bounds/origin describe cell centers; VTK ImageData origin is the lower point boundary.
    origin = tuple(value - 0.5 * delta for value, delta in zip(origin, spacing))
    spacing = tuple(value * stride for value in spacing)
    return step, time, names, sampled_data, origin, spacing


def _sanitize_name(name: str, fallback: str) -> str:
    name = name.strip() or fallback
    name = re.sub(r"[^0-9A-Za-z_]+", "_", name)
    if name[0].isdigit():
        name = "var_" + name
    return name


def write_vti(
    out_path: Path,
    names: list[str],
    data: np.ndarray,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    time: float,
    step: int,
) -> None:
    """結合済みPointDataをVTK ImageData形式で書き出す。"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nx, ny, nz, nvar = data.shape
    x0, y0, z0 = origin
    dx, dy, dz = spacing
    whole_extent = f"0 {nx} 0 {ny} 0 {nz}"

    vtk = ET.Element("VTKFile", {"type": "ImageData", "version": "1.0", "byte_order": "LittleEndian", "header_type": "UInt64"})
    image = ET.SubElement(vtk, "ImageData", {
        "WholeExtent": whole_extent,
        "Origin": f"{x0:.17g} {y0:.17g} {z0:.17g}",
        "Spacing": f"{dx:.17g} {dy:.17g} {dz:.17g}",
    })
    piece = ET.SubElement(image, "Piece", {"Extent": whole_extent})

    field_data = ET.SubElement(piece, "FieldData")
    ET.SubElement(field_data, "DataArray", {
        "type": "Float64", "Name": "TimeValue", "NumberOfTuples": "1", "format": "ascii",
    }).text = f"{time:.17e}"
    ET.SubElement(field_data, "DataArray", {
        "type": "Int32", "Name": "Step", "NumberOfTuples": "1", "format": "ascii",
    }).text = str(int(step))

    container = ET.SubElement(piece, "CellData", {"Scalars": "density" if "density" in names else ""})
    ET.SubElement(piece, "PointData")

    array_nbytes = nx * ny * nz * np.dtype("<f8").itemsize
    offset = 0
    for ivar in range(nvar):
        ET.SubElement(container, "DataArray", {
            "type": "Float64",
            "Name": _sanitize_name(names[ivar] if ivar < len(names) else "", f"var{ivar + 1}"),
            "NumberOfComponents": "1",
            "format": "appended",
            "offset": str(offset),
        })
        offset += 8 + array_nbytes

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
        for ivar in range(nvar):
            flat = np.ravel(np.asarray(data[..., ivar], dtype="<f8"), order="F")
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


def discover_files(
    input_path: Path,
    case_name: str | None = None,
    layout: str = "auto",
    meta: dict | None = None,
) -> list[Path]:
    if input_path.is_dir():
        if case_name:
            files = sorted(Path(p) for p in glob.glob(str(input_path / f"{case_name}_field_*.slf")))
            if files:
                return files
        files = sorted(Path(p) for p in glob.glob(str(input_path / "*.slf")))
        rank_files = [path for path in files if re.search(r"rank\d+", path.stem)]
        global_files = [path for path in files if path not in rank_files]

        if layout == "rank":
            return rank_files
        if layout == "global":
            return global_files
        if not rank_files or not global_files:
            return rank_files if rank_files else global_files

        parallel = (meta or {}).get("parallel", {})
        decomposition = str(parallel.get("decomposition", "")).strip().lower()
        mpi_enabled = bool(parallel.get("mpi_enabled", False))
        if decomposition == "serial-global" or not mpi_enabled:
            selected, selected_name = global_files, "global"
        else:
            selected, selected_name = rank_files, "rank-wise"
        print(
            f"INFO: found both global and rank-wise SLF files; meta.json selects {selected_name} files. "
            "Use --layout global or --layout rank to override.",
            file=sys.stderr,
        )
        return selected
    return [input_path]


def group_by_step(files: list[Path]) -> dict[int, list[Path]]:
    groups: dict[int, list[Path]] = {}
    for p in files:
        step, _rank = parse_step_rank(p)
        groups.setdefault(step, []).append(p)
    return dict(sorted(groups.items()))


def complete_step_groups(
    groups: dict[int, list[Path]],
    meta: dict,
) -> dict[int, list[Path]]:
    """実行中にまだ全rankが書き終わっていないステップを除外する。"""
    if all(
        parse_step_rank(path)[1] is None
        for files in groups.values()
        for path in files
    ):
        return groups
    parallel = meta.get("parallel", {}) if isinstance(meta, dict) else {}
    expected = int(parallel.get("mpi_nprocs", 1))
    if expected <= 1:
        return groups

    complete: dict[int, list[Path]] = {}
    for step, files in groups.items():
        ranks = {parse_step_rank(path)[1] for path in files}
        ranks.discard(None)
        if len(files) == expected and len(ranks) == expected:
            complete[step] = files
        else:
            print(
                f"WARNING: skipping incomplete step {step}: "
                f"found {len(files)} files/{len(ranks)} ranks, expected {expected}",
                file=sys.stderr,
            )
    return complete


def select_step_groups(
    groups: dict[int, list[Path]],
    specification: str,
) -> dict[int, list[Path]]:
    """all、latest、カンマ区切り番号、start:stop:strideを解釈する。"""
    if not groups:
        return {}
    spec = specification.strip().lower()
    if spec == "all":
        return groups
    if spec == "latest":
        latest = max(groups)
        return {latest: groups[latest]}

    available = sorted(groups)
    selected: set[int] = set()
    for token in (item.strip() for item in spec.split(",")):
        if not token:
            continue
        if ":" not in token:
            selected.add(int(token))
            continue

        parts = token.split(":")
        if len(parts) not in {2, 3}:
            raise ValueError(f"Invalid step range: {token!r}")
        start = int(parts[0]) if parts[0] else available[0]
        stop = int(parts[1]) if parts[1] else available[-1]
        stride = int(parts[2]) if len(parts) == 3 and parts[2] else 1
        if stride <= 0:
            raise ValueError(f"Step stride must be positive: {token!r}")
        selected.update(
            step
            for step in available
            if start <= step <= stop and (step - start) % stride == 0
        )

    missing = sorted(step for step in selected if step not in groups)
    if missing:
        raise ValueError(
            f"Requested steps are not available: {missing}; available={available}"
        )
    result = {step: groups[step] for step in available if step in selected}
    if not result:
        raise ValueError(
            f"No steps matched {specification!r}; available={available}"
        )
    return result


def print_inspection(
    input_path: Path,
    meta_path: Path | None,
    files: list[Path],
    groups: dict[int, list[Path]],
    first: SLFData,
    meta: dict,
    derive_mode: str,
    requested_fields: list[str] | None,
    stride: int,
    layout: str,
) -> None:
    grid = get_global_grid(meta, first)
    steps = list(groups)
    rank_wise = any(parse_step_rank(path)[1] is not None for path in files)
    print(f"[INFO] SLF ParaView converter version: {TOOL_VERSION}")
    print(f"[INFO] input: {input_path.resolve()}")
    print(f"[INFO] meta: {meta_path.resolve() if meta_path is not None else 'not found'}")
    print(f"[INFO] layout: {layout} ({'rank-wise' if rank_wise else 'global'})")
    print(f"[INFO] SLF files: {len(files)}")
    print(f"[INFO] selected steps: {steps}")
    print(f"[INFO] source fields: {first.names}")
    print(f"[INFO] derive mode: {derive_mode}")
    print(f"[INFO] requested fields: {requested_fields or 'all'}")
    print(f"[INFO] global grid: {grid}; spatial stride: {stride}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Merge SolverLibrary .slf files into full-domain ParaView .vti files")
    parser.add_argument("--version", action="version", version=f"%(prog)s {TOOL_VERSION}")
    parser.add_argument("input", help="Input .slf file or directory containing .slf files")
    parser.add_argument("-o", "--output-dir", default="paraview", help="Output directory")
    parser.add_argument("--meta", default=None, help="Path to meta.json. Default: input_dir/meta.json or parent/meta.json if found")
    parser.add_argument("--derive", choices=["auto", "none", "nse", "gpe"], default="auto")
    parser.add_argument(
        "--layout",
        choices=["auto", "global", "rank"],
        default="auto",
        help="Select global or rank-wise SLFs when both exist. Default: infer from meta.json",
    )
    parser.add_argument(
        "--fields",
        default="all",
        help="Comma-separated output fields, for example density,phase. Default: all primary and derived fields",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Keep every Nth grid cell along each axis. Use 2 for a lighter preview. Default: 1",
    )
    parser.add_argument(
        "--steps",
        default="all",
        help=(
            "Steps to convert: all, latest, comma-separated values, or "
            "inclusive start:stop:stride. Example: 0:1000:100"
        ),
    )
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--pvd-name", default="collection.pvd")
    parser.add_argument(
        "--inspect-only",
        action="store_true",
        help="Validate input and show the selected SLF files/steps without writing VTI/PVD files",
    )
    args = parser.parse_args(argv)

    if args.stride < 1:
        parser.error("--stride must be at least 1")
    requested_fields = None
    if args.fields.strip().lower() != "all":
        requested_fields = [name.strip() for name in args.fields.split(",") if name.strip()]
        if not requested_fields:
            parser.error("--fields must be 'all' or a comma-separated list")

    in_path = Path(args.input)
    if not in_path.exists():
        parser.error(f"Input path does not exist: {in_path}")
    if args.meta:
        meta_path = Path(args.meta)
        if not meta_path.exists():
            parser.error(f"meta.json does not exist: {meta_path}")
    else:
        candidates = [in_path / "meta.json", in_path.parent / "meta.json"] if in_path.is_dir() else [
            in_path.parent / "meta.json",
            in_path.parent.parent / "meta.json",
        ]
        meta_path = next((p for p in candidates if p.exists()), None)

    meta = load_meta(meta_path) if meta_path is not None else {}
    case_name = str(meta.get("case_name", "")).strip() or None
    files = discover_files(in_path, case_name=case_name, layout=args.layout, meta=meta)
    if not files:
        print(f"No .slf files found: {in_path}", file=sys.stderr)
        if case_name:
            print(f"Looked for {case_name}_field_*.slf, then *.slf", file=sys.stderr)
        return 1
    try:
        first = read_slf(files[0], allow_nonuniform=True)
    except (OSError, EOFError, ValueError) as exc:
        parser.error(str(exc))
    derive_mode = infer_derive_mode(args.derive, meta, first)
    rank_lookup = build_rank_range_lookup(meta)

    rank_wise = any(parse_step_rank(path)[1] is not None for path in files)
    if first.grid_header is not None and args.stride != 1:
        parser.error('Nonuniform SLF requires --stride 1')
    if rank_wise and not rank_lookup and first.grid_header is None:
        parser.error(
            "Rank-wise SLF conversion requires meta.json with parallel.rank_ranges. "
            "Pass --meta <output/meta.json>."
        )

    all_groups = group_by_step(files)
    groups = all_groups if first.grid_header is not None else complete_step_groups(all_groups, meta)
    if not groups:
        parser.error(
            "No complete time steps are available. The simulation may still be writing rank files, "
            "or meta.json may contain the wrong mpi_nprocs value."
        )
    try:
        groups = select_step_groups(groups, args.steps)
    except ValueError as exc:
        parser.error(str(exc))
    print_inspection(
        in_path,
        meta_path,
        files,
        groups,
        first,
        meta,
        derive_mode,
        requested_fields,
        args.stride,
        args.layout,
    )
    if args.inspect_only:
        if first.grid_header is not None:
            try:
                for step_files in groups.values():
                    snapshot, _ = merge_mapped(step_files)
                    mapped_axes(snapshot)
            except (ValueError, EOFError, OSError) as exc:
                parser.error(str(exc))
        print("[OK] Inspection completed; no VTI/PVD files were written.")
        return 0
    out_dir = Path(args.output_dir)
    datasets: list[tuple[float, Path]] = []

    for step, step_files in groups.items():
        if first.grid_header is not None:
            try:
                snapshot, data = merge_mapped(step_files)
                out = out_dir / f'field_{snapshot.grid_header[1]:06d}.vtr'
                write_mapped_vtr(out, snapshot, data, derive_mode, requested_fields, args.stride)
            except (ValueError, EOFError, OSError) as exc:
                parser.error(str(exc))
            datasets.append((snapshot.time, out.relative_to(out_dir)))
            print(f'Wrote mapped-grid VTR: {out}')
            continue
        step_i, time, names, data, origin, spacing = merge_step(
            step_files,
            meta,
            rank_lookup,
            derive_mode,
            args.gamma,
            requested_fields,
            args.stride,
        )
        out = out_dir / f"field_{step_i:06d}.vti"
        write_vti(out, names, data, origin, spacing, time, step_i)
        datasets.append((time, out.relative_to(out_dir)))
        print(
            f"Wrote merged full-domain VTI: {out}  shape={data.shape}  "
            f"derive={derive_mode}  stride={args.stride}  fields={names}"
        )

    pvd = out_dir / args.pvd_name
    write_pvd(pvd, datasets)
    print(f"Wrote: {pvd}")
    print("Open collection.pvd in ParaView to visualize the full computational domain.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
