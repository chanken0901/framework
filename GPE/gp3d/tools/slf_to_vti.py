#!/usr/bin/env python3
"""単一のSolverLibrary SLFをVTK ImageData (.vti)へ変換する。

GPEの複素波動関数はpsi_realとpsi_imagの2変数として保存される。本ツールは
それらをVTIへ格納し、density=|psi|^2、abs_psi、phaseも追加する。
MPI rankファイルの全領域結合には別のmerged変換ツールを使用する。
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import os
import struct
import sys
import zlib
from array import array
from pathlib import Path


HEADER_INT = "<i"
HEADER_DOUBLE = "<d"


def read_exact(handle, nbytes: int) -> bytes:
    data = handle.read(nbytes)
    if len(data) != nbytes:
        raise EOFError(f"unexpected end of file while reading {nbytes} bytes")
    return data


def read_slf(path: Path) -> dict:
    """一つのSLFからヘッダー、変数名、float64データを読み込む。"""
    with path.open("rb") as handle:
        magic = read_exact(handle, 8)
        if magic != b"SLF1\0\0\0\0":
            raise ValueError(f"{path} is not an SLF1 file")

        version = struct.unpack(HEADER_INT, read_exact(handle, 4))[0]
        dtype_code = struct.unpack(HEADER_INT, read_exact(handle, 4))[0]
        ndim = struct.unpack(HEADER_INT, read_exact(handle, 4))[0]
        shape = struct.unpack("<4i", read_exact(handle, 16))
        meta = struct.unpack("<8i", read_exact(handle, 32))
        time = struct.unpack(HEADER_DOUBLE, read_exact(handle, 8))[0]
        bounds = struct.unpack("<6d", read_exact(handle, 48))
        nvar = struct.unpack(HEADER_INT, read_exact(handle, 4))[0]
        names = [
            read_exact(handle, 32).decode("ascii", errors="ignore").rstrip("\0 ").strip()
            for _ in range(nvar)
        ]

        if version != 1:
            raise ValueError(f"unsupported SLF version: {version}")
        if dtype_code != 2:
            raise ValueError(f"unsupported SLF dtype_code {dtype_code}; expected 2=float64")
        if ndim != 4:
            raise ValueError(f"unsupported SLF ndim {ndim}; expected 4")
        if shape[3] != nvar:
            raise ValueError(f"header mismatch: shape nvar={shape[3]}, nvar={nvar}")

        count = shape[0] * shape[1] * shape[2] * nvar
        raw = read_exact(handle, count * 8)

    values = array("d")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()

    return {
        "path": path,
        "shape": shape,
        "meta": meta,
        "time": time,
        "bounds": bounds,
        "names": names,
        "values": values,
    }


def variable_values(values: array, nx: int, ny: int, nz: int, nvar: int, ivar: int) -> array:
    """Extract one Fortran-ordered variable from field(nx,ny,nz,nvar)."""
    ncell = nx * ny * nz
    start = ivar * ncell
    out = array("d", values[start : start + ncell])
    return out


def derived_gpe_arrays(names: list[str], values: array, nx: int, ny: int, nz: int, nvar: int) -> dict[str, array]:
    """psi_real/psi_imagがある場合に密度、絶対値、位相を計算する。"""
    lowered = [name.lower() for name in names]
    result: dict[str, array] = {}

    for ivar, name in enumerate(names):
        result[name or f"var{ivar + 1}"] = variable_values(values, nx, ny, nz, nvar, ivar)

    if "psi_real" not in lowered or "psi_imag" not in lowered:
        return result

    real_part = result[names[lowered.index("psi_real")]]
    imag_part = result[names[lowered.index("psi_imag")]]
    density = array("d")
    phase = array("d")

    for re, im in zip(real_part, imag_part):
        density.append(re * re + im * im)
        phase.append(math.atan2(im, re))

    result["density"] = density
    result["phase"] = phase
    return result


def xml_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace('"', "&quot;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def encoded_data(values: array, compress: bool) -> str:
    raw = values.tobytes()
    if sys.byteorder != "little":
        tmp = array("d", values)
        tmp.byteswap()
        raw = tmp.tobytes()

    payload = zlib.compress(raw) if compress else raw
    header = struct.pack("<Q", len(payload))
    return base64.b64encode(header + payload).decode("ascii")


def write_vti(slf: dict, output_path: Path, compress: bool = False) -> None:
    """SLFの格子と変数をVTIのPointDataとして書き出す。"""
    nx, ny, nz, nvar = slf["shape"]
    x_min, x_max, y_min, y_max, z_min, z_max = slf["bounds"]

    spacing = (
        (x_max - x_min) / max(nx - 1, 1),
        (y_max - y_min) / max(ny - 1, 1),
        (z_max - z_min) / max(nz - 1, 1),
    )
    arrays = derived_gpe_arrays(slf["names"], slf["values"], nx, ny, nz, nvar)

    compressor_attr = ' compressor="vtkZLibDataCompressor"' if compress else ""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write(f'<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"{compressor_attr}>\n')
        out.write(
            f'  <ImageData WholeExtent="0 {nx - 1} 0 {ny - 1} 0 {nz - 1}" '
            f'Origin="{x_min:.17g} {y_min:.17g} {z_min:.17g}" '
            f'Spacing="{spacing[0]:.17g} {spacing[1]:.17g} {spacing[2]:.17g}">\n'
        )
        out.write(f'    <Piece Extent="0 {nx - 1} 0 {ny - 1} 0 {nz - 1}">\n')
        out.write('      <PointData Scalars="density">\n')
        for name, arr in arrays.items():
            out.write(
                f'        <DataArray type="Float64" Name="{xml_escape(name)}" '
                f'NumberOfComponents="1" format="binary">\n'
            )
            out.write(f"          {encoded_data(arr, compress)}\n")
            out.write("        </DataArray>\n")
        out.write("      </PointData>\n")
        out.write("      <CellData>\n")
        out.write("      </CellData>\n")
        out.write("    </Piece>\n")
        out.write("  </ImageData>\n")
        out.write("</VTKFile>\n")


def default_output_path(input_path: Path, output_dir: Path | None) -> Path:
    if output_dir is None:
        return input_path.with_suffix(".vti")
    return output_dir / f"{input_path.stem}.vti"


def load_meta(path: Path | None) -> dict | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert SLF1 files to VTI.")
    parser.add_argument("inputs", nargs="+", type=Path, help="SLF files or directories to convert")
    parser.add_argument("-o", "--output", type=Path, help="Output .vti path for one input, or output directory for many")
    parser.add_argument("--meta", type=Path, help="Optional meta.json path; parsed for validation only")
    parser.add_argument("--compress", action="store_true", help="Compress VTI binary arrays with zlib")
    args = parser.parse_args(argv)

    expanded_inputs: list[Path] = []
    meta = load_meta(args.meta)
    case_name = str(meta.get("case_name", "")).strip() if meta else ""
    for input_path in args.inputs:
        if input_path.is_dir():
            if case_name:
                case_files = sorted(input_path.glob(f"{case_name}_field_*.slf"))
                if case_files:
                    expanded_inputs.extend(case_files)
                    continue
            expanded_inputs.extend(sorted(input_path.glob("*.slf")))
        else:
            expanded_inputs.append(input_path)

    if not expanded_inputs:
        parser.error("no SLF files matched the input")

    if args.output and len(expanded_inputs) > 1 and args.output.suffix:
        parser.error("--output must be a directory when converting multiple files")

    if meta is not None and meta.get("equation", "").lower() != "gpe":
        print(f"warning: meta equation is {meta.get('equation')!r}, not 'GPE'", file=sys.stderr)

    output_dir = args.output if args.output and not args.output.suffix else None

    for input_path in expanded_inputs:
        slf = read_slf(input_path)
        output_path = args.output if args.output and len(expanded_inputs) == 1 and args.output.suffix else default_output_path(input_path, output_dir)
        write_vti(slf, output_path, compress=args.compress)
        print(f"{input_path} -> {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
