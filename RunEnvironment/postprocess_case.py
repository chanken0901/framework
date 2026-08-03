#!/usr/bin/env python3
"""Convert the generated case's SLF output into ParaView VTI/PVD files."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shlex
import subprocess
import sys
from pathlib import Path


class PostprocessError(RuntimeError):
    """Raised when the generated run environment is incomplete."""


def _load_json(path: Path) -> dict:
    if not path.is_file():
        raise PostprocessError(f"required file was not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise PostprocessError(f"expected a JSON object: {path}")
    return value


def _resolve(root: Path, value: str | None, default: Path) -> Path:
    path = Path(value) if value else default
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def _converter_path(root: Path, model: str) -> Path:
    candidates = {
        "nse": [
            root / "SolverLibrary" / "NSE" / "tools"
            / "slf_to_paraview_merged_cropghost.py",
        ],
        "gpe": [
            root / "SolverLibrary" / "GPE" / "gp3d" / "tools"
            / "slf_to_paraview_merged_cropghost.py",
        ],
    }
    for path in candidates.get(model, []):
        if path.is_file():
            return path
    raise PostprocessError(
        f"ParaView converter is not available for model {model!r}. "
        "Regenerate the execution environment from the current FrameWork."
    )


def _default_fields(model: str) -> str:
    if model == "nse":
        return "rho,u,v,w,p"
    if model == "gpe":
        return "density,phase"
    return "all"


def _display_command(command: list[str]) -> str:
    if sys.platform == "win32":
        return subprocess.list2cmdline(command)
    return shlex.join(command)


def build_command(args: argparse.Namespace) -> tuple[Path, list[str]]:
    root = Path(__file__).resolve().parents[1]
    lock = _load_json(root / "environment.lock.json")
    model = str(lock.get("model", "")).strip().lower()
    if model not in {"nse", "gpe"}:
        raise PostprocessError(f"unsupported model in environment.lock.json: {model}")

    case_directory_value = str(lock.get("case_directory", "")).strip()
    if case_directory_value:
        case_directory = Path(case_directory_value)
    else:
        case_id = str(lock.get("case_id", "")).strip()
        if not case_id:
            raise PostprocessError("case_directory/case_id is missing from environment.lock.json")
        case_directory = Path("cases") / case_id
    case_root = _resolve(root, args.case_directory, case_directory)

    input_dir = _resolve(root, args.input_dir, case_root / "output")
    output_dir = _resolve(root, args.output_dir, case_root / "paraview")
    meta = _resolve(root, args.meta, input_dir / "meta.json")
    converter = _converter_path(root, model)
    fields = args.fields if args.fields is not None else _default_fields(model)

    command = [
        sys.executable,
        str(converter),
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--meta",
        str(meta),
        "--derive",
        args.derive,
        "--layout",
        args.layout,
        "--fields",
        fields,
        "--stride",
        str(args.stride),
        "--steps",
        args.steps,
        "--gamma",
        str(args.gamma),
        "--pvd-name",
        args.pvd_name,
    ]
    return output_dir, command


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the current generated case's rank-wise SLF output into "
            "merged ParaView VTI/PVD files."
        )
    )
    parser.add_argument("--case-directory", help="Case directory relative to the run root")
    parser.add_argument("--input-dir", help="SLF directory; default: <case>/output")
    parser.add_argument("--output-dir", help="VTI directory; default: <case>/paraview")
    parser.add_argument("--meta", help="meta.json path; default: <input-dir>/meta.json")
    parser.add_argument(
        "--steps",
        default="latest",
        help="all, latest, comma-separated steps, or inclusive start:stop:stride",
    )
    parser.add_argument(
        "--fields",
        help="Comma-separated fields; defaults to rho,u,v,w,p for NSE",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=2,
        help="Spatial preview stride. Default: 2; use 1 for full resolution",
    )
    parser.add_argument(
        "--derive",
        choices=["auto", "none", "nse", "gpe"],
        default="auto",
    )
    parser.add_argument(
        "--layout",
        choices=["auto", "global", "rank"],
        default="auto",
    )
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--pvd-name", default="collection.pvd")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.stride < 1:
        parser.error("--stride must be at least 1")
    if importlib.util.find_spec("numpy") is None:
        parser.error(
            "NumPy is required. Install it into this Python environment with "
            f"{sys.executable} -m pip install numpy"
        )

    try:
        output_dir, command = build_command(args)
    except (OSError, ValueError, json.JSONDecodeError, PostprocessError) as exc:
        parser.error(str(exc))

    print(f"[INFO] ParaView output: {output_dir}")
    print(f"[CMD] {_display_command(command)}")
    if args.dry_run:
        return 0
    completed = subprocess.run(command, check=False)
    return int(completed.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
