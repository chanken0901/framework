#!/usr/bin/env python3
"""Run ParaView conversion or NSE turbulence statistics for a generated case."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shlex
import subprocess
import sys
from pathlib import Path

from case_configuration import resolve_case_configuration
from case_input import CaseInputError, derive_nse_hit_transport
from yaml_support import YamlFormatError, load_yaml


TOOL_VERSION = "2.0.0"


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


def _statistics_path(root: Path, model: str) -> Path:
    if model != "nse":
        raise PostprocessError("Turbulence statistics are currently available only for NSE")
    path = (
        root
        / "SolverLibrary"
        / "NSE"
        / "tools"
        / "nse_turbulence_statistics.py"
    )
    if path.is_file():
        return path
    raise PostprocessError(
        "NSE turbulence-statistics tool is not available. "
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


def _nested(mapping: dict, dotted_path: str, default=None):
    value = mapping
    for part in dotted_path.split("."):
        if not isinstance(value, dict) or part not in value:
            return default
        value = value[part]
    return value


def _resolved_nse_reynolds(case_config: dict) -> float | None:
    flow_type = str(_nested(case_config, "flow.type", "")).strip().lower()
    flow_type = "_".join(flow_type.replace("-", " ").split())
    if flow_type in {
        "hit",
        "hit_spectral",
        "homogeneous_isotropic_turbulence",
    }:
        hit = _nested(case_config, "flow.hit", {})
        if isinstance(hit, dict):
            mach = hit.get("turbulent_mach_number")
            re_lambda = hit.get(
                "turbulent_reynolds_number",
                hit.get("taylor_reynolds_number"),
            )
            spectrum = hit.get("spectrum")
            length = hit.get("integral_length")
            if isinstance(spectrum, dict):
                spectrum_type = str(spectrum.get("type", "")).strip().lower()
                spectrum_type = "_".join(
                    spectrum_type.replace("-", " ").split()
                )
                selected = spectrum.get(spectrum_type, {})
                if isinstance(selected, dict):
                    if spectrum_type in {"johnsen", "k4_gaussian"}:
                        length = selected.get("characteristic_length")
                    elif spectrum_type == "pope":
                        length = selected.get("integral_length")
            if mach not in {None, ""} and re_lambda not in {None, ""}:
                if length in {None, ""}:
                    raise PostprocessError(
                        "HIT Reynolds resolution requires a spectrum length"
                    )
                try:
                    return derive_nse_hit_transport(
                        mach, re_lambda, length
                    )["solver_reynolds"]
                except CaseInputError as exc:
                    raise PostprocessError(str(exc)) from exc

    value = _nested(case_config, "physics.nse.reynolds_number")
    return None if value is None else float(value)


def _case_context(args: argparse.Namespace) -> dict:
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
    meta = _resolve(root, args.meta, input_dir / "meta.json")
    case_yaml = case_root / "case.yaml"
    if not case_yaml.is_file():
        raise PostprocessError(f"case.yaml was not found: {case_yaml}")
    case_config = resolve_case_configuration(case_yaml).document
    gamma = args.gamma
    if gamma is None:
        gamma = float(_nested(case_config, "physics.nse.gamma", 1.4))
    reynolds = args.reynolds
    if reynolds is None and model == "nse":
        reynolds = _resolved_nse_reynolds(case_config)
    return {
        "root": root,
        "model": model,
        "case_root": case_root,
        "input_dir": input_dir,
        "meta": meta,
        "gamma": gamma,
        "reynolds": reynolds,
    }


def build_command(args: argparse.Namespace) -> tuple[Path, list[str]]:
    context = _case_context(args)
    root = context["root"]
    model = context["model"]
    output_dir = _resolve(root, args.output_dir, context["case_root"] / "paraview")
    converter = _converter_path(root, model)
    fields = args.fields if args.fields is not None else _default_fields(model)
    command = [
        sys.executable,
        str(converter),
        str(context["input_dir"]),
        "--output-dir",
        str(output_dir),
        "--meta",
        str(context["meta"]),
        "--derive",
        args.derive,
        "--layout",
        args.layout,
        "--fields",
        fields,
        "--stride",
        str(args.stride),
        "--steps",
        args.steps or "latest",
        "--gamma",
        str(context["gamma"]),
        "--pvd-name",
        args.pvd_name,
    ]
    if getattr(args, "inspect_only", False):
        command.append("--inspect-only")
    return output_dir, command


def build_statistics_command(args: argparse.Namespace) -> tuple[Path, list[str]]:
    context = _case_context(args)
    root = context["root"]
    model = context["model"]
    statistics = _statistics_path(root, model)
    reynolds = context["reynolds"]
    if reynolds is None or reynolds <= 0.0:
        raise PostprocessError(
            "A positive Reynolds number is required. Set "
            "flow.hit turbulent targets, physics.nse.reynolds_number, "
            "or use --reynolds."
        )
    output = _resolve(
        root,
        args.statistics_output,
        context["case_root"] / "statistics" / "turbulence_statistics.csv",
    )
    command = [
        sys.executable,
        str(statistics),
        str(context["input_dir"]),
        "--output",
        str(output),
        "--meta",
        str(context["meta"]),
        "--layout",
        args.layout,
        "--steps",
        args.steps or "all",
        "--gamma",
        str(context["gamma"]),
        "--reynolds",
        str(reynolds),
        "--density-floor",
        str(args.density_floor),
        "--pressure-floor",
        str(args.pressure_floor),
    ]
    return output, command


def build_commands(args: argparse.Namespace) -> list[tuple[str, Path, list[str]]]:
    commands: list[tuple[str, Path, list[str]]] = []
    if args.task in {"paraview", "all"}:
        output, command = build_command(args)
        commands.append(("ParaView", output, command))
    if args.task in {"statistics", "all"}:
        output, command = build_statistics_command(args)
        commands.append(("Turbulence statistics", output, command))
    return commands


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run ParaView conversion and NSE turbulence-statistics "
            "postprocessing for the current generated case."
        )
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {TOOL_VERSION}")
    parser.add_argument(
        "--task",
        choices=["paraview", "statistics", "all"],
        default="paraview",
        help="Postprocessing task. Default: paraview",
    )
    parser.add_argument("--case-directory", help="Case directory relative to the run root")
    parser.add_argument("--input-dir", help="SLF directory; default: <case>/output")
    parser.add_argument("--output-dir", help="VTI directory; default: <case>/paraview")
    parser.add_argument(
        "--statistics-output",
        help="Statistics CSV path; default: <case>/statistics/turbulence_statistics.csv",
    )
    parser.add_argument("--meta", help="meta.json path; default: <input-dir>/meta.json")
    parser.add_argument(
        "--steps",
        default=None,
        help=(
            "all, latest, comma-separated steps, or inclusive start:stop:stride. "
            "Defaults: latest for ParaView and all for statistics"
        ),
    )
    parser.add_argument(
        "--fields",
        help=(
            "Comma-separated fields. Defaults: rho,u,v,w,p for NSE and "
            "density,phase for GPE"
        ),
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
    parser.add_argument(
        "--gamma", type=float, help="Override physics.nse.gamma from case.yaml"
    )
    parser.add_argument(
        "--reynolds",
        type=float,
        help="Override physics.nse.reynolds_number for turbulence statistics",
    )
    parser.add_argument("--density-floor", type=float, default=1.0e-12)
    parser.add_argument("--pressure-floor", type=float, default=1.0e-12)
    parser.add_argument("--pvd-name", default="collection.pvd")
    parser.add_argument(
        "--inspect-only",
        action="store_true",
        help=(
            "Read SLF headers/meta.json and show the selected steps without "
            "writing VTI/PVD files (ParaView task only)"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show the child command only; do not read SLF files or write output",
    )
    args = parser.parse_args(argv)

    if args.stride < 1:
        parser.error("--stride must be at least 1")
    if args.inspect_only and args.task != "paraview":
        parser.error("--inspect-only can be used only with --task paraview")
    if importlib.util.find_spec("numpy") is None:
        parser.error(
            "NumPy is required. Install it into this Python environment with "
            f"{sys.executable} -m pip install numpy"
        )

    try:
        commands = build_commands(args)
    except (
        OSError,
        ValueError,
        json.JSONDecodeError,
        YamlFormatError,
        PostprocessError,
    ) as exc:
        parser.error(str(exc))

    print(f"[INFO] postprocess_case version: {TOOL_VERSION}")
    for label, output, command in commands:
        print(f"[INFO] {label} output: {output}")
        print(f"[CMD] {_display_command(command)}")
        if args.dry_run:
            continue
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            return int(completed.returncode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
