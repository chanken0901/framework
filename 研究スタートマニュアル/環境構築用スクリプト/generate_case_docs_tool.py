#!/usr/bin/env python3
"""
generate_case_docs_tool_v3.py

case.yaml から input.dat, README.md, meta.json を一括生成するスクリプト。
さらに、case_index.csv を用いて既存ケースとの重複チェックを行う。

基本方針
--------
- 人間が編集する正本: case.yaml
- ソルバーが読む: input.dat
- 人間が読む: README.md
- Python等の管理スクリプトが読む: meta.json
- 既存ケース一覧: cases/case_index.csv

重複チェック
------------
デフォルトでは、以下のキーが既存ケースと一致すると重複と判定する。

  physics_model
  mach_number
  reynolds_number
  reynolds_lambda
  turbulent_mach_number
  nx
  ny
  nz
  scheme
  reconstruction
  time_integration

同一条件が見つかった場合は、input.dat / README.md / meta.json を生成せず中止する。

使い方
------
通常:

    python scripts/generate_case_docs_tool_v3.py --case cases/case0001 --overwrite

重複チェックを無効化:

    python scripts/generate_case_docs_tool_v3.py --case cases/case0001 --overwrite --skip-duplicate-check

比較キーを指定:

    python scripts/generate_case_docs_tool_v3.py --case cases/case0001 --overwrite \
      --duplicate-keys physics_model mach_number nx ny nz scheme time_integration

個別生成:

    python scripts/generate_case_docs_tool_v3.py --case cases/case0001 --only input --overwrite

確認だけ:

    python scripts/generate_case_docs_tool_v3.py --case cases/case0001 --dry-run
"""

"""
generate_case_docs_tool_v4.py

case.yaml から input.dat, README.md, meta.json を一括生成するスクリプト。
case_index.csv を用いて既存ケースとの重複チェックも行う。

v4 の主な変更点
----------------
- case.yaml の新設計に対応:
    physics.model
    physics.nse.*
    flow.type
    solver.type
- input.dat 生成レイアウトを共通 / 物理モデル固有 / flow 固有に分離。
- NSE と TGV をまず正式対応。
- 将来 GPE を追加しやすい構造に変更。
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple


# ============================================================
# input.dat layout definitions
# ============================================================

COMMON_INPUT_LAYOUT = [
    ("case", [
        ("case_id", "case_id"),
        ("case_label", "case_label"),
        ("physics_model", "physics.model"),
        ("flow_type", "flow.type"),
        ("solver_type", "solver.type"),
    ]),
    ("grid", [
        ("nx", "grid.nx"),
        ("ny", "grid.ny"),
        ("nz", "grid.nz"),
        ("x_min", "grid.x_min"),
        ("x_max", "grid.x_max"),
        ("y_min", "grid.y_min"),
        ("y_max", "grid.y_max"),
        ("z_min", "grid.z_min"),
        ("z_max", "grid.z_max"),
    ]),
    ("time", [
        ("cfl", "time.cfl"),
        ("t_max", "time.t_max"),
        ("dt", "time.dt"),
        ("output_frequency", "time.output_frequency"),
        ("save_interval", "time.save_interval"),
    ]),
    ("numerics", [
        ("flux", "numerics.flux"),
        ("reconstruction", "numerics.reconstruction"),
        ("time_integration", "numerics.time_integration"),
    ]),
    ("storage", [
        ("output_location", "storage.output_location"),
        ("restart_location", "storage.restart_location"),
    ]),
]

NSE_INPUT_LAYOUT = [
    ("physics_nse", [
        ("gamma", "physics.nse.gamma"),
        ("rho_0", "physics.nse.rho_0"),
        ("mach_number", "physics.nse.mach_number"),
        ("reynolds_number", "physics.nse.reynolds_number"),
        ("prandtl_number", "physics.nse.prandtl_number"),
    ]),
]

GPE_INPUT_LAYOUT = [
    ("physics_gpe", [
        ("hbar", "physics.gpe.hbar"),
        ("mass", "physics.gpe.mass"),
        ("interaction_strength", "physics.gpe.interaction_strength"),
        ("density0", "physics.gpe.density0"),
        ("chemical_potential", "physics.gpe.chemical_potential"),
        ("damping", "physics.gpe.damping"),
    ]),
]

TGV_INPUT_LAYOUT = [
    ("flow_tgv", [
        ("tgv_amplitude", "flow.tgv.amplitude"),
        ("tgv_mode_x", "flow.tgv.mode_x"),
        ("tgv_mode_y", "flow.tgv.mode_y"),
        ("tgv_mode_z", "flow.tgv.mode_z"),
    ]),
]

HIT_INPUT_LAYOUT = [
    ("flow_hit", [
        ("hit_reynolds_lambda", "flow.hit.reynolds_lambda"),
        ("hit_turbulent_mach_number", "flow.hit.turbulent_mach_number"),
        ("hit_integral_length", "flow.hit.integral_length"),
        ("hit_random_seed", "flow.hit.random_seed"),
    ]),
]


# ============================================================
# case_index.csv mapping
# ============================================================

INDEX_FIELD_MAP = {
    "case_id": "case_id",
    "case_label": "case_label",
    "status": "status",
    "description": "description",

    "physics_model": "physics.model",
    "flow_type": "flow.type",
    "solver_type": "solver.type",

    "nx": "grid.nx",
    "ny": "grid.ny",
    "nz": "grid.nz",

    "mach_number": "physics.nse.mach_number",
    "reynolds_number": "physics.nse.reynolds_number",
    "prandtl_number": "physics.nse.prandtl_number",

    "flux": "numerics.flux",
    "reconstruction": "numerics.reconstruction",
    "time_integration": "numerics.time_integration",

    "raw_data_location": "storage.raw_data_location",

    # TGV固有
    "tgv_amplitude": "flow.tgv.amplitude",
    "tgv_mode_x": "flow.tgv.mode_x",
    "tgv_mode_y": "flow.tgv.mode_y",
    "tgv_mode_z": "flow.tgv.mode_z",

    # HIT固有：必要になったら使う
    "hit_reynolds_lambda": "flow.hit.reynolds_lambda",
    "hit_turbulent_mach_number": "flow.hit.turbulent_mach_number",
}

COMMON_DUPLICATE_KEYS = [
    "physics_model",
    "flow_type",
    "solver_type",
    "nx",
    "ny",
    "nz",
    "mach_number",
    "reynolds_number",
    "prandtl_number",
    "flux",
    "reconstruction",
    "time_integration",
]

DUPLICATE_KEYS_BY_FLOW = {
    "TGV": COMMON_DUPLICATE_KEYS + [
        "tgv_amplitude",
        "tgv_mode_x",
        "tgv_mode_y",
        "tgv_mode_z",
    ],
    "HIT": COMMON_DUPLICATE_KEYS + [
        "hit_reynolds_lambda",
        "hit_turbulent_mach_number",
    ],
}


# ============================================================
# small YAML parser
# ============================================================

def strip_comment(line: str) -> str:
    if "#" in line:
        return line.split("#", 1)[0]
    return line


def parse_scalar(text: str) -> Any:
    text = text.strip()
    if text == "":
        return ""
    if (text.startswith('"') and text.endswith('"')) or (text.startswith("'") and text.endswith("'")):
        return text[1:-1]

    lower = text.lower()
    if lower in {"true", "yes", "on"}:
        return True
    if lower in {"false", "no", "off"}:
        return False
    if lower in {"null", "none", "~"}:
        return ""

    try:
        if re.search(r"[.eE]", text):
            return float(text)
        return int(text)
    except ValueError:
        return text


def parse_simple_yaml(path: Path) -> Dict[str, Any]:
    root: Dict[str, Any] = {}
    stack: List[Tuple[int, Any]] = [(-1, root)]
    lines = path.read_text(encoding="utf-8").splitlines()

    for i, raw in enumerate(lines):
        line = strip_comment(raw).rstrip()
        if not line.strip():
            continue

        indent = len(line) - len(line.lstrip(" "))
        content = line.strip()

        while stack and indent <= stack[-1][0]:
            stack.pop()

        parent = stack[-1][1]

        if content.startswith("- "):
            if not isinstance(parent, list):
                raise ValueError(f"List item found but parent is not a list near line {i+1}")
            parent.append(parse_scalar(content[2:].strip()))
            continue

        if ":" not in content:
            raise ValueError(f"Invalid YAML line {i+1}: {raw}")

        key, value = content.split(":", 1)
        key = key.strip()
        value = value.strip()

        if not isinstance(parent, dict):
            raise ValueError(f"Parent is not dictionary near line {i+1}")

        if value == "":
            next_content = ""
            next_indent = -1
            for later in lines[i + 1:]:
                cleaned = strip_comment(later).rstrip()
                if not cleaned.strip():
                    continue
                next_indent = len(cleaned) - len(cleaned.lstrip(" "))
                next_content = cleaned.strip()
                break

            if next_content and next_indent > indent:
                child: Any = [] if next_content.startswith("- ") else {}
                parent[key] = child
                stack.append((indent, child))
            else:
                parent[key] = ""
        else:
            parent[key] = parse_scalar(value)

    return root


def get(data: Dict[str, Any], path: str, default: Any = "") -> Any:
    cur: Any = data
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def input_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        return ", ".join(input_value(v) for v in value)
    return str(value)


def normalize_for_compare(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple, set)) and len(value) == 0:
        return ""

    text = str(value).strip()
    if text == "":
        return ""

    try:
        num = float(text)
        return f"{num:.15g}"
    except ValueError:
        return text.lower()


# ============================================================
# git helpers
# ============================================================

def git_commit_hash() -> str:
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True, check=True)
        return r.stdout.strip()
    except Exception:
        return ""


def git_is_dirty() -> bool:
    try:
        r = subprocess.run(["git", "status", "--porcelain"], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True, check=True)
        return bool(r.stdout.strip())
    except Exception:
        return False


# ============================================================
# layout selector
# ============================================================

def normalized_name(value: Any) -> str:
    return str(value).strip().lower()


def build_input_layout(data: Dict[str, Any]) -> List[Tuple[str, List[Tuple[str, str]]]]:
    physics_model = normalized_name(get(data, "physics.model"))
    flow_type = normalized_name(get(data, "flow.type"))

    layout: List[Tuple[str, List[Tuple[str, str]]]] = []
    layout.extend(COMMON_INPUT_LAYOUT)

    if physics_model in {"nse", "navier_stokes", "compressible_nse", "compressible_euler"}:
        layout.extend(NSE_INPUT_LAYOUT)
    elif physics_model == "gpe":
        layout.extend(GPE_INPUT_LAYOUT)
    else:
        # 未対応モデルでも共通部分だけは出力し、後段で気付けるようにする。
        print(f"[WARNING] Unknown physics.model: {get(data, 'physics.model')}")

    if flow_type == "tgv":
        layout.extend(TGV_INPUT_LAYOUT)
    elif flow_type == "hit":
        layout.extend(HIT_INPUT_LAYOUT)

    return layout


# ============================================================
# case_index.csv helpers
# ============================================================

def case_index_path_from_case_dir(case_dir: Path) -> Path:
    return case_dir.parent / "case_index.csv"


def read_case_index(case_index: Path) -> List[Dict[str, str]]:
    if not case_index.exists():
        return []

    with case_index.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def yaml_to_index_values(data: Dict[str, Any]) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for index_key, yaml_path in INDEX_FIELD_MAP.items():
        values[index_key] = input_value(get(data, yaml_path, ""))
    return values


def update_current_case_in_index(data: Dict[str, Any], case_index: Path) -> None:
    rows = read_case_index(case_index)
    current = yaml_to_index_values(data)
    current_case_id = current.get("case_id", "").strip()

    if not current_case_id:
        raise ValueError("case_id is empty in case.yaml")
    if not case_index.exists():
        raise FileNotFoundError(f"case_index.csv not found: {case_index}")

    with case_index.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

    if not fieldnames:
        raise ValueError(f"case_index.csv has no header: {case_index}")

    now = datetime.now().isoformat(timespec="seconds")
    found = False
    updated_rows: List[Dict[str, str]] = []

    for row in rows:
        new_row = dict(row)
        if row.get("case_id", "").strip() == current_case_id:
            for key, value in current.items():
                if key in fieldnames:
                    new_row[key] = value
            if "updated_at" in fieldnames:
                new_row["updated_at"] = now
            found = True
        updated_rows.append(new_row)

    if not found:
        new_row = {name: "" for name in fieldnames}
        for key, value in current.items():
            if key in fieldnames:
                new_row[key] = value
        if "created_at" in fieldnames:
            new_row["created_at"] = now
        if "updated_at" in fieldnames:
            new_row["updated_at"] = now
        updated_rows.append(new_row)

    temporary = case_index.with_suffix(case_index.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in updated_rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
    temporary.replace(case_index)


def check_duplicate_case(data: Dict[str, Any], case_index: Path, duplicate_keys: List[str]) -> List[Dict[str, str]]:
    rows = read_case_index(case_index)
    if not rows:
        return []

    current = yaml_to_index_values(data)
    current_case_id = current.get("case_id", "")

    duplicates: List[Dict[str, str]] = []
    for row in rows:
        row_case_id = row.get("case_id", "")
        if current_case_id and row_case_id == current_case_id:
            continue

        is_same = True
        for key in duplicate_keys:
            current_value = normalize_for_compare(current.get(key, ""))
            row_value = normalize_for_compare(row.get(key, ""))
            if current_value != row_value:
                is_same = False
                break

        if is_same:
            duplicates.append(row)

    return duplicates


def print_duplicate_error(duplicates: List[Dict[str, str]], duplicate_keys: List[str], data: Dict[str, Any]) -> None:
    current = yaml_to_index_values(data)

    print("[ERROR] Duplicate case detected.")
    print("        The current case has the same condition as an existing case.")
    print("")
    print("Comparison keys:")
    for key in duplicate_keys:
        print(f"  - {key}: {current.get(key, '')}")
    print("")
    print("Matched existing case(s):")

    for row in duplicates:
        print(f"  - case_id: {row.get('case_id', '')}")
        print(f"    case_label: {row.get('case_label', '')}")
        print(f"    status: {row.get('status', '')}")
        print(f"    description: {row.get('description', '')}")
        print(f"    raw_data_location: {row.get('raw_data_location', '')}")
        print("")

    print("Generation was stopped.")
    print("")
    print("To allow this intentionally, use:")
    print("  --skip-duplicate-check")
    print("")
    print("Or change the comparison keys using:")
    print("  --duplicate-keys physics_model flow_type solver_type mach_number nx ny nz flux time_integration")


# ============================================================
# document generators
# ============================================================

def generate_input_dat(data: Dict[str, Any], source: Path) -> str:
    lines: List[str] = [
        "# input.dat",
        "# Automatically generated from case.yaml.",
        "# Do not edit this file directly.",
        f"# Source: {source}",
        f"# Generated at: {datetime.now().isoformat(timespec='seconds')}",
        "",
    ]

    for section, items in build_input_layout(data):
        section_lines = []
        for output_key, yaml_path in items:
            value = get(data, yaml_path, "")
            if value == "":
                continue
            section_lines.append(f"{output_key} = {input_value(value)}")

        if section_lines:
            lines += [
                "# =========================",
                f"# {section}",
                "# =========================",
            ]
            lines += section_lines
            lines.append("")

    return "\n".join(lines)


def generate_readme(data: Dict[str, Any]) -> str:
    case_id = get(data, "case_id")
    case_label = get(data, "case_label")
    status = get(data, "status")
    description = get(data, "description")

    used_in = get(data, "outputs.used_in", [])
    if isinstance(used_in, list) and used_in:
        used_text = "\n".join(f"- {v}" for v in used_in)
    elif used_in:
        used_text = f"- {used_in}"
    else:
        used_text = "- "

    physics_model = get(data, "physics.model")
    flow_type = get(data, "flow.type")
    solver_type = get(data, "solver.type")

    return f"""# {case_id}: {case_label}

## 1. Case overview

- Case ID: `{case_id}`
- Case label: `{case_label}`
- Status: `{status}`
- Description: {description}
- Project: {get(data, "project.name")}
- Author: {get(data, "project.author")}

## 2. Model / flow / solver

- Physics model: `{physics_model}`
- Flow type: `{flow_type}`
- Solver type: `{solver_type}`

## 3. NSE parameters

- Gamma: {get(data, "physics.nse.gamma")}
- Reference density: {get(data, "physics.nse.rho_0")}
- Mach number: {get(data, "physics.nse.mach_number")}
- Reynolds number: {get(data, "physics.nse.reynolds_number")}
- Prandtl number: {get(data, "physics.nse.prandtl_number")}

## 4. Flow parameters

- TGV amplitude: {get(data, "flow.tgv.amplitude")}
- TGV mode x: {get(data, "flow.tgv.mode_x")}
- TGV mode y: {get(data, "flow.tgv.mode_y")}
- TGV mode z: {get(data, "flow.tgv.mode_z")}
- HIT Reynolds lambda: {get(data, "flow.hit.reynolds_lambda")}
- HIT turbulent Mach number: {get(data, "flow.hit.turbulent_mach_number")}

## 5. Grid

- nx: {get(data, "grid.nx")}
- ny: {get(data, "grid.ny")}
- nz: {get(data, "grid.nz")}
- x range: {get(data, "grid.x_min")} to {get(data, "grid.x_max")}
- y range: {get(data, "grid.y_min")} to {get(data, "grid.y_max")}
- z range: {get(data, "grid.z_min")} to {get(data, "grid.z_max")}

## 6. Numerical method

- Flux: `{get(data, "numerics.flux")}`
- Reconstruction: `{get(data, "numerics.reconstruction")}`
- Time integration: `{get(data, "numerics.time_integration")}`

## 7. Time settings

- CFL: {get(data, "time.cfl")}
- dt: {get(data, "time.dt")}
- t_max: {get(data, "time.t_max")}
- output_frequency: {get(data, "time.output_frequency")}

## 8. Data location

```text
{get(data, "storage.raw_data_location")}
```

## 9. Used in

{used_text}

## 10. Notes

Detailed notes should be written in `notes.md`.

---

This README was automatically generated from `case.yaml`.
"""


def generate_meta_json(data: Dict[str, Any], source: Path) -> str:
    payload = {
        "schema_version": "2.0",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "source": str(source),
        "case": {
            "case_id": get(data, "case_id"),
            "case_label": get(data, "case_label"),
            "status": get(data, "status"),
            "description": get(data, "description"),
            "created_at": get(data, "created_at"),
        },
        "project": get(data, "project", {}),
        "physics": get(data, "physics", {}),
        "flow": get(data, "flow", {}),
        "solver": get(data, "solver", {}),
        "grid": get(data, "grid", {}),
        "time": get(data, "time", {}),
        "numerics": get(data, "numerics", {}),
        "storage": get(data, "storage", {}),
        "outputs": get(data, "outputs", {}),
        "git": {
            "commit_hash": git_commit_hash(),
            "dirty": git_is_dirty(),
        },
    }
    return json.dumps(payload, ensure_ascii=False, indent=2)


def write(path: Path, text: str, overwrite: bool, dry_run: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Use --overwrite.")

    if dry_run:
        print(f"\n[DRY-RUN] {path}")
        print("-" * 60)
        print(text)
        print("-" * 60)
        return

    path.write_text(text, encoding="utf-8")
    print(f"[OK] Wrote {path}")


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate input.dat, README.md, meta.json from case.yaml with duplicate check."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--case", help="Case directory, e.g. cases/case0001")
    group.add_argument("--yaml", help="Path to case.yaml")

    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--only", choices=["all", "input", "readme", "meta"], default="all")

    parser.add_argument(
        "--case-index",
        default=None,
        help="Path to case_index.csv. Default: inferred from case directory.",
    )
    parser.add_argument(
        "--skip-duplicate-check",
        action="store_true",
        help="Skip duplicate condition check.",
    )
    parser.add_argument(
        "--duplicate-keys",
        nargs="+",
        default=None,
        help="Keys used for duplicate check. Default: automatically selected by flow.type.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.case:
        case_dir = Path(args.case)
        yaml_path = case_dir / "case.yaml"
    else:
        yaml_path = Path(args.yaml)
        case_dir = yaml_path.parent

    if not yaml_path.exists():
        print(f"[ERROR] case.yaml not found: {yaml_path}")
        return 1

    try:
        data = parse_simple_yaml(yaml_path)
    except Exception as exc:
        print(f"[ERROR] Failed to parse case.yaml: {exc}")
        return 1

    case_index = Path(args.case_index) if args.case_index else case_index_path_from_case_dir(case_dir)

    if not args.skip_duplicate_check:
        if not case_index.exists():
            print(f"[ERROR] case_index.csv not found: {case_index}")
            print("        Duplicate safety check cannot be performed, so generation was stopped.")
            print("        Specify the correct file with --case-index, or intentionally use")
            print("        --skip-duplicate-check.")
            return 1

        print(f"[INFO] Duplicate check index: {case_index.resolve()}")

        invalid_keys = [k for k in duplicate_keys if k not in INDEX_FIELD_MAP]
        if invalid_keys:
            print("[ERROR] Invalid duplicate key(s):")
            for k in invalid_keys:
                print(f"  - {k}")
            print("")
            print("Allowed keys:")
            for k in sorted(INDEX_FIELD_MAP):
                print(f"  - {k}")
            return 1

        flow_type = str(get(data, "flow.type", "")).strip()

        if args.duplicate_keys:
            duplicate_keys = args.duplicate_keys
        else:
            duplicate_keys = DUPLICATE_KEYS_BY_FLOW.get(flow_type, COMMON_DUPLICATE_KEYS)

        duplicates = check_duplicate_case(data, case_index, duplicate_keys)
        if duplicates:
            print_duplicate_error(duplicates, duplicate_keys, data)
            return 1

    try:
        if args.only in {"all", "input"}:
            write(case_dir / "input.dat", generate_input_dat(data, yaml_path), args.overwrite, args.dry_run)
        if args.only in {"all", "readme"}:
            write(case_dir / "README.md", generate_readme(data), args.overwrite, args.dry_run)
        if args.only in {"all", "meta"}:
            write(case_dir / "meta.json", generate_meta_json(data, yaml_path), args.overwrite, args.dry_run)
    except FileExistsError as exc:
        print(f"[ERROR] {exc}")
        return 1

    if not args.dry_run:
        try:
            update_current_case_in_index(data, case_index)
            print(f"[OK] Updated {case_index}")
        except Exception as exc:
            print(f"[ERROR] Failed to update case_index.csv: {exc}")
            return 1

    print("[DONE] Case documents generated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
