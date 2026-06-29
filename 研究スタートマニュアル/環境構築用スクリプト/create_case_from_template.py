#!/usr/bin/env python3
"""
create_case_from_template.py

テンプレートファイル templates/case_template.yaml を使って
新しいケースディレクトリ cases/caseXXXX/ を作成するスクリプト。

設計方針
--------
- スクリプト本体は基本的に編集しない。
- デフォルト値を変えたい場合は templates/case_template.yaml を編集する。
- 各ケースの条件を変えたい場合は生成された cases/caseXXXX/case.yaml を編集する。
- input.dat, README.md, meta.json は case.yaml から別スクリプトで自動生成する。

想定ディレクトリ
----------------
ProjectName/
├── scripts/
│   └── create_case_from_template.py
├── templates/
│   └── case_template.yaml
└── cases/
    ├── case_index.csv
    └── case0001/
        ├── case.yaml
        └── notes.md

使い方
------
通常:

    python scripts/create_case_from_template.py --label baseline

ケース番号を明示:

    python scripts/create_case_from_template.py --case-id case0010 --label test

テンプレートを指定:

    python scripts/create_case_from_template.py --template templates/case_template.yaml --label M05

実行せず確認:

    python scripts/create_case_from_template.py --label baseline --dry-run
"""

from __future__ import annotations

import argparse
import csv
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List


CASE_INDEX_COLUMNS = [
    "case_id",
    "case_label",
    "status",
    "description",
    "physics_model",
    "mach_number",
    "reynolds_number",
    "reynolds_lambda",
    "turbulent_mach_number",
    "nx",
    "ny",
    "nz",
    "scheme",
    "reconstruction",
    "time_integration",
    "raw_data_location",
    "used_in",
    "created_at",
    "updated_at",
    "git_commit_hash",
    "notes",
]


NOTES_TEMPLATE = """# Notes for {case_id}

## Purpose

Write the purpose of this case here.

## Changes from previous case

Write what was changed from the previous case.

## Results

Write important findings here.

## Problems / Warnings

Write problems, warnings, or failed attempts here.
"""


def ensure_case_index(case_index: Path) -> None:
    case_index.parent.mkdir(parents=True, exist_ok=True)
    if case_index.exists():
        return

    with case_index.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CASE_INDEX_COLUMNS)
        writer.writeheader()


def read_case_index(case_index: Path) -> List[Dict[str, str]]:
    ensure_case_index(case_index)
    with case_index.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_case_index(case_index: Path, rows: List[Dict[str, str]]) -> None:
    with case_index.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CASE_INDEX_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in CASE_INDEX_COLUMNS})


def next_case_id(cases_root: Path, case_index: Path) -> str:
    max_id = 0

    if cases_root.exists():
        for path in cases_root.iterdir():
            if path.is_dir():
                m = re.fullmatch(r"case(\d{4,})", path.name)
                if m:
                    max_id = max(max_id, int(m.group(1)))

    for row in read_case_index(case_index):
        cid = row.get("case_id", "")
        m = re.fullmatch(r"case(\d{4,})", cid)
        if m:
            max_id = max(max_id, int(m.group(1)))

    return f"case{max_id + 1:04d}"


def replace_placeholders(text: str, mapping: Dict[str, str]) -> str:
    for key, value in mapping.items():
        text = text.replace("{{" + key + "}}", value)
    return text


def parse_simple_value_from_yaml(text: str, dotted_key: str) -> str:
    """
    case_index.csv へ最低限の情報を転記するための簡易抽出関数。
    複雑なYAML構文には対応しない。
    """
    parts = dotted_key.split(".")
    if len(parts) == 1:
        pattern = re.compile(rf"^\s*{re.escape(parts[0])}\s*:\s*(.*?)\s*$")
        for line in text.splitlines():
            m = pattern.match(line)
            if m:
                return m.group(1).strip().strip('"').strip("'")
        return ""

    parent, child = parts[0], parts[1]
    in_parent = False
    parent_indent = None

    for line in text.splitlines():
        if not line.strip() or line.strip().startswith("#"):
            continue

        indent = len(line) - len(line.lstrip(" "))

        if re.match(rf"^\s*{re.escape(parent)}\s*:\s*$", line):
            in_parent = True
            parent_indent = indent
            continue

        if in_parent:
            if indent <= (parent_indent or 0):
                in_parent = False
                parent_indent = None
                continue

            m = re.match(rf"^\s*{re.escape(child)}\s*:\s*(.*?)\s*$", line)
            if m:
                return m.group(1).strip().strip('"').strip("'")

    return ""


def append_case_index(case_index: Path, case_yaml_text: str) -> None:
    rows = read_case_index(case_index)

    case_id = parse_simple_value_from_yaml(case_yaml_text, "case_id")
    case_label = parse_simple_value_from_yaml(case_yaml_text, "case_label")

    if not case_id:
        raise ValueError("case_id could not be read from case.yaml")

    for row in rows:
        if row.get("case_id") == case_id:
            raise ValueError(f"case_id already exists in case_index.csv: {case_id}")

    now = datetime.now().isoformat(timespec="seconds")

    row = {
        "case_id": case_id,
        "case_label": case_label,
        "status": parse_simple_value_from_yaml(case_yaml_text, "status"),
        "description": parse_simple_value_from_yaml(case_yaml_text, "description"),
        "physics_model": parse_simple_value_from_yaml(case_yaml_text, "physics.model"),
        "mach_number": parse_simple_value_from_yaml(case_yaml_text, "physics.mach_number"),
        "reynolds_number": parse_simple_value_from_yaml(case_yaml_text, "physics.reynolds_number"),
        "reynolds_lambda": parse_simple_value_from_yaml(case_yaml_text, "physics.reynolds_lambda"),
        "turbulent_mach_number": parse_simple_value_from_yaml(case_yaml_text, "physics.turbulent_mach_number"),
        "nx": parse_simple_value_from_yaml(case_yaml_text, "grid.nx"),
        "ny": parse_simple_value_from_yaml(case_yaml_text, "grid.ny"),
        "nz": parse_simple_value_from_yaml(case_yaml_text, "grid.nz"),
        "scheme": parse_simple_value_from_yaml(case_yaml_text, "numerics.scheme"),
        "reconstruction": parse_simple_value_from_yaml(case_yaml_text, "numerics.reconstruction"),
        "time_integration": parse_simple_value_from_yaml(case_yaml_text, "numerics.time_integration"),
        "raw_data_location": parse_simple_value_from_yaml(case_yaml_text, "storage.raw_data_location"),
        "created_at": now,
        "updated_at": now,
    }

    rows.append(row)
    write_case_index(case_index, rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a new case directory from templates/case_template.yaml."
    )

    parser.add_argument("--cases-root", default="cases")
    parser.add_argument("--case-index", default=None)
    parser.add_argument("--template", default="templates/case_template.yaml")
    parser.add_argument("--case-id", default=None)
    parser.add_argument("--label", default="baseline")
    parser.add_argument("--description", default="")
    parser.add_argument("--project-name", default=None)
    parser.add_argument("--data-root", default="/mnt/data_nas")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    project_root = Path.cwd()
    project_name = args.project_name or project_root.name

    cases_root = Path(args.cases_root)
    case_index = Path(args.case_index) if args.case_index else cases_root / "case_index.csv"
    template_path = Path(args.template)

    ensure_case_index(case_index)

    if not template_path.exists():
        print(f"[ERROR] Template not found: {template_path}")
        return 1

    case_id = args.case_id or next_case_id(cases_root, case_index)
    case_dir = cases_root / case_id
    raw_data_location = f"{args.data_root.rstrip('/')}/{project_name}/{case_id}"

    mapping = {
        "case_id": case_id,
        "case_label": args.label,
        "description": args.description,
        "project_name": project_name,
        "raw_data_location": raw_data_location,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }

    template_text = template_path.read_text(encoding="utf-8")
    case_yaml_text = replace_placeholders(template_text, mapping)

    print("[INFO] Create case from template")
    print(f"  case_id:           {case_id}")
    print(f"  case_label:        {args.label}")
    print(f"  template:          {template_path}")
    print(f"  case_dir:          {case_dir}")
    print(f"  raw_data_location: {raw_data_location}")

    if args.dry_run:
        print("[DRY-RUN] No files were created.")
        print("")
        print(case_yaml_text)
        return 0

    if case_dir.exists() and not args.overwrite:
        print(f"[ERROR] Case directory already exists: {case_dir}")
        print("        Use --overwrite if intentional.")
        return 1

    case_dir.mkdir(parents=True, exist_ok=True)

    case_yaml_path = case_dir / "case.yaml"
    notes_path = case_dir / "notes.md"

    if case_yaml_path.exists() and not args.overwrite:
        print(f"[ERROR] case.yaml already exists: {case_yaml_path}")
        return 1

    case_yaml_path.write_text(case_yaml_text, encoding="utf-8")

    if not notes_path.exists() or args.overwrite:
        notes_path.write_text(NOTES_TEMPLATE.format(case_id=case_id), encoding="utf-8")

    append_case_index(case_index, case_yaml_text)

    print("[OK] Case created.")
    print(f"  case.yaml: {case_yaml_path}")
    print(f"  notes.md:  {notes_path}")
    print(f"  index:     {case_index}")
    print("")
    print("Next steps:")
    print(f"  1. Edit {case_yaml_path}")
    print("  2. Generate input.dat from case.yaml")
    print(f"     python scripts/generate_case_docs_tool.py --case {case_dir} --overwrite")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
