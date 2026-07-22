#!/usr/bin/env python3
"""Create ``cases/caseXXXX/case.yaml`` from the STEP1-generated schema."""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

from case_registry import next_case_id, set_nested, sync_case_index
from yaml_support import YamlFormatError, dump_yaml, load_yaml


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a case from templates/case_schema.yaml")
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--template", default="templates/case_schema.yaml")
    parser.add_argument("--index-schema", default="templates/case_index_schema.yaml")
    parser.add_argument("--case-id")
    parser.add_argument("--label", default="baseline")
    parser.add_argument("--description", default="")
    parser.add_argument("--data-root")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def create_case(args: argparse.Namespace) -> Path:
    root = Path(args.project_root).resolve()
    template_path = (root / args.template).resolve()
    index_schema_path = (root / args.index_schema).resolve()
    cases_root = root / "cases"
    case_index = cases_root / "case_index.csv"
    case_id = args.case_id or next_case_id(cases_root, case_index)
    case_dir = cases_root / case_id
    case_path = case_dir / "case.yaml"
    if case_dir.exists() and not args.overwrite:
        raise FileExistsError(f"case directory already exists: {case_dir}")

    template = load_yaml(template_path)
    if not isinstance(template, dict):
        raise ValueError("case template must be a mapping")
    case = copy.deepcopy(template)
    set_nested(case, "case_id", case_id)
    set_nested(case, "case_label", args.label)
    if args.description:
        set_nested(case, "description", args.description)
    if args.data_root:
        project_name = str(case.get("project", {}).get("name", root.name))
        data_path = Path(args.data_root) / project_name / case_id
        set_nested(case, "storage.raw_data_location", data_path.as_posix())

    text = dump_yaml(case)
    if args.dry_run:
        print(text)
        return case_path
    case_dir.mkdir(parents=True, exist_ok=True)
    case_path.write_text(text, encoding="utf-8")
    (case_dir / "notes.md").write_text(f"# Notes for {case_id}\n", encoding="utf-8")
    index_schema = load_yaml(index_schema_path)
    if not isinstance(index_schema, dict):
        raise ValueError("case index schema must be a mapping")
    sync_case_index(case, index_schema, case_index)
    print(f"[OK] Created {case_path}")
    return case_path


def main(argv: list[str] | None = None) -> int:
    try:
        create_case(parse_args(argv))
    except (FileExistsError, OSError, ValueError, YamlFormatError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
