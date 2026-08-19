#!/usr/bin/env python3
"""Schema-driven case registry helpers shared by case-generation tools."""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any


def nested(data: dict[str, Any], path: str, default: Any = "") -> Any:
    current: Any = data
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def set_nested(data: dict[str, Any], path: str, value: Any) -> None:
    parts = path.split(".")
    current = data
    for key in parts[:-1]:
        child = current.get(key)
        if not isinstance(child, dict):
            child = {}
            current[key] = child
        current = child
    current[parts[-1]] = value


def csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return str(value)


def normalized(value: Any) -> str:
    text = csv_value(value).strip()
    if not text:
        return ""
    try:
        return f"{float(text):.15g}"
    except ValueError:
        return text.casefold()


def next_case_id(cases_root: Path, case_index: Path) -> str:
    maximum = 0
    if cases_root.exists():
        for path in cases_root.iterdir():
            match = re.fullmatch(r"case(\d{4,})", path.name) if path.is_dir() else None
            if match:
                maximum = max(maximum, int(match.group(1)))
    if case_index.is_file():
        with case_index.open("r", newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                match = re.fullmatch(r"case(\d{4,})", row.get("case_id", ""))
                if match:
                    maximum = max(maximum, int(match.group(1)))
    return f"case{maximum + 1:04d}"


def row_from_case(case: dict[str, Any], schema: dict[str, Any]) -> dict[str, str]:
    paths = schema.get("column_paths", {}) or {}
    row = {str(column): csv_value(nested(case, str(path))) for column, path in paths.items()}
    now = datetime.now().isoformat(timespec="seconds")
    row.setdefault("created_at", now)
    row["updated_at"] = now
    return row


def find_duplicates(
    case: dict[str, Any], schema: dict[str, Any], case_index: Path
) -> list[dict[str, str]]:
    if not case_index.is_file():
        return []
    current = row_from_case(case, schema)
    current_id = current.get("case_id", "")
    keys = [str(key) for key in schema.get("duplicate_check_keys", []) or []]
    if not keys:
        return []
    with case_index.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return [
        row
        for row in rows
        if row.get("case_id", "") != current_id
        and all(normalized(row.get(key, "")) == normalized(current.get(key, "")) for key in keys)
    ]


def sync_case_index(
    case: dict[str, Any], schema: dict[str, Any], case_index: Path
) -> None:
    columns = [str(column) for column in schema.get("columns", []) or []]
    if not columns:
        raise ValueError("case index schema has no columns")
    case_index.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    existing_columns: list[str] = []
    if case_index.is_file():
        with case_index.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            existing_columns = list(reader.fieldnames or [])
            rows = list(reader)

    options = schema.get("options", {}) or {}
    if options.get("keep_unknown_existing_columns", True):
        columns.extend(column for column in existing_columns if column not in columns)
    if not options.get("auto_add_missing_columns", True) and existing_columns:
        columns = existing_columns

    current = row_from_case(case, schema)
    case_id = current.get("case_id", "")
    if not case_id:
        raise ValueError("case_id is required")
    found = False
    for row in rows:
        if row.get("case_id", "") == case_id:
            created = row.get("created_at", "")
            row.update(current)
            if created:
                row["created_at"] = created
            found = True
    if not found:
        rows.append(current)

    temporary = case_index.with_suffix(case_index.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    temporary.replace(case_index)
