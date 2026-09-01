#!/usr/bin/env python3
"""Maintain one flattened case index outside generated run environments."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from yaml_support import YamlFormatError, load_yaml


BASE_COLUMNS = (
    "case_key",
    "model",
    "case_id",
    "case_label",
    "status",
    "description",
    "profile",
    "processes",
    "environment",
    "case_path",
    "registered_at_utc",
    "updated_at_utc",
)
SUMMARY_PATHS = {
    "case_id": "case_id",
    "case_label": "case_label",
    "status": "status",
    "description": "description",
}
SUMMARY_FLAT_KEYS = {
    "case_id",
    "case_label",
    "status",
    "description",
    "physics.model",
    "solver.profile",
    "solver.mpi_processes",
}
INDEX_REPLACE_TIMEOUT_SECONDS = 2.0
INDEX_REPLACE_RETRY_SECONDS = 0.1


class GlobalCaseIndexError(RuntimeError):
    """Raised when the shared case index cannot be updated safely."""


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise GlobalCaseIndexError(f"{label} must be a mapping")
    return value


def _nested(data: dict[str, Any], path: str, default: Any = "") -> Any:
    current: Any = data
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def _csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return str(value)


def _flatten(value: Any, prefix: str, output: dict[str, str]) -> None:
    if isinstance(value, dict):
        if not value and prefix:
            output[prefix] = "{}"
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            _flatten(child, child_prefix, output)
        return
    output[prefix] = _csv_value(value)


def _relative_case_path(case_path: Path, index_path: Path) -> str:
    try:
        relative = case_path.resolve().relative_to(index_path.parent.resolve())
        return relative.as_posix()
    except ValueError:
        return str(case_path.resolve())


def row_from_case(
    case: dict[str, Any],
    *,
    model: str,
    profile: str,
    processes: int,
    environment: str,
    case_path: Path,
    index_path: Path,
    registered_at_utc: str = "",
) -> dict[str, str]:
    case_id = str(_nested(case, "case_id", "")).strip()
    if not case_id:
        raise GlobalCaseIndexError(f"case_id is missing: {case_path}")
    model_name = str(model or _nested(case, "physics.model", "")).strip().lower()
    if not model_name:
        raise GlobalCaseIndexError(f"physics model is missing: {case_path}")

    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    row = {
        "case_key": f"{model_name}:{case_id}",
        "model": model_name,
        "case_id": case_id,
        "case_label": _csv_value(_nested(case, SUMMARY_PATHS["case_label"])),
        "status": _csv_value(_nested(case, SUMMARY_PATHS["status"])),
        "description": _csv_value(_nested(case, SUMMARY_PATHS["description"])),
        "profile": str(profile or _nested(case, "solver.profile", "")),
        "processes": str(
            processes or _nested(case, "solver.mpi_processes", 1)
        ),
        "environment": environment,
        "case_path": _relative_case_path(case_path, index_path),
        "registered_at_utc": registered_at_utc or now,
        "updated_at_utc": now,
    }
    flattened: dict[str, str] = {}
    _flatten(case, "", flattened)
    for key in SUMMARY_FLAT_KEYS:
        flattened.pop(key, None)
    row.update(flattened)
    return row


def _read_index(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        return [], []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = list(reader.fieldnames or [])
        if not columns:
            raise GlobalCaseIndexError(f"case index has no header: {path}")
        return columns, list(reader)


def _columns_for(
    existing: list[str], rows: list[dict[str, str]]
) -> list[str]:
    columns = list(BASE_COLUMNS)
    columns.extend(column for column in existing if column not in columns)
    discovered = sorted(
        {
            key
            for row in rows
            for key in row
            if key not in columns
        }
    )
    columns.extend(discovered)
    return columns


def _write_index(
    path: Path,
    rows: list[dict[str, str]],
    existing_columns: list[str] | None = None,
    *,
    replace_timeout_seconds: float = INDEX_REPLACE_TIMEOUT_SECONDS,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = _columns_for(existing_columns or [], rows)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
        deadline = time.monotonic() + max(0.0, replace_timeout_seconds)
        while True:
            try:
                temporary.replace(path)
                break
            except PermissionError as exc:
                if time.monotonic() >= deadline:
                    raise GlobalCaseIndexError(
                        "cannot update the shared case index because Windows "
                        f"is preventing replacement of {path}. Close Excel, "
                        "an editor, or another program that has case_index.csv "
                        "open; the index will be synchronized on the next run"
                    ) from exc
                time.sleep(INDEX_REPLACE_RETRY_SECONDS)
    finally:
        if temporary.exists():
            temporary.unlink()


@contextmanager
def _index_lock(path: Path, timeout_seconds: float = 15.0) -> Iterator[None]:
    lock_path = path.with_name(f"{path.name}.lock")
    deadline = time.monotonic() + timeout_seconds
    descriptor: int | None = None
    while descriptor is None:
        try:
            descriptor = os.open(
                lock_path,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            )
            os.write(
                descriptor,
                f"pid={os.getpid()}\ncreated={time.time()}\n".encode("ascii"),
            )
        except FileExistsError:
            try:
                stale = time.time() - lock_path.stat().st_mtime > 300.0
            except FileNotFoundError:
                continue
            if stale:
                try:
                    lock_path.unlink()
                except FileNotFoundError:
                    pass
                continue
            if time.monotonic() >= deadline:
                raise GlobalCaseIndexError(
                    f"timed out waiting for case index lock: {lock_path}"
                )
            time.sleep(0.05)
    try:
        yield
    finally:
        os.close(descriptor)
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass


def sync_case_document(
    case: dict[str, Any],
    *,
    index_path: Path,
    model: str,
    profile: str,
    processes: int,
    environment: str,
    case_path: Path,
) -> None:
    index_path = index_path.resolve()
    index_path.parent.mkdir(parents=True, exist_ok=True)
    with _index_lock(index_path):
        existing_columns, rows = _read_index(index_path)
        case_id = str(_nested(case, "case_id", "")).strip()
        model_name = str(model or _nested(case, "physics.model", "")).strip().lower()
        case_key = f"{model_name}:{case_id}"
        registered = ""
        retained = []
        for row in rows:
            if row.get("case_key") == case_key:
                registered = row.get("registered_at_utc", "")
            else:
                retained.append(row)
        retained.append(
            row_from_case(
                case,
                model=model_name,
                profile=profile,
                processes=processes,
                environment=environment,
                case_path=case_path,
                index_path=index_path,
                registered_at_utc=registered,
            )
        )
        retained.sort(key=lambda row: (row.get("model", ""), row.get("case_id", "")))
        _write_index(index_path, retained, existing_columns)


def sync_environment_case(
    environment_root: Path,
    lock: dict[str, Any],
    index_path: Path,
) -> None:
    environment_root = environment_root.resolve()
    case_path = environment_root / str(lock["case_directory"]) / "case.yaml"
    if not case_path.is_file():
        raise GlobalCaseIndexError(f"case design not found: {case_path}")
    case = _mapping(load_yaml(case_path), "case YAML")
    try:
        mpi_processes = int(_nested(case, "solver.mpi_processes", 1))
    except (TypeError, ValueError) as exc:
        raise GlobalCaseIndexError(
            "case solver.mpi_processes must be an integer"
        ) from exc
    sync_case_document(
        case,
        index_path=index_path,
        model=str(lock.get("model", "")),
        profile=str(lock.get("profile", "")),
        processes=mpi_processes,
        environment=environment_root.name,
        case_path=case_path,
    )


def case_index_path(environment_root: Path, lock: dict[str, Any]) -> Path:
    value = str(lock.get("case_index_path") or "../case_index.csv")
    path = Path(value)
    if not path.is_absolute():
        path = environment_root / path
    return path.resolve()


def _environment_records(root: Path) -> list[tuple[Path, dict[str, Any]]]:
    records: list[tuple[Path, dict[str, Any]]] = []
    for environment in sorted(root.iterdir()):
        lock_path = environment / "environment.lock.json"
        if not environment.is_dir() or not lock_path.is_file():
            continue
        try:
            lock = json.loads(lock_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise GlobalCaseIndexError(f"invalid environment lock: {lock_path}") from exc
        if not isinstance(lock, dict) or lock.get("schema_version") != 1:
            raise GlobalCaseIndexError(f"invalid environment lock: {lock_path}")
        records.append((environment, lock))
    return records


def rebuild_case_index(root: Path, index_path: Path, dry_run: bool = False) -> int:
    root = root.resolve()
    index_path = index_path.resolve()
    records = _environment_records(root)
    if dry_run:
        print(f"[DRY-RUN] index: {index_path}")
        for environment, lock in records:
            print(
                f"  {lock.get('model', '')}:{lock.get('case_id', '')} "
                f"<- {environment}"
            )
        return len(records)

    rows: list[dict[str, str]] = []
    existing_registered: dict[str, str] = {}
    _, existing_rows = _read_index(index_path)
    for row in existing_rows:
        key = row.get("case_key", "")
        if key:
            existing_registered[key] = row.get("registered_at_utc", "")

    for environment, lock in records:
        case_path = environment / str(lock["case_directory"]) / "case.yaml"
        case = _mapping(load_yaml(case_path), f"case YAML {case_path}")
        case_id = str(_nested(case, "case_id", "")).strip()
        model = str(lock.get("model") or _nested(case, "physics.model", "")).lower()
        key = f"{model}:{case_id}"
        rows.append(
            row_from_case(
                case,
                model=model,
                profile=str(lock.get("profile", "")),
                processes=int(lock.get("processes", 1)),
                environment=environment.name,
                case_path=case_path,
                index_path=index_path,
                registered_at_utc=existing_registered.get(key, ""),
            )
        )
    rows.sort(key=lambda row: (row.get("model", ""), row.get("case_id", "")))
    with _index_lock(index_path):
        _write_index(index_path, rows)
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build one shared, flattened case_index.csv for all run environments."
    )
    parser.add_argument("--root", required=True, help="ResearchRuns root")
    parser.add_argument("--index", help="Output CSV; default: <root>/case_index.csv")
    parser.add_argument("--rebuild", action="store_true", help="Rebuild from all cases")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        root = Path(args.root).resolve()
        index = Path(args.index).resolve() if args.index else root / "case_index.csv"
        if not args.rebuild:
            raise GlobalCaseIndexError("--rebuild is required")
        count = rebuild_case_index(root, index, args.dry_run)
        if not args.dry_run:
            print(f"[OK] Global case index: {index}")
            print(f"     cases: {count}")
        return 0
    except (
        GlobalCaseIndexError,
        KeyError,
        OSError,
        ValueError,
        YamlFormatError,
    ) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
