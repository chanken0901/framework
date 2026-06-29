#!/usr/bin/env python3
"""
init_research_project_schema_driven.py

STEP1: project_schema.yaml を Single Source of Truth として読み込み，
研究プロジェクトの標準ディレクトリ構造と，各種テンプレートを生成するスクリプト。

生成される主なファイル:
  - README.md
  - .gitignore
  - Makefile
  - cases/case_index.csv
  - templates/case_schema.yaml
  - templates/case_index_schema.yaml
  - templates/input_schema.yaml
  - templates/meta_schema.yaml
  - templates/readme_schema.yaml
  - config/machine.yaml
  - config/compiler.yaml
  - config/paths.yaml

使い方:
  新しいプロジェクトを作る場合:
      python init_research_project_schema_driven.py --schema project_schema.yaml --project MyProject

  現在のディレクトリに構造を作る場合:
      python init_research_project_schema_driven.py --schema project_schema.yaml --here

  既存テンプレート等を上書きしたい場合:
      python init_research_project_schema_driven.py --schema project_schema.yaml --project MyProject --overwrite

設計方針:
  - CASE_INDEX_COLUMNS や CASE_YAML_TEMPLATE を Python 内に固定しない。
  - ディレクトリ，case_index.csv，case_schema.yaml などは project_schema.yaml から生成する。
  - STEP2以降のツールは project_schema.yaml ではなく templates/*.yaml を読む想定にする。
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyYAML が必要です。次を実行してください: python -m pip install pyyaml"
    ) from exc


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Schema root must be a mapping: {path}")
    return data


def write_text_if_needed(path: Path, text: str, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_yaml_if_needed(path: Path, data: dict[str, Any], overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, allow_unicode=True, sort_keys=False)


def flatten_directory_groups(directories: dict[str, Any]) -> list[str]:
    paths: list[str] = []
    for _, group_paths in directories.items():
        if isinstance(group_paths, list):
            paths.extend(str(p) for p in group_paths)
    return paths


def enabled_items(schema: dict[str, Any], section: str, collection_key: str) -> dict[str, Any]:
    sec = schema.get(section, {}) or {}
    enabled = sec.get("enabled", []) or []
    all_items = sec.get(collection_key, {}) or {}
    return {name: all_items[name] for name in enabled if name in all_items}


def parameter_defaults(parameters: dict[str, Any], output_key: bool = False) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, spec in (parameters or {}).items():
        if not isinstance(spec, dict):
            out[key] = spec
            continue
        name = spec.get("output_key", key) if output_key else key
        out[name] = spec.get("default")
    return out


def parameter_descriptions(parameters: dict[str, Any], output_key: bool = False) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, spec in (parameters or {}).items():
        if not isinstance(spec, dict):
            out[key] = {"default": spec, "description": ""}
            continue
        name = spec.get("output_key", key) if output_key else key
        out[name] = {
            "default": spec.get("default"),
            "description": spec.get("description", ""),
        }
    return out


def unique(seq: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


# -----------------------------------------------------------------------------
# Directory / top-level files
# -----------------------------------------------------------------------------


def make_dirs(root: Path, schema: dict[str, Any]) -> None:
    for rel in flatten_directory_groups(schema.get("directories", {}) or {}):
        p = root / rel
        p.mkdir(parents=True, exist_ok=True)
        gitkeep = p / ".gitkeep"
        if not gitkeep.exists():
            gitkeep.write_text("", encoding="utf-8")


def create_readme(root: Path, schema: dict[str, Any], overwrite: bool) -> None:
    project = schema.get("project", {}) or {}
    name = project.get("name", "ResearchProject")
    description = project.get("description", "Numerical simulation research framework")

    directory_lines = []
    for group, paths in (schema.get("directories", {}) or {}).items():
        if isinstance(paths, list):
            directory_lines.append(f"### {group}")
            directory_lines.extend(f"- `{p}/`" for p in paths)
            directory_lines.append("")

    text = f"""# {name}

{description}

このプロジェクトは `project_schema.yaml` を Single Source of Truth として初期化されています。

## 基本方針

- ディレクトリ構造，ケース管理表，テンプレートは `project_schema.yaml` から生成する。
- STEP2以降のツールは `templates/*.yaml` を参照してケース生成・重複確認・入力ファイル生成を行う。
- ソースコード，解析コード，論文，計算条件は Git/GitHub で管理する。
- VTK, HDF5, restart, checkpoint などの巨大データは Git 管理しない。
- 巨大データの保存場所は `case.yaml`, `input.dat`, `meta.json`, `case_index.csv` に記録する。

## 生成された主要ディレクトリ

{chr(10).join(directory_lines).rstrip()}

## 初期化後の推奨操作

```bash
git init
git add .
git commit -m "Initial schema-driven project structure"
```
"""
    write_text_if_needed(root / "README.md", text, overwrite)


def create_gitignore(root: Path, overwrite: bool) -> None:
    text = """# Build products
*.o
*.mod
*.smod
*.a
*.so
*.dll
*.dylib
*.exe
*.out
bin/
build/

# Python
__pycache__/
*.pyc
.venv/
venv/
.env

# LaTeX temporary files
*.aux
*.bbl
*.bcf
*.blg
*.fdb_latexmk
*.fls
*.log
*.run.xml
*.synctex.gz
*.toc

# Large simulation data
*.h5
*.hdf5
*.nc
*.vtk
*.vtu
*.plt
*.raw
*.bin
*.dat.bin

# Case outputs
cases/*/run_env/
cases/*/output/
cases/*/logs/
cases/*/restart/
cases/*/checkpoint/
cases/*/tmp/

# Local machine settings
config/local.yaml
config/machine.local.yaml
config/machine.local.mk
config/machine.local.sh

# OS/editor files
.DS_Store
Thumbs.db
.vscode/
.idea/
"""
    write_text_if_needed(root / ".gitignore", text, overwrite)


def create_makefile(root: Path, overwrite: bool) -> None:
    text = """# Project-level Makefile

.PHONY: help case-list clean

help:
	@echo "Available targets:"
	@echo "  make case-list   Show registered cases"
	@echo "  make clean       Remove common temporary files"

case-list:
	@if [ -f cases/case_index.csv ]; then \
		cat cases/case_index.csv; \
	else \
		echo "cases/case_index.csv not found"; \
	fi

clean:
	find . -name "__pycache__" -type d -prune -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
"""
    write_text_if_needed(root / "Makefile", text, overwrite)


# -----------------------------------------------------------------------------
# Schema-driven generation
# -----------------------------------------------------------------------------


def build_case_index_columns(schema: dict[str, Any]) -> list[str]:
    columns = schema.get("case_index", {}).get("columns", {}) or {}
    physics_enabled = schema.get("physics", {}).get("enabled", []) or []
    flow_enabled = schema.get("flow", {}).get("enabled", []) or []

    result: list[str] = []

    for key in ["case", "model"]:
        result.extend(columns.get(key, []) or [])

    physics_columns = columns.get("physics", {}) or {}
    for model in physics_enabled:
        result.extend(physics_columns.get(model, []) or [])

    flow_columns = columns.get("flow", {}) or {}
    for flow in flow_enabled:
        result.extend(flow_columns.get(flow, []) or [])

    for key in ["grid", "numerics", "storage", "metadata"]:
        result.extend(columns.get(key, []) or [])

    return unique(str(c) for c in result)


def create_case_index(root: Path, schema: dict[str, Any], overwrite: bool) -> None:
    options = schema.get("case_index", {}).get("options", {}) or {}
    if not options.get("auto_create_case_index", True):
        return

    path = root / "cases" / "case_index.csv"
    if path.exists() and not overwrite:
        return

    columns = build_case_index_columns(schema)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()


def build_case_schema(schema: dict[str, Any]) -> dict[str, Any]:
    project = schema.get("project", {}) or {}
    physics_enabled = schema.get("physics", {}).get("enabled", []) or []
    flow_enabled = schema.get("flow", {}).get("enabled", []) or []
    solver_enabled = schema.get("solver", {}).get("enabled", []) or []

    physics_models = enabled_items(schema, "physics", "models")
    flow_types = enabled_items(schema, "flow", "types")
    solver_types = enabled_items(schema, "solver", "types")

    default_physics = physics_enabled[0] if physics_enabled else None
    default_flow = flow_enabled[0] if flow_enabled else None
    default_solver = solver_enabled[0] if solver_enabled else None

    physics_block: dict[str, Any] = {"model": default_physics}
    for name, spec in physics_models.items():
        case_key = spec.get("case_key", name.lower())
        physics_block[case_key] = parameter_defaults(spec.get("parameters", {}) or {})

    flow_block: dict[str, Any] = {"type": default_flow}
    for name, spec in flow_types.items():
        case_key = spec.get("case_key", name.lower())
        flow_block[case_key] = parameter_defaults(spec.get("parameters", {}) or {}, output_key=False)

    return {
        "case_id": "case0001",
        "case_label": project.get("default_case_label", "baseline"),
        "status": "planned",
        "description": "baseline case",
        "project": {
            "name": project.get("name", "ResearchProject"),
            "author": project.get("author", ""),
        },
        "physics": physics_block,
        "flow": flow_block,
        "solver": {"type": default_solver, "available": list(solver_types.keys())},
        "grid": parameter_defaults(schema.get("grid", {}).get("parameters", {}) or {}),
        "time": parameter_defaults(schema.get("time", {}).get("parameters", {}) or {}),
        "numerics": parameter_defaults(schema.get("numerics", {}).get("parameters", {}) or {}),
        "storage": parameter_defaults(schema.get("storage", {}).get("parameters", {}) or {}),
        "outputs": {"used_in": []},
        "notes": "",
    }


def build_case_index_schema(schema: dict[str, Any]) -> dict[str, Any]:
    physics_models = enabled_items(schema, "physics", "models")
    flow_types = enabled_items(schema, "flow", "types")

    duplicate_common = schema.get("case_index", {}).get("duplicate_check", {}).get("common", []) or []
    duplicate_keys = list(duplicate_common)

    for _, spec in physics_models.items():
        duplicate_keys.extend(spec.get("duplicate_keys", []) or [])
    for _, spec in flow_types.items():
        duplicate_keys.extend(spec.get("duplicate_keys", []) or [])
    duplicate_keys.extend(schema.get("grid", {}).get("duplicate_keys", []) or [])
    duplicate_keys.extend(schema.get("numerics", {}).get("duplicate_keys", []) or [])

    return {
        "columns": build_case_index_columns(schema),
        "duplicate_check_keys": unique(str(k) for k in duplicate_keys),
        "options": schema.get("case_index", {}).get("options", {}) or {},
    }


def build_input_schema(schema: dict[str, Any]) -> dict[str, Any]:
    physics_models = enabled_items(schema, "physics", "models")
    flow_types = enabled_items(schema, "flow", "types")

    sections: dict[str, Any] = {
        "case": ["case_id", "case_label", "status", "description"],
        "model": ["physics_model", "flow_type", "solver_type"],
        "grid": list((schema.get("grid", {}).get("parameters", {}) or {}).keys()),
        "time": list((schema.get("time", {}).get("parameters", {}) or {}).keys()),
        "numerics": list((schema.get("numerics", {}).get("parameters", {}) or {}).keys()),
        "storage": list((schema.get("storage", {}).get("parameters", {}) or {}).keys()),
    }

    for name, spec in physics_models.items():
        sections[spec.get("input_section", f"physics_{name.lower()}")] = list(
            (spec.get("parameters", {}) or {}).keys()
        )
    for name, spec in flow_types.items():
        sections[spec.get("input_section", f"flow_{name.lower()}")] = list(
            (spec.get("parameters", {}) or {}).keys()
        )

    return {
        "format": "key_value",
        "extension": ".dat",
        "comment_prefix": "#",
        "sections": sections,
        "note": "generate_case_docs_tool.py が case.yaml から input.dat を生成するときに参照する。",
    }


def build_meta_schema(schema: dict[str, Any]) -> dict[str, Any]:
    return {
        "case": ["case_id", "case_label", "status", "description"],
        "project": ["name", "author"],
        "model": ["physics_model", "flow_type", "solver_type"],
        "physics": list(enabled_items(schema, "physics", "models").keys()),
        "flow": list(enabled_items(schema, "flow", "types").keys()),
        "grid": list((schema.get("grid", {}).get("parameters", {}) or {}).keys()),
        "time": list((schema.get("time", {}).get("parameters", {}) or {}).keys()),
        "numerics": list((schema.get("numerics", {}).get("parameters", {}) or {}).keys()),
        "storage": list((schema.get("storage", {}).get("parameters", {}) or {}).keys()),
        "git": {
            "record_commit_hash": schema.get("git", {}).get("record_commit_hash", True),
            "warn_if_dirty": schema.get("git", {}).get("warn_if_dirty", True),
        },
    }


def build_readme_schema(schema: dict[str, Any]) -> dict[str, Any]:
    return {
        "title": "{case_id}",
        "sections": [
            "Case Information",
            "Physical model",
            "Flow/problem setting",
            "Grid",
            "Numerical methods",
            "Storage",
            "Outputs",
            "Notes",
        ],
        "fields": {
            "case": ["case_id", "case_label", "status", "description"],
            "model": ["physics_model", "flow_type", "solver_type"],
            "grid": ["nx", "ny", "nz"],
            "storage": ["raw_data_location", "output_location", "restart_location"],
        },
    }


def build_parameter_catalog(schema: dict[str, Any]) -> dict[str, Any]:
    catalog: dict[str, Any] = {}

    for name, spec in enabled_items(schema, "physics", "models").items():
        catalog[f"physics.{name}"] = parameter_descriptions(spec.get("parameters", {}) or {})

    for name, spec in enabled_items(schema, "flow", "types").items():
        catalog[f"flow.{name}"] = parameter_descriptions(
            spec.get("parameters", {}) or {}, output_key=True
        )

    for section in ["grid", "time", "numerics", "storage"]:
        catalog[section] = parameter_descriptions(
            schema.get(section, {}).get("parameters", {}) or {}
        )

    return catalog


def create_templates(root: Path, schema: dict[str, Any], overwrite: bool) -> None:
    template_dir = root / (schema.get("templates", {}).get("output_directory", "templates"))
    template_dir.mkdir(parents=True, exist_ok=True)

    generated = schema.get("templates", {}).get("generate", []) or []
    builders = {
        "case_schema.yaml": build_case_schema,
        "case_index_schema.yaml": build_case_index_schema,
        "input_schema.yaml": build_input_schema,
        "meta_schema.yaml": build_meta_schema,
        "readme_schema.yaml": build_readme_schema,
    }

    for filename in generated:
        builder = builders.get(filename)
        if builder is not None:
            write_yaml_if_needed(template_dir / filename, builder(schema), overwrite)

    write_yaml_if_needed(template_dir / "parameter_catalog.yaml", build_parameter_catalog(schema), overwrite)


def create_config(root: Path, schema: dict[str, Any], overwrite: bool) -> None:
    project_name = (schema.get("project", {}) or {}).get("name", "ResearchProject")
    storage_defaults = parameter_defaults(schema.get("storage", {}).get("parameters", {}) or {})

    machine = {
        "machine_name": "local_linux_server",
        "paths": {
            "repo_root": f"/mnt/repo_nas/Kento/ResearchRepo/{project_name}",
            "data_root": str(storage_defaults.get("raw_data_location", f"/mnt/data_nas/{project_name}"))
            .replace("/case0001", ""),
            "scratch_root": f"/scratch/kento/{project_name}",
        },
        "compiler": {
            "fortran": "gfortran",
            "mpi_fortran": "mpif90",
            "cuda_fortran": "nvfortran",
        },
        "scheduler": {"type": "slurm"},
    }

    write_yaml_if_needed(root / "config" / "machine.yaml", machine, overwrite)
    write_yaml_if_needed(root / "config" / "compiler.yaml", {"compiler": machine["compiler"]}, overwrite)
    write_yaml_if_needed(root / "config" / "paths.yaml", {"paths": machine["paths"]}, overwrite)


def copy_schema(root: Path, schema_path: Path, overwrite: bool) -> None:
    dst = root / "project_schema.yaml"
    if dst.exists() and not overwrite:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(schema_path, dst)


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def init_project(root: Path, schema_path: Path, overwrite: bool = False) -> None:
    schema = read_yaml(schema_path)

    root.mkdir(parents=True, exist_ok=True)
    make_dirs(root, schema)

    copy_schema(root, schema_path, overwrite)
    create_readme(root, schema, overwrite)
    create_gitignore(root, overwrite)
    create_makefile(root, overwrite)
    write_text_if_needed(root / "LICENSE", "All rights reserved.\n", overwrite)

    create_case_index(root, schema, overwrite)
    create_templates(root, schema, overwrite)
    create_config(root, schema, overwrite)

    stamp = root / "docs" / "project_initialized.txt"
    if overwrite or not stamp.exists():
        stamp.write_text(
            f"Project initialized at {datetime.now().isoformat(timespec='seconds')}\n"
            f"Schema: {schema_path.resolve()}\n",
            encoding="utf-8",
        )

    manifest = {
        "project_root": str(root.resolve()),
        "schema": str(schema_path.resolve()),
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "generated_files": [
            "README.md",
            ".gitignore",
            "Makefile",
            "LICENSE",
            "project_schema.yaml",
            "cases/case_index.csv",
            "templates/case_schema.yaml",
            "templates/case_index_schema.yaml",
            "templates/input_schema.yaml",
            "templates/meta_schema.yaml",
            "templates/readme_schema.yaml",
            "templates/parameter_catalog.yaml",
            "config/machine.yaml",
            "config/compiler.yaml",
            "config/paths.yaml",
            "docs/project_initialized.txt",
        ],
    }
    write_text_if_needed(
        root / "docs" / "generated_manifest.json",
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        overwrite,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a schema-driven research project directory structure."
    )

    parser.add_argument(
        "--schema",
        default="project_schema.yaml",
        help="Path to project_schema.yaml. Default: project_schema.yaml",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--project",
        help="Project directory name/path to create, e.g. SphericalShockProject",
    )
    group.add_argument(
        "--here",
        action="store_true",
        help="Create the structure in the current directory.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing generated files.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    schema_path = Path(args.schema)
    if not schema_path.exists():
        raise FileNotFoundError(f"Schema file not found: {schema_path}")

    root = Path.cwd() if args.here else Path(args.project)
    init_project(root, schema_path=schema_path, overwrite=args.overwrite)

    print("[OK] Schema-driven research project directory structure created.")
    print(f"  Root:   {root.resolve()}")
    print(f"  Schema: {schema_path.resolve()}")
    print("")
    print("Generated from project_schema.yaml:")
    print("  directories")
    print("  cases/case_index.csv")
    print("  templates/case_schema.yaml")
    print("  templates/case_index_schema.yaml")
    print("  templates/input_schema.yaml")
    print("  templates/meta_schema.yaml")
    print("  templates/readme_schema.yaml")
    print("")
    print("Next steps:")
    print(f"  cd {root}")
    print("  git init")
    print("  git add .")
    print('  git commit -m "Initial schema-driven project structure"')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
