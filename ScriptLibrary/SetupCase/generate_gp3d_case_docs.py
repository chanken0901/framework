#!/usr/bin/env python3
"""ケースYAMLから人間向け文書とGP3D入力を生成する。

入力は ``cases/<case_id>/case.yaml``、出力は同じディレクトリの
``input.nml``、``README.md``、``case_meta.json`` である。計算そのものや
実行環境のコピーは行わず、ケース内容の確認と記録だけを担当する。
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from case_registry import find_duplicates, sync_case_index
from gp3d_case import CaseValidationError, case_digest, generate_input_namelist, nested
from run_env_builder import BuildEnvironmentError, resolve_solver
from yaml_support import YamlFormatError, load_yaml


FRAMEWORK_ROOT = Path(__file__).resolve().parents[2]


def _write(path: Path, text: str, overwrite: bool, dry_run: bool) -> None:
    """上書き規則とDryRunを共通化してテキスト成果物を書き出す。"""
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists; use --overwrite")
    if dry_run:
        print(f"\n--- {path} ---\n{text}")
        return
    path.write_text(text, encoding="utf-8")
    print(f"[OK] Wrote {path}")


def _readme(case: dict, profile: str) -> str:
    """ケースの主要条件を一覧できる短いMarkdownを組み立てる。"""
    return f"""# {case.get('case_id')}: {case.get('case_label', '')}

{case.get('description', '')}

## 計算条件

- 物理モデル: `{nested(case, 'physics.model')}`
- 初期条件: `{nested(case, 'flow.type')}`
- ソルバープロファイル: `{profile}`
- 格子数: `{nested(case, 'grid.nx')} x {nested(case, 'grid.ny')} x {nested(case, 'grid.nz')}`
- 時間刻み: `{nested(case, 'time.dt')}`
- 実時間発展ステップ数: `{nested(case, 'time.nsteps')}`

`input.nml`と`case_meta.json`は`case.yaml`から自動生成されます。
これらを独立した設定原本として直接編集しないでください。
"""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate GP3D case documents from case.yaml")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--case", help="Case directory containing case.yaml")
    group.add_argument("--yaml", help="Path to case.yaml")
    parser.add_argument(
        "--project",
        default="project_schema.yaml",
        help="Path to project_schema.yaml",
    )
    parser.add_argument("--machine", help="Override machine YAML")
    parser.add_argument("--solver-root", help="Override the GP3D package root")
    parser.add_argument("--solver-library", help="Override the SolverLibrary root")
    parser.add_argument("--index-schema", default="templates/case_index_schema.yaml")
    parser.add_argument("--skip-duplicate-check", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """ケース、プロジェクト、ソルバーマニフェストを読み、3成果物を生成する。"""
    args = parse_args(argv)
    case_path = Path(args.yaml).resolve() if args.yaml else (Path(args.case) / "case.yaml").resolve()
    case_dir = case_path.parent
    try:
        case = load_yaml(case_path)
        project_path = Path(args.project).resolve()
        project = load_yaml(project_path)
        machine_relative = nested(project, "paths.default_machine", "config/machine.local.yaml")
        machine_path = Path(args.machine).resolve() if args.machine else (
            project_path.parent / str(machine_relative)
        ).resolve()
        machine = load_yaml(machine_path)
        solver_id = str(nested(case, "solver.implementation", "gp3d"))
        solver_root, manifest_path = resolve_solver(
            project_path,
            project,
            machine,
            solver_id,
            Path(args.solver_root) if args.solver_root else None,
            Path(args.solver_library) if args.solver_library else None,
        )
        manifest = load_yaml(manifest_path)
        profile = str(nested(case, "solver.profile", ""))
        input_text = generate_input_namelist(case, manifest, profile)
        index_schema_path = project_path.parent / args.index_schema
        index_schema = load_yaml(index_schema_path)
        case_index = project_path.parent / "cases" / "case_index.csv"
        duplicates = find_duplicates(case, index_schema, case_index)
        if duplicates and not args.skip_duplicate_check:
            ids = ", ".join(row.get("case_id", "") for row in duplicates)
            raise ValueError(f"duplicate case conditions already exist: {ids}")
        meta = {
            "schema_version": 1,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "case_id": case.get("case_id"),
            "case_sha256": case_digest(case),
            "solver_id": solver_id,
            "solver_profile": profile,
            "physics": case.get("physics", {}),
            "flow": case.get("flow", {}),
            "grid": case.get("grid", {}),
            "time": case.get("time", {}),
            "output": case.get("output", {}),
        }
        _write(case_dir / "input.nml", input_text, args.overwrite, args.dry_run)
        _write(case_dir / "README.md", _readme(case, profile), args.overwrite, args.dry_run)
        _write(
            case_dir / "case_meta.json",
            json.dumps(meta, ensure_ascii=False, indent=2) + "\n",
            args.overwrite,
            args.dry_run,
        )
        if not args.dry_run:
            sync_case_index(case, index_schema, case_index)
    except (
        BuildEnvironmentError,
        CaseValidationError,
        YamlFormatError,
        FileExistsError,
        KeyError,
        OSError,
        ValueError,
    ) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
