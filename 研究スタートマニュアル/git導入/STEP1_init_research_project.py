#!/usr/bin/env python3
"""
init_research_project.py

研究プロジェクト用の標準ディレクトリ構造を生成するスクリプト。

想定用途:
  - CFD / LES / DNS / SBLI
  - 量子流体 / Gross-Pitaevskii equation
  - MPI / CUDA Fortran / HPC
  - Git/GitHubによるコード管理
  - NASによる巨大データ管理

方針:
  - OS依存を避けるため、Python標準ライブラリのみを使用する。
  - Gitで管理する軽量ファイル群と、NASに置く巨大データを分離する。
  - cases/case_index.csv を最初から作成する。
  - 空ディレクトリもGitで追跡できるように .gitkeep を置く。

使い方:
  新しいプロジェクトを作る場合:
      python init_research_project.py --project MyProject

  現在のディレクトリに構造を作る場合:
      python init_research_project.py --here

  既存ファイルを上書きしたい場合:
      python init_research_project.py --project MyProject --overwrite
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path


CASE_INDEX_COLUMNS = [
    "case_id",
    "case_label",
    "status",
    "description",

    "physics_model",
    "flow_type",
    "solver_type",

    "mach_number",
    "reynolds_number",
    "prandtl_number",

    "nx",
    "ny",
    "nz",

    "flux",
    "reconstruction",
    "time_integration",
    "raw_data_location",
    "used_in",
    "created_at",
    "updated_at",
    "git_commit_hash",
    "notes",

    "tgv_amplitude",
    "tgv_mode_x",
    "tgv_mode_y",
    "tgv_mode_z",
    "hit_reynolds_lambda",
    "hit_turbulent_mach_number",

]


README_TEXT = """# Research Project

このディレクトリは、数値シミュレーション研究のための標準プロジェクト構造です。

## 基本方針

- ソースコード、解析コード、論文、計算条件は Git/GitHub で管理する。
- VTK, HDF5, restart, checkpoint などの巨大データは Git 管理しない。
- 巨大データは NAS や外部ストレージに保存し、その場所を `input.dat`, `case.yaml`, `meta.json` に記録する。
- 計算ケースは `cases/case0001/`, `cases/case0002/` のように安定した番号で管理する。
- 物理パラメータはディレクトリ名ではなく、`case.yaml` や `case_index.csv` に記録する。

## 主要ディレクトリ

```text
src/           ソースコード
scripts/       自動化スクリプト
config/        環境・コンパイラ・パス設定
templates/     ケース生成テンプレート
cases/         計算ケース
analysis/      解析コード
paper/         論文
presentation/  発表資料
docs/          設計メモ・研究ノート
tests/         テスト
```
"""


GITIGNORE_TEXT = """# Build products
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


MAKEFILE_TEXT = """# Project-level Makefile

.PHONY: help case-list clean

help:
\t@echo "Available targets:"
\t@echo "  make case-list   Show registered cases"
\t@echo "  make clean       Remove common temporary files"

case-list:
\t@if [ -f cases/case_index.csv ]; then \\
\t\tcat cases/case_index.csv; \\
\telse \\
\t\techo "cases/case_index.csv not found"; \\
\tfi

clean:
\tfind . -name "__pycache__" -type d -prune -exec rm -rf {} + 2>/dev/null || true
\tfind . -name "*.pyc" -delete 2>/dev/null || true
"""


CASE_YAML_TEMPLATE = """# case.yaml
# このファイルをケースの正本として使う想定です。
# generate_case_docs.py などで input.dat, README.md, meta.json を自動生成します。

case_id: case0001
case_label: baseline
status: planned
description: baseline case

project:
  name: ResearchProject
  author: Kento Tanaka

physics:
  model: NSE
  nse:
    gamma: 1.4
    mach_number: 0.5
    reynolds_number: 1.0e5
    prandtl_number: 0.72

flow:
  type: TGV
  tgv:
    amplitude: 1.0
    mode_x: 1
    mode_y: 1
    mode_z: 1

solver:
  type: nse_fvm

grid:
  nx: 512
  ny: 512
  nz: 512

numerics:
  flux: Roe
  reconstruction: WENOZ
  time_integration: TVD_RK3

storage:
  raw_data_location: /mnt/data_nas/ResearchProject/case0001

outputs:
  used_in:
    - Figure 1
"""


INPUT_TEMPLATE = """# input.dat
# ソルバーが実行時に読み込む計算条件ファイル。
# 通常は case.yaml から自動生成することを推奨。

case_id = case0001
case_label = baseline
physics_model = NSE
flow_type = TGV
solver_type = nse_fvm

gamma = 1.4
mach_number = 0.5
reynolds_number = 1.0e5
prandtl_number = 0.72

nx = 512
ny = 512
nz = 512

flux = Roe
reconstruction = WENOZ
time_integration = TVD_RK3

raw_data_location = /mnt/data_nas/ResearchProject/case0001
"""


README_TEMPLATE = """# {case_id}

このREADMEは、各ケースの説明用テンプレートです。

## Case Information

- Case ID:
- Case label:
- Status:
- Description:

## Physical parameters

- Mach number:
- Reynolds number:
- Grid:

## Data location

```text
/mnt/data_nas/ProjectName/case0001
```

## Notes

自由記述は `notes.md` に記入します。
"""


META_TEMPLATE = """{
  "case": {
    "case_id": "case0001",
    "case_label": "baseline",
    "status": "planned"
  },
  "physics": {
    "model": "NSE",
    "mach_number": 0.5,
    "reynolds_number": 100000.0,
    "prandtl_number": 0.72
  },
  "flow": {
      "type": "TGV"
  },
  "solver": {
      "type": "nse_fvm"
  },
  "grid": {
    "nx": 512,
    "ny": 512,
    "nz": 512
  },
  "storage": {
    "raw_data_location": "/mnt/data_nas/ResearchProject/case0001"
  },
  "git": {
    "commit_hash": ""
  }
}
"""


CONFIG_MACHINE = """# machine.yaml
# 実行環境ごとの設定を書くファイル。
# 例: Linux共有サーバー、HPC、ローカルPCなど。

machine_name: local_linux_server

paths:
  repo_root: /mnt/repo_nas/Kento/ResearchRepo/ProjectName
  data_root: /mnt/data_nas/ProjectName
  scratch_root: /scratch/kento/ProjectName

compiler:
  fortran: gfortran
  mpi_fortran: mpif90
  cuda_fortran: nvfortran

scheduler:
  type: slurm
"""


def write_text_if_needed(path: Path, text: str, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        return
    path.write_text(text, encoding="utf-8")


def make_dirs(root: Path) -> None:
    directories = [
        "docs",
        "docs/references",
        "docs/meeting_notes",

        "src",
        "src/core",
        "src/physics",
        "src/solver",
        "src/numerics",
        "src/mpi",
        "src/io",
        "src/analysis",

        "scripts",
        "config",
        "templates",
        "cases",
        "analysis",
        "analysis/pod",
        "analysis/dmd",
        "analysis/information_theory",
        "analysis/visualization",
        "paper",
        "paper/manuscript",
        "paper/figures",
        "paper/tables",
        "paper/bibliography",
        "presentation",
        "presentation/conference",
        "presentation/seminar",
        "tests",
        "tests/unit",
        "tests/regression",
        "tests/benchmark",
    ]

    for d in directories:
        p = root / d
        p.mkdir(parents=True, exist_ok=True)
        gitkeep = p / ".gitkeep"
        if not gitkeep.exists():
            gitkeep.write_text("", encoding="utf-8")


def create_case_index(root: Path, overwrite: bool = False) -> None:
    path = root / "cases" / "case_index.csv"
    if path.exists() and not overwrite:
        return

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CASE_INDEX_COLUMNS)
        writer.writeheader()


def create_templates(root: Path, overwrite: bool = False) -> None:
    write_text_if_needed(root / "templates" / "case.yaml", CASE_YAML_TEMPLATE, overwrite)
    write_text_if_needed(root / "templates" / "input_template.dat", INPUT_TEMPLATE, overwrite)
    write_text_if_needed(root / "templates" / "README_template.md", README_TEMPLATE, overwrite)
    write_text_if_needed(root / "templates" / "meta_template.json", META_TEMPLATE, overwrite)


def create_config(root: Path, overwrite: bool = False) -> None:
    write_text_if_needed(root / "config" / "machine.yaml", CONFIG_MACHINE, overwrite)
    write_text_if_needed(root / "config" / "compiler.yaml", "# compiler settings\n", overwrite)
    write_text_if_needed(root / "config" / "paths.yaml", "# path settings\n", overwrite)


def init_project(root: Path, overwrite: bool = False) -> None:
    root.mkdir(parents=True, exist_ok=True)

    make_dirs(root)

    write_text_if_needed(root / "README.md", README_TEXT, overwrite)
    write_text_if_needed(root / ".gitignore", GITIGNORE_TEXT, overwrite)
    write_text_if_needed(root / "Makefile", MAKEFILE_TEXT, overwrite)
    write_text_if_needed(root / "LICENSE", "All rights reserved.\n", overwrite)

    create_case_index(root, overwrite)
    create_templates(root, overwrite)
    create_config(root, overwrite)

    stamp = root / "docs" / "project_initialized.txt"
    if overwrite or not stamp.exists():
        stamp.write_text(
            f"Project initialized at {datetime.now().isoformat(timespec='seconds')}\n",
            encoding="utf-8",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a standard research project directory structure."
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
        help="Overwrite existing template files, README.md, .gitignore, Makefile, etc.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.here:
        root = Path.cwd()
    else:
        root = Path(args.project)

    init_project(root, overwrite=args.overwrite)

    print("[OK] Research project directory structure created.")
    print(f"  Root: {root.resolve()}")
    print("")
    print("Main directories:")
    print("  src/           source code")
    print("  cases/         simulation cases")
    print("  analysis/      post-processing")
    print("  paper/         manuscript")
    print("  scripts/       automation scripts")
    print("  docs/          documentation")
    print("")
    print("Next steps:")
    print(f"  cd {root}")
    print("  git init")
    print("  git add .")
    print('  git commit -m "Initial project structure"')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
