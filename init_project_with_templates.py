#!/usr/bin/env python3
"""
init_project.py

Research project initializer for simulation-based research.

This script creates the standard directory structure for CFD, quantum fluid,
GPE, HPC, and LaTeX-based paper-writing workflows.

Design policy
-------------
1. Case directories use stable IDs:
     cases/case0001/
     cases/case0002/

2. Physical/numerical parameters are NOT encoded in directory names.
   They are recorded in:
     cases/case_index.csv
     cases/case0001/input.dat
     cases/case0001/meta.json

3. Human-readable labels are recorded as:
     case_label = M2_Re1e5_N512

4. Large data are NOT stored in this project directory.
   They should be stored in external storage such as NAS/HDD/HPC storage.

5. Case templates are stored in:
     cases/templates/

Created structure
-----------------
project_root/
├── src/
│   ├── solver/
│   ├── physics/
│   ├── io/
│   └── mpi/
├── post/
│   ├── common/
│   ├── statistics/
│   └── visualization/
├── paper/
│   └── sections/
├── fig/
│   ├── draft/
│   └── final/
├── cases/
│   ├── templates/
│   │   ├── input.dat.template
│   │   ├── slurm.template
│   │   └── README_template.md
│   └── case_index.csv
├── scripts/
├── config/
├── docs/
├── tests/
├── data_links/
├── README.md
├── .gitignore
└── Makefile

Usage
-----
Create a new project directory:

    python init_project.py --project spherical_shock_project

Initialize the current directory:

    python init_project.py --here

Overwrite generated project files if needed:

    python init_project.py --here --overwrite
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path


CASE_INDEX_COLUMNS = [
    "case_id",
    "case_label",
    "case_dir",
    "status",
    "project_name",
    "description",
    "mach_number",
    "reynolds_number",
    "reynolds_lambda",
    "turbulent_mach_number",
    "prandtl_number",
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
    "git_branch",
    "git_commit_hash",
    "notes",
]


README_TEMPLATE = """# Research Project

This repository is a structured research project for simulation-based research.

## Basic Policy

This project separates:

1. source code,
2. simulation case settings,
3. post-processing scripts,
4. figures,
5. paper manuscripts,
6. large raw data.

Large raw data should not be stored directly in this repository.
Instead, store large data in external storage such as NAS, HDD, or HPC storage,
and record its location in each case's `input.dat` and `meta.json`.

## Directory Structure

```text
project-root/
├── src/          # simulation source code
├── post/         # post-processing and analysis scripts
├── paper/        # LaTeX manuscript and bibliography
├── fig/          # figures for papers and reports
│   ├── draft/    # working figures
│   └── final/    # final figures used in manuscripts
├── cases/        # simulation case directories and case registry
│   ├── templates/
│   └── case_index.csv
├── scripts/      # automation scripts
├── config/       # machine-dependent configuration files
├── docs/         # research notes and documentation
├── tests/        # verification and small test cases
└── data_links/   # optional symbolic links or text files pointing to large data
```

## Case Directory Policy

Case directories use stable IDs:

```text
cases/
├── case0001/
├── case0002/
└── case_index.csv
```

Do not encode physical parameters directly in the directory name.
Instead, record them in:

```text
cases/case_index.csv
cases/case0001/input.dat
cases/case0001/meta.json
```

For human readability, use `case_label`:

```text
case_id    = case0001
case_label = M2_Re1e5_N512
```

This makes it easy to change label rules later without renaming directories.

## Case Templates

The directory

```text
cases/templates/
```

contains templates used by case-generation scripts:

```text
input.dat.template
slurm.template
README_template.md
```

These files define the default structure of newly created cases.

## Case Registry

The file

```text
cases/case_index.csv
```

is the project-level case registry.

It stores a summary of all cases, including:

- case ID,
- case label,
- physical parameters,
- numerical parameters,
- raw data location,
- paper/figure usage,
- Git commit hash.

## Generated Case Files

Each case directory is expected to contain:

```text
input.dat      # single source of truth for the case
README.md      # generated human-readable case description
meta.json      # generated machine-readable metadata
notes.md       # free-form human notes
build.sh       # generated build script
run.sh         # generated run script
job_slurm.sh   # generated Slurm job script
```

## Data Policy

Git-managed:

```text
src/
post/
paper/
fig/
cases/
scripts/
config/
docs/
tests/
README.md
Makefile
.gitignore
```

Not Git-managed:

```text
large raw data
HPC output files
temporary files
binary field data
```

Large data should be stored externally and referenced by path.
"""


GITIGNORE_TEMPLATE = """# =========================
# Build products
# =========================
*.o
*.mod
*.smod
*.a
*.so
*.dylib
*.dll
*.exe
*.out
bin/
build/

# =========================
# Python
# =========================
__pycache__/
*.pyc
.venv/
venv/
.env

# =========================
# LaTeX build products
# =========================
*.aux
*.bbl
*.bcf
*.blg
*.fdb_latexmk
*.fls
*.log
*.out
*.run.xml
*.synctex.gz
*.toc

# Keep source PDF figures, but ignore compiled paper PDFs if desired.
# Uncomment the next line if you do not want to track PDFs.
# *.pdf

# =========================
# Large simulation data
# =========================
*.h5
*.hdf5
*.nc
*.vtk
*.vtu
*.plt
*.dat.bin
*.raw

# =========================
# Case outputs
# =========================
cases/*/output/
cases/*/tmp/
cases/*/logs/
cases/*/restart/
cases/*/checkpoint/

# =========================
# Local machine settings
# =========================
config/machine.local.mk
config/machine.local.sh

# =========================
# OS/editor files
# =========================
.DS_Store
Thumbs.db
.vscode/
.idea/
"""


MAKEFILE_TEMPLATE = """# Project-level Makefile

.PHONY: help case-list clean

help:
\t@echo "Available targets:"
\t@echo "  make case-list   Show registered cases"
\t@echo "  make clean       Remove common temporary files"

case-list:
\t@if [ -f cases/case_index.csv ]; then \\
\t\tcolumn -s, -t < cases/case_index.csv | less -S; \\
\telse \\
\t\techo "cases/case_index.csv not found"; \\
\tfi

clean:
\tfind . -name "__pycache__" -type d -prune -exec rm -rf {} +
\tfind . -name "*.pyc" -delete
"""


INPUT_DAT_TEMPLATE = """# input.dat.template
# This file is the template for each simulation case.
# create_case.py copies this file to cases/caseXXXX/input.dat.
#
# Policy:
#   - input.dat is the single source of truth for each case.
#   - README.md and meta.json are generated from input.dat.
#   - Do not encode parameters in the directory name.
#   - Use case_label for human readability.

# =========================
# Case identity
# =========================
case_id = {case_id}
case_label = {case_label}
case_name = {case_id}
project_name = {project_name}
status = planned
description = short description of this case

# =========================
# Physics
# =========================
mach_number = 2.0
reynolds_number = 1.0e5
reynolds_lambda =
turbulent_mach_number =
prandtl_number = 0.72

# =========================
# Grid
# =========================
nx = 512
ny = 512
nz = 512

# =========================
# Numerics
# =========================
scheme = Roe
reconstruction = WENOZ
time_integration = TVD_RK3

# =========================
# Storage
# =========================
# Large data must not be placed in cases/caseXXXX/.
# Store it externally and record the path below.
raw_data_location = /mnt/storage/{project_name}/{case_id}

# =========================
# Paper / figure usage
# =========================
used_in =

# =========================
# Duplicate check
# =========================
duplicate_check_keys = mach_number, reynolds_number, nx, ny, nz, scheme, reconstruction
"""


SLURM_TEMPLATE = """#!/usr/bin/env bash
#SBATCH --job-name={case_id}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$CASE_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

echo "[INFO] Running on host: $(hostname)"
echo "[INFO] Case directory: $CASE_DIR"

# Edit module commands for your HPC environment.
# module purge
# module load gcc
# module load openmpi
# module load cuda

bash "$CASE_DIR/build.sh"
bash "$CASE_DIR/run.sh"
"""


CASE_README_TEMPLATE = """# {case_id}

> This README is generated from `input.dat`.
> Do not edit this file directly. Edit `input.dat` and rerun the case-generation script.
> Free-form human notes should be written in `notes.md`.

## Case Identity

| Item | Value |
|---|---|
| Case ID | {case_id} |
| Case label | {case_label} |
| Status | {status} |

## Project

{project_name}

## Description

{description}

## Physical Parameters

| Parameter | Value |
|---|---:|
| Mach number | {mach_number} |
| Reynolds number | {reynolds_number} |
| Reynolds lambda | {reynolds_lambda} |
| Turbulent Mach number | {turbulent_mach_number} |
| Prandtl number | {prandtl_number} |

## Grid

| Direction | Number |
|---|---:|
| Nx | {nx} |
| Ny | {ny} |
| Nz | {nz} |

## Numerical Method

| Item | Value |
|---|---|
| Scheme | {scheme} |
| Reconstruction | {reconstruction} |
| Time integration | {time_integration} |

## Raw Data Location

```text
{raw_data_location}
```

## Used In

{used_in}
"""


def create_case_index(path: Path, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        return

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CASE_INDEX_COLUMNS)
        writer.writeheader()


def write_if_needed(path: Path, text: str, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        return
    path.write_text(text, encoding="utf-8")


def initialize_project(root: Path, overwrite: bool = False) -> None:
    directories = [
        "src",
        "src/solver",
        "src/physics",
        "src/io",
        "src/mpi",
        "post",
        "post/common",
        "post/visualization",
        "post/statistics",
        "paper",
        "paper/sections",
        "fig",
        "fig/draft",
        "fig/final",
        "cases",
        "cases/templates",
        "scripts",
        "config",
        "docs",
        "tests",
        "data_links",
    ]

    for directory in directories:
        (root / directory).mkdir(parents=True, exist_ok=True)

    write_if_needed(root / "README.md", README_TEMPLATE, overwrite=overwrite)
    write_if_needed(root / ".gitignore", GITIGNORE_TEMPLATE, overwrite=overwrite)
    write_if_needed(root / "Makefile", MAKEFILE_TEMPLATE, overwrite=overwrite)

    create_case_index(root / "cases" / "case_index.csv", overwrite=overwrite)

    write_if_needed(
        root / "cases" / "templates" / "input.dat.template",
        INPUT_DAT_TEMPLATE,
        overwrite=overwrite,
    )

    write_if_needed(
        root / "cases" / "templates" / "slurm.template",
        SLURM_TEMPLATE,
        overwrite=overwrite,
    )

    write_if_needed(
        root / "cases" / "templates" / "README_template.md",
        CASE_README_TEMPLATE,
        overwrite=overwrite,
    )

    created_file = root / "docs" / "project_initialized.txt"
    if overwrite or not created_file.exists():
        created_file.write_text(
            f"Project initialized at {datetime.now().isoformat(timespec='seconds')}\n",
            encoding="utf-8",
        )

    for directory in [
        "src",
        "post",
        "paper",
        "fig/draft",
        "fig/final",
        "cases",
        "cases/templates",
        "scripts",
        "config",
        "docs",
        "tests",
        "data_links",
    ]:
        keep = root / directory / ".gitkeep"
        if not keep.exists():
            keep.write_text("", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Initialize a standard research project directory structure."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--project",
        help="Project directory to create, e.g. spherical_shock_project",
    )
    group.add_argument(
        "--here",
        action="store_true",
        help="Initialize the current directory.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite README.md, .gitignore, Makefile, case_index.csv, and templates if they already exist.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.here:
        root = Path.cwd()
    else:
        root = Path(args.project)
        root.mkdir(parents=True, exist_ok=True)

    initialize_project(root, overwrite=args.overwrite)

    print("[OK] Project structure initialized.")
    print(f"  Root:       {root.resolve()}")
    print(f"  Case index: {(root / 'cases' / 'case_index.csv').resolve()}")
    print(f"  Templates:  {(root / 'cases' / 'templates').resolve()}")
    print("")
    print("Recommended case naming policy:")
    print("  cases/case0001/")
    print("  cases/case0002/")
    print("  Use case_label in input.dat and case_index.csv for readability.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
