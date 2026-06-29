#!/usr/bin/env python3
"""
init_project.py

Research project initializer.

This script creates the standard directory structure for a simulation-based
research project. It is designed for CFD, quantum fluid, GPE, HPC, and
LaTeX-based paper-writing workflows.

Design policy
-------------
- Case directories use stable IDs:
    cases/case0001/
    cases/case0002/

- Physical/numerical conditions are NOT encoded in directory names.
  They are recorded in:
    cases/case_index.csv
    cases/case0001/input.dat
    cases/case0001/meta.json

- Large data are NOT stored in this project directory.
  They should be stored in an external storage location such as NAS/HDD/HPC storage.

Created structure
-----------------
project_root/
├── src/
├── post/
├── paper/
├── fig/
│   ├── draft/
│   └── final/
├── cases/
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

Overwrite README, .gitignore, Makefile, or case_index.csv if needed:

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

This makes it easy to change the label rule later without renaming directories.

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

.PHONY: help init-list case-list clean

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

    created_file = root / "docs" / "project_initialized.txt"
    if overwrite or not created_file.exists():
        created_file.write_text(
            f"Project initialized at {datetime.now().isoformat(timespec='seconds')}\n",
            encoding="utf-8",
        )

    # Placeholders to keep empty directories visible in Git if desired.
    for directory in [
        "src",
        "post",
        "paper",
        "fig/draft",
        "fig/final",
        "cases",
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
        help="Overwrite README.md, .gitignore, Makefile, and case_index.csv if they already exist.",
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
    print("")
    print("Next steps:")
    print("  1. Put this script in scripts/ if desired.")
    print("  2. Create a case with create_case.py.")
    print("  3. Initialize Git in the project root.")
    print("")
    print("Recommended case naming policy:")
    print("  cases/case0001/")
    print("  cases/case0002/")
    print("  Use case_label in input.dat and case_index.csv for readability.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
