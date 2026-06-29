#!/usr/bin/env python3
"""
check_ubuntu_research_env.py

Ubuntu research server software checklist generator.

This script checks whether recommended software for a simulation-based
research server is installed, and generates checklist reports.

Target environment:
  - Ubuntu / Debian-like Linux server
  - CFD, quantum fluid, GPE, MPI, GPU, LaTeX, Git/GitHub workflow

Generated files:
  - software_checklist.md
  - software_checklist.csv
  - software_checklist.json

Usage:
  python check_ubuntu_research_env.py

Specify output directory:
  python check_ubuntu_research_env.py --out reports

Notes:
  - This script does not install anything.
  - It only checks availability and versions.
  - It uses only Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional


SOFTWARE = [
    {
        "name": "Git",
        "category": "Version control",
        "commands": ["git"],
        "version_cmd": ["git", "--version"],
        "apt": "sudo apt install git",
        "usage": "git status / git commit / git push",
        "description": "ソースコード・論文・解析スクリプトの履歴管理に使う。",
        "priority": "必須",
    },
    {
        "name": "OpenSSH Server",
        "category": "Remote access",
        "commands": ["sshd", "ssh"],
        "version_cmd": ["ssh", "-V"],
        "apt": "sudo apt install openssh-server",
        "usage": "ssh user@server",
        "description": "Windowsなどの外部PCからUbuntu共有サーバーへ接続するために使う。",
        "priority": "必須",
    },
    {
        "name": "tmux",
        "category": "Terminal",
        "commands": ["tmux"],
        "version_cmd": ["tmux", "-V"],
        "apt": "sudo apt install tmux",
        "usage": "tmux new -s work",
        "description": "SSH切断後も解析やジョブを継続するために使う。",
        "priority": "必須",
    },
    {
        "name": "htop",
        "category": "Monitoring",
        "commands": ["htop"],
        "version_cmd": ["htop", "--version"],
        "apt": "sudo apt install htop",
        "usage": "htop",
        "description": "CPU・メモリ使用量を確認する。",
        "priority": "推奨",
    },
    {
        "name": "tree",
        "category": "Utility",
        "commands": ["tree"],
        "version_cmd": ["tree", "--version"],
        "apt": "sudo apt install tree",
        "usage": "tree cases",
        "description": "ディレクトリ構造を見やすく表示する。",
        "priority": "推奨",
    },
    {
        "name": "rsync",
        "category": "File transfer",
        "commands": ["rsync"],
        "version_cmd": ["rsync", "--version"],
        "apt": "sudo apt install rsync",
        "usage": "rsync -av src/ dst/",
        "description": "共有サーバー・NAS・HPC間でファイル同期する。",
        "priority": "必須",
    },
    {
        "name": "Python",
        "category": "Python",
        "commands": ["python3"],
        "version_cmd": ["python3", "--version"],
        "apt": "sudo apt install python3",
        "usage": "python3 script.py",
        "description": "解析スクリプトや自動化スクリプトの実行に使う。",
        "priority": "必須",
    },
    {
        "name": "pip",
        "category": "Python",
        "commands": ["pip3"],
        "version_cmd": ["pip3", "--version"],
        "apt": "sudo apt install python3-pip",
        "usage": "pip3 install package",
        "description": "Pythonパッケージ管理に使う。",
        "priority": "推奨",
    },
    {
        "name": "conda / mamba",
        "category": "Python",
        "commands": ["conda", "mamba"],
        "version_cmd": None,
        "apt": "Miniforge installer recommended",
        "usage": "conda create -n research python=3.11",
        "description": "Python解析環境を分離して管理する。Miniforge推奨。",
        "priority": "推奨",
    },
    {
        "name": "JupyterLab",
        "category": "Python",
        "commands": ["jupyter"],
        "version_cmd": ["jupyter", "lab", "--version"],
        "apt": "conda install jupyterlab",
        "usage": "jupyter lab",
        "description": "対話的な解析・可視化に使う。",
        "priority": "任意",
    },
    {
        "name": "GCC",
        "category": "Compiler",
        "commands": ["gcc", "g++"],
        "version_cmd": ["gcc", "--version"],
        "apt": "sudo apt install build-essential",
        "usage": "gcc --version",
        "description": "C/C++コンパイラ。ライブラリビルドにも必要。",
        "priority": "必須",
    },
    {
        "name": "gfortran",
        "category": "Compiler",
        "commands": ["gfortran"],
        "version_cmd": ["gfortran", "--version"],
        "apt": "sudo apt install gfortran",
        "usage": "gfortran main.f90",
        "description": "Fortranコードのコンパイルに使う。",
        "priority": "必須",
    },
    {
        "name": "OpenMPI",
        "category": "MPI",
        "commands": ["mpirun", "mpif90"],
        "version_cmd": ["mpirun", "--version"],
        "apt": "sudo apt install openmpi-bin libopenmpi-dev",
        "usage": "mpirun -np 16 ./solver input.dat",
        "description": "MPI並列計算に使う。",
        "priority": "必須",
    },
    {
        "name": "FFTW",
        "category": "Numerical library",
        "commands": ["fftw-wisdom"],
        "version_cmd": ["fftw-wisdom", "--version"],
        "apt": "sudo apt install libfftw3-dev",
        "usage": "link with -lfftw3",
        "description": "高速FFTライブラリ。擬スペクトル法やGPE計算で重要。",
        "priority": "推奨",
    },
    {
        "name": "HDF5 tools",
        "category": "Data",
        "commands": ["h5dump"],
        "version_cmd": ["h5dump", "-V"],
        "apt": "sudo apt install hdf5-tools libhdf5-dev",
        "usage": "h5dump file.h5",
        "description": "HDF5ファイルの確認・大規模データ保存に使う。",
        "priority": "推奨",
    },
    {
        "name": "ParaView",
        "category": "Visualization",
        "commands": ["paraview"],
        "version_cmd": ["paraview", "--version"],
        "apt": "sudo apt install paraview",
        "usage": "paraview",
        "description": "VTK/HDF5などの3次元可視化に使う。",
        "priority": "任意",
    },
    {
        "name": "ffmpeg",
        "category": "Visualization",
        "commands": ["ffmpeg"],
        "version_cmd": ["ffmpeg", "-version"],
        "apt": "sudo apt install ffmpeg",
        "usage": "ffmpeg -i img%04d.png movie.mp4",
        "description": "画像列から動画を作る。",
        "priority": "推奨",
    },
    {
        "name": "TeX Live / pdflatex",
        "category": "LaTeX",
        "commands": ["pdflatex"],
        "version_cmd": ["pdflatex", "--version"],
        "apt": "sudo apt install texlive-full",
        "usage": "pdflatex main.tex",
        "description": "LaTeX論文のコンパイルに使う。",
        "priority": "推奨",
    },
    {
        "name": "latexmk",
        "category": "LaTeX",
        "commands": ["latexmk"],
        "version_cmd": ["latexmk", "--version"],
        "apt": "sudo apt install latexmk",
        "usage": "latexmk -pdf main.tex",
        "description": "LaTeXの自動コンパイルに使う。",
        "priority": "推奨",
    },
    {
        "name": "biber",
        "category": "LaTeX",
        "commands": ["biber"],
        "version_cmd": ["biber", "--version"],
        "apt": "sudo apt install biber",
        "usage": "biber main",
        "description": "biblatexの参考文献処理に使う。",
        "priority": "任意",
    },
    {
        "name": "GitHub CLI",
        "category": "GitHub",
        "commands": ["gh"],
        "version_cmd": ["gh", "--version"],
        "apt": "sudo apt install gh",
        "usage": "gh auth login",
        "description": "GitHub操作をコマンドラインから行う。",
        "priority": "任意",
    },
    {
        "name": "CUDA Toolkit",
        "category": "GPU",
        "commands": ["nvcc"],
        "version_cmd": ["nvcc", "--version"],
        "apt": "Install from NVIDIA official repository",
        "usage": "nvcc --version",
        "description": "NVIDIA GPU向け開発環境。",
        "priority": "GPU使用時必須",
    },
    {
        "name": "NVIDIA SMI",
        "category": "GPU",
        "commands": ["nvidia-smi"],
        "version_cmd": ["nvidia-smi"],
        "apt": "Install NVIDIA driver",
        "usage": "nvidia-smi",
        "description": "GPU状態・ドライバ・メモリ使用量を確認する。",
        "priority": "GPU使用時必須",
    },
    {
        "name": "Nsight Systems",
        "category": "GPU",
        "commands": ["nsys"],
        "version_cmd": ["nsys", "--version"],
        "apt": "Install from NVIDIA",
        "usage": "nsys profile ./solver",
        "description": "GPU/CPU性能解析に使う。",
        "priority": "任意",
    },
    {
        "name": "NFS client",
        "category": "NAS",
        "commands": ["mount.nfs"],
        "version_cmd": None,
        "apt": "sudo apt install nfs-common",
        "usage": "sudo mount -t nfs server:/share /mnt/nas",
        "description": "NASをNFSでマウントするために使う。",
        "priority": "NAS使用時必須",
    },
    {
        "name": "Apptainer",
        "category": "Container",
        "commands": ["apptainer"],
        "version_cmd": ["apptainer", "--version"],
        "apt": "sudo apt install apptainer",
        "usage": "apptainer run image.sif",
        "description": "HPC向けコンテナ環境。再現性向上に有効。",
        "priority": "発展",
    },
]


def run_cmd(cmd: List[str]) -> str:
    try:
        out = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5,
        )
        return out.strip().splitlines()[0] if out.strip() else ""
    except Exception as exc:
        return ""


def command_exists(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def check_item(item: Dict) -> Dict:
    commands = item["commands"]
    found_commands = [cmd for cmd in commands if command_exists(cmd)]
    missing_commands = [cmd for cmd in commands if not command_exists(cmd)]

    installed = len(found_commands) > 0

    version = ""
    if installed:
        version_cmd = item.get("version_cmd")
        if version_cmd:
            version = run_cmd(version_cmd)
        else:
            version = f"found: {', '.join(found_commands)}"

    return {
        "name": item["name"],
        "category": item["category"],
        "priority": item["priority"],
        "installed": installed,
        "status": "OK" if installed else "MISSING",
        "found_commands": ", ".join(found_commands),
        "missing_commands": ", ".join(missing_commands),
        "version": version,
        "description": item["description"],
        "install": item["apt"],
        "usage": item["usage"],
    }


def collect_system_info() -> Dict:
    return {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "python": platform.python_version(),
    }


def write_markdown(path: Path, system_info: Dict, rows: List[Dict]) -> None:
    ok = sum(1 for r in rows if r["installed"])
    total = len(rows)

    lines = []
    lines.append("# Ubuntu Research Server Software Checklist")
    lines.append("")
    lines.append("## System Information")
    lines.append("")
    for k, v in system_info.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Installed: {ok} / {total}")
    lines.append(f"- Missing: {total - ok} / {total}")
    lines.append("")
    lines.append("## Checklist")
    lines.append("")
    lines.append("| Status | Software | Priority | Category | Version / Found | Install | Usage |")
    lines.append("|---|---|---|---|---|---|---|")

    for r in rows:
        mark = "OK" if r["installed"] else "MISSING"
        version = r["version"].replace("|", "\\|")
        lines.append(
            f"| {mark} | {r['name']} | {r['priority']} | {r['category']} | "
            f"{version} | `{r['install']}` | `{r['usage']}` |"
        )

    lines.append("")
    lines.append("## Details")
    lines.append("")

    for r in rows:
        lines.append(f"### {r['name']}")
        lines.append("")
        lines.append(f"- Status: {r['status']}")
        lines.append(f"- Priority: {r['priority']}")
        lines.append(f"- Category: {r['category']}")
        lines.append(f"- Description: {r['description']}")
        lines.append(f"- Found commands: {r['found_commands']}")
        lines.append(f"- Missing commands: {r['missing_commands']}")
        lines.append(f"- Version: {r['version']}")
        lines.append(f"- Install: `{r['install']}`")
        lines.append(f"- Usage: `{r['usage']}`")
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def write_csv(path: Path, rows: List[Dict]) -> None:
    fields = [
        "status",
        "name",
        "priority",
        "category",
        "version",
        "description",
        "install",
        "usage",
        "found_commands",
        "missing_commands",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({field: r.get(field, "") for field in fields})


def write_json(path: Path, system_info: Dict, rows: List[Dict]) -> None:
    payload = {
        "system": system_info,
        "software": rows,
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check recommended software for an Ubuntu research server."
    )
    parser.add_argument(
        "--out",
        default="software_check",
        help="Output directory. Default: software_check",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    system_info = collect_system_info()
    rows = [check_item(item) for item in SOFTWARE]

    write_markdown(out / "software_checklist.md", system_info, rows)
    write_csv(out / "software_checklist.csv", rows)
    write_json(out / "software_checklist.json", system_info, rows)

    ok = sum(1 for r in rows if r["installed"])
    total = len(rows)

    print("[OK] Software checklist generated.")
    print(f"  Installed: {ok} / {total}")
    print(f"  Missing:   {total - ok} / {total}")
    print(f"  Markdown:  {out / 'software_checklist.md'}")
    print(f"  CSV:       {out / 'software_checklist.csv'}")
    print(f"  JSON:      {out / 'software_checklist.json'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
