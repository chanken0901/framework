# Project Name

（例: Spherical Shock Wave Turbulence Interaction）

---

## Overview

本プロジェクトは〇〇に関する数値シミュレーションおよび解析を行うための研究プロジェクトである。

### Research Objective

研究目的を簡潔に記載する。

例：

* 球面衝撃波と等方性乱流の相互作用機構の解明
* 量子渦と非線形波動の相互生成機構の解明
* 衝撃波面変形の統計則の導出

---

## Directory Structure

```text
project-root/
├── src/
├── post/
├── cases/
├── fig/
├── paper/
├── scripts/
├── config/
├── docs/
├── tests/
├── README.md
├── meta.json
└── .gitignore
```

### src/

数値計算コード本体

例

```text
src/
├── solver/
├── physics/
├── io/
├── mpi/
└── main.f90
```

---

### post/

解析コード

例

```text
post/
├── spectrum/
├── statistics/
├── visualization/
└── common/
```

---

### cases/

計算条件

例

```text
cases/
└── case0001_M2_Re1e5/
    ├── input.dat
    ├── meta.json
    └── README.md
```

注意：

巨大データは保存しない。

---

### fig/

論文図

```text
fig/
├── draft/
└── final/
```

draft:
作業中

final:
論文採用版

---

### paper/

論文原稿

```text
paper/
├── main.tex
├── refs.bib
├── sections/
└── figures/
```

---

## Software Requirements

### Compiler

```text
gfortran >= 12
```

または

```text
nvfortran >= 24
```

### MPI

```text
OpenMPI
```

または

```text
MPICH
```

### Python

```text
Python >= 3.11
```

主要ライブラリ

```text
numpy
scipy
matplotlib
h5py
pandas
pyvista
```

---

## Build

### CPU version

```bash
make cpu
```

### MPI version

```bash
make mpi
```

### GPU version

```bash
make gpu
```

---

## Run

### Local execution

```bash
./solver input.dat
```

### MPI execution

```bash
mpirun -np 16 ./solver input.dat
```

### Slurm execution

```bash
sbatch run.sh
```

---

## Post-processing

解析実行

```bash
python post/main.py
```

または

```bash
make post
```

---

## Figure Generation

図作成

```bash
python post/plot_spectrum.py
```

出力先

```text
fig/draft/
```

論文採用時

```text
fig/final/
```

へ移動する。

---

## Data Policy

### Git管理対象

```text
src/
post/
paper/
fig/
cases/
scripts/
config/
```

### Git管理対象外

```text
*.h5
*.hdf5
*.vtk
*.vtu
output/
data_raw/
```

---

## External Data Storage

巨大データ保存場所

例

```text
/mnt/storage/project_name/
```

または

```text
NAS/project_name/
```

GitHubには保存しない。

---

## Reproducibility

再現に必要な情報

* 入力条件
* コンパイル条件
* Git commit hash
* Python環境
* 使用ライブラリ

を保存すること。

---

## Publications

関連論文

1. Author, Journal, Year
2. Author, Journal, Year

---

## Authors

Principal Investigator

```text
Kento Tanaka
```

Collaborators

```text
Name
Affiliation
```

---

## License

研究室内利用

または

```text
MIT License
```

---

## Contact

担当者

```text
email@example.com
```
