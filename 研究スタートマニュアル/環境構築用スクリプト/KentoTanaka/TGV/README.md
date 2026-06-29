# ResearchProject

Numerical simulation research framework

このプロジェクトは `project_schema.yaml` を Single Source of Truth として初期化されています。

## 基本方針

- ディレクトリ構造，ケース管理表，テンプレートは `project_schema.yaml` から生成する。
- STEP2以降のツールは `templates/*.yaml` を参照してケース生成・重複確認・入力ファイル生成を行う。
- ソースコード，解析コード，論文，計算条件は Git/GitHub で管理する。
- VTK, HDF5, restart, checkpoint などの巨大データは Git 管理しない。
- 巨大データの保存場所は `case.yaml`, `input.dat`, `meta.json`, `case_index.csv` に記録する。

## 生成された主要ディレクトリ

### docs
- `docs/`
- `docs/references/`
- `docs/meeting_notes/`

### source
- `src/`
- `src/core/`
- `src/physics/`
- `src/solver/`
- `src/numerics/`
- `src/mpi/`
- `src/io/`
- `src/analysis/`

### workflow
- `scripts/`
- `config/`
- `templates/`
- `cases/`

### analysis
- `analysis/`
- `analysis/pod/`
- `analysis/dmd/`
- `analysis/information_theory/`
- `analysis/visualization/`

### writing
- `paper/`
- `paper/manuscript/`
- `paper/figures/`
- `paper/tables/`
- `paper/bibliography/`
- `presentation/`
- `presentation/conference/`
- `presentation/seminar/`

### tests
- `tests/`
- `tests/unit/`
- `tests/regression/`
- `tests/benchmark/`

## 初期化後の推奨操作

```bash
git init
git add .
git commit -m "Initial schema-driven project structure"
```
