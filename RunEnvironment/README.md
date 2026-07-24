# NAS外部実行環境ジェネレーター

NASは`SolverLibrary`と`ScriptLibrary`の正本だけを保持します。計算時は
`environment.yaml`に従って、必要なソースとスクリプトだけをワークステーションの
ローカルSSD、またはスパコンのscratchへコピーします。生成後のビルドと実行では
NASを参照しません。

## 処理の流れ

1. モデルとプロファイルを選択する。
2. `solver_manifest.yaml`から必要コンポーネントを解決する。
3. Fortranモジュール依存関係を検査する。
4. NAS外へソース、CMake、共通ビルダー、machine設定をコピーする。
5. `SetupCase/create_case_from_template.py`でcaseを生成する。
6. `case.yaml`から`input.dat`または`input.nml`を生成する。
7. 搬送元、SHA-256、Gitコミットを`provenance.json`へ記録する。
8. ローカルコピーだけを使ってビルド・実行する。

## ワークステーション

NAS上の[`environment.yaml`](environment.yaml)はテンプレートとして保ち、計算機側へ
コピーした設計書の`source.framework_root`と`destination.root`を変更します。

```powershell
$tool = "\\Mozart\share\研究フレームワーク構築\ScriptLibrary\RunEnvironment"
New-Item -ItemType Directory -Force C:\ResearchDesigns | Out-Null
Copy-Item "$tool\environment.yaml" C:\ResearchDesigns\gpe_case0001.yaml

python "$tool\prepare_environment.py" `
  C:\ResearchDesigns\gpe_case0001.yaml --dry-run
python "$tool\prepare_environment.py" `
  C:\ResearchDesigns\gpe_case0001.yaml
```

生成先へ移動して実行します。

```powershell
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

`--all`はビルドから実行までを連続して行います。`model.include_tests: true`の場合は
CTestも含みます。実行時のカレントディレクトリは`cases/<case_id>`なので、相対指定
された`output`はcaseフォルダ内へ生成されます。

## モデルの変更

NSE:

```yaml
model:
  name: nse
  profile: cpu_mpi
case:
  template: ScriptLibrary/RunEnvironment/case_templates/nse_taylor_green.yaml
```

GPE:

```yaml
model:
  name: gpe
  profile: cpu_mpi_fftw
case:
  template: ScriptLibrary/RunEnvironment/case_templates/gpe_quantum_taylor_green.yaml
```

GPEでは`cpu_serial_dft`、`cpu_serial_fftw`、`cpu_mpi_dft`、
`cpu_mpi_fftw`、`cuda_single`、`cuda_mpi_cufftmp`を選択できます。

## スパコン

NASがスパコンからマウントされている場合は
[`environment.hpc.yaml`](environment.hpc.yaml)を編集してログインノードで生成します。

```bash
python3 prepare_environment.py environment.hpc.yaml --dry-run
python3 prepare_environment.py environment.hpc.yaml
cd "$SCRATCH/gpe_case0001"
python3 tools/run_case.py --validate-only
sbatch submit.slurm
```

標準の`submit.slurm`は計算ノード上でビルドしてから実行します。ビルドをログインノードで
済ませる運用に変更する場合は、`submit.slurm`の`--build`を削除してください。

NASをマウントできない場合は、ワークステーションで持ち運び用アーカイブを作ります。

```powershell
python .\prepare_environment.py .\environment.hpc.yaml `
  --framework-root "\\Mozart\share\研究フレームワーク構築" `
  --output "C:\ResearchRuns\gpe_case0001" `
  --archive --archive-format gztar
```

生成された`tar.gz`だけをスパコンへ転送してscratchで展開します。スパコンの
コンパイラ、MPI、FFTW、CUDA、cuFFTMpの場所は、転送前または展開後に
`config/machine.yaml`へ記述します。

## 更新と再生成

生成環境は正本ではありません。ソルバーやcase条件を変更した場合はNAS側を更新し、
`--overwrite`で再生成します。削除対象がこのツールの生成物であることを示す
`.generated_run_environment.json`がないディレクトリは上書きしません。
