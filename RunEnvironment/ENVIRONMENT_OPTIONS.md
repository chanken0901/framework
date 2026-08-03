# 実行環境設計書の選択肢一覧

## NSE CUDA

単一GPU版は次の組み合わせで選択します。

```yaml
select:
  model: nse_cuda_single
  execution: release

parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

`nse_cuda_single`はsolver profile `cuda_single`を選択します。CPU/MPI版と同じ
`environment.nse.yaml`をひな型として使い、`select.model`と`parallel`を
切り替えます。

```powershell
$tool = "C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\RunEnvironment"
$design = "$env:USERPROFILE\ResearchRuns\Designs\nse_cuda.yaml"

Copy-Item "$tool\environment.nse.yaml" $design
```

コピーした設計書を次のように変更します。

```yaml
select:
  model: nse_cuda_single
  execution: release

parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

```powershell
python "$tool\prepare_environment.py" $design --dry-run
python "$tool\prepare_environment.py" $design
```

生成後は`run_case.py --prepare`、`--validate-only`、`--build`、`--run`の順に
実行します。単一GPU版は2次・6次KEEP、SSPRK3、周期境界、6次粘性項に対応します。

`environment.gpe.yaml`または`environment.nse.yaml`では、各項目を値のまとまりを
表すIDで選択できます。
標準候補の正本は`environment_options.yaml`です。この文書は人が確認するための
一覧であり、スクリプトが実際に読むのはYAMLカタログです。

## 一覧表示コマンド

標準候補と設計書から追加された候補をまとめて表示します。

```powershell
python .\prepare_environment.py .\environment.gpe.yaml --list-options
```

特定の分類だけを表示できます。

```powershell
python .\prepare_environment.py .\environment.gpe.yaml --list-options source
python .\prepare_environment.py .\environment.gpe.yaml --list-options model
python .\prepare_environment.py .\environment.gpe.yaml --list-options destination
```

外部ツールで処理するときはJSONで出力できます。

```powershell
python .\prepare_environment.py .\environment.gpe.yaml `
  --list-options model --options-format json
```

## 標準候補

### source

| ID | 用途 |
|---|---|
| `framework_relative` | 設計書を`RunEnvironment`内に置く場合にローカルFrameWorkを相対参照 |
| `local_framework` | `C:\Users\Owner\Documents\Codex\FrameWork`を参照 |
| `mozart_nas` | 同期先の`\\Mozart\share\FrameWork`を参照 |
| `hpc_mounted_framework` | 特別な運用でLinuxの`/mnt/framework/FrameWork`を参照 |

### destination

| ID | 用途 |
|---|---|
| `windows_research_runs` | `${USERPROFILE}/ResearchRuns`へ自動採番で生成し、直下の`case_index.csv`へ登録 |
| `local_generated` | 設計書と同じ場所の`generated`へ生成し、同階層の共通台帳へ登録 |
| `hpc_scratch` | `${SCRATCH}`へ生成し、`${SCRATCH}/case_index.csv`へ登録 |

### model

| ID | モデルとソルバープロファイル |
|---|---|
| `nse_cpu_mpi` | NSE / `cpu_mpi` |
| `nse_cpu_mpi_2decomp_fftw` | NSE / `cpu_mpi_2decomp_fftw`、2DECOMP&FFT分散初期化・Forcing対応 |
| `nse_cuda_single` | NSE / `cuda_single`、単一GPU・cuFFT Forcing対応 |
| `gpe_cpu_serial_dft` | GPE / `cpu_serial_dft` |
| `gpe_cpu_serial_fftw` | GPE / `cpu_serial_fftw` |
| `gpe_cpu_mpi_dft` | GPE / `cpu_mpi_dft` |
| `gpe_cpu_mpi_fftw` | GPE / `cpu_mpi_fftw` |
| `gpe_cuda_single` | GPE / `cuda_single` |
| `gpe_cuda_mpi_cufftmp` | GPE / `cuda_mpi_cufftmp` |

### target

| ID | 用途 |
|---|---|
| `windows_gnu_msmpi` | Windows、GNU Fortran、Microsoft MPI |
| `linux_gnu_mpi` | Linux、GNU Fortran、システムMPI |
| `linux_hpc_slurm` | Linuxスパコン、Slurm |

### case

| ID | 用途 |
|---|---|
| `nse_case` | NSE共通case。流れ場は生成後の`case.yaml`で選択 |
| `gpe_quantum_taylor_green` | GPE量子Taylor-Green渦 |

case候補には対応モデルの条件があります。例えば、
`gpe_quantum_taylor_green`と`nse_cpu_mpi`を同時に選ぶと生成前にエラーになります。

### execution

| ID | ビルド構成 |
|---|---|
| `debug` | Debug |
| `release` | Release |
| `relwithdebinfo` | RelWithDebInfo |

MPI、OpenMP、CUDAの使用可否は選択肢IDではなく、設計書の`parallel`で指定します。
MPIプロセス数とOpenMPスレッド数はここには記述しません。

```yaml
parallel:
  use_mpi: true
  use_openmp: true
  use_cuda: false
```

### scheduler

| ID | 用途 |
|---|---|
| `disabled` | ジョブ投入スクリプトを生成しない |
| `slurm` | 並列数を`sbatch`実行時に指定する汎用Slurmスクリプト |

### archive

| ID | 用途 |
|---|---|
| `none` | アーカイブを作成しない |
| `zip` | ZIPを作成 |
| `gztar` | tar.gzを作成 |

## 並列設定の分担

並列機能の有無と実行時並列数は、次のように分けて管理します。

| 設定キー | 設定場所 | 役割 |
|---|---|---|
| `parallel.use_mpi` | 実行環境設計書 | MPI対応ソースと実行方式を選ぶ |
| `parallel.use_openmp` | 実行環境設計書 | OpenMP局所ループを有効化する |
| `parallel.use_cuda` | 実行環境設計書 | CUDAバックエンドを選ぶ |
| `solver.mpi_processes` | 生成後の`case.yaml` | 実行時のMPIプロセス数を指定する |
| `solver.omp_threads` | 生成後の`case.yaml` | MPIランク当たりのOpenMPスレッド数を指定する |

MPI/CUDAの使用有無を変更した場合は、実行環境を再生成してビルドし直します。
`solver.mpi_processes`と`solver.omp_threads`だけを変更する場合は、再ビルドは不要です。

## 設計書での選択

```yaml
schema_version: 1
environment_id: gpe_case0001
option_catalogs: []

select:
  source: local_framework
  destination: windows_research_runs
  model: gpe_cpu_mpi_fftw
  target: windows_gnu_msmpi
  case: gpe_quantum_taylor_green
  execution: release
  scheduler: disabled
  archive: none

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: false
```

候補の一部だけを変更するときは、同名セクションへ上書き値を書きます。
次の例では、生成先と同時コンパイル数だけが選択候補から変更されます。

```yaml
destination:
  root: D:/ResearchRuns/case0042
  case_index: D:/ResearchRuns/case_index.csv

execution:
  parallel_jobs: 16
```

実行時の並列数は、生成された`cases/<case_id>/case.yaml`で指定します。

```yaml
solver:
  mpi_processes: 8
  omp_threads: 2
```

通常は`python .\tools\run_case.py --run`でこの値を使います。一時的な上書きは
`--processes 16 --omp-threads 4`で指定できます。

Slurmではcaseの値ではなく、投入時の資源指定を優先します。

```bash
sbatch --nodes=1 --ntasks-per-node=8 --cpus-per-task=2 submit.slurm
```

`submit.slurm`は`SLURM_NTASKS`と`SLURM_CPUS_PER_TASK`を`run_case.py`へ渡します。

解決順序は「標準カタログ、追加カタログ、`select`、設計書内の上書き、
コマンドライン上書き」です。生成環境には元の`environment.source.yaml`に加えて、
最終値を記録した`environment.resolved.yaml`が保存されます。

## 候補の追加

`environment_options.local.example.yaml`を参考に、設計書と同じフォルダへ
独自カタログを作ります。

```yaml
schema_version: 1
catalog_id: laboratory_local
allow_override: false

options:
  destination:
    workstation_nvme:
      label: Workstation NVMe
      description: 高速なローカル作業領域
      values:
        root: E:/ResearchRuns/gpe_case0001
        case_index: E:/ResearchRuns/case_index.csv

  execution:
    release_large_build:
      label: Release large build
      description: 16並列でコンパイルするRelease構成
      values:
        configuration: Release
        parallel_jobs: 16
```

設計書から追加カタログを読み込み、追加したIDを選択します。相対パスは
設計書が置かれたフォルダを基準に解決されます。

```yaml
option_catalogs:
  - environment_options.local.yaml

select:
  destination: workstation_nvme
  execution: release_large_build
```

既存の分類へ新しいIDを追加するだけなら、Pythonコードの変更は不要です。
同じ分類とIDの重複は誤上書きを防ぐためエラーになります。意図的に標準候補を
置き換えるカタログだけ、`allow_override: true`を指定してください。

新しい分類そのものを追加する場合は、生成処理にも意味を持たせる必要があるため、
`environment_options.py`の`OPTION_CATEGORIES`と
`prepare_environment.py`の利用箇所を併せて拡張します。
