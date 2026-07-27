# 実行環境設計書の選択肢一覧

`environment.yaml`では、各項目を値のまとまりを表すIDで選択できます。
標準候補の正本は`environment_options.yaml`です。この文書は人が確認するための
一覧であり、スクリプトが実際に読むのはYAMLカタログです。

## 一覧表示コマンド

標準候補と設計書から追加された候補をまとめて表示します。

```powershell
python .\prepare_environment.py .\environment.yaml --list-options
```

特定の分類だけを表示できます。

```powershell
python .\prepare_environment.py .\environment.yaml --list-options source
python .\prepare_environment.py .\environment.yaml --list-options model
python .\prepare_environment.py .\environment.yaml --list-options destination
```

外部ツールで処理するときはJSONで出力できます。

```powershell
python .\prepare_environment.py .\environment.yaml `
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

| ID | 初期条件 |
|---|---|
| `nse_taylor_green` | NSE Taylor-Green渦 |
| `gpe_quantum_taylor_green` | GPE量子Taylor-Green渦 |

case候補には対応モデルの条件があります。例えば、
`gpe_quantum_taylor_green`と`nse_cpu_mpi`を同時に選ぶと生成前にエラーになります。

### execution

| ID | 構成 |
|---|---|
| `debug_single` | Debug、1プロセス、1 OpenMPスレッド |
| `serial_release` | Release、1プロセス、1 OpenMPスレッド |
| `mpi4_release` | Release、4 MPIプロセス |
| `mpi8_release` | Release、8 MPIプロセス |

### scheduler

| ID | 用途 |
|---|---|
| `disabled` | ジョブ投入スクリプトを生成しない |
| `slurm_cpu_8` | 1ノード、8 CPUタスクのSlurmスクリプト |

### archive

| ID | 用途 |
|---|---|
| `none` | アーカイブを作成しない |
| `zip` | ZIPを作成 |
| `gztar` | tar.gzを作成 |

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
  execution: mpi4_release
  scheduler: disabled
  archive: none
```

候補の一部だけを変更するときは、同名セクションへ上書き値を書きます。
次の例では、生成先とMPIプロセス数だけが選択候補から変更されます。

```yaml
destination:
  root: D:/ResearchRuns/case0042
  case_index: D:/ResearchRuns/case_index.csv

execution:
  processes: 6
```

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
    mpi16_release:
      label: Release sixteen MPI processes
      description: 16 MPIプロセスで実行
      values:
        configuration: Release
        processes: 16
        omp_threads: 1
        parallel_jobs: 16
```

設計書から追加カタログを読み込み、追加したIDを選択します。相対パスは
設計書が置かれたフォルダを基準に解決されます。

```yaml
option_catalogs:
  - environment_options.local.yaml

select:
  destination: workstation_nvme
  execution: mpi16_release
```

既存の分類へ新しいIDを追加するだけなら、Pythonコードの変更は不要です。
同じ分類とIDの重複は誤上書きを防ぐためエラーになります。意図的に標準候補を
置き換えるカタログだけ、`allow_override: true`を指定してください。

新しい分類そのものを追加する場合は、生成処理にも意味を持たせる必要があるため、
`environment_options.py`の`OPTION_CATEGORIES`と
`prepare_environment.py`の利用箇所を併せて拡張します。
