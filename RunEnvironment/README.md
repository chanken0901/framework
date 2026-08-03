# 外部実行環境ジェネレーター

## NSE単一GPU

NSEのCPU/MPI版と単一GPU版は、どちらも`environment.nse.yaml`を設計書の
ひな型にします。GPEと同様に、同じ設計書の`select.model`と
`parallel`でMPI、OpenMP、CUDAの使用可否を指定し、`select.execution`では
Debug/Releaseだけを選びます。

```powershell
$tool = "C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"
$design = "$designs\nse_cuda.yaml"

New-Item -ItemType Directory -Force $designs | Out-Null
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

CPU MPI + OpenMP版では次を選択します。

```yaml
select:
  model: nse_cpu_mpi
  execution: release

parallel:
  use_mpi: true
  use_openmp: true
  use_cuda: false
```

実行環境生成後は、GPEおよびNSE CPU版と同じコマンドを使います。

```powershell
python "$tool\prepare_environment.py" $design --dry-run
python "$tool\prepare_environment.py" $design

# 生成ログに表示されたnse_caseNNNNへ移動する
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

実行前に`cases\caseNNNN\case.yaml`で並列数を指定します。

```yaml
solver:
  mpi_processes: 4
  omp_threads: 4
```

## GPEとNSEを同じ手順で実行する

モデル別の標準設計書は次の2ファイルです。

| モデル | 標準設計書 | 生成先の名前 |
|---|---|---|
| GPE | `environment.gpe.yaml` | `gpe_caseNNNN` |
| NSE | `environment.nse.yaml` | `nse_caseNNNN` |

どちらも実行環境を生成した後は、同じ`tools/run_case.py`を使用します。

### GPE

```powershell
$framework = "C:\Users\Owner\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"

New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item "$tool\environment.gpe.yaml" "$designs\gpe_qtgv.yaml"
python "$tool\prepare_environment.py" "$designs\gpe_qtgv.yaml"
```

### NSE

```powershell
$framework = "C:\Users\Owner\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"

New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item "$tool\environment.nse.yaml" "$designs\nse_tgv.yaml"

python "$tool\prepare_environment.py" "$designs\nse_tgv.yaml"

# 生成ログに表示された自動採番後のディレクトリへ移動する
Set-Location "$env:USERPROFILE\ResearchRuns\nse_caseNNNN"

python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

計算条件は`cases\caseNNNN\case.yaml`を編集します。編集後の`input.dat`は
`--prepare`、`--validate-only`、`--build`、`--run`のいずれでも必要に応じて
自動再生成されます。結果は`cases\caseNNNN\output`へ出力されます。
現在のNSE MPI分割は4プロセス以上を必要とします。
`solver.mpi_processes`と`solver.omp_threads`は実行時設定であり、実行環境を
作り直さずに変更できます。

## ParaView可視化

生成した実行環境では、GPE/NSEとも同じコマンドでSLFをVTI/PVDへ変換できます。
既定では最新の完全な時刻だけを空間方向に2点おきで変換するため、まず可視化を
確認するときに向いています。

```powershell
python .\tools\postprocess_case.py
```

NSEの全解像度データについて、0、500、1000ステップだけを変換する例です。

```powershell
python .\tools\postprocess_case.py `
  --steps 0,500,1000 `
  --stride 1 `
  --fields rho,u,v,w,p
```

結果は`cases\caseNNNN\paraview`へ生成されます。ParaViewでは
`collection.pvd`を開きます。計算中に実行した場合、全MPI rankの書き込みが
完了していない最新ステップは自動的に除外されます。

ローカルFrameWorkは`SolverLibrary`と`ScriptLibrary`の編集元です。計算時は
モデル別の`environment.<model>.yaml`に従って、必要なソースとスクリプトだけをワークステーションの
ローカルSSD、またはスパコンのscratchへコピーします。NASは同期ミラーとして扱い、
直接編集、ビルド、実行には使いません。

## 処理の流れ

1. モデルとプロファイルを選択する。
2. `solver_manifest.yaml`から必要コンポーネントを解決する。
3. Fortranモジュール依存関係を検査する。
4. FrameWork外へソース、CMake、共通ビルダー、machine設定をコピーする。
5. `SetupCase/create_case_from_template.py`でcaseを生成する。
6. `case.yaml`から`input.dat`または`input.nml`を生成する。
7. FrameWork外側の共通`case_index.csv`へcase条件を登録する。
8. 搬送元、SHA-256、Gitコミットを`provenance.json`へ記録する。
9. ローカルコピーだけを使ってビルド・実行する。

## 設計書の選択肢

`source`、`destination`、`model`、`target`、`case`、`execution`、
`scheduler`、`archive`は、`environment_options.yaml`に登録されたIDで
選択できます。現在の全候補、上書き規則、候補の追加方法は
[`ENVIRONMENT_OPTIONS.md`](ENVIRONMENT_OPTIONS.md)にまとめています。

```powershell
python .\prepare_environment.py .\environment.gpe.yaml --list-options
python .\prepare_environment.py .\environment.nse.yaml --list-options model
```

設計書では次のように選びます。

```yaml
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

研究室や計算機固有の候補は別のYAMLカタログへ追加できます。既存の分類に候補を
追加する場合、Pythonコードの変更は不要です。

## ワークステーション

ローカルFrameWork上の[`environment.gpe.yaml`](environment.gpe.yaml)と
[`environment.nse.yaml`](environment.nse.yaml)をモデル別テンプレートとして保ち、
計算機側へコピーした設計書の`select.source`と`select.destination`を選びます。
通常のWindows運用では`select.source: local_framework`を使います。NASは同期ミラーとして
扱い、直接編集やビルドには使いません。候補にないパスは`source.framework_root`
または`destination.root`で上書きできます。

```powershell
$framework = "C:\Users\Owner\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"
New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item "$tool\environment.gpe.yaml" "$designs\gpe_case0001.yaml"

python "$tool\prepare_environment.py" `
  "$designs\gpe_case0001.yaml" --dry-run
python "$tool\prepare_environment.py" `
  "$designs\gpe_case0001.yaml"
```

`windows_research_runs`は`${USERPROFILE}/ResearchRuns`を親ディレクトリとして、
既存の実行環境とアーカイブを調べ、`gpe_case0001`、`gpe_case0002`のように
最初の空き番号を自動採番します。生成環境内の`case.id`も同じ番号になります。
NSEは`nse_case0001`から始まり、モデルごとに独立して採番されます。

すべてのモデルとcaseの共通台帳は、次の1ファイルです。

```text
C:\Users\Owner\ResearchRuns\case_index.csv
```

各実行環境の`cases`には個別の`case_index.csv`を置きません。共通台帳は
`case.yaml`を階層名付きの列へ展開するため、GPE・NSE・流れ場固有の条件が増えても
自動的に新しい列が追加されます。`case_key`は`gpe:case0002`のようにモデル名と
case番号を組み合わせ、モデル間で同じcase番号を使用しても区別します。

`--dry-run`はディレクトリを作成しないため、続けて実際の生成を行うと同じ番号が
使用されます。`--output`でパスを明示した場合は自動採番を行わず、そのパスを
そのまま使用します。

生成先へ移動して実行します。

```powershell
# case.yamlからinput.nml/input.datを明示的に再生成
python .\tools\run_case.py --prepare

# case.yamlの方が新しい場合は、以下のコマンドでも実行前に自動再生成
python .\tools\run_case.py --validate-only

# CMake構成とビルドだけを実行
python .\tools\run_case.py --build

# 既存の実行ファイルを使って計算だけを実行
python .\tools\run_case.py --run
```

`--run`は`case.yaml`の`solver.mpi_processes`と`solver.omp_threads`を読みます。
一度だけ別の数で実行する場合は、次のようにコマンドラインで上書きできます。

```powershell
python .\tools\run_case.py --run --processes 8 --omp-threads 2
```

`--run`はCMakeキャッシュ生成、Visual Studio初期化、configure、buildを行いません。
実行ファイルが存在しない場合は、先に`--build`を実行するようエラーで案内します。

`--all`はビルドから実行までを連続して行います。毎回ビルドしたくない場合は
`--all`ではなく`--run`を使用します。`model.include_tests: true`の場合は
CTestも含みます。実行時のカレントディレクトリは`cases/<case_id>`なので、相対指定
された`output`はcaseフォルダ内へ生成されます。

`run_case.py`は`environment.lock.json`からモデル、プロファイル、caseを取得します。
`case.yaml`が入力ファイルより新しい場合は、`case_input.py`を自動実行してから
検証、ビルド、実行へ進みます。また、`--prepare`、検証、ビルド、テスト、実行の
開始時に、現在の`case.yaml`から共通`case_index.csv`の該当行を更新します。

既存環境をまとめて再走査し、共通台帳を再構築する場合は次を使います。

```powershell
python "$tool\global_case_index.py" `
  --root "$env:USERPROFILE\ResearchRuns" `
  --rebuild
```

`case_input.py`を直接使う場合、`--profile`は`case.yaml`および
`environment.lock.json`のプロファイルと一致させる必要があります。CPU/MPI版と
CUDA版の入力が混在しないよう、不一致はエラーになります。通常は直接呼び出さず、
`python .\tools\run_case.py --prepare`を使用してください。

## モデルの変更

NSE:

```yaml
select:
  model: nse_cpu_mpi
  case: nse_case
```

NSEの流れ場はenvironment設計書では分けません。生成後の
`cases/caseNNNN/case.yaml`にある`flow.type`を変更します。

```yaml
flow:
  type: taylor_green
```

分散FFTでHIT初期条件を生成する場合は、同じ`environment.nse.yaml`で
`model: nse_cpu_mpi_2decomp_fftw`を選び、生成後の`case.yaml`を
次のように変更します。

```yaml
flow:
  type: hit_spectral
  hit:
    spectrum: johnsen
    random_seed: 13579
    rms_velocity: 0.1
    peak_wavenumber: 4.0
```

GPE:

```yaml
select:
  model: gpe_cpu_mpi_fftw
  case: gpe_quantum_taylor_green
```

GPEでは`cpu_serial_dft`、`cpu_serial_fftw`、`cpu_mpi_dft`、
`cpu_mpi_fftw`、`cuda_single`、`cuda_mpi_cufftmp`を選択できます。

## スパコン

[`environment.hpc.yaml`](environment.hpc.yaml)をローカルの設計書置場へコピーし、
ワークステーションで持ち運び用アーカイブを作ります。

```powershell
python "$tool\prepare_environment.py" "$designs\gpe_hpc.yaml" `
  --dry-run

python "$tool\prepare_environment.py" "$designs\gpe_hpc.yaml" `
  --archive --archive-format gztar
```

`windows_research_runs`が空いている`gpe_caseNNNN`を自動採番します。
生成された`tar.gz`だけをスパコンへ転送してscratchで展開します。スパコンの
コンパイラ、MPI、FFTW、CUDA、cuFFTMpの場所は、転送前または展開後に
`config/machine.yaml`へ記述します。

```bash
mkdir -p "$SCRATCH/gpe_case0001"
tar -xzf gpe_case0001.tar.gz -C "$SCRATCH/gpe_case0001" --strip-components=1
cd "$SCRATCH/gpe_case0001"
python3 tools/run_case.py --prepare
python3 tools/run_case.py --validate-only
python3 tools/run_case.py --build
sbatch --nodes=1 --ntasks-per-node=8 --cpus-per-task=2 submit.slurm
```

ビルドは利用機関が許可するログインノード、ビルドノード、またはインタラクティブジョブで
1回実行します。標準の`submit.slurm`は`python3 tools/run_case.py --run`だけを実行し、
Slurmジョブごとの再ビルドを避けます。MPIプロセス数は`SLURM_NTASKS`、OpenMP
スレッド数は`SLURM_CPUS_PER_TASK`から取得するため、`sbatch`実行時に指定します。

## 更新と再生成

生成環境は編集元ではありません。ソルバーや共通スクリプトを変更した場合はローカル
FrameWorkを更新し、Gitへ記録してから`--overwrite`で再生成します。NASへの同期は
テストとGit更新の後に行います。削除対象がこのツールの生成物であることを示す
`.generated_run_environment.json`がないディレクトリは上書きしません。
