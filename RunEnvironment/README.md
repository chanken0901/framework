# 外部実行環境ジェネレーター

ローカルFrameWorkは`SolverLibrary`と`ScriptLibrary`の編集元です。計算時は
`environment.yaml`に従って、必要なソースとスクリプトだけをワークステーションの
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
python .\prepare_environment.py .\environment.yaml --list-options
python .\prepare_environment.py .\environment.yaml --list-options model
```

設計書では次のように選びます。

```yaml
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

研究室や計算機固有の候補は別のYAMLカタログへ追加できます。既存の分類に候補を
追加する場合、Pythonコードの変更は不要です。

## ワークステーション

ローカルFrameWork上の[`environment.yaml`](environment.yaml)をテンプレートとして保ち、
計算機側へコピーした設計書の`select.source`と`select.destination`を選びます。
通常のWindows運用では`select.source: local_framework`を使います。NASは同期ミラーとして
扱い、直接編集やビルドには使いません。候補にないパスは`source.framework_root`
または`destination.root`で上書きできます。

```powershell
$framework = "C:\Users\Owner\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"
New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item "$tool\environment.yaml" "$designs\gpe_case0001.yaml"

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
  case: nse_taylor_green
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
sbatch submit.slurm
```

ビルドは利用機関が許可するログインノード、ビルドノード、またはインタラクティブジョブで
1回実行します。標準の`submit.slurm`は`python3 tools/run_case.py --run`だけを実行し、
Slurmジョブごとの再ビルドを避けます。

## 更新と再生成

生成環境は編集元ではありません。ソルバーや共通スクリプトを変更した場合はローカル
FrameWorkを更新し、Gitへ記録してから`--overwrite`で再生成します。NASへの同期は
テストとGit更新の後に行います。削除対象がこのツールの生成物であることを示す
`.generated_run_environment.json`がないディレクトリは上書きしません。
