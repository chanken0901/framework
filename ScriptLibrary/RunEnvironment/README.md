# 外部実行環境ジェネレーター

Stage 8の多成分MPIペンシル/OpenMP環境は`environment.nse_multicomponent.parallel.yaml`、
OpenMP単独は`environment.nse_multicomponent.openmp.yaml`を使用する。
設定とWindows/Linux手順は[多成分並列計算](../../SolverLibrary/NSE/docs/NSE_MULTICOMPONENT_PARALLEL.md)を参照。

ケース設定の分割形式、従来形式との互換性、`resolved_case.yaml`については
[`CASE_CONFIGURATION.md`](CASE_CONFIGURATION.md)を参照してください。新しい多成分ケースは
必要な拡張だけを`config/*.yaml`として生成します。

多成分・反応流拡張のStage 0検証環境には
`environment.nse_multicomponent.yaml`を使用します。この環境は現行単成分NSEとは
別モデル`nse_multicomponent`を生成し、まだ流体時間発展を行いません。

Stage 1の周期パッシブスカラー移流には
`environment.nse_multicomponent.passive_scalar.yaml`を使用します。一定速度場、
一次風上法、SSPRK3で部分密度を保存形式により移流し、最終場をCSV出力します。
Stage 1はCPU逐次専用で、運動量・全エネルギーは更新しません。

Stage 2の非反応・非粘性多成分Euler計算には
`environment.nse_multicomponent.inviscid.yaml`を使用します。全species部分密度、
運動量、全エネルギーをRusanov流束とSSPRK3で連成更新します。物理拡散、粘性、
熱伝導、化学反応はまだ含みません。

Stage 3の温度・組成依存熱力学には
`environment.nse_multicomponent.thermally_perfect.yaml`を使用します。化学種ごとの
分子量とNASA-7係数から温度、圧力、比熱比、音速を計算し、最終CSVへ温度`T`も
出力します。物性値と流れ場は整合するSI単位で指定してください。Stage 3もCPU逐次、
非粘性、周期境界、非反応です。

Stage 4の混合平均輸送には
`environment.nse_multicomponent.viscous.yaml`を使用します。設計書から
`case_templates/nse_multicomponent_viscous.yaml`を展開し、一定のspecies拡散係数、
基準粘性係数、Prandtl数と周期species波を設定します。生成される
`cpu_serial_viscous`環境はCPU逐次で、拡散・粘性・熱伝導を含む非反応
多成分Navier--Stokes計算をビルド・試験できます。

Stage 5の0次元有限反応速度化学には
`environment.nse_multicomponent.reactor.yaml`を使用します。設計書から反応器の
`case.yaml`と、多成分・NASA-7熱力学・一段Arrhenius反応の3つの拡張YAMLを生成します。
`cpu_serial_reactor`環境は断熱・定容の均質反応器をSSPRK3で積分し、時刻、温度、
圧力、全species質量分率をCSVへ出力します。Stage 5はCPU逐次の化学反応単体検証で、
流体輸送を行いません。

Stage 6の反応性多成分Navier--Stokes計算には
`environment.nse_multicomponent.reactive.yaml`を使用します。設計書から周期species波の
`case.yaml`と、多成分・NASA-7熱力学・混合平均輸送・一段Arrhenius反応の4つの拡張YAMLを
生成します。`cpu_serial_reactive`環境は、Stage 4の対流・species拡散・Newton粘性・
Fourier熱伝導とStage 5の化学反応をStrang分割で結合します。全体時間刻みは対流・拡散・
化学制約から決まり、化学半stepは必要に応じてSSPRK3でsubcycleします。現在はCPU逐次、
直交等間隔格子、全方向周期境界、一次Rusanov流束に限定しています。

Stage 7の反応衝撃波管と面別物理境界には
`environment.nse_multicomponent.reactive_boundaries.yaml`を使用します。設計書から
`case_templates/nse_multicomponent_reactive_shock_tube.yaml`を展開し、x方向の二状態反応流、
6面別の`periodic`／`reflective`／`dirichlet`／`non_reflecting`、境界参照状態、
途中スナップショットと積分履歴を設定します。solver profileは
`cpu_serial_reactive_boundaries`です。現在はCPU逐次、直交等間隔格子、一次Rusanov流束、
Strang分割に限定されます。境界と出力の詳細は
[`../../SolverLibrary/NSE/docs/NSE_MULTICOMPONENT_REACTIVE_BOUNDARIES.md`](../../SolverLibrary/NSE/docs/NSE_MULTICOMPONENT_REACTIVE_BOUNDARIES.md)
を参照してください。

> OS固有の操作はWindows（PowerShell）とLinux（bash）を併記します。共通のパス、Python、CMake、MPIの対応表は[`../../docs/WINDOWS_LINUX_COMMANDS.md`](../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

NSEの`case.yaml`では`numerics.convective_scheme`に`keep2`、`keep6`、
`weno5z_roe`、`hybrid`を指定できる。ハイブリッド構成例は次のとおり。

```yaml
numerics:
  convective_scheme: hybrid
  hybrid:
    smooth_scheme: keep6
    shock_scheme: weno5z_roe
    sensor: ducros_pressure
    sensor_onset: 0.01
    sensor_full: 0.10
```

非ハイブリッド計算では`convective_scheme`を`keep2`、`keep6`、または
`weno5z_roe`にする。`hybrid`以下の項目はその場合使用されない。
対流流束の選択項目は`numerics.convective_scheme`だけであり、旧
`numerics.flux`と未使用だった`numerics.reconstruction`は受け付けない。
旧項目が残っている場合は、設定を黙って無視せず移行エラーを表示する。

## NSE CUDA（単一GPU／MPI＋マルチGPU）

NSEのCPU/MPI版、単一GPU版、MPI＋CUDA版は、いずれも`environment.nse.yaml`を設計書の
ひな型にします。`select.model`は物理モデルの`nse`だけを指定し、
`parallel`でMPI、OpenMP、CUDAの使用可否を指定します。solver profileは
この組合せと`case.yaml`が要求するFFT機能から自動選択され、
`select.execution`ではDebug/Releaseだけを選びます。

```powershell
$tool = "C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"
$design = "$designs\nse_cuda.yaml"

New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item "$tool\environment.nse.yaml" $design
```

Linux（bash）:

```bash
tool="$HOME/Research/FrameWork/ScriptLibrary/RunEnvironment"
designs="$HOME/ResearchRuns/Designs"
design="$designs/nse_cuda.yaml"

mkdir -p "$designs"
cp "$tool/environment.nse.yaml" "$design"
```

Linuxでは、下記のmodel／parallel設定に加えて`select.source`、
`select.destination`、`select.target`をそれぞれLinux用へ変更します。具体例は
[`../../docs/WINDOWS_LINUX_COMMANDS.md`](../../docs/WINDOWS_LINUX_COMMANDS.md)の第6章を参照してください。

コピーした設計書を次のように変更します。

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

CPU MPI + OpenMP版では次を選択します。

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: true
  use_openmp: true
  use_cuda: false
```

MPI＋CUDAマルチGPU版では次を選択します。

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: true
```

この組合せでは通常の時間発展に`cuda_mpi` profileが選択されます。
`solver.mpi_processes`をGPU総数にし、通常は1 MPI rankを1 GPUへ割り当てます。
Taylor–Greenと保存済み乱流場は`cuda_mpi`、分散HIT初期化または
Petersen–Livescu forcingはLinux用`cuda_mpi_cufftmp`を使用します。

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

Linux（bash）:

```bash
python3 "$tool/prepare_environment.py" "$design" --dry-run
python3 "$tool/prepare_environment.py" "$design"

# 生成ログに表示されたnse_caseNNNNへ移動する
cd "$HOME/ResearchRuns/nse_caseNNNN"
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --run
```

実行前に`cases\caseNNNN\case.yaml`で並列数を指定します。

```yaml
solver:
  mpi_processes: 4
  omp_threads: 4
```

CPU逐次・MPI・OpenMP・単一GPU CUDA・MPI＋CUDA版では、`case.yaml`の対流流束を次の4種類から選べます。

```yaml
numerics:
  # keep2, keep6, weno5z_roe, hybrid
  convective_scheme: weno5z_roe
```

`weno5z_roe`は特性空間の5次精度WENO-Z再構築とRoe流束です。`hybrid`は
滑らかな領域のKEEPと衝撃波領域のWENO5-Z/Roeをセンサーで連続的に混合します。
どちらも単一GPU CUDA版とMPI＋CUDA版に対応しており、`cuda_single`／`cuda_mpi`プロファイルでも同じ
`case.yaml`設定を使用できます。

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

Linux（bash）:

```bash
framework="$HOME/Research/FrameWork"
tool="$framework/ScriptLibrary/RunEnvironment"
designs="$HOME/ResearchRuns/Designs"

mkdir -p "$designs"
cp "$tool/environment.gpe.yaml" "$designs/gpe_qtgv.yaml"
# コピー後、source/destination/targetをLinux用へ変更する
python3 "$tool/prepare_environment.py" "$designs/gpe_qtgv.yaml"
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

Linux（bash）:

```bash
framework="$HOME/Research/FrameWork"
tool="$framework/ScriptLibrary/RunEnvironment"
designs="$HOME/ResearchRuns/Designs"

mkdir -p "$designs"
cp "$tool/environment.nse.yaml" "$designs/nse_tgv.yaml"
# コピー後、source/destination/targetをLinux用へ変更する
python3 "$tool/prepare_environment.py" "$designs/nse_tgv.yaml"

# 生成ログに表示された自動採番後のディレクトリへ移動する
cd "$HOME/ResearchRuns/nse_caseNNNN"
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --run
```

計算条件は`cases\caseNNNN\case.yaml`を編集します。編集後の`input.dat`は
`--prepare`、`--validate-only`、`--build`、`--run`のいずれでも必要に応じて
自動再生成されます。結果は`cases\caseNNNN\output`へ出力されます。
現在のNSE MPI分割は2プロセス以上を必要とします。MPI＋CUDA版ではさらに
各y/z局所ブロックが`nghost`セル以上になる格子数とrank数を選びます。
`solver.mpi_processes`と`solver.omp_threads`は実行時設定であり、実行環境を
作り直さずに変更できます。

## ParaView可視化

生成した実行環境では、GPE/NSEとも`tools\postprocess_case.py`を共通入口として
SLFをVTI/PVDへ変換します。以下のコマンドは、`environment.lock.json`がある
**実行環境ルート**で実行してください。

まず、後処理ツールが現行版か確認します。

```powershell
Test-Path .\environment.lock.json
python .\tools\postprocess_case.py --version
```

`True`と`postprocess_case.py 2.0.0`以上が表示されれば、作業場所と版は正しい状態です。
次に、実行される子コマンドだけを表示する場合は`--dry-run`を使用します。この段階では
SLFを読み込まず、VTI/PVDも生成しません。

```powershell
python .\tools\postprocess_case.py --dry-run
```

入力ファイル、`meta.json`、選択時刻、変数、格子数を検査し、まだ変換しない場合は
`--inspect-only`を使用します。

```powershell
python .\tools\postprocess_case.py --inspect-only --steps latest
```

既定の変換は、最新の完全な時刻だけを空間方向に2点おきで処理します。最初の
可視化確認向けです。

```powershell
python .\tools\postprocess_case.py
```

全解像度の0、500、1000ステップだけを変換する例です。NSEでは次のようにします。

```powershell
python .\tools\postprocess_case.py `
  --steps 0,500,1000 `
  --stride 1 `
  --fields rho,u,v,w,p
```

GPEの全解像度では`--fields density,phase`を指定します。`--steps`には`all`、
`latest`、`0,500,1000`、または終端を含む`0:1000:100`を指定できます。

結果は`cases\caseNNNN\paraview`へ生成されます。ParaViewでは
`collection.pvd`を開きます。計算中に実行した場合、全MPI rankの書き込みが
完了していない最新ステップは自動的に除外されます。

生成済み実行環境はFrameWorkの更新を自動追従しません。`--version`が使えない、
2.0.0未満、または`unrecognized arguments: --steps`が出る場合は、現在の
FrameWorkから実行環境を再生成してください。詳細はFrameWorkルートの
`後処理_ParaView可視化手順書.md`または同名のWord版を参照します。

## NSE乱流統計の後処理

保存済みSLFから全時刻の乱流統計を計算する場合は、次を実行します。

```powershell
python .\tools\postprocess_case.py --task statistics
```

結果は`cases\caseNNNN\statistics\turbulence_statistics.csv`へ時系列で出力され、
各統計量の定義は`turbulence_statistics_metadata.json`へ保存されます。
方向別速度変動RMS、積分スケール、Taylorマイクロスケール、Kolmogorovスケール、
`Re_L`、`Re_lambda`、乱流Mach数、散逸率、Reynolds応力と等方性誤差を含みます。
`gamma`と基準Reynolds数は`case.yaml`の`physics.nse`から自動取得します。

ParaView変換と乱流統計を続けて実行する場合は次の形です。

```powershell
python .\tools\postprocess_case.py --task all
```

`--steps 0,500,1000`や`--steps 0:1000:100`を付けると、統計処理する時刻を
限定できます。統計計算ではFFTを行うため空間間引きは使用せず、SLFの全格子を
処理します。

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
  model: gpe
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
`--all`ではなく`--run`を使用します。`solver.include_tests: true`の場合は
CTestも含みます。実行時のカレントディレクトリは`cases/<case_id>`なので、相対指定
された`output`はcaseフォルダ内へ生成されます。

`run_case.py`は`environment.lock.json`からモデル、プロファイル、caseを取得します。
`case.yaml`が入力ファイルより新しい場合は、`case_input.py`を自動実行してから
検証、ビルド、実行へ進みます。また、`--prepare`、検証、ビルド、テスト、実行の
開始時に、現在の`case.yaml`から共通`case_index.csv`の該当行を更新します。
Windowsで`case_index.csv`をExcelやエディタが開いていて置換できない場合は、警告を
表示して入力生成、ビルド、実行を継続します。共通台帳は`case.yaml`から再構築できる
補助ファイルであり、ファイルを閉じた後の次回実行時に自動で再同期されます。

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
  model: nse
  case: nse_case

parallel:
  use_mpi: true
  use_openmp: true
  use_cuda: false
```

NSEの流れ場はenvironment設計書では分けません。生成後の
`cases/caseNNNN/case.yaml`にある`flow.type`を変更します。

```yaml
flow:
  type: taylor_green
```

分散FFTでHIT初期条件を生成する場合も、environment設計書へ特殊profileを
追加する必要はありません。CPU/MPI環境には通常版と2DECOMP&FFT版が同梱され、
`run_case.py`が`case.yaml`から必要なprofileを自動選択します。生成後の
`case.yaml`だけを次のように変更します。

```yaml
flow:
  type: hit
  hit:
    turbulent_mach_number: 0.5
    turbulent_reynolds_number: 30.0
    random_seed: 13579
    dealias_fraction: 0.6666666666666667
    isotropy_mode: projected_shell
    isotropy_k_cutoff: 2.5
    isotropy_tolerance: 1.0e-8
    isotropy_max_iterations: 80
    spectrum:
      type: pope  # johnsen または pope
      johnsen:
        characteristic_length: 1.0
        length_scale_ratio: 2.0
      pope:
        integral_length: 1.0
        energy_constant: 1.5
        large_scale_constant: 6.78
        dissipation_constant: 0.40
        large_scale_exponent: 2.0
        dissipation_exponent: 5.2
```

`turbulent_reynolds_number`は元コードの`Re_ini`に対応する
テイラー・マイクロスケールReynolds数`Re_lambda`です。ソルバーは
選択したスペクトルの代表長さ`L`とともに、次の関係から初期条件と
粘性係数を一貫して決めます。

```text
u'       = M_t / sqrt(3)
Re_L     = 3 Re_lambda^2 / 20
nu       = sqrt(3/2) u' L / Re_L
Re_solver = 1 / nu
lambda   = L sqrt(10 / Re_L)
eta      = L Re_L^(-3/4)
```

初期速度場は逆FFT後に成分RMSが`u'`へ一致するよう正規化されます。
HITでは`physics.nse.mach_number`と`reynolds_number`より、この導出値が
優先されます。`prandtl_number`は独立に指定し、熱拡散係数へ
`1/(Re_solver Pr)`の組合せで反映されます。

Johnsenでは`length_scale_ratio`が元コードの`L_lambda`に対応し、
`k_peak = 2 length_scale_ratio / characteristic_length`です。Popeでは
`K_C`、`C_L`、`C_eta`、`p_0`、`beta`に対応する値を`pope`ブロックだけで
指定します。`spectrum.type`で選ばれていないブロックは使用されません。

`projected_shell`は低波数シェルの二次統計を等方化し、Helmholtz射影で
発散ゼロ条件を維持します。`isotropy_k_cutoff`は、低波数Forcingを使う場合は
通常`forcing.petersen_livescu.k_cutoff`と同じ値にします。

保存済みNSE乱流場を長いx領域へ配置する場合は、元ケースのrank別SLFを
ghostなしの可搬SLFへ変換します。実行環境ルートで次を実行します。

元計算と読込み先はCPU MPI／CUDAのどちらでもよく、CPU MPI→CPU MPI、
CPU MPI→CUDA、CUDA→CPU MPI、CUDA→CUDAのすべてに対応します。元計算と
読込み先のMPIプロセス数を一致させる必要はありません。

```powershell
python .\SolverLibrary\NSE\tools\nse_prepare_imported_turbulence.py `
  .\previous_case\output `
  --step latest `
  --output .\cases\caseNNNN\initial_data\turbulence.slf
```

`case.yaml`では次のように指定します。

```yaml
flow:
  type: imported_turbulence
  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: embed  # embed または tile
    x_start: 2.0
    blend_cells: 8
    velocity_offset: [0.5, 0.0, 0.0]
    background:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143
```

`embed`は一つの乱流ブロックを背景場へ滑らかに接続し、`tile`はx方向へ
周期反復します。CPU/CUDA別のコピー可能なPowerShell手順、格子互換条件、
正常終了の確認方法は
[`NSE_IMPORTED_TURBULENCE.md`](../../SolverLibrary/NSE/docs/NSE_IMPORTED_TURBULENCE.md)
を参照してください。

CPU/MPI/OpenMP、単一GPU CUDA、MPI＋CUDAの各版で、各物理面を`periodic`、
`non_reflecting`、`reflective`から選べます。無反射面には`reference_state`を指定し、
`boundary.reference_states`に密度、速度、圧力を定義します。鏡像面は基準状態を持たず、
法線運動量だけを反転する自由滑り・断熱条件です。詳細は
[`NSE_BOUNDARY_CONDITIONS.md`](../../SolverLibrary/NSE/docs/NSE_BOUNDARY_CONDITIONS.md)
を参照してください。

Forcingも流れ場と同様に`type`で選択し、方式固有の設定を同名のブロックへ
まとめます。無効化するときは`type: none`のままにします。

```yaml
forcing:
  type: petersen_livescu
  petersen_livescu:
    spectrum: low_wavenumber
    fft_backend: auto
    k_cutoff: 2.5
    target_dissipation: 0.1
    dilatational_ratio: 0.0
    denominator_floor: 1.0e-14
    max_coefficient: 0.0
    report_interval: 100
```

未登録の`type`、方式固有ブロック内の未知のキー、新旧形式の混在は
`python tools/run_case.py --prepare`でエラーになります。

GPE:

```yaml
select:
  model: gpe
  case: gpe_quantum_taylor_green

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: false
  fft_decomposition: slab
```

GPEでは`parallel.use_mpi`と`parallel.use_cuda`から、逐次CPU、MPI CPU、
単一GPU、複数GPUを選択します。CPU版は通常FFTW profileを自動選択します。
MPI CPU版とMPI/cuFFTMp版では`parallel.fft_decomposition`に`slab`（既定）または
`pencil`を指定できます。CPUで`pencil`を選ぶと`cpu_mpi_pencil_fftw`、
CUDA併用で選ぶと`cuda_mpi_cufftmp_pencil`が自動選択されます。
後者にはcuFFTMp 11.4.0（NVIDIA HPC SDK 25.3）以降が必要です。
検証用DFTが必要な場合だけ、`solver.profile`へ`cpu_serial_dft`または
`cpu_mpi_dft`を明示します。

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

`windows_research_runs`のように自動採番が有効な設計では、上書きするcase番号を
必ず`--case-id`で指定します。対象指定のない`--overwrite`は次の空き番号を
生成せず、誤操作防止のエラーになります。

```powershell
$tool = "C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\RunEnvironment"
$design = "$env:USERPROFILE\ResearchRuns\Designs\nse.yaml"

# case0015を置換する計画だけを確認
python "$tool\prepare_environment.py" $design `
  --case-id case0015 --overwrite --dry-run

# nse_case0015を削除して、現在のFrameWorkから再生成
python "$tool\prepare_environment.py" $design `
  --case-id case0015 --overwrite
```

上書きでは`build`、`cases`、計算結果を含む生成環境全体を置換します。残す必要がある
`case.yaml`や`output`は実行前に別の場所へ退避してください。生成先を完全なパスで
指定する従来の`--output <path> --overwrite`も引き続き利用できます。

## NSEの平面衝撃波–乱流干渉

保存乱流の外側へ平面衝撃波を配置し、衝撃波背後状態をDirichlet境界から供給する
NSEケースは`flow.type: shock_turbulence_interaction`で生成できます。設計書の完全な
YAML例、衝撃波前後状態の定義、配置制約、実行手順は
[`NSE_SHOCK_TURBULENCE_INTERACTION.md`](../../SolverLibrary/NSE/docs/NSE_SHOCK_TURBULENCE_INTERACTION.md)
を参照してください。`case.yaml`変更後は`python .\tools\run_case.py --prepare`を再実行します。

## NSEの有限高圧室–衝撃波–乱流干渉

鏡像閉端を持つ有限高圧室から衝撃波と膨張波を発生させ、局所乱流へ入射させる場合は
`flow.type: shock_tube_turbulence_interaction`を使用します。高圧・低圧状態、隔膜位置、
境界条件、実行手順は
[`NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md`](../../SolverLibrary/NSE/docs/NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)
を参照してください。
Stage 9の平面ノズルは `environment.nse_multicomponent.nozzle.yaml` を使用します。
形状は `config/geometry.yaml` に分離されます。
手順と制約は [一般座標手順書](../../SolverLibrary/NSE/docs/NSE_MULTICOMPONENT_GEOMETRY.md) を参照してください。
