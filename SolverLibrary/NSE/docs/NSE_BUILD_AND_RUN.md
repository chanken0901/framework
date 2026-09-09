# NSE 生成・ビルド・実行手順

> OS固有の操作はWindows（PowerShell）とLinux（bash）を併記します。共通のパス、Python、CMake、MPIの対応表は[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

更新日: 2026-09-02

## 1. 推奨フロー

通常の計算では、ローカル`FrameWork`を編集元とし、
`ScriptLibrary/RunEnvironment`から`ResearchRuns`へ外部実行環境を生成する。
生成環境には選択したソース、ビルド設計、マシン設定、case、実行ツールが入る。

```text
FrameWork (edit source)
  -> prepare_environment.py
ResearchRuns/nse_caseNNNN (build and run)
  -> run_case.py --prepare
  -> run_case.py --validate-only
  -> run_case.py --build
  -> run_case.py --run
```

NASは同期ミラーとし、NAS上で直接ビルドまたは実行しない。

## 2. 前提ソフトウェア

Windows CPU/MPI版:

- Python 3
- CMake 3.24以上
- Ninja
- GNU Fortran
- Microsoft MPI RuntimeとSDK

Windows CUDA版（単一GPU／MPI＋CUDA）では、上記にCUDA Toolkit、NVIDIAドライバ、
Visual Studio Build Toolsのx64 C++環境が必要となる。

```powershell
python --version
cmake --version
ninja --version
gfortran --version
mpiexec -help
nvcc --version
```

Linux/スパコンでは、同等のPython、CMake、Fortranコンパイラ、
MPI、選択したFFT/CUDAバックエンドをmoduleなどで読み込む。

```bash
python3 --version
cmake --version
ninja --version
gfortran --version
mpirun --version
nvcc --version
```

## 3. 新規clone

`2decomp-fft`はGitサブモジュールである。FrameWorkモノレポの新規取得時は次を使う。

Windows（PowerShell）:

```powershell
git clone --recurse-submodules https://github.com/chanken0901/framework.git FrameWork
```

Linux（bash）:

```bash
git clone --recurse-submodules https://github.com/chanken0901/framework.git FrameWork
```

既存cloneへ取得する場合:

```console
git submodule update --init --recursive
```

このサブモジュール更新コマンドはWindows／Linux共通である。

## 4. 外部実行環境を生成する

PowerShellで設計書のコピーを作る。

```powershell
$framework = "$env:USERPROFILE\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"
$design = "$designs\nse.yaml"

New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item "$tool\environment.nse.yaml" $design
```

Linux（bash）では次のようにコピーする。

```bash
framework="$HOME/Research/FrameWork"
tool="$framework/ScriptLibrary/RunEnvironment"
designs="$HOME/ResearchRuns/Designs"
design="$designs/nse.yaml"

mkdir -p "$designs"
cp "$tool/environment.nse.yaml" "$design"
```

Linuxではコピーした設計書の`source.framework_root`、`destination`、
`select.target`をLinux用へ変更する。具体例は
[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)の第6章を参照する。

CPU MPI + OpenMP:

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: true
  use_openmp: true
  use_cuda: false
```

CPU MPIでHIT初期化またはPetersen-Livescu forcingを選ぶ場合も、上と同じ
environment設計を使う。互換profileは生成時に同梱され、`case.yaml`の内容に応じて
`run_case.py`が2DECOMP&FFT版へ自動切替する。

単一GPU CUDA:

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

MPI＋CUDAマルチGPU:

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: true
```

この組合せでは通常の時間発展に`cuda_mpi` profileが選択される。`case.yaml`の
`solver.mpi_processes`は使用するGPU総数にし、原則として1 MPI rankを1 GPUへ
割り当てる。Taylor–Greenと保存済み乱流場は`cuda_mpi`で実行する。
分散HIT初期化またはPetersen–Livescu forcingにはLinux上で
`cuda_mpi_cufftmp`を選択する。cuFFTMp版のビルド方法は
[`NSE_CUFFTMP.md`](NSE_CUFFTMP.md)を参照する。

生成前に選択肢と展開内容を確認する。

```powershell
python "$tool\prepare_environment.py" $design --list-options model
python "$tool\prepare_environment.py" $design --dry-run
python "$tool\prepare_environment.py" $design
```

Linux（bash）:

```bash
python3 "$tool/prepare_environment.py" "$design" --list-options model
python3 "$tool/prepare_environment.py" "$design" --dry-run
python3 "$tool/prepare_environment.py" "$design"
```

生成ログに表示された`nse_caseNNNN`へ移動する。

```powershell
Set-Location "$env:USERPROFILE\ResearchRuns\nse_caseNNNN"
```

Linux（bash）:

```bash
cd "$HOME/ResearchRuns/nse_caseNNNN"
```

## 5. case設計書

計算条件は`cases\caseNNNN\case.yaml`で管理する。`input.dat`を正本として
直接編集しない。主な項目は次のとおり。

- `grid`: 格子数、領域、ghostセル数
- `time`: CFL、時間刻み、最大ステップ、出力間隔
- `physics.nse`: 比熱比、Mach数、Reynolds数、Prandtl数
- `flow`: Taylor-Green、HIT、または保存済み乱流場の初期条件
- `forcing`: Forcingの有無と方式
- `boundary`: 6物理面の周期／無反射／鏡像条件と無反射基準状態
- `numerics`: 対流流束、粘性項、時間積分
- `solver`: MPIランク数、OpenMPスレッド数、CUDAデバイス

対流流束は実行時に4方式から選択できる。

```yaml
numerics:
  # keep2, keep6, weno5z_roe, hybrid
  convective_scheme: hybrid
  hybrid:
    smooth_scheme: keep6
    shock_scheme: weno5z_roe
    sensor: ducros_pressure
    sensor_onset: 0.01
    sensor_full: 0.10
  viscous_scheme: central6
  time_integrator: ssprk3
```

`hybrid`以外の場合、`hybrid`以下の項目は使用されない。設定と混合則は
[`NSE_HYBRID_FLUX.md`](NSE_HYBRID_FLUX.md)を参照する。

CPU版でx方向だけを無反射、y-z方向を周期にする例を次に示す。

```yaml
boundary:
  faces:
    x_min: {type: non_reflecting, reference_state: far_field}
    x_max: {type: non_reflecting, reference_state: far_field}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
  reference_states:
    far_field:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143
  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto
```

周期面は方向ごとの対で指定する。この面別設定はCPU/MPI/OpenMP、単一GPU、MPI＋CUDAの
全profileで共通に使用できる。詳細は
[`NSE_BOUNDARY_CONDITIONS.md`](NSE_BOUNDARY_CONDITIONS.md)を参照する。

自由滑り・断熱の鏡像面にする場合は、基準状態を付けずに`reflective`を指定する。

```yaml
boundary:
  faces:
    x_min: {type: reflective}
    x_max: {type: reflective}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
  reference_states: {}
```

`reflective`では密度、接線運動量、全エネルギーを偶対称、法線運動量を奇対称に
ghostへ写す。粘着壁（no-slip）や規定温度壁が必要な場合は別の境界条件を実装する。

保存済み乱流を長いx領域へ配置する場合は、先にrank別SLFを可搬SLFへ変換する。
元計算と読込み先の組合せは、CPU MPI→CPU MPI、CPU MPI→CUDA、CUDA→CPU MPI、
CUDA→CUDAのすべてに対応する。元計算と読込み先でMPIプロセス数を一致させる
必要はない。

```powershell
python .\SolverLibrary\NSE\tools\nse_prepare_imported_turbulence.py `
  .\previous_case\output `
  --step latest `
  --output .\cases\caseNNNN\initial_data\turbulence.slf
```

Linux（bash）:

```bash
python3 ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py \
  ./previous_case/output \
  --step latest \
  --output ./cases/caseNNNN/initial_data/turbulence.slf
```

続いて`case.yaml`の`flow.type`を`imported_turbulence`へ変更する。`embed`、`tile`、`periodic_embed`の
入力例、CPU/CUDA別のWindows／Linux実行手順、格子互換条件、平均速度の追加方法は
[`NSE_IMPORTED_TURBULENCE.md`](NSE_IMPORTED_TURBULENCE.md)を参照する。

全領域でWENO5-Z/Roeを使う場合は`convective_scheme: weno5z_roe`、全領域で
KEEP6を使う場合は`convective_scheme: keep6`とする。後者ではWENO計算を行わず、
生成後の`input.dat`にも`hybrid_*`は出力されない。対流方式の選択項目は
`numerics.convective_scheme`だけであり、旧`numerics.flux`と未使用だった
`numerics.reconstruction`が残っている場合は移行エラーとなる。

## 6. 入力生成から実行まで

```powershell
# case.yamlからinput.datを生成
python .\tools\run_case.py --prepare

# 設計書、プロファイル、バックエンドの組合せを検証
python .\tools\run_case.py --validate-only

# CMake configureとビルド
python .\tools\run_case.py --build

# 既存実行ファイルで計算
python .\tools\run_case.py --run
```

Linux（bash）:

```bash
# case.yamlからinput.datを生成
python3 ./tools/run_case.py --prepare

# 設計書、プロファイル、バックエンドの組合せを検証
python3 ./tools/run_case.py --validate-only

# CMake configureとビルド
python3 ./tools/run_case.py --build

# 既存実行ファイルで計算
python3 ./tools/run_case.py --run
```

一時的にMPIランク数とOpenMPスレッド数を上書きする場合:

```powershell
python .\tools\run_case.py --run --processes 8 --omp-threads 2
```

```bash
python3 ./tools/run_case.py --run --processes 8 --omp-threads 2
```

テストを含む実行環境では次を使う。

```powershell
python .\tools\run_case.py --test
```

```bash
python3 ./tools/run_case.py --test
```

`--run`はconfigureやbuildを行わない。再ビルドが必要な場合は先に`--build`を実行する。

## 7. 終了と出力の確認

出力は`cases\caseNNNN\output`へ書き込まれる。最終保存ステップのファイル出力後、
標準出力に次の完了メッセージが表示されれば正常終了である。

```text
NSE calculation completed successfully: step=<final step>, time=<final time>
```

`IEEE_UNDERFLOW_FLAG`や`IEEE_DENORMAL`は、単独では発散を意味しない。ただし現行版は
正常終了時に上記メッセージを出すため、正常性は完了メッセージ、最終出力、
ステップログの有限値を合わせて判定する。

ParaView変換:

```powershell
python .\tools\postprocess_case.py
```

```bash
python3 ./tools/postprocess_case.py
```

最新の完全な保存ステップが変換され、
`cases\caseNNNN\paraview\collection.pvd`が生成される。

## 8. SolverLibrary上で直接ビルドする

開発時のテストでは共通ランナー`ScriptLibrary/BuildSolver`を使う。

```powershell
Set-Location "$env:USERPROFILE\Documents\Codex\FrameWork\ScriptLibrary\BuildSolver"

python .\build_model.py .\build.yaml --model nse --list-profiles
python .\build_model.py .\build.yaml `
  --model nse --profile cpu_mpi --test
python .\build_model.py .\build.yaml `
  --model nse --profile cuda_single --test
python .\build_model.py .\build.yaml `
  --model nse --profile cuda_mpi --build
```

Linux（bash）:

```bash
cd "$HOME/Research/FrameWork/ScriptLibrary/BuildSolver"

python3 ./build_model.py ./build.yaml --model nse --list-profiles
python3 ./build_model.py ./build.yaml \
  --model nse --profile cpu_mpi --test
python3 ./build_model.py ./build.yaml \
  --model nse --profile cuda_single --test
python3 ./build_model.py ./build.yaml \
  --model nse --profile cuda_mpi --build
```

主なプロファイル:

| profile | 用途 |
|---|---|
| `cpu_mpi` | MPI + OpenMP、Taylor-Greenなど |
| `cpu_mpi_2decomp_fftw` | 2DECOMP&FFTを使うHIT初期化・Forcing |
| `cuda_single` | 単一GPU CUDA |
| `cuda_mpi` | MPI＋CUDA、1 rank＝1 GPUのマルチGPU時間発展 |

ビルド生成物は既定で次へ置く。

```text
%LOCALAPPDATA%\SolverLibraryBuild\nse-<profile>-<configuration>
```

## 9. 互換・低レベル手順

NSE単体の`config/build.yaml`を読む互換ランナー:

```powershell
Set-Location "$env:USERPROFILE\Documents\Codex\FrameWork\SolverLibrary\NSE"
python .\tools\build_from_yaml.py .\config\build.yaml --validate-only
python .\tools\build_from_yaml.py .\config\build.yaml
```

Linux（bash）:

```bash
cd "$HOME/Research/FrameWork/SolverLibrary/NSE"
python3 ./tools/build_from_yaml.py ./config/build.yaml --validate-only
python3 ./tools/build_from_yaml.py ./config/build.yaml
```

CMakeを直接操作する必要がある場合:

```powershell
cmake --preset windows-msmpi-release
cmake --build --preset windows-msmpi-release
```

Linux:

```bash
cmake --preset linux-mpi-release -DCMAKE_Fortran_COMPILER=mpifort
cmake --build --preset linux-mpi-release
```

通常はこれらの低レベル手順より、外部実行環境または`BuildSolver`を使う。
生成された`InitialCache.cmake`や`resolved_build.json`は直接編集しない。

## 10. 更新時の検証

変更内容に応じて次を実行する。

```powershell
# CPU/MPI全テスト
python .\build_model.py .\build.yaml `
  --model nse --profile cpu_mpi --test

# 単一GPU全テスト
python .\build_model.py .\build.yaml `
  --model nse --profile cuda_single --test

# MPI＋CUDAビルド（実機smoke testは複数GPUノードで実行）
python .\build_model.py .\build.yaml `
  --model nse --profile cuda_mpi --build

# 外部実行環境ツール
python -m unittest discover ..\RunEnvironment\tests
```

Linux（bash）:

```bash
# CPU/MPI全テスト
python3 ./build_model.py ./build.yaml \
  --model nse --profile cpu_mpi --test

# 単一GPU全テスト
python3 ./build_model.py ./build.yaml \
  --model nse --profile cuda_single --test

# MPI＋CUDAビルド（実機smoke testは複数GPUノードで実行）
python3 ./build_model.py ./build.yaml \
  --model nse --profile cuda_mpi --build

# 外部実行環境ツール
python3 -m unittest discover ../RunEnvironment/tests
```

ハイブリッド流束を変更した場合は、CPUの保存性・切替えテストと
CPU/CUDA一致テストの両方を確認する。

## 11. 衝撃波–乱流干渉ケース

保存乱流SLFの作成後、`case.yaml`で`flow.type: shock_turbulence_interaction`を選ぶ。
Dirichlet駆動面を含む完全な設定例と実行前チェックは
[`NSE_SHOCK_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TURBULENCE_INTERACTION.md)を参照する。
caseを編集した後は必ず`--prepare`を再実行し、Rankine–Hugoniot関係から得た背後状態を
`input.dat`へ反映してからビルド・実行する。

有限高圧室から衝撃波と膨張波を同時に発生させる場合は、
`flow.type: shock_tube_turbulence_interaction`を選ぶ。鏡像閉端、隔膜、高圧・低圧状態、
局所乱流の完全な設定例は
[`NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)を参照する。
