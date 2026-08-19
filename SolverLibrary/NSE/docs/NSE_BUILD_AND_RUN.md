# NSE 生成・ビルド・実行手順

更新日: 2026-08-18

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

Windows単一GPU版では、上記にCUDA Toolkit、NVIDIAドライバ、
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

## 3. 新規clone

`2decomp-fft`はGitサブモジュールである。新規取得時は次を使う。

```powershell
git clone --recurse-submodules <SolverLibrary URL>
```

既存cloneへ取得する場合:

```powershell
git submodule update --init --recursive
```

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

CPU MPI + 2DECOMP&FFTのHIT初期化・Forcing:

```yaml
select:
  model: nse
  execution: release

solver:
  profile: cpu_mpi_2decomp_fftw

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: false
```

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

生成前に選択肢と展開内容を確認する。

```powershell
python "$tool\prepare_environment.py" $design --list-options model
python "$tool\prepare_environment.py" $design --dry-run
python "$tool\prepare_environment.py" $design
```

生成ログに表示された`nse_caseNNNN`へ移動する。

```powershell
Set-Location "$env:USERPROFILE\ResearchRuns\nse_caseNNNN"
```

## 5. case設計書

計算条件は`cases\caseNNNN\case.yaml`で管理する。`input.dat`を正本として
直接編集しない。主な項目は次のとおり。

- `grid`: 格子数、領域、ghostセル数
- `time`: CFL、時間刻み、最大ステップ、出力間隔
- `physics.nse`: 比熱比、Mach数、Reynolds数、Prandtl数
- `flow`: Taylor-GreenまたはHIT初期条件
- `forcing`: Forcingの有無と方式
- `numerics`: 対流流束、粘性項、境界条件、時間積分
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
  boundary_condition: periodic
  time_integrator: ssprk3
```

`hybrid`以外の場合、`hybrid`以下の項目は使用されない。設定と混合則は
[`NSE_HYBRID_FLUX.md`](NSE_HYBRID_FLUX.md)を参照する。

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

一時的にMPIランク数とOpenMPスレッド数を上書きする場合:

```powershell
python .\tools\run_case.py --run --processes 8 --omp-threads 2
```

テストを含む実行環境では次を使う。

```powershell
python .\tools\run_case.py --test
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
```

主なプロファイル:

| profile | 用途 |
|---|---|
| `cpu_mpi` | MPI + OpenMP、Taylor-Greenなど |
| `cpu_mpi_2decomp_fftw` | 2DECOMP&FFTを使うHIT初期化・Forcing |
| `cuda_single` | 単一GPU CUDA |

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

# 外部実行環境ツール
python -m unittest discover ..\RunEnvironment\tests
```

ハイブリッド流束を変更した場合は、CPUの保存性・切替えテストと
CPU/CUDA一致テストの両方を確認する。
