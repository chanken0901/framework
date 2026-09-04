# NSE・GPE共通ビルドランナー

> Windows（PowerShell）とLinux（bash）のコマンド対応は[`../../docs/WINDOWS_LINUX_COMMANDS.md`](../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。本書のGit/CMakeオプションは、OS固有のパスを除いて両環境で共通です。

`build.yaml`の`selected_model`でNSEまたはGPEを選び、同じコマンドで検証、CMake設定生成、
ビルド、テスト、実行を行います。SolverLibrary内の各`solver_manifest.yaml`が、使用可能な
ソルバープロファイルとCMake設定を公開します。

## 1. モデルの選択

[`build.yaml`](build.yaml)を編集します。

```yaml
selected_model: nse
```

GPEへ切り替える場合は次のようにします。

```yaml
selected_model: gpe
```

モデルごとのプロファイルは同じ設計書の`models`で指定します。

```yaml
models:
  nse:
    profile: cpu_mpi

  gpe:
    profile: cpu_mpi_dft
```

利用可能なモデルとプロファイルはコマンドでも確認できます。

```powershell
python .\build_model.py .\build.yaml --list-models
python .\build_model.py .\build.yaml --model gpe --list-profiles
```

Linux（bash）:

```bash
python3 ./build_model.py ./build.yaml --list-models
python3 ./build_model.py ./build.yaml --model gpe --list-profiles
```

## 2. 対応プロファイル

### NSE

| profile | 内容 |
|---|---|
| `cpu_mpi` | MPI領域分割とOpenMP局所ループ |
| `cpu_mpi_2decomp_fftw` | MPI/OpenMPと2DECOMP&FFTによるHIT初期化・Forcing |
| `cuda_single` | 単一GPUのKEEP/WENO/ハイブリッド、`central6`粘性項、SSPRK3 |

`RunEnvironment`から生成したCPU/MPI環境では、通常版と2DECOMP&FFT版が
同梱される。`case.yaml`がHIT初期条件またはPetersen-Livescu forcingを使う場合、
`run_case.py`が`cpu_mpi_2decomp_fftw`を自動選択するため、利用者がprofileを
重複指定する必要はない。

CUDA版のビルドとテスト:

```powershell
python .\build_model.py .\build.yaml `
  --model nse --profile cuda_single --test
```

Linux（bash）:

```bash
python3 ./build_model.py ./build.yaml \
  --model nse --profile cuda_single --test
```

### NSE多成分・反応流拡張

`nse_multicomponent`は現行`nse`と別manifest、別実行ファイルで管理します。
Stage 0は状態レイアウトとprovider契約、Stage 1は一定速度場による保存形式の
周期パッシブスカラー移流、Stage 2は共通γの非反応・非粘性多成分Euler方程式、
Stage 3はNASA-7物性による温度・組成依存熱力学、Stage 4は混合平均拡散、
Newton粘性、Fourier熱伝導、Stage 5は断熱・定容0次元反応器の一段不可逆
Arrhenius反応、Stage 6はStrang分割による流体・反応結合を検証します。

| プロファイル | 内容 |
|---|---|
| `cpu_serial_foundation` | 1成分極限、熱力学・輸送・反応provider契約 |
| `cpu_serial_passive_scalar` | 一次風上法、SSPRK3、周期パッシブスカラー移流 |
| `cpu_serial_inviscid` | 共通γ理想気体、Rusanov流束、非粘性多成分Sod問題 |
| `cpu_serial_thermally_perfect` | NASA-7熱力学、温度反転、温度・組成依存の混合比熱比 |
| `cpu_serial_viscous` | 混合平均species拡散、Newton粘性、Fourier熱伝導 |
| `cpu_serial_reactor` | NASA-7生成エネルギー、一段Arrhenius反応、SSPRK3均質反応器 |
| `cpu_serial_reactive` | Stage 4流体輸送とStage 5化学反応のStrang分割結合 |

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_foundation `
  --validate-only
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_foundation \
  --validate-only
```

Stage 1のビルドと回帰試験:

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_passive_scalar `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_passive_scalar \
  --test
```

Stage 2のビルドと回帰試験:

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_inviscid `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_inviscid \
  --test
```

Stage 3のビルドと回帰試験:

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_thermally_perfect `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_thermally_perfect \
  --test
```

Stage 4のビルドと回帰試験:

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_viscous `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_viscous \
  --test
```

Stage 5のビルドと回帰試験:

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_reactor `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_reactor \
  --test
```

Stage 6のビルドと回帰試験:

```powershell
python .\build_model.py .\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_reactive `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_reactive \
  --test
```

### GPE

| プロファイル | 内容 |
|---|---|
| `cpu_serial_dft` | 検証用DFTによるCPU逐次実行 |
| `cpu_serial_fftw` | FFTWによるCPU逐次実行 |
| `cpu_mpi_dft` | MPI分散FFTと局所DFT |
| `cpu_mpi_fftw` | MPI分散FFTと局所FFTW |
| `cpu_mpi_pencil_fftw` | 2DECOMP&FFTによるMPIペンシル分割FFT |
| `cuda_single` | 単一GPUのCUDA・cuFFT |
| `cuda_mpi_cufftmp` | 1ランク1GPUのMPI・cuFFTMp |
| `cuda_mpi_cufftmp_pencil` | 1ランク1GPUのMPI・cuFFTMpペンシル分割 |

## 3. ビルド

WindowsではローカルFrameWorkのBuildSolverから実行します。

```powershell
Set-Location "$env:USERPROFILE\Documents\Codex\FrameWork\ScriptLibrary\BuildSolver"
```

Linux（bash）:

```bash
cd "$HOME/Research/FrameWork/ScriptLibrary/BuildSolver"
```

NASの`\\Mozart\share\FrameWork`は同期ミラーであり、BuildSolverの直接実行場所には
しません。通常の計算では、RunEnvironmentで生成した外部実行環境内の
`ScriptLibrary\BuildSolver`を`run_case.py`から呼び出します。

設計書の検証だけを行います。

```powershell
python .\build_model.py .\build.yaml --validate-only
```

```bash
python3 ./build_model.py ./build.yaml --validate-only
```

検証、CMake configure、buildをまとめて実行します。

```powershell
python .\build_model.py .\build.yaml
```

```bash
python3 ./build_model.py ./build.yaml
```

コマンドラインで一時的にモデルとプロファイルを変更することもできます。

```powershell
python .\build_model.py .\build.yaml `
  --model gpe `
  --profile cpu_mpi_dft `
  --build
```

```bash
python3 ./build_model.py ./build.yaml \
  --model gpe \
  --profile cpu_mpi_dft \
  --build
```

生成物は既定で次に置かれます。NASへオブジェクトファイルは作りません。

```text
%LOCALAPPDATA%\SolverLibraryBuild\nse-cpu_mpi-release
%LOCALAPPDATA%\SolverLibraryBuild\gpe-cpu_mpi_dft-release
```

## 4. 実行

既存の実行ファイルを再ビルドせずに起動する場合は、`--build`を付けず
`--run`だけを指定します。`--run`ではCMake configure、build、Visual Studio
ビルド環境の初期化を行いません。

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cuda_single --run `
  --input-file "C:\path\to\case\input.nml"
```

```bash
python3 ./build_model.py ./build.yaml \
  --model gpe --profile cuda_single --run \
  --input-file "$HOME/path/to/case/input.nml"
```

### NSE

```powershell
python .\build_model.py .\build.yaml `
  --model nse --build --run `
  --input-file "C:\path\to\case\input.dat" `
  --processes 4 --omp-threads 8
```

```bash
python3 ./build_model.py ./build.yaml \
  --model nse --build --run \
  --input-file "$HOME/path/to/case/input.dat" \
  --processes 4 --omp-threads 8
```

NSEは入力を実行ディレクトリの`input.dat`として読みます。ランナーがローカルの
`build/run/input.dat`へ配置してから`mpiexec`を起動します。

### GPE逐次

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_serial_dft `
  --build --run `
  --input-file "C:\path\to\case\input.nml"
```

```bash
python3 ./build_model.py ./build.yaml \
  --model gpe --profile cpu_serial_dft \
  --build --run \
  --input-file "$HOME/path/to/case/input.nml"
```

### GPE MPI

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_mpi_dft `
  --build --run `
  --input-file "C:\path\to\case\input.nml" `
  --processes 4
```

```bash
python3 ./build_model.py ./build.yaml \
  --model gpe --profile cpu_mpi_dft \
  --build --run \
  --input-file "$HOME/path/to/case/input.nml" \
  --processes 4
```

入力パスを設計書へ保存する場合は、モデルごとの`input_file`へ記述します。

```yaml
models:
  gpe:
    profile: cpu_mpi_dft
    input_file: C:/Research/cases/gpe001/input.nml
```

## 5. FFTW・CUDA・cuFFTMp

外部ライブラリの場所はケース設計書ではなく、
`machine_profiles/windows_gnu_msmpi.yaml`または`linux_gnu_mpi.yaml`へ記述します。

```yaml
libraries:
  fftw_root: C:/fftw
  decomp2d_root: C:/path/to/2decomp-fft
  cuda_compiler: C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.0/bin/nvcc.exe
  cuda_toolkit_root: C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.0
  cuda_architectures: 86
```

cuFFTMpはLinux用マシンプロファイルへ次を設定します。

```yaml
libraries:
  cufftmp_root: /path/to/nvhpc/math_libs
  nvshmem_root: /path/to/nvhpc/comm_libs/nvshmem
  cufftmp_api: auto
```

`cuda_mpi_cufftmp_pencil`は`cufftMpMakePlanDecomposition`を使用するため、
cuFFTMp 11.4.0（NVIDIA HPC SDK 25.3）以降が必要です。

Windowsの`cuda_single`では、ランナーがVisual Studio x64 C++環境を検出してから
CMakeを起動します。`cuda_mpi_cufftmp`はWindowsでは明示的に拒否されます。

## 6. Debug・テスト

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_serial_dft `
  --configuration Debug --build
```

```bash
python3 ./build_model.py ./build.yaml \
  --model gpe --profile cpu_serial_dft \
  --configuration Debug --build
```

CTestを含めて実行します。

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_serial_dft `
  --test
```

```bash
python3 ./build_model.py ./build.yaml \
  --model gpe --profile cpu_serial_dft \
  --test
```

`--dry-run`ではファイル生成や外部コマンド実行を行わず、予定された処理だけを表示します。

## 7. 設定ファイルの役割

| ファイル | 変更する場面 |
|---|---|
| `build.yaml` | モデル、プロファイル、Debug・Release、実行ランク数 |
| `model_catalog.yaml` | 新しい物理モデルをSolverLibraryへ登録するとき |
| `machine_profiles/*.yaml` | 計算機、コンパイラ、MPI、FFTW、CUDAを変更するとき |
| `SolverLibrary/*/solver_manifest.yaml` | ソルバー構成やバックエンドを追加するとき |
| NSE/GPEの入力ファイル | 格子、時間刻み、初期条件、物理パラメータを変更するとき |

ビルド設計と物理条件を分離するため、流れ場や原子種などの条件を`build.yaml`へは
記述しません。

### 条件付きFortran依存関係

Fortranソースの`#if`、`#ifdef`、`#ifndef`で有効になる`use`文がある場合は、
そのプロファイルでCMakeが定義するマクロを`solver_manifest.yaml`にも登録します。

```yaml
profiles:
  cuda_mpi_cufftmp:
    fortran_preprocessor_defines:
      - NSE_INIT_CUFFTMP
```

依存関係検査はこの一覧に従って有効な分岐だけを解析します。値はCMakeの
`target_compile_definitions`と一致させてください。未指定のプロファイルでは
マクロなしとして検査します。

## 8. 再現性情報

各ビルドディレクトリの`generated`へ次を保存します。

- `InitialCache.cmake`: YAMLから生成したCMake初期キャッシュ
- `resolved_build.json`: 選択モデル、プロファイル、コンポーネント、各設計書のSHA-256

生成ファイルを直接編集せず、YAMLまたはソルバーマニフェストを更新して再生成します。
