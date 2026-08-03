# NSE・GPE共通ビルドランナー

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

## 2. 対応プロファイル

### NSE

| profile | 内容 |
|---|---|
| `cuda_single` | 単一GPUのCUDA KEEP・SSPRK3 |

CUDA版のビルドとテスト:

```powershell
python .\build_model.py .\build.yaml `
  --model nse --profile cuda_single --test
```

| プロファイル | 内容 |
|---|---|
| `cpu_mpi` | MPI領域分割とOpenMP局所ループ |

### GPE

| プロファイル | 内容 |
|---|---|
| `cpu_serial_dft` | 検証用DFTによるCPU逐次実行 |
| `cpu_serial_fftw` | FFTWによるCPU逐次実行 |
| `cpu_mpi_dft` | MPI分散FFTと局所DFT |
| `cpu_mpi_fftw` | MPI分散FFTと局所FFTW |
| `cuda_single` | 単一GPUのCUDA・cuFFT |
| `cuda_mpi_cufftmp` | 1ランク1GPUのMPI・cuFFTMp |

## 3. ビルド

WindowsではローカルFrameWorkのBuildSolverから実行します。

```powershell
Set-Location "$env:USERPROFILE\Documents\Codex\FrameWork\ScriptLibrary\BuildSolver"
```

NASの`\\Mozart\share\FrameWork`は同期ミラーであり、BuildSolverの直接実行場所には
しません。通常の計算では、RunEnvironmentで生成した外部実行環境内の
`ScriptLibrary\BuildSolver`を`run_case.py`から呼び出します。

設計書の検証だけを行います。

```powershell
python .\build_model.py .\build.yaml --validate-only
```

検証、CMake configure、buildをまとめて実行します。

```powershell
python .\build_model.py .\build.yaml
```

コマンドラインで一時的にモデルとプロファイルを変更することもできます。

```powershell
python .\build_model.py .\build.yaml `
  --model gpe `
  --profile cpu_mpi_dft `
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

### NSE

```powershell
python .\build_model.py .\build.yaml `
  --model nse --build --run `
  --input-file "C:\path\to\case\input.dat" `
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

### GPE MPI

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_mpi_dft `
  --build --run `
  --input-file "C:\path\to\case\input.nml" `
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

Windowsの`cuda_single`では、ランナーがVisual Studio x64 C++環境を検出してから
CMakeを起動します。`cuda_mpi_cufftmp`はWindowsでは明示的に拒否されます。

## 6. Debug・テスト

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_serial_dft `
  --configuration Debug --build
```

CTestを含めて実行します。

```powershell
python .\build_model.py .\build.yaml `
  --model gpe --profile cpu_serial_dft `
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

## 8. 再現性情報

各ビルドディレクトリの`generated`へ次を保存します。

- `InitialCache.cmake`: YAMLから生成したCMake初期キャッシュ
- `resolved_build.json`: 選択モデル、プロファイル、コンポーネント、各設計書のSHA-256

生成ファイルを直接編集せず、YAMLまたはソルバーマニフェストを更新して再生成します。
