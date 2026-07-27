# NSE CMakeビルド手順

> 通常のビルド条件変更にはYAML方式を使用してください。`config/build.yaml`から
> CMake設定を生成する手順は[`YAML_BUILD.md`](YAML_BUILD.md)に記載しています。
> この文書はCMakeを直接操作する場合の手順です。

## 1. 構成

CMakeでは、整理済みの`src`配下だけをビルドします。

- `nse_core`: common、grid、field、io、mpiモジュールをまとめた静的ライブラリ
- `nse_solver`: `src/main/main_nse.f90`から作る実行ターゲット
- Windows出力: `%LOCALAPPDATA%/SolverLibraryBuild/NSE/<preset>/bin/solver.exe`
- Linux出力: `build/<preset>/bin/solver`

従来の`Makefile`は互換確認用として残しています。ルート直下の旧Fortranファイルと
`src/sample`はCMakeターゲットへ含めません。

## 2. Windowsで必要なもの

- CMake 3.24以上
- Ninja
- gfortran
- Microsoft MPI Runtime
- Microsoft MPI SDK

PowerShellで次を確認します。

```powershell
cmake --version
ninja --version
gfortran --version
mpiexec -help
```

`gfortran`がPATHにない場合は、configure時にコンパイラを指定します。

```powershell
cmake --preset windows-msmpi-release `
  -DCMAKE_Fortran_COMPILER=C:/msys64/ucrt64/bin/gfortran.exe
```

## 3. Windows Releaseビルド

NSEディレクトリで実行します。

```powershell
cmake --preset windows-msmpi-release
cmake --build --preset windows-msmpi-release
```

実行例:

```powershell
$env:OMP_NUM_THREADS = "8"
$solver = "$env:LOCALAPPDATA\SolverLibraryBuild\NSE\windows-msmpi-release\bin\solver.exe"
mpiexec -n 4 $solver
```

入力ファイルを相対パスで読む場合は、従来と同じ実行ディレクトリから`solver.exe`を
起動してください。

## 4. Windows Debugビルド

```powershell
cmake --preset windows-msmpi-debug
cmake --build --preset windows-msmpi-debug
```

Debugでは`-O0 -g -fcheck=all -fbacktrace`を有効にします。

## 5. Linux・スパコン

MPIラッパーとNinjaをmoduleで読み込んでから実行します。

```bash
module load gcc
module load openmpi
module load cmake
module load ninja

cmake --preset linux-mpi-release \
  -DCMAKE_Fortran_COMPILER=mpifort
cmake --build --preset linux-mpi-release

export OMP_NUM_THREADS=8
mpiexec -n 4 ./build/linux-mpi-release/bin/solver
```

計算機側にNinjaがない場合は、プリセットを使わずUnix Makefilesを指定できます。

```bash
cmake -S . -B build/linux-make \
  -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=Release \
  -DNSE_MPI_PROVIDER=SYSTEM \
  -DCMAKE_Fortran_COMPILER=mpifort
cmake --build build/linux-make -j 8
```

## 6. 主なCMakeオプション

| オプション | 既定値 | 内容 |
|---|---:|---|
| `NSE_MPI_PROVIDER` | `AUTO` | Windowsでは`MSMPI`、それ以外では`SYSTEM` |
| `MSMPI_ROOT` | Microsoft MPI SDK標準パス | `Include`と`Lib/x64`を持つSDKルート |
| `NSE_ENABLE_WARNINGS` | `ON` | GNU Fortranの`-Wall -Wextra`を有効化 |
| `CMAKE_BUILD_TYPE` | プリセットで指定 | `Release`、`Debug`、`RelWithDebInfo` |

Microsoft MPI構成では、従来のMakefileと同じ`mpif.h`、`msmpifec.lib`、`msmpi.lib`を
使用します。GNU Fortranでは`-fallow-invalid-boz`と`-fallow-argument-mismatch`も
自動的に追加します。

## 7. クリーンビルド

CMakeはソース外のビルドディレクトリへ生成します。Windowsでは、NASや深い日本語パスに
よる260文字制限を避けるため、`%LOCALAPPDATA%/SolverLibraryBuild/NSE`を使用します。
設定を完全に作り直す場合は、対象プリセットのビルドディレクトリを削除してから再度
configureしてください。ソースファイルやNAS上の計算結果は削除しません。
