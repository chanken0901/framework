# NSE YAMLビルド手順

## 1. 目的

`config/build.yaml`をビルド設計書として読み込み、次の処理を自動化します。

1. モジュールと依存関係の検証
2. マシンプロファイルの読み込み
3. CMake初期キャッシュの生成
4. CMake configureとbuild
5. 必要に応じたMPI・OpenMP実行

生成ファイルとオブジェクトはSolverLibraryには書き込まず、既定では
`%LOCALAPPDATA%/SolverLibraryBuild/NSE`へ置きます。NAS上のNSEを参照する場合も、
`cmd.exe`の作業ディレクトリがUNCパスになりません。

## 2. 設定ファイル

### `config/build.yaml`

ケースまたはビルド単位で変更する設計書です。

| 項目 | 内容 |
|---|---|
| `machine_profile` | 使用するマシンプロファイル |
| `solver.module_set` | 使用するモジュール集合 |
| `solver.executable_name` | 実行ファイル名 |
| `build.configuration` | `Release`、`Debug`、`RelWithDebInfo` |
| `build.output_root` | ローカルビルド出力先 |
| `build.parallel_jobs` | コンパイルの並列数 |
| `run.mpi_processes` | MPIランク数 |
| `run.omp_threads` | 1ランク当たりのOpenMPスレッド数 |

### `config/build_profiles/*.yaml`

計算機やコンパイラに依存する設定です。

- `windows_msmpi.yaml`: Windows、gfortran、Microsoft MPI
- `linux_system_mpi.yaml`: Linux・スパコン、`mpifort`、システムMPI

コンパイラの場所、最適化フラグ、MPI実装を変更する場合はこちらを編集します。

### `config/module_catalog.yaml`

NSEモジュールのソースパス、役割、依存関係を定義します。新しいFortranモジュールを
追加したときに更新します。通常のDebug・Release切り替えでは編集しません。

## 3. Windowsでのビルド

NSEディレクトリで実行します。

```powershell
python .\tools\build_from_yaml.py .\config\build.yaml
```

オプションを指定しない場合は、検証、生成、configure、buildまで行います。
実行ファイルは次に生成されます。

```text
%LOCALAPPDATA%\SolverLibraryBuild\NSE\windows-msmpi-release\bin\solver.exe
```

設計書の検証だけを行う場合は次です。

```powershell
python .\tools\build_from_yaml.py .\config\build.yaml --validate-only
```

生成されたCMake設定だけを確認する場合は次です。

```powershell
python .\tools\build_from_yaml.py .\config\build.yaml --generate-only
```

## 4. Debugビルド

設計書を変更せず、一時的にDebugへ切り替えられます。

```powershell
python .\tools\build_from_yaml.py .\config\build.yaml `
  --configuration Debug `
  --build-dir "$env:LOCALAPPDATA\SolverLibraryBuild\NSE\windows-msmpi-debug" `
  --build
```

恒久的に変更する場合は`config/build.yaml`の`build.configuration`と
`build.directory_name`を編集します。

## 5. MPI実行

NSEは実行ディレクトリにある`input.dat`を読みます。スクリプトは指定された入力を
ローカルの`build/run/input.dat`へコピーしてから起動します。

```powershell
python .\tools\build_from_yaml.py .\config\build.yaml `
  --build --run `
  --input-file "C:\path\to\case\input.dat" `
  --processes 4 `
  --omp-threads 8
```

`--build --run`では、YAML生成、configure、build、入力配置、`mpiexec`実行までを
一つのコマンドで行います。

## 6. ローカルFrameWorkで使う場合

PowerShellからローカルのNSEソルバーへ移動します。

```powershell
Set-Location "$env:USERPROFILE\Documents\Codex\FrameWork\SolverLibrary\NSE"
python .\tools\build_from_yaml.py .\config\build.yaml
```

現在の標準運用では、RunEnvironmentで`${USERPROFILE}\ResearchRuns`へ外部実行環境を
生成し、その中でビルド・実行します。NASの`\\Mozart\share\FrameWork`は同期ミラーとし、
UNCパスを作業ディレクトリにしてビルドしません。

## 7. Linux・スパコン

`config/build.yaml`をケース側へ複製し、次のようにプロファイルを変更します。

```yaml
machine_profile: build_profiles/linux_system_mpi.yaml

build:
  configuration: Release
  output_root: ${HOME}/.cache/SolverLibraryBuild/NSE
  directory_name: linux-mpi-release
  warnings: true
  parallel_jobs: 16
  configure_fresh: true
```

実行例です。

```bash
python3 tools/build_from_yaml.py config/build.yaml --build --run \
  --input-file /path/to/case/input.dat \
  --processes 16 --omp-threads 4
```

ジョブスケジューラを使用する場合は、ビルド後の`bin/solver`を`mpiexec`または
`srun`から起動します。`tools.mpi_launcher`はマシンプロファイルで変更できます。

## 8. 生成物

ビルドディレクトリの`generated`には次が保存されます。

- `NSEInitialCache.cmake`: YAMLから生成したCMake設定
- `resolved_build.json`: 設計書、プロファイル、選択モジュール、ハッシュ値

これらは再現性確認用です。直接編集せず、YAMLを変更して再生成してください。
