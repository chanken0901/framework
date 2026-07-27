# NSE SolverLibrary

MPI・OpenMP対応の三次元Navier-Stokes方程式ソルバーです。

## 推奨ビルド

通常はCMakeファイルを直接編集せず、`config/build.yaml`を設計書としてビルドします。
スクリプトがモジュール依存関係を解決し、CMake設定をローカルビルドディレクトリへ
生成します。

```powershell
python .\tools\build_from_yaml.py .\config\build.yaml
```

NSEとGPEを同じ設計書から選択するフレームワークでは、
`ScriptLibrary/BuildSolver/build_model.py`を使用します。NSE単独ランナーは既存ケースとの
互換性のため残しています。

検証、Debugビルド、MPI実行、NAS上での利用方法は
[`YAML_BUILD.md`](YAML_BUILD.md)を参照してください。

主な設定場所は次の通りです。

| 変更内容 | 編集ファイル |
|---|---|
| Debug・Release、警告、出力先、並列ビルド数 | `config/build.yaml` |
| コンパイラ、MPI、OpenMP、最適化フラグ | `config/build_profiles/*.yaml` |
| ソースモジュールと依存関係 | `config/module_catalog.yaml` |
| CMakeターゲットそのもの | `CMakeLists.txt` |

## CMake直接実行

従来のCMakeプリセットも引き続き利用できます。

```powershell
cmake --preset windows-msmpi-release
cmake --build --preset windows-msmpi-release
```

必要なソフトウェアとCMake直接実行の詳細は[`CMAKE_BUILD.md`](CMAKE_BUILD.md)を
参照してください。

従来の`Makefile`とMS-MPI用バッチファイルは、既存環境の比較・移行確認用として
残しています。

## ソース構成

| ディレクトリ | 内容 |
|---|---|
| `src/common` | 精度、定数、共通設定、モデル設定 |
| `src/grid` | FVM格子 |
| `src/field` | NSE保存変数と作業配列 |
| `src/io` | 入力条件とSLF出力 |
| `src/mpi` | MPI領域分割と通信 |
| `src/main` | NSEメインプログラム |
| `config` | YAMLビルド設計書、モジュールカタログ、マシンプロファイル |
| `tools` | YAML検証、CMake設定生成、ビルド、実行スクリプト |

CMakeは上記の整理済みソースだけを使用します。ルート直下の旧ソースと`src/sample`は
ビルド対象外です。
