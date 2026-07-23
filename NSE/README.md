# NSE SolverLibrary

MPI・OpenMP対応の三次元Navier-Stokes方程式ソルバーです。

## ビルド

新しい標準ビルドはCMakeです。Windowsのgfortran＋Microsoft MPI、Linux・HPCの
システムMPIに対応しています。

```powershell
cmake --preset windows-msmpi-release
cmake --build --preset windows-msmpi-release
```

必要なソフトウェア、実行方法、Debugビルド、HPC向け設定は
[`CMAKE_BUILD.md`](CMAKE_BUILD.md)を参照してください。

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

CMakeは上記の整理済みソースだけを使用します。ルート直下の旧ソースと`src/sample`は
ビルド対象外です。
