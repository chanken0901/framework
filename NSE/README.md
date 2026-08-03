# NSE SolverLibrary

統計的定常乱流のPetersen-Livescu Forcing、CPU分散FFT、単一GPU cuFFT、
および将来のcuFFTMp拡張境界については
[`docs/NSE_FORCING.md`](docs/NSE_FORCING.md)を参照してください。

## 単一GPU CUDA版

`cuda_single` profileでは、周期境界、KEEP対流流束、CFL評価、SSPRK3を
単一GPUで実行します。保存変数は時間発展中GPUへ常駐し、SLF出力時だけ
CPUへ転送します。

```powershell
python ..\..\ScriptLibrary\BuildSolver\build_model.py `
  ..\..\ScriptLibrary\BuildSolver\build.yaml `
  --model nse --profile cuda_single --test
```

対応範囲、実行環境の生成、GPUメモリ量は
[`docs/NSE_CUDA.md`](docs/NSE_CUDA.md)を参照してください。

初期条件、移流スキーム、粘性項、境界条件の追加方法は
[`docs/NSE_MODULE_DESIGN.md`](docs/NSE_MODULE_DESIGN.md)を参照してください。

## 選択可能なKEEP対流項

CPU/MPI/OpenMP版と単一GPU CUDA版のKEEP対流項は、2次精度と6次精度を
実行時に選択できます。既定値は6次精度です。どちらも対称な二点KEEP流束を
保存形に合成し、KEEPの運動エネルギー・内部エネルギー保存構造を維持します。

```yaml
numerics:
  convective_scheme: keep6  # keep2 または keep6
```

次数はスキーム名に含めて指定します。`keep` 単独の指定は使用できません。

現在の周期境界実装との共通化のため、どちらの精度でもゴーストセル数は3です。

## 6次精度粘性項

CPU/MPI/OpenMP版と単一GPU CUDA版の両方で、`central6`粘性項を利用できます。
一定粘性係数のNewton流体、Stokesの仮定、Fourier熱伝導を用い、一階・
二階・混合微分を陽的6次精度中心差分で評価します。必要なゴーストセル数は
3です。

`case.yaml`では次のように有効化します。

```yaml
physics:
  nse:
    reynolds_number: 100.0
    prandtl_number: 0.72

numerics:
  viscous_scheme: central6
```

`viscous_scheme: none`を指定すると、同じCPU/GPUプロファイルのまま粘性項を
無効化できます。自動時間刻みでは、対流CFL条件に加えて6次差分の拡散安定
条件も適用されます。

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

## ParaView可視化

外部実行環境では、caseのSLF出力を次の共通コマンドで変換します。

```powershell
python .\tools\postprocess_case.py
```

既定値は最新の完全なステップ、`rho,u,v,w,p`、空間間引き2です。全解像度で
指定ステップを変換する場合は、次のように指定します。

```powershell
python .\tools\postprocess_case.py --steps 0,500,1000 --stride 1
```

生成された`cases\<case_id>\paraview\collection.pvd`をParaViewで開きます。

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
| `tools` | YAML検証、CMake設定生成、ビルド、実行、SLFからVTIへの変換 |

CMakeは上記の整理済みソースだけを使用します。ルート直下の旧ソースと`src/sample`は
ビルド対象外です。
