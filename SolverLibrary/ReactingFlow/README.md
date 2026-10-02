# ReactingFlow：独立反応流ライブラリ

2026-10-02追加：[2D非一様初期条件](docs/INITIAL2D.md)。
衝撃波管用split、局所領域box、任意分布profileと初期場書出しを実装。実行検証は保留中。

2026-10-02追加：[2D反応Euler実行入口](docs/FLOW2D.md)を追加。
時間積分・反応分割・境界・CSV出力を接続。ビルドのみ確認、実行テストは保留。
以下の格子追加時点の記録と異なり時間発展は実装されたが、分子輸送・高次化・物理検証は未完了。

2026-10-02追加：[2D構造格子・ノズル生成](docs/GRID2D.md)をFortranに実装。
`rf_mesh2d` で壁面座標から可視化用格子を出力する。流体時間発展は未接続、実行検証は保留中。

2026-10-02：[多次元有限体積基盤](docs/FINITE_VOLUME.md)の先行実装を追加。
任意方向面流束・すべり壁・静止メッシュ残差を提供する。実行検証は保留中で、
多次元反応流・ノズルの実行ソルバーはまだ未完成。既存1Dは変更しない。

2026-10-01：CFDの[試行段の範囲違反回復](docs/OPERATOR_RECOVERY.md)と
[旧NSE実コードとの非反応／一段反応比較](docs/LEGACY_COMPARISON.md)を追加。
[6 µs反応波検証](docs/ACOUSTIC_RETURN_VALIDATION.md)は定常波保持の基準に未達。
短時間合格とは区別し、R4全体は未完了としている。

定常分布の移流・反応釣合いを調べる[Fortran離散残差診断](docs/FLOW_BALANCE.md)を追加。
化学種MUSCLを修正し、1.0 µsの過駆動波試験は512セルで合格。

追加検証：[1.0 µsの過駆動反応波・格子比較と成果物保存](docs/EXTENDED_WAVE_VALIDATION.md)。
同一格子幅での[下流境界距離比較](docs/BOUNDARY_DISTANCE_VALIDATION.md)も明示的に実行できる。

長時間反応波検証用の任意出力：[波面位置・速度履歴](docs/WAVE_HISTORY.md)（Fortran、1D）。

反応流入境界・任意初期分布・CFD反応波比較：[入力と検証範囲](docs/REACTIVE_BOUNDARIES.md)。
CFD反応波の[格子・時間・化学誤差の検証結果](docs/CFD_WAVE_VALIDATION.md)：一次精度・MUSCLとも短時間試験合格。MUSCLの物性下限付近の停止は修正済み。

実在機構の輸送データ生成・Fortran評価：[テーブル輸送の手順](docs/TABULATED_TRANSPORT.md)。

化学平衡・CJ速度探索とCJ連動ZND：[入力・検証範囲](docs/CJ_REFERENCE.md)。
指定速度ZNDの従来入力も維持：[基本仕様](docs/ZND_REFERENCE.md)。R4全体は未完了。

凍結組成の衝撃波基準計算を追加：[入力・実行・制限](docs/SHOCK_REFERENCE.md)。CJ/ZNDではない。

状態：**R1〜R3、R4aの1次元反応Euler、R4bの定係数輸送、R4cのMUSCL、R4dの再試行、R4eの種別定係数拡散をFortranで実装済み。R4全体は進行中。**
従来NSEのStage 0〜10とは別の移行段階で管理する。

単成分NSE（熱揺らぎ拡張を含む）に依存しないライブラリとして再構築する。
既存の`NSE/src/extensions/multicomponent`は検証用・移行元として保持する。
現時点では既存ソルバーの物理的な移動、モデル名変更、実行環境の自動切替は行っていない。
NSE配下に同じ実装を複製して同期する方式にはしない。

## 構成

```text
ReactingFlow/
  CMakeLists.txt      Fortranビルド
  src/fortran/        計算本体（NSE/Python依存なし）
  reference/python/   旧Python版：入力変換・比較検証用として保全
  tools/              機構入力変換・オフライン検証のみ
  examples/           構造検証例とFortran反応器入力例
  tests/              回帰テスト
  docs/               移行計画・各段階の完了条件
```

R0では化学種・元素組成・反応物／生成物の係数を不変データに変換し、
元素保存を有理数で検査する。分子量は元素質量からkg/molで導出する。
未知のキー、誤った単位、重複、NaN、無効な係数は拒否する。
中性気体のみ。電離・表面反応は対象外。原子量は入力責任とし、自動補正しない。

R0の内部YAMLは**Cantera YAMLではない**。構造だけを検証する形式で、速度評価には使わない。
反応速度・NASA物性・第三体等の未対応項目を黙って読み捨てることはしない。
R1の外部入力は別の`import_cantera`を使用する。対応範囲は下記のR1仕様を参照。

## 実行・検証

通常の計算は[Fortran版のビルド・実行手順](docs/REACTORS.md)を使用する。
以下は入力構造の検査のみ。以前のPython版は`reference/python/`へ移し、
Fortran計算時には呼び出さない。旧`run_reactor.py`は比較専用の
`run_reference_reactor.py`へ名称変更した。

FrameWorkルートから、Windows（PowerShell）：

```powershell
python -m pip install PyYAML
python SolverLibrary/ReactingFlow/tools/validate_mechanism.py SolverLibrary/ReactingFlow/examples/topology_demo.yaml
python -m unittest discover -s SolverLibrary/ReactingFlow/tests
```

Linux（bash、利用中の仮想環境内）：

```bash
python3 -m pip install PyYAML
python3 SolverLibrary/ReactingFlow/tools/validate_mechanism.py SolverLibrary/ReactingFlow/examples/topology_demo.yaml
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests
```

`[OK] topology`は元素保存等の検査が通った意味であり、燃焼計算の検証済みを意味しない。
既存のケースは従来の`nse_multicomponent`で実行する。
新ライブラリの実行用manifestは、CFD入口を移植してから登録する。

詳細は[移行順序](docs/MIGRATION.md)を参照。

## R1：外部入力・熱力学

[R1仕様と実行手順](docs/THERMODYNAMICS.md)を追加した。
NASA-7/9物性、混合気体EOS、温度復元、基準量変換はFortran実装。
以前のPython版はオフライン照合用に残す。
入力アダプターと照合テストにのみCantera 3.2.0を使用する。
R2の反応速度評価は[速度仕様](docs/KINETICS.md)を参照。
R3の定容・定圧断熱反応器、保存検査、着火比較は[反応器仕様と実行手順](docs/REACTORS.md)を参照。
R4a〜eの[1次元流体・詳細反応・輸送・MUSCL・再試行の仕様と手順](docs/FLOW1D.md)を追加した。
温度・組成依存のテーブル輸送は実装済み（上記手順参照）。多次元・ノズル、
長時間デトネーション検証、MPI・GPUは独立ライブラリでは未完了。
> R4f update: fixed-state (Dirichlet) boundaries are implemented in the Fortran 1D solver,
> including conservative convective and transport boundary fluxes. See [FLOW1D](docs/FLOW1D.md).
> R4 remains in progress; nonreflecting boundaries, molecular transport and detonation validation remain.

> R4g update: optional common-exponent power-law temperature scaling for transport is implemented.
> Constant models remain available. This is a simplified model, not molecular mixture transport.
> See [FLOW1D](docs/FLOW1D.md) and `examples/flow1d_h2_power_law.in`.

> R4h update: species Sutherland viscosity with Wilke mixture viscosity is available.
> This adds composition-dependent shear viscosity, not a full molecular transport package.
> Parameters are supplied explicitly; see [FLOW1D](docs/FLOW1D.md) and `examples/flow1d_h2_wilke.in`.

> R4i update: optional original-Eucken species conductivity with WMS mixing, NASA7/9 heat capacity,
> and conservative interval-based timestep bounds. See [FLOW1D](docs/FLOW1D.md).
> This remains approximate transport, not validated detailed molecular heat/mass transport.

> R4j: named-species transport files, SI mass checks and input-relative paths are supported.
> See [remaining work](docs/R4_REMAINING.md) for R4 and subsequent R5-R7 scope.

> R4k: composition-dependent mixture-averaged diffusion from prescribed constant binary coefficients.
> See [input, equations and limitations](docs/MIXTURE_DIFFUSION.md). Detailed temperature/pressure-dependent
> binary transport properties and full multicomponent diffusion remain future work.

> R4l: named-species binary diffusion files are supported, with molar-mass checks and complete pair validation.
> See [format and example](docs/MIXTURE_DIFFUSION.md) and `examples/flow1d_h2_binary_file.in`.

> R4m: optional binary diffusion T-power/p scaling, independent of viscosity and conductivity.
> See [parameters and limitations](docs/MIXTURE_DIFFUSION.md). This is not a collision-integral model.

> R4n: named binary diffusion temperature tables, log interpolation and inverse-pressure correction.
> Out-of-range temperatures are rejected. See the same diffusion guide for the file format.

> R4o: recoverable cell-chemistry failures now reject the complete split step and halve dt.
> See [retry scope and limits](docs/CHEMISTRY_RETRY.md). RHS/domain errors remain fatal.

> R4 consolidation 1: local acoustic characteristic outlet and acoustic/conservation tests.
> See [scope, limitations and example](docs/CHARACTERISTIC_OUTLET.md). R4 is not complete.
