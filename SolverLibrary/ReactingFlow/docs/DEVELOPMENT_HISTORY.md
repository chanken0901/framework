# ReactingFlow：開発履歴（保管）

2026-10-14にREADME冒頭から分離。日付ごとの記述は当時の状態です。
現在の制約は[ライブラリ入口](../README.md)、[R4残作業](R4_REMAINING.md)と各機能文書を参照してください。

2026-10-02追加：[2D非一様初期条件](INITIAL2D.md)。
衝撃波管用split、局所領域box、任意分布profileと初期場書出しを実装。実行検証は保留中。

2026-10-02追加：[2D反応Euler実行入口](FLOW2D.md)を追加。
時間積分・反応分割・境界・CSV出力を接続。ビルドのみ確認、実行テストは保留。
以下の格子追加時点の記録と異なり時間発展は実装されたが、分子輸送・高次化・物理検証は未完了。

2026-10-02追加：[2D構造格子・ノズル生成](GRID2D.md)をFortranに実装。
`rf_mesh2d` で壁面座標から可視化用格子を出力する。流体時間発展は未接続、実行検証は保留中。

2026-10-02：[多次元有限体積基盤](FINITE_VOLUME.md)の先行実装を追加。
任意方向面流束・すべり壁・静止メッシュ残差を提供する。実行検証は保留中で、
多次元反応流・ノズルの実行ソルバーはまだ未完成。既存1Dは変更しない。

2026-10-01：CFDの[試行段の範囲違反回復](OPERATOR_RECOVERY.md)と
[旧NSE実コードとの非反応／一段反応比較](LEGACY_COMPARISON.md)を追加。
[6 µs反応波検証](ACOUSTIC_RETURN_VALIDATION.md)は定常波保持の基準に未達。
短時間合格とは区別し、R4全体は未完了としている。

定常分布の移流・反応釣合いを調べる[Fortran離散残差診断](FLOW_BALANCE.md)を追加。
化学種MUSCLを修正し、1.0 µsの過駆動波試験は512セルで合格。

追加検証：[1.0 µsの過駆動反応波・格子比較と成果物保存](EXTENDED_WAVE_VALIDATION.md)。
同一格子幅での[下流境界距離比較](BOUNDARY_DISTANCE_VALIDATION.md)も明示的に実行できる。

長時間反応波検証用の任意出力：[波面位置・速度履歴](WAVE_HISTORY.md)（Fortran、1D）。

反応流入境界・任意初期分布・CFD反応波比較：[入力と検証範囲](REACTIVE_BOUNDARIES.md)。
CFD反応波の[格子・時間・化学誤差の検証結果](CFD_WAVE_VALIDATION.md)：一次精度・MUSCLとも短時間試験合格。MUSCLの物性下限付近の停止は修正済み。

実在機構の輸送データ生成・Fortran評価：[テーブル輸送の手順](TABULATED_TRANSPORT.md)。

化学平衡・CJ速度探索とCJ連動ZND：[入力・検証範囲](CJ_REFERENCE.md)。
指定速度ZNDの従来入力も維持：[基本仕様](ZND_REFERENCE.md)。R4全体は未完了。

凍結組成の衝撃波基準計算を追加：[入力・実行・制限](SHOCK_REFERENCE.md)。CJ/ZNDではない。
