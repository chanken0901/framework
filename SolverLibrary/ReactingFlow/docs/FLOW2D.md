# 実験的2次元反応Euler実行入口

2026-10-02。CPU逐次Fortranの `rf_flow2d` を追加した。
**ビルド確認のみ。実行テストはユーザー指示により保留中。**
6 µs反応波の未達は未解決であり、実用燃焼・デトネーション計算の検証完了ではない。

## 今回接続したもの

- [ノズル格子](GRID2D.md)と[面流束](FINITE_VOLUME.md)。2次元平面、静止格子。
- 一次精度Rusanov空間離散化、非分割SSPRK3流体時間積分。
- 化学半段→流体全段→化学半段のStrang分割。時間結合は2次精度を意図するが未検証。
- 1Dと共通のDVODE化学積分・元素検査・内部エネルギー整合性検査。
  化学中は両運動量成分と生成エネルギー込み全エネルギーを保持し、反応熱は別加算しない。
- 各段の状態適合性・CFL確認。回復可能な失敗時は全段を棄却し、元の状態からdtを半減して再試行。
  最大30試行。採用段の境界流束のみ収支に加える。クリッピングはしない。
- 境界は `outflow`（ゼロ勾配）、`reflecting`（すべり鏡像）、`dirichlet`（固定ghost）。
  ゼロ勾配は無反射境界を意味しない。

状態は `[rhoY...,rho*u,rho*v,rhoE]`。化学の再利用は運動量の大きさを1D化学APIへ渡す
アダプターで行い、2Dの運動エネルギーを落とさない。
既存1Dの数値処理は変更せず、`chemistry_cells` の公開範囲のみ拡張した。

## 入力

`examples/flow2d_h2air.in` を参照。`&flow2d` の主な項目：

- `mesh_profile`, `ny`: 格子。相対パスは入力ファイルのあるフォルダ基準。
- `temperature`, `pressure`, `velocity(2)`, `mass_fractions(ns)`: 一様初期状態（SI）。
- `boundary_kind(4)`: 順序はi最小、i最大、j最小、j最大。
- 固定面のみ `boundary_temperature(4)`, `boundary_pressure(4)`,
  `boundary_velocity(2,4)`, `boundary_y(ns,4)` を与える。未指定の組成は拒否する。
- `chemistry`, `chemistry_rtol`, `chemistry_atol_species`, `chemistry_atol_temperature`,
  `chemistry_max_steps`: 詳細反応の使用と積分制御。
- `end_time`, `max_dt`, `cfl`, `max_steps`, `write_every`: 時間・出力制御。

組成順序は機構の化学種順序。入力例はCantera h2o2の
`H2,H,O,O2,OH,H2O,HO2,H2O2,AR,N2` 専用。他の機構にそのまま使用しない。
質量分率の自動正規化は行わない。対応していない境界名やnamelist項目は拒否する。

## ビルド・実行

FrameWorkルートから実行。既存の機構変換で作成した `mechanism.rf` を用意する。
入力例は動作検証用の短時間・粗格子設定で、物理検証済み条件ではない。

Windows PowerShell：

```powershell
cmake --build build/reactingflow-fortran-release --target rf_flow2d
.\build\reactingflow-fortran-release\rf_flow2d.exe mechanism.rf SolverLibrary/ReactingFlow/examples/flow2d_h2air.in build/flow2d.csv
```

Linux：

```bash
cmake --build build/reactingflow-fortran-release --target rf_flow2d
./build/reactingflow-fortran-release/rf_flow2d mechanism.rf SolverLibrary/ReactingFlow/examples/flow2d_h2air.in build/flow2d.csv
```

出力先の親フォルダを用意すること。既存CSVは上書きしない。
初期状態・指定間隔・最終時刻をCSVへ出力し、正常終了メッセージは最終書込み後に表示する。
CSVは各セルの中心座標・面積・密度・速度2成分・温度・圧力・質量分率を持つ。
末尾には境界流束を差し引いた質量・運動量・全エネルギーの**絶対**収支残差と棄却回数を出力する。
化学種別質量は反応で変化するので個々の収支をゼロとは判定しない。
`# SUCCESS` は正常終了のみを意味し、物理的妥当性や保存誤差の合格判定ではない。

## 未実装・未検証

多次元MUSCL、分子輸送（粘性・熱伝導・種拡散）、特性境界・周期接続、非一様初期状態、
再始動、MPI/OpenMP/CUDA、3D、軸対称、case.yaml/実行環境の運用統合は未実装。
従ってNavier–Stokesのノズル燃焼実用ソルバーは未完成。
反応速度内部の致命的エラーや浮動小数点オーバーフローの全てを再試行で回復するわけではない。
現時点では各段で幾何検査と一時配列確保を行う基準実装であり、大規模計算用の高速化は未着手。
検証コード `rf_flow2d_unit` は追加・ビルドのみで、合格とは扱わない。
