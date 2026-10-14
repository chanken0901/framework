# NSE：過去の冒頭説明（保管）

2026-10-14にREADMEから分離。以下は追加当時の記録であり、現状は[資料一覧](../README.md)を優先してください。
旧多成分系列の機能は削除していません。新規の独立反応流開発とは区別してください。

2026-10-06：[不等間隔格子の共通基盤](../numerics/NSE_NONUNIFORM_GRID.md)を追加。
現在は座標・微分係数生成のみで、NSE本計算への接続は未完了です。
既存の等間隔計算とCPU/MPI/CUDA/MPI+CUDAの選択機能は維持します。

> Windows（PowerShell）とLinux（bash）の共通コマンド対応は[`../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

> 多成分・反応流の追加開発は、現行単成分NSEと分離した
> [`NSE_MULTICOMPONENT_ROADMAP.md`](../legacy-multicomponent/NSE_MULTICOMPONENT_ROADMAP.md)
> および`solver_manifest_multicomponent.yaml`で管理します。
> Stage 1ではCPU逐次の周期パッシブスカラー移流まで実装済みです。
> Stage 4ではNASA-7物性に加え、混合平均species拡散、Newton粘性、
> Fourier熱伝導まで実装済みです。Stage 5では断熱・定容0次元反応器、
> 一段不可逆Arrhenius反応、SSPRK3時間積分まで実装済みです。
> Stage 6ではStage 4の流体輸送とStage 5の化学反応をStrang分割で結合し、
> 化学反応の自動subcyclingを実装済みです。
> Stage 7では反応流に面別の周期・鏡像・Dirichlet・特性緩和無反射境界、
> x方向反応衝撃波管、途中スナップショットと積分履歴を追加しています。
> Stage 2の共通γモデルも独立profileとして維持しています。
> Stage 8はMPIペンシル分割とOpenMPに対応します。
> Stage 9は平面ノズルの静止一般座標に対応します。設定と制約は
> [一般座標手順書](../legacy-multicomponent/NSE_MULTICOMPONENT_GEOMETRY.md)を参照してください。
> 設定・Windows/Linux手順は[多成分並列計算](../legacy-multicomponent/NSE_MULTICOMPONENT_PARALLEL.md)を参照してください。
