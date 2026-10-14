# NSE資料：用途別一覧

単成分NSEの現行仕様を先に示す。旧多成分系列と開発履歴は末尾に分ける。
記載の機能を自由に組み合わせられるとは限らないため、各文書の対応範囲を確認する。

## 計算の開始・継続

- [全体仕様・無次元化](NSE_SOLVER_SPECIFICATION.md)
- [ビルドと実行](NSE_BUILD_AND_RUN.md)
- [SLF再スタート](NSE_RESTART.md)
- [境界条件](NSE_BOUNDARY_CONDITIONS.md)
- [不等間隔格子](NSE_NONUNIFORM_GRID.md)

## 物理設定・数値計算

- [揺らぎ拡張：周期境界に限定](NSE_FLUCTUATING_HYDRODYNAMICS.md)
- [HITと分散FFT](NSE_HIT_DISTRIBUTED_FFT.md)
- [Forcing](NSE_FORCING.md)／[ターゲットの指定](NSE_FORCING_TARGET.md)
- [保存乱流の読み込み・配置](NSE_IMPORTED_TURBULENCE.md)
- [衝撃波と乱流の干渉](NSE_SHOCK_TURBULENCE_INTERACTION.md)
- [有限高圧室による衝撃波管・乱流干渉](NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)
- [ハイブリッド流束](NSE_HYBRID_FLUX.md)／[WENO5Z＋Roe](NSE_WENO5Z_ROE.md)
- [正値性と停止判定](NSE_POSITIVITY.md)

## GPU・通信

- [CUDA](NSE_CUDA.md)
- [CUDA-aware MPI](NSE_CUDA_AWARE_MPI.md)
- [分散FFT：cuFFTMp](NSE_CUFFTMP.md)

## 後処理

- [乱流統計](NSE_TURBULENCE_STATISTICS.md)
- [SLFのFFT解析](NSE_SLF_FFT.md)
- [衝撃波位置・速度・前後統計](NSE_SHOCK_ANALYSIS.md)
- [不等間隔SLFの可視化・統計](NSE_NONUNIFORM_POSTPROCESS.md)

## 開発・互換系列

- [モジュール設計](NSE_MODULE_DESIGN.md)
- [旧NSE多成分系列のロードマップ](NSE_MULTICOMPONENT_ROADMAP.md)
- 旧系列：[並列化](NSE_MULTICOMPONENT_PARALLEL.md)／[CUDA](NSE_MULTICOMPONENT_CUDA.md)／[一般座標](NSE_MULTICOMPONENT_GEOMETRY.md)／[反応境界](NSE_MULTICOMPONENT_REACTIVE_BOUNDARIES.md)
- 詳細反応機構の新規開発先：[独立ReactingFlow](../../ReactingFlow/README.md)
- [過去のREADME冒頭記録](DEVELOPMENT_HISTORY.md)

Stage 0〜10と独立ReactingFlowのR0〜R4は別系列であり、対応状況を読み替えない。
