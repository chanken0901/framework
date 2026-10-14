# NSE資料：用途別一覧

単成分NSEの現行仕様を先に示す。旧多成分系列と開発履歴は末尾に分ける。
記載の機能を自由に組み合わせられるとは限らないため、各文書の対応範囲を確認する。

2026-10-14に本文を内容別フォルダへ移動した。`overview/`は全体仕様、`usage/`は操作、
`numerics/`は数値計算法・格子・境界、`physics/`は物理設定・初期条件、
`parallel/`はGPU・通信、`postprocess/`は解析、`legacy-multicomponent/`は旧系列、
`history/`は過去の記録。以下のリンクから該当文書を開く。

## 計算の開始・継続

- [全体仕様・無次元化](overview/NSE_SOLVER_SPECIFICATION.md)
- [ビルドと実行](usage/NSE_BUILD_AND_RUN.md)
- [SLF再スタート](usage/NSE_RESTART.md)
- [境界条件](numerics/NSE_BOUNDARY_CONDITIONS.md)
- [不等間隔格子](numerics/NSE_NONUNIFORM_GRID.md)

## 物理設定・数値計算

- [揺らぎ拡張：周期境界に限定](physics/NSE_FLUCTUATING_HYDRODYNAMICS.md)
- [HITと分散FFT](physics/NSE_HIT_DISTRIBUTED_FFT.md)
- [Forcing](physics/NSE_FORCING.md)／[ターゲットの指定](physics/NSE_FORCING_TARGET.md)
- [保存乱流の読み込み・配置](physics/NSE_IMPORTED_TURBULENCE.md)
- [衝撃波と乱流の干渉](physics/NSE_SHOCK_TURBULENCE_INTERACTION.md)
- [有限高圧室による衝撃波管・乱流干渉](physics/NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)
- [ハイブリッド流束](numerics/NSE_HYBRID_FLUX.md)／[WENO5Z＋Roe](numerics/NSE_WENO5Z_ROE.md)
- [正値性と停止判定](numerics/NSE_POSITIVITY.md)

## GPU・通信

- [CUDA](parallel/NSE_CUDA.md)
- [CUDA-aware MPI](parallel/NSE_CUDA_AWARE_MPI.md)
- [分散FFT：cuFFTMp](parallel/NSE_CUFFTMP.md)

## 後処理

- [乱流統計](postprocess/NSE_TURBULENCE_STATISTICS.md)
- [SLFのFFT解析](postprocess/NSE_SLF_FFT.md)
- [衝撃波位置・速度・前後統計](postprocess/NSE_SHOCK_ANALYSIS.md)
- [不等間隔SLFの可視化・統計](postprocess/NSE_NONUNIFORM_POSTPROCESS.md)

## 開発・互換系列

- [モジュール設計](overview/NSE_MODULE_DESIGN.md)
- [旧NSE多成分系列のロードマップ](legacy-multicomponent/NSE_MULTICOMPONENT_ROADMAP.md)
- 旧系列：[並列化](legacy-multicomponent/NSE_MULTICOMPONENT_PARALLEL.md)／[CUDA](legacy-multicomponent/NSE_MULTICOMPONENT_CUDA.md)／[一般座標](legacy-multicomponent/NSE_MULTICOMPONENT_GEOMETRY.md)／[反応境界](legacy-multicomponent/NSE_MULTICOMPONENT_REACTIVE_BOUNDARIES.md)
- 詳細反応機構の新規開発先：[独立ReactingFlow](../../ReactingFlow/README.md)
- [過去のREADME冒頭記録](history/DEVELOPMENT_HISTORY.md)

Stage 0〜10と独立ReactingFlowのR0〜R4は別系列であり、対応状況を読み替えない。
