# ReactingFlow資料：用途別一覧

独立Fortran反応流ライブラリの資料。旧NSE多成分Stage系列とは区別する。
**R4全体は未完了**。実装済み・短時間試験合格・長時間検証済みを同じ意味に扱わない。

2026-10-14に本文を内容別フォルダへ移動した。`physics/`は物性・反応・輸送、
`numerics/`は離散化・積分、`flow/`はCFD・格子・境界、`reference/`は基準計算、
`validation/`は検証と診断、`development/`は計画・進捗・履歴を扱う。

## 入口・進捗

- [ライブラリ概要と実行手順](../README.md)
- [移行方針](development/MIGRATION.md)
- [R4の残作業・完了条件](development/R4_REMAINING.md)
- [過去の追加記録](development/DEVELOPMENT_HISTORY.md)

## 物性・化学反応・輸送

- [熱力学](physics/THERMODYNAMICS.md)／[反応速度](physics/KINETICS.md)
- [反応器](physics/REACTORS.md)／[反応器の参照検証](validation/REFERENCE_REACTORS.md)
- [化学積分の再試行](numerics/CHEMISTRY_RETRY.md)
- [混合気体の拡散](physics/MIXTURE_DIFFUSION.md)
- [テーブル輸送](physics/TABULATED_TRANSPORT.md)

## CFD・格子・境界

- [1D流体計算](flow/FLOW1D.md)
- [2D流体計算](flow/FLOW2D.md)／[2D格子](flow/GRID2D.md)／[2D初期条件](flow/INITIAL2D.md)
- [有限体積基盤](numerics/FINITE_VOLUME.md)
- [反応境界・初期分布](flow/REACTIVE_BOUNDARIES.md)
- [特性流出境界](flow/CHARACTERISTIC_OUTLET.md)
- [試行段の範囲違反回復](numerics/OPERATOR_RECOVERY.md)

## 衝撃波・デトネーションの基準計算

- [凍結組成の衝撃波](reference/SHOCK_REFERENCE.md)
- [CJ速度](reference/CJ_REFERENCE.md)／[ZND](reference/ZND_REFERENCE.md)
- [波面履歴](validation/WAVE_HISTORY.md)
- [流れの離散残差診断](validation/FLOW_BALANCE.md)

## 検証記録

- [CFD反応波検証](validation/CFD_WAVE_VALIDATION.md)
- [1 µs拡張試験](validation/EXTENDED_WAVE_VALIDATION.md)
- [6 µs試験：基準未達の記録](validation/ACOUSTIC_RETURN_VALIDATION.md)
- [境界距離比較](validation/BOUNDARY_DISTANCE_VALIDATION.md)
- [旧NSE系列との比較](validation/LEGACY_COMPARISON.md)

2D文書に記載されたビルド確認だけで、物理検証や本番運用が完了したとは判断しない。
