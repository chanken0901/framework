# WENO5-Z/Roe対流流束

## 概要

CPU逐次、MPI、OpenMP、およびMPI+OpenMP版では、圧縮性Navier-Stokes方程式の
対流項として`weno5z_roe`を選択できます。セル中心の保存変数から特性空間で
左右状態を5次精度WENO-Z再構築し、Roe近似Riemann solverで面流束を求めます。

```yaml
numerics:
  convective_scheme: weno5z_roe
```

単一GPU CUDA版でも`weno5z_roe`を選択できます。CUDA版は保存変数をGPUに
常駐させ、各セルで必要な左右6面のWENO/Roe流束を直接評価します。

## モジュール構成

| 層 | ファイル | 役割 |
|---|---|---|
| 再構築 | `src/numerics/reconstruction/mod_reconstruction_weno5z.f90` | 左右のWENO5-Z再構築 |
| Riemann solver | `src/numerics/riemann/mod_riemann_roe.f90` | Roe平均、固有分解、entropy fix、数値流束 |
| 複合流束 | `src/numerics/convective/mod_convective_weno5z_roe.f90` | 特性射影、三方向の面流束計算 |
| CUDA流束 | `src/gpu/nse_cuda_weno5z_roe.cuh` | GPU上の特性再構築、Roe流束、三方向RHS |
| KEEP | `src/numerics/convective/mod_convective_keep.f90` | KEEP2/KEEP6流束 |
| 選択 | `src/numerics/convective/mod_convective_dispatch.f90` | 実行時スキーム選択 |

空間演算、SSPRK3、MPI領域分割は`mod_convective_scheme`の共通APIだけを呼びます。
そのため、WENO/Roe固有処理は時間積分やMPIコードへ漏れません。

## 数値手法

- WENO-Zの大域滑らかさ指標は`tau5 = abs(beta0-beta2)`です。
- 非線形重みは指数2、既定のepsilonは`1.0e-20`です。
- 左右再構築には鏡映した最適重み`(0.1,0.6,0.3)`と
  `(0.3,0.6,0.1)`を使用します。
- 各面のRoe平均から左右固有ベクトルを作り、5点ステンシルを特性空間へ
  射影してから再構築します。
- Roe固有値にはHarten-Hyman型entropy fixを適用します。
- x、y、z方向は保存変数と流束を法線座標へ回転して同じ一次元solverを使います。
- 必要なゴーストセル数は3です。

参照コードの右側候補多項式の並びに対して、最適重みを左側と同じ
`(0.1,0.6,0.3)`にすると滑らかな領域でも右再構築が5次精度になりません。
本実装では候補の並びに対応する鏡映重みへ修正しています。

## 検証

`tests/test_weno5z_roe.f90`で次を検査します。

1. 左右WENO-Z再構築の5次収束
2. Roe左右固有ベクトルの逆行列関係
3. 左右同一状態で数値流束と物理流束が一致すること
4. x、y、z方向の流束回転
5. 滑らかなentropy waveの5次収束と保存性
6. Sod不連続面で流束が有限値を保つこと

CUDAビルドでは`nse_cuda_weno5z_roe_compare`により、同じ三次元初期場を
CPU版とCUDA版でSSPRK3の一ステップだけ進め、保存変数が一致することを検査します。

## 実装済みKEEP/Roeハイブリッドとの連携

`convective_scheme: hybrid`では、Ducros-pressureセンサーにより面ごとに
KEEP流束とWENO5-Z/Roe流束を連続混合します。現在の依存方向は次です。

```text
mod_convective_dispatch
  -> mod_convective_hybrid
       -> mod_convective_leaf_registry
            -> mod_convective_keep
            -> mod_convective_weno5z_roe
```

ハイブリッド層はセンサー評価と混合率の計算を担当し、構成流束の実計算は
leaf registryを介して既存のKEEPおよびWENO5-Z/Roe実装を再利用します。
混合率が0と1の間では両方の面流束を計算して線形混合し、0または1では必要な
構成流束だけを計算します。

## 制約

- 現在はpositivity-preserving limiterを実装していません。極端な強衝撃波では、
  CFLを下げるだけでなく密度・圧力正値性を保証する制限法の追加が必要です。
- 現在の時間積分は陽的SSPRK3です。
- 単一GPU CUDA版は実装済みです。MPI+CUDA版は未実装です。
