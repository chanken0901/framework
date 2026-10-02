# 多次元有限体積基盤（実装先行・未検証）

2026-10-02。長時間反応波試験を保留し、R4の多次元・形状対応の基盤を先行実装する。
6 µsの不合格を解決済みとはしない。これは実用2D/3D反応流ソルバーの完成ではない。

## 今回の範囲

`src/fortran/mod_rf_finite_volume.f90` はCPU逐次のFortranモジュール。
既存の `rf_flow1d` の入力・実行・数値方式は変更しない。
1〜3次元の状態変換、任意方向のRusanov面流束、すべり鏡像壁、
静止メッシュの一次精度非粘性残差と保守的なCFL刻み評価を提供する。
NASA熱力学・温度反転・許容範囲は既存1Dの実装を共有する。

保存変数は `q=[rho*Y(1:ns),rho*u(1:nd),rho*E]`。
全エネルギーは生成エネルギーと全速度成分の運動エネルギーを含む。単位系はSI。
化学種数は機構から決まり固定上限を追加しない。

## 面・セルの契約

`rf_face_mesh` は `owner(nface)`, `neighbor(nface)`,
`area_vector(nd,nface)`, `volume(ncell)` を所有する。
面積ベクトルSはownerからneighborへ向く。境界面はneighbor=0、Sは外向き。
内部面は一度だけ登録する。セル体積は正、各セルの外向き面積ベクトルの和はゼロ。
検査は形状、有限値、接続添字、正の面積・体積と
`norm(sum(S))/sum(norm(S)) <= 1e-12` を確認する。
幾何形状の自己交差や重複面などを完全に検出するメッシュ品質検査ではない。
1DのSは断面積、2Dでは単位奥行きの面長、3Dでは面積として扱う。

面積込み物理流束は次のとおり。

```text
F_species = rho*Y * (u dot S)
F_momentum = rho*u * (u dot S) + p*S
F_energy = (rho*E+p) * (u dot S)
lambda = abs(u dot S) + sound_speed*norm(S)
Fhat = (F_left + F_right - max(lambda_left,lambda_right)*(q_right-q_left))/2
```

`finite_volume_rhs` はownerからFhatを引き、neighborへ同じFhatを加え、
最後に各セル体積で割る。`boundary_rate` は外向き境界流束の負の総和。
従って `sum(volume*dq)=boundary_rate` となる。
入力ghostは `(ns+nd+1,nface)`、内部面の列は参照しない。
境界状態は呼出側で設定する。鏡像壁には `reflect_normal` を使用できる。
この壁は法線運動量のみ反転し、接線運動量・化学種・全エネルギーを保持する。
粘着壁・熱伝達壁・触媒壁ではない。

CFL刻みは `cfl*min(volume/sum_face(lambda))`、`0<cfl<=1`。
化学・分子輸送の刻み制約や高次時間積分の安定性を保証しない。
1Dの既存CFL評価とは定義が異なり、一様格子ではより保守的になる。

## API

- `primitive_nd` / `conserved_nd`: 状態変換。後者の任意 `ok` で不適合状態を返せる。
- `normal_flux` / `rusanov_normal_flux`: 面積込み流束と面積込み特性速度。
- `reflect_normal`: すべり鏡像ghost生成。
- `validate_face_mesh`: 静止メッシュ契約検査。
- `finite_volume_rhs`: セル残差、境界収支、流体CFL刻み。状態は変更しない。

## 次に実装する部分

1. [2D構造格子・ノズル生成](GRID2D.md)と境界面分類は追加済み（ビルド確認のみ）。
2. [2D時間積分・段棄却・反応結合と実行入口](FLOW2D.md)を追加済み（ビルド確認のみ）。
3. 多次元MUSCL、分子輸送、特性流入出境界。
4. 実行検証後のMPIペンシル・GPU常駐移行（R5/R6）。

曲面形状を面積ベクトルで扱う土台に、2D時間発展の実験的入口を追加した。
粘性・高次化を含む実用一般座標ソルバーはまだ未完成。ノズル形状生成と時間発展の追加範囲は
上記リンクを参照。移動格子は対象外。
検証コード `rf_finite_volume_unit` を追加したが、今回は実行せずビルド確認のみ。
長時間試験・通常回帰も再開していない。
