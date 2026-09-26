# 1D特性流出境界（R4 集約作業1）

`left_bc` または `right_bc` に `characteristic` を指定する。
同じ側の `*_boundary_temperature`、`*_boundary_pressure`、`*_boundary_velocity`、
全種の `*_boundary_y` で外部参照状態を指定する。既存のDirichlet入力と同じ形式。
反対側は独立に選べるが、周期境界との片側だけの組合せは禁止。

これは局所音響特性を用いた流出専用の外部状態構成であり、完全な反応性NSCBCではない。
特性波を使う境界条件の背景は
[Poinsot・Lele (1992)](https://ntrs.nasa.gov/citations/19920065566) を参照。
本実装が同論文の方式全体を再現するという意味ではない。

## 状態の構成

外向き法線nを左端−1、右端+1、セル側状態をρ,u,p、凍結組成音速をaとする。
外向き速度un=n*uについて、un>=aなら全波が流出するためセル状態をそのまま使う。
0<=un<aでは外部参照の圧力・法線速度から

```
A = ((pref-p) - rho*a*(un_ref-un))/2
pb = p + A
rhob = rho + A/a**2
ub = u - n*A/(rho*a)
Yb = Y
Tb = T*(pb/p)*(rho/rhob)
```

を構成し、既存のRusanov流束へ渡す。内部からの出ていく音響成分と
組成を保ち、入ってくる音響成分を外部状態で指定する局所線形化である。
参照温度・組成は参照状態の妥当性検証に使い、流出組成を上書きしない。
境界の拡散流束は0（既存outflowと同じ）、境界隣接セルのMUSCL傾きは片側外挿により0。
生成した外部状態の波速も時間刻み判定に含める。

## 制限

- セル側または生成した境界状態が逆流する場合は明示的に停止。流入への自動切替はしない。
- 非正圧力・非正密度・NASA物性範囲外は停止。外部状態をクリップしない。
- 小振幅音響の検証までであり、強い反応波、亜音速／超音速遷移、デトネーションの通過、
  境界近傍の強い粘性・熱伝導・種拡散について実用検証済みではない。
- 「厳密無反射」や、従来境界に常に勝ることを意味しない。

## 検証と実行例

`rf_flow_unit` に、一様状態、超音速流出、64セルの小振幅音響パルス流出、
境界流束を含めた全保存量収支の試験を追加した。
入射圧力振幅10 Paに対し、終了時の流入向き成分のRMSは約1.36e-4 Pa。
比較した固定外部状態＋Rusanovも約9.80e-5 Paで、この試験では本方式の優位性は示されていない。
合格基準は流入成分RMSが1e-3 Pa未満、保存量相対収支誤差1e-10未満。
左右の入力・一様流出、参照状態欠落、逆流拒否は実行ファイルの回帰試験で確認する。

`examples/flow1d_h2_characteristic.in` はh2o2の10種順の一様窒素流の例。
機構をh2o2.rfへ変換済みなら、リポジトリ直下で以下を実行する。出力名は未作成のものを使う。

```powershell
.\build\reactingflow-fortran-release\rf_flow1d.exe h2o2.rf SolverLibrary\ReactingFlow\examples\flow1d_h2_characteristic.in characteristic.csv
```

```bash
./build/reactingflow-fortran-release/rf_flow1d h2o2.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_characteristic.in characteristic.csv
```
