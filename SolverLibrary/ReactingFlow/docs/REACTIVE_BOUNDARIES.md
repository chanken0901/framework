# 反応流入境界とCFD反応波検証

最新の格子・時間・化学許容誤差比較と停止診断は[CFD反応波検証結果](CFD_WAVE_VALIDATION.md)を参照。

独立ReactingFlowのCPU逐次・1D用。計算本体・初期分布読込みはFortran。
Python/Cantera/SciPyは入力変換と独立比較試験にのみ使用する。

## 流入の設定

`left_bc='reacting_inlet'`（または`right_bc`）を指定し、その側の
`*_boundary_temperature`, `*_boundary_pressure`, `*_boundary_velocity`,
`*_boundary_y`を全て与える。単位はK、Pa、m/s、質量分率。
種の順序は機構ファイルと一致させる。既存境界は変更していない。

外向き法線nは左−1、右+1。内部状態をρ,u,p、凍結音速をaとする。
亜音速流入では、外部の温度・組成・速度を指定し、圧力は外へ伝わる
凍結音響特性との局所線形適合条件から求める：

```
Tb = Tref, Yb = Yref, ub = uref
pb = p + rho*a*n*(u-uref)
```

ρbは混合気体の状態方程式で求め、Rusanov流束に渡す。
超音速流入では参照状態全体を与える。亜音速時の参照圧力は
参照状態・音速の検証用であり、pbを固定する値ではない。
外部リザーバは一定で、反応積分は内部セルに行う。
輸送を有効にした場合は半セル距離の温度・速度・組成勾配で境界拡散流束を評価し、
境界流束込みの保存量収支に含める。境界波速・拡散係数を時間刻み制約にも含める。

これは完全な反応性NSCBCではない。逆流、参照／内部の亜音速・超音速不一致、
非正圧力、物性範囲外は停止する。自動的な流入／流出切替や無反射性能は保証しない。
流出側には既存の[局所特性流出](CHARACTERISTIC_OUTLET.md)を組み合わせられる。

## 初期分布ファイル

`&flow1d`内で`initial_profile='initial.rf'`を指定する。
相対パスは入力namelistファイルのディレクトリ基準。
指定すると従来の左右二状態初期化を置き換えるが、namelist後の左右組成2行は
従来入力との互換性のため引き続き必要。未指定なら従来通り。

```
RF_FLOW_PROFILE_V1
機構のcanonical_sha256
化学種数 セル数 領域長さ[m]
T[K] p[Pa] u[m/s] Y1 Y2 ... Yns
（左から順に全セルの値）
```

一様格子のセル中心で評価した初期値を指定する。コメントは入れない。
機構ハッシュ・種数・セル数・領域長さ・状態の妥当性・不足行・余分な記録を検査する。
異なる機構や格子のファイルを暗黙に補間・流用しない。

## 実行

`h2o2.rf`は[反応器手順](REACTORS.md)でCantera h2o2.yamlから変換する。
入力例は燃焼流入の短時間動作確認であり、検証済みデトネーション条件ではない。
出力は新しいファイル名を使う。

Windows PowerShell（FrameWorkルート）：

```powershell
.\build\reactingflow-fortran-release\rf_flow1d.exe h2o2.rf SolverLibrary\ReactingFlow\examples\flow1d_h2_reacting_inlet.in inlet.csv
$env:RF_FORTRAN_BUILD = (Resolve-Path build/reactingflow-fortran-release).Path
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -v
```

Linux bash：

```bash
./build/reactingflow-fortran-release/rf_flow1d h2o2.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_reacting_inlet.in inlet.csv
export RF_FORTRAN_BUILD="$PWD/build/reactingflow-fortran-release"
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -v
```

比較試験にはCantera、NumPy、SciPy、PyYAMLが必要。通常のFortran計算には不要。

## 検証範囲と残課題

- 左右・亜音速／超音速の一様流、定係数輸送あり／なし、逆流や参照入力不足の拒否。
- 単体試験で音響適合式の非一様状態を直接検査。
- 発熱のない一次反応A→Bと移流の解析解：40/80セルで組成誤差減少を検査。
- H2/空気・上流300 K、1 atm、Mach 5の過駆動波：独立Cantera/Radau ZNDを
  衝撃波固定座標でCFD初期値にし、短時間後の圧力分布・波面位置・保存量を比較する。
  これは実験室座標で自走するデトネーションの検証とは異なる。
- NASA温度下限で運動エネルギーを差し引く際の丸め誤差だけを許容する修正を追加。
  実質的な物性範囲外や負圧のクリッピングは行わない。

修正前の記録：粗い32/64セルでは反応帯圧力誤差が収束せず、64/128セルでは減少したものの
設定した3%基準に未達だった。128セルの圧力L1誤差/max(p_ref)は約0.03997。
256セルでは`Flow step transport/chemistry retry limit exceeded (30 attempts)`で停止した。
追加診断では物性下限直下の内部エネルギーによる棄却を確認し、
ゼロ傾き面の保存量保持・SSPRK増分形・温度復元精度を修正した。
修正後は同じMUSCL反応波試験と全86件の回帰が合格。256セルの圧力誤差は約2.12%。
一次精度の格子・時間・化学感度試験も合格。物性範囲の拡張や保存量クリッピングは行っていない。
詳細は上記の最新報告を参照。
Debug/ReleaseのFortran単体試験は各4件合格、燃焼流入入力例は正常終了。
細分化試験だけから時間・化学積分・分割誤差が
十分小さいとは判断しない。長時間の波面速度、自己維持波、CJ極限、火炎境界、
完全NSCBC、多次元については引き続き未検証／未実装であり、R4全体の完了ではない。
