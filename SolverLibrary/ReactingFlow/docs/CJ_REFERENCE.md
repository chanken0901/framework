# 化学平衡とCJ速度の基準計算

`rf_cj` はFortranで理想中性気体の平衡CJ速度・下流状態・組成を求める。
反応速度を長時間積分して平衡を近似するのではなく、NASA熱物性と元素保存で平衡を解く。
Python/Canteraは機構変換と独立比較試験だけに用い、計算本体には不要。

## 方法・仮定

温度T・密度rho固定で、元素ポテンシャルlambdaからモル量n_i [mol/kg]を求める。

```
n_i = p_ref_i/(rho R_u T) exp[-g_i^0/(R_u T) + sum_e a_ei lambda_e]
sum_i a_ei n_i = b_e
Y_i = W_i n_i
```

存在しない元素を含む種は候補から除外する。減衰Newton法で元素残差を解く。
微量種による悪条件化にはNewton行列の正則化を使うが、収束判定は元の元素保存残差
（最大相対値1e-11未満）で行う。従属元素などで行列が特異なら停止し、未検証の解を返さない。

各密度比r=rho/rho0について、平衡Hugoniotの高圧枝を温度二分法で求める。

```
e - e0 = (p+p0)(1/rho0 - 1/rho)/2
D^2 = (p-p0) / [rho0 (1-1/r)]
```

指定密度比区間を21点で走査し、内部最小と単峰性を確認して黄金分割探索する。
平衡Hugoniot上のDの最小値をCJ候補とし、エネルギー保存、エントロピー増大、
上流超音速、`D/r = a_equilibrium` を確認する。平衡音速は平衡状態のT・rho微分から
`a_eq² = p_rho - p_T s_rho/s_T` として計算し、凍結音速と混同しない。
NASA範囲外・探索区間端の最小・非単峰区間は拒否する。

方法の参考：[Caltech SDToolboxのCJ/Hugoniot解説](https://shepherd.caltech.edu/EDL/PublicResources/sdt/nb/sdt_intro.slides.html)。
上流コードを組み込んだものではなく、比較もSDToolboxそのものではなくCantera平衡＋SciPyで行っている。

平衡候補は機構内の全許容気相種。固体・液体・電離・非理想EOSは対象外。
速度機構の不可逆性や切断された反応経路で、その熱力学平衡に到達できるとは限らない。
任意機構・任意条件での大域的な最小値や、実験の非理想デトネーション速度を保証しない。

## 入力と実行

`&cj` に上流 `temperature` [K]、`pressure` [Pa]、`ratio_min` と `ratio_max` を指定し、
次行に機構の種順で質量分率を書く。密度比の既定区間は1.2〜2.5。
区間を広げるとNASA範囲外になる場合があるので、探索拒否を無理に回避せず状態を確認する。

[REACTORS.md](REACTORS.md) の手順でビルドとh2o2機構変換を行った後、FrameWorkルートで実行する。

Windows PowerShell:

```powershell
cmake --build build/reactingflow --parallel 4
.\build\reactingflow\rf_cj.exe build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/cj_h2_air.in build/reactingflow/cj.csv
.\build\reactingflow\rf_znd.exe build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/znd_cj_overdriven_h2_air.in build/reactingflow/znd_full.csv
```

Linux bash:

```bash
cmake --build build/reactingflow --parallel 4
./build/reactingflow/rf_cj build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/cj_h2_air.in build/reactingflow/cj.csv
./build/reactingflow/rf_znd build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/znd_cj_overdriven_h2_air.in build/reactingflow/znd_full.csv
```

既存出力は上書きしない。CJ出力は速度[m/s]、温度[K]、圧力[Pa]、密度[kg/m3]、全種Y。
機構ハッシュ、上流状態・組成、保存残差も記録する。

## ZNDとの接続

`speed_mode='mach'` は従来の直接Mach指定を維持する。
`speed_mode='cj'` では同じ機構・上流条件でCJを解き、`D=overdrive*D_CJ` とする。
ここでoverdriveは**速度比**であり、文献で使うことがある二乗比ではない。
1以上を指定する。CJモードではmach入力は使用せず計算で置き換える。
`ratio_min/max` はCJ探索へ渡す。

`equilibrium_tolerance>0` を指定すると終端で同じ元素量のTV平衡を別途解き、
`max |Y_end-Y_equilibrium(T_end,rho_end)|` が許容値以下か検査する。
0（既定）ならこの検査をしない。非平衡なら途中CSVを残してエラーとなりSUCCESSを書かない。
これは指定組成許容誤差の判定であり、全状態・全化学種反応速度の完全停止の証明ではない。

`overdrive=1` のCJ極限は特異点付近の積分が難しいため、全条件での安定到達は未保証。
まず例の1.05（過駆動）で検証する。音速点拒否は維持し、到達を成功と誤認しない。
最大発熱率とその滞留時間・距離を受理ステップから記録する。連続補間による極値ではないため、
誘導時間・長さとして使う場合はmax_stepと許容誤差を変えて収束確認すること。

## 検証済み範囲

- 平衡TV：H2/O2/N2、H2/O2、CH4/O2/N2、500/1500/3000K、0.1/3 kg/m3をCanteraと比較。
- CJ：H2:O2:N2=2:1:3.76、300K、101325/202650PaをCantera平衡＋SciPy根探索・最小化と比較。
  平衡等エントロピー差分による独立な音速条件も確認。
- ZND：Mach5、100usまでの反応帯をCantera速度＋Radauの別形式のODEと比較。
  rtol=1e-8/1e-10と最大刻みの半減で比較し、終端平衡組成の許容誤差1e-6を検査。
- CJ連動D/D_CJ=1.05の終端平衡検査、探索区間不適合、未平衡の短時間停止、上書き拒否。

R4のCFD格子収束、反応境界、多次元・一般座標・ノズル、詳細輸送の検証は別作業であり、
この基準計算の成功だけではR4全体を完了としない。

## 添付例の実行確認値

Cantera3.2.0付属h2o2.yamlを変換した機構と、添付例の丸めた質量分率での確認値：

- CJ速度：約1976.008 m/s、下流温度：約2964.365 K。
- D/D_CJ=1.05、100usのZND：最大エネルギー相対誤差4.63e-10未満、元素相対誤差1.64e-15未満。
- 終端の平衡組成最大絶対差：約2.23e-16以下。

これらは添付例の動作確認値であり、任意条件の誤差保証ではない。
Debug/Releaseで77件の比較・回帰試験と4件のFortran試験に合格した。
