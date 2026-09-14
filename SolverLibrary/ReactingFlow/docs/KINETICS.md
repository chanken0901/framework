# R2 詳細反応の瞬時速度評価

独立Python実装。Canteraは入力変換とテストにのみ使用し、速度評価時には呼び出さない。
この段階では反応ODEの時間積分、着火、CFD、MPI/OpenMP、GPUを実装していない。

## 対応範囲

| 形式 | 実装 |
|---|---|
| Arrhenius素反応 | A T^b exp(-Ea/RT)、可逆／不可逆、複数反応 |
| 可逆反応 | NASAの標準Gibbsエネルギーと標準圧から濃度平衡定数を導出 |
| 第三体 | 化学種ごとの衝突効率と既定効率を使用 |
| Falloff | Lindemann、Troe（3/4パラメータ）、SRI（3/5パラメータ） |
| PLOG | 対数圧力補間、同一圧力の速度和、範囲外は最近接端点 |
| Chebyshev | 温度・圧力の二次元展開、フィット範囲外は拒否 |
| 非標準反応次数 | 不可逆反応の非負次数。非反応物への次数指定も対応 |
| duplicate反応 | Canteraの入力検証後、各反応を別々に評価して生成率へ加算 |

未対応：chemically-activated、Tsang、Blowers–Masel、Linear-Burke、表面反応、
電離・プラズマ、負の反応次数、可逆反応の任意次数、PLOG以外の負A。
Troeは0<=A<=1、正のT1/T3に限定。対応外は読み込み時に反応番号・式付きで拒否し、
勝手に省略しない。R1の読み込みにもこの検査が適用される。
PLOGの負A項はグループの和が正の場合だけ許容し、評価温度で非正ならエラー。

## 単位と保存

内部濃度はmol/m³。Canteraのkmol基準係数から、反応次数nに対して
`A_mol = A_kmol * 1000^(1-n)`へ変換する。第三体の低圧係数は次数が一つ増える。
活性化エネルギーはJ/molへ変換する。Chebyshevも係数の定数項で単位を補正する。

`kinetics.evaluate(T, rho, Y)`はSI温度・密度と質量分率を受け取り、次を返す。

- `forward`, `reverse`, `net`：反応進行速度mol/(m³ s)。第三体効果込み。
- `molar_production`：化学種生成率mol/(m³ s)。
- `mass_production`：化学種生成率kg/(m³ s)。
- `heat_release`：`-sum(h_s * omega_s)`、W/m³、診断用。

圧力は同じ状態から理想気体EOSで計算する。別の矛盾する圧力入力は受け付けない。
保存エネルギーに生成エンタルピーを含めるため、heat_releaseを全エネルギー方程式へ
そのまま追加すると反応熱を二重計上する。R3/R4でエネルギー規約を維持する。

速度・質量作用則は対数で評価し、正逆差の相殺にはexpm1を使用。
ゼロ濃度、ゼロ有効第三体、不可逆反応を明示処理する。
オーバーフローは拒否し、非有限値や負密度・不正組成を補正しない。
float64以下の極小速度はアンダーフローでゼロになり得る。
Chebyshevの範囲端点はEOS丸め誤差（相対2e-14）のみ端点に合わせる。

## 検証

Cantera 3.2.0を基準として以下を実施。データの複製配布はせず、同梱資産を参照する。

- h2o2.yaml：10種29反応、400〜2800 K、10³〜10⁷ Pa。
- gri30.yaml：53種325反応、500/1200/2500 K、10⁴/10⁵/10⁷ Pa。
- 合成の各falloff形式：低圧・高圧極限を含め10⁻²〜10¹² Pa。
- PLOGの同圧力複数項・負A項・範囲外、Chebyshevの交差項・端点。
- 非標準次数、ゼロ濃度・第三体、質量／元素保存、平衡での詳細釣合い。
- Canteraなしで解析解と比較する二分子反応テスト。

速度一致は着火遅れ・燃焼速度・デトネーションの検証とは異なる。
現在の実装は基準CPU評価器であり、大規模CFDの高速化は後段。

## 実行（FrameWorkルート）

Windows PowerShell（Linux bashではpythonをpython3へ置き換える）：

```powershell
python -m pip install cantera==3.2.0 PyYAML
python -m unittest discover -s SolverLibrary/ReactingFlow/tests
python SolverLibrary/ReactingFlow/tools/inspect_thermo.py path/to/mechanism.yaml --temperature 1000 --pressure 101325 --mass-fractions '{"H2":0.1,"O2":0.9}' --rates
```

`--rates`は画面JSONに瞬時速度を追加するだけで、ファイル変更・時間積分は行わない。
機構内の名前と指定組成を一致させる。`path/to/mechanism.yaml`は実在ファイルへ変更する。

APIでは`import_cantera(path).kinetics.evaluate(T, rho, Y)`を使用する。
Yの順序は`imported.gas.names`。機構は毎セル読み直さず、一度だけコンパイルする。

式・外部入力の参照：[Cantera 3.2 rate constants](https://cantera.org/stable/reference/kinetics/rate-constants.html)。
