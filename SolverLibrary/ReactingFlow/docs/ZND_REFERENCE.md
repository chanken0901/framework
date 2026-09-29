# 指定速度ZND基準計算

`rf_znd` は平面・定常・非粘性・断熱の反応帯をFortranで積分する独立窓口。
上流静止、右向き衝撃波を仮定する。`mach` から衝撃波速度Dを決め、
凍結衝撃波直後の状態を初期値に詳細反応を積分する。
これは時間発展CFDではなく、CJ速度の自動探索でもない。
追加の `speed_mode='cj'` と終端平衡検査は [CJ_REFERENCE.md](CJ_REFERENCE.md) を参照。
既定の `speed_mode='mach'` では、指定速度がCJ以上かの判定はしない。

## 方程式と座標

独立変数tauは衝撃波を通過した流体粒子の滞留時間[s]。
xは衝撃波から下流へ測った距離[m]で、実験室系の右向き座標とは異なる。
uは正の衝撃波固定系流速。実験室系流速はD-u。

```
j = rho0 D                  (一定質量流束)
B = D + p0/j
rho = j/u
p = j(B-u)
T = u(B-u)/R(Y)
dYk/dtau = omega_k/rho
dR/dtau = sum(Rk dYk/dtau)
Hdot = sum(hk dYk/dtau)
a_f^2 = (cp/cv) R T
eta = 1 - u^2/a_f^2
du/dtau = u [ (dR/dtau)/R - Hdot/(cp T) ] / eta
dx/dtau = u
```

hkは生成エンタルピーを含む種別比エンタルピー[J/kg]。
質量・運動量は代数式で保存し、全エンタルピーh+u²/2と元素量は各受理ステップで検査する。
最大質量分率の種を従属変数として質量分率和を1に保つ。
積分は既存DVODEのBDF・数値差分密Jacobian。
Newton試行の負組成はRHS評価だけ非負化・再正規化するが、受理解はクリップせず検査する。

物理モデルの参考：[Caltech SDToolbox](https://shepherd.caltech.edu/EDL/PublicResources/sdt/)。
SDToolboxコードは組み込んでおらず、今回の比較先もSDToolboxそのものではない。

## 入力・失敗時の扱い

`&znd` と、その後に機構の種順の質量分率を書く。単位はSI。
`examples/znd_h2_air.in` はh2o2.yaml種順専用。

- `temperature`, `pressure`, `mach`：上流条件と指定速度。
- `end_time`, `max_step`, `max_steps`：滞留時間上限・内部刻み上限・受理ステップ上限。
- `rtol`, `atol_species`, `atol_velocity`, `atol_distance`：積分許容誤差。
- `sonic_margin`：etaの下限。既定1e-6、許容1e-8以上1未満。

音速点付近、NASA範囲外、負組成、積分失敗、保存誤差過大では停止する。
Newton試行が範囲を外れた場合も停止し、汎用的な回復処理はまだない。
途中CSVが残ることがあるため末尾のSUCCESS行を確認すること。
既定のSUCCESSは指定滞留時間への到達だけを意味し、平衡・CJへの到達ではない。
`equilibrium_tolerance>0` の場合だけ、追加の終端平衡組成検査にも合格したことを示す。
既存出力は上書きしない。粘性・熱伝導・拡散・曲率・ノズルはこの基準モデルに含めない。

## 実行（FrameWorkルート）

[REACTORS.md](REACTORS.md) のビルド・機構変換を先に行う。

Windows PowerShell:

```powershell
cmake --build build/reactingflow --parallel 4
.\build\reactingflow\rf_znd.exe build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/znd_h2_air.in build/reactingflow/znd.csv
```

Linux bash:

```bash
cmake --build build/reactingflow --parallel 4
./build/reactingflow/rf_znd build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/znd_h2_air.in build/reactingflow/znd.csv
```

CSVは滞留時間・距離・温度・圧力・密度・衝撃波固定系流速・eta・全種質量分率。
機構ハッシュ、指定速度、積分設定、エネルギー／元素保存誤差も記録する。

## 検証範囲と残作業

Cantera3.2の熱物性・反応速度とSciPy Radauを使う独立比較試験を追加した。
比較側は圧力・密度を微分する形式で、Fortran側のuによる代数消去形式と区別した。
純N2、Mach3の非反応状態保持、およびH2:O2:N2=2:1:3.76、300K、101325Pa、
Mach5、滞留時間1usの反応帯を比較する。rtol=1e-8と1e-10で照合し、
有限区間の状態、保存量、音速限界の拒否、上書き防止を検査する。
入力例の組成は丸めた質量分率であり、比較試験の厳密なモル比とはわずかに異なる。

その後、限定条件で終端平衡までの比較とCJ速度探索を追加した（上記CJ_REFERENCE参照）。
誘導長の十分な収束・CJ近傍特異点・CFDとの比較はまだ未検証。
任意条件のデトネーションを計算可能とする保証ではなく、R4は継続中。
