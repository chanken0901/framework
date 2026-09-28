# 凍結組成の垂直衝撃波基準計算

`rf_normal_shock` はFortranの独立した基準計算窓口。上流静止・右向き衝撃波について、
NASA熱物性と質量・運動量・エネルギー保存から下流状態を求める。
化学反応は進めず、上下流の質量分率Yを固定する。**CJ速度計算、ZND反応帯、
時間発展CFDではない。R4全体の完了を示すものでもない。**

## 入力と出力

入力はSI単位。`&shock` の `temperature` [K]、`pressure` [Pa]、`mach` に続けて、
機構ファイルの化学種順の全質量分率を書く。例は `examples/shock_h2_air.in`。
Machは1.001以上。NASA範囲外は外挿せず停止する。根の探索経路が範囲外に出る場合も
拒否するため、最終状態だけが範囲内であっても強い衝撃波を拒否することがある。

CSVには衝撃波速度 [m/s]、下流温度 [K]、圧力 [Pa]、密度 [kg/m3]、
実験室系の下流速度 [m/s] を出力する。組成、機構ハッシュ、保存残差とSUCCESS行も記録する。
既存出力は上書きしない。入力組成と機構の種順を必ず合わせること。

## 保存則

上流音速をa0、衝撃波速度をD=M a0、密度比をrとすると、

```
rho1 = r rho0
p1 = p0 + rho0 D^2 (1 - 1/r)
T1 = p1 / (rho1 Rmix)
h(T1,Y) + (D/r)^2 / 2 = h(T0,Y) + D^2 / 2
u1_lab = D (1 - 1/r)
```

r=1の自明解を除き、圧縮解を探索・二分法で解く。保存残差、圧縮、エントロピー、
衝撃波固定系の下流亜音速条件を検査する。

## 実行（FrameWorkルート）

[REACTORS.md](REACTORS.md) の手順でビルドし、h2o2機構を `build/reactingflow/h2o2.rf` に変換する。
Pythonは機構の事前変換・比較テストのみで、衝撃波計算本体には不要。

Windows PowerShell:

```powershell
cmake --build build/reactingflow --parallel 4
.\build\reactingflow\rf_normal_shock.exe build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/shock_h2_air.in build/reactingflow/shock.csv
```

Linux bash:

```bash
cmake --build build/reactingflow --parallel 4
./build/reactingflow/rf_normal_shock build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/shock_h2_air.in build/reactingflow/shock.csv
```

## 回帰検証

- 定比熱gamma=1.4の解析解：Mach 1.001、1.1、2、5の温度・圧力・密度・速度と保存残差。
- Cantera熱物性とSciPyの独立根探索：H2/O2/N2混合気、Mach 1.1、2、3。
- 亜音速入力、NASA範囲超過、出力上書きの拒否。

上記は基準状態の検証。CFDの衝撃波伝播精度、CJ/ZND、格子収束の検証は別途必要。
