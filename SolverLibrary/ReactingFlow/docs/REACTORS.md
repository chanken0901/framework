# R3：硬い反応系の0次元積分

独立ReactingFlowの**Python CPU基準実装**。空間格子、流入出、壁熱伝達、
衝撃波、MPI、GPUは扱わない。既存NSEの実行環境やcase.yamlを変更しない。
Canteraは機構の入力と照合テストに使用し、時間発展中の反応速度・物性は
R1/R2の独立実装で計算する。対応機構の範囲は[KINETICS.md](KINETICS.md)と同じ。

## 方程式と単位

均一・閉鎖・断熱・中性理想混合気体。全量SI、質量分率Y、時刻s、温度K、
圧力Pa、密度kg/m³。初期T/P/Yから密度を求める。

- 共通：`dY_s/dt = omega_s / rho`（omegaはkg/m³/s）。
- `constant_volume`：rho固定、`dT/dt = -sum(e_s omega_s)/(rho cv)`。
  圧力はEOSで変化する。比内部エネルギーeを保存検査する。
- `constant_pressure`：p固定、rhoはEOSで変化し、
  `dT/dt = -sum(h_s omega_s)/(rho cp)`。比エンタルピーhを保存検査する。
  閉鎖系の質量一定に対し体積が変化する理想化（流通反応器ではない）。

`e_s=(h_molar,s-RT)/M_s`、`h_s=h_molar,s/M_s`。
生成エネルギーを含むため、R2のheat_releaseを別途加算しない。
CFDの全エネルギー式と、この温度ODEは同じ形式ではない。

## 積分・失敗時の扱い

SciPy 1.17.1の可変次数・可変刻みBDF、密な数値差分Jacobianと陰的Newton解法。
解析Jacobian、疎行列最適化、自作Fortran積分器は未実装。
状態は全化学種YとT。局所誤差の尺度は`atol_i + rtol*abs(state_i)`。
Yの絶対許容誤差とTの絶対許容誤差を分ける。

Newtonと数値Jacobianの**試行値のみ**、負のYを0として和で割った組成で
RHSを評価する。これは物理領域外の数値的延長であり、解ベクトルへ書き戻さない。
負の試行値を使った評価回数を`negative_trial_rhs_calls`に記録する。
受理ステップは元のYのまま検査し、負値（微小値も含む）、非有限値、
組成和の誤差>1e-12、NASA適用範囲外なら停止する。黙ったクリップや正規化はしない。
**BDFは正値性保証法ではない**。検査で停止した場合は許容誤差・max_stepを小さくして
再検証する。任意の機構・初期条件で必ず計算できるという保証ではない。
NASA範囲外のNewton試行も停止し、外挿しない。

各受理ステップの保存残差を出力する：

- `mass_sum_error`：max|sum(Y)-1|。
- `element_relative_error`：元素モル量/mol kg⁻¹の初期値との差をmax(1,|初期値|)で割った最大値。
- `energy_relative_error`：max|e-e0|またはmax|h-h0|を
  `max(1 J/kg, |初期保存エネルギー|, cp0*T0)`で割った値。

元素残差>1e-8、エネルギー残差>max(100*rtol,1e-7)でも停止する。
これは異常検出であり、許容誤差への収束検証を代替しない。
ステップ数上限・BDF失敗を成功扱いしない。出力ファイルは完走した場合だけ新規作成する。

## 実行手順

FrameWorkルートで実行する。開発用`run_reactor.py`を入口とする。
R7で共通run/postprocessへ統合予定であり、まだResearchRunsへ自動コピーしない。
機構ファイルはCantera同梱データを直接参照できる。機構ファイルの再配布は不要。

Windows（PowerShell、使用するPython仮想環境内）：

```powershell
python -m pip install -r SolverLibrary/ReactingFlow/requirements-reference.txt
$mechanism = python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')"
python SolverLibrary/ReactingFlow/tools/run_reactor.py "$mechanism" --temperature 1000 --pressure 101325 --mass-fractions '{"H2":0.0285,"O2":0.2264,"N2":0.7451}' --end-time 0.001 --mode constant_volume --output reactor_cv.json
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -v
```

Linux（bash、使用するPython仮想環境内）：

```bash
python3 -m pip install -r SolverLibrary/ReactingFlow/requirements-reference.txt
mechanism=$(python3 -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')")
python3 SolverLibrary/ReactingFlow/tools/run_reactor.py "$mechanism" --temperature 1000 --pressure 101325 --mass-fractions '{"H2":0.0285,"O2":0.2264,"N2":0.7451}' --end-time 0.001 --mode constant_volume --output reactor_cv.json
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests -v
```

定圧は`--mode constant_pressure --output reactor_cp.json`に変更する。
出力が既にある場合は上書きしないため、新しい名前を指定する。
組成は**質量分率**（省略種は0、合計1必須）。モル比2:1:3.76をそのまま入力しない。
数値設定の既定値はrtol=1e-7、atol-species=1e-14、atol-temperature=1e-6 K。
`--max-step`はs、既定は終了時刻、`--max-steps`は既定100000。
`--rtol 1e-9 --atol-species 1e-16 --atol-temperature 1e-8`で収束を比較できる。

JSONには機構ハッシュ、設定、種順序、全受理時刻とT/P/rho/Y、保存残差、
RHS/Jacobian/LU回数、着火遅れを格納する。等間隔出力ではない。
着火遅れは**初めてT0+400 Kを上向きに越える時刻**。
BDFの補間多項式で求根する。`--ignition-rise`で温度上昇を変更できる。
未到達はnull。最大dT/dt時刻やOHピーク時刻とは異なるので比較条件を揃える。

## 検証範囲と注意点

- Cantera 3.2.0同梱h2o2.yaml、水素/空気H2:O2:N2=2:1:3.76（モル比）、
  1000 K・101325 Pa・1 ms：定容/定圧のT/P/rho/Y履歴と着火時刻を照合。
- rtol=1e-5から1e-9への厳格化で着火時刻の誤差減少を検査。
- 質量・元素・エネルギー、無反応不変状態、解析解を持つ高速一次反応、
  不正入力・ステップ上限・負の受理状態の拒否を検査。
- 同機構のNASA切替点1000 Kでは初期組成のhに約-0.139543 J/kgの不連続があり、
  相対エネルギー残差が約9e-8で頭打ちになる。機構の係数は変更しない。
  1100 K開始の滑らかな範囲でも別途エネルギー誤差の収束を確認する。
- 大規模炭化水素機構の着火、長時間平衡、全圧力依存形式の時間発展、
  デトネーションの検証はこのR3テストの範囲外。速度単体のR2検証とは区別する。

理論・積分器の参照：
[Cantera定圧理想気体反応器](https://www.cantera.org/3.2/reference/reactors/ideal-gas-constant-pressure-reactor.html)、
[SciPy BDF](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.BDF.html)。
