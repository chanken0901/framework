# 過駆動反応波の下流境界距離依存性

## 検証結果（2026-10-01）

Release版の独立実行で合格（870.3秒）。判定基準の緩和なし。

| 共通の0〜2 mmにおける定常基準との差 | 2 mm/512セル | 3 mm/768セル |
|---|---:|---:|
| 圧力誤差 | 3.78300% | 3.78235% |
| 温度誤差 | 2.56845% | 2.56824% |
| H2O質量分率の平均絶対差 | 0.00502753 | 0.00502758 |

両解間の圧力差／基準最大圧力は1.60660e-5（0.00160660%）、
温度差／基準最大温度は4.26303e-6（0.000426303%）、
H2O質量分率の平均絶対差は5.63329e-8。
後半の波面移動速度は両者とも約−3.93636 m/s。
位置・検出・保存量・正圧を含む全基準に合格した。
今回の条件・時刻では境界距離感度が小さいが、下記の音響伝播時間に関する制限がある。

成果物：`build/reactingflow-boundary-distance-wave/`。
通常回帰90件中88件合格・長時間系2件スキップ（本試験は別実行）、
Fortran CTestはRelease/Debugそれぞれ5/5合格。

## 目的と比較条件

Fortran `rf_flow1d` のMUSCL反応波計算について、流出境界を遠ざけた影響を調べる。
検証用Pythonは独立Cantera/Radau基準の生成・実行・比較だけを担当する。
ソルバーの支配方程式や数値方式は変更しない。

| 条件 | 基準領域 | 下流延長領域 |
|---|---:|---:|
| 領域長 | 2 mm | 3 mm |
| セル数 | 512 | 768 |
| 格子幅 | 3.90625 µm | 3.90625 µm |
| 初期衝撃波位置 | 0.5 mm | 0.5 mm |
| 終了時刻 | 1.0 µs | 1.0 µs |

H2/空気、上流300 K・1 atm、上流Mach 5の衝撃波固定座標系。
左は同じ反応流入、右は各境界位置での定常ZND参照状態を用いる局所線形特性流出。
最大刻み2 ns、化学相対許容誤差1e-9、種絶対許容誤差1e-16、温度絶対許容誤差1e-8。
初期分布の共通部分が一致し、格子点が重なることも検査する。

## 判定

比較対象は両ケースの **共通の0〜2 mm**。領域を広げたことで平均誤差が
見かけ上小さくなるのを避けるため、3 mm全体での誤差を合否に使わない。

- 各解と定常基準の圧力・温度平均絶対差／基準最大値：それぞれ5%未満。
- 各解と定常基準のH2O質量分率の平均絶対差：0.02未満。
- 両解間の圧力・温度平均絶対差／基準最大値：それぞれ1%未満。
- 両解間のH2O質量分率の平均絶対差：0.005未満。
- 両ケースで波面を全記録時刻で検出し、最終位置偏差は0.1 mm以下。
- 後半の波面移動速度の絶対値／上流速度：5%未満。
- 正圧、正常終了、境界流束込みの質量・運動量・エネルギー・元素保存誤差1e-8未満。

これは上記条件での境界距離感度の検査であり、無反射性の証明ではない。
過駆動波の有限時間検証であり、自走CJ波・任意機構・長時間安定性を保証しない。
時間刻みはCFL制限も受けるため、同じ最大刻みでも両領域でステップ列が厳密に一致するとは限らない。

特に、2 mm領域の初期定常分布から凍結音速を用いて評価した
上流向き音響特性の速さ `a-u` は約276〜548 m/s。
右境界から初期波面までの移動時間を `Σ dx/(a-u)` で見積もると約4.60 µsとなる。
したがって1.0 µsの一致だけでは、境界で生じた擾乱が波面に戻った後の影響を検証できない。
この時間は初期分布に基づく概算で、反応・波形変化を含む厳密な伝播時間ではない。
今後は音響伝播時間を越える計算で再確認する必要がある。

## 再実行

FrameWorkルートで、Release版をビルドしNumPy/SciPy/Cantera/PyYAMLが利用できるPythonを使う。
通常回帰ではスキップし、明示的に有効化する。
成果物保存先は未使用の名前を指定する（既存結果を上書きしない）。

PowerShell:

```powershell
$env:RF_FORTRAN_BUILD = (Resolve-Path build/reactingflow-fortran-release).Path
$env:RF_BOUNDARY_WAVE_TESTS = '1'
$env:RF_WAVE_ARTIFACTS = "$PWD/build/boundary-distance-run01"
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -k boundary_distance -v
Remove-Item Env:RF_BOUNDARY_WAVE_TESTS
Remove-Item Env:RF_WAVE_ARTIFACTS
```

Linux:

```bash
RF_FORTRAN_BUILD="$PWD/build/reactingflow-fortran-release" \
RF_BOUNDARY_WAVE_TESTS=1 \
RF_WAVE_ARTIFACTS="$PWD/build/boundary-distance-run01" \
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -k boundary_distance -v
```

各ケースの機構・入力・初期分布・全場出力・波面履歴・ソルバーログを保存する。
