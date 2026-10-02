# 6 µs反応波試験：未達と診断

## 2026-10-02の追加結果

元のMach 5条件について、次の試験はいずれも6 µsまで正常終了したが定常基準に不合格。
基準変更やMach 6への置き換えによって合格扱いにはしない。

| 条件 | 圧力相対平均誤差 | 温度相対平均誤差 | H2O L1誤差 | 棄却回数 |
| --- | ---: | ---: | ---: | ---: |
| 1024セル、max_dt=2e-9 s | 0.363041 | 0.527572 | 0.138882 | 449 |
| 512セル、max_dt=3e-10 s | 0.363264 | 0.526780 | 0.138883 | 76,342 |

前者の実測最大刻みは3.18534e-10 sで、max_dtではなくCFLで制限される。
証跡は `build/reactingflow-wave-cause-fine/` と
`build/reactingflow-wave-cause-small-step/` に保存。
倍の格子数、小さい刻みともに解消しないが、これだけで物理的不安定性と断定できない。

### 棄却理由の短時間再現

小刻み試験の途中の基本変数スナップショットを初期条件として50 ns再計算した。
`build/reactingflow-retry-replay-long/10000/` と `20000/` は各520回の棄却で正常終了。
ログでは輸送中のセル123において内部エネルギーがNASA下限を
約5.59958e-8 J/kg下回り、刻み半減を繰り返すことを確認した。
これは丸め誤差近傍の下限判定を調べる具体的な再現例であり、元の76,342回全ての分類でも、
波の弱化の根本原因の確定でもない。基本変数からの再構成なのでビット単位の再始動でもない。

`&flow1d` に `diagnose_retries=.true.` を追加すると、各棄却の化学半段・輸送・
最終状態判定の情報と試行番号・刻みを標準出力へ出す。既定値は `.false.`。
ログ量が増えるため通常計算には不要。30回目の失敗診断は従来どおり残る。
成功試行を棄却として出さないよう、ログは採用判定の後に置く。
診断の有無で保存変数・境界流束・採用刻み・棄却回数が完全一致する単体試験を追加した。
更新後のRelease/DebugのFortran単体試験は各5件合格。全Python回帰は今回未再実行。

2026-10-01。[境界距離比較](BOUNDARY_DISTANCE_VALIDATION.md)を同じ格子幅で6 µsへ延長した。
**ソルバーは正常終了したが、定常反応波の検証基準には不合格。**
短時間の合格でこの結果を上書きせず、R4の未解決事項として扱う。

## 確認結果

誤差の評価区間は両ケース共通の0〜2 mm。

| 条件 | 圧力誤差 | 温度誤差 | H2O質量分率L1差 | 全領域の最終衝撃波候補 |
|---|---:|---:|---:|---:|
| 2 mm / 512セル | 34.914% | 52.082% | 0.138883 | 1.37891 mm |
| 3 mm / 768セル | 34.080% | 51.936% | 0.138883 | 1.38281 mm |

初期位置は0.5 mm。両ケースとも約3.8063 µsに探索範囲の上端0.8 mmを越えた。
最終候補は最終全場の最大隣接圧力差から別途調べた位置であり、追跡履歴の継続値ではない。
両解間の共通区間の圧力差は基準最大値の1.06483%、温度差は0.179767%、
H2O平均絶対差は4.93907e-8。圧力差の1%基準も未達。
保存誤差は最大約1.23e-14、最終正圧・終了時刻6 µs・正常終了を確認した。

同じ波面移動が下流を延長した領域でも発生しており、2 mmの右境界だけを原因とは断定できない。
格子解像度、数値流束／化学分割誤差、境界の定常波保持、物理的な反応波非定常性の
切り分けが必要。現段階で「物理的不安定だから問題ない」「実装バグが原因」とは断定しない。

## 検証コードの修正

未検出レコードの `position=0` は無効値だが、旧検証コードはこれも直線回帰に含めたため、
後半の速度を約−280 m/sと誤って報告した。これは物理的な波面速度ではない。
今回、対象区間に未検出が1件でもあれば速度を無効（NaN）とするよう修正し、単体試験を追加した。
元から全時刻での検出が合否条件なので、旧処理でもこの試験が合格することはない。
Fortranの検出器と探索範囲、物理的な合格基準は変更していない。

## 再実行

通常回帰ではスキップする。必要パッケージはNumPy/SciPy/Cantera/PyYAML。
各ケースに最大14400秒の実行タイムアウトを設ける。細分化した格子は実行に時間がかかる。
タイムアウト時も、出力途中のファイル・入力・標準出力を保存してから試験を失敗させる。
長時間試験は1000ステップごとに全場も保存する。途中出力は正常終了を意味しない。
既存結果を上書きしないため、保存先は新しい名前を使用する。
保存先は計算開始前に確保し、重複していれば長時間計算を起動せず拒否する。
保存を指定した試験は入力も保存先へコピーし、そのディレクトリで結果を生成するため、
実行途中の出力も残る。Pythonのプライベート一時フォルダの存続・読み取り権限には依存しない。

PowerShell（FrameWorkルート）:

```powershell
$env:RF_FORTRAN_BUILD = (Resolve-Path build/reactingflow-fortran-release).Path
$env:RF_ACOUSTIC_WAVE_TESTS = '1'
$env:RF_WAVE_ARTIFACTS = "$PWD/build/acoustic-return-run01"
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -k acoustic_return -v
Remove-Item Env:RF_ACOUSTIC_WAVE_TESTS
Remove-Item Env:RF_WAVE_ARTIFACTS
```

Linux:

```bash
RF_FORTRAN_BUILD="$PWD/build/reactingflow-fortran-release" \
RF_ACOUSTIC_WAVE_TESTS=1 RF_WAVE_ARTIFACTS="$PWD/build/acoustic-return-run01" \
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -k acoustic_return -v
```

現在の記録条件では不合格の再現を期待する試験であり、実用計算例として推奨しない。

## 記録の版

結果は `build/reactingflow-r4-acoustic-wave/` に保存。
実行バイナリは `build/reactingflow-r4-acoustic-bin/rf_flow1d.exe`、SHA256:
`B0561C2430010BE2EC441C2349CB6894F0773BF1625C8C35AEB827F0FF31E13D`。
温度反転Newton化後、試行段の範囲違反回復処理追加前のバイナリである。
この6 µs結果を回復処理追加後の検証済み結果と混同しない。

## 原因切り分けの追加結果（2026-10-01）

回復処理追加後の現行Releaseバイナリでも、Mach 5・2 mm/512セルを再実行した。
今回の対照試験に用いた `rf_flow1d.exe` のSHA256は
`3B65F9DFCD107E7836E5BE91DDDE9F018DF9769B72E257394F1DBF1AA5EF4E0F`。
圧力誤差0.34914120941743165、温度誤差0.5208161097180227、
H2O平均絶対差0.1388830615627946、初回未検出時刻3.8063017087126695 µsであり、
従来の不合格が再現した。したがって、範囲違反時の回復処理はこの問題の解決ではない。

| 切り分け条件 | 6 µs時点の結果 | 解釈の範囲 |
|---|---|---|
| Mach 5、反応なし、凍結衝撃波、128セル | 圧力誤差0.4953%、波面0.484375 mm、合格 | 流体・境界処理だけでは同規模の波面流出は再現しない |
| Mach 6、反応あり、2 mm/512セル | 圧力誤差0.3374%、温度誤差0.3776%、H2O差0.00130119、波面0.5 mm、合格 | より大きいオーバードライブでの独立した対照例。Mach 5の合格を意味しない |
| Mach 6、反応あり、3 mm/768セル、共通0〜2 mm区間で評価 | 圧力誤差0.3608%、温度誤差0.2874%、H2O差0.000563092、波面0.496094 mm、合格 | 同じ格子幅で下流だけ延長した対照例 |

Mach 6の後半区間の回帰波面速度は−0.02529 m/s、全時刻で波面を検出した。
3 mm領域でも全時刻検出し、後半の回帰速度は−0.07739 m/s。
両計算を共通0〜2 mmで比較した差は、圧力0.53657%、温度0.50695%、
H2O平均絶対差0.00147702であり、元の境界距離比較基準（1%、1%、0.005）を満たした。
初期場による2 mm領域の下流境界から波面への音響伝播時間の概算は2.6234 µsで、
6 µsはその時間を越える。ただし、任意の波形・条件に対する無反射性の保証ではない。
結果はそれぞれ `build/reactingflow-wave-cause-frozen/`、
`build/reactingflow-wave-cause-overdrive/`、現行Mach 5再現は
`build/reactingflow-wave-cause-coarse-current/` に保存した。
Mach 6・3 mm領域の結果は `build/reactingflow-wave-cause-overdrive-far/` に保存した。

同じ機構・上流組成のCJ速度は1976.318978 m/s。
Mach 5の速度2043.870865 m/sに対応するオーバードライブは
`f=(D/DCJ)^2 ≈ 1.0695`。ZND定常解が存在することと、その時間発展の安定性は別問題である。
低オーバードライブでの水素・空気の脈動・離脱は報告されているが、
文献の機構・圧力条件は本試験と同一ではなく、文献だけで今回の原因を断定しない。
参考：[Daimon and Matsuo (2007)](https://doi.org/10.1063/1.2801478)。

初期反応帯全セルでFortran反応速度とCanteraを比較した代表最大値規格化差は、
512/1024セルでそれぞれ3.17e-14/3.25e-14。
初期物理流束の相対変動は質量約1e-15、運動量約7e-13、エネルギー約4e-12だった。
初期解の反応速度の単位や保存流束の不整合を示す結果ではない。
一方、衝撃波直近と端点を除くH2離散残差の絶対値積分／化学源絶対値積分は
512セルで1.437%、1024セルで0.213%であり、離散化誤差の解像度依存性は存在する。
この初期残差の改善だけで6 µsの格子収束を証明したことにはならない。

### 衝撃波弱化後の着火遅れ

Mach 5・512セルの最終場では最大温度1284.47 K、最大H2O質量分率3.42e-6となり、
計算領域内の反応生成物がほぼ失われていた。圧力勾配の探索失敗だけでは説明できない。
初期の衝撃波直後と、6 µs時点の波面の少し下流から温度・圧力を採り、
新気組成の定容断熱反応を独立に比較した。

| 代表状態 | 温度 K | 圧力 Pa | Fortranによる200 K昇温時間 |
|---|---:|---:|---:|
| Mach 5・初期衝撃波後 | 1620.9502 | 3002967.01 | 0.348136 µs |
| Mach 5・弱化後 | 1279.79952 | 2171654.72 | 23.2998 µs |
| Mach 6・初期衝撃波後 | 2154.55077 | 4358041.49 | 0.0642809 µs |

`test_postshock_ignition_histories` でFortranとCanteraの温度・圧力・密度・全化学種の
時系列を比較し合格した。Fortranが記録した昇温時刻におけるCantera温度も、
初期温度+200 Kに0.02 K以内で一致することを検査する。
弱化後の着火遅れは初期の約67倍であり、反応が波面を支えにくくなる方向と整合する。
ただしこれは代表状態の**0次元定容試験**であり、CFD内の着火位置・遅れを直接測定した値ではない。
この一致は孤立した化学計算の確認であり、CFDの分割誤差や最初の波面弱化の原因を除外するものではない。

### 対照試験の再実行

`test_frozen_stationary_shock` は通常回帰に追加した。
元のMach 5試験は条件・合格基準を変更せず残している。
細分化Mach 5と高オーバードライブMach 6は、独立した任意試験として実行する。
以下はMach 6の2 mm/512セルと3 mm/768セルを同じ0〜2 mm区間で比較する例。
再実行時は保存先名を変更する。

```powershell
$env:RF_FORTRAN_BUILD = (Resolve-Path build/reactingflow-fortran-release).Path
$env:RF_WAVE_CAUSE_TESTS = '1'
$env:RF_WAVE_ARTIFACTS = "$PWD/build/wave-cause-run01"
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -k znd_overdrive_acoustic_return -v
Remove-Item Env:RF_WAVE_CAUSE_TESTS
Remove-Item Env:RF_WAVE_ARTIFACTS
```

```bash
RF_FORTRAN_BUILD="$PWD/build/reactingflow-fortran-release" \
RF_WAVE_CAUSE_TESTS=1 RF_WAVE_ARTIFACTS="$PWD/build/wave-cause-run01" \
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests -p test_reactive_boundaries.py -k znd_overdrive_acoustic_return -v
```

Mach 5・1024セルを調べる場合は、`-k` の値を `znd_fine_acoustic_return` に変更する。
試験のタイムアウト延長や診断の追加は、数値解法自体の修正ではない。
