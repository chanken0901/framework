# R1〜R3：Fortran計算本体と0次元反応器

**計算本体はFortran**。NASA物性、混合気体EOS、温度反転、基準量変換、
詳細反応速度、定容／定圧の時間積分を`src/fortran/`で実装する。
NSE、GPE、Pythonインタープリターへの実行時依存はない。
空間離散化、流入出、壁熱伝達、衝撃波、MPI、CUDAはまだ未実装。

## 構成と依存

- `reactingflow`：Fortran静的ライブラリ。種数・反応数は可変長。
- `rf_reactor`：Fortranの0次元計算実行ファイル。
- `rf_probe`：Fortran物性／速度照合用、`rf_unit`：Fortran単体検証。
- `tools/export_mechanism.py`：Cantera YAMLをSI機構データへ変換する前処理のみ。
  NSE/GPEの環境生成と同様にPythonを使用するが、計算時には不要。
- `reference/python/`：以前のPython実装を削除せず検証用に分離。
  `run_reference_reactor.py`は比較用で、通常の計算入口ではない。

時間積分は[Modern Fortran DVODE](https://github.com/jacobwilliams/dvode)の
固定commit `59d2eedd1410a4bd6c196c24d08f95a7c20565c0`。
可変刻み・次数1〜5のBDF（MF=22、数値差分の密Jacobian、陰的Newton）。
BLAS/LINPACKのFortran実装も含む。Python/SciPyを呼び出すラッパーではない。
初回CMake構成時に取得し、実行時にはネットワーク不要。
オフラインでは同commitのソースを用意して
`-DRF_DVODE_SOURCE_DIR=/absolute/path/to/dvode`を指定する。
取得先の`LICENSE.md`はBSD-3-Clause。バイナリ配布時もライセンス文書を添付する。
上流は開発中と明記しているため固定版と回帰テストを維持し、無検証で更新しない。

## ビルド・実行（FrameWorkルート）

Windows（PowerShell）：gfortran、CMake、NinjaをPATHへ登録する。
Pythonは機構変換だけに必要。初回はGitとネットワーク接続も必要。

```powershell
cmake -S SolverLibrary/ReactingFlow -B build/reactingflow -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/reactingflow --parallel 4
ctest --test-dir build/reactingflow --output-on-failure
python -m pip install Cantera==3.2.0 PyYAML
$mechanism = python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')"
python SolverLibrary/ReactingFlow/tools/export_mechanism.py "$mechanism" build/reactingflow/h2o2.rf
.\build\reactingflow\rf_reactor.exe build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/reactor_h2_air.in build/reactingflow/reactor_cv.csv
```

Linux（bash）：gfortran、CMake、Ninja、Gitをインストールし、Pythonは仮想環境を使う。

```bash
cmake -S SolverLibrary/ReactingFlow -B build/reactingflow -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/reactingflow --parallel 4
ctest --test-dir build/reactingflow --output-on-failure
python3 -m venv .venv
source .venv/bin/activate
python -m pip install Cantera==3.2.0 PyYAML
mechanism=$(python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')")
python SolverLibrary/ReactingFlow/tools/export_mechanism.py "$mechanism" build/reactingflow/h2o2.rf
./build/reactingflow/rf_reactor build/reactingflow/h2o2.rf SolverLibrary/ReactingFlow/examples/reactor_h2_air.in build/reactingflow/reactor_cv.csv
```

機構変換と計算結果は既存ファイルを上書きしない。再実行時は新しい出力名を指定する。
既存の別構成のビルドディレクトリは流用せず、専用ビルド先を使用する。
WindowsでgfortranのDLLが見つからない場合はコンパイラーのbinをPATHへ追加する。
既存ResearchRunsへの自動コピーやrun/postprocess統合はR7で行い、今回は切り替えない。

## 入力と方程式

第1引数は機構、第2引数はFortran namelistと組成、第3引数は出力CSV。
入力例は`examples/reactor_h2_air.in`。
namelist `/reactor/`の後に**機構の種順序に従う全質量分率**を1行で書く。
サンプルはCantera 3.2.0のh2o2.yaml専用。別機構へ同じ組成行を流用しない。
合計1、非負、種数一致が必要。一般Cantera YAMLをFortranへ直接読ませない。
変換済み機構は`RFMECH1`識別子、入力／canonicalのSHA256、種と反応の順序を保持する。

| namelist変数 | 既定値 | 意味 |
|---|---|---|
| mode | 'constant_volume' | 定圧は'constant_pressure' |
| temperature / pressure | 1000 / 101325 | 初期K / Pa |
| end_time | 0.001 | 終了時刻s |
| rtol | 1e-7 | 相対誤差（1e-12〜1e-2） |
| atol_species / atol_temperature | 1e-14 / 1e-6 | 質量分率 / Kの絶対誤差 |
| max_step / max_steps | 0.001 / 100000 | 最大刻みs / 最大ステップ数 |
| ignition_rise | 400 | 着火判定の初期温度からの上昇K |

閉鎖・断熱の理想混合気体、全量SI。
`dY_s/dt=omega_s/rho`、定容はrho固定で
`dT/dt=-sum(e_s*omega_s)/(rho*cv)`、定圧はp固定で
`dT/dt=-sum(h_s*omega_s)/(rho*cp)`。
`h_s=h_molar,s/M_s`、`e_s=(h_molar,s-RT)/M_s`。
生成エネルギーを含むため、診断用heat_releaseを別に加算しない。
定圧は閉鎖系の体積が変化する理想化で、流通反応器ではない。

初期最大質量分率の1種を従属変数として、`Y_dep=1-sum(Y_other)`で復元する。
積分状態はTと残りN-1種。これは出力の正規化ではなく、質量保存拘束の消去。
Newton試行のみ負値を0にして正規化した組成でRHSを評価する。
受理解をクリップせず、微小値も含む負のY、NASA範囲外、非有限値はエラー。
BDFは正値性保証法ではない。任意の条件で完走できる保証はなく、失敗時は
許容誤差とmax_stepを小さくして再検証する。試行温度の外挿もしない。

## 出力・保存検査

CSVは全受理時刻（不等間隔）の`time,temperature,pressure,density,Y_...`。
コメント行に機構ハッシュ、許容誤差、着火時刻、保存残差、積分統計を記録する。
着火時刻はT0+ignition_riseへの初回上向き交差をDVODE補間で求める。
未到達は`-1`。最大dT/dt時刻とは異なる。

- 質量残差：max|sum(Y)-1|、許容1e-12。
- 元素残差：元素モル量/mol kg⁻¹の初期値との差をmax(1,|初期値|)で割った最大値、許容1e-8。
- エネルギー残差：定容eまたは定圧hの初期値との差を
  max(1 J/kg,|初期保存量|,cp0*T0)で割った最大値、許容max(100*rtol,1e-7)。

計算中にCSVへ順次出力するため、異常終了では部分結果が残る。
**末尾の`# SUCCESS`と正常終了コードがあるものだけ完走結果**として扱う。
画面の`[OK] Fortran reactor completed`も完走時のみ。
閾値判定は収束検証の代替ではない。h2o2.yamlの1000 K切替には約0.14 J/kgの
エンタルピー不連続があり、相対エネルギー誤差は約9e-8で頭打ちになる。

## 比較テスト

FortranのみのCTestはPythonなしで動く。Cantera照合を含める場合は追加で以下を実行する。

Windows：

```powershell
python -m pip install -r SolverLibrary/ReactingFlow/requirements-reference.txt
$env:RF_FORTRAN_BUILD = (Resolve-Path build/reactingflow).Path
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -v
```

Linux：

```bash
python -m pip install -r SolverLibrary/ReactingFlow/requirements-reference.txt
export RF_FORTRAN_BUILD="$PWD/build/reactingflow"
python -m unittest discover -s SolverLibrary/ReactingFlow/tests -v
```

`RF_FORTRAN_BUILD`未指定ではFortran照合がskipされる。skipを検証済みと扱わない。
水素/空気の着火履歴、保存・許容誤差収束、NASA7/9、GRI30の速度、第三体、
Lindemann/Troe/SRI、PLOG、Chebyshevを照合する。
今回の実機検証はWindows/gfortran。Linuxコマンドは併記しているがLinux実行は未検証。
MPI/CUDA、反応流CFD、デトネーションは次段階であり、今回の完成範囲に含まない。
