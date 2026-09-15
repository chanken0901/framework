# R4a/b/c：1次元流体・詳細反応・輸送の結合

R4aのEuler基盤、R4bの定係数輸送、R4cのMUSCL再構築を実装した。
**Fortran、CPU逐次実行、1次元・一様直交格子の基準計算**。
R4全体の完成ではない。温度・組成依存の分子輸送、2D/3D、ノズル、MPI/OpenMP/CUDA、
CJ/ZNDデトネーション検証、旧反応流CFDとの全体同値検証は残っている。
現在のコードを実用的なデトネーションソルバーと見なさない。

## できること

- 左右のT/P/u/Yを指定した衝撃波管・接触面・一様場の初期化。
- `chemistry=.false.`による非反応多成分Euler計算。
- `chemistry=.true.`による、R2で対応する詳細機構を使う反応Euler計算。
- `transport_model='constant'`による粘性・熱伝導・種拡散。反応の有無とは独立。
- `reconstruction='muscl'`による対流の空間2次再構築。既定は従来の1次精度。
- 両端周期、左右別々の鏡像壁／ゼロ勾配流出。
- 初期状態、指定ステップ間隔、最終時刻のCSV出力と保存量監視。

## 数値仕様

保存変数は`[rho*Y_1,...,rho*Y_N,rho*u,rho*E]`。
全密度は種密度の和。`E=e(T,Y)+u²/2`で生成エネルギーを含む。
一般的な定比熱比式ではなくNASA熱力学からT・p・音速を復元する。

流束は旧`mod_mc_euler_flux`と同じRusanov式を独立状態配列／EOSへ移した。
空間は既定で区分一定の**1次精度**。MUSCLを選択すると滑らかな領域で2次精度。
流体時間積分はSSPRK3。
化学は各セル定容DVODE BDF、Strang分割（化学dt/2→流体dt→化学dt/2）。
この分割の形式的時間精度は2次だが、衝撃波や硬い反応では次数低下があり得る。
時間・空間・化学許容誤差への収束を個別に確認する。

化学過程ではrho・rho*u・rho*Eを変えず、種を更新した後に保存エネルギーからTを復元する。
ODE温度との差が`max(1e-3 K,100*chemistry_rtol*T)`を超えたら停止。
診断用反応熱の別加算はしない。各セル元素残差も検査する。
圧縮性多成分の接触面では保存型混合による圧力誤差が生じ得る。
MUSCL以外の高次再構築、低散逸流束、圧力平衡保持処理はまだない。

輸送無効時の刻みは`min(max_dt,終了までの時間,CFL*dx/max(|u|+c))`。
第1化学半ステップ後と各流体段の音速を再確認し、CFL超過時は元の状態から
dtを半分にしてやり直す（最大30回）。これで反応後の音速増大も考慮する。
輸送有効時は、後述の拡散制限も毎段で併用する。
負の種密度、不正な内部エネルギー、NASA範囲外はクリップせず停止する。
**任意の強い衝撃波に対する正値性保証や汎用リカバリーは未実装**。

境界の`outflow`はゼロ勾配外挿であり、無反射境界ではない。
`reflecting`は法線運動量だけ符号反転する滑り鏡像壁。
`periodic`は両端同時指定のみ。Dirichletと無反射は未移植。

## R4b：粘性・熱伝導・化学種拡散

定係数の基準輸送モデルを独立した`mod_rf_transport`で実装した。
Fortran APIの`rf_transport`に係数をまとめ、既存の流体APIでは末尾の任意引数として受け取る。
以前の引数・入力はそのまま使用できる。係数を全て0にした結果はEuler版と一致する。
将来の温度・組成依存モデルはこの輸送モジュールに追加する。

| namelist設定 | 単位・既定値 | 内容 |
|---|---|---|
| transport_model | 'none' | 'constant'で明示指定した係数を有効化 |
| viscosity | Pa s、0 | 動粘度ではなく粘性係数mu |
| bulk_viscosity | Pa s、0 | 体積粘性係数zeta（0ならStokes仮説） |
| thermal_conductivity | W/(m K)、0 | 熱伝導率kappa |
| mass_diffusivity | m²/s、0 | 全化学種に共通の質量分率勾配ベース拡散係数D |

全係数は有限・非負。noneのまま非ゼロ係数を指定した場合は、黙って無視せずエラー。
constantでも個別係数を0にできるため、熱伝導だけなどの切り分けが可能。
機構YAMLの輸送パラメーターから自動算出しない。例の係数は検証用であり実在気体の推奨値ではない。
**化学種ごとのD、温度依存粘性、混合平均／多成分輸送、Soret、Dufour、圧力拡散は未実装**。
これらが必要な燃焼速度・火炎構造を定量評価できる完成段階とは見なさない。

全流束は`F=F_Euler+F_diff`として同じセル面で保存的に差分する。

```text
J_s^raw = -rho_face * D * (Y_s,R - Y_s,L)/dx
J_s = J_s^raw - Y_s,face * sum(J^raw)
tau = (4*mu/3 + zeta) * (u_R - u_L)/dx
F_diff(species) = J_s
F_diff(momentum) = -tau
F_diff(energy) = -u_face*tau - kappa*(T_R-T_L)/dx + sum(h_s(T_face)*J_s)
```

面のrho/u/T/Yは左右算術平均、勾配は隣接セルの差分。輸送の空間離散化は2次。
最後の1種で丸め誤差のみを閉じてsum(J)=0とする（濃度のクリップではない）。
h_sには生成エンタルピーを含み、種拡散に伴うエネルギー輸送を省略しない。
対流の既定値は空間1次のRusanov。滑らかな解の空間2次計算には下記のMUSCLを選ぶ。

輸送は流体と同じSSPRK3の各段に入る。化学とのStrang分割は変更しない。
刻みの保守的な目安は次の通りで、max_dtと終了時刻でも制限する。

```text
nu_bound = (4*mu/3+zeta)/min(rho)
          + kappa/min(rho*cv) + D*max(rho)/min(rho)
dt <= CFL / (max(|u|+c)/dx + 2*nu_bound/dx²)
```

圧縮性のエネルギー式に合わせ熱拡散の制限にはcpではなくcvを用いる。
これは非線形系全体の正値性保証ではない。格子を細かくすると拡散制限はdx²に比例して厳しくなる。
周期境界では面流束を共有する。鏡像壁はu_wall=0、断熱・種の不透過を課し、
壁の粘性応力は運動量の境界収支に入る。壁面エネルギー・種流束は0。
outflowでは全原始変数の法線勾配0として、追加の輸送流束は0。

輸送付きサンプルは`examples/flow1d_h2_transport.in`。
下記の準備済み機構を使い、Windowsでは：

```powershell
.\build\reactingflow\rf_flow1d.exe build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_transport.in build/reactingflow/h2_transport.csv
```

Linuxでは：

```bash
./build/reactingflow/rf_flow1d build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_transport.in build/reactingflow/h2_transport.csv
```

出力コメントにモデルと全輸送係数を記録する。境界補正した保存検査には輸送流束も含む。
参考：[Canteraの種拡散流束・補正速度の説明](https://www.cantera.org/stable/reference/onedim/governing-equations.html)。
本実装の共通定係数Fickモデルを、Canteraの混合平均輸送モデルと同一とは扱わない。

## R4c：MUSCL再構築

`&flow1d`に次を追加する。省略すれば従来の結果を維持する。

```fortran
 reconstruction='muscl', ! 'first_order' (default) or 'muscl'
```

Riemann流束はどちらもRusanov。`muscl`はT・p・u・全質量分率Yを原始変数として、
MC（monotonized central）リミター付きの左右面状態を作る。
最大質量分率の種を従属傾斜とし、種の傾斜和を0にする。
その後、全種共通の縮小係数で左右面のYを隣接セルの最小・最大の範囲に制限する。
セル平均の濃度をクリップ・再正規化する方式ではない。
面上の保存変数はNASA熱力学とEOSで生成し、単一の共有面流束でセルを更新する。
温度・圧力もMC制限するが、多成分系全体のTVD性・更新後の正値性を保証するものではない。
衝撃波・極値付近ではリミターが作用し、局所的に次数が下がる。

対流CFLの波速はセル中心だけでなく左右再構築面も含めて評価し、SSPRK3各段で再確認する。
鏡像境界では内側の再構築面状態の運動量を反転、周期境界では周期隣接セルを参照する。
輸送項は引き続きセル中心間の2次差分であり、対流の再構築とは独立。
反応の有無、輸送の有無、既存の3種類の境界と組み合わせられる。
WENO/KEEPの移植やCJ/ZND検証が済んだという意味ではない。

Fortran APIでは末尾の任意引数`reconstruction`で指定する。
`flow_timestep`を単独でMUSCLに使う場合は`left_bc`と`right_bc`も渡す（省略時はoutflow）。
`advance_flow`と実行プログラムは指定された境界を刻み評価にも渡す。

反応＋輸送＋MUSCLの例は`examples/flow1d_h2_muscl.in`。
下記の準備済み機構を使い、Windows PowerShellでは：

```powershell
.\build\reactingflow\rf_flow1d.exe build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_muscl.in build/reactingflow/h2_muscl.csv
```

Linux bashでは：

```bash
./build/reactingflow/rf_flow1d build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_muscl.in build/reactingflow/h2_muscl.csv
```

CSVコメントの`reconstruction`に採用方式を記録する。未知の方式名は拒否する。

## ビルド・実行

機構変換やビルドの準備は[R1〜R3手順](REACTORS.md)と共通。
以前生成した実行環境・case.yamlは変更しない。まだ`run_case.py`には統合していない。

Windows PowerShell、FrameWorkルート：

```powershell
cmake -S SolverLibrary/ReactingFlow -B build/reactingflow -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/reactingflow --parallel 4
ctest --test-dir build/reactingflow --output-on-failure
python -m pip install Cantera==3.2.0 PyYAML
$mechanism = python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')"
python SolverLibrary/ReactingFlow/tools/export_mechanism.py "$mechanism" build/reactingflow/h2o2_flow.rf
.\build\reactingflow\rf_flow1d.exe build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_n2_shock.in build/reactingflow/n2_shock.csv
.\build\reactingflow\rf_flow1d.exe build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_reactive.in build/reactingflow/h2_reactive.csv
```

Linux bash、FrameWorkルート（Python依存は仮想環境内へ導入）：

```bash
cmake -S SolverLibrary/ReactingFlow -B build/reactingflow -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/reactingflow --parallel 4
ctest --test-dir build/reactingflow --output-on-failure
python3 -m venv .venv
source .venv/bin/activate
python -m pip install Cantera==3.2.0 PyYAML
mechanism=$(python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')")
python SolverLibrary/ReactingFlow/tools/export_mechanism.py "$mechanism" build/reactingflow/h2o2_flow.rf
./build/reactingflow/rf_flow1d build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_n2_shock.in build/reactingflow/n2_shock.csv
./build/reactingflow/rf_flow1d build/reactingflow/h2o2_flow.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_reactive.in build/reactingflow/h2_reactive.csv
```

全て出力の上書きは禁止。再実行は新しい名前を指定する。
Windows/gfortranで検証。Linuxコマンドは併記したが、この環境でのLinux実行は未検証。

## 入力・出力

Fortran namelist `&flow1d ... /`の後に左右の全質量分率を各1行、機構の種順序で書く。
サンプルはCantera 3.2.0 h2o2.yamlの順序専用。

| 設定 | 意味・制約 |
|---|---|
| nx, length, interface_x | セル数>=2、領域[0,length] m、左右状態の切替位置m |
| left/right_temperature, pressure, velocity | 各側のK、Pa、m/s |
| left_bc, right_bc | periodic / reflecting / outflow |
| chemistry | .true.で反応、既定.false. |
| reconstruction | first_order（既定）/ muscl（MC制限付き空間2次） |
| end_time, max_dt, cfl | 終了s、最大刻みs、0<CFL<=0.5 |
| max_steps, write_every | 受理ステップ上限、書き出し間隔（正の整数） |
| chemistry_rtol | 化学相対誤差、既定1e-9 |
| chemistry_atol_species, chemistry_atol_temperature | 化学絶対誤差、既定1e-16 / 1e-8 K |
| chemistry_max_steps | 1セルの1化学半ステップに許すBDFステップ数、既定100000 |

CSV列は`step,time,x,density,velocity,temperature,pressure,Y_...`。
SLFではないため既存のNSE用FFT・ParaView変換へ直接渡さない。
最終時刻はwrite_everyの倍数でなくても必ず書く。
末尾`# SUCCESS`と終了コード0がある結果だけが完走。途中停止は部分CSVが残る。

質量・運動量・全エネルギー・元素について、境界流束の時間積分を差し引いて保存残差を検査する。
SSPRK3の各段の境界寄与は1/6、1/6、2/3の重みで集計する。
残差はmax(1,|初期積分量|)で割る。元素はmolの積分量（単位断面積当たり）。
いずれかが1e-7以上なら停止。鏡像壁からの運動量変化は境界力として含める。

## 検証・残作業

実施：Fortran流束単体、周期一様場、静止鏡像壁、無反応化学、
非反応接触面の種保存、定比熱理想気体Sodの厳密解に対する格子収束、
一様水素/空気とCantera定容反応器の比較、非一様反応場の刻み細分化と保存検査。
R4bでは粘性応力・粘性仕事・Fourier熱流束・種エンタルピー流束、断熱壁、
周期熱／種拡散の正弦波減衰と空間2次収束、拡散による刻み制限、
反応＋輸送の保存、輸送係数0でのEuler版との一致、不正係数の拒否を検証した。
正弦波の解析解比較は輸送演算子を単独で検証し、Rusanovの数値拡散を混ぜていない。
R4cでは滑らかな周期組成波・密度波移流のセル平均厳密解に対する32/64/128セルの2次収束、
Sod厳密解に対する同一格子での1次精度との比較、面の種範囲・総和、周期／壁の保存、
反応＋輸送＋MUSCLの保存、既定／明示first_orderの完全一致、不正方式の拒否を検証した。
Sod参考：[Clawpack Euler Riemann問題](https://www.clawpack.org/riemann_book/html/Euler_approximate.html)。
同梱窒素ショック例はNASAの温度依存比熱であり、定比熱Sod厳密解とは別物。

R4残作業：旧CFDとの全体同値検証、分子輸送モデルと境界の拡充、より高次・頑健性検証、
ノズル／一般座標、1D反応波のCJ速度・ZND構造・格子／分割誤差の検証。
この確認が済むまではR5（並列化）完了へ進めない。
