# 単成分NSEの不等間隔格子：実装状況と接続方針

## 本計算の接続（2026-10-07・最新）

限定した組合せで、格子生成→初期化→時間発展→座標付き出力を実行可能にした。
**全機能の不等間隔対応や、噴流の本番条件の検証が完了したという意味ではない。**
後述の開発履歴にある「共通停止ガード維持」は当時の状態であり、現在は以下の
許可リストに置き換わる。等間隔の既存設定・SLF出力は変更しない。

| 項目 | 不等間隔で現在選択できる内容 |
|---|---|
| 格子 | 直交・方向分離 `sinh`、各方向のstrength=0ならその方向は等間隔 |
| 初期条件 | `taylor_green`（物理セル中心の座標で評価）、`uniform_flow`（2026-10-08追加） |
| 対流 | `keep2`、`weno5z_roe`、`hybrid`（KEEP2＋WENO5Z_Roe） |
| 粘性・熱伝導 | `fv2`、または `none`。FV2にはRe>0、Pr>0が必要 |
| 時間積分 | 既存SSPRK3、固定dt／CFL自動dt |
| 境界 | 各面の周期・特性緩和無反射・鏡像・固定値。周期は方向ごとの両端ペア |
| 実行構成 | CPU/OpenMP、MPI、CUDA、MPI＋CUDA。従来の環境・実行窓口を使用 |
| 出力 | `vtr`。各rankのVTRと全rankを参照するPVTR。保存量5変数、物理セル辺座標、ghost除外 |

KEEP6は計量重みとFV物理体積の整合が未解決のため、本計算では引き続き拒否する。
CENTRAL6、HIT、読み込み乱流、restart、forcing、揺らぎもこの不等間隔モードでは
未対応として拒否する。これらは**等間隔では従来どおり**。別スキームへの自動変更はしない。
WENOという名前だけで全体5次精度を保証しない。粘性FV2は滑らかな格子で2次を意図する方式。

### 実行方法

#### 一様流・入口出口の確認（2026-10-08追加）

`flow.type: uniform_flow` では `flow.uniform_state` に密度、3方向速度、圧力を
明示する。単位は既存NSEと同じ無次元量であり、Mach数として再解釈しない。
指定がない場合や密度・圧力がsolverの下限以下の場合は生成時／読込み時に停止する。
等間隔格子でも使用できるが、既存の初期条件・既定値の意味は変更しない。

```yaml
flow:
  type: uniform_flow
  uniform_state:
    density: 1.2
    velocity: [0.13, -0.04, 0.02]
    pressure: 0.9
```

既存の `examples/nonuniform_tgv.case.yaml` のflow節を上記へ置き換えれば
全周期の一様流計算になる。固定値入口＋無反射出口の場合は、既存の
`boundary.faces.x_min` を `type: DIRICHLET`、`x_max` を `type: NON_REFLECTING`
とし、各 `reference_state` に上記と同じdensity/velocity/pressureを指定する。
y/zは両端周期とする。これは一定の入口状態であり、噴流の半径方向分布ではない。
境界の参照状態は初期値から自動で補わず、従来どおり明示指定する。

Fortranの直接入力では `initial_condition='uniform_flow'` と
`&nse uniform_state=1.2,0.13,-0.04,0.02,0.9 /` を指定する。
保存量は `[rho,rho*u,rho*v,rho*w,p/(gamma-1)+rho*|u|^2/2]` で初期化する。

検証スクリプトに `--uniform-flow` を追加すると、12³伸長格子で移動一様流の
保持（全保存量の各セル絶対差2e-11未満）を検査する。3対流×2粘性×
周期／固定値入口＋無反射出口×CPU・4 MPI・CUDA・4 MPI＋CUDAの48実行が成功。
MPI＋CUDAは1GPU共有の試験であり、実際の複数GPU通信の確認ではない。
既存のCUDA hybrid衝撃波回帰試験も維持する。今回CUDAカーネルの変更はない。
従来のTaylor–Greenについても同じ4構成・48実行を再確認した。
CPU CTest 35件、CUDAの対流／衝撃波／positivity試験7件、入力生成119件が成功。
今回の更新はFrameWork内のみで、既存ResearchRuns環境には自動反映しない。

#### 既存の環境生成・実行手順

完全な入力例は `examples/nonuniform_tgv.case.yaml`。環境生成時のMPI/CUDAの選択は
変更不要。生成済みの `cases/<case-id>/case.yaml` ではcase_id、出力先、solver設定を
自分の環境のまま維持し、例のflow/grid/numerics/output設定を反映する。
特に `grid.mapping` だけ変更せず、`viscous_scheme: fv2` と `output.format: vtr` も設定する。
既存の拡張config参照でforcing/FHが有効なら、その不等間隔ケースでは無効化が必要。

新しいソースで生成した実行環境のルートで、Windows PowerShell:

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --configure --build --run
```

Linux:

```bash
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --configure --build --run
```

古い実行環境はソースとtoolsのコピーを持つため、FrameWorkだけ更新しても変わらない。
従来の環境再生成手順で更新し、計算結果は新しい出力先へ保存する。
手動CMakeビルドの粘性バックエンドは `NSE_VISCOUS_SCHEME=central6` のままでよい。
このバックエンドに実行時選択のFV2も含む。`NSE_VISCOUS_SCHEME=none` のビルドではFV2は使えない。
MPIのy/z各局所ブロックはnghost以上のセル数を確保する。

ParaViewでは `output_nonuniform/field_*.pvtr` の連番を開く。
VTRは現在ASCII形式のためSLFより容量・書込時間が大きい。出力頻度に注意する。
不等間隔ではoutput_frequency>0なら、出力間隔の途中で終了した場合も最終状態を保存する。
output_frequency=0では時間発展後の出力を行わず、write_initialは独立に適用する。
`meta.json` のspacingはnullであり、最小セル幅を一様格子間隔と偽って記録しない。
座標はVTR内のセル辺から取得する。既存SLF用postprocess/FFT/import/restartは
VTRに未対応なので、この出力には使用しない。非等間隔データをそのまま通常FFTへ渡すことも不可。

### CFL・実用上の注意

自動dtの対流制限は `CFL / [max(|u_d|+c) * sum_d(1/hmin_d)]`。
全領域の方向別最小物理幅を使用し、CPU/CUDAで同じ定義とする。
粘性制限も最小幅と最小密度を使用する。固定dtは自動調整されないため、まず
低いCFLの自動dtで確認する。例の0.15は検証用初期値で、安定性の保証値ではない。
極端なstretchは最小幅を非常に小さくし、精度・安定性・計算量を悪化させる。

周期接続部や無反射境界で全体の高次精度を保証しない。
また、無反射は既存の特性緩和であって反射を厳密にゼロにはしない。
噴流入口の半径方向速度・温度分布の設定は別途未実装。
長時間噴流DNSへ移る前に、その初期条件・入口分布の追加と格子収束・反射率検証が必要。

### 自動検証

`tests/check_nonuniform_production.py` は12³、20ステップで、各3対流方式×2粘性方式×
全周期／x無反射・yz周期を実行する。後者は自動dt、前者は固定dt。
正の密度・圧力、5保存量の物理体積積分（周期）、実座標・MPI領域被覆・初期場、
CPU/MPI/CUDA間の場と時刻の一致を検査する。実行にはPython標準ライブラリのみ使用する。
例（FrameWorkルート、実行ファイルの場所は環境に合わせる）:

```powershell
python SolverLibrary/NSE/tests/check_nonuniform_production.py --cpu build/nse-positivity-cpu/bin/solver.exe --cuda build/nse-positivity-cuda/bin/nse_cuda.exe --work-dir build
```

```bash
python3 SolverLibrary/NSE/tests/check_nonuniform_production.py --cpu build/nse-cpu/bin/solver --cuda build/nse-cuda/bin/nse_cuda --work-dir build
```

`--mpi <exe> --mpiexec <launcher>` で4プロセスCPU比較、
`--mpi-cuda <exe>` で4プロセスCUDA比較を追加できる。
`--shared-gpu` は検証目的で全rankを1GPUへ割当てる。複数GPU間通信や性能の検証ではない。
許容誤差: CPU/GPU保存量の各セル値2e-10、時刻2e-12、周期体積積分2e-11×max(1,初期積分絶対値)。

今回の確認結果:

- 上記12条件×CPU/OpenMP・4 MPI・CUDA・4 MPI＋CUDAの計48実行が成功。
- MPI＋CUDAは1GPU共有・host-staged通信。複数GPU実機・CUDA-aware通信の再検証は未実施。
- CPUの12条件は `--steps 200` でも成功。長時間DNSの保証ではない。
- 既存CPU CTest 35件、case入力テスト118件が成功。
- 環境生成テストは51件中50件成功。残る1件は既存のforcingテンプレートに対し、
  廃止済み `target_dissipation` の文字列を期待するテストの不一致。
  今回そのforcing仕様やテストは変更していない。
- VTR/PVTRはXMLの読み戻し、セル数・座標・領域被覆・データを検証済み。
  ParaView GUIでの表示確認は未実施。

## この接続作業のコミット用コマンド

Gitの操作は未実施。この接続は以前の不等間隔演算・境界実装に依存するため、
それらが未コミットなら先に変更を確認して含めること。
下記は今回編集したファイルをまとめてステージするPowerShell用。
同じファイル内に別作業の変更がある場合はそれも含まれるため、diffを確認する。
今回のツール環境ではGitのworktree認識が失敗し、既存のステージ状態は確認できなかった。

```powershell
Set-Location 'C:\Users\Owner\Documents\Codex\FrameWork'
$nse = 'SolverLibrary/NSE'
$runenv = 'ScriptLibrary/RunEnvironment'
$files = @(
  "$nse/src/common/mod_common_config.f90",
  "$nse/src/grid/mod_grid_axis.f90",
  "$nse/src/io/mod_input_reader.f90",
  "$nse/src/io/mod_slf_output.f90",
  "$nse/src/time/mod_nse_time_integration.f90",
  "$nse/src/gpu/nse_cuda_bridge.cu",
  "$nse/src/main/main_nse.f90",
  "$nse/src/main/main_nse_cuda.f90",
  "$nse/src/main/main_nse_mpi_cuda.f90",
  "$nse/solver_manifest.yaml",
  "$nse/docs/NSE_NONUNIFORM_GRID.md",
  "$nse/docs/NSE_BUILD_AND_RUN.md",
  "$nse/examples/nonuniform_tgv.case.yaml",
  "$nse/tests/input_nonuniform_production.dat",
  "$nse/tests/check_nonuniform_production.py",
  "$runenv/case_input.py",
  "$runenv/case_templates/nse.yaml",
  "$runenv/tests/test_case_input.py"
)
git add -- $files
git diff --cached --stat
git diff --cached --check
# 内容を確認してから実行:
git commit -m "Enable validated nonuniform NSE runs with physical-coordinate VTK output"
```

Linuxの場合はリポジトリのルートへ移動し、同じファイル一覧を `git add --` に渡す。
commitコマンドは共通。pushは現在のブランチ・送信先を確認して別途行う。

## 開発履歴：周期・無反射の境界処理（本計算接続前）

CPU、CUDA、MPI+CUDAの格子生成に各方向の周期指定を接続した。
周期方向は反対側セルの幅を使い、座標を領域長だけ平行移動してghostへ延長する。
非周期方向は従来どおり端の幅列を鏡映する。周期は両端をセットで指定し、
片側だけ周期、反対側を無反射という設定は不可。各方向を独立に選択できる。

無反射は既存の `characteristic_relaxation` を継続使用する。
緩和係数を `alpha=1-exp(-strength*distance/L)` とし、非等間隔では
`distance=abs(x_ghost-x_interior_boundary_cell)` を使用する。
等間隔の `ghost_layer*dx` に一致する定義。
`length_scale: auto` のLはその方向の**全領域長**であり、最小幅×セル数や
MPI局所領域長ではない。CUDAにも全領域長と局所座標を転送する。
FV2用にも境界の角を含めた段階的halo交換を適用する。

検証済み:

- 4 MPIプロセス（y/z=2×2）＋OpenMP 2スレッドで、伸長格子の
  x無反射・y/z周期、および全方向無反射の境界単体試験成功。
- 単一GPUで、伸長格子のCPU/CUDA境界値比較成功（全ghost・角を含む）。
  比較許容誤差5e-12。等間隔の既存CPU/MPI・CUDA試験も成功。
- 格子生成への周期フラグ受渡し、周期ghost座標の単体試験を追加・維持。
- 続きの検証では、全6面・3層それぞれについてsinh写像の解析座標と
  特性緩和式から期待保存量を独立計算し、誤差2e-12以下で一致。
  `length_scale:auto`と明示長さ0.37を、等間隔・不等間隔の両方で確認した。
  GPU比較にはyのみ無反射、zのみ無反射、全周期、明示緩和長の組合せも追加。

**境界処理の接続であり、不等間隔本計算の共通停止ガードはまだ維持する。**
KEEP6の非周期高次閉包、周期sinh接続部の滑らかさ、粘性・積分重み・出力との
整合は別途残る。境界値比較は音波の反射率や噴流の長時間安定性を保証しない。
この特性緩和境界は厳密な無反射でも、完全な多次元粘性NSCBCでもない。

### 噴流向けの方向別設定例（境界設定の例）

```yaml
boundary:
  faces:
    x_min: {type: non_reflecting, reference_state: ambient}
    x_max: {type: non_reflecting, reference_state: ambient}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
  reference_states:
    ambient:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.7142857142857143
  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto
```

数値はgamma=1.4の周囲静止場の例で、実際の無次元化・周囲状態に合わせる。
この例は噴流を注入する設定ではない。入口で噴流を駆動する場合は
入口速度・密度・温度の分布指定が別途必要。横方向を周期にすると、
孤立した無限遠の噴流ではなく周期配列の計算になる。

## 2026-10-07：写像KEEP6空間演算

**KEEP6単独の非粘性空間演算をCPU/CUDAへ接続した。本計算の共通停止ガードは維持する。**
対象は滑らかな直交分離写像 x(i), y(j), z(k)。一般の曲がった格子ではない。
既存の等間隔KEEP6と同じ対称2点流束・計算座標上の6次係数を使用し、
分母に同じ差分から計算した正の計量を用いる。

```text
h_i = 3/4*(x[i+1]-x[i-1]) - 3/20*(x[i+2]-x[i-2])
      + 1/60*(x[i+3]-x[i-3])
R_i = -2/h_i * sum_s d_s*(F#(q_i,q_{i+s})-F#(q_{i-s},q_i))
d_s = [3/4,-3/20,1/60]
```

`h_i`は物理セル幅ではなく計算座標間隔を含む離散ヤコビアン。
3Dの積分重みは `h_x(i)*h_y(j)*h_z(k)`。
保存性・運動エネルギー対流収支は**この重み**に対して成立する。
物理セル幅積の既存 `vol` と同一ではないため、物理体積での厳密保存とは主張しない。
状態は写像格子点の点値として扱い、WENOのセル平均との整合は別途必要。
CPUではKEEP6単独の対流残差のみ計量を使う分岐を追加。
CUDAでは局所中心座標から同じ計量を初期化時に作成しGPUへ常駐させる。
等間隔時は既存演算のまま。非正・非有限計量はKEEP6で拒否する。

検証結果:

- 滑らかな周期写像 `x=s+0.2*sin(s)` の24/48/96点で、5保存量の微分の
  L2収束次数が **5.61276、5.89596**。
- 計量重み付き保存誤差 **5.55e-17**（同試験）。
- sinh伸長の固定物理区間 `0.2<=x<=0.8` で32→64点の質量微分誤差比
  **61.8499**。境界から固定セル数を除くだけでは比較区間が変わるため、
  この試験では同じ物理区間を比較する。
- 一様流保持、一定圧力の運動エネルギー対流収支、既存等間隔KEEP6回帰も検証。
- 関連単体・回帰10件成功。CPU/OpenMP、MPI、CUDA、MPI+CUDAでビルド成功。

残る制限:

- 非周期境界の6次閉包と、sinh周期接続部の写像の滑らかさは未解決。
  上記の内部・滑らかな周期写像の収束結果を境界までの6次精度と解釈しない。
- WENOとのハイブリッドは積分重みの統一が必要なため、非等間隔KEEP6を含む場合は
  CPU/CUDAで拒否する。KEEP2へ自動降格しない。
- 粘性・揺らぎ・forcing・SLF/統計・再スタートも計量重みと整合させる必要がある。
  共通停止ガードを手動で解除して本計算してはならない。
- CUDAはビルド確認まで。GPU数値一致・MPI分割数依存性は未検証。

写像計量と保存形を同時に扱う背景は
[curvilinear KEEP論文](https://www.sciencedirect.com/science/article/pii/S0021999121003776)
を参照。本実装は直交分離写像の限定実装であり、同論文の一般曲線座標方式を
全て実装したものではない。

## 2026-10-07：KEEP2の保存的流束

以下はKEEP2接続時点の履歴。KEEP6の最新状況は上の節を参照。

KEEP2の対称2点流束を非等間隔直交格子の面積・体積と組み合わせる。
CPUの既存空間演算は `(Aplus*Fplus-Aminus*Fminus)/V` を使用するため、
2点流束そのものに距離重みを入れない。CUDAのKEEP残差にはGPU常駐の局所セル幅を
接続した。WENOとのハイブリッドでもKEEP2側の流束は共用する。

**不等間隔KEEP6は未実装**。等間隔の6次係数をそのまま局所幅で割っても
非等間隔で6次精度にはならない。CPUの面流束入口とCUDAの格子アップロードで
この組合せを拒否する（ハイブリッドに含まれるKEEP6も対象）。
KEEP6をKEEP2へ自動降格しない。等間隔KEEP6は従来経路を維持する。
高次の保存的計量演算子とその境界閉包の実装は残る。

隣接セルの質量流束をm、算術平均速度をubarとすると、運動量流束の
対流部分は `m*ubar`。この関係を保持し、体積重み付き運動エネルギーの
対流収支を整合させる。圧縮性流では圧力仕事があるので、一般の流れで
運動エネルギー全体が一定という意味ではない。時間積分の厳密保存も保証しない。
理論背景は [Jamesonの有限体積KEP論文](https://citeseerx.ist.psu.edu/document?doi=1b9c3025db3dcc62ec592b8d66598069c5ba9866&repid=rep1&type=pdf) を参照。

`nse_test_nonuniform_keep`で周期データの体積重み付き5保存量収支、一定圧力下の
運動エネルギー対流収支、一様流保持、滑らかな伸長の内部2次収束を確認。
非等間隔KEEP6拒否、既存等間隔KEEP6精度試験も成功。急変格子・境界閉包・
GPU数値一致・本計算全体は未検証。本ソルバーの共通停止ガードは維持する。

## 2026-10-07：粘性・GPU常駐幾何・ハイブリッド検出

現在も本計算の停止ガードは維持する。以下は実装済みの構成要素であり、
不等間隔の全機能が使用可能になったという意味ではない。

- CPUに保存的な面流束型の粘性・熱伝導 `fv2` を追加済み。
  同じ面流束の差をセル幅で割る。従来の `central6` とは別方式で、
  CENTRAL6を黙って低次の方式へ変更しない。滑らかな伸長上の低次方式であり、
  急激な格子幅変化に対する2次収束を保証しない。
- CUDAにもFV2カーネルを追加済み。中心座標・セル幅・WENO面係数を各MPI領域に
  切り出して初期化時に一度転送し、GPU上に保持する。毎ステップの幾何転送は不要。
- CUDA WENOの再構築と流束差の局所セル幅補正を接続済み。
- CPU/CUDAのDucros検出に3点不等間隔微分を追加。
  圧力曲率を距離で重み付けし、物理座標で線形の圧力を格子伸長だけで検出しない。
  CUDAハイブリッドの流束差も局所セル幅を使用する。
  **KEEP2側は上記の最新節を参照。KEEP6側は未対応で、本計算はまだ使用不可。**

中心から左右の中心までの距離をhm,hpとすると、速度勾配は
`[hp*(u0-um)/hm + hm*(up-u0)/hp]/(hm+hp)`。
圧力検出は
`abs(hm*(pp-p0)-hp*(p0-pm))/(hm*pp+(hm+hp)*p0+hp*pm+epsilon*(hm+hp))`。
等間隔時は既存の処理経路を維持する。

今回のCPU/OpenMP単体試験は3件成功：非等間隔WENO係数・等間隔極限、
既存ハイブリッド回帰、不等間隔演算（線形圧力検出、FV2一様流保持・
線形速度場の粘性仕事、圧力ジャンプ検出、局所GPU幾何パックの添字）。
CPU/OpenMP・MPI・CUDA・MPI+CUDAのビルド成功。MPI+CUDAリンク時には既存の
`corrupt .drectve`警告が出るため、ビルド成功だけで実機動作確認済みとはしない。
長時間流体試験・GPU演算の数値比較は未実施。
FV2の格子収束、境界流束、安定時間刻み、MPI/GPU間の数値一致は未検証。

Windows/Linux共通の検証コマンド（CPU構成をビルド済みの場合）:

```text
cmake --build build/nse-positivity-cpu --target nse_test_nonuniform_operators nse_test_nonuniform_weno nse_test_hybrid_flux
ctest --test-dir build/nse-positivity-cpu -R "nse_(nonuniform_operators|nonuniform_weno|hybrid_flux)$" --output-on-failure
```

## 2026-10-06：CPU側WENO再構築への接続

`mod_reconstruction_nonuniform` を追加し、CPUのWENO5Z/Roe面流束に
不等間隔用の分岐を接続した。共通の特性変換・正値性制限・Roe流束は維持する。
等間隔時は既存関数を使用するため、今回の再構築へ自動的に置き換わらない。
本ソルバーの停止ガードは維持し、単独の不等間隔本計算はまだ許可しない。

5セルの境界からセル平均モーメント行列を作り、3個の二次候補多項式と
四次の最適多項式を計算する。面位置を0、左隣接セル幅を1に規格化した座標で
滑らかさは `(a1-a2)^2 + (13/3)*a2^2` とする。
最適線形重みは四次多項式と一致するよう生成し、非正の重みや不整合を拒否する。
右状態は座標と値を鏡映して同じ処理を使用する。
各面の係数は格子生成時に保存し、時間発展中の逆行列計算を避ける。

非線形重みは従来と同じZ型の `tau=abs(beta0-beta2)` を用い、
重み計算を対数表現にして大きな比のオーバーフローを避ける。
任意伸長における非線形収束次数・衝撃波挙動は**未検証**であり、
四次多項式の線形再現性だけで不等間隔WENOの5次収束を保証しない。
KEEP6・ハイブリッド全体は未対応。KEEP2・粘性とCUDA構成要素は上記の最新節を参照。

CPU/OpenMP構成でビルド成功。等間隔極限、一定値保持、セル平均多項式再現の
検証コード `nse_test_nonuniform_weno` を追加し、2026-10-07に実行成功。

2026-10-06。対象は単成分NSE。CPU、MPI、OpenMP、CUDA、MPI+CUDAの既存選択機能を
維持して直交不等間隔格子に対応させる作業。**複数のCPU/GPU構成要素を接続。本計算は未対応。**
等間隔格子の生成・数値処理・環境選択は維持する。

## 入力・プレビュー（追加実装）

case.yamlで次を指定すると、環境の種類とは独立してinput.datへ変換される。
既存のnx/ny/nz、領域範囲などと同じ `grid` 内に記述する。

```yaml
grid:
  # nx/ny/nz、領域範囲などの既存項目も必要
  mapping:
    type: sinh
    strength: [0.0, 2.0, 2.0]
```

これは現時点では**幾何プレビュー専用**。本ソルバーに渡すと、
CPU/MPI/CUDA/MPI+CUDA共通の入力段階で「演算未接続」と明示して停止する。
solver側からこの停止を解除する入力スイッチはない。
`mapping` を省略すれば以前と同じinput.datを生成する。
`type: uniform` のstrengthは全て0に限る。
未知キー、3要素以外、非有限値、負値、20超、真偽値を拒否する。
GPEへsinhを転用する設定は拒否する。

伸長用に全域の軸情報から局所範囲のセル中心・面積・体積を生成する経路を追加した。
MPI通信を行う処理ではなく、各領域が同じ全域座標から切り出す構成。
非周期ゴースト座標延長では物理端の幅を鏡映する。周期境界との接続は最新節を参照。
sim%dx/dy/dzは伸長時には全域最小幅となるため、平均間隔とは区別する必要がある。

プレビューは各軸の左右境界・中心・幅をCSVへ出す（3D配列は作らない）。
ソルバー、MPIランチャー、GPU計算を起動せず、1プロセスで実行する。

Windows PowerShell、既存CPUビルドを使う例：

```powershell
cmake --build build/nse-positivity-cpu --target nse_grid_preview
.\build\nse-positivity-cpu\bin\nse_grid_preview.exe SolverLibrary/NSE/tests/input_grid_preview.dat build/grid_preview.csv
```

Linux、同じ名前でCPUビルドを構成済みの場合：

```bash
cmake --build build/nse-positivity-cpu --target nse_grid_preview
./build/nse-positivity-cpu/bin/nse_grid_preview SolverLibrary/NSE/tests/input_grid_preview.dat build/grid_preview.csv
```

生成環境ではサンプルの代わりに自分のinput.datを指定できる。
出力先の親フォルダが必要。既存CSVは上書きしない。
今回の確認はCPU構成のビルドのみ。実行例、Python回帰、MPI/CUDA実行は未検証。

## 今回追加したもの

`src/grid/mod_grid_axis.f90`（Fortran、MPI/CUDA依存なし）：

- `build_axis(edges,ng,periodic,axis)`：明示した昇順のセル境界座標から1方向の幾何を生成。
- `build_sinh_axis(n,ng,lower,upper,strength,periodic,axis)`：中心を細かくする対称伸長。
  強度0で等間隔。強度の範囲0〜20。極端な座標で境界が重なる場合は拒否。
- 物理境界0〜n、セル1〜nと、ng層のゴースト座標・セル幅・セル中心。
  周期は反対側の幅、非周期は端のセル幅列を鏡映して座標を延長する。
  これは幾何の延長のみで、物理境界条件を設定する処理ではない。
- 全域の最小セル幅。
- 7点の点値に対する1階・2階微分係数。CPU/GPUに共通の係数表を渡すためのデータ。
  係数は `(-3:3,global_cell)` の順に保持する。この点微分表のGPU利用は未接続。
  別のWENO面係数と中心座標・セル幅のGPU転送・常駐は接続済み。

伸長関数は、s=i/nとして
`x=lower+(upper-lower)*(1+sinh(strength*(2*s-1))/sinh(strength))/2`。
噴流中心が領域中央にない場合や片側だけの伸長は、明示座標で表す想定。
APIにはn>=ng>=3の制限がある。

微分係数は6次までの多項式点値の微分を再現するもの。
**有限体積のセル平均再構築係数でもKEEPの保存的流束係数でもない。**
任意の不等間隔格子で全演算の6次精度を保証するものではなく、
現行KEEP6/WENO5Z/CENTRAL6への単純置換は禁止する。
座標だけ変更して既存の等間隔用係数を流用する実装にはしない。

## 残作業（完了条件）

1. 共通入力：case.yaml→input.dat、テンプレートコメント、既定uniform、値検査は追加済み。
   実行環境生成を通した回帰確認は残る。
2. 格子と入出力：既存 `mod_grid_fvm` への接続、MPI局所範囲への切出し、
   SLF座標とメタデータ、初期化と境界の距離依存処理。伸長用の局所幾何生成までは追加済み。
3. 保存的演算：KEEP2/KEEP6、WENO5Z、ハイブリッド、粘性・熱流束・CFLを
   不等間隔用に整合させる。形式的な点微分の精度と保存性を別々に確認する。
4. GPU：同じ座標・係数をGPUへ常駐させ、CUDA/MPI+CUDAの演算を接続。
   既存CUDA-aware MPIとホスト経由の選択も保持する。
5. 拡張機能：揺らぎのセル体積・散逸離散化、FFTベースのHIT初期化・フォーシング、
   読込乱流の座標対応、後処理との整合を確認する。
   通常FFTは非等間隔の物理座標をそのまま扱えないため、格子写像や補間の定義が必要。
6. 検証：等間隔回帰、保存量、一様流保持、格子収束、境界、MPI分割数、CPU/GPU比較。

いずれかの環境を廃止したり、不等間隔指定を黙って等間隔へ戻したりしない。
上記が未完了である間は、不等間隔計算が可能になったとは案内しない。
入力のsinh指定はプレビューにのみ使用でき、未接続の計算へ入ることはない。

## 確認状況

- ソースを共通CMakeと環境生成用マニフェストに登録。
- 既存CMakeキャッシュでもソース一覧へ追加されるよう対応。
- CPU/OpenMP構成で係数生成の単体テスト実行ファイルまでビルド成功。
- テストコードは等間隔極限、中心伸長、周期ゴースト、多項式微分、重複座標拒否を含む。
- 最新の単体試験結果は冒頭を参照。MPI/CUDAの本計算による検証は未実施。
