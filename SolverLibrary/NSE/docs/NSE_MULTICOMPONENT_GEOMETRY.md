# Stage 9: 平面ノズルの一般座標・境界適合格子

## 実装範囲

対象は `nse_multicomponent` の静止・単一ブロック構造格子。
`cartesian`（既定）と `planar_nozzle` を選択できる。
既存の単成分 `nse`、Stage 0--8の設定・直交格子は変更しない。

平面ノズルは縮小拡大断面をz方向に押し出した形状で、状態は3次元で計算する。
任意CAD/外部格子ファイルの読込み、円形ノズルの軸対称方程式、軸上の特異点処理、
マルチブロック、移動格子/ALE、CUDAは未対応。
境界の `reflective` は自由滑り・断熱・非透過の静止壁である。
非滑り壁・指定壁温・壁面反応は未実装であり、実機の壁面境界層の予測には使用しない。

## 設定

`case.yaml`に以下を追加する。既存の4拡張参照はそのまま維持する。

```yaml
extensions:
  multicomponent: config/multicomponent.yaml
  thermodynamics: config/thermodynamics.yaml
  transport: config/transport.yaml
  chemistry: config/chemistry.yaml
  geometry: config/geometry.yaml
```

`config/geometry.yaml`:

```yaml
schema_version: 1
extension: geometry
config:
  type: planar_nozzle
  inlet_half_height: 0.20
  throat_half_height: 0.10
  exit_half_height: 0.25
  throat_x: 0.50
```

長さは他のソルバー入力と同じ単位系を用いる。半高さなので全流路高さは2倍。
格子のx,z座標は物理座標、y座標は計算座標ηであり、全領域で `y_min: -1`、
`y_max: 1` とする。y両面には通常 `reflective` を設定する。y周期は拒否する。
x周期を使う場合は入口・出口の半高さが一致する必要がある。
側壁があるため、壁の法線方向に速度を持つ一様流は定常解ではない。

現在のノズル生成器は、入口・スロート間、スロート・出口間をそれぞれ二次曲線で補間する。

- 左側: h(x) = ht + (hin-ht) ((x-xt)/(xmin-xt))²
- 右側: h(x) = ht + (hout-ht) ((x-xt)/(xmax-xt))²
- 写像: X=x, Y=η h(x), Z=z

実際の離散格子はxの各節点にhを与え、その間を直線で結ぶ。
半高さは正値、スロートはx領域内部でなければならない。
入力生成では現時点で `reactive_navier_stokes` モードを対象とする。
標準の3成分・一段反応機構は結合検証用であり、実在燃料の燃焼機構ではない。

## 数値処理

- 各共有面に一つの面積ベクトルを定義し、面流束の差を物理セル体積で除す。
  保存変数は従来どおり物理体積あたりの部分密度・Cartesian運動量・全エネルギー。
- Rusanov流束は物理面法線方向に投影する。対流精度は既存どおり一次。
- セル内の勾配を逆写像Jacobianで物理空間へ変換する。
  輸送面の勾配には隣接セル間の差に基づく非直交補正を施す。
- 粘性応力、熱伝導、化学種拡散を物理法線へ投影する。化学種拡散流束の総和はゼロ。
- 境界の状態と参照速度を面法線座標へ回転して既存の境界処理を適用し、物理座標へ戻す。
  無反射は従来のcharacteristic-relaxation近似であり厳密な変比熱NSCBCではない。
- 対流dtは各面の面積・法線速度・音速とセル体積から計算する。
  輸送dtには格子の縮小・非直交性を含む保守的な制限を加える。
  極端に細い/歪んだ格子を安全に自動修復する機能ではない。
- 化学反応は従来の局所Strang分割を使用する。
- 保存量/履歴は物理セル体積で重み付けする。CSVのx,y,zは物理空間のセル中心位置。
  CSVは点データであり、任意メッシュ読込みやVTK構造格子出力を追加したわけではない。

幾何を `mod_mc_geometry`、一般座標流束を `mod_mc_mapped_flux`、
物性依存の輸送を選択式providerへ分離した。
別の写像を追加する際は、面積ベクトルの閉包、一様場保持、体積正値性と
勾配変換の整合性を同時に検証する。

## MPI / OpenMP

既存profileをそのまま使用し、幾何のためのprofileは増やさない。

| 実行 | profile |
|---|---|
| 逐次 | cpu_serial_reactive_boundaries |
| OpenMP | cpu_openmp_reactive |
| MPI / MPI＋OpenMP | cpu_mpi_reactive_pencil |

`use_mpi`、`use_openmp`、`process_grid: [Py,Pz]` はStage 8と同じ。
xを保持し計算座標のy,zを分割する。MPI haloに対しても同一の写像を適用する。
化学反応は所有セルのみ、保存量は物理体積で集約する。
一般座標側は検証を優先した実装で、輸送勾配の一時配列を評価ごとに確保する。
従来直交格子の高速経路とworkspaceは維持する。

## 実行手順

FrameWorkルートから実行環境を生成する。

Windows PowerShell:

```powershell
python .\ScriptLibrary\RunEnvironment\prepare_environment.py .\ScriptLibrary\RunEnvironment\environment.nse_multicomponent.nozzle.yaml
```

Linux bash:

設計書の `select.target` を `linux_gnu_mpi`、
`select.destination` を `local_generated` に変更する。
必要なら `destination.root` と `destination.case_index` にLinux上の保存先を指定する。

```bash
python3 ./ScriptLibrary/RunEnvironment/prepare_environment.py ./ScriptLibrary/RunEnvironment/environment.nse_multicomponent.nozzle.yaml
```

生成先はコマンド末尾に表示される。そのディレクトリへ移動して次を実行する。

| 操作 | PowerShell | bash |
|---|---|---|
| 設定再生成 | `python .\tools\run_case.py --prepare` | `python3 ./tools/run_case.py --prepare` |
| ビルド | `python .\tools\run_case.py --build` | `python3 ./tools/run_case.py --build` |
| 計算 | `python .\tools\run_case.py --run` | `python3 ./tools/run_case.py --run` |

同梱例は圧力差を持つ反応流の短時間テストであり、定常ノズル性能計算ではない。
まず10ステップの動作を確認し、その後に境界参照状態、物性、反応機構、
格子解像度と時間条件を目的に合わせて設定する。

## 検証と残る確認

Windows/GNU Fortran＋Microsoft MPIで次を確認した。

- 離散的な面積ベクトルの閉包と正のセル体積。
- ノズル壁に接する一様な奥行き方向流の対流・輸送RHSがゼロ。
- 閉じた自由滑り断熱壁での質量・全エネルギー保存（反応Strang更新を含む）。
- MPI 1/2/3/4プロセスと逐次計算のRHS・反応流更新の一致。2スレッド併用。
- 既存の直交Euler/粘性/反応流テスト。
- YAML拡張から生成した環境のビルドと平面ノズル10ステップ実行。

Linux実機、長時間定常ノズル、格子収束、実在反応機構、実験との比較は未検証。
数値回帰の合格は、燃焼器性能やデトネーションの予測精度を保証しない。

