# 多成分反応流の面別境界と時系列出力（Stage 7）

更新日: 2026-09-04

## 対応範囲

Stage 7はStage 6の熱的完全混合気体、混合平均輸送、一段不可逆Arrhenius反応、
Strang分割を維持したまま、次を追加するCPU逐次検証段階である。

- 6物理面ごとの`periodic`、`reflective`、`dirichlet`、`non_reflecting`
- 密度、3方向速度、圧力、全species質量分率を含む参照状態
- x方向の二状態初期条件`reactive_shock_tube_x`
- step 0、指定間隔、最終stepの場スナップショット
- species質量、全質量、運動量、全エネルギー、最小物理量の時系列

使用するsolver profileは`cpu_serial_reactive_boundaries`である。Stage 6の
`cpu_serial_reactive`は全方向周期の既存仕様として残り、同じ入力を継続利用できる。

## 境界条件

境界は`x_min`、`x_max`、`y_min`、`y_max`、`z_min`、`z_max`を個別に指定する。
片側だけを周期境界にはできず、同じ方向の両面を必ず`periodic`にする。

| `type` | 処理 |
|---|---|
| `periodic` | 反対側セルを使用する |
| `reflective` | 法線運動量だけを反転する自由滑り・断熱鏡像境界 |
| `dirichlet` | 指定した参照状態を境界外状態として使用する |
| `non_reflecting` | 外向き特性を内部から、流入特性を参照状態へ緩和して構成する |

`non_reflecting`は内部状態と参照状態の平均比熱比を固定した
characteristic-relaxation近似である。`relaxation_strength`は非負値、
`length_scale`は正値または`auto`を指定する。`auto`では境界法線方向の領域長を使う。
これはStage 7の検証用境界であり、変比熱の厳密なNSCBCや反応性流入・流出モデルではない。

```yaml
boundary:
  faces:
    x_min: {type: dirichlet, reference_state: driver}
    x_max: {type: non_reflecting, reference_state: far_field}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: reflective}
    z_max: {type: reflective}
  reference_states:
    driver:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 356334.11220656743
      mass_fractions: {fuel: 0.45, oxidizer: 0.45, product: 0.10}
    far_field:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 267250.5841549256
      mass_fractions: {fuel: 0.49, oxidizer: 0.49, product: 0.02}
  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto
```

`dirichlet`と`non_reflecting`には`reference_state`が必須である。参照状態の
species名は多成分設定と完全一致し、質量分率は非負かつ総和1でなければならない。
対流流束と粘性・拡散・熱伝導の両方が同じ面別境界状態を使用する。

## 反応衝撃波管初期条件

`reactive_shock_tube_x`は`interface_location`の左右へ、密度、速度、圧力、
全species質量分率から構成した二つの熱力学的状態を配置する。

```yaml
flow:
  type: reactive_shock_tube
  reactive_shock_tube:
    initial_condition: reactive_shock_tube_x
    interface_location: 0.35
    left:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 356334.11220656743
      mass_fractions: {fuel: 0.45, oxidizer: 0.45, product: 0.10}
    right:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 267250.5841549256
      mass_fractions: {fuel: 0.49, oxidizer: 0.49, product: 0.02}
```

## 出力

```yaml
output:
  write_final: true
  filename: multicomponent_reactive_shock_tube_final.csv
  write_snapshots: true
  write_history: true
  output_every: 5
  snapshot_prefix: multicomponent_reactive_shock_tube
  history_filename: multicomponent_reactive_shock_tube_history.csv
```

`write_snapshots: true`ではstep 0、`output_every`の倍数、最終stepへ
`<snapshot_prefix>_########.csv`を出力する。最終stepが指定間隔と一致しても重複しない。
`write_final`の最終場ファイルは従来どおり別に出力する。

履歴CSVにはstep、時刻、`dt`、species別質量、全質量、3方向運動量、全エネルギー、
最小species部分密度、最小混合密度、最小圧力、最小温度を記録する。非周期境界では
領域内保存量は境界流束によって変化するため、終了時の値は保存誤差ではなく
初期値に対する最大相対領域総量変化として表示する。

## ビルドと実行

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\BuildSolver\build_model.py `
  .\ScriptLibrary\BuildSolver\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_reactive_boundaries `
  --test

python .\ScriptLibrary\RunEnvironment\prepare_environment.py `
  .\ScriptLibrary\RunEnvironment\environment.nse_multicomponent.reactive_boundaries.yaml
```

Linux（bash）:

```bash
python3 ./ScriptLibrary/BuildSolver/build_model.py \
  ./ScriptLibrary/BuildSolver/build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_reactive_boundaries \
  --test

python3 ./ScriptLibrary/RunEnvironment/prepare_environment.py \
  ./ScriptLibrary/RunEnvironment/environment.nse_multicomponent.reactive_boundaries.yaml
```

生成されたケースでは`tools/run_case.py --prepare`、`--build`、`--run`の順で実行する。

## 現在の制約

### 2026-09-07: 逐次RHSの高速化

対流・輸送RHSは、各セルの密度・速度・質量分率・温度・圧力・音速を
RHS評価ごとに一度生成して共有する。状態を更新した後は必ず再計算するため、
SSPRK各段や化学反応後の状態に古いキャッシュを使用しない。
種密度・混合密度・圧力の検査をこの生成処理に統合し、非有限値も検査する。

Euler流束は方向別の格子列を走査し、直前の面流束を保持して内部面と周期接続面を
一度だけ評価する。全格子サイズの面流束配列は追加しない。
輸送項は既存RHSへ直接加算し、独立した`transport_rhs`を廃止した。
作業配列は実行中再利用し、格子サイズ・化学種数が変わると再確保する。

追加のYAML設定は不要である。既存のRHS・時間積分APIも、末尾の任意引数
`workspace`を省略して呼べる。省略時は呼出し内で一時workspaceを使用する。
実行ドライバはworkspaceを保持して全ステップで再利用する。
温度逆算の二分法、数値スキーム、出力形式、並列化設定は従来どおりである。

`nse_multicomponent_rhs_workspace`で旧六面評価との一致、全方向の境界条件、
単一セル軸、状態更新後の再利用、格子サイズ変更、輸送項の直接加算を検証する。
同テストの時間表示は小規模Euler RHSの比較であり、計算全体の速度倍率ではない。

### 対応範囲

- CPU逐次、三次元直交等間隔格子だけに対応する。
- 対流流束は一次Rusanov、時間積分はSSPRK3、結合はStrang分割である。
- 反応機構は一段不可逆Arrhenius反応だけである。
- `reflective`はno-slip壁ではない。
- MPI、OpenMP、CUDA、および一般座標はStage 7の対象外である。
