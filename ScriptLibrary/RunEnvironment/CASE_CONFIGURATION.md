# ケース設定と拡張YAML

更新日: 2026-09-04

## 1. 方針

通常の格子、初期条件、時間積分、数値解法、並列数、出力設定は、従来どおり
`case.yaml`へ記述する。多成分、熱力学、輸送、将来の化学反応など、追加モデルだけを
`config/*.yaml`へ分離する。

旧来の一体型`case.yaml`もそのまま読み込める。既存ケースを移行しなくてもよい。
新しい分割形式を使う場合だけ、`case.yaml`へ`schema_version: 2`と`extensions`を記述する。

## 2. 生成される構成

Stage 4の多成分・粘性ケースは次の構成になる。

```text
cases/case0001/
├── case.yaml
├── config/
│   ├── multicomponent.yaml
│   ├── thermodynamics.yaml
│   └── transport.yaml
├── resolved_case.yaml
├── input.dat
└── notes.md
```

`case.yaml`は利用する拡張だけを明示的に参照する。

```yaml
schema_version: 2

physics:
  model: nse_multicomponent

extensions:
  multicomponent: config/multicomponent.yaml
  thermodynamics: config/thermodynamics.yaml
  transport: config/transport.yaml

grid:
  nx: 64
  ny: 4
  nz: 4
```

拡張ファイルは、種類と設定本体を次の形式で宣言する。

```yaml
schema_version: 1
extension: transport
config:
  model: mixture_averaged
  reference_dynamic_viscosity: 1.8e-5
  prandtl_number: 0.72
```

現在登録されている拡張と、合成後の配置先は次のとおり。

| `extensions`のキー | 合成後の配置先 |
|---|---|
| `multicomponent` | `physics.multicomponent` |
| `thermodynamics` | `thermodynamics` |
| `transport` | `transport` |
| `chemistry` | `chemistry` |

ファイルは自動検索しない。`case.yaml`の`extensions`に書かれたファイルだけを読む。
このため、別の計算用設定が同じディレクトリにあっても誤って混入しない。

## 3. Stageごとに生成されるファイル

| 計算 | 生成する拡張YAML |
|---|---|
| Stage 0 基盤 | `multicomponent.yaml` |
| Stage 1 パッシブスカラー | `multicomponent.yaml` |
| Stage 2 非粘性多成分 | `multicomponent.yaml`, `thermodynamics.yaml` |
| Stage 3 温度依存熱力学 | `multicomponent.yaml`, `thermodynamics.yaml` |
| Stage 4 粘性・拡散・熱伝導 | `multicomponent.yaml`, `thermodynamics.yaml`, `transport.yaml` |
| Stage 5 0次元有限反応速度化学 | `multicomponent.yaml`, `thermodynamics.yaml`, `chemistry.yaml` |
| Stage 6 反応性多成分流 | `multicomponent.yaml`, `thermodynamics.yaml`, `transport.yaml`, `chemistry.yaml` |
| Stage 7 反応衝撃波管・面別境界 | `multicomponent.yaml`, `thermodynamics.yaml`, `transport.yaml`, `chemistry.yaml` |

非反応計算では`chemistry.model`の既定値が`none`なので、空の化学反応設定ファイルは
生成しない。Stage 5からStage 7が`config/chemistry.yaml`を生成する。Stage 5は格子輸送を
行わないため`transport.yaml`を生成せず、Stage 6とStage 7は流体輸送と反応を結合するため両方を生成する。

Stage 5の反応式はspecies名をキーにして記述する。`orders`を省略すると反応物の
量論係数を反応次数として使う。入力生成時にspecies名、係数の正値性、分子量を
含む反応式の質量保存を検査する。

```yaml
schema_version: 1
extension: chemistry
config:
  model: one_step_arrhenius
  reaction:
    reactants: {fuel: 1.0, oxidizer: 1.0}
    products: {product: 2.0}
    pre_exponential_factor: 1000.0
    temperature_exponent: 0.0
    activation_temperature: 2000.0
```

0次元反応器の初期状態、時間刻み、出力条件は通常の実行条件なので`case.yaml`に
残す。質量分率もspecies名で指定する。

```yaml
flow:
  type: homogeneous_reactor
  homogeneous_reactor:
    density: 1.0
    temperature: 1200.0
    mass_fractions: {fuel: 0.5, oxidizer: 0.5, product: 0.0}

time:
  dt: 0.0
  maximum_dt: 1.0e-3
  chemistry_cfl: 0.1
  nsteps: 100
```

Stage 6では流体と化学の刻み制約、結合方式を`case.yaml`へ記述する。
`dt: 0.0`では対流・拡散・化学制約の最小値を使い、固定`dt`では化学半stepだけを
`maximum_chemistry_substeps`の範囲内で自動subcycleする。

```yaml
time:
  cfl: 0.2
  diffusion_cfl: 0.4
  chemistry_cfl: 0.1
  maximum_chemistry_substeps: 10000
  dt: 0.0
  nsteps: 20

numerics:
  convective_scheme: rusanov1
  boundary_condition: periodic
  time_integration: ssprk3
  coupling_scheme: strang
  chemistry_time_integration: ssprk3_subcycled
```

Stage 7では6物理面を`boundary.faces`で個別に指定する。`dirichlet`と
`non_reflecting`は、密度、3方向速度、圧力、全species質量分率を持つ
`reference_state`を参照する。周期境界は同じ方向の両面を対で指定する。

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

boundary:
  faces:
    x_min: {type: dirichlet, reference_state: driver}
    x_max: {type: non_reflecting, reference_state: far_field}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
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

output:
  write_final: true
  filename: multicomponent_reactive_shock_tube_final.csv
  write_snapshots: true
  write_history: true
  output_every: 5
  snapshot_prefix: multicomponent_reactive_shock_tube
  history_filename: multicomponent_reactive_shock_tube_history.csv
```

完全なひな型は`case_templates/nse_multicomponent_reactive_shock_tube.yaml`にある。
Stage 6 profileへStage 7専用の初期条件または面別境界を指定した場合は、対応profileへの
変更を促すエラーにする。

## 4. 編集と入力再生成

利用者が編集するのは`case.yaml`と、そこから参照される`config/*.yaml`である。
`resolved_case.yaml`と`input.dat`は生成物なので直接編集しない。

Windows（PowerShell）:

```powershell
python .\tools\run_case.py --prepare
```

Linux（bash）:

```bash
python3 ./tools/run_case.py --prepare
```

`run_case.py`は`case.yaml`だけでなく、参照中の拡張YAMLの更新時刻も確認する。
いずれかが`input.dat`より新しければ、自動的に設定を再合成して入力を再生成する。

`case_input.py`を直接使った場合も、入力生成に成功すると同じケースディレクトリへ
`resolved_case.yaml`を出力する。このファイルは全拡張を展開した単独で読める設定で、
実行条件の確認と再現性記録に使用できる。

`created_at`などの日時・日付は、`resolved_case.yaml`ではISO形式の文字列として
保存する。YAML読込み時に日時型へ変換される環境でも生成でき、タイムゾーンと
小数秒を保持する。元の`case.yaml`や拡張YAMLを変更する必要はない。

## 5. エラーにする条件

曖昧な設定や可搬性を損なう参照は、暗黙の優先順位を付けず入力生成前に拒否する。

- 同じ項目を`case.yaml`本体と拡張YAMLの両方に定義した
- 未登録の拡張名を指定した
- `case.yaml`の`extensions`と拡張ファイル内の`extension`が一致しない
- 拡張ファイルの`schema_version`が1ではない
- 絶対パス、`..`、シンボリックリンク、junctionを参照した
- 同じ拡張ファイルを複数の拡張名から参照した
- 拡張ファイルに`schema_version`、`extension`、`config`以外の最上位キーがある

拡張パスは常に`case.yaml`が置かれたディレクトリを基準に解決する。したがって、
Windowsで生成したケースをLinuxへ移しても、相対ディレクトリ構成を保てば変更は不要である。

## 6. 旧形式との互換性

次のような従来の単一ファイル形式は引き続き有効である。

```yaml
physics:
  model: nse_multicomponent
  multicomponent:
    mode: inviscid_euler
    species: [species_a, species_b]

thermodynamics:
  model: calorically_perfect
  gamma: 1.4
```

旧形式には`extensions`を書かない。分割形式と旧形式の定義を混在させると重複エラーに
なるため、移行するときは対象セクション全体を拡張ファイルへ移す。

## 7. 開発者が拡張種類を追加する場合

新しい種類を追加するときは、次の順で明示的に登録する。

1. `case_configuration.py`の`EXTENSION_TARGETS`へ名前と合成先を追加する。
2. `case_templates/extensions/`へ拡張テンプレートを追加する。
3. 必要な計算だけ、`environment_options.yaml`の`extension_templates`へ登録する。
4. `environment.schema.json`の拡張名一覧を更新する。
5. 合成先、重複拒否、ケース外参照拒否、入力再生成を自動テストへ追加する。

この登録手順により、新機能を追加しても既存の単成分NSEや不要な計算設定へ項目を
増やさずに済む。
