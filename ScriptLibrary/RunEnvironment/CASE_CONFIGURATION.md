# ケース設定と拡張YAML

更新日: 2026-09-03

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

非反応計算では`chemistry.model`の既定値が`none`なので、空の化学反応設定ファイルは
生成しない。Stage 5だけが`config/chemistry.yaml`を生成する。Stage 5は格子輸送を
行わないため、`transport.yaml`は生成しない。

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
