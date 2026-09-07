# NSE多成分・反応流拡張仕様

更新日: 2026-09-04

## 1. 目的と実装状態

現行の単成分NSEを基準実装として維持しながら、混合性多成分気体と有限反応速度
反応流を段階的に追加する。本拡張は現行NSEと別manifest、別profileで構成し、
既存ソースへ多成分条件分岐を追加しない。

| モデル | profile | 状態 |
|---|---|---|
| `nse` | 既存profile | 現行単成分NSE |
| `nse_multicomponent` | `cpu_serial_foundation` | Stage 0基盤 |
| `nse_multicomponent` | `cpu_serial_passive_scalar` | Stage 1パッシブスカラー |
| `nse_multicomponent` | `cpu_serial_inviscid` | Stage 2非反応・非粘性多成分流 |
| `nse_multicomponent` | `cpu_serial_thermally_perfect` | Stage 3温度・組成依存熱力学 |
| `nse_multicomponent` | `cpu_serial_viscous` | Stage 4拡散・粘性・熱伝導 |
| `nse_multicomponent` | `cpu_serial_reactor` | Stage 5 0次元有限反応速度化学 |
| `nse_multicomponent` | `cpu_serial_reactive` | Stage 6反応性多成分Navier--Stokes |
| `nse_multicomponent` | `cpu_serial_reactive_boundaries` | Stage 7反応流境界・初期条件・時系列出力 |

多成分profileの実行ファイル名はすべて`nse_multicomponent`である。profileごとに
main programとコンパイル対象を切り替えるため、各Stageの実行内容は混在しない。

## 2. Stage 0: 拡張契約

Stage 0は流体計算を進めず、状態レイアウト、species一覧、provider選択、
1成分極限で5保存変数になることを検証する。

```yaml
physics:
  multicomponent:
    mode: foundation
```

利用するproviderは次のとおりである。

- 熱力学: `calorically_perfect`
- 輸送: `none`
- 化学反応: `none`

## 3. Stage 1: パッシブスカラー移流

Stage 1は一定の搬送速度`U_c`により各部分密度を保存形式で移流する。

```text
∂rho_s/∂t + div(rho_s U_c) = 0
```

運動量と全エネルギーは固定し、speciesから流体へのフィードバックを行わない。
全方向周期境界、一次風上有限体積法、SSPRK3を使用する。物理拡散はなく、
`transport.model`は`none`である。一次風上法に由来する数値拡散は存在する。

## 4. Stage 2: 非反応・非粘性多成分流

Stage 2は全成分部分密度、運動量、全エネルギーを連成して進める。

```text
∂rho_s/∂t       + div(rho_s u)                 = 0
∂(rho u)/∂t     + div(rho u⊗u + pI)            = 0
∂(rho E)/∂t     + div((rho E + p)u)            = 0
rho             = sum_s rho_s
```

化学反応源項、粘性応力、熱伝導、species拡散は含まない。したがってStage 2は
Navier-Stokes方程式ではなく、多成分Euler方程式である。

### 4.1 熱力学

Stage 2の`calorically_perfect` providerは全speciesで共通の比熱比`gamma`を使用する。

```text
p = (gamma - 1) [rho E - |rho u|²/(2 rho)]
c = sqrt(gamma p / rho)
```

組成依存比熱、温度依存比熱、分子量からの混合気体定数はStage 3以降で追加する。

### 4.2 数値流束

セル界面には一次精度Rusanov流束を使用する。

```text
F* = 0.5(F_L + F_R) - 0.5 a_max (Q_R - Q_L)
a_max = max(|u_n,L| + c_L, |u_n,R| + c_R)
```

各species流束の総和は全質量流束と一致する。全方向で同じ流束関数を使用し、
周期面でも同じ保存的な面流束を共有する。

### 4.3 時間積分とCFL

時間積分はSSPRK3である。`dt: 0.0`の場合は毎step、現在の状態から次式で算出する。

```text
dt = CFL / max_cells[
  (|u|+c)/dx + (|v|+c)/dy + (|w|+c)/dz
]
```

固定`dt`の場合も多次元Courant数が1以下であることを検査する。各SSPRK段の前後で
部分密度、混合密度、圧力の正値性を検査し、破れた場合は直ちに停止する。

### 4.4 初期条件と出力

初期条件はx方向の多成分Sod問題`multispecies_sod_x`である。左右それぞれに
密度、速度、圧力、全speciesの質量分率を指定する。質量分率は非負かつ総和1を要求する。

最終CSVには座標、全species部分密度、混合密度、速度、圧力、全エネルギーを出力する。

## 5. Stage 3: 温度・組成依存熱力学

Stage 3はStage 2と同じ非反応・非粘性Euler方程式を解き、熱力学providerを
`thermally_perfect`へ交換する。各speciesについて分子量とNASA-7係数の低温・高温域を
`config/thermodynamics.yaml`へ指定する。係数配列は`physics.multicomponent.species`の名前をキーにするため、
化学種の並び順を変更しても別speciesの物性を誤って割り当てない。

```text
R_s       = R_u / W_s
R_mix     = sum_s Y_s R_s
p         = rho R_mix T
cp_mix(T) = sum_s Y_s cp_s(T)
gamma_mix = cp_mix / (cp_mix - R_mix)
c         = sqrt(gamma_mix p / rho)
```

NASA-7の`a1`～`a5`から`cp_s(T)`、`a1`～`a6`から基準生成エンタルピーを含む
内部エネルギーを評価する。`a7`も将来の反応平衡・可逆反応に備えて保持する。
保存変数の内部エネルギーから温度を求める際は、設定温度範囲内の単調な二分法を使う。
範囲外エネルギー、`cv <= 0`、中間温度で不連続な係数は計算開始前に拒否する。

Stage 3のSI単位規約は次のとおり。

- 分子量: kg/kmol
- 普遍気体定数: J/(kmol K)
- 温度: K
- 密度、速度、圧力、全エネルギーも整合するSI単位系

最終CSVにはStage 2の変数に加えて温度`T`を出力する。Stage 2の共通`gamma`モデルと
既存出力形式は`cpu_serial_inviscid`で引き続き利用できる。
NASA-7の`a6`が定める生成エンタルピーの基準によっては`rhoE`が負になり得るが、
温度・圧力が正で設定範囲内なら異常ではない。Stage 3の正値性判定も`rhoE`の符号ではなく、
部分密度、混合密度、温度、圧力に対して行う。

## 6. Stage 4: 混合平均輸送

Stage 4はStage 3の熱力学に、化学種拡散、Newton粘性およびFourier熱伝導を加えた
非反応多成分Navier--Stokes方程式を解く。対流流束はStage 2/3と同じRusanov法、
輸送流束はセル中心勾配を面へ補間する二次精度中心差分で評価する。

```text
J_s^0 = -rho D_s grad(Y_s)
J_s   = J_s^0 - Y_s sum_r J_r^0
sum_s J_s = 0

tau = mu [grad(u) + grad(u)^T - (2/3) div(u) I]
k   = mu cp_mix / Pr
```

保存式へ加える輸送項は次のとおりである。

```text
species:  -div(J_s)
momentum:  div(tau)
energy:    div(tau.u + k grad(T) - sum_s h_s J_s)
```

`mixture_averaged` providerは現段階では一定の基準粘性係数、一定の各成分拡散係数、
一定Prandtl数を受け取る。熱伝導率は局所温度・組成から求めた`cp_mix`を通して変化する。
拡散補正により各面の全成分拡散質量流束を厳密にゼロとし、全質量を保存する。
エネルギー流束にはNASA-7から得た各成分エンタルピーを含める。

自動時間刻みは対流CFL制約と次の陽的拡散制約の小さい方を使う。

```text
alpha_max = max_s(D_s, mu/rho, k/[rho cv_mix])
dt_diff   = diffusion_cfl /
            {2 alpha_max (1/dx^2 + 1/dy^2 + 1/dz^2)}
```

検証用初期条件`periodic_species_wave_x`は、一定密度・温度・速度の周期場で2成分の
質量分率に逆符号の正弦波を加える。これは拡散による振幅減衰、成分流束和ゼロ、
周期領域の保存性を分離して確認するための初期条件である。

## 7. Stage 5: 0次元有限反応速度化学

Stage 5は流体輸送と化学反応を分離して検証する、断熱・定容の0次元均質反応器である。
保存状態はStage 0--4と同じで、化学反応源だけをSSPRK3で時間積分する。

```text
d(rho_s)/dt = omega_s
d(rho u)/dt = 0
d(rho E)/dt = 0
```

最初の反応providerは一段不可逆Arrhenius反応である。`C_s=rho_s/W_s`を
kmol/m3単位のモル濃度とすると、進行速度と質量生成速度は次式になる。

```text
q = A T^beta exp(-Ta/T) product_s(C_s^alpha_s)
omega_s = W_s (nu_product,s - nu_reactant,s) q
```

`A`の単位は総反応次数に依存する。`Ta`は活性化温度[K]であり、活性化エネルギーを
指定する場合は`Ta=Ea/Ru`へ変換する。指数計算はlog空間で評価し、非常に小さい
反応速度は明示的にゼロへ丸める。これにより正常な低反応速度で
`IEEE_UNDERFLOW_FLAG`や`IEEE_DENORMAL`を発生させない。

反応式はspecies名をキーにして入力する。入力生成時とFortran初期化時の両方で、
speciesの存在、反応物と生成物の分離、係数の非負性、分子量を含む質量保存を検査する。
`orders`を省略すると反応物の量論係数を反応次数として使う。

自動時間刻みは、消費されるspeciesが負にならないためのdepletion時間を使う。

```text
dt = min(maximum_dt, chemistry_cfl * min_s[rho_s/(-omega_s)])
```

`dt: 0.0`でこの自動刻みを選ぶ。正の固定`dt`も指定できるが、各stepで上式の制約を
満たさない場合は停止する。各SSPRK段で部分密度、混合密度、温度、圧力を検査する。

反応源の全エネルギー成分はゼロである。NASA-7の`a6`を含む生成エネルギーが
保存状態に含まれるため、発熱反応では組成変化と全エネルギー保存から温度が上昇する。
履歴CSVには初期状態（step 0）、指定間隔、必ず最終stepを出力し、step、時刻、dt、
密度、温度、圧力、全species質量分率を記録する。

## 8. Stage 6: 流体・反応結合

Stage 6はStage 4の多成分Navier--Stokes輸送とStage 5の有限反応速度化学を、
Strang分割で結合する。

```text
Q^(1) = C(dt/2) Q^n
Q^(2) = F(dt)   Q^(1)
Q^(n+1) = C(dt/2) Q^(2)
```

`C`は各セル独立の化学反応、`F`は対流、species拡散、Newton粘性、Fourier熱伝導を
含む流体更新である。流体更新にはStage 4のSSPRK3をそのまま使用する。化学更新も
Stage 5と同じSSPRK3を共通integratorから呼び出し、各半stepをspecies depletion時間で
必要な回数だけsubcycleする。固定`dt`でも化学反応だけを自動分割できる。

`dt: 0.0`の場合、流体の対流・拡散制約と化学反応のdepletion制約の最小値を採用する。

```text
dt = min(dt_convective, dt_diffusive, dt_chemistry)
```

化学半step後にも流体CFLを再検査し、部分密度、混合密度、圧力、温度を各演算子の
更新中に検証する。周期領域では総質量、三方向運動量、全エネルギーを保存する。
最終CSVには座標、全species部分密度、混合密度、速度、圧力、温度、全エネルギーを
出力する。

初期Stage 6は、周期直交格子上の反応species波を検証問題とする。一様反応場では
Stage 5と同じ時間発展になり、反応速度ゼロではStage 4と同じ流体更新になることを
回帰試験で確認する。

## 9. Stage 7: 反応流境界・初期条件・出力

Stage 7はStage 6の流体・反応結合を維持し、6物理面を独立指定できる境界APIを追加する。
各面の選択肢は`periodic`、`reflective`、`dirichlet`、`non_reflecting`である。
`reflective`は法線運動量だけを反転する自由滑り・断熱鏡像境界、`dirichlet`は指定した
参照状態を境界外状態にする。`non_reflecting`は内部と参照状態の平均比熱比を固定し、
外向き特性を内部から、流入特性を参照状態へ緩和して構成するcharacteristic-relaxation
近似である。参照状態には密度、3方向速度、圧力、全species質量分率を指定する。

境界条件は対流流束と輸送流束の両方へ適用する。周期境界は同じ座標方向の両面を対で
指定しなければならない。初期条件`reactive_shock_tube_x`はx方向の界面位置を境に、
左右それぞれの密度、速度、圧力、組成から二状態反応流を生成する。

場のCSVはstep 0、指定間隔、必ず最終stepへ出力できる。積分履歴にはspecies別質量、
全質量、3方向運動量、全エネルギー、最小species部分密度、最小混合密度、最小圧力、
最小温度を記録する。非周期境界では領域内保存量が境界流束により変化するため、
終了時には保存誤差ではなく初期値に対する最大相対領域総量変化を表示する。

設定、式、出力名の詳細は`docs/NSE_MULTICOMPONENT_REACTIVE_BOUNDARIES.md`を参照する。

## 10. 保存変数

全Stageで状態レイアウトを共有する。

```text
Q = [rho_1, ..., rho_Ns, rho*u, rho*v, rho*w, rho*E]
rho = sum_s rho_s
nvar = Ns + 4
```

`Ns=1`では5保存変数になる。変数番号は`mod_mc_state_layout`だけが決定する。

## 11. Provider境界

熱力学、輸送、化学反応は同一APIを持つ代替Fortran moduleとして実装し、manifestが
ビルド時に各1個を選択する。

```text
mod_mc_thermodynamics_provider
mod_mc_transport_provider
mod_mc_chemistry_provider
```

Stage 2では熱力学providerが混合密度、圧力、音速を提供する。Stage 3では同じAPIに
温度、混合気体定数、混合比熱比、primitive状態からの全エネルギー生成を追加した。
Stage 4では輸送providerを`mixture_averaged`へ交換し、Stage 0--3では`none`を維持する。
Stage 5では輸送providerを`none`に戻し、化学反応providerを
`one_step_arrhenius`へ交換する。Stage 0--4の化学反応providerは引き続き`none`である。
Stage 6とStage 7では`thermally_perfect`、`mixture_averaged`、`one_step_arrhenius`を同時に選択する。
入力名とコンパイル済みproviderが異なる場合は開始前に停止する。

## 12. 設計書

新しく生成する多成分ケースは、共通条件を`case.yaml`、拡張固有条件を
`config/*.yaml`へ分ける。詳細な規約は
`ScriptLibrary/RunEnvironment/CASE_CONFIGURATION.md`を参照する。

Stage 2の`case.yaml`は利用する拡張を明示的に参照する。

```yaml
schema_version: 2

physics:
  model: nse_multicomponent

extensions:
  multicomponent: config/multicomponent.yaml
  thermodynamics: config/thermodynamics.yaml

flow:
  type: multispecies_sod
  multispecies_sod:
    initial_condition: multispecies_sod_x
    interface_location: 0.5
    left:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 1.0
      mass_fractions: [0.8, 0.2]
    right:
      density: 0.125
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.1
      mass_fractions: [0.2, 0.8]

numerics:
  convective_scheme: rusanov1
  boundary_condition: periodic
  time_integration: ssprk3
```

`config/multicomponent.yaml`には拡張の識別情報と多成分設定を書く。

```yaml
schema_version: 1
extension: multicomponent
config:
  mode: inviscid_euler
  species: [species_a, species_b]
```

Stage 2の`config/thermodynamics.yaml`は次のとおりである。

```yaml
schema_version: 1
extension: thermodynamics
config:
  model: calorically_perfect
  gamma: 1.4
```

Stage 3では`case_templates/nse_multicomponent_thermally_perfect.yaml`をひな型にする。
物性は`config/thermodynamics.yaml`でspecies名をキーにして指定する。

```yaml
schema_version: 1
extension: thermodynamics
config:
  model: thermally_perfect
  universal_gas_constant: 8314.46261815324
  temperature_min: 200.0
  temperature_max: 6000.0
  species_data:
    N2:
      molecular_weight: 28.0134
      temperature_midpoint: 1000.0
      nasa7_low: [3.53100528, -1.23660987e-4, -5.02999433e-7, 2.43530612e-9, -1.40881235e-12, -1046.97628, 2.96747468]
      nasa7_high: [2.95257626, 1.39690040e-3, -4.92631603e-7, 7.86010367e-11, -4.60755321e-15, -923.948645, 5.87188762]
```

Stage 4では`case_templates/nse_multicomponent_viscous.yaml`をひな型にする。
輸送物性は`config/transport.yaml`へ分離する。species名をキーにするため、species一覧との
過不足を入力生成時に拒否する。

```yaml
schema_version: 1
extension: transport
config:
  model: mixture_averaged
  reference_dynamic_viscosity: 1.8e-5
  prandtl_number: 0.72
  species_data:
    N2: {diffusivity: 2.0e-5}
    O2: {diffusivity: 2.0e-5}
```

格子、初期条件、時間条件はStage 4でも`case.yaml`に残す。

```yaml
flow:
  type: periodic_species_wave
  periodic_species_wave:
    initial_condition: periodic_species_wave_x
    density: 1.0
    temperature: 300.0
    velocity: [0.0, 0.0, 0.0]
    mean_mass_fractions: [0.5, 0.5]
    positive_species: N2
    negative_species: O2
    amplitude: 0.1
    wavenumber: 1

time:
  cfl: 0.2
  diffusion_cfl: 0.4
  dt: 0.0
```

Stage 0--4は非反応なので`chemistry`ファイルを生成せず、既定値`none`を使う。
Stage 5では`case_templates/nse_multicomponent_reactor.yaml`をひな型にし、
`config/chemistry.yaml`を追加する。

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

初期状態と時間・出力条件は`case.yaml`へ記述する。

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

output:
  write_history: true
  output_every: 1
  filename: homogeneous_reactor.csv
```

Stage 6は`case_templates/nse_multicomponent_reactive.yaml`をひな型にし、
`multicomponent`、`thermodynamics`、`transport`、`chemistry`の4拡張を参照する。
結合方式と化学積分方式は`case.yaml`で明示する。

```yaml
extensions:
  multicomponent: config/multicomponent.yaml
  thermodynamics: config/thermodynamics.yaml
  transport: config/transport.yaml
  chemistry: config/chemistry.yaml

time:
  cfl: 0.2
  diffusion_cfl: 0.4
  chemistry_cfl: 0.1
  maximum_chemistry_substeps: 10000
  dt: 0.0

numerics:
  convective_scheme: rusanov1
  boundary_condition: periodic
  time_integration: ssprk3
  coupling_scheme: strang
  chemistry_time_integration: ssprk3_subcycled
```

Stage 7は`case_templates/nse_multicomponent_reactive_shock_tube.yaml`をひな型にする。
4拡張はStage 6と同じであり、初期条件、面別境界、参照状態、時系列出力を
`case.yaml`に追加する。

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

旧来の一体型`case.yaml`も引き続き読み込めるが、同じ設定を本体と拡張ファイルへ
重複定義することはできない。入力生成時には全設定を展開した`resolved_case.yaml`も
保存する。

## 13. ビルドと実行環境

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\BuildSolver\build_model.py `
  .\ScriptLibrary\BuildSolver\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_viscous `
  --test
```

Linux（bash）:

```bash
python3 ./ScriptLibrary/BuildSolver/build_model.py \
  ./ScriptLibrary/BuildSolver/build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_viscous \
  --test
```

実行環境設計書は
`ScriptLibrary/RunEnvironment/environment.nse_multicomponent.viscous.yaml`
を使用する。

Stage 5はprofileを`cpu_serial_reactor`へ、実行環境設計書を
`ScriptLibrary/RunEnvironment/environment.nse_multicomponent.reactor.yaml`へ変更する。

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\BuildSolver\build_model.py `
  .\ScriptLibrary\BuildSolver\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_reactor `
  --test
```

Linux（bash）:

```bash
python3 ./ScriptLibrary/BuildSolver/build_model.py \
  ./ScriptLibrary/BuildSolver/build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_reactor \
  --test
```

Stage 6はprofile `cpu_serial_reactive`と実行環境設計書
`ScriptLibrary/RunEnvironment/environment.nse_multicomponent.reactive.yaml`を使用する。

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\BuildSolver\build_model.py `
  .\ScriptLibrary\BuildSolver\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_reactive `
  --test
```

Linux（bash）:

```bash
python3 ./ScriptLibrary/BuildSolver/build_model.py \
  ./ScriptLibrary/BuildSolver/build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_reactive \
  --test
```

Stage 7はprofile `cpu_serial_reactive_boundaries`と実行環境設計書
`ScriptLibrary/RunEnvironment/environment.nse_multicomponent.reactive_boundaries.yaml`を使用する。

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\BuildSolver\build_model.py `
  .\ScriptLibrary\BuildSolver\build.yaml `
  --model nse_multicomponent `
  --profile cpu_serial_reactive_boundaries `
  --test
```

Linux（bash）:

```bash
python3 ./ScriptLibrary/BuildSolver/build_model.py \
  ./ScriptLibrary/BuildSolver/build.yaml \
  --model nse_multicomponent \
  --profile cpu_serial_reactive_boundaries \
  --test
```

## 14. 制約

### 14.1 Stage 4

- CPU逐次実行のみ
- 三次元直交等間隔格子
- 全方向周期境界のみ
- 一次Rusanov流束のみ
- 理想混合気体のみ（実在気体効果なし）
- NASA-7係数の設定温度範囲内のみ
- 化学反応なし
- 粘性係数、Prandtl数、species拡散係数は入力値で一定
- Soret効果、Dufour効果、圧力拡散、実在気体輸送なし
- x方向の周期species波初期条件のみ
- CSVは最終時刻だけ出力

### 14.2 Stage 5

- CPU逐次の0次元均質反応器のみ
- 断熱・定容のみ
- 一段不可逆反応を1本だけ使用可能
- 圧力依存反応、第三体、Falloff、可逆反応、平衡定数は未実装
- 陽的SSPRK3のみで、stiff chemistry用陰解法は未実装
- 空間輸送、粘性、拡散との結合は未実装
- 標準テンプレートの3species物性・反応係数は数値検証用であり、実在反応機構ではない

### 14.3 Stage 6

- CPU逐次実行のみ
- 三次元直交等間隔格子、全方向周期境界のみ
- 対流は一次Rusanov流束、流体と化学の結合はStrang分割のみ
- 一段不可逆反応を1本だけ使用可能
- 化学反応は陽的SSPRK3 subcyclingで、stiff chemistry用陰解法は未実装
- 初期条件はx方向の周期species波のみ
- 標準テンプレートの物性・反応係数は結合検証用であり、実在反応機構ではない
- MPI、OpenMP、CUDAは未対応

### 14.4 Stage 7

- CPU逐次、三次元直交等間隔格子のみ
- 対流は一次Rusanov流束、時間積分はSSPRK3、流体・化学結合はStrang分割のみ
- 一段不可逆反応を1本だけ使用可能
- `non_reflecting`は固定平均比熱比のcharacteristic-relaxation近似であり、厳密な変比熱NSCBCではない
- `reflective`は自由滑り・断熱鏡像境界であり、no-slip壁ではない
- MPI、OpenMP、CUDA、一般座標は未対応
- 標準テンプレートの物性・反応係数は数値結合検証用であり、実在反応機構ではない

未実装の選択肢は黙って別方式として扱わず、入力生成時または計算開始前に拒否する。

## 15. 後続Stage

8. MPI/OpenMP最適化
9. 一般座標
10. CUDA

## 16. 必須回帰条件

- 現行`nse`のmanifest、実行ファイル名、既定profileを維持する。
- Stage 0/1 profileを独立してビルド・実行できる。
- `Ns=1`で5保存変数になる。
- 同一左右状態のRusanov流束が物理流束と一致する。
- species流束の総和が全質量流束と一致する。
- 一様周期場の離散右辺がゼロになる。
- 周期計算で全species質量、運動量、全エネルギーを保存する。
- 部分密度、混合密度、圧力の正値性を維持する。
- primitive状態から生成した全エネルギーから温度と圧力を復元できる。
- 混合比熱比が温度と組成に応じて変化する。
- NASA-7の温度範囲外やspecies物性の不一致を開始前に拒否する。
- 反応追加後も成分生成速度の総和をゼロにする。
- 0次元反応器で反応物が減少し、生成物が増加する。
- 発熱反応で全エネルギーを保存しながら温度が上昇する。
- 化学反応の自動時間刻みがspecies depletion制約を満たす。
- 成分拡散流束の総和を各セルでゼロにする。
- 一様場の輸送右辺をゼロにする。
- 周期species波と周期速度波を拡散・粘性で減衰させる。
- 輸送を含む周期計算でも全species質量、運動量、全エネルギーを保存する。
- Stage 6の周期反応流で総質量、運動量、全エネルギーを保存する。
- 一様反応場のStage 6結果がStage 5の化学更新と一致する。
- 反応速度ゼロのStage 6結果がStage 4の流体更新と一致する。
- Dirichlet境界が指定した密度、速度、圧力、全species質量分率を再現する。
- 一様参照状態の無反射境界が一様場を変化させない。
- 鏡像境界が法線運動量だけを反転し、全エネルギーを保存する。
- 非周期境界を含む粘性・反応流の右辺と時間発展が有限かつ正値である。
- Stage 7のstep 0、指定間隔、最終stepスナップショットと積分履歴を生成する。
