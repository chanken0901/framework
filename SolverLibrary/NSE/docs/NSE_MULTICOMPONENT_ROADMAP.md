# NSE多成分・反応流拡張仕様

更新日: 2026-09-02

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
`case.yaml`へ指定する。係数配列は`physics.multicomponent.species`の名前をキーにするため、
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

## 7. 保存変数

全Stageで状態レイアウトを共有する。

```text
Q = [rho_1, ..., rho_Ns, rho*u, rho*v, rho*w, rho*E]
rho = sum_s rho_s
nvar = Ns + 4
```

`Ns=1`では5保存変数になる。変数番号は`mod_mc_state_layout`だけが決定する。

## 8. Provider境界

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
化学反応providerは全Stageで`none`である。入力名とコンパイル済みproviderが異なる場合は
開始前に停止する。

## 9. 設計書

完全なひな型は
`ScriptLibrary/RunEnvironment/case_templates/nse_multicomponent_inviscid.yaml`
に置く。

```yaml
physics:
  model: nse_multicomponent
  multicomponent:
    mode: inviscid_euler
    species: [species_a, species_b]

thermodynamics:
  model: calorically_perfect
  gamma: 1.4

transport:
  model: none

chemistry:
  model: none

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

Stage 3では`case_templates/nse_multicomponent_thermally_perfect.yaml`をひな型にする。
物性はspecies名をキーにして指定する。

```yaml
physics:
  multicomponent:
    mode: thermally_perfect_euler
    species: [N2, O2]

thermodynamics:
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
輸送物性もspecies名をキーにするため、species一覧との過不足を入力生成時に拒否する。

```yaml
physics:
  multicomponent:
    mode: viscous_navier_stokes
    species: [N2, O2]

transport:
  model: mixture_averaged
  reference_dynamic_viscosity: 1.8e-5
  prandtl_number: 0.72
  species_data:
    N2: {diffusivity: 2.0e-5}
    O2: {diffusivity: 2.0e-5}

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

## 10. ビルドと実行環境

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

## 11. Stage 4の制約

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

未実装の選択肢は黙って別方式として扱わず、入力生成時または計算開始前に拒否する。

## 12. 後続Stage

5. 0次元有限反応速度化学
6. Strang分割による流体・反応結合
7. 反応流境界条件、初期条件、出力
8. MPI/OpenMP最適化
9. 一般座標
10. CUDA

## 13. 必須回帰条件

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
- 成分拡散流束の総和を各セルでゼロにする。
- 一様場の輸送右辺をゼロにする。
- 周期species波と周期速度波を拡散・粘性で減衰させる。
- 輸送を含む周期計算でも全species質量、運動量、全エネルギーを保存する。
