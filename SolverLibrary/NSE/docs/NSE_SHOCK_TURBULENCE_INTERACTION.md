# 平面衝撃波–乱流干渉の初期条件と流入境界

> Windows（PowerShell）とLinux（bash）の実行コマンド対応は[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

更新日: 2026-09-02
対応: CPU/MPI/OpenMP、単一 GPU CUDA、MPI＋CUDA

## 1. 概要

`flow.type: shock_turbulence_interaction` は、保存済みの乱流場を長い x 方向領域へ埋め込み、乱流の外側に平面衝撃波を置く初期条件です。衝撃波は時間発展により乱流へ入射します。衝撃波の後方にある物理量は、進行方向上流側の x 境界に Dirichlet 条件として継続的に与えられます。

正の x 方向へ伝播させる場合の配置は次のとおりです。

```text
x_min                                                      x_max
  |  衝撃波背後状態  | 衝撃波 | 一様な衝撃波前方状態 | 乱流ブロック | 前方状態 |
  |<-- Dirichlet ---|  x_s  |                    x_start              |
                               衝撃波の進行方向 --->
```

負の x 方向では左右が逆になり、`x_max` が Dirichlet 流入面になります。

## 2. そのまま使える case.yaml 設定例

次の例は Mach 1.5 の衝撃波を正の x 方向へ進めます。`position` と `x_start` は実際の格子に合わせて変更してください。

```yaml
flow:
  type: shock_turbulence_interaction

  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: embed
    x_start: 4.0
    blend_cells: 0
    # 元SLFが変動速度だけを持つ場合に加える平均速度。
    # 省略時は planar_shock.upstream.velocity と同じ値になる。
    velocity_offset: [0.0, 0.0, 0.0]
    # 省略時は planar_shock.upstream と同じ状態になる。
    background:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.7142857142857143

  planar_shock:
    position: 2.0
    propagation_direction: positive_x
    mach_number: 1.5
    upstream:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.7142857142857143

boundary:
  faces:
    x_min:
      type: dirichlet
      reference_state: shock_downstream
    x_max:
      type: non_reflecting
      reference_state: shock_upstream
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}

  reference_states:
    shock_downstream:
      source: planar_shock.downstream
    shock_upstream:
      source: planar_shock.upstream

  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto

forcing:
  type: none

numerics:
  convective_scheme: hybrid
  hybrid:
    smooth_scheme: keep6
    shock_scheme: weno5z_roe
    sensor: ducros_pressure
    sensor_onset: 0.01
    sensor_full: 0.10
  viscous_scheme: central6
  time_integrator: ssprk3

time:
  cfl: 0.25
  use_fixed_dt: false
```

`source: planar_shock.downstream` を使うと、Rankine–Hugoniot 関係から求めた衝撃波背後状態が Dirichlet 境界へ自動的に渡されます。同じ値を二重入力しないため、初期条件と境界条件の不整合を防げます。

## 3. 衝撃波背後状態の決定

### 3.1 衝撃 Mach 数から自動計算する方法

`mach_number > 1`、比熱比 `gamma`、衝撃波前方の密度・速度・圧力を与えると、完全気体の一次元垂直衝撃波関係から背後状態を計算します。

```text
r = rho2/rho1
  = ((gamma + 1) M_s^2) / ((gamma - 1) M_s^2 + 2)

p2/p1
  = 1 + (2 gamma/(gamma + 1)) (M_s^2 - 1)

c1 = sqrt(gamma p1/rho1)
u2 = u1 + s M_s c1 (1 - 1/r)
```

ここで `s=+1` は `positive_x`、`s=-1` は `negative_x` です。接線方向速度は変えません。`M_s` は衝撃波前方流体に対する衝撃波 Mach 数です。

### 3.2 背後状態を直接指定する方法

実験値や別の Riemann 問題から得た状態を使用する場合は、`mach_number` を削除し、次のように指定します。

```yaml
  planar_shock:
    position: 2.0
    propagation_direction: positive_x
    upstream:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.7142857142857143
    downstream:
      density: 1.8620689655172413
      velocity: [0.6944444444444444, 0.0, 0.0]
      pressure: 1.755952380952381
```

`mach_number` と `downstream` の同時指定、および両方の省略はエラーです。

## 4. 配置条件

- `planar_shock.position` は `x_min < position < x_max` を満たし、対象格子のセル境界 `x_min + n*dx` 上になければなりません。
- `positive_x` では衝撃波を乱流ブロックの左外側へ置き、`x_min` を Dirichlet にします。
- `negative_x` では衝撃波を乱流ブロックの右外側へ置き、`x_max` を Dirichlet にします。
- 乱流読込みモードは `embed` または `periodic_embed` です。後者は`x_length`で区間長を指定します（[配置手順](NSE_IMPORTED_TURBULENCE.md)）。`tile` は使用できません。
- 読込み SLF と対象格子の `ny`、`nz`、`dx` は一致させます。実行時の MPI プロセス数や CPU/CUDA の種類は、SLF 作成時と一致させる必要はありません。
- `imported_turbulence.background` は衝撃波前方状態と一致させます。省略すれば自動的に一致します。
- `blend_cells: 0` は保存乱流をそのまま配置します。滑らかな接続が必要なら正の値を使用できますが、乱流ブロック端の変動が減衰します。

配置に重なりや不整合がある場合、`input.dat` 生成時または初期化時に停止し、問題の面・状態・位置を表示します。

## 5. Dirichlet 境界の数値仕様

`type: dirichlet` は、指定した原始変数 `(rho,u,v,w,p)` を保存変数 `(rho,rho*u,rho*v,rho*w,E)` に変換し、該当する物理境界の全 3 ghost 層へ毎回固定します。

```text
E = p/(gamma - 1) + rho (u^2 + v^2 + w^2)/2
```

これは面上の補間値だけを制約する条件ではなく、境界外側に一定状態の流体リザーバを置く条件です。そのため、衝撃波背後状態を領域へ供給し続ける本用途に適しています。CPU、CUDA、MPI 内部 halo を含む面・辺・角の処理順序は共通です。

## 6. 数値流束の選択

衝撃波は不連続なので、`keep6` 単独は推奨しません。次のいずれかを使用してください。

- `weno5z_roe`: 全領域で衝撃波捕獲を優先する。
- `hybrid`: 滑らかな乱流領域では KEEP6、衝撃波近傍では WENO5-Z/Roe を使う。

最初の動作確認では `cfl: 0.2`～`0.3` 程度から開始し、密度・圧力の最小値と衝撃波位置を確認してから増やしてください。現在は positivity-preserving limiter を実装していないため、強い衝撃波では格子解像度、CFL、流束選択の影響を受けます。

## 7. 作成から実行まで

まず、既存の rank 別出力から ghost を除いた可搬 SLF を作ります。

```powershell
python .\SolverLibrary\NSE\tools\nse_prepare_imported_turbulence.py `
  .\previous_case\output `
  --step latest `
  --output .\cases\caseNNNN\initial_data\turbulence.slf
```

Linux（bash）:

```bash
python3 ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py \
  ./previous_case/output \
  --step latest \
  --output ./cases/caseNNNN/initial_data/turbulence.slf
```

生成済み実行環境では `case.yaml` を編集した後、次の順で再生成・検証・ビルド・実行します。

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

Linux（bash）:

```bash
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --run
```

正常終了時は最終 step の出力後に次のメッセージが表示されます。

```text
NSE calculation completed successfully: step=<final step>, time=<final time>
```

## 8. 実行前チェック

1. `turbulence.slf` が case ディレクトリ基準のパスに存在する。
2. 衝撃波位置がセル境界上で、乱流ブロックの外側にある。
3. 進行方向が `positive_x` なら `x_min`、`negative_x` なら `x_max` が Dirichlet である。
4. Dirichlet の参照状態に `source: planar_shock.downstream` を使用している。
5. 反対側の x 面は通常、前方状態を参照する `non_reflecting` である。
6. y、z の周期境界は min/max を対で指定している。
7. `forcing.type: none` である。FFT forcing は全 6 面周期を必要とするため併用できない。
8. `convective_scheme` が `weno5z_roe` または `hybrid` である。

## 9. 実装と試験

主な実装は次のファイルにあります。

- `src/init/mod_init_shock_turbulence.f90`: 乱流配置、衝撃波前後状態の配置、整合性検査
- `src/boundary/mod_boundary_runtime.f90`: CPU/MPI/OpenMP Dirichlet 境界
- `src/gpu/nse_cuda_bridge.cu`: CUDA Dirichlet 境界
- `ScriptLibrary/RunEnvironment/case_input.py`: YAML 解決、Rankine–Hugoniot 計算、早期検査
- `tests/test_imported_turbulence.f90`: 衝撃波と乱流の初期配置
- `tests/test_boundary_dirichlet.f90`: MPI の全 ghost 層、面・辺・角
- `tests/test_cuda_boundary_dirichlet.f90`: CPU/CUDA 一致
- `tests/input_tgv_dirichlet_small.dat`: CPU、CUDA、MPI＋CUDA smoke test
