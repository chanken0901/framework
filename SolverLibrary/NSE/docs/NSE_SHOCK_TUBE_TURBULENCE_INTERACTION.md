# 有限高圧室から発生する衝撃波–乱流干渉

> Windows（PowerShell）とLinux（bash）の実行コマンド対応は[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

更新日: 2026-09-02
対応: CPU/MPI/OpenMP、単一 GPU CUDA、MPI＋CUDA

## 1. 概要

`flow.type: shock_tube_turbulence_interaction`は、閉じた有限長の高圧室と低圧室の
隔膜を時刻0で取り除く衝撃波管初期条件です。低圧域の下流には保存済みの局所乱流を
配置します。

```text
x_min                                                           x_max
  | 高圧driver | 隔膜 | 低圧driven | 局所乱流 | 低圧driven |
  | reflective | x_d  |                         non_reflecting |
                  衝撃波・接触面 --->
              <--- 膨張波（x_minで反射後、右へ進む）
```

これは`低圧－高圧－低圧－乱流－低圧`という配置から左側低圧域を省き、高圧室内の
切断面より右側だけを計算する構成です。`x_min`は高圧室の閉端または鏡像面です。
省略した左側を波が通過する問題ではなく、`x_min`で反射する問題になります。高圧状態を
Dirichlet境界から供給し続けないため、膨張波が反射して衝撃波へ追いつく有限ドライバー
問題を計算できます。隔膜除去直後には衝撃波と膨張波に加えて接触不連続面も生じます。

## 2. そのまま使える設定例

次の部分をNSEの`case.yaml`へ設定します。座標は対象格子のセル境界に合わせてください。

```yaml
flow:
  type: shock_tube_turbulence_interaction

  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: embed
    x_start: 4.0
    blend_cells: 0
    # 省略時はshock_tube.driven.velocityと同じ値。
    velocity_offset: [0.0, 0.0, 0.0]
    # 省略時はshock_tube.drivenと同じ状態。
    background:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.7142857142857143

  shock_tube:
    diaphragm_position: 2.0
    driver:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 3.5714285714285716
    driven:
      density: 1.0
      velocity: [0.0, 0.0, 0.0]
      pressure: 0.7142857142857143

boundary:
  faces:
    x_min: {type: reflective}
    x_max:
      type: non_reflecting
      reference_state: driven_state
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}

  reference_states:
    driven_state:
      source: shock_tube.driven

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

## 3. 状態と配置の制約

- `driver.pressure`は`driven.pressure`より大きくします。衝撃波Mach数はこのRiemann
  問題から決まり、平面衝撃波ケースのように直接指定しません。
- `driver.velocity`のx成分は0にします。これは静止した閉端高圧室と`reflective`
  境界を整合させるためです。
- `diaphragm_position`は`x_min < x_d < x_max`のxセル境界上に置きます。
- 乱流は`mode: embed`で、`diaphragm_position <= x_start`とします。
- `imported_turbulence.background`は`shock_tube.driven`と同じ状態にします。省略すれば
  自動的に同じ値になります。
- `x_min`は`reflective`、`x_max`は`non_reflecting`にし、後者の参照状態には
  `source: shock_tube.driven`を使います。
- y、zを周期にする場合、min/maxを必ず対で指定します。
- SLFと対象格子の`ny`、`nz`、`dx`、y-z範囲を一致させます。SLF作成時と実行時の
  MPIプロセス数、CPU/CUDA方式は一致させる必要がありません。

これらは`input.dat`生成時とFortran初期化時の両方で検査されます。

## 4. 膨張波が衝撃波へ追いつく条件

高圧室長さ`L_driver = diaphragm_position - x_min`を短くすると、左向き膨張波が
`x_min`へ到達して反射するまでの時間が短くなります。その結果、反射膨張波が
右向き衝撃波へ早く追いつきます。圧力比、比熱比、隔膜から乱流までの距離も到達位置を
変えるため、まず粗い格子で密度・圧力のx方向平均を追跡し、乱流へ到達する前後の
衝撃波強度を確認してください。

`x_min`を鏡像面と解釈する場合、全領域には反対側にも同じ構造が鏡映されます。物理的な
剛体閉端として解釈する場合、現行`reflective`は自由滑り・断熱壁です。粘性no-slip壁では
ありません。

## 5. 数値設定

衝撃波、接触面、膨張波を含むため、`weno5z_roe`または`hybrid`を使用してください。
滑らかな乱流の精度を保ちながら不連続を捕獲する用途では、例のKEEP6/WENO5-Z–Roe
ハイブリッドを推奨します。初回は`cfl: 0.2`～`0.3`から開始します。現行ソルバーには
positivity-preserving limiterがないため、強い圧力比では格子解像度とCFLを段階的に
調整してください。

## 6. SLF作成から実行まで

元の乱流計算のrank別出力を、ghostなしの可搬SLFへ変換します。

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

生成済みケースでは`case.yaml`を編集した後、必ず入力を再生成します。

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

正常終了時は最終stepの出力後に終了メッセージが表示されます。

## 7. 実装箇所

- `src/init/mod_init_shock_tube_turbulence.f90`: 低圧背景・乱流・有限高圧室の配置と検査
- `src/init/mod_init_imported_turbulence.f90`: バックエンド非依存の可搬SLF読込み
- `src/boundary/mod_boundary_runtime.f90`: CPU/MPI/OpenMPの鏡像・無反射境界
- `src/gpu/nse_cuda_bridge.cu`: CUDAの鏡像・無反射境界
- `ScriptLibrary/RunEnvironment/case_input.py`: YAML解決と早期検査
- `tests/test_imported_turbulence.f90`: 初期配置の回帰試験
- `ScriptLibrary/RunEnvironment/tests/test_case_input.py`: 正常入力と不整合拒否の試験
