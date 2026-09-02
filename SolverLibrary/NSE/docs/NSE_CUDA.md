# NSE CUDA版（単一GPU／MPI＋マルチGPU）

> Windows（PowerShell）とLinux（bash）の共通コマンド対応は[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。cuFFTMpはLinux専用です。

## 対応範囲

`cuda_single` profileと`cuda_mpi` profileは、CPU版と同じ保存変数
`rho, rho_u, rho_v, rho_w, rho_E`を使い、次の組み合わせをGPUで計算します。

- 対流項: `keep2`、`keep6`、特性空間`weno5z_roe`、またはKEEP/WENO `hybrid`
- 粘性項: 6次精度中心差分 `central6`、または `none`
- 境界条件: 6面別の周期境界、特性無反射境界、鏡像境界、またはDirichlet固定状態境界
- 時間積分: SSPRK3
- 精度: float64
- `cuda_single`: GPU 1台、MPI/OpenMPなし
- `cuda_mpi`: MPI 1 rankにつきGPU 1台、OpenMPなし

`central6`は一定粘性係数のNewton流体、Stokesの仮定、Fourier熱伝導を
CUDAカーネルで計算します。対流方式、粘性方式、面別境界、SSPRK3は
単一GPU版とマルチGPU版で共通です。

## GPU化される処理

初期条件はCPUで生成し、その後、保存変数をGPUへ一度転送します。時間発展中は
次のデータをGPUメモリへ常駐させます。

- 現在値 `Q`
- SSPRK3の開始値 `Q0`
- 右辺 `RHS`
- 粘性計算用の基本変数 `u, v, w, T`
- CFL評価用の局所最大波速度

各時間ステップでは、周期／無反射／鏡像／Dirichlet ghost更新、CFL・粘性時間刻み評価、選択した対流流束、
三方向の流束発散、粘性・熱伝導項、SSPRK3の3段階更新をCUDAカーネルで
実行します。SLF出力時だけ`Q`をCPUへ戻します。

GPU版は面流束配列を保存せず、各セルで必要な左右6面の流束を直接計算します。
これにより、三方向分の大きな流束配列をGPUメモリへ保持しません。

## MPI＋CUDAの実装

`cuda_mpi`は、既存CPU版と同じくx方向を各rankが全域保持し、y-z面を
二次元MPI分割します。各rankは自分の局所領域だけをGPUに保持します。

SSPRK3の各段で、次の順序でghostセルを完成させます。

1. CUDAカーネルでx面、非分割方向、分割方向の物理端に必要な周期／無反射／鏡像／Dirichlet ghostを更新
2. y面をGPUからhost staging bufferへpackし、内部隣接rankまたは周期反対端rankとMPI交換してGPUへunpack
3. z面を同様に交換。無反射／鏡像物理端では通信せずGPU生成値を保持し、辺・角ghostを完成
4. 局所RHSとSSPRK3段更新をGPUで実行

CFL時間刻みは各GPUの局所値を計算した後、`MPI_Allreduce(MIN)`で全rankを
同じ値にします。SLF出力はCPU版と同じrank別形式で、`meta.json`には全rankの
局所範囲が記録されます。通信はhost staging方式なのでCUDA-aware MPIは必須では
ありません。

各y/z局所ブロックは少なくとも`nghost`セル必要です。GPU割当ての既定は
MPI shared-memory communicatorから得たノード内rank番号です。Slurmなどが
各rankへ`CUDA_VISIBLE_DEVICES`を1台だけ公開する場合は、そのrankから見える
device 0を自動選択します。1台のGPUを複数rankで共有するデバッグ時だけ、
`NSE_CUDA_DEVICE_POLICY=fixed`を設定できます。

`cuda_mpi`はTaylor–Green初期条件、`imported_turbulence`、
`shock_turbulence_interaction`、`shock_tube_turbulence_interaction`を利用します。
初期条件はCPUで共通実装から生成した後、各GPUへ転送されるため、有限高圧室ケースも
CPU版と同じ配置になります。入力は
[`NSE_SHOCK_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TURBULENCE_INTERACTION.md)および
[`NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)を参照してください。
分散HIT初期化またはPetersen–Livescu forcingが必要なLinux計算機では、
`cuda_mpi_cufftmp`を選択します。このプロファイルは既存のY-Z領域分割を
`cufftMpMakePlanDecomposition`へ直接渡します。必要環境とビルド方法は
[`NSE_CUFFTMP.md`](NSE_CUFFTMP.md)を参照してください。

## 実行環境の生成

CPU/MPI版と共通のNSE設計書テンプレートをコピーします。

```powershell
$framework = "C:\Users\Owner\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"

New-Item -ItemType Directory -Force $designs | Out-Null
Copy-Item `
  "$tool\environment.nse.yaml" `
  "$designs\nse_tgv_cuda.yaml"
```

Linux（bash）:

```bash
framework="$HOME/Research/FrameWork"
tool="$framework/ScriptLibrary/RunEnvironment"
designs="$HOME/ResearchRuns/Designs"

mkdir -p "$designs"
cp "$tool/environment.nse.yaml" "$designs/nse_tgv_cuda.yaml"
```

LinuxではCUDAのparallel設定に加えて、`source.framework_root`、`destination`、
`select.target`をLinux用へ変更します。cuFFTMpを使う場合は
[`NSE_CUFFTMP.md`](NSE_CUFFTMP.md)の専用profileと要件も確認してください。

コピーした`nse_tgv_cuda.yaml`の選択を次のように変更します。

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

MPI＋CUDA（マルチGPU）では、同じ設計書を次のようにします。

```yaml
select:
  model: nse
  execution: release

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: true
```

```powershell
python "$tool\prepare_environment.py" `
  "$designs\nse_tgv_cuda.yaml" --dry-run
python "$tool\prepare_environment.py" `
  "$designs\nse_tgv_cuda.yaml"
```

Linux（bash）:

```bash
python3 "$tool/prepare_environment.py" \
  "$designs/nse_tgv_cuda.yaml" --dry-run
python3 "$tool/prepare_environment.py" \
  "$designs/nse_tgv_cuda.yaml"
```

生成された`nse_caseNNNN`へ移動し、通常の共通フローを実行します。

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

単一GPU版ではMPIとOpenMPを無効にし、マルチGPU版ではMPIだけを有効にします。
profile、MPI、CUDAはenvironment設計から自動決定されるため、`case.yaml`には
重ねて記述しません。単一GPUのcase設定例は次のとおりです。

```yaml
solver:
  use_openmp: false
  mpi_processes: 1
  omp_threads: 1
  cuda_device: 0
```

マルチGPUでは`mpi_processes`をGPU総数に合わせます。通常は1 rank＝1 GPUです。

```yaml
solver:
  use_openmp: false
  mpi_processes: 4
  omp_threads: 1
  cuda_device: 0
```

対流流束はCPU版と共通のケース設定で選びます。

```yaml
numerics:
  convective_scheme: weno5z_roe  # keep2, keep6, weno5z_roe, hybrid
```

CPU版とCUDA版のどちらも`keep`単独の指定は使用できません。

特定方向だけを無反射にする場合もCPU版と同じ`boundary`設定を使います。例えば
x両端を無反射、y-z両端を周期にする設定は次のとおりです。

```yaml
boundary:
  faces:
    x_min: {type: non_reflecting, reference_state: far_field}
    x_max: {type: non_reflecting, reference_state: far_field}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
  reference_states:
    far_field:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143
  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto
```

周期面は同じ方向の両端を対で指定します。Petersen–Livescu FFT forcingは物理上の
前提から全6面周期を要求するため、無反射面または鏡像面との併用は入力生成時に拒否されます。

鏡像境界もCPU版と同じ`type: reflective`で指定します。密度、接線運動量、全エネルギーは
対応する内点からそのまま鏡映し、面法線方向の運動量だけ符号を反転します。これは
自由滑り・断熱条件であり、no-slip壁ではありません。MPI＋CUDAの物理端では通信せず、
GPU上で鏡像ghostを作ります。

## ライブラリ上で直接ビルドする場合

BuildSolverを使うと、WindowsのVisual Studio x64環境も自動設定されます。

```powershell
Set-Location `
  "C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\BuildSolver"

python .\build_model.py .\build.yaml `
  --model nse `
  --profile cuda_single `
  --build
```

Linux（bash）:

```bash
cd "$HOME/Research/FrameWork/ScriptLibrary/BuildSolver"
python3 ./build_model.py ./build.yaml \
  --model nse \
  --profile cuda_single \
  --build
```

テストも含める場合:

```powershell
python .\build_model.py .\build.yaml `
  --model nse `
  --profile cuda_single `
  --test
```

Linux（bash）:

```bash
python3 ./build_model.py ./build.yaml \
  --model nse \
  --profile cuda_single \
  --test
```

実行ファイル名は`nse_cuda.exe`です。

MPI＋CUDA版:

```powershell
python .\build_model.py .\build.yaml `
  --model nse `
  --profile cuda_mpi `
  --build
```

Linux（bash）:

```bash
python3 ./build_model.py ./build.yaml \
  --model nse \
  --profile cuda_mpi \
  --build
```

実行ファイル名は`nse_mpi_cuda.exe`です。例えば4 GPUでは次のように実行します。

```powershell
mpiexec -n 4 .\nse_mpi_cuda.exe .\input.dat
```

Linux（bash）:

```bash
mpirun -np 4 ./nse_mpi_cuda ./input.dat
```

Slurmでは施設のMPI構成に従い、`mpirun`の代わりに`srun`を使う場合があります。

## GPUアーキテクチャ

使用GPUに合わせて、machine profileの`cuda_architectures`を設定します。

```yaml
libraries:
  cuda_compiler: C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.3/bin/nvcc.exe
  cuda_toolkit_root: C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.3
  cuda_architectures: 86
```

値はGPUのcompute capabilityに対応します。実機と異なる値を指定すると、
実行時に「unsupported toolchain」や「no kernel image」のエラーになる場合が
あります。

## メモリ使用量

主なGPUメモリ量は、おおよそ次式です。

```text
3 * Ntotal * 5 * 8 byte
+ Ntotal * 4 * 8 byte
+ reduction workspace

Ntotal = (nx+2*nghost) * (local_ny+2*nghost) * (local_nz+2*nghost)
```

単一GPUでは`local_ny=ny`、`local_nz=nz`です。MPI＋CUDAでは各GPUの局所y-z
範囲を使います。このほかにCUDA reduction領域と、最大y/z面1組分のGPU・host
staging bufferが必要です。

## ParaView

CUDA版もCPU版と同じSLF形式を出力します。

```powershell
python .\tools\postprocess_case.py
```

Linux（bash）:

```bash
python3 ./tools/postprocess_case.py
```

生成された`cases\<case_id>\paraview\collection.pvd`をParaViewで開きます。

## 検証内容

実装時に次を確認しています。

1. CUDA 13.3とGNU FortranによるWindowsビルド
2. 一定流がSSPRK3の1ステップ後も保存されること
3. 周期ghostセルが反対側の物理セルと一致すること
4. GPUのCFL時間刻みがCPU計算と一致すること
5. KEEP2、KEEP6、WENO5-Z/Roe、KEEP/WENOハイブリッドの三次元SSPRK3結果がCPU参照実装と一致すること
6. 非一様な三次元場の`central6`粘性項と時間刻みがCPU参照実装と一致すること
7. 8立方Taylor-Green渦のKEEP、WENO5-Z/Roe、ハイブリッド 2ステップsmoke test
8. 従来の`cpu_mpi` profileが引き続きビルドできること
9. `cuda_mpi`のWindows MPI＋CUDAビルドと2 rank／4 rank smoke test
10. 2×2 y-z分割の全物理セル2,560値が単一GPU結果と完全一致すること
11. 全6面無反射およびx無反射・y-z周期で、辺・角を含むCPU/CUDA ghost値が許容誤差内で一致すること
12. 単一GPUと2 rank MPI＋CUDAで無反射境界を含む時間発展が正常終了すること
13. 全6面鏡像で、面・辺・角のCPU/CUDA ghost値が許容誤差内で一致すること
14. 単一GPUと2 rank MPI＋CUDAで鏡像境界を含む時間発展が正常終了すること
