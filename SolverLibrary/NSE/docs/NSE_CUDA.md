# NSE 単一GPU CUDA版

## 対応範囲

`cuda_single` profileは、現在のCPU版と同じ保存変数
`rho, rho_u, rho_v, rho_w, rho_E`を使い、次の組み合わせをGPUで計算します。

- 対流項: `keep2`、`keep6`、特性空間`weno5z_roe`、またはKEEP/WENO `hybrid`
- 粘性項: 6次精度中心差分 `central6`、または `none`
- 境界条件: 三方向周期境界
- 時間積分: SSPRK3
- 精度: float64
- GPU数: 1
- MPI/OpenMP: 使用しない

`central6`は一定粘性係数のNewton流体、Stokesの仮定、Fourier熱伝導を
CUDAカーネルで計算します。現在の`cuda_single`は単一GPU専用です。

## GPU化される処理

初期条件はCPUで生成し、その後、保存変数をGPUへ一度転送します。時間発展中は
次のデータをGPUメモリへ常駐させます。

- 現在値 `Q`
- SSPRK3の開始値 `Q0`
- 右辺 `RHS`
- 粘性計算用の基本変数 `u, v, w, T`
- CFL評価用の局所最大波速度

各時間ステップでは、周期ghost更新、CFL・粘性時間刻み評価、選択した対流流束、
三方向の流束発散、粘性・熱伝導項、SSPRK3の3段階更新をCUDAカーネルで
実行します。SLF出力時だけ`Q`をCPUへ戻します。

GPU版は面流束配列を保存せず、各セルで必要な左右6面の流束を直接計算します。
これにより、三方向分の大きな流束配列をGPUメモリへ保持しません。

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

コピーした`nse_tgv_cuda.yaml`の選択を次のように変更します。

```yaml
select:
  model: nse_cuda_single
  execution: release

parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

```powershell
python "$tool\prepare_environment.py" `
  "$designs\nse_tgv_cuda.yaml" --dry-run
python "$tool\prepare_environment.py" `
  "$designs\nse_tgv_cuda.yaml"
```

生成された`nse_caseNNNN`へ移動し、通常の共通フローを実行します。

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

CUDA版ではMPIとOpenMPを無効にします。GPU番号は`case.yaml`の次の項目で
選択できます。

```yaml
solver:
  profile: cuda_single
  use_mpi: false
  use_openmp: false
  use_cuda: true
  mpi_processes: 1
  omp_threads: 1
  cuda_device: 0
```

対流流束はCPU版と共通のケース設定で選びます。

```yaml
numerics:
  convective_scheme: weno5z_roe  # keep2, keep6, weno5z_roe, hybrid
```

CPU版とCUDA版のどちらも`keep`単独の指定は使用できません。

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

テストも含める場合:

```powershell
python .\build_model.py .\build.yaml `
  --model nse `
  --profile cuda_single `
  --test
```

実行ファイル名は`nse_cuda.exe`です。

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

Ntotal = (nx+2*nghost) * (ny+2*nghost) * (nz+2*nghost)
```

256立方格子、`nghost=3`では、主要配列だけで約2.60 GiBです。このほかに
CUDA reduction用の小さな作業領域が必要です。

## ParaView

CUDA版もCPU版と同じSLF形式を出力します。

```powershell
python .\tools\postprocess_case.py
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
