# Stage 10: GPU常駐型の多成分・反応流

更新日: 2026-09-08

## 対象と実装状態

`nse_multicomponent` のStage 6--9の反応性Navier--StokesをCUDAで実行する。
単成分 `nse` のCUDAソース・計算法・通信方式は変更していない。

| profile | GPU | MPI | OpenMP |
|---|---|---|---|
| `cuda_single_reactive` | 1 GPU | なし | なし |
| `cuda_mpi_reactive_pencil` | 各rankのGPU | xを保持しy/zを分割 | なし |

既存CPU profileは引き続き選択可能。CUDAは `reactive_navier_stokes` モード用であり、
旧Stage 0--5用profileを置換しない。CUDA機器や実行条件が不適合ならエラーにし、
CPU版への暗黙のフォールバックはしない。

## GPUに常駐するもの

- species部分密度、Cartesian運動量3成分、全エネルギー（`Ns+4`変数）。
- RK開始時の保存変数、右辺、面流束、primitive値、物理勾配、halo用作業領域。
- NASA7物性、輸送係数、反応係数、境界参照状態と格子パラメータ。

初期化時に場を一度アップロードした後、熱力学復元、一次Rusanov対流流束、
粘性・熱伝導・species拡散、境界条件、SSPRK3、Strang分割の化学反応subcyclingを
GPU上で計算する。ノズルの計量もGPU上で評価する。
毎段階の全配列ダウンロード／再アップロードは行わない。
NASA7・輸送・反応データはCPU providerが検証した係数を渡すため、GPU用の別物性入力は不要。

CPU側には実行制御、ファイル入出力、MPI呼び出しが残る。
時間刻み、エラー状態、substep数、保存量・最小値はGPUで集約し、少量の値だけCPUへ戻す。
MPI版では全rankの最小時間刻みを使う。全状態をCPUへ戻すのはsnapshot／履歴／最終出力時。
出力用のCPU配列も保持するため、「GPU常駐」はCPUメモリが不要という意味ではない。

## 維持する物理・数値仕様

- NASA7温度・組成依存理想気体混合物。
- 既存の定粘性・定Prandtl・species別定拡散係数の輸送モデル。
- 質量保存補正付きspecies拡散と一段不可逆Arrhenius反応。
- 化学反応の陽的SSPRK3 subcycling、流体とのStrang分割。
- `cartesian` と静止押出し平面 `planar_nozzle`。
- 面ごとの `periodic`、`reflective`、`dirichlet`、`non_reflecting`。
- 流体CFL・拡散制約、化学反応のdepletion制約、NASA温度範囲・有限値・正値性検査。
- 周期領域の最終質量・運動量・全エネルギー保存検査。

数値モデル自体はCPU版と同じである。多成分側の対流は一次Rusanovであり、
単成分側のKEEP6/WENOを新たに移植したものではない。
鏡像壁は自由滑り・断熱、無反射は既存の固定平均比熱比によるcharacteristic-relaxation近似。
任意外部格子、軸対称、非滑り壁、指定壁温、壁面反応、複数反応・stiff陰解法は含まない。

## 設計書の指定

単一GPU用の完成例は
`ScriptLibrary/RunEnvironment/environment.nse_multicomponent.cuda.yaml`。
既存の `case.yaml` と `config/*.yaml` の分離仕様は変わらない。

```yaml
parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
solver:
  profile: cuda_single_reactive
  include_tests: true
```

MPI＋CUDAでは上記を次へ変更して環境を生成する。

```yaml
parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: true
solver:
  profile: cuda_mpi_reactive_pencil
  include_tests: true
```

生成後の `cases/caseXXXX/case.yaml` のsolver節では、例えば4 rankなら以下を設定する。
`use_mpi` は従来どおり維持する。profileとの矛盾は入力生成時に拒否する。

```yaml
solver:
  profile: cuda_mpi_reactive_pencil
  use_mpi: true
  use_openmp: false
  use_cuda: true
  mpi_processes: 4
  omp_threads: 1
  decomposition: pencil
  process_grid: [2, 2]
```

`process_grid: [0, 0]` ならMPIが分割数を決める。
分割する各方向で、各rankの所有セルは少なくとも2個必要。

## Windowsで生成・実行

GNU Fortran、CMake、Ninja、CUDA Toolkit、対応GPUドライバ、Visual Studio C++ Build Toolsを用意する。
MPI版にはMicrosoft MPI SDKとランタイムも必要。
machine profileのCUDAコンパイラ・Toolkit・architectureは実機に合わせる。

FrameWorkのルートで実行する。既存の実行環境を上書きする指定は付けない。

```powershell
python .\ScriptLibrary\RunEnvironment\prepare_environment.py `
  .\ScriptLibrary\RunEnvironment\environment.nse_multicomponent.cuda.yaml
```

表示された生成先へ移動し、必要ならcase/configを編集して実行する。

```powershell
Set-Location '表示された生成先の絶対パス'
python .\tools\run_case.py --prepare
python .\tools\run_case.py --build --run
```

## Linuxで生成・実行

Linux用の設計書を作るときは、上記完成例をコピーし、`select.target` を
`linux_gnu_mpi`、`select.destination` を `local_generated` にする。
トップレベルに生成先と索引のパスを明記する（`YOUR_USER` は実際のユーザー名へ変更）。

```yaml
destination:
  root: /home/YOUR_USER/ResearchRuns/mc_cuda_case0001
  case_index: /home/YOUR_USER/ResearchRuns/case_index.csv
```

CUDA、GNU Fortran、MPI、CMake、Ninjaの環境を利用先の手順に従って読み込む。
Linux machine profileの `compiler.fortran` は既定で `mpifort`。
`libraries.cuda_architectures` 等も対象GPUに合わせる。
FrameWorkルートで、自分の設計書のパスを指定して生成する。

```bash
python3 ./ScriptLibrary/RunEnvironment/prepare_environment.py ./my_cuda_environment.yaml \
  --framework-root "$PWD"
cd /home/YOUR_USER/ResearchRuns/mc_cuda_case0001
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --build --run
```

バッチ計算ではサイト指定のGPU割当て・MPI起動方法を使う。
この変更で特定スパコンのジョブ設定を自動認定したわけではない。

## GPU割当てと通信

単一GPU版は可視GPUのdevice 0、MPI版はノード内MPI rankをdevice番号として選ぶ。
標準の運用は1 MPI rank / GPU。可視GPUが不足すれば開始時にエラーにする。
スケジューラが各rankに異なるGPUを1台だけ見せる設定では、各rankの可視device 0を指定する。

```powershell
$env:NSE_MC_CUDA_DEVICE = '0'
python .\tools\run_case.py --run
```

```bash
export NSE_MC_CUDA_DEVICE=0
python3 ./tools/run_case.py --run
```

全rankから同じGPU集合が見える状態でこの指定を行うと、全rankが同じGPUを共有する。
複数GPUへ自動分散する指定ではない。GPUの可視化とrank割当てを必ず確認する。

今回はGPU間の直接MPI通信は実装していない。
GPU上で2層haloをpackし、halo部分だけCPUバッファ経由でMPI通信してGPUにunpackする。
y方向の通信後にコーナーを含めたz方向の通信を行う。
通信領域以外の状態配列はGPU上に残る。CUDA-aware MPI、GPUDirect RDMA、非同期重畳は後続課題。

## メモリと性能

`Nv=Ns+4`、所有セル数を `N`、両側2層のghostを含むセル数を `Ng` とすると、
主要GPU配列の概算は `8 * [Ng * (10*Nv+3) + N*Nv]` byte。
これにhaloバッファ、少量の係数・集約領域、CUDAの実行資源が加わる。
MPIでは各rankがこのメモリを必要とする。
CPUには所有領域の入出力配列とhalo、rank 0には出力時の全領域収集メモリも必要。

常駐化は全配列の反復転送を除去するための実装であり、すべての格子サイズでCPUより
速いことを保証しない。小規模格子、多数species、細いペンシル、頻繁なCSV出力は不利。
長時間性能、メモリ削減、通信重畳などの最適化は別途測定して行う。

## 検証

Windows、GNU Fortran、CUDA 13.3、Microsoft MPI、1台のGPUで確認。

- CPU/GPUの時間刻み・流体右辺・4ステップStrang更新を6種類の境界／格子条件で比較。
- Ns=2、4、64と、1セル厚の方向を含む小格子でCPU/GPUの反応流更新を比較。
- 周期／閉壁領域の質量・全エネルギー保存。
- GPUが存在しない指定、負のspecies密度、NASA範囲外、subcycling上限の拒否。
- 初期uploadと明示出力以外に全配列転送が増えないことをカウンタで検査。
- MPI 1/2/3/4 rankの非一様場・割り切れない分割・コーナー・反応更新をCPU逐次版と比較。
  この回帰は全rankがdevice 0を共有する機能検証で、複数GPU実機試験ではない。
- YAMLから環境を生成し、64×16×4の反応ノズルを10ステップ実行。最終出力後に正常終了。
  同じ入力を2 MPI rank（device 0共有）でも10ステップ実行した。

複数GPU／複数ノード実機、Linux実機、長時間性能は未検証。
Compute Sanitizerはこの環境で対象アプリを起動できず、メモリ検査は未完了。
本テストの合格を実在燃焼器・デトネーションの精度検証の代わりにはしない。
