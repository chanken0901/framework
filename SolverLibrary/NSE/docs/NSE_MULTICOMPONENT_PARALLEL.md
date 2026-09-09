# Stage 8: 多成分・反応流のMPI/OpenMP

更新日: 2026-09-07

Stage 7の物理モデルと境界条件を保ち、OpenMPおよびMPIペンシル分割を追加した。
既存の逐次profileは引き続き使用できる。

| 実行方式 | profile | use_mpi | use_openmp |
|---|---|---|---|
| 逐次 | cpu_serial_reactive_boundaries | false | false |
| OpenMP | cpu_openmp_reactive | false | true |
| MPI | cpu_mpi_reactive_pencil | true | false |
| MPI＋OpenMP | cpu_mpi_reactive_pencil | true | true |

MPI profileはOpenMP対応でビルドする。実行時の`use_openmp: false`では1スレッドを使用する。
この場合は`omp_threads: 1`も指定する。
`use_mpi`の仕様は維持し、選択profileと矛盾する指定は入力生成時に拒否する。
生成済み実行環境のprofileを変更する際は、環境を再生成して再ビルドする。

## ペンシル分割

x方向は分割せず、y・z方向を`Py × Pz`のCartesianプロセス格子で分割する
（x-pencil）。FFTライブラリは使用しない。保存変数数は`Ns+4`で可変である。
分割数で割り切れない格子は各rankへ余りを配分する。

対流・輸送RHSの各SSPRK段で2層のhaloを交換する。y方向の交換後にyのhaloを含む
z方向の交換を行うため、横方向の勾配に必要なコーナーも正しい隣接rankから取得する。
物理境界は大域領域の外面に適用する。無反射境界の自動長さ尺度は大域領域の長さを使う。
化学反応は各rankの所有セルだけで進める。

時間刻みは全rankの最小値、化学substep数は最大値、保存量は総和で集約する。
OpenMPではセルごとの熱力学・化学反応・時間刻み・勾配を並列化する。
面流束は独立した格子列を各スレッドへ割り当て、隣接セル更新の競合を防ぐ。
MPIはOpenMP並列領域の外で呼び、`MPI_THREAD_FUNNELED`を要求する。

## case.yaml

例: 4 MPIプロセス × 2 OpenMPスレッド。

```yaml
solver:
  profile: cpu_mpi_reactive_pencil
  use_mpi: true
  use_openmp: true
  use_cuda: false
  mpi_processes: 4
  omp_threads: 2
  decomposition: pencil
  process_grid: [2, 2]   # [Py, Pz]; product = mpi_processes
```

`process_grid`を省略するか`[0, 0]`にすると`MPI_Dims_create`で自動決定する。
`[2, 0]`のように一方を固定することもできる。分割する方向は各rankが最低2セルを
所有する必要がある。例えば24 MPIなら`[6,4]`に対し`ny >= 12, nz >= 8`が必要である。
偏った格子で自動分割がこの条件を満たさない場合は、格子に合う分割数を明示する。
`Py=1`または`Pz=1`も同じペンシル実装で扱い、1 MPIプロセスでも実行できる。
スラブを選択する設定は提供しない。

通常の物理・時間・境界設定および`config/*.yaml`の4拡張はStage 7と同じ。
`input.dat`には必要なMPI profileの場合だけ次のグループを生成する。

```fortran
&multicomponent_parallel
  decomposition = "pencil"
  process_grid = 2, 2
/
```

## ビルド・検証

以下はFrameWorkルートで実行する。既存のテスト入力を使って最初に動作確認できる。

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\BuildSolver\build_model.py --model nse_multicomponent --profile cpu_mpi_reactive_pencil --build-dir .\build\stage8 --omp-threads 2 --test
python .\ScriptLibrary\BuildSolver\build_model.py --model nse_multicomponent --profile cpu_mpi_reactive_pencil --build-dir .\build\stage8 --processes 4 --omp-threads 2 --input-file .\SolverLibrary\NSE\tests\input_multicomponent_reactive_boundaries.dat --run
```

Linux（bash、GNU Fortran＋システムMPI）:

```bash
python3 ./ScriptLibrary/BuildSolver/build_model.py --model nse_multicomponent --profile cpu_mpi_reactive_pencil --machine-profile ./ScriptLibrary/BuildSolver/machine_profiles/linux_gnu_mpi.yaml --build-dir ./build/stage8 --omp-threads 2 --test
python3 ./ScriptLibrary/BuildSolver/build_model.py --model nse_multicomponent --profile cpu_mpi_reactive_pencil --machine-profile ./ScriptLibrary/BuildSolver/machine_profiles/linux_gnu_mpi.yaml --build-dir ./build/stage8 --processes 4 --omp-threads 2 --input-file ./SolverLibrary/NSE/tests/input_multicomponent_reactive_boundaries.dat --run
```

OpenMP単独ではprofileを`cpu_openmp_reactive`、ビルド先を別名、プロセス数を1にする。
MPIのみの比較では`--omp-threads 1`を使用する。
直接実行時は`OMP_NUM_THREADS`を設定してから`mpiexec -n 4 <実行ファイル> <input.dat>`
を実行する。ハイブリッドでは`MPIプロセス数 × OMP_NUM_THREADS`が割当CPU数を
超えないようにする。

## 実行環境生成

MPIには`ScriptLibrary/RunEnvironment/environment.nse_multicomponent.parallel.yaml`、
OpenMP単独には`environment.nse_multicomponent.openmp.yaml`を使用する。
生成時の既定はMPI 4プロセス、OpenMP 1スレッドである。
生成された`case.yaml`の`solver.omp_threads`などを上記の例に合わせて編集する。

Windows（PowerShell）:

```powershell
python .\ScriptLibrary\RunEnvironment\prepare_environment.py .\ScriptLibrary\RunEnvironment\environment.nse_multicomponent.parallel.yaml
```

Linuxでは同設計書の`select.target`を`linux_gnu_mpi`、
`select.destination`を`local_generated`に変更し、必要なら
`destination.root`と`destination.case_index`を希望するLinuxパスに指定する。

```bash
python3 ./ScriptLibrary/RunEnvironment/prepare_environment.py ./ScriptLibrary/RunEnvironment/environment.nse_multicomponent.parallel.yaml
```

生成先に移動して実行する。

| 処理 | PowerShell | bash |
|---|---|---|
| 入力再生成 | `python .\tools\run_case.py --prepare` | `python3 ./tools/run_case.py --prepare` |
| ビルド | `python .\tools\run_case.py --build` | `python3 ./tools/run_case.py --build` |
| 実行 | `python .\tools\run_case.py --run` | `python3 ./tools/run_case.py --run` |

## 出力・制約・検証範囲

CSV、時系列snapshot、履歴と終了メッセージはrank 0が書く。場の出力時だけ
全領域をrank 0へ集めるため、CSV形式と座標・変数順は逐次版と同じになる。
rank 0には全領域の場と収集バッファのメモリが必要で、大規模出力は今後の改善対象である。
MPI版の終了時に、出力を含む最遅rankの経過秒を表示する。

領域分割の方向は現段階でy・z固定。一般座標はStage 9、
GPU常駐CUDAは [Stage 10](NSE_MULTICOMPONENT_CUDA.md) に追加した。MPI-IO、
通信と計算の非同期重畳、化学反応量に応じたMPI負荷分散は未実装である。
細いペンシルではhalo計算の割合が増えるため、小規模格子で速度向上は保証されない。

回帰試験は1・2・3・4 MPI、OpenMP併用、割り切れない格子、全方向の物理境界、
y/z方向に非一様な場、コーナー通信、二段のStrang更新を逐次計算と比較する。
テストの速度比だけで大規模・多ノード性能を評価しない。
