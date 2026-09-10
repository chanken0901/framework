# 単成分NSE：CUDA-aware MPIによるGPUバッファ直接通信

単成分NSEの`cuda_mpi`／`cuda_mpi_cufftmp`に、GPU上のhaloを直接MPIへ渡す経路を追加した。
LLNS揺らぎ拡張にも共通で適用する。多成分・反応流とGPEは今回の変更対象外。
CPU版、単一GPU版、物理条件やMPIのy-z分割方法は変更しない。

## 通信仕様

- GPUでpack → GPUバッファを`MPI_Sendrecv`へ渡す → GPUでunpack。
- y交換後にz交換し、辺・角のghostを引き継ぐ。非周期端の`MPI_PROC_NULL`では受信領域を上書きしない。
- pack完了をCUDA同期で保証してからMPIを呼ぶ。受信はblockingで完了を待ち、unpack完了後にバッファを再使用する。
- 通信バッファはGPU上で再使用する。直接通信時にはFortran側のホストhalo配列を確保しない。
- 従来のGPU送信作業領域に加え、最大方向halo要素数×8 byteの受信領域がGPUに必要。
- MPIと計算の非同期重畳、stream-aware MPIは今回未実装。CFL・異常判定等のスカラーMPI通信、出力時のCPU転送は残る。

アプリケーションがホストへ中継しない仕様であり、実際のNVLink／PCIe P2P／GPUDirect RDMAの
利用はMPI・UCX・GPU・NIC・ドライバ・配置に依存する。`device`ログだけではRDMA使用の証明にならない。
[Open MPI CUDAガイド](https://docs.open-mpi.org/en/v5.0.x/tuning-apps/networking/cuda.html)も参照。

## 必要環境と選択

現実装の能力確認はOpen MPIの`MPIX_Query_cuda_support()`を利用する。
CUDA対応Open MPI、または同APIを提供するOpen MPIベースの配布環境を使用する。
MPI初期化後に全rankの能力と設定を照合する。検出できないMPIへGPUポインタを試しに送ることはしない。
[能力確認APIの仕様](https://docs.open-mpi.org/en/main/man-openmpi/man3/MPIX_Query_cuda_support.3.html)。

ビルド時：`NSE_ENABLE_CUDA_AWARE_MPI=ON`（既定OFF）。`NSE_USE_MPI=ON`、
`NSE_GPU_BACKEND=cuda`、`NSE_MPI_PROVIDER=SYSTEM`が必要。
CとFortranで**同じMPI実装・バージョン**をリンクすること。
MS-MPI用ビルドでONを指定すると構成時にエラーにする。

実行時の環境変数`NSE_CUDA_MPI_TRANSPORT`は全rankで同じ値を指定する。

| 値 | 動作 |
| --- | --- |
| `auto`（既定） | 全rankでビルド・実行時対応を確認できればdevice、確認できなければstaged |
| `device` | 直接通信を必須にする。非対応・未検出・ビルド無効なら通信開始前にMPI全体を異常終了 |
| `staged` | 従来のCPUバッファ経由を使用。比較・互換運用向け |

起動時にroot rankが`NSE CUDA MPI halo transport: device ...`または`staged ...`を表示する。
直接通信が必要な本番実行・検証では`device`を明示し、意図しないフォールバックを防ぐ。
1rank／GPUを基本とし、GPU選択は従来の`NSE_CUDA_DEVICE_POLICY`とジョブ配置を使用する。

## Linuxでのビルド・実行

生成環境では`parallel.use_mpi: true`、`parallel.use_cuda: true`、`parallel.use_openmp: false`で生成し、
生成された`ScriptLibrary/BuildSolver/build.local.yaml`の既存`models.nse`へ以下を追加する。
同階層の`profile`等を消さず、既存`cmake_overrides`があればそこへ追記する。

```yaml
models:
  nse:
    cmake_overrides:
      NSE_ENABLE_CUDA_AWARE_MPI: true
```

Linux用machine profileとCUDA対応MPIを用意したうえで、生成環境のルートで実行する。

```bash
python3 ./tools/run_case.py --build
export NSE_CUDA_MPI_TRANSPORT=device
python3 ./tools/run_case.py --run
```

FrameWork直下から直接ビルドする場合（CUDA Toolkitと互換ホストコンパイラがPATH上にあること）：

```bash
cmake -S SolverLibrary/NSE -B build/nse-device-mpi \
  -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc \
  -DNSE_MPI_PROVIDER=SYSTEM -DNSE_USE_MPI=ON -DNSE_ENABLE_OPENMP=OFF \
  -DNSE_GPU_BACKEND=cuda -DNSE_ENABLE_CUDA_AWARE_MPI=ON \
  -DNSE_INIT_FFT_BACKEND=none -DNSE_FORCING_FFT_BACKEND=none \
  -DNSE_VISCOUS_SCHEME=central6 -DNSE_BOUNDARY_SCHEME=runtime \
  -DBUILD_TESTING=ON -DNSE_ENABLE_MULTI_GPU_TESTS=ON
cmake --build build/nse-device-mpi -j 8
export NSE_CUDA_MPI_TRANSPORT=device
# 2/4rankのテストに必要なGPUをジョブスケジューラで確保して実行する。
ctest --test-dir build/nse-device-mpi -R 'nse_cuda_fh_mpi|nse_mpi_cuda_' --output-on-failure
```

cuFFTMpを使うケースは既存`cuda_mpi_cufftmp`のFFT依存設定を維持してONを追加する。
この選択はNSEのhalo通信のみを変更し、cuFFTMp内部の通信方法を変更しない。

## Windows / PowerShell

現在のMS-MPI環境では直接通信を使わず、ビルドオプションはOFFのままにする。

```powershell
$env:NSE_CUDA_MPI_TRANSPORT = 'staged'
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

Windows側からLinux向け実行環境を生成・転送する場合は、Linux側のビルドでONにする。
既存ResearchRunsのソースや実行ファイルは自動更新されない。

## 検証状況と残作業

開発機では直接通信部分のCUDA/MPI APIコンパイル、従来経路の回帰、未対応時の`device`拒否と
不正設定拒否を検証した。MS-MPI＋1GPUのためCUDA-aware Open MPI経路の実通信、
複数GPU・複数ノードでの一致と性能は未検証。
対応実機では上記試験に加え、同一入力・seed・dtで`staged`と`device`の結果を比較する。
LLNSの2/4rank試験はCUDA-awareビルドで`device`を強制し、単一領域CUDAとの一致と保存性を確認する。
Nsight SystemsやMPI/UCXの診断で通信経路を確認した後に性能比較を行うこと。
