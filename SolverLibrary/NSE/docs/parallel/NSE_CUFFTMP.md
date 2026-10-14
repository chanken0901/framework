# NSE MPI＋CUDA分散FFT（cuFFTMp）

## 対応機能

`cuda_mpi_cufftmp`プロファイルは、既存のMPI＋CUDA時間発展を維持したまま、
次の2機能をcuFFTMpで実行します。

- `flow.type: hit_spectral`の分散HIT初期化
- `forcing.type: petersen_livescu`の分散Helmholtz分解

NSEの領域分割は、各rankがx方向全域とy-z方向の局所ブロックを保持する方式です。
cuFFTMpには`cufftMpMakePlanDecomposition`で、このY-Z分割を入力・出力の両方に
指定します。root rankへの場の集約は行いません。

HIT初期化では、既存の2DECOMP&FFT版と同じ乱数生成、低波数等方化、RMS調整、
圧力Poisson計算を使用します。FFT時だけ局所ホスト配列とcuFFTMp記述子の間を
転送します。初期化終了後、その一時的なFFT planと記述子は解放されます。

Petersen–Livescu forcingでは、保存変数、密度重み付き速度、スペクトル、RHSを
GPU上に保持します。3速度成分を分散FFTし、局所スペクトル上でHelmholtz分解を
行います。分母と圧力膨張相関のスカラー値だけを`MPI_Allreduce`し、逆FFT後に
各rankの運動量RHSへ加算します。

## 必要環境

- Linux
- NVIDIA GPU（MPI rankごとに1台）
- NVIDIA HPC SDKまたは互換性のあるCUDA C++ツールチェーン
- 同じHPC SDKリリースに含まれるcuFFTMpとNVSHMEM
- NVSHMEM bootstrapと互換性のあるMPI
- cuFFTMp 11.4.0以降

Windowsでは既存の`cuda_mpi`プロファイルを利用できますが、cuFFTMp版の構成は
CMakeが明示的に拒否します。

## ビルド

HPC SDKに合わせて次の環境変数を設定します。cuFFTMp、NVSHMEM、MPIを異なる
HPC SDKリリースから混在させないでください。

```bash
export CUFFTMP_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/math_libs
export NVSHMEM_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/comm_libs/nvshmem
export MPI_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/comm_libs/hpcx/latest/ompi

cmake -S SolverLibrary/NSE -B build/nse-cufftmp \
  -DCMAKE_BUILD_TYPE=Release \
  -DNSE_USE_MPI=ON \
  -DNSE_ENABLE_OPENMP=OFF \
  -DNSE_GPU_BACKEND=cuda \
  -DNSE_INIT_FFT_BACKEND=cufftmp \
  -DNSE_FORCING_FFT_BACKEND=cufftmp
cmake --build build/nse-cufftmp -j
```

環境によって標準配置が異なる場合は、従来の変数`CUFFT_INC`、`CUFFT_LIB`、
`NVSHMEM_INC`、`NVSHMEM_LIB`も使用できます。

## case.yaml

HITとforcingを同時に使う例です。

```yaml
solver:
  profile: cuda_mpi_cufftmp
  mpi_processes: 4

flow:
  type: hit_spectral

forcing:
  type: petersen_livescu
  petersen_livescu:
    spectrum: low_wavenumber
    fft_backend: cufftmp
    k_cutoff: 2.5
    target_dissipation: 0.1
    dilatational_ratio: 0.0
    denominator_floor: 1.0e-14
    max_coefficient: 0.0
    report_interval: 100
```

HITだけを使う場合も`cuda_mpi_cufftmp`を選び、forcingを`none`にします。
forcingだけを使う場合は、Taylor–Greenまたは読み込み初期条件と組み合わせられます。

## 実行

スケジューラがrankごとにGPUを1台だけ公開する場合、各rankからはdevice 0として
見えるため、そのまま実行できます。

```bash
mpirun -n 4 build/nse-cufftmp/bin/nse_mpi_cuda cases/case0001/input.dat
```

cuFFTMp/NVSHMEMの初期化APIは全MPI rankが同じ順序で呼びます。planは
`MPI_Finalize`より前に破棄されます。実行時のbootstrapエラーは、MPIとNVSHMEMの
組合せ、`NVSHMEM_BOOTSTRAP`、ジョブランチャー、GPU割当を確認してください。

## 検証範囲

Windows開発機では、従来のMPI＋CUDA時間発展、単一GPU、CPU/MPIの回帰テストと、
cuFFTMp用Fortran/C++インターフェースの契約検証を行います。実際の分散FFTの数値・
性能検証は、cuFFTMpと複数GPUを備えたLinux計算機上で実行してください。
