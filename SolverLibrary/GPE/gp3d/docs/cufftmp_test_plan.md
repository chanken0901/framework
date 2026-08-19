# cuFFTMpスパコン検証手順

この手順は、最初の移植結果を再現可能な形で採取し、次の修正と性能最適化に使うためのものです。
最初は同一ノード内の2 GPUで正しさを確認し、その後に1、2、4 GPUの比較を行います。

## 1. 実行環境の記録

サイトのモジュール名に合わせてHPC SDKとMPIをロードし、次の出力を保存してください。

```bash
module list 2>&1
which mpifort mpicxx nvcc cmake
mpifort --version
mpicxx --version
nvcc --version
mpirun --version
nvidia-smi
nvidia-smi topo -m
echo "CUFFTMP_HOME=$CUFFTMP_HOME"
echo "NVSHMEM_HOME=$NVSHMEM_HOME"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
```

cuFFTMpでは、同じノード内の全rankから同一の`CUDA_VISIBLE_DEVICES`が見える必要があります。
rankごとにGPUを一つだけ見せる設定にはしません。コードがノード内rank番号からGPU番号を選びます。
NVSHMEMの通信用スレッドがあるため、各rankには最低2 CPU coreを割り当てます。

HPC SDK付属の`nvshmem_bootstrap_mpi.so`は、同梱されるHPC-X向けです。
まずHPC SDK/HPC-Xの`mpifort`と`mpicxx`を同じビルドで使用してください。
サイト標準MPIを使う場合は、そのMPIに対して構築されたNVSHMEM bootstrap pluginが必要になることがあります。
`nvshmem_bootstrap_mpi.so`をロードできない場合は、NVSHMEMの`lib`を`LD_LIBRARY_PATH`へ追加し、
MPI ABIが異なる場合はサイト管理者に対応pluginの有無を確認します。

## 2. ビルド

```bash
cmake -S . -B build/cufftmp -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_MPI=ON \
  -DFFT_BACKEND=dft \
  -DGPU_BACKEND=cufftmp \
  -DCUFFTMP_ROOT="$CUFFTMP_HOME" \
  -DNVSHMEM_ROOT="$NVSHMEM_HOME" \
  -DBUILD_TESTING=ON
cmake --build build/cufftmp --parallel 8
```

HPC SDK 25.3より前で`cufftMpMakePlan3d`が見つからない場合は、
`-DCUFFTMP_API=legacy`を追加して再構成します。

## 3. 2 GPU正しさ確認

```bash
export NVSHMEM_SYMMETRIC_SIZE=8G
ctest --test-dir build/cufftmp -V -R cufftmp
```

次の二つが成功すれば、実時間・虚時間split-step、診断量、ARGLEの最小経路を確認できています。

- `cufftmp_splitstep_np2`
- `cufftmp_argle_smoke_np2`

CTestがサイトのGPU割り当てと合わない場合は、確保した計算ノード上で直接実行します。

```bash
mpirun -np 2 build/cufftmp/gp3d_test_cufftmp_splitstep
mpirun -np 2 build/cufftmp/gp3d_cufftmp tests/cuda_smoke.nml
```

## 4. 性能比較

同じノード、同じGPU型、同じ入力でrank数だけを変えます。

```bash
for np in 1 2 4; do
  mpirun -np "$np" build/cufftmp/gp3d_cufftmp tests/cufftmp_scaling_256.nml \
    | tee "cufftmp_np${np}.log"
done
```

`split-step timing`表の次の行を比較します。

- `nonlinear_local`
- `fft_forward_inverse`
- `kinetic_spectral`
- `allocation_and_other`
- `step_total`

初回は`256^3`で通信を含む基本動作を確認します。これで複数GPUが不利でも直ちに異常とは限りません。
その後、利用可能メモリに応じて`512^3`以上へ上げると分散化の効果を判断しやすくなります。

## 5. Slurm例

サイト固有のMPI起動方法が優先です。次は1ノード2 GPUの最小例です。

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00

module purge
# module load <site MPI and NVIDIA HPC SDK modules>

export OMP_NUM_THREADS=1
export NVSHMEM_SYMMETRIC_SIZE=8G
export CUDA_DEVICE_ORDER=PCI_BUS_ID

srun --gpu-bind=none build/cufftmp/gp3d_test_cufftmp_splitstep
srun --gpu-bind=none build/cufftmp/gp3d_cufftmp tests/cufftmp_scaling_256.nml
```

## 6. フィードバックに含めるもの

次の内容があれば、API差異、NVSHMEM設定、数値誤差、通信時間のどこを直すべきか判断できます。

1. 環境記録の全出力
2. CMake configureとbuildの末尾
3. `ctest -V -R cufftmp`の全出力
4. `cufftmp_np1.log`、`cufftmp_np2.log`、可能なら`cufftmp_np4.log`
5. ジョブスクリプトと割り当てたnode、rank、GPU、CPU core数
6. エラー時は最初のエラーだけでなく、その前後を含む標準出力と標準エラー
