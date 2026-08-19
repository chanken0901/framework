# GP3Dソルバーパッケージ

このディレクトリは、SolverLibraryで管理するGP3Dパッケージの原本です。
数値計算モジュールと、薄い実行プログラムから構成されています。
ケースごとに生成される実行環境には、選択したプロファイルで必要なファイルだけがコピーされます。

## コード全体仕様書

物理モデル、数値解法、モジュールAPI、CPU/MPI/CUDA/cuFFTMpバックエンド、
入力、SLF出力、再スタート、ビルド、テストおよび拡張規約の詳細は、
[GPE/GP3Dコード全体仕様書（Markdown）](../../../GPE_GP3Dコード全体仕様書.md)を参照してください。
配布・閲覧用の[Word版](../../../GPE_GP3Dコード全体仕様書.docx)も同じ場所にあります。

## ディレクトリ構成

| パス | 役割 |
| --- | --- |
| `src/common` | 共通の数値型 |
| `src/grid` | 周期境界を持つ直交格子 |
| `src/init` | Taylor-Green/ARGLEを含むGPE初期条件 |
| `src/fft` | 逐次および分散CPU FFTバックエンド |
| `src/gpu` | CUDA/cuFFTおよびMPI/cuFFTMpブリッジとFortranバインディング |
| `src/solver` | Split-operator法による時間積分 |
| `src/io` | 入力、SLF出力、再スタート読み込み |
| `src/mpi` | 実MPIモジュールと逐次実行用スタブ |
| `src/main` | CPU/MPI、単一GPU、複数GPU用エントリーポイント |
| `tools` | ビルド・実行ワークフローとParaView変換ツール |
| `tests` | バックエンドと再スタート機能のテスト |

`solver_manifest.yaml`が、このパッケージを組み立てるための公開された定義です。
各プロファイルには、コンポーネント、CMake定義、実行ファイル、追加テストが記述されます。
フレームワークの生成スクリプトは、この定義から必要ファイルを解決し、
Fortranモジュールの依存関係を検査してからコピーします。

NSEとGPEを同じ入口からビルド・実行する場合は、フレームワーク側の
`ScriptLibrary/BuildSolver/build_model.py`を使用します。共通`build.yaml`の
`selected_model: gpe`と、本マニフェストのプロファイル名を指定します。

ケースに依存するパラメータを、このパッケージ内へ追加しないでください。
初期条件の種類とそのパラメータは`cases/<case_id>/case.yaml`へ記述し、
フレームワークの入力アダプターを通してGP3Dのnamelistへ変換します。

## 対応している構成

- 検証用DFTまたはFFTWを使用するCPU逐次実行
- z方向スラブ分割と、局所DFTまたはFFTWを使用するMPI分散実行
- MPI rank内のCPU処理をOpenMPで並列化するMPI+OpenMPハイブリッド実行
- CUDAとcuFFTを使用するNVIDIA GPU 1台での実行
- 1 MPI rankにつき1 GPUを使用する、CUDAとcuFFTMpによる分散実行

複数GPU版は`cuda_mpi_cufftmp`プロファイルとして独立しているため、
既存のCPU版および単一GPU版で使用しているケースYAMLの形式を維持できます。

## MPI+OpenMPハイブリッド実行

`cpu_mpi_dft`と`cpu_mpi_fftw`はOpenMP対応でビルドされます。同じ実行ファイルのまま、
ケースごとに`solver.use_openmp`を切り替えられるため、ON/OFFのたびに再ビルドする必要はありません。

```yaml
solver:
  profile: cpu_mpi_fftw
  processes: 4
  use_openmp: true
  omp_threads: 4
```

`case.yaml`を編集した後、入力を再生成して実行します。

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --run
```

起動時の`# OpenMP compiled=... active=... threads_per_rank=...`で実際の設定を確認できます。
OpenMPは各rank内の局所項、運動エネルギー係数、ARGLE、診断量、初期条件、および
分散FFTを構成する独立な1次元変換に使われます。MPI転置通信はOpenMP並列領域の外で実行されます。

## cuFFTMp版のビルドと実行

cuFFTMp版はLinux専用です。NVIDIA HPC SDK、対応するNVSHMEM、MPI、CMake、
CUDA対応C++コンパイラ、Fortranコンパイラが必要です。計算ノードでは
各rankから同じ`CUDA_VISIBLE_DEVICES`が見えるようにし、1 rankを1 GPUへ割り当てます。

HPC SDKの環境をロードした後、通常は次のように構成します。

```bash
cmake -S . -B build/cufftmp -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_MPI=ON \
  -DFFT_BACKEND=dft \
  -DGPU_BACKEND=cufftmp \
  -DCUFFTMP_ROOT="$CUFFTMP_HOME" \
  -DNVSHMEM_ROOT="$NVSHMEM_HOME"
cmake --build build/cufftmp --parallel 8
```

HPC SDK 25.3以降のAPIは自動検出されます。古いHPC SDKで自動検出が合わない場合は
`-DCUFFTMP_API=legacy`、新しいAPIを明示する場合は`-DCUFFTMP_API=modern`を追加します。

2 GPUでの比較テストと実行例は次のとおりです。

```bash
ctest --test-dir build/cufftmp --output-on-failure -R cufftmp
export NVSHMEM_SYMMETRIC_SIZE=8G
mpirun -np 2 build/cufftmp/gp3d_cufftmp tests/cuda_smoke.nml
```

JSON設計書からまとめて実行する場合は、サンプルを計算条件に合わせて編集してから
Linux用ランナーへ渡します。

```bash
python3 tools/run_workflow_linux.py workflow.cufftmp.example.json
```

実空間の波動関数はzスラブで保持され、順FFT後はcuFFTMpの組み込みyスラブ配置を
運動項カーネルが直接処理します。SLFは従来どおりrank別zスラブとして出力されるため、
既存の再スタート機能とParaView変換ツールを利用できます。

SLFからVTI/PVDへの変換と、ParaViewでの密度・位相・量子渦表示の標準手順は、
[`後処理_ParaView可視化手順書.md`](../../../後処理_ParaView可視化手順書.md)を参照してください。

スパコンでの環境記録、2 GPU検証、1/2/4 GPU比較、Slurm例は
`docs/cufftmp_test_plan.md`にまとめています。
