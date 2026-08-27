# NSE SolverLibrary

MPI＋CUDAで分散HIT初期化またはPetersen–Livescu forcingを使う場合は、
[`docs/NSE_CUFFTMP.md`](docs/NSE_CUFFTMP.md)の`cuda_mpi_cufftmp`手順を参照してください。

MPI/OpenMP、単一GPU CUDA、またはMPI＋CUDAマルチGPUで実行する、5保存変数の三次元圧縮性
Navier-Stokesソルバーです。周期境界、SSPRK3、KEEP/WENO系の対流流束、
6次精度粘性項、Taylor-Green/HIT初期条件、Petersen-Livescu Forcingをモジュール化しています。

## 最初に読む文書

- ソルバー全体仕様・無次元化: [`docs/NSE_SOLVER_SPECIFICATION.md`](docs/NSE_SOLVER_SPECIFICATION.md)
- 生成・ビルド・実行: [`docs/NSE_BUILD_AND_RUN.md`](docs/NSE_BUILD_AND_RUN.md)
- 対流ハイブリッド: [`docs/NSE_HYBRID_FLUX.md`](docs/NSE_HYBRID_FLUX.md)
- CUDA対応範囲: [`docs/NSE_CUDA.md`](docs/NSE_CUDA.md)
- モジュール設計と拡張: [`docs/NSE_MODULE_DESIGN.md`](docs/NSE_MODULE_DESIGN.md)
- WENO5-Z/Roe: [`docs/NSE_WENO5Z_ROE.md`](docs/NSE_WENO5Z_ROE.md)
- HITと分散FFT: [`docs/NSE_HIT_DISTRIBUTED_FFT.md`](docs/NSE_HIT_DISTRIBUTED_FFT.md)
- Forcing: [`docs/NSE_FORCING.md`](docs/NSE_FORCING.md)
- 乱流統計: [`docs/NSE_TURBULENCE_STATISTICS.md`](docs/NSE_TURBULENCE_STATISTICS.md)
- 保存済み乱流場の配置: [`docs/NSE_IMPORTED_TURBULENCE.md`](docs/NSE_IMPORTED_TURBULENCE.md)

## 対応プロファイル

| profile | 実行方式 | 主な用途 |
|---|---|---|
| `cpu_mpi` | MPI + OpenMP | Taylor-Green、汎用CPU計算 |
| `cpu_mpi_2decomp_fftw` | MPI + 2DECOMP&FFT | HIT初期化、分散FFT Forcing |
| `cuda_single` | 単一GPU | CUDAによる時間発展、cuFFT HIT/Forcing |
| `cuda_mpi` | MPI + CUDA | 1 rank＝1 GPUのy-z分割時間発展 |
| `cuda_mpi_cufftmp` | MPI + CUDA + cuFFTMp | 分散HIT初期化、分散FFT Forcing |

新規cloneでは`2decomp-fft`サブモジュールも取得します。

```powershell
git clone --recurse-submodules <SolverLibrary URL>
```

## 対流流束

CPU/MPI/OpenMP版、単一GPU CUDA版、MPI＋CUDA版で、次の4方式を実行時に選択できます。

| 設定値 | 内容 |
|---|---|
| `keep2` | 2次精度KEEP |
| `keep6` | 6次精度KEEP（既定） |
| `weno5z_roe` | 特性空間の5次精度WENO-Z再構築 + Roe流束 |
| `hybrid` | Ducros-pressureセンサーでKEEP/WENOを連続混合 |

全領域でWENO5-Z/Roeを使う場合:

```yaml
numerics:
  convective_scheme: weno5z_roe
```

滑らかな領域をKEEP6、衝撃波領域をWENO5-Z/Roeにする場合:

```yaml
numerics:
  convective_scheme: hybrid
  hybrid:
    smooth_scheme: keep6
    shock_scheme: weno5z_roe
    sensor: ducros_pressure
    sensor_onset: 0.01
    sensor_full: 0.10
```

対流流束の選択項目は`numerics.convective_scheme`だけです。
`convective_scheme: keep6`ではWENO計算を行いません。非ハイブリッド方式では
`hybrid`以下の値は計算にも生成後の`input.dat`にも使用されません。
旧`numerics.flux`、旧`numerics.reconstruction`、`keep`や`weno`単独の名前は
入力できません。各方式は3層のghostセルを使います。

## 粘性項と時間積分

`central6`は、一定粘性係数のNewton流体、Stokesの仮定、Fourier熱伝導を
6次精度中心差分で評価します。`none`で粘性項を無効化できます。
どちらもCPU版、単一GPU版、MPI＋CUDA版に対応します。

```yaml
physics:
  nse:
    reynolds_number: 100.0
    prandtl_number: 0.72

numerics:
  viscous_scheme: central6
  boundary_condition: periodic
  time_integrator: ssprk3
```

自動時間刻みは対流CFL条件と拡散安定条件を併用します。

## 最短実行例

通常は`ScriptLibrary/RunEnvironment`から生成した実行環境で次を実行します。

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

最終出力後に次が表示されれば正常終了です。

```text
NSE calculation completed successfully: step=<final step>, time=<final time>
```

SLFをParaView用に変換する場合:

```powershell
python .\tools\postprocess_case.py
```

## ソース構成

| ディレクトリ | 内容 |
|---|---|
| `src/common` | 精度、定数、共通設定、モデル設定 |
| `src/grid`, `src/field` | FVM格子と保存変数 |
| `src/init`, `src/forcing` | 初期条件とForcing |
| `src/boundary` | 境界条件 |
| `src/numerics` | 対流、再構築、Riemann solver、粘性、空間演算 |
| `src/time` | 時間積分 |
| `src/gpu` | CUDAバックエンド |
| `src/io`, `src/mpi`, `src/main` | 入出力、MPI、実行制御 |
| `config` | NSE単体YAMLビルド設計 |
| `tests` | CPU、MPI、CUDA回帰テスト |

CMake、`solver_manifest.yaml`、`config/module_catalog.yaml`が同じ現行ソース構成を参照します。
