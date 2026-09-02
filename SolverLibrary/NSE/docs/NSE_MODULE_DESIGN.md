# NSEモジュール設計

> テスト・生成コマンドのWindows（PowerShell）／Linux（bash）対応は[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

実装済みのKEEP/WENOハイブリッド流束と、ハイブリッド／非ハイブリッドの
入力設定は[`NSE_HYBRID_FLUX.md`](NSE_HYBRID_FLUX.md)に定義する。
面別周期／無反射／鏡像境界の仕様は
[`NSE_BOUNDARY_CONDITIONS.md`](NSE_BOUNDARY_CONDITIONS.md)に定義する。

## 目的

NSEソルバーの初期条件、移流スキーム、粘性項、境界条件を、mainプログラムを変更せずに交換・追加できる構造にする。
現在は5変数の圧縮性Navier-Stokes方程式を扱い、2次・6次精度KEEP対流項、
5次精度WENO-Z/Roe対流項、KEEP/WENOハイブリッド対流項、非粘性計算または
6次精度中心差分による粘性計算を選択できる。

## 依存方向

依存関係は次の一方向とする。

```text
main
  -> time
     -> numerics/spatial_operator
        -> boundary
        -> numerics/convective
        -> numerics/viscous
  -> init
  -> io

各機能モジュール
  -> field/grid/config/mpi
  -> common
```

`common`、`grid`、`field`、`mpi`から上位の`main`や個別スキームを参照してはならない。
mainは計算手順だけを制御し、流束式、初期条件式、境界処理を持たない。

## 現在の構成

| 分類 | ファイル | 公開APIまたは役割 |
|---|---|---|
| 実行制御 | `src/main/main_nse.f90` | 入力、MPI初期化、時間ループ、出力 |
| 単一GPU実行制御 | `src/main/main_nse_cuda.f90` | 単一GPUの時間ループと出力 |
| MPI＋CUDA実行制御 | `src/main/main_nse_mpi_cuda.f90` | MPI分割、GPU割当て、段別SSPRK3、全rank CFL同期 |
| CUDA API | `src/gpu/mod_nse_gpu_cuda.f90` | FortranからCUDA context、時間積分、halo pack/unpackを操作 |
| MPI＋CUDA halo | `src/gpu/mod_nse_gpu_mpi.f90` | host stagingによるy→z面交換と辺・角ghost伝播 |
| 初期条件選択 | `src/init/mod_nse_initial_conditions.f90` | `initialize_nse_state` |
| Taylor–Green | `src/init/mod_init_taylor_green.f90` | `initialize_taylor_green` |
| 分散FFT HIT | `src/init/mod_init_hit_spectral_2decomp.f90` | `initialize_hit_spectral` |
| 保存乱流場 | `src/init/mod_init_imported_turbulence.f90` | 可搬SLFの`embed`/`tile`配置 |
| 平面衝撃波–乱流 | `src/init/mod_init_shock_turbulence.f90` | 形成済み衝撃波とDirichletドライバーの配置 |
| 有限高圧室–乱流 | `src/init/mod_init_shock_tube_turbulence.f90` | 隔膜、高圧室、低圧域、局所乱流の配置 |
| CPU面別境界 | `src/boundary/mod_boundary_runtime.f90` | 周期／特性無反射／鏡像／Dirichletの実行時dispatch |
| CUDA面別境界 | `src/gpu/mod_nse_gpu_cuda.f90`、`nse_cuda_bridge.cu` | GPU上の周期／特性無反射／鏡像／Dirichlet ghost処理 |
| 対流項ディスパッチ | `src/numerics/convective/mod_convective_dispatch.f90` | 単一方式とハイブリッドの実行時選択 |
| KEEP流束 | `src/numerics/convective/mod_convective_keep.f90` | 2次・6次精度KEEP |
| WENO-Z再構築 | `src/numerics/reconstruction/mod_reconstruction_weno5z.f90` | 左右5次精度再構築 |
| Roe流束 | `src/numerics/riemann/mod_riemann_roe.f90` | 固有分解とentropy fix |
| WENO-Z/Roe複合流束 | `src/numerics/convective/mod_convective_weno5z_roe.f90` | 特性再構築と三方向面流束 |
| 構成流束レジストリ | `src/numerics/convective/mod_convective_leaf_registry.f90` | 方式名から保存形面流束を選択 |
| ハイブリッド流束 | `src/numerics/convective/mod_convective_hybrid.f90` | センサー、混合率、共有面流束の合成 |
| 粘性項なし | `src/numerics/viscous/mod_viscous_none.f90` | 非粘性計算用の空実装 |
| 6次精度粘性項 | `src/numerics/viscous/mod_viscous_central6.f90` | `add_viscous_rhs`、粘性時間刻み評価 |
| 空間演算 | `src/numerics/mod_nse_spatial_operator.f90` | 境界、流束、発散、粘性項を合成 |
| SSPRK3 | `src/time/mod_nse_time_integration.f90` | CFL時間刻みと3段Runge–Kutta |

## 交換可能モジュールの規約

移流、粘性、境界の各分類は、上位層へ同じFortranモジュール名と公開APIを提供する。
対流項はCPU版で複数の流束モジュールを同時にコンパイルし、ディスパッチ層が
`case.yaml`の指定から実行時に選択する。この方式により、空間演算とmainは
個別スキーム名を知らなくてよい。

### 移流スキーム

Fortranモジュール名は`mod_convective_scheme`とし、次を公開する。

```fortran
compute_convective_flux
validate_convective_scheme
convective_required_ghost_cells
```

現在のKEEP実装は、`convective_scheme='keep2'`では隣接点対と係数`[1/2]`を、
`convective_scheme='keep6'`では距離1～3の点対と係数`[3/4,-3/20,1/60]`を使い、
対称な二点KEEP流束を保存形の面流束へ合成する。
CPU/MPI/OpenMP版とCUDA版は同じ離散式を使用し、3層のゴーストセルを必要とする。

`convective_scheme='weno5z_roe'`では、Roe平均で得た固有ベクトルを各面の
5点ステンシルへ適用し、特性空間で左右状態をWENO-Z再構築する。再構築状態から
Harten-Hyman型entropy fix付きRoe流束を計算する。CPU逐次、MPI、OpenMP、
MPI+OpenMP、単一GPU CUDA版、MPI＋CUDA版に対応し、3層のゴーストセルを必要とする。

### 粘性項

Fortranモジュール名は`mod_viscous_scheme`とし、次を公開する。

```fortran
add_viscous_rhs
validate_viscous_scheme
viscous_required_ghost_cells
viscous_scheme_name
viscous_dt_limit
```

### 境界条件

現行の周期境界backendと、将来の面別runtime dispatcherは、上位層に対して
Fortranモジュール名`mod_nse_boundary`と次のAPIを維持する。

```fortran
apply_nse_boundary
validate_boundary_scheme
boundary_required_ghost_cells
boundary_scheme_name
```

面別境界の実装後、`apply_nse_boundary`は6面の設定を実行時に選択し、一回の呼出しで
MPI内部halo、物理面、辺、角を完成させる。`boundary_scheme_name`は個々の面種別ではなく
コンパイル済みdispatcher backend名`runtime`を返す。`validate_boundary_scheme`は
6面の組合せ、周期面の対、基準状態、実行profileとの互換性を検証する。

## 対流スキームを追加する手順

1. 再構築、Riemann solver、複合流束を責務ごとのモジュールへ分ける。
2. 複合流束は、一つの共有面に一つの保存形面流束を返すAPIを公開する。
3. ハイブリッドの構成流束に使う場合は`mod_convective_leaf_registry.f90`へ方式名と呼び出し先を登録する。
4. 独立した主方式として入力選択させる場合は`mod_convective_dispatch.f90`の検証と分岐も追加する。
5. CMake、`solver_manifest.yaml`、`config/module_catalog.yaml`へ依存順を登録する。
6. `case_input.py`へ許可する入力名と実行バックエンド制約を登録する。
7. CUDAで使用する方式はGPU側の面流束APIとCPU/CUDA一致テストも追加する。
8. 一様場保存、周期移流、衝撃波管、格子収束の順で検証する。

KEEP/WENOハイブリッドではDucros-pressureセンサーを独立させ、ディスパッチ層の下に
ハイブリッド複合流束を実装している。KEEPとWENO5-Z/Roeの内部式は変更せず、面ごとの
混合係数だけをハイブリッド層が担当する。詳細は
`docs/NSE_HYBRID_FLUX.md`を参照する。

## 初期条件を追加する手順

1. `src/init/mod_init_<name>.f90`を追加する。
2. 保存変数`Q = [rho, rho*u, rho*v, rho*w, rho*E]`を局所内部セルへ設定する手続きを公開する。
3. `mod_nse_initial_conditions.f90`へ`use`と`select case`を追加する。
4. `solver_manifest.yaml`と`config/module_catalog.yaml`へ登録する。
5. ケース入力の`initial_condition`で名前を選ぶ。

流れ場固有パラメータが増えた場合は、式の中へ固定値として埋め込まず、`nse_config`と`/nse/`入力へ追加する。

## 境界条件を追加する手順

面ごとの組合せを可能にするため、境界条件ごとに同名の`mod_nse_boundary`を交換する
従来方式から、1つのdispatcherへ複数の一次元面演算子を登録する方式へ移行する。

1. `src/boundary/operators/mod_boundary_<name>.f90`へ、一つの法線方向と面向きを受け取る
   境界演算子を追加する。
2. `mod_nse_boundary`のdispatcherへ入力名、検証、CPU呼出し先を登録する。
3. 物理境界とMPI内部境界を区別し、ghostセルを一度のAPI呼出しで完成させる。
4. x→y→zの演算子合成順を守り、混合境界の辺・角を完成させる。
5. CUDA対応では同じ演算子のdevice実装を追加し、CPUへ暗黙fallbackしない。
6. CMake、マニフェスト、モジュールカタログの境界backendを`runtime`へ更新する。
7. 一様場、保存量、反射率、辺・角、MPI分割位置、CPU/CUDA一致テストを追加する。

個別方式を追加するたびに新しい`mod_nse_boundary`を作ってコンパイル時に一つだけ選ぶ
設計へ戻してはならない。面別設定の正規入力と無反射／鏡像境界の数値契約は
`NSE_BOUNDARY_CONDITIONS.md`に従う。

## 入力とビルドの選択

生成後のFortran入力`/nse/`では、対流・粘性方式に加えて6面の境界種別と
無反射面の具体的な基準状態を指定する。

```fortran
convective_scheme = 'weno5z_roe'
viscous_scheme = 'central6'
boundary_x_min = 'non_reflecting'
boundary_x_max = 'non_reflecting'
boundary_y_min = 'periodic'
boundary_y_max = 'periodic'
boundary_z_min = 'periodic'
boundary_z_max = 'periodic'
time_integrator = 'ssprk3'
```

対流項は`keep2`、`keep6`、`weno5z_roe`、`hybrid`で指定する。`keep`や`weno`単独の指定は
使用できない。4方式はCPU逐次、MPI、OpenMP、MPI+OpenMP、単一GPU CUDA版、
MPI＋CUDA版で利用できる。
`hybrid`の構成流束とセンサーは`NSE_HYBRID_FLUX.md`の設定で選択する。

CMakeでは対応するバックエンドを選択する。

```powershell
cmake -S . -B build `
  -DNSE_CONVECTIVE_BACKEND=runtime `
  -DNSE_VISCOUS_SCHEME=central6 `
  -DNSE_BOUNDARY_SCHEME=runtime
```

Linux（bash）:

```bash
cmake -S . -B build \
  -DNSE_CONVECTIVE_BACKEND=runtime \
  -DNSE_VISCOUS_SCHEME=central6 \
  -DNSE_BOUNDARY_SCHEME=runtime
```

入力とコンパイル済みバックエンドが一致しない場合、実行開始時に停止する。

case YAMLの`boundary.faces`と`reference_states`を
`case_input.py`が検証し、名前付き基準状態を各面の具体的な密度、速度、圧力へ解決して
Fortran入力へ渡す。Fortran側は面ごとの種別・密度・速度・圧力配列として保持し、YAMLの名前解決を
重複実装しない。ビルド時の境界backendは上記の`runtime`へ統一する。

移行期間は旧`boundary_condition = 'periodic'`を6面周期へ展開する。面別指定と旧指定の
同時使用は入力エラーとし、暗黙の優先順位を設けない。CPUおよびCUDA profileはともに
`runtime` backendを使い、実行時に6面の設定を解決する。

## 現段階の制約

- 旧MPIハロー交換との互換性により、面別境界は`nghost = 3`、`nv = 5`を要求する。
- 面別周期／無反射／鏡像／Dirichlet境界はCPU/MPI/OpenMP、単一GPU CUDA、MPI＋CUDA版で利用できる。
- `viscous_scheme = 'central6'`は一定粘性係数、Stokesの仮定、Fourier熱伝導を用いる。`reynolds_number`と`prandtl_number`は粘性・熱伝導項へ反映される。
- `viscous_scheme = 'none'`を選ぶと、同じ実行プロファイルで非粘性計算を行える。
- `central6`は3層のghostセルを必要とし、CPU/MPI/OpenMP版、単一GPU CUDA版、MPI＋CUDA版で利用できる。
- MPI分割は既存のy-z二次元分割を維持している。
- `default`は互換性のためTaylor–Greenへ対応付けている。
- `hit_spectral`は2DECOMP&FFT版を選んだMPIビルドで利用できる。
- `imported_turbulence`は通常CPU版、2DECOMP&FFT版、単一GPU版、MPI＋CUDA版で利用できる。入力SLFの作成方法は`NSE_IMPORTED_TURBULENCE.md`に従う。
- `shock_turbulence_interaction`は`imported_turbulence`を再利用し、`mod_init_shock_turbulence`が平面衝撃波の配置とDirichlet駆動面の整合性を検査する。
- `shock_tube_turbulence_interaction`は同じ可搬SLF読込みを再利用し、`mod_init_shock_tube_turbulence`が有限高圧室、鏡像閉端、低圧無反射端の整合性を検査する。
- WENO5-Z/RoeはCPU版、単一GPU CUDA版、MPI＋CUDA版で利用できる。positivity-preserving limiterは未実装である。
- MPI＋CUDA時間発展のhalo通信はCUDA-aware MPIを要求しないhost staging方式である。`cuda_mpi_cufftmp`では、同じY-Z分割をcuFFTMpへ渡して分散HIT初期化とFFT forcingを実行する。
- SSPRK3とCFL時間刻みのみを実装している。

この制約は各モジュールの検証手続きで明示的に検査し、未対応条件を黙って計算しない。
