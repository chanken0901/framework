# NSEモジュール設計

## 目的

NSEソルバーの初期条件、移流スキーム、粘性項、境界条件を、mainプログラムを変更せずに交換・追加できる構造にする。
現在は5変数の圧縮性Navier-Stokes方程式を扱い、2次・6次精度KEEP対流項と、
非粘性計算または6次精度中心差分による粘性計算を選択できる。

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
| 初期条件選択 | `src/init/mod_nse_initial_conditions.f90` | `initialize_nse_state` |
| Taylor–Green | `src/init/mod_init_taylor_green.f90` | `initialize_taylor_green` |
| 分散FFT HIT | `src/init/mod_init_hit_spectral_2decomp.f90` | `initialize_hit_spectral` |
| 周期境界 | `src/boundary/mod_boundary_periodic.f90` | `apply_nse_boundary`ほか |
| 選択可能KEEP流束 | `src/numerics/convective/mod_convective_keep.f90` | 2次・6次精度の`compute_convective_flux` |
| 粘性項なし | `src/numerics/viscous/mod_viscous_none.f90` | 非粘性計算用の空実装 |
| 6次精度粘性項 | `src/numerics/viscous/mod_viscous_central6.f90` | `add_viscous_rhs`、粘性時間刻み評価 |
| 空間演算 | `src/numerics/mod_nse_spatial_operator.f90` | 境界、流束、発散、粘性項を合成 |
| SSPRK3 | `src/time/mod_nse_time_integration.f90` | CFL時間刻みと3段Runge–Kutta |

## 交換可能モジュールの規約

移流、粘性、境界の各バックエンドは、分類ごとに同じFortranモジュール名と公開APIを実装する。
CMakeは各分類から1つだけを選択してコンパイルする。この方式により、空間演算とmainは個別スキーム名を知らなくてよい。

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

Fortranモジュール名は`mod_nse_boundary`とし、次を公開する。

```fortran
apply_nse_boundary
validate_boundary_scheme
boundary_required_ghost_cells
boundary_scheme_name
```

## WENOを追加する手順

1. `src/numerics/convective/mod_convective_weno.f90`を追加する。
2. `mod_convective_scheme`の公開APIをすべて実装する。
3. 必要なゴーストセル数を`convective_required_ghost_cells`から返す。
4. `CMakeLists.txt`の`NSE_CONVECTIVE_BACKEND`候補とソース選択へ`weno_family`を追加する。
5. `solver_manifest.yaml`と`config/module_catalog.yaml`へモジュールを登録する。
6. 入力の`convective_scheme = 'weno'`とビルド設定を一致させる。
7. 一様場保存、周期移流、衝撃波管、格子収束の順で検証する。

WENOの再構築、特性分解、数値流束はWENOモジュール内に閉じ込める。
時間積分、MPI分割、出力コードへWENO固有処理を置かない。

## 初期条件を追加する手順

1. `src/init/mod_init_<name>.f90`を追加する。
2. 保存変数`Q = [rho, rho*u, rho*v, rho*w, rho*E]`を局所内部セルへ設定する手続きを公開する。
3. `mod_nse_initial_conditions.f90`へ`use`と`select case`を追加する。
4. `solver_manifest.yaml`と`config/module_catalog.yaml`へ登録する。
5. ケース入力の`initial_condition`で名前を選ぶ。

流れ場固有パラメータが増えた場合は、式の中へ固定値として埋め込まず、`nse_config`と`/nse/`入力へ追加する。

## 境界条件を追加する手順

1. `src/boundary/mod_boundary_<name>.f90`を追加する。
2. `mod_nse_boundary`の公開APIをすべて実装する。
3. 物理境界とMPI内部境界を区別し、ゴーストセルを一度のAPI呼び出しで完成させる。
4. CMake、マニフェスト、モジュールカタログへ候補を登録する。
5. 保存量、反射条件、MPI分割位置に関する単体・並列テストを追加する。

## 入力とビルドの選択

入力`/nse/`では次を指定する。

```fortran
convective_scheme = 'keep6'
viscous_scheme = 'central6'
boundary_condition = 'periodic'
time_integrator = 'ssprk3'
```

対流項は`keep2`または`keep6`で指定します。`keep`単独の指定は使用できません。

CMakeでは対応するバックエンドを選択する。

```powershell
cmake -S . -B build `
  -DNSE_CONVECTIVE_BACKEND=keep_family `
  -DNSE_VISCOUS_SCHEME=central6 `
  -DNSE_BOUNDARY_SCHEME=periodic
```

入力とコンパイル済みバックエンドが一致しない場合、実行開始時に停止する。

## 現段階の制約

- 旧MPIハロー交換との互換性により、周期境界は`nghost = 3`、`nv = 5`を要求する。
- `viscous_scheme = 'central6'`は一定粘性係数、Stokesの仮定、Fourier熱伝導を用いる。`reynolds_number`と`prandtl_number`は粘性・熱伝導項へ反映される。
- `viscous_scheme = 'none'`を選ぶと、同じ実行プロファイルで非粘性計算を行える。
- `central6`は3層のghostセルを必要とし、CPU/MPI/OpenMP版と単一GPU CUDA版で利用できる。
- MPI分割は既存のy-z二次元分割を維持している。
- `default`は互換性のためTaylor–Greenへ対応付けている。
- `hit_spectral`は2DECOMP&FFT版を選んだMPIビルドで利用できる。
- SSPRK3とCFL時間刻みのみを実装している。

この制約は各モジュールの検証手続きで明示的に検査し、未対応条件を黙って計算しない。
