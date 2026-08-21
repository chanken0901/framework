# GPE/GP3Dコード全体仕様書

**版:** 1.1

**更新日:** 2026-08-21
**対象:** SolverLibrary GPE/gp3d  
**正本:** `C:\Users\Owner\Documents\Codex\FrameWork`

---

## 1. 文書の目的

本書は、三次元Gross–Pitaevskii方程式ソルバーGP3Dの現行実装を、物理モデル、数値解法、ソフトウェア構造、入力、出力、並列化、GPU化、ビルド、テスト、制約の観点から定義するコード仕様書です。

主な用途は次のとおりです。

| 用途 | 内容 |
|---|---|
| 実装理解 | 各モジュールの責任とデータの流れを確認する |
| case設計 | `case.yaml` の値が計算へどう反映されるか確認する |
| 保守 | 変更対象と影響範囲を判断する |
| 移植 | Windows、Linux、MPI、CUDA、cuFFTMpの要件を確認する |
| 検証 | 数値結果とバックエンドの同等性を確認する |
| Codex連携 | 本書を渡して現行構造を維持したコード変更を依頼する |

本書は2026-08-21時点の次のソースを直接確認して更新しています。

```text
C:\Users\Owner\Documents\Codex\FrameWork\SolverLibrary\GPE\gp3d
```

GP3DはFrameWorkモノレポ内でGit管理されています。厳密な版はリポジトリで
`git rev-parse HEAD`と`git status --short`を実行して確認します。

---

## 2. 適用範囲

### 2.1 実装済み

| 分類 | 実装 |
|---|---|
| 方程式 | 零温度の三次元Gross–Pitaevskii方程式 |
| 時間発展 | 実時間および固定ノルム虚時間のStrang split-operator法 |
| 初期緩和 | Taylor–Green速度場を用いる無拘束・半陰的ARGLE |
| 空間離散 | 周期直交格子と擬スペクトル法 |
| CPU逐次 | 参照DFT、FFTW3 |
| CPU並列 | MPI zスラブ分割、または2DECOMP&FFTによるペンシル分割FFT |
| 単一GPU | CUDAカーネルとcuFFT |
| 複数GPU | MPI、CUDA、cuFFTMp、1 MPI rankにつき1 GPU |
| 初期条件 | Gaussian、直線渦、Thomas–Fermi渦、渦輪、渦タングル、量子Taylor–Green |
| 入出力 | Fortran namelist、SLF1、meta.json |
| 再スタート | global SLFおよびrank分割SLF、MPI分割数変更対応 |
| 後処理 | SLFからVTI/PVDへの変換 |

### 2.2 適用外

現行コードには次のモデルは含まれません。

| 適用外 | 備考 |
|---|---|
| 有限温度GPE | 散逸GPE、ZNG、SPGPEなどは未実装 |
| 多成分GPE | 波動関数は1成分のみ |
| 回転座標系 | 明示的な角運動量項は未実装 |
| 吸収境界 | 周期境界のみ |
| 適応格子 | 一様直交格子のみ |
| 非周期FFT | 未実装 |
| 自動時間刻み | `dt` は固定 |

---

## 3. システム概要

GP3Dは、共通の物理・入出力モジュールと、同じ公開APIを持つ交換可能なFFT、MPI、GPU実装から構成されます。

```text
case.yaml
   |
   v
ScriptLibrary input adapter
   |
   v
input.nml
   |
   v
main program
   |
   +-- gp3d_input -------- grid / params / state / output config
   +-- gp3d_initial_conditions
   +-- gp3d_restart
   +-- gp3d_solver ------- Split-step / ARGLE / diagnostics
   |      |
   |      +-- gp3d_fft --- DFT / FFTW / MPI slab / MPI pencil
   |      +-- gp3d_mpi --- stub / real MPI
   |
   +-- gp3d_gpu ---------- CUDA / cuFFTMp backend
   |
   +-- gp3d_io ----------- SLF / meta.json
```

バックエンド差し替えは、CMakeが同じFortranモジュール名の実装から一つを選ぶことで行います。

| 公開モジュール | 選択可能な実装 |
|---|---|
| `gp3d_mpi` | 逐次スタブ、実MPI |
| `gp3d_fft` | 逐次DFT、逐次FFTW、MPIスラブFFT、MPIペンシルFFT |
| `gp3d_local_fft` | 1次元DFT、1次元FFTW |
| `gp3d_gpu` | 単一GPU CUDA、MPI/cuFFTMp |

---

## 4. 物理モデル

### 4.1 Gross–Pitaevskii方程式

現行コードが解く方程式は次です。

```text
i hbar dpsi/dt =
  [ -(hbar^2 / 2m) nabla^2 + V(x,y,z) + g abs(psi)^2 ] psi
```

| 記号 | 意味 |
|---|---|
| `psi` | 複素波動関数 |
| `rho = abs(psi)^2` | 粒子数密度または無次元密度 |
| `hbar` | 換算Planck定数 |
| `m` | 粒子質量 |
| `g` | 接触相互作用係数 |
| `V` | 外部ポテンシャル |

原子種の違いを物理単位で扱う場合は、主に `mass`、`g`、`hbar`、トラップ係数、規格化 `norm` へ反映します。散乱長から `g` を作る単位変換はコード内で自動化されていません。

### 4.2 無次元パラメータ

`use_dimensionless_parameters = true` の場合、内部パラメータは次へ置き換えられます。

```text
hbar = 1
mass = 1 / (2 alpha)
g = beta
```

したがって時間発展方程式は次の形になります。

```text
i dpsi/dt =
  [ -alpha nabla^2 + V + beta abs(psi)^2 ] psi
```

無次元モードでは、`&gpe` の `hbar`、`mass`、時間発展用の `g` は使用されません。`alpha` と `beta` が優先されます。

### 4.3 外部ポテンシャル

全caseで、状態初期化時に次の調和ポテンシャルを構築します。

```text
V(x,y,z) = 0.5 [ (wx x)^2 + (wy y)^2 + (wz z)^2 ]
```

これは現行コードの式です。標準的な物理式に現れる質量係数 `m` は、このポテンシャル式には自動で掛かりません。物理単位で厳密な調和トラップを設定する場合は、入力係数の定義を揃えるか、将来ポテンシャルモジュールを拡張する必要があります。

`wx = wy = wz = 0` とすれば、一様な周期領域になります。

### 4.4 ノルム

波動関数のノルムは次で計算します。

```text
N = integral abs(psi)^2 dV
  approximately sum(abs(psi)^2) dx dy dz
```

`norm` は初期状態の目標ノルム、および通常の虚時間発展で各ステップ後に戻す目標ノルムです。

### 4.5 エネルギー

診断量として次を計算します。

```text
E = integral [
      (hbar^2 / 2m) abs(grad psi)^2
      + V abs(psi)^2
      + (g / 2) abs(psi)^4
    ] dV
```

運動エネルギーはFourier空間で、ポテンシャル・相互作用エネルギーは実空間で評価します。

---

## 5. 空間格子と境界条件

### 5.1 格子

格子は一様なセル中心直交格子です。

```text
dx = (x_max - x_min) / nx
x(i) = x_min + (i - 0.5) dx
```

`y`、`z` も同様です。端点は重複しません。

### 5.2 境界条件

FFTを使用するため、全方向が周期境界です。

| 要件 | 内容 |
|---|---|
| 周期性 | `psi` とその空間構造が対向境界で接続すること |
| 端点 | `x_max` 自体の格子点は持たない |
| トラップ | 周期境界と調和ポテンシャルを併用する場合、境界で密度が十分小さいことが望ましい |
| Taylor–Green | 周期領域と整合するよう角度座標へ写像する |

境界から波が発生して見える場合、非周期な初期条件、境界で無視できない密度、急な初期緩和、格子不足などを確認します。

### 5.3 波数

FFTのモードは次です。

```text
k = 2 pi mode / L
mode = 0, 1, ..., floor(n/2), negative modes
```

偶数格子では、インデックス `n/2 + 1` が正のNyquistモードとして定義されます。

### 5.4 MPI分割

実空間の所有配置とSLF出力は、FFT方式によらずzスラブです。CPU/MPIと
MPI/cuFFTMpのFFT内部配置は
`FFT_DECOMPOSITION=slab`または`pencil`から選択します。

| 項目 | 内容 |
|---|---|
| ローカル形状 | `(nx, ny, local_nz)` |
| 分割 | `nz` をrank数でほぼ均等分割 |
| 余り | 小さいrankから1面ずつ追加 |
| 制約 | `nprocs <= nz` |
| スラブFFT追加制約 | `nprocs <= ny` |
| CPUペンシルFFT追加制約 | 2DECOMP&FFTの`p_row <= min(nx,ny)`、`p_col <= min(ny,nz)` |
| cuFFTMpペンシル追加制約 | `p_y <= ny`、`p_z <= nz`、`p_y * p_z = nprocs` |

したがってスラブCPU分散FFTとスラブcuFFTMpでは`nprocs <= min(ny,nz)`、
CPU／cuFFTMpペンシルFFTではzスラブ保存配列に由来する`nprocs <= nz`が必要です。
ペンシルのプロセス格子は各バックエンドがrank数から自動決定します。

---

## 6. 状態と共通データ型

型定義は `src/common/gp3d_types.f90` にあります。

### 6.1 `gp3d_grid_t`

| メンバー | 内容 |
|---|---|
| `nx, ny, nz` | 全体格子数 |
| `local_nz` | 現rankのz方向格子数 |
| `k_start, k_end` | 現rankが所有する全体zインデックス |
| `rank, nprocs` | 並列配置 |
| `lx, ly, lz` | 周期長 |
| `dx, dy, dz` | 格子幅 |
| `x, y, z` | 全体座標配列 |
| `kx, ky, kz` | 全体波数配列 |

各rankは `z` と `kz` も全体配列として保持し、局所ループでは `k_start` から全体インデックスへ変換します。

### 6.2 `gp3d_params_t`

時間ループが直接参照する導出済みパラメータです。

| メンバー | 内容 |
|---|---|
| `dt` | 固定時間刻み |
| `mass, hbar, g` | 実際に使用する物理係数 |
| `norm` | 目標ノルム |
| `nsteps` | 今回追加で進めるステップ数 |
| `output_every` | 出力間隔 |
| `imaginary_time` | 通常Split-stepを虚時間化するフラグ |

### 6.3 `gp3d_state_t`

| 配列 | 型と形状 | 内容 |
|---|---|---|
| `psi` | complex float64 `(nx, ny, local_nz)` | 波動関数 |
| `potential` | float64 `(nx, ny, local_nz)` | 外部ポテンシャル |

---

## 7. Split-operator時間積分

### 7.1 演算子分割

Hamiltonianを次へ分けます。

```text
T = -(hbar^2 / 2m) nabla^2
W(psi) = V + g abs(psi)^2
```

1ステップはStrang分割で進めます。

```text
psi(t + dt) approximately
  exp[-i dt W / (2 hbar)]
  exp[-i dt T / hbar]
  exp[-i dt W / (2 hbar)]
  psi(t)
```

実装順序は次です。

1. 実空間で局所項を半ステップ適用する。
2. 順FFTする。
3. Fourier空間で運動項を1ステップ適用する。
4. 逆FFTする。
5. 実空間で局所項を半ステップ適用する。
6. 虚時間の場合だけノルムを規格化する。

### 7.2 局所半ステップ

各格子点で独立に次を掛けます。

実時間:

```text
psi <- psi exp[-i dt (V + g abs(psi)^2) / (2 hbar)]
```

虚時間:

```text
psi <- psi exp[-dt (V + g abs(psi)^2) / (2 hbar)]
```

有限差分・有限体積の空間フラックス計算は行いません。

### 7.3 運動項

Fourier係数に対して次を掛けます。

実時間:

```text
psi_k <- psi_k exp[
  -i dt (hbar^2 / 2m) (kx^2 + ky^2 + kz^2) / hbar
]
```

虚時間:

```text
psi_k <- psi_k exp[
  -dt (hbar^2 / 2m) (kx^2 + ky^2 + kz^2) / hbar
]
```

### 7.4 精度

| 項目 | 仕様 |
|---|---|
| 時間分割 | 対称Strang分割 |
| 形式精度 | 滑らかな解に対し時間2次 |
| 空間精度 | 周期的で滑らかな解に対する擬スペクトル精度 |
| 時間刻み | 固定 |
| 非線形評価 | 各局所半ステップ開始時の現在の `psi` |
| dealiasing | 未実装 |

### 7.5 虚時間発展

`imaginary_time = true` の通常Split-stepでは、指数因子から虚数単位を外し、各ステップ後に `norm` へ規格化します。

これは固定ノルムの基底状態・準安定状態探索です。実時間の物理的な時刻とは異なります。

ARGLEはこの通常虚時間法とは別の処理です。

---

## 8. ARGLE仕様

### 8.1 目的

ARGLEは、量子Taylor–Green初期条件から強い音波成分を減らし、指定したTaylor–Green速度場に対応する波動関数を実時間計算前に緩和するために使用します。

`argle_enabled = true` かつ再スタートでない場合、通常の時間ループより前に実行されます。

### 8.2 緩和方程式

`alpha = hbar / (2 mass)` とします。

無次元モード:

```text
dpsi/dtau =
  alpha nabla^2 psi
  - i u dot grad(psi)
  + [ beta (1 - abs(psi)^2) - V - abs(u)^2 / (4 alpha) ] psi
```

物理パラメータモード:

```text
dpsi/dtau =
  alpha nabla^2 psi
  - i u dot grad(psi)
  + [ (mu - V - g abs(psi)^2) / hbar
      - abs(u)^2 / (4 alpha) ] psi
```

### 8.3 Taylor–Green速度場

周期角度座標を使い、次を与えます。

```text
ux = U sin(x) cos(y) cos(z)
uy = -U cos(x) sin(y) cos(z)
uz = 0
```

`U` は `tg_velocity_amplitude` です。物理座標は各領域長を使って `[-pi, pi)` の角度へ写像します。

### 8.4 半陰的更新

Laplacian項をCrank–Nicolson型に扱い、反応・移流項を陽的に扱います。

```text
psi_k(new) =
  [ (1 - 0.5 dtau alpha k^2) psi_k(old)
    + dtau RHS_k(old) ]
  / [1 + 0.5 dtau alpha k^2]
```

ARGLEは無拘束法であり、各反復後のノルム規格化を行いません。

### 8.5 収束判定

```text
max_delta_rate = max(abs(psi_new - psi_old)) / argle_dtau
```

`argle_tolerance > 0` かつ `max_delta_rate < argle_tolerance` で終了します。`argle_tolerance = 0` は収束による早期終了を無効にします。

### 8.6 FFT回数

ARGLE 1反復では概ね次の5回の三次元FFTを使用します。

1. `psi_old` の順FFT
2. x微分用の逆FFT
3. y微分用の逆FFT
4. RHSの順FFT
5. 更新後 `psi` の逆FFT

### 8.7 出力に関する注意

| 設定 | 実際の動作 |
|---|---|
| `argle_write_seed` | ARGLE前の解析的seedを `output/argle_seed` へ保存 |
| `argle_output_every` | 反復ログの表示間隔 |
| ARGLE中間場 | 現行コードは保存しない |
| 再スタート | ARGLEをスキップする |

`argle_write_seed` は緩和後の状態ではなく、緩和前の比較用seedを保存します。

---

## 9. 初期条件

初期条件の選択は `&simulation initial_condition`、フレームワークでは `flow.type` で行います。

### 9.1 対応一覧

| 正式名 | エイリアス | 背景と構造 |
|---|---|---|
| `gaussian` | 大文字表記も一部対応 | 三次元Gaussian |
| `uniform_vortex` | `vortex_uniform` | 一様背景中のz方向直線渦 |
| `tf_vortex` | `thomas_fermi_vortex`, `vortex` | Thomas–Fermi背景中のz方向直線渦 |
| `vortex_ring` | `ring` | Gaussian背景中のz軸法線渦輪 |
| `vortex_tangle` | `random_vortices`, `quantum_turbulence` | Gaussian背景中のランダムz方向直線渦群 |
| `ring_tangle` | `vortex_ring_tangle`, `random_rings` | 一様背景中のランダム配向渦輪群 |
| `quantum_taylor_green` | `taylor_green`, `tg` | 周期Taylor–Green量子渦配置 |
| `restart_slf` | `restart`, `slf` | SLFから復元 |

小文字で指定することを標準とします。

### 9.2 Gaussian

```text
psi = exp[
  -0.5 ((x/sigma0)^2 + (y/sigma0)^2 + (z/sigma0)^2)
]
```

通常は生成後に `norm` へ規格化されます。

### 9.3 直線渦

z方向へ伸びる渦を、既存の背景波動関数へ乗算します。

```text
phase = charge atan2(y - y0, x - x0)
core = [ s / sqrt(1 + s^2) ]^abs(charge)
s = radius / healing_length
```

`uniform_vortex` の背景密度は入力から変更できず、現在は `density0 = 1` です。

### 9.4 Thomas–Fermi渦

背景密度は次です。

```text
rho = max((mu - V) / g, 0)
```

その上へ直線渦を刻印します。

無次元モードでも、この初期背景計算は `cfg%g` を参照します。一方、時間発展は `beta` を `g` として使います。無次元Thomas–Fermi渦を使う場合は、この差を理解して `g` と `beta` を整合させる必要があります。

### 9.5 渦輪

円筒座標でリングからの距離を作り、直線渦と同様の芯と位相を刻印します。

`vortex_ring` は原点中心、z軸法線です。内部API `gp3d_imprint_vortex_ring_oriented` は任意中心・任意法線に対応します。

### 9.6 ランダム直線渦群

`vortex_tangle` は、xy平面上のランダム位置にz方向直線渦を配置します。符号は正負交互です。

初期状態は等方的な三次元渦タングルではなく、平行な直線渦群です。その後の時間発展で変形・相互作用させる設計です。

### 9.7 ランダム渦輪群

`ring_tangle` は、中心、半径、法線を乱数で選ぶ三次元渦輪群です。符号は正負交互です。

ランダム配向を含むため、多数の絡み合った量子渦を作る初期条件としては `vortex_tangle` より直接的です。

### 9.8 量子Taylor–Green

周期Taylor–Green対称性を持つ4本組の素渦を、変換座標 `lambda`、`mu` 上で合成します。

素渦は次の形です。

```text
psi_e = tanh(r / (sqrt(2) xi)) (lambda + i mu) / r
```

4個の積を `winding` 乗し、背景密度を掛けます。

無次元モードでは次を導出します。

```text
healing_length = sqrt(alpha / beta)
sound_speed = sqrt(2 alpha beta)
Mach = tg_velocity_amplitude / sound_speed
```

`tg_auto_winding = true` の場合:

```text
winding = floor(tg_velocity_amplitude / (2 pi alpha))
```

結果が1未満の場合はエラーとします。

### 9.9 位相ノイズ

`phase_noise > 0` の場合、各格子点へ `[-phase_noise, phase_noise]` の決定論的位相を加えます。

全体格子インデックスと `random_seed` から値を作るため、MPI分割数を変えても同じ全体格子点には同じノイズが入ります。

---

## 10. モジュール構成

### 10.1 ディレクトリ

| パス | 責任 |
|---|---|
| `src/common` | 精度、定数、共有型 |
| `src/grid` | 周期格子、波数、zスラブ分割 |
| `src/init` | 初期波動関数 |
| `src/solver` | Split-step、ARGLE、診断量、タイミング |
| `src/fft` | 逐次・MPI FFTバックエンド |
| `src/mpi` | MPI抽象化 |
| `src/gpu` | CUDA/cuFFT、cuFFTMp、Fortran C binding |
| `src/io` | namelist、SLF、meta.json、再スタート |
| `src/main` | CPU/MPI、単一GPU、複数GPUの制御 |
| `tests` | 数値・入出力テスト |
| `tools` | 実行ワークフローと後処理 |
| `docs` | スパコン検証手順 |

### 10.2 主要モジュール

| モジュール | 主要責任 |
|---|---|
| `gp3d_types` | 共有データ構造 |
| `gp3d_grid` | 格子と波数の生成 |
| `gp3d_input` | namelist読み込みと問題構築 |
| `gp3d_initial_conditions` | 初期条件生成 |
| `gp3d_solver` | CPU時間積分と診断 |
| `gp3d_fft` | 三次元FFT共通API |
| `gp3d_local_fft` | MPI FFTの局所1次元変換 |
| `gp3d_mpi` | MPIコンテキストと集団通信 |
| `gp3d_gpu` | GPU共通API |
| `gp3d_io` | SLFとmeta.json |
| `gp3d_restart` | SLF再開 |

### 10.3 公開CPUソルバーAPI

| 手続き | 役割 |
|---|---|
| `gp3d_state_allocate` | `psi` と `potential` を確保 |
| `gp3d_set_harmonic_potential` | 調和ポテンシャル設定 |
| `gp3d_set_gaussian_initial_state` | Gaussian初期化 |
| `gp3d_normalize` | ノルム規格化 |
| `gp3d_density_norm` | 全体ノルム |
| `gp3d_energy` | 全体エネルギー |
| `gp3d_step_split_operator` | Split-step 1ステップ |
| `gp3d_relax_taylor_green_argle` | CPU ARGLE |
| `gp3d_report_step_timing` | タイミング集約 |

### 10.4 FFT共通API

すべてのFFT実装は次を公開します。

| 手続き | 契約 |
|---|---|
| `gp3d_fft_init` | plan生成 |
| `gp3d_fft_forward` | 無規格化の順変換 |
| `gp3d_fft_inverse` | `1 / (nx ny nz)` 規格化付き逆変換 |
| `gp3d_fft_finalize` | plan解放 |

### 10.5 GPU共通API

| 手続き | 契約 |
|---|---|
| `gp3d_gpu_init` | GPU、FFT plan、常駐配列を確保 |
| `gp3d_gpu_upload` | 初期状態と係数を転送 |
| `gp3d_gpu_download` | 出力用に `psi` を取得 |
| `gp3d_gpu_step` | GPU常駐Split-step |
| `gp3d_gpu_relax_taylor_green_argle` | GPU常駐ARGLE |
| `gp3d_gpu_diagnostics` | normとenergy |
| `gp3d_gpu_finalize` | GPU資源解放 |

---

## 11. 実行フロー

### 11.1 CPU逐次・CPU MPI

`src/main/main.f90` の処理順は次です。

1. MPI実装または逐次スタブを初期化する。
2. `input.nml` を読む。
3. 実際のMPI rank数を設定へ反映する。
4. 格子、状態、ポテンシャル、初期条件または再スタートを構築する。
5. FFT planを生成する。
6. 必要ならARGLE前seedを保存する。
7. 再スタートでなければARGLEを実行する。
8. `meta.json` と初期SLFを保存する。
9. 初期normとenergyを表示する。
10. `nsteps` 回の時間発展を行う。
11. 出力ステップで診断とSLF保存を行う。
12. タイミングを表示する。
13. FFTとMPIを終了する。

### 11.2 単一GPU

`src/main/main_cuda.f90` は、ホストで初期化後に次を行います。

1. CUDA device、cuFFT plan、GPU配列を確保する。
2. `psi`、`potential`、波数をGPUへ転送する。
3. ARGLEをGPUで実行する。
4. 時間ループをGPU常駐のまま進める。
5. 出力時だけ `psi` をホストへ戻す。
6. 診断量はGPU上で集約する。

### 11.3 複数GPU

`src/main/main_cufftmp.f90` は、各MPI rankが一つのGPUとzスラブを所有します。

1. rank内ローカル番号からGPUを選択する。
2. cuFFTMp分散planを作る。
3. 実空間zスラブをGPUへ転送する。
4. 順FFT後のcuFFTMp yスラブまたはxペンシル配置をカーネルが直接処理する。
5. 出力時だけローカルzスラブをダウンロードする。
6. 診断量と規格化にはMPI Allreduceを使う。

---

## 12. CPU FFTバックエンド

### 12.1 逐次DFT

`src/fft/gp3d_fft.f90` は外部ライブラリを使わない直接三次元DFTです。

| 項目 | 仕様 |
|---|---|
| 用途 | 小格子の参照解、テスト |
| 計算量 | 全格子点数を `N` として概ね `O(N^2)` |
| 実用規模 | 非常に小さい格子のみ |
| 外部依存 | なし |

### 12.2 逐次FFTW

`src/fft/gp3d_fft_fftw.f90` はFFTW3のC APIを使用します。

| 項目 | 仕様 |
|---|---|
| plan | `FFTW_ESTIMATE` |
| 配列 | complex float64 |
| 次元順 | Fortran列優先に合わせて `(nz, ny, nx)` をFFTWへ渡す |
| 逆変換 | 実行後に全格子点数で除算 |
| 外部依存 | `libfftw3` |

### 12.3 MPIスラブ分散FFT

`src/fft/gp3d_fft_mpi.f90` の処理は次です。

1. zスラブ上でx方向局所FFT
2. zスラブ上でy方向局所FFT
3. `MPI_Alltoallv` でzスラブからyスラブへ転置
4. yスラブ上でz方向局所FFT
5. `MPI_Alltoallv` でyスラブからzスラブへ戻す

1回の三次元変換につき2回のAlltoallvがあります。Split-step 1ステップでは順変換と逆変換があるため、合計4回のAlltoallvを使用します。

局所1次元変換は、参照DFTまたはFFTWから選択します。

### 12.4 MPIペンシル分散FFT

`src/fft/gp3d_fft_pencil_2decomp.f90`は2DECOMP&FFTとFFTWを使用します。
公開APIの入力・出力は従来互換のzスラブですが、FFT内部は2次元プロセス格子です。

順変換は次の順です。

1. GP3D zスラブから2DECOMP Xペンシルへ`MPI_Alltoallv`で再分配する。
2. 2DECOMP&FFTでXペンシルからYペンシル、Zペンシルへ転置しながら三次元FFTする。
3. ZペンシルからGP3D zスラブへ`MPI_Alltoallv`で戻す。

逆変換は逆順に処理し、最後に`1 / (nx ny nz)`で規格化します。物理場とSLFを
zスラブのまま残すため、既存の再スタート・後処理と互換です。一方、1回の変換で
外側2回と2DECOMP内部2回の計4回の転置通信を行うため、少ないrankではスラブ版より
高速とは限りません。現段階では保存配置の制約により`nprocs <= nz`も残ります。

### 12.5 CPUメモリ確保

現行実装は次を時間ループ中に確保・解放します。

| 場所 | 一時配列 |
|---|---|
| Split-step | `psi_k` を毎ステップallocate/deallocate |
| MPI三次元FFT | zスラブ、yスラブ、line bufferを変換ごとにallocate/deallocate |
| MPI転置 | send/recv bufferを転置ごとにallocate/deallocate |
| ペンシル三次元FFT | X/Zペンシルと外側再分配bufferを変換ごとにallocate/deallocate |
| energy | `psi_k` を診断ごとにallocate/deallocate |

正しさを優先した現行仕様であり、性能最適化ではplanまたはworkspaceへ常設する余地があります。

---

## 13. 単一GPU CUDA仕様

### 13.1 構成

| 層 | ファイル |
|---|---|
| Fortran API | `src/gpu/gp3d_gpu_cuda.f90` |
| C ABI/CUDA | `src/gpu/gp3d_cuda_bridge.cu` |
| main | `src/main/main_cuda.f90` |

### 13.2 GPU常駐配列

通常時間発展では、少なくとも次をGPUに保持します。

| 配列 | 内容 |
|---|---|
| `psi` | 実空間波動関数 |
| `spectral` | Fourier空間work |
| `potential` | 外部ポテンシャル |
| `kx, ky, kz` | 波数 |
| `reductions` | norm、energy、ARGLE指標 |

ARGLE使用時は `old_psi`、`grad_x`、`grad_y`、`rhs` を追加確保します。

### 13.3 cuFFT

| 項目 | 仕様 |
|---|---|
| plan | `cufftPlan3d(nz, ny, nx, CUFFT_Z2Z)` |
| 精度 | complex float64 |
| 逆変換規格化 | CUDA kernelで全格子点数の逆数を乗算 |
| データ転送 | 初期upload、出力時download |

### 13.4 GPUタイミング

CUDA eventを使って次を測定します。

| 区分 | 内容 |
|---|---|
| `nonlinear_local` | 2回の局所半ステップ |
| `fft_forward_inverse` | 順cuFFTと逆cuFFT |
| `kinetic_spectral` | 運動項kernel |
| `allocation_and_other` | 逆FFT規格化、虚時間規格化など |
| `step_total` | GPU上の1ステップ全体 |

---

## 14. 複数GPU cuFFTMp仕様

### 14.1 基本要件

| 項目 | 要件 |
|---|---|
| OS | Linux |
| GPU | NVIDIA CUDA対応GPU |
| 配置 | 1 MPI rankにつき1 GPU |
| FFT | cuFFTMp |
| 通信 | MPI、NVSHMEM |
| コンパイラ | Fortran、MPI C++、nvcc |
| 推奨環境 | NVIDIA HPC SDKと整合するMPI/NVSHMEM |
| ペンシル版API | cuFFTMp 11.4.0（NVIDIA HPC SDK 25.3）以降 |

### 14.2 GPU選択

ノード内の `MPI_COMM_TYPE_SHARED` communicatorを作り、ノード内rank番号をCUDA device番号として使用します。

全rankから同じ `CUDA_VISIBLE_DEVICES` の一覧が見える必要があります。rankごとに1 GPUだけを可視化する運用とは組み合わせません。

### 14.3 データ配置

| 段階 | スラブ版 | ペンシル版 |
|---|---|---|
| 実空間 | `[local_z][ny][nx]` | `[local_z][ny][nx]` |
| 順FFT後 | `[nz][local_y][nx]` | `[local_spectral_z][local_y][nx]` |
| 逆FFT後 | 実空間zスラブ | 実空間zスラブ |

ペンシル版は`cufftMpMakePlanDecomposition`へ入力zスラブboxと出力xペンシルboxを渡します。
`p_y × p_z`はrank数の因数のうち正方形に近く、y・z方向に空領域を作らない組を選びます。
運動項、スペクトル微分、スペクトル診断は順FFT後の局所波数領域を直接処理し、
ホストへ集約しません。

### 14.4 規格化と診断

各GPUで局所和を計算し、MPI Allreduceで全体値へします。

| 量 | MPI演算 |
|---|---|
| ノルム | SUM |
| エネルギー | SUM |
| ARGLE平均密度 | SUM |
| ARGLE最小密度 | MIN |
| ARGLE最大変化率 | MAX |

### 14.5 エラー処理

一つのrankでcuFFTMp処理が失敗した場合、Fortran wrapperはエラー文字列を出力し、`MPI_Abort` で全rankを停止します。

### 14.6 タイミング

測定時は開始前にMPI Barrierを行い、区間ごとにCUDA deviceを同期して `MPI_Wtime` を記録します。測定自体が同期コストを増やすため、本番性能測定以外では無効化できます。

---

## 15. 入力ファイル仕様

入力はFortran namelist形式で、`&simulation` と `&gpe` の2節を使います。

実行ファイルの第1引数で入力ファイルを指定します。

```powershell
gp3d_sequential.exe input.nml
```

引数を省略すると `input.nml` を使用します。

### 15.1 `&simulation`

#### 識別・初期条件

| キー | 既定値 | 実際の用途 |
|---|---|---|
| `equation` | `GPE` | meta.jsonの識別。方程式選択には使わない |
| `case_name` | `gp3d` | 出力メタデータのcase名 |
| `input_file` | `input.nml` | 設定へ保持するが時間発展では未使用 |
| `initial_condition` | `gaussian` | 初期条件選択 |
| `restart_file` | 空 | 指定時は初期条件生成を上書きして再スタート |

#### 格子

| キー | 既定値 | 実際の用途 |
|---|---:|---|
| `nx, ny, nz` | `16, 16, 16` | 全体格子数 |
| `nghost` | `0` | SLF/metaへ記録。計算配列へghostを追加しない |
| `x_min, x_max` | `-6, 6` | x領域 |
| `y_min, y_max` | `-6, 6` | y領域 |
| `z_min, z_max` | `-6, 6` | z領域 |

#### 時間

| キー | 既定値 | 実際の用途 |
|---|---:|---|
| `dt` | `2.5e-4` | 常に使う固定時間刻み |
| `t_max` | `0` | `nsteps <= 0` の場合だけ `ceil(t_max/dt)` を導出 |
| `nsteps` | `20` | 今回追加で進めるステップ数 |
| `cfl` | `0` | 現行GPE計算では未使用 |
| `use_fixed_dt` | true | 現行コードは常に固定dt。フラグ自体は未使用 |

`nsteps > 0` の場合、`t_max` は終了判定に使用されません。

#### 出力

| キー | 既定値 | 実際の用途 |
|---|---|---|
| `output_frequency` | `10` | 正値なら、その倍数stepで出力。0以下なら途中出力なし |
| `output_dir` | `output` | SLFとmeta.jsonの出力先 |
| `output_format` | `slf` | metaへ記録。writerは現時点で常にSLF |
| `precision_name` | `float64` | metaへ記録。writerは現時点で常にfloat64 |
| `write_initial` | true | 時間ループ開始前の状態を出力 |
| `write_meta` | true | root rankがmeta.jsonを出力 |
| `timing_enabled` | false | Split-step区間計測 |

最終stepが `output_frequency` の倍数でない場合、最終状態は自動保存されません。

#### バックエンド・並列

| キー | 既定値 | 実際の用途 |
|---|---|---|
| `backend` | `dft` | 情報用。実際のバックエンドはビルドプロファイルで決まる |
| `use_mpi` | false | mainが実MPI状態で上書き |
| `use_openmp` | false | metaへ記録。計算ループのOpenMP化は未実装 |
| `use_cuda` | false | GPU mainが上書き |
| `cuda_device` | `0` | 単一GPUで使用。cuFFTMpはノード内rankから自動選択 |
| `rank` | `0` | mainが実行時rankで上書き |
| `nprocs` | `1` | mainが実行時process数で上書き |

### 15.2 `&gpe`

#### 方程式係数

| キー | 既定値 | 用途 |
|---|---:|---|
| `use_dimensionless_parameters` | false | `alpha, beta` モードを選択 |
| `alpha` | `0.05` | 無次元運動項係数 |
| `beta` | `40` | 無次元非線形係数 |
| `g` | `1` | 物理モードの相互作用係数 |
| `hbar` | `1` | 物理モード |
| `mass` | `1` | 物理モード |
| `norm` | `1` | 目標ノルム |
| `imaginary_time` | true | 通常Split-stepの実時間・虚時間選択 |

#### ポテンシャルと背景

| キー | 既定値 | 用途 |
|---|---:|---|
| `sigma0` | `1` | Gaussian幅 |
| `wx, wy, wz` | `1, 1, 1` | 調和ポテンシャル係数 |
| `mu` | `8` | Thomas–Fermi背景、物理モードARGLE |
| `healing_length` | `0.25` | 渦芯幅。無次元Taylor–Greenでは導出値を使用 |

#### 渦

| キー | 既定値 | 用途 |
|---|---:|---|
| `vortex_charge` | `1` | 単一渦・渦輪の巻き数 |
| `vortex_x0, vortex_y0` | `0, 0` | 直線渦中心 |
| `ring_radius` | `2` | 単一渦輪半径 |
| `ring_z0` | `0` | 単一渦輪z位置 |
| `tangle_nlines` | `8` | ランダム直線渦数 |
| `tangle_nrings` | `8` | ランダム渦輪数 |
| `ring_radius_min` | `0.8` | ランダム渦輪最小半径 |
| `ring_radius_max` | `1.8` | ランダム渦輪最大半径 |
| `phase_noise` | `0` | 位相ノイズ振幅 |
| `random_seed` | `12345` | 再現用seed |

#### Taylor–GreenとARGLE

| キー | 既定値 | 用途 |
|---|---:|---|
| `tg_velocity_amplitude` | `1` | Taylor–Green速度振幅 |
| `tg_winding` | `1` | 手動巻き数 |
| `tg_auto_winding` | true | 巻き数を自動導出 |
| `argle_enabled` | false | 実時間前のARGLE |
| `argle_write_seed` | true | 緩和前seedを保存 |
| `argle_steps` | `0` | 最大反復数 |
| `argle_output_every` | `100` | 進捗ログ間隔 |
| `argle_dtau` | `1e-3` | 擬時間刻み |
| `argle_tolerance` | `1e-8` | 最大変化率の収束閾値。0で無効 |

`argle_enabled` はコード上、初期条件名を検査しません。ただし右辺は常にTaylor–Green速度場を使用するため、標準運用では `quantum_taylor_green` と組み合わせます。

### 15.3 namelist読み込み上の注意

| 状況 | 現行動作 |
|---|---|
| 入力ファイルがない | 警告後に既定値で続行 |
| `&simulation` 読み込み失敗 | エラー停止 |
| `&gpe` 読み込み失敗 | 既定値を保持して続行 |

`&gpe` の綴り違いを見逃さないため、通常は `case.yaml` から自動生成し、直接編集しません。

---

## 16. `case.yaml`との契約

フレームワークでは `case.yaml` を正本とし、`ScriptLibrary` のinput adapterが `input.nml` を生成します。

### 16.1 主要対応

| `case.yaml` | `input.nml` |
|---|---|
| `case_id` | `simulation.case_name` |
| `flow.type` | `simulation.initial_condition` |
| `restart.file` | `simulation.restart_file` |
| `grid.*` | `simulation` の格子項目 |
| `time.*` | `simulation` の時間項目 |
| `output.*` | `simulation` の出力項目 |
| `solver.profile` | `solver_manifest.yaml` のprofile |
| `solver.processes` | 実行process数 |
| `physics.gpe.*` | `gpe` namelist |

### 16.2 GPEキーの互換名

| 旧名 | 正式名 |
|---|---|
| `interaction_strength` | `g` |
| `chemical_potential` | `mu` |

`density0` と `damping` は現行GP3D namelistに対応しないため、SetupCase validatorでは拒否します。

### 16.3 生成原則

1. `case.yaml` を編集する。
2. `case_input.py` または `run_case.py --prepare` で `input.nml` を再生成する。
3. `input.nml` を直接編集しない。
4. `solver.profile` と実行環境のprofileを一致させる。
5. 非MPI profileでは `solver.processes = 1` とする。

---

## 17. 出力仕様

### 17.1 ファイル

| 実行方式 | ファイル名 |
|---|---|
| 逐次CPU・単一GPU | `field_000000.slf` |
| CPU MPI・cuFFTMp | `field_000000_rank00000.slf` |
| メタデータ | `meta.json` |

SLFには `psi_real` と `psi_imag` の2変数だけを保存します。`density`、`phase`、`abs_psi` は後処理で導出します。

### 17.2 SLF1バイナリ構造

SLFはFortran stream unformattedで、次の順に書きます。

| 順序 | 型 | 内容 |
|---:|---|---|
| 1 | char 8 byte | `"SLF1"` とゼロ埋め |
| 2 | int32 | version = 1 |
| 3 | int32 | dtype_code = 2、float64 |
| 4 | int32 | ndim = 4 |
| 5 | int32 x 4 | local shape `(nx, ny, local_nz, nvar)` |
| 6 | int32 x 8 | step、rank、global shape、nghost、nprocs、k_start |
| 7 | float64 | time |
| 8 | float64 x 6 | セル中心座標の最小・最大 |
| 9 | int32 | nvar |
| 10 | char 32 x nvar | 変数名 |
| 11 | float64配列 | Fortran列優先field data |

データ配列は `(i, j, k, variable)` で、`i` が最速です。

writerはCPUネイティブの整数・浮動小数点表現を使用します。現行のWindows/Linux x86_64ではlittle-endianですが、異なるendian間の可搬性は保証していません。

### 17.3 SLF meta配列

| 要素 | 内容 |
|---:|---|
| 1 | step |
| 2 | rank。global出力は0 |
| 3 | global nx |
| 4 | global ny |
| 5 | global nz |
| 6 | nghost |
| 7 | source nprocs |
| 8 | local k_start |

### 17.4 `meta.json`

`meta.json` は次を記録します。

| 分類 | 内容 |
|---|---|
| case | equation、case_name |
| grid | 全体格子数、領域長、origin、spacing |
| data | precision、format、primary variables |
| parallel | MPI、OpenMP、CUDAフラグ |
| decomposition | `serial-global` または `z-slab-distributed-fft` |
| rank ranges | 各rankのi、j、k担当範囲 |

`origin` は最初のセル中心です。ParaView変換ツールはVTIのセル境界原点へ補正します。

### 17.5 出力タイミング

| 条件 | 出力 |
|---|---|
| `write_initial = true` | ARGLE後またはrestart読込後の開始状態 |
| `output_frequency > 0` | stepが頻度の倍数 |
| `output_frequency <= 0` | 時間ループ中の出力なし |
| 最終stepが頻度の倍数でない | 最終状態は保存しない |

---

## 18. 再スタート仕様

### 18.1 対応入力

| 入力 | 対応 |
|---|---|
| global SLF | 対応 |
| rank分割SLF family | 対応 |
| 保存時と異なるMPI分割数 | 対応 |
| float32 SLF | 非対応 |
| 格子数変更 | 非対応 |
| 領域境界変更 | 非対応 |

### 18.2 rank familyの指定

`restart_file` に次のどちらかを指定できます。

```text
output/field_000100_rank00000.slf
output/field_000100.slf
```

後者が存在せず、`field_000100_rank00000.slf` が存在する場合は自動解決します。

### 18.3 分割数変更

再スタートreaderは、保存SLFのz範囲と現在rankのz範囲の重なりだけを読みます。全source rankファイルを走査するため、保存時と再開時のMPI process数が異なっても、全体格子が同じなら復元できます。

### 18.4 検証条件

readerは次を検証します。

1. SLF1、version 1、float64、rank-4 dataである。
2. global格子数が一致する。
3. x-y local shapeが一致する。
4. z rangeが有効である。
5. `psi_real` と `psi_imag` がある。
6. セル中心領域境界が許容誤差内で一致する。
7. rank familyのstep、time、source nprocsが一致する。
8. 現rankのzスラブ全体がsource filesで覆われる。

### 18.5 再開後の時間

| 値 | 動作 |
|---|---|
| `start_step` | SLF headerのstep |
| `start_time` | SLF headerのtime |
| `nsteps` | startから追加で進める回数 |
| 出力step | `start_step + advance_step` |
| ARGLE | 実行しない |
| 初期条件生成 | 実行しない |
| 再規格化 | 読み込み直後には行わない |

---

## 19. 診断量とタイミング

### 19.1 標準出力

開始時と出力stepで次を表示します。

```text
# step norm energy
```

MPI・cuFFTMpではroot rankだけが表示します。

### 19.2 タイミング項目

| ラベル | 内容 |
|---|---|
| `nonlinear_local` | 2回の局所半ステップ |
| `fft_forward_inverse` | 1回の順FFTと1回の逆FFT |
| `kinetic_spectral` | Fourier空間運動項 |
| `allocation_and_other` | allocation、逆FFT規格化、虚時間規格化など |
| `step_total` | 1ステップ全体 |

出力列は次です。

| 列 | 内容 |
|---|---|
| `rank_mean_s` | rank間平均累積時間 |
| `rank_max_s` | 最も遅いrankの累積時間 |
| `max_s_per_step` | rank最大時間をstep数で割った値 |
| `mean_percent` | rank平均の `step_total` に対する割合 |

### 19.3 バックエンドごとの時計

| バックエンド | 計測 |
|---|---|
| CPU | `system_clock` のwall time |
| 単一GPU | CUDA event |
| cuFFTMp | CUDA同期とMPI Wtime |

GPU診断量を出力stepで計算する際の追加FFTは、Split-step timingには含みません。

---

## 20. ビルドシステム

### 20.1 CMake

`CMakeLists.txt` が正本です。Fortran 2008を要求します。

主要cache変数は次です。

| 変数 | 値 | 用途 |
|---|---|---|
| `USE_MPI` | ON/OFF | 実MPIを選択 |
| `FFT_BACKEND` | `dft`, `fftw` | CPU FFT |
| `FFT_DECOMPOSITION` | `slab`, `pencil` | CPU MPI／MPI cuFFTMp FFTの分割方式 |
| `GPU_BACKEND` | `none`, `cuda`, `cufftmp` | GPU実装 |
| `GP3D_CUDA_ARCHITECTURES` | 例 `86` | CUDA architecture |
| `FFTW_ROOT` | path | FFTW検索 |
| `GP3D_2DECOMP_ROOT` | path | ペンシル版の2DECOMP&FFT検索 |
| `CUFFTMP_ROOT` | path | cuFFTMp検索 |
| `NVSHMEM_ROOT` | path | NVSHMEM検索 |
| `CUFFTMP_API` | `auto`, `modern`, `legacy` | cuFFTMp API世代 |
| `BUILD_TESTING` | ON/OFF | CTest |

### 20.2 SolverLibrary profiles

`solver_manifest.yaml` がフレームワーク公開契約です。

| profile | 実行ファイル | MPI | CPU FFT | GPU |
|---|---|---:|---|---|
| `cpu_serial_dft` | `gp3d_sequential` | OFF | 参照DFT | なし |
| `cpu_serial_fftw` | `gp3d_sequential` | OFF | FFTW | なし |
| `cpu_mpi_dft` | `gp3d_mpi` | ON | 分散FFT + 局所DFT | なし |
| `cpu_mpi_fftw` | `gp3d_mpi` | ON | 分散FFT + 局所FFTW | なし |
| `cpu_mpi_pencil_fftw` | `gp3d_mpi` | ON | 2DECOMP&FFTペンシル + FFTW | なし |
| `cuda_single` | `gp3d_cuda` | OFF | テスト用DFTを同梱 | cuFFT |
| `cuda_mpi_cufftmp` | `gp3d_cufftmp` | ON | 比較テスト用分散DFTを同梱 | cuFFTMp |
| `cuda_mpi_cufftmp_pencil` | `gp3d_cufftmp` | ON | 比較テスト用分散DFTを同梱 | cuFFTMpペンシル |

通常はprofileを直接指定せず、実行環境設計書で分割方式を選びます。

```yaml
parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: false
  fft_decomposition: pencil  # slab または pencil
```

CPUでは`slab`が`cpu_mpi_fftw`、`pencil`が`cpu_mpi_pencil_fftw`を選びます。
`use_cuda: true`を併用すると、それぞれ`cuda_mpi_cufftmp`と
`cuda_mpi_cufftmp_pencil`を選びます。
この値はビルド構成なので、変更時は実行環境の再生成と再ビルドが必要です。

### 20.3 不正な組み合わせ

| 組み合わせ | 動作 |
|---|---|
| `GPU_BACKEND=cuda` かつ `USE_MPI=ON` | CMakeエラー |
| `GPU_BACKEND=cufftmp` かつ `USE_MPI=OFF` | CMakeエラー |
| Windowsで `GPU_BACKEND=cufftmp` | CMakeエラー |
| `FFT_DECOMPOSITION=pencil`かつMPI無効 | CMakeエラー |
| `FFT_DECOMPOSITION=pencil`かつ単一GPU CUDA | CMakeエラー |
| CPUで`FFT_DECOMPOSITION=pencil`かつ`FFT_BACKEND!=fftw` | CMakeエラー |
| cuFFTMpペンシルかつ`cufftMpMakePlanDecomposition`なし | CMakeエラー |
| 未知のFFT/GPU backend | CMakeエラー |

### 20.4 Windows MPI

CMakeの標準MPI検出に失敗した場合、Windowsでは次のMicrosoft MPIをfallbackとして使用します。

```text
C:\Program Files (x86)\Microsoft SDKs\MPI
```

link libraryは `msmpifec` と `msmpi` です。

CPUペンシル版は、同じFortranコンパイラ・MPI ABI・float64・FFTWでビルドした
2DECOMP&FFTを追加で必要とします。マシンプロファイルの
`libraries.decomp2d_root`へインストールprefixを指定し、BuildSolverが
`GP3D_2DECOMP_ROOT`へ変換します。

### 20.5 Windows CUDA

単一GPU版は次を組み合わせます。

| 要素 | 用途 |
|---|---|
| GNU Fortran | Fortran本体と最終実行ファイル |
| MSVC | nvccのhost C++ compiler |
| CUDA Toolkit | nvcc、cuFFT、CUDA runtime |
| bridge DLL | MSVC側CUDAコードとGNU FortranのC ABI接続 |

CUDA ToolkitがGPU driverの対応範囲より新しい場合、ビルドに成功しても実行時に「unsupported PTX toolchain」となることがあります。`GP3D_CUDA_ARCHITECTURES` とdriver/toolkitの整合を確認します。

### 20.6 Linux cuFFTMp

必要な主要ファイルは次です。

```text
cufftMp.h
libcufftMp.so
libnvshmem_host.so
```

スラブ版はHPC SDKの世代に応じ、`cufftMpMakePlan3d` またはlegacy
`cufftMpAttachComm`を選択します。ペンシル版は`cufftMpMakePlanDecomposition`を
使用するため、cuFFTMp 11.4.0（NVIDIA HPC SDK 25.3）以降が必要です。

---

## 21. フレームワークでのビルド・実行

### 21.1 設計情報

| 正本 | 役割 |
|---|---|
| `case.yaml` | 物理・格子・時間・出力条件 |
| `solver_manifest.yaml` | ソルバーprofileと必要ファイル |
| `build.local.yaml` など | 計算機・toolchain |
| `input.nml` | 自動生成されたソルバー入力 |
| `resolved_build.json` | 解決済みビルド計画 |

### 21.2 標準フロー

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --test
python .\tools\run_case.py --run
```

`--run` は既存実行ファイルを使い、毎回ビルドしません。`--all` はprepareを除き、構成に応じてビルド・テスト・実行をまとめるため、反復実行では使い分けます。

### 21.3 直接CMake

CPU逐次DFTの最小例:

```powershell
cmake --fresh `
  -S . `
  -B build\cpu-serial-dft `
  -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DUSE_MPI=OFF `
  -DFFT_BACKEND=dft `
  -DGPU_BACKEND=none `
  -DBUILD_TESTING=ON

cmake --build build\cpu-serial-dft --parallel 8
ctest --test-dir build\cpu-serial-dft --output-on-failure
```

通常はフレームワークの設計書を正本とし、直接CMakeは移植・デバッグ時に使用します。

---

## 22. テスト仕様

### 22.1 自動テスト

| テスト | 検証内容 |
|---|---|
| `restart_slf_roundtrip` | global/rank SLF書込・再読込、分割変更 |
| `distributed_fft_np2` | 2 rank Fourier modeと往復誤差 |
| `distributed_fft_np3` | 不均等3 rank分割 |
| `distributed_fft_pencil_np4` | 非立方・不均等格子、2×2プロセス格子のペンシルFFTと往復誤差 |
| `taylor_green_argle_serial` | 逐次ARGLE最小経路 |
| `taylor_green_argle_np2` | MPI ARGLE |
| `serial_restart_main` | mainでの再スタートstep/time |
| `mpi_restart_main_np2` | MPI mainでの再スタート |
| `cuda_splitstep_compare` | CPUと単一GPUの実時間・虚時間、norm、energy |
| `cuda_smoke` | GPU main、ARGLE、timing |
| `cufftmp_splitstep_np2` | CPU MPIと2 GPUのfield、norm、energy |
| `cufftmp_argle_smoke_np2` | cuFFTMp mainとARGLE |
| `cufftmp_pencil_splitstep_np4` | 非立方・不均等格子、cuFFTMp 2×2ペンシルとCPU MPIの比較 |

### 22.2 数値許容誤差

| 比較 | 許容値 |
|---|---:|
| 分散FFT Fourier mode | `1e-10` |
| 分散FFT roundtrip | `1e-10` |
| 単一GPU field/scalar | `5e-10` |
| cuFFTMp field/scalar | `2e-9` |
| restart field | `1e-13` |

### 22.3 2026-08-21の確認結果

ペンシル分割追加時に、現行ソースからクリーンビルドして次を確認しました。

| profile相当 | compiler/runtime | 結果 |
|---|---|---|
| `cpu_serial_dft` | GNU Fortran 16.1.0 | 3件中3件成功 |
| `cpu_mpi_dft` | GNU Fortran 16.1.0、Microsoft MPI | 5件中5件成功 |
| `cpu_serial_fftw` | FFTW未使用 | 今回未実施 |
| `cpu_mpi_fftw` | GNU Fortran 16.1.0、Microsoft MPI、FFTW | 6件中6件成功 |
| `cpu_mpi_pencil_fftw` | GNU Fortran 16.1.0、Microsoft MPI、2DECOMP&FFT 2.1.0、FFTW | 7件中7件成功 |
| `cuda_single` | GPU環境依存 | 今回未実施 |
| `cuda_mpi_cufftmp` | Linuxスパコンが必要 | 今回未実施 |
| `cuda_mpi_cufftmp_pencil` | Linuxスパコン、cuFFTMp 11.4.0以降が必要 | 今回未実施 |

Microsoft MPIの `mpif.h` からBOZ literal warning、link時に `.drectve` warningが出ましたが、対象テストはすべて成功しました。

---

## 23. エラー条件と不変条件

### 23.1 格子

| 条件 | 動作 |
|---|---|
| 格子数が0以下 | error stop |
| 領域長が0以下 | error stop |
| `nprocs > nz` | error stop |
| 分散FFTで `nprocs > ny` | error stop |
| ペンシルFFTで`nprocs > nz` | zスラブ保存配置のためerror stop |

### 23.2 物理係数

| 条件 | 動作 |
|---|---|
| `hbar <= 0` または `mass <= 0` | error stop |
| dimensionlessで `alpha <= 0` または `beta <= 0` | error stop |
| TF背景で `g <= 0` | error stop |
| healing lengthが0以下 | error stop |
| Gaussian幅が0以下 | error stop |
| target norm計算時にzero wavefunction | error stop |

### 23.3 ARGLE

| 条件 | 動作 |
|---|---|
| enabledでstepsが0以下 | error stop |
| `argle_dtau <= 0` | error stop |
| toleranceが負 | error stop |
| auto windingが0 | error stop |

### 23.4 GPU

| 条件 | 動作 |
|---|---|
| 単一GPUでlocal zが全体zと異なる | error stop |
| cuFFTMpでMPI無効 | error stop |
| node-local rank数がvisible GPU数を超える | backend error |
| Fortran分割とcuFFTMp natural layout不一致 | backend error |

---

## 24. 現行制約と既知の注意点

### 24.1 数値モデル

| 制約 | 影響 |
|---|---|
| 周期境界のみ | トラップ外縁や非周期問題には大領域・低境界密度が必要 |
| dealiasingなし | 強い非線形・高波数成分でaliasingの可能性 |
| 固定dt | CFLや誤差推定による自動制御なし |
| 最終step強制出力なし | output frequencyとの整合が必要 |
| 単一成分・零温度 | 有限温度、多成分現象は対象外 |

### 24.2 設定

| 設定 | 現状 |
|---|---|
| `cfl` | 未使用 |
| `use_fixed_dt` | 未使用。常に固定dt |
| `nghost` | 出力metadataのみ |
| `use_openmp` | metadataのみ |
| `output_format` | 実装はSLFのみ |
| `precision_name` | 実装はfloat64のみ |
| `backend` | 実行backend選択には使わない |
| `equation` | GPE以外を指定しても計算式はGPE |

### 24.3 初期条件

| 注意 | 内容 |
|---|---|
| `vortex_tangle` | 初期状態はz方向平行直線渦群 |
| `ring_tangle` | 三次元ランダム配向 |
| `uniform_vortex` | 背景密度は1固定 |
| dimensionless TF | 初期背景は `g`、時間発展は `beta` |
| ARGLE | フラグがtrueなら初期条件名に関係なくTG速度場を使用 |
| phase outside fluid | 低密度領域の位相は物理的意味が弱い |

### 24.4 性能

| 制約 | 内容 |
|---|---|
| CPU `psi_k` | 毎step allocate |
| MPI work | FFT・転置ごとにallocate |
| CPU MPI分割 | スラブ、またはFFT内部ペンシル |
| スラブMPI通信 | Split-step 1 stepにつきAlltoallv 4回 |
| ペンシルMPI通信 | Split-step 1 stepにつき外側再分配4回と2DECOMP内部転置4回 |
| ARGLE | 1反復につき三次元FFT約5回 |
| 出力 | host SLFであり、GPUでは出力時downloadが必要 |

### 24.5 移植

| 制約 | 内容 |
|---|---|
| SLF endian | native endian |
| Windows CUDA | MSVC host compilerが必要 |
| cuFFTMp | Linuxのみ |
| NVSHMEM | MPI ABIとbootstrap pluginの整合が必要 |
| GPE Git版 | FrameWorkモノレポ内で追跡する |

---

## 25. 拡張方法

### 25.1 初期条件を追加する

変更箇所:

1. `gp3d_model_config_t` に必要なパラメータを追加する。
2. `gp3d_input.f90` の `&gpe` に追加する。
3. `gp3d_initial_conditions.f90` に生成手続きを追加する。
4. `set_initial_condition` に名前を追加する。
5. `ScriptLibrary/SetupCase/gp3d_case.py` の許可名・GPEキーを更新する。
6. RunEnvironment input adapterとcase templateを更新する。
7. 小格子テストを追加する。
8. 本仕様書を更新する。

### 25.2 新しいFFTバックエンドを追加する

新実装は `module gp3d_fft` と共通APIを維持します。

```text
gp3d_fft_plan_t
gp3d_fft_init
gp3d_fft_forward
gp3d_fft_inverse
gp3d_fft_finalize
```

CMakeのsource選択、manifest component、profile、比較テストを追加します。

### 25.3 新しいGPUバックエンドを追加する

`module gp3d_gpu` の共通APIを維持し、GPU mainから差し替えられる設計にします。C ABIを使う場合は、Fortran側とC/CUDA側のdouble complex配列契約を維持します。

### 25.4 ポテンシャルを追加する

現行は `gp3d_set_harmonic_potential` を直接呼びます。複数ポテンシャルへ拡張する場合は、`src/potential` または同等のモジュールを設け、`case.yaml` の選択値から構築する形が望まれます。

### 25.5 出力形式を追加する

`output_format` を実際に有効化するには、次が必要です。

1. `gp3d_io` にwriterを追加する。
2. mainでformat dispatchする。
3. restart対応範囲を定義する。
4. meta契約を定義する。
5. 後処理ツールを追加する。
6. roundtrip testを追加する。

### 25.6 有限温度モデル

有限温度化はFFTだけの交換では完了しません。方程式、状態変数、ノイズ、散逸、診断、時間積分、入力、テストを独立モデルまたはsolver variantとして追加します。

---

## 26. 保守時の影響範囲

| 変更 | 主な影響先 |
|---|---|
| 物理係数 | types、input、solver、GPU bridges、tests |
| 時間積分式 | solver、CUDA kernel、cuFFTMp kernel、比較tests |
| grid layout | grid、FFT、GPU、IO、restart、postprocess |
| SLF header | IO、restart、Python converters、tests |
| case.yaml key | SetupCase、RunEnvironment adapter、schema、manual |
| profile | manifest、CMake、BuildSolver catalog、environment options |
| GPU memory layout | Fortran binding、C ABI、bridge、tests |
| ARGLE | CPU solver、2 GPU bridges、input、tests |

CPUとGPUで同じ数式を維持する変更では、必ずCPU referenceとの比較テストを更新します。

---

## 27. 再現性要件

本計算では次を保存します。

| 保存対象 | 理由 |
|---|---|
| `case.yaml` | 物理・数値条件の正本 |
| `input.nml` | 実際にソルバーへ渡した入力 |
| `solver_manifest.yaml` | profile構成 |
| `resolved_build.json` | 解決済みビルド条件 |
| compiler/CUDA/MPI version | 数値・ABI再現 |
| SolverLibrary Git commit | ソース版 |
| dirty status | 未コミット変更の有無 |
| `meta.json` | 格子・分割・変数 |
| SLF | 波動関数 |
| stdout/stderr | norm、energy、ARGLE、timing |
| ParaView state | 可視化条件 |

GP3DがGit未追跡の状態では、commit hashだけでソースを再現できません。正式計算前に追跡・コミットします。

---

## 28. 後処理

標準変換ツール:

```text
tools/slf_to_paraview_merged_cropghost.py
```

導出量:

```text
density = psi_real^2 + psi_imag^2
phase = atan2(psi_imag, psi_real)
abs_psi = sqrt(density)
```

global出力とrank分割出力の両方に対応します。

詳細は次を参照します。

```text
C:\Users\Owner\Documents\Codex\FrameWork\後処理_ParaView可視化手順書.md
```

---

## 29. ファイル責任一覧

| ファイル | 責任 |
|---|---|
| `CMakeLists.txt` | source/backend/target/test選択 |
| `solver_manifest.yaml` | フレームワーク公開profile |
| `src/common/gp3d_types.f90` | 型と既定値 |
| `src/grid/gp3d_grid.f90` | 周期格子と分割 |
| `src/init/gp3d_initial_conditions.f90` | 初期条件 |
| `src/solver/gp3d_solver.f90` | CPU数値計算 |
| `src/fft/gp3d_fft.f90` | 逐次参照DFT |
| `src/fft/gp3d_fft_fftw.f90` | 逐次FFTW |
| `src/fft/gp3d_fft_mpi.f90` | MPI分散FFT |
| `src/fft/gp3d_fft_pencil_2decomp.f90` | 2DECOMP&FFTペンシル分散FFT |
| `src/fft/gp3d_local_fft_dft.f90` | 局所参照DFT |
| `src/fft/gp3d_local_fft_fftw.f90` | 局所FFTW |
| `src/mpi/gp3d_mpi_stub.f90` | 逐次互換MPI |
| `src/mpi/gp3d_mpi_real.f90` | 実MPI |
| `src/gpu/gp3d_gpu_cuda.f90` | 単一GPU Fortran API |
| `src/gpu/gp3d_cuda_bridge.cu` | 単一GPU CUDA/cuFFT |
| `src/gpu/gp3d_gpu_cufftmp.f90` | cuFFTMp Fortran API |
| `src/gpu/gp3d_cufftmp_bridge.cu` | MPI/cuFFTMp |
| `src/io/gp3d_input.f90` | namelistと問題構築 |
| `src/io/gp3d_io.f90` | SLF/meta writer |
| `src/io/gp3d_restart.f90` | SLF restart reader |
| `src/main/main.f90` | CPU/MPI main |
| `src/main/main_cuda.f90` | 単一GPU main |
| `src/main/main_cufftmp.f90` | 複数GPU main |
| `tools/slf_to_paraview_merged_cropghost.py` | ParaView時系列変換 |
| `tests/*` | 回帰・比較テスト |

---

## 30. 要求仕様一覧

### 30.1 物理・数値

| ID | 要求 |
|---|---|
| PHY-001 | 単一成分三次元GPEをcomplex float64で解く |
| PHY-002 | 実時間と固定ノルム虚時間を選択できる |
| PHY-003 | 物理パラメータと無次元alpha/betaを選択できる |
| NUM-001 | 周期一様格子と擬スペクトル法を使用する |
| NUM-002 | Strang split-operator法を使用する |
| NUM-003 | 逆FFTは全格子点数で規格化する |
| NUM-004 | ARGLEは無拘束・半陰的更新とする |

### 30.2 並列・GPU

| ID | 要求 |
|---|---|
| PAR-001 | CPU MPI版は全体場をrootへ集約せず分散FFTする |
| PAR-002 | 実空間所有配置はzスラブとする |
| PAR-003 | CPU MPI FFTはスラブとペンシルを選択でき、既定値をスラブとする |
| GPU-001 | 単一GPU版は時間ループ中に `psi` をGPU常駐させる |
| GPU-002 | cuFFTMp版は1 rank 1 GPUとする |
| GPU-003 | cuFFTMp順FFT後の局所波数空間layoutを直接処理する |
| GPU-004 | cuFFTMp FFTはスラブとペンシルを選択でき、実空間とSLFはzスラブを維持する |

### 30.3 入出力

| ID | 要求 |
|---|---|
| IO-001 | 波動関数を `psi_real` と `psi_imag` でSLFへ保存する |
| IO-002 | MPI出力はrank別SLFとrank rangeを保存する |
| IO-003 | global/rank SLFから再スタートできる |
| IO-004 | 保存時と再開時のMPI分割数変更を許容する |
| IO-005 | case.yamlからinput.nmlを自動生成できる |

### 30.4 保守

| ID | 要求 |
|---|---|
| SW-001 | FFT実装は共通 `gp3d_fft` APIを維持する |
| SW-002 | MPI実装は共通 `gp3d_mpi` APIを維持する |
| SW-003 | GPU実装は共通 `gp3d_gpu` APIを維持する |
| SW-004 | profileは `solver_manifest.yaml` に定義する |
| SW-005 | CPU/GPU変更には比較テストを伴わせる |

---

## 31. 用語

| 用語 | 意味 |
|---|---|
| GPE | Gross–Pitaevskii equation |
| GP3D | 本三次元GPEソルバーパッケージ |
| SLF | Solver Library Format |
| Strang分割 | 半step、全step、半stepの対称演算子分割 |
| 擬スペクトル法 | 微分演算をFourier空間で処理する方法 |
| zスラブ | z方向の一部を各rankが所有する分割 |
| yスラブ | z方向FFTのためy方向を分割した一時配置 |
| Xペンシル | x方向を全保持し、y・z方向を2次元プロセス格子で分割する配置 |
| Zペンシル | z方向を全保持し、x・y方向を2次元プロセス格子で分割する配置 |
| 2DECOMP&FFT | ペンシル分割の転置と分散FFTを提供するMPIライブラリ |
| ARGLE | Advective Real Ginzburg–Landau Equationによる緩和 |
| cuFFT | NVIDIAの単一GPU FFT |
| cuFFTMp | NVIDIAのMPI対応分散GPU FFT |
| NVSHMEM | GPU間通信・対称メモリ基盤 |

---

## 32. 更新規則

次の変更では、本書を同じ変更単位で更新します。

1. 方程式または時間積分式を変更した。
2. input key、既定値、case.yaml mappingを変更した。
3. 初期条件を追加・削除した。
4. SLFまたはmeta.jsonの形式を変更した。
5. backend、profile、実行ファイルを変更した。
6. GPU memory layoutまたはMPI分割を変更した。
7. テスト基準を変更した。
8. 既知の制約を解消または追加した。

Markdown版を内容の正本とし、Word版はMarkdown版から同時に再生成します。
