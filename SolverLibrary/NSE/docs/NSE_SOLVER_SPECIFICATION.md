# NSEソルバー総合仕様書

更新日: 2026-08-25
仕様区分: 現行実装準拠（as implemented）
対象: `FrameWork/SolverLibrary/NSE` および `ScriptLibrary/RunEnvironment`

## 1. 目的と適用範囲

本書は、FrameWorkに含まれる三次元圧縮性Navier–Stokes（NSE）ソルバーの
支配方程式、無次元化、数値解法、初期条件、境界条件、Forcing、並列化、
入出力および既知の制約を一つにまとめた総合仕様書である。

本書では理想的な設計ではなく、2026-08-25時点のローカルソースが実際に行う
計算を仕様とする。既存文書と実装が一致しない箇所は「現行実装上の注意事項」へ
明記する。

## 2. ソルバー概要

| 項目 | 現行仕様 |
|---|---|
| 方程式 | 三次元圧縮性Euler／Navier–Stokes方程式 |
| 気体モデル | 比熱比一定の熱量的完全気体 |
| 保存変数 | 5変数 `[rho, rho*u, rho*v, rho*w, rho*E]` |
| 格子 | 一様直交Cartesian、セル中心有限体積法 |
| 境界条件 | x、y、z全方向の周期境界のみ |
| 対流流束 | KEEP2、KEEP6、WENO5-Z/Roe、KEEP/WENOハイブリッド |
| 粘性項 | 無効、または一定輸送係数の6次精度中心差分 |
| 時間積分 | 3段3次SSPRK（SSPRK3） |
| 初期条件 | Taylor–Green渦、スペクトルHIT、保存済み乱流場のx方向配置 |
| 外力 | なし、またはPetersen–Livescu線形Forcing |
| CPU並列 | MPIによるy-z分割＋OpenMP |
| GPU | 単一NVIDIA GPU、またはMPI＋複数NVIDIA GPU、CUDA、float64 |
| 出力 | rank別SLF保存変数、JSONメタデータ |
| 再スタート | 時刻・stepを継続する再スタートは未実装。保存済み乱流場の初期条件読込みは対応 |

## 3. 支配方程式

### 3.1 保存形

無次元保存変数を

```text
U = [rho, rho*u, rho*v, rho*w, rho*E]^T
```

とし、ソルバーは次の保存形を解く。

```text
partial(U)/partial(t) + div(Fc) = div(Fv) + Sforcing
```

法線方向をxとしたときの非粘性流束は

```text
Fc_x = [
  rho*u,
  rho*u^2 + p,
  rho*u*v,
  rho*u*w,
  u*(rho*E + p)
]^T
```

であり、y、z方向には速度・運動量成分を回転して同じ式を用いる。

### 3.2 状態方程式

```text
p = (gamma - 1) * [rho*E - 0.5*rho*(u^2 + v^2 + w^2)]
a = sqrt(gamma*p/rho)
theta = p/rho
```

`theta`は粘性熱伝導項で使用する無次元温度に相当する量である。ソルバーは
気体定数を独立な入力として持たない。

### 3.3 粘性応力と熱伝導

`viscous_scheme: central6`では、一定粘性係数、Stokesの仮定、Fourier熱伝導を
用いる。粘性応力の無次元形を

```text
tau_ij = partial(u_i)/partial(x_j)
       + partial(u_j)/partial(x_i)
       - (2/3)*delta_ij*div(u)
```

とすると、運動量と全エネルギーの粘性右辺は

```text
momentum_i:
  (1/Re) * partial(tau_ij)/partial(x_j)

energy:
  (1/Re) * [
    u_i*partial(tau_ij)/partial(x_j)
    + tau_ij*partial(u_i)/partial(x_j)
    + gamma/((gamma-1)*Pr) * laplacian(theta)
  ]
```

である。密度方程式には粘性項を加えない。輸送係数の温度依存性、Sutherland則、
バルク粘性は実装していない。

### 3.4 Forcing源項

Petersen–Livescu Forcingを有効にした場合の保存変数源項は、現行実装では

```text
Sforcing = [0, f_x, f_y, f_z, 0]^T
```

である。運動量RHSにはForcingを加えるが、保存全エネルギーRHSへ直接の
Forcing項は加えない。これは現行実装の明示的なモデル選択である。

## 4. 無次元化

### 4.1 採用しているスケーリング

本ソルバーは、速度を基準音速で無次元化する音響スケーリングを採用している。
基準長さ、基準密度、基準音速をそれぞれ `L_ref`、`rho_ref`、`a_ref` とすると、

| 物理量 | 無次元量 |
|---|---|
| 座標 | `x* = x/L_ref` |
| 時間 | `t* = t*a_ref/L_ref` |
| 速度 | `u* = u/a_ref` |
| 密度 | `rho* = rho/rho_ref` |
| 圧力 | `p* = p/(rho_ref*a_ref^2)` |
| 比全エネルギー | `E* = E/a_ref^2` |
| 温度相当量 | `theta* = p*/rho* = R*T/a_ref^2` |

このため、支配方程式には `1/M^2` のようなMach数係数が現れない。
`mach_number`は方程式係数ではなく、初期速度の大きさを与える入力である。

### 4.2 Reynolds数とPrandtl数

入力Reynolds数は

```text
Re = rho_ref*a_ref*L_ref/mu_ref
Pr = mu_ref*cp/kappa_ref
```

に対応し、動粘性係数は局所密度を用いて

```text
nu* = 1/(Re*rho*)
```

となる。したがって `reynolds_number` は基準音速に基づくReynolds数であり、
初期流速に基づくReynolds数とは一般に一致しない。

例えば `rho0=1`、代表長さを1、Taylor–Green速度振幅を `U0*=M` とすると、

```text
Re_TGV = U0*/nu* = M * Re_input
```

である。目標とする渦Reynolds数が `Re_TGV` の場合は、音響スケーリング上の入力を

```text
Re_input = Re_TGV/M
```

とする必要がある。

### 4.3 基準静止状態

音速を1にする基準状態は

```text
p0 = rho0/gamma
a0 = sqrt(gamma*p0/rho0) = 1
```

である。既定値 `rho0=1`、`gamma=1.4` では `p0=1/gamma` となる。

### 4.4 有次元量への戻し方

解析結果を有次元化するときは、計算前に定めた `L_ref`、`rho_ref`、`a_ref` を用いる。

```text
x_dim   = L_ref * x_solver
t_dim   = (L_ref/a_ref) * t_solver
u_dim   = a_ref * u_solver
rho_dim = rho_ref * rho_solver
p_dim   = rho_ref*a_ref^2 * p_solver
E_dim   = a_ref^2 * E_solver
```

ソルバーはSI単位、基準長さ、基準音速、基準密度を入力として保持しない。
入力YAMLの座標、時間、速度、Forcing、散逸率などは、利用者があらかじめ
無次元化した値である。

### 4.5 HIT目標値からの輸送係数導出

HITでは `turbulent_reynolds_number` をTaylorマイクロスケールReynolds数
`Re_lambda` と解釈する。目標乱流Mach数を `M_t`、スペクトル代表長さを `L` とし、

```text
u'        = M_t/sqrt(3)                  # 1成分RMS
Re_L      = 3*Re_lambda^2/20
U_L       = sqrt(3/2)*u'
nu        = U_L*L/Re_L
Re_solver = 1/nu
lambda    = L*sqrt(10/Re_L)
eta       = L*Re_L^(-3/4)
```

を用いる。HITを選択すると、`physics.nse.mach_number` と
`physics.nse.reynolds_number` はこれらの導出値で上書きされる。

この導出は `rho0=1` を前提にすると粘性実装 `nu=1/(Re*rho)` と一致する。
`rho0`を1以外にする場合、指定した `Re_lambda` と実際の粘性係数の対応は変わる。

## 5. 保存変数と数値的下限値

保存配列の成分は次のとおりである。

| 添字 | 変数 |
|---:|---|
| 1 | `rho` |
| 2 | `rho*u` |
| 3 | `rho*v` |
| 4 | `rho*w` |
| 5 | `rho*E` |

基本変数を復元するときは

```text
rho_eval = max(rho, small_rho)
p_eval   = max((gamma-1)*(rhoE-kinetic_energy), small_p)
```

を使用する。ただし、この処理は保存変数そのものを正値へ修正するものではない。
WENO5-Z/Roeにもpositivity-preserving limiterは実装されていない。したがって
`small_rho`と`small_p`は評価時の除算・平方根保護であり、発散防止を保証しない。

## 6. 格子と有限体積離散

### 6.1 格子

- 一様直交Cartesian格子
- セル中心保存変数
- `dx=(x_max-x_min)/nx`、y、zも同様
- 非一様格子、曲線座標、埋め込み境界は未実装
- 現行境界・全数値方式はghostセルをちょうど3層要求する

セル `i,j,k` の半離散式は

```text
dU_ijk/dt =
  -(Ax*F_{i+1/2}-Ax*F_{i-1/2})/V
  -(Ay*G_{j+1/2}-Ay*G_{j-1/2})/V
  -(Az*H_{k+1/2}-Az*H_{k-1/2})/V
  + Rviscous + Rforcing
```

である。共有面には一つの面流束を定義し、隣接セルが逆符号で使う。

### 6.2 KEEP2／KEEP6

対称二点KEEP流束を保存形の面流束へ合成する。

| 方式 | 正側一次微分係数 | 形式精度 |
|---|---|---:|
| `keep2` | `[1/2]` | 2次 |
| `keep6` | `[3/4, -3/20, 1/60]` | 6次 |

KEEP6は距離1～3の点対を使用する。非ハイブリッドで `keep2` または `keep6` を
選択した場合、WENO計算は行わない。

### 6.3 WENO5-Z／Roe

`weno5z_roe`は次の処理を行う。

1. 面中央の左右状態からRoe平均と固有ベクトルを計算する。
2. 保存変数を局所法線方向へ回転する。
3. Roe左固有ベクトルで5点stencilを特性空間へ射影する。
4. WENO5-Zで左右状態を再構築する。
5. Harten–Hyman型entropy fix付きRoe流束を計算する。
6. 流束を全体座標へ戻す。

WENO-Zの既定epsilonは `1.0e-20`、最適重みは左再構築で
`[0.1, 0.6, 0.3]`、右再構築でその鏡像である。形式精度は滑らかな領域で5次。

### 6.4 KEEP/WENOハイブリッド

`convective_scheme: hybrid`では、面ごとに

```text
F_hybrid = (1-alpha)*F_smooth + alpha*F_shock
```

を計算する。既定構成は `F_smooth=KEEP6`、`F_shock=WENO5-Z/Roe` である。
構成流束には `keep2`、`keep6`、`weno5z_roe`を指定できる。

Ducros-pressureセンサーは、速度発散、渦度および方向別圧力曲率から0～1の値を
作る。セルセンサーの概略は

```text
D = compression^2/(div(u)^2 + |omega|^2 + epsilon)
P = |p_plus - 2*p_center + p_minus|
    /(p_plus + 2*p_center + p_minus + epsilon)
S_cell = clamp(D*P, 0, 1)
S_face = max(S_left, S_right)
```

である。実装では極小勾配によるunderflowを避けるため、発散と渦度を最大絶対値で
正規化してから二乗する。二つの閾値間は

```text
x = clamp((S_face-onset)/(full-onset), 0, 1)
alpha = x^2*(3-2*x)
```

で連続化する。

### 6.5 粘性項 `central6`

速度3成分と `theta=p/rho` をghostセルを含めて復元し、一次、二次、混合微分を
6次精度中心差分で評価する。一次微分は

```text
(-f_{i-3}+9f_{i-2}-45f_{i-1}+45f_{i+1}-9f_{i+2}+f_{i+3})/(60*dx)
```

二次微分は

```text
(2f_{i-3}-27f_{i-2}+270f_{i-1}-490f_i
 +270f_{i+1}-27f_{i+2}+2f_{i+3})/(180*dx^2)
```

である。`viscous_scheme: none`では粘性・熱伝導項を加えない。

## 7. 時間積分と時間刻み

### 7.1 SSPRK3

時間積分はSSPRK3のみである。

```text
U1     = Un + dt*L(Un)
U2     = 3/4*Un + 1/4*(U1 + dt*L(U1))
U(n+1) = 1/3*Un + 2/3*(U2 + dt*L(U2))
```

境界条件、対流項、粘性項、Forcingは各段のRHS評価で再計算される。
Petersen–Livescu Forcingも1時間ステップにつき3回評価される。

### 7.2 固定時間刻み

`time.use_fixed_dt: true`では `time.dt` を使用する。ただし最終ステップで
`t+dt>t_max`となる場合は、`dt=t_max-t`へ切り詰める。

### 7.3 自動時間刻み

`time.use_fixed_dt: false`では、全領域の最大特性速度

```text
lambda_max = max(|u|+a, |v|+a, |w|+a)
dt_conv = CFL * min(dx,dy,dz)/lambda_max
```

をMPI全rankで評価する。粘性が有効な場合はさらに

```text
alpha_max = max(4/3, gamma/Pr)/(Re*rho_min)
dt_diff = 2 /
  [(272/45)*alpha_max*(1/dx^2 + 1/dy^2 + 1/dz^2)]
dt = min(dt_conv, dt_diff)
```

とする。

### 7.4 終了条件

時間ループは

```text
t < t_max かつ step < nsteps
```

の間だけ継続する。どちらかへ到達すると正常終了し、

```text
NSE calculation completed successfully: step=<step>, time=<time>
```

を表示する。

## 8. 境界条件

現行実装は三方向周期境界のみである。

- x方向は各MPI rankが全x範囲を保持するためローカルコピー
- y、z方向はMPI halo交換と周期端rank間通信
- 6次混合微分に必要な辺・角ghostを埋めるため、段階的にhaloを伝播
- ghostセル数はちょうど3
- 保存変数数はちょうど5

壁面、流入流出、対称、非反射境界は未実装である。

## 9. 初期条件

### 9.1 Taylor–Green渦

現行式は

```text
rho = rho0
u = M*sin(x)*cos(y)*cos(z)
v = -M*cos(x)*sin(y)*cos(z)
w = 0
p = 1/gamma
  + rho0*M^2/16 * (cos(2x)+cos(2y))*(cos(2z)+2)
rhoE = p/(gamma-1) + 0.5*rho*(u^2+v^2+w^2)
```

である。標準的には `[0,2*pi]^3` を用いる。

重要な実装仕様:

- 速度振幅は `physics.nse.mach_number` であり、`flow.taylor_green.amplitude`は未使用。
- 波数は1へ固定され、`mode_x`、`mode_y`、`mode_z`は未使用。
- 平均圧力は常に `1/gamma`。したがって音速基準との厳密な整合は `rho0=1` のときだけ成立する。

### 9.2 スペクトルHIT

HIT初期化は周期領域上で次を行う。

1. JohnsenまたはPope型エネルギースペクトルからFourier係数振幅を作る。
2. seedと全体波数番号からMPI分割数に依存しない位相を生成する。
3. 波数に直交する二つの基底でsolenoidal速度を直接生成する。
4. Hermitian対称性を課し、逆FFT後の実数性を保証する。
5. `dealias_fraction`を各軸Nyquist波数へ適用する。
6. 必要なら低波数shellのReynolds応力を等方化する。
7. 1成分RMSを `M_t/sqrt(3)` へ全体正規化する。
8. 速度勾配から周期Poisson方程式を解いて初期圧力変動を求める。
9. 密度、速度、圧力から保存変数を構築する。

Johnsenスペクトル形状は

```text
r = k/k_peak
Eshape(k) = r^4*exp(-2*r^2)
k_peak = 2*length_scale_ratio/integral_length
```

Popeスペクトル形状は

```text
Eshape(k) = C*K^(-5/3)*f_L(kL)*f_eta(k*eta)
```

であり、定数と指数はcaseで指定できる。最終的な速度RMSはスペクトル形状とは
独立に目標値へ再正規化される。

初期圧力は非圧縮速度場を仮定した

```text
laplacian(p) = -rho0 *
  [partial(u_i)/partial(x_j)]*[partial(u_j)/partial(x_i)]
```

をスペクトル法で解く。ゼロ波数は

```text
mean(p) = rho0/gamma
```

に固定する。初期圧力最小値が `small_p` 以下の場合は計算を開始せず停止する。

### 9.3 HIT FFTバックエンド

| 実行方式 | 初期化FFT |
|---|---|
| CPU MPI/OpenMP | 2DECOMP&FFT + FFTW3 |
| 単一GPU | cuFFT |
| MPI + CUDA | cuFFTMp（`cuda_mpi_cufftmp`） |

生成環境では、`flow.type: hit`を検出すると対応profileを自動選択する。

### 9.4 保存済み乱流場

`flow.type: imported_turbulence`では、付属ツールでghostを除去して単一ファイルへ
集約した保存量SLFを読み込む。`embed`は長いx領域の一部へ乱流ブロックを置き、
左右端を指定セル数のraised-cosineで背景保存量と混合する。`tile`は元のxセル列を
周期反復する。読み込んだ速度へ一定の`velocity_offset`を加える場合は、圧力を保って
運動量と全エネルギーを再構築する。

CPU版は全x範囲を各rankが保持する既存分割を利用し、自rankのy-z範囲だけをSLFから
直接読む。このため初期データ作成時と本計算時のMPI並列数は独立である。
CUDA版はCPU上で同じ可搬SLFから各rankの初期場を構築した後、完成した保存変数場を
そのrankのGPUへ転送する。したがってCPU MPI、単一GPU、MPI＋CUDAのどの出力からでも、
CPU MPI、単一GPU、MPI＋CUDAのいずれへ読み込める。元計算のバックエンドやMPI
プロセス数を読込み先へ引き継がない。

| 元計算 | 読込み先 | データ受渡し |
|---|---|---|
| CPU MPI | CPU MPI | rank別出力→可搬SLF→読込み先MPI分割で局所読込み |
| CPU MPI | 単一GPU | rank別出力→可搬SLF→CPU初期化→GPU転送 |
| 単一GPU | CPU MPI | CUDA単一rank出力→可搬SLF→CPU MPI局所読込み |
| 単一GPU | 単一GPU | CUDA単一rank出力→可搬SLF→CPU初期化→GPU転送 |
| CPU MPI／単一GPU／MPI＋CUDA | MPI＋CUDA | 可搬SLF→読込み先MPI分割で局所読込み→各GPUへ転送 |

格子補間は行わず、`dx,dy,dz`、`ny,nz`、y-z領域の一致を要求する。`tile`ではさらに
対象`nx`が元データ`nx`の整数倍でなければならない。全セルについて有限値、密度、
圧力を検査してから計算を開始する。変換と入力の詳細は
[`NSE_IMPORTED_TURBULENCE.md`](NSE_IMPORTED_TURBULENCE.md)に定義する。

## 10. Petersen–Livescu Forcing

密度重み付き速度

```text
w_i = sqrt(rho)*u_i
```

をFFTし、Helmholtz分解でsolenoidal成分 `w_s` とdilatational成分 `w_d` に分ける。
`low_wavenumber`では `0<|k|<k_cutoff` だけを保持し、`full_spectrum`では
ゼロ波数以外を保持する。

目標注入率 `epsilon` と比 `r=epsilon_d/epsilon_s` から

```text
epsilon_s = epsilon/(1+r)
epsilon_d = epsilon-epsilon_s
PD = <p*div(u)>
c_s = epsilon_s/<w_s.w>
c_d = (epsilon_d-PD)/<w_d.w>
```

を求め、運動量へ

```text
f = sqrt(rho)*(c_s*w_s + c_d*w_d)
```

を加える。分母下限と任意の係数上限を持つ。周期境界が必須である。

| 実行方式 | Forcing FFT |
|---|---|
| CPU MPI/OpenMP | 2DECOMP&FFT + FFTW3 |
| 単一GPU | cuFFT |
| MPI + CUDA | cuFFTMp（`cuda_mpi_cufftmp`） |

## 11. 並列化と実行profile

### 11.1 CPU MPI/OpenMP

CPU版はx方向を各rankが全域保持し、y-z平面を二次元MPI分割する。

```text
local domain = full x × local y block × local z block
```

MPIプロセス格子は因数対のうち和が小さい組合せを選ぶ。2または3 rankでは
1×P、4 rank以上では可能な限り正方形に近いy-z分割を選ぶ。MPI実行は2 rank以上を
要求し、`run_case.py`も1 rankを拒否する。

OpenMPはセルループ、流束計算、SSPRK更新などへ適用される。MPI rank数と
rankあたりOpenMPスレッド数は実行時に変更できる。

### 11.2 単一GPU CUDA

- GPU数は1
- MPI/OpenMPは無効
- 保存変数、RHS、SSPRK作業配列をGPUへ常駐
- 周期ghost、CFL、対流、粘性、Forcing、SSPRKをGPUで実行
- 出力時のみ保存変数をCPUへ戻す
- float64固定

### 11.3 MPI＋CUDAマルチGPU

- x全域×局所y×局所zを1 MPI rankのGPU 1台へ保持
- SSPRK3の各段でy面、続いてz面をhost staging方式で交換
- y交換後にz交換することで面だけでなく辺・角ghostも完成
- CFL時間刻みは全rankの最小値へ同期
- GPU番号はノード内rankから自動選択。1 GPUだけがrankへ公開された場合はdevice 0を使用
- 各局所y/zブロックは`nghost`セル以上
- CUDA-aware MPIは不要
- `cuda_mpi_cufftmp`でcuFFTMpによる分散HIT初期化とFFT forcingを実装

### 11.4 profile

| profile | 用途 | 主な能力 |
|---|---|---|
| `cpu_mpi` | 通常CPU計算 | MPI/OpenMP、Taylor–Green |
| `cpu_mpi_2decomp_fftw` | スペクトルCPU計算 | HIT、Petersen–Livescu Forcing |
| `cuda_single` | 単一GPU | CUDA、cuFFT HIT/Forcing |
| `cuda_mpi` | 複数GPU | MPI＋CUDA、Taylor–Green、保存済み乱流場 |
| `cuda_mpi_cufftmp` | 複数GPU分散FFT | MPI＋CUDA、cuFFTMp HIT/Forcing |

NSEのCPU/MPI生成環境には互換profileを同梱し、`case.yaml`の要求機能を満たす
最小profileを `run_case.py` が選択する。`solver.profile`を明示した場合は固定指定となり、
互換profileへ暗黙に変更しない。

## 12. case.yaml入力仕様

### 12.1 主要項目

| YAMLパス | 意味・制約 |
|---|---|
| `physics.nse.nv` | 5固定 |
| `physics.nse.gamma` | 比熱比。通常は1より大きい値 |
| `physics.nse.small_rho` | 密度評価下限。正値 |
| `physics.nse.small_p` | 圧力評価下限。正値 |
| `physics.nse.rho0` | 無次元初期密度。現行初期条件では1を推奨 |
| `physics.nse.mach_number` | 非HIT初期速度振幅。方程式係数ではない |
| `physics.nse.reynolds_number` | 音響スケーリングReynolds数 |
| `physics.nse.prandtl_number` | Prandtl数 |
| `flow.type` | `taylor_green`, `hit`, `imported_turbulence` |
| `flow.imported_turbulence.file` | ghostなしの可搬NSE SLF |
| `flow.imported_turbulence.mode` | `embed`または`tile` |
| `flow.imported_turbulence.x_start` | 元乱流セル列を開始するxセル境界座標 |
| `flow.imported_turbulence.blend_cells` | `embed`両端の混合セル数 |
| `flow.imported_turbulence.velocity_offset` | 読込み速度へ加える一定速度3成分 |
| `flow.imported_turbulence.background` | `embed`外側の密度、速度、圧力 |
| `flow.hit.turbulent_mach_number` | HIT目標乱流Mach数 `M_t` |
| `flow.hit.turbulent_reynolds_number` | HIT目標 `Re_lambda` |
| `flow.hit.random_seed` | MPI分割数に依存しない乱数seed |
| `flow.hit.dealias_fraction` | `(0,1]`。既定は2/3 |
| `flow.hit.isotropy_mode` | `none` または `projected_shell` |
| `flow.hit.spectrum.type` | `johnsen` または `pope` |
| `forcing.type` | `none` または `petersen_livescu` |
| `forcing.petersen_livescu.spectrum` | `full_spectrum` または `low_wavenumber` |
| `forcing.petersen_livescu.fft_backend` | 通常は `auto` |
| `forcing.petersen_livescu.k_cutoff` | 低波数Forcingの上限波数 |
| `forcing.petersen_livescu.target_dissipation` | 正の無次元目標注入率 |
| `forcing.petersen_livescu.dilatational_ratio` | 非負の `epsilon_d/epsilon_s` |
| `grid.nx,ny,nz` | 物理セル数 |
| `grid.nghost` | 3固定 |
| `grid.*_min,*_max` | 無次元周期領域 |
| `time.cfl` | 自動時間刻みCFL |
| `time.dt` | 固定刻み、または初期値 |
| `time.t_max` | 無次元終了時刻 |
| `time.nsteps` | 最大ステップ数 |
| `time.use_fixed_dt` | 固定／自動時間刻み |
| `time.output_frequency` | ステップ出力間隔 |
| `numerics.convective_scheme` | `keep2`, `keep6`, `weno5z_roe`, `hybrid` |
| `numerics.hybrid.smooth_scheme` | ハイブリッド滑らか側の構成流束 |
| `numerics.hybrid.shock_scheme` | ハイブリッド衝撃波側の構成流束 |
| `numerics.hybrid.sensor_onset` | 混合開始閾値。0以上 |
| `numerics.hybrid.sensor_full` | 完全切替閾値。onsetより大きい値 |
| `numerics.viscous_scheme` | `none` または `central6` |
| `numerics.boundary_condition` | `periodic`のみ |
| `numerics.time_integrator` | `ssprk3`のみ |
| `solver.mpi_processes` | MPI版は2以上。`cuda_mpi`では通常GPU総数と同じ |
| `solver.use_openmp` | case単位のOpenMP使用可否 |
| `solver.omp_threads` | rankあたりスレッド数 |
| `solver.cuda_device` | 単一GPU番号、または各ノードのGPU番号割当ての開始値 |
| `output.directory` | 出力先。caseディレクトリからの相対指定を推奨 |
| `output.write_initial` | step 0の保存可否 |
| `output.write_meta` | メタデータ指定。CPU版の注意は15.3節参照 |

### 12.2 最小例

```yaml
physics:
  model: nse
  nse:
    nv: 5
    gamma: 1.4
    small_rho: 1.0e-12
    small_p: 1.0e-12
    rho0: 1.0
    mach_number: 0.5
    reynolds_number: 200.0
    prandtl_number: 0.72

flow:
  type: taylor_green

grid:
  nx: 64
  ny: 64
  nz: 64
  nghost: 3
  x_min: 0.0
  x_max: 6.283185307179586
  y_min: 0.0
  y_max: 6.283185307179586
  z_min: 0.0
  z_max: 6.283185307179586

time:
  cfl: 0.5
  dt: 1.0e-4
  t_max: 0.1
  nsteps: 1000
  output_frequency: 100
  use_fixed_dt: false

numerics:
  convective_scheme: keep6
  viscous_scheme: central6
  boundary_condition: periodic
  time_integrator: ssprk3

solver:
  use_openmp: true
  mpi_processes: 4
  omp_threads: 1
```

`case.yaml`が正本であり、生成されたFortran namelist `input.dat`を直接編集しない。

## 13. 出力と後処理

### 13.1 SLF出力

各rankは

```text
field_<step:6桁>_rank<rank:5桁>.slf
```

へ保存変数をfloat64で出力する。SLFにはstep、rank、全体格子数、ghost数、時刻、
領域、変数名を含む。CPU版のrank別配列にはghost領域も含まれ、後処理ツールが
`meta.json`のrank範囲を用いて物理領域をcrop・統合する。

保存される一次変数は

```text
rho, rho_u, rho_v, rho_w, rho_E
```

のみである。ParaView変換時に `u,v,w,p` を導出できる。

### 13.2 出力時刻

- `write_initial: true`ならstep 0を保存
- step>0は `step mod output_frequency == 0` のとき保存
- 最終stepを無条件で保存する処理はない

したがって最終状態を必ず保存したい場合は、終了stepを
`output_frequency`の倍数にする必要がある。

### 13.3 時間計測

CPU版は各ステップについて

```text
step, time, dt, step_wall_seconds, total_wall_seconds
```

をroot rankから出力し、`MPI_Wtime()`で計測する。CUDA版も同じ内容を
`system_clock`で出力する。時間にはステップ内の出力処理も含まれる。

### 13.4 乱流統計

後処理では、平均値を差し引いた速度変動からRMS、乱流運動エネルギー、
Reynolds応力、等方性誤差、積分スケール、散逸率、Taylor長、Kolmogorov長、
`Re_L`、`Re_lambda`、乱流Mach数を計算できる。

## 14. 検証項目

現行テスト群は主に次を確認する。

- KEEP2／KEEP6の形式精度と周期保存性
- WENO5-Z再構築、Roe流束、entropy fix、衝撃波管
- ハイブリッドのKEEP/WENO切替え、連続混合、保存性
- 6次精度粘性項と拡散時間刻み
- 周期境界の面、辺、角ghost
- HITの等方化数学、目標Mach／Reynolds導出
- Forcing係数とCPU／CUDAの整合
- CPU参照実装とCUDA SSPRK3結果の一致
- SLF変換と乱流統計後処理

## 15. 現行実装上の注意事項と制約

### 15.1 無次元化に直接関係する事項

1. `mach_number`は方程式に現れるMach数係数ではなく、初期速度振幅である。
2. `reynolds_number`は基準音速に基づくReynolds数である。
3. Taylor–Greenの平均圧力は `1/gamma` 固定で、`rho0/gamma`ではない。
   したがって `rho0=1`を使用すること。
4. HITの目標Reynolds数からの導出も、粘性実装との整合上 `rho0=1`を前提とする。
5. ソルバーは有次元基準量を保存しない。再有次元化に使う基準量はケース管理側で別途記録すること。

### 15.2 YAML上に存在するが現行計算へ反映されない項目

- `flow.taylor_green.amplitude`
- `flow.taylor_green.mode_x`
- `flow.taylor_green.mode_y`
- `flow.taylor_green.mode_z`
- `physics.nse.nghost`（`grid.nghost`が実際の値）
- `numerics.time_integration`（`numerics.time_integrator`が正規項目）

### 15.3 入出力上の制約

- `output.format`を変えてもソルバー本体はSLFを書き出す。
- 計算精度はfloat64固定で、`output.precision`は実質的にメタデータである。
- CPU版は現状 `output.write_meta`にかかわらず `meta.json`を書き出す。
- 最終stepは無条件保存されない。
- 時刻とstepを継続する再スタート読込みは未実装。保存済み乱流をstep 0の初期条件として読む機能とは区別する。
- Linuxで出力先が未作成の場合は、実行前にディレクトリを作成するのが安全である。

### 15.4 数値・物理モデル上の制約

- 一様直交格子、三方向周期境界のみ。
- 5保存変数の完全気体のみ。
- 粘性係数と熱伝導係数は一定。
- 化学反応、多成分、LES/RANS、重力、一般物体力は未実装。
- positivity-preserving limiterは未実装。
- `small_rho`、`small_p`は評価保護であり、保存状態の修復ではない。
- CUDA版は単一GPUとMPI＋CUDAマルチGPUに対応する。MPI＋CUDAのhalo通信は
  host staging方式であり、分散HIT初期化とFFT forcingにはLinux上のcuFFTMpが必要である。
- CPU版は最低4 MPIプロセスを要求し、CPU逐次profileはない。

## 16. 実装ファイル対応表

| 仕様 | 主実装 |
|---|---|
| 設定・HIT輸送係数 | `src/common/mod_model_config.f90` |
| 入力namelist | `src/io/mod_input_reader.f90` |
| 格子 | `src/grid/mod_grid_fvm.f90` |
| 保存変数配列 | `src/field/mod_field_nse.f90` |
| Taylor–Green | `src/init/mod_init_taylor_green.f90` |
| HIT CPU | `src/init/mod_init_hit_spectral_2decomp.f90` |
| HIT GPU | `src/init/mod_init_hit_spectral_cufft.f90` |
| 周期境界 | `src/boundary/mod_boundary_periodic.f90` |
| KEEP | `src/numerics/convective/mod_convective_keep.f90` |
| WENO5-Z | `src/numerics/reconstruction/mod_reconstruction_weno5z.f90` |
| Roe | `src/numerics/riemann/mod_riemann_roe.f90` |
| ハイブリッド | `src/numerics/convective/mod_convective_hybrid.f90` |
| 粘性・熱伝導 | `src/numerics/viscous/mod_viscous_central6.f90` |
| 空間演算 | `src/numerics/mod_nse_spatial_operator.f90` |
| SSPRK3・CFL | `src/time/mod_nse_time_integration.f90` |
| Forcing | `src/forcing/mod_nse_forcing_2decomp.f90` |
| CPU実行制御 | `src/main/main_nse.f90` |
| CUDA実行制御 | `src/main/main_nse_cuda.f90` |
| SLF出力 | `src/io/mod_slf_output.f90` |
| case.yaml変換 | `ScriptLibrary/RunEnvironment/case_input.py` |
| profile選択 | `ScriptLibrary/RunEnvironment/profile_selection.py` |

## 17. 仕様変更時の管理方針

次のいずれかを変更した場合は本書も同時に更新する。

- 支配方程式、状態方程式、無次元化
- 初期条件式またはHIT目標値の定義
- 対流、粘性、時間積分、境界条件
- Forcingのエネルギー取扱い
- profile、MPI分割、CUDA対応範囲
- case.yamlの正規キー、既定値、入力制約
- SLF仕様、出力タイミング、後処理定義

特に、Taylor–Greenの `rho0`対応、未使用の振幅・mode指定、出力selectorを
将来修正した場合は、15章の注意事項を削除または更新し、回帰テストを追加すること。
