# 分散FFTによる一様等方性乱流初期条件

## 概要

`initial_condition = 'hit_spectral'`を指定すると、2DECOMP&FFTとFFTW3を使って
一様等方性乱流（HIT）の初期速度場を生成する。

初期化処理でも計算本体と同じx-pencil分割を使用する。rootプロセスへの三次元場の
集約は行わない。

処理手順は次の通り。

1. z-pencil上の各波数について、Johnsen型またはPope型の`E(k)`を評価する。
2. seedと全体波数番号から、2組の再現可能なランダム位相を生成する。
3. 波数に直交する2本の基底へランダム複素係数を展開し、発散のない
   `u_hat`、`v_hat`、`w_hat`を直接構築する。
4. `u_hat(-k) = conjg(u_hat(k))`を明示的に課して、実数速度場を保証する。
5. 2/3則を既定値として高波数成分を除く。
6. 2DECOMP&FFTとFFTW3による分散逆FFTでx-pencilへ戻す。
7. x-pencilを計算本体の領域分割へMPI再分配する。
8. 全ランクの速度二乗和を集約し、1成分RMSを指定値へ正規化する。
9. 一様密度・一様圧力と組み合わせて保存変数を構築する。

実空間乱数場を先に生成して順FFTする処理は行わない。元コードと同じく、
エネルギースペクトルと乱数から波数空間速度を直接生成し、逆FFTだけを
ライブラリへ任せる。

## ビルド

PowerShellでは、インストール済み2DECOMP&FFTの場所を指定して構成する。

```powershell
$env:DECOMP2D_ROOT = "C:\Users\Owner\Documents\Codex\ThirdParty\install\2decomp-fft-v2.1.0"

cmake --fresh `
  -S . `
  -B build\hit-2decomp `
  -G Ninja `
  -DCMAKE_Fortran_COMPILER=C:/msys64/ucrt64/bin/gfortran.exe `
  -DNSE_USE_MPI=ON `
  -DNSE_MPI_PROVIDER=MSMPI `
  -DNSE_INIT_FFT_BACKEND=2decomp_fftw `
  -DNSE_2DECOMP_ROOT="$env:DECOMP2D_ROOT"

cmake --build build\hit-2decomp --parallel 8
```

YAMLビルドでは`config/build.hit.yaml`を使用する。

```powershell
python .\tools\build_from_yaml.py .\config\build.hit.yaml --build
```

統合フレームワークではsolver profileに`cpu_mpi_2decomp_fftw`を指定する。

流れ場ごとの`environment.*.yaml`は作成しない。
NSE共通の`ScriptLibrary/RunEnvironment/environment.nse.yaml`で
次の選択肢を指定する。

```yaml
select:
  model: nse_cpu_mpi_2decomp_fftw
  case: nse_case
  execution: release

parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: false
```

ここで`model`は分散FFTを利用できるビルド構成を選び、
`case`はNSE共通の`case.yaml`入力雛形を選ぶだけである。
HITのスペクトル、乱数seed、RMS速度などの流れ場条件は、
生成後の`cases/caseNNNN/case.yaml`で管理する。
ビルド時と実行時にケーステンプレートを直接渡す必要はない。

`ScriptLibrary/BuildSolver/build.local.yaml`では次のように指定できる。

```yaml
models:
  nse:
    profile: cpu_mpi_2decomp_fftw
    cmake_overrides:
      NSE_2DECOMP_ROOT: C:/Users/Owner/Documents/Codex/ThirdParty/install/2decomp-fft-v2.1.0
```

または、ビルド前に`DECOMP2D_ROOT`環境変数を設定する。

## 入力条件

`&simulation`で初期条件を選ぶ。

```fortran
initial_condition = 'hit_spectral'
```

`&nse`でスペクトル条件を指定する。

```fortran
hit_spectrum = 'johnsen'
hit_seed = 13579
hit_rms_velocity = 0.1
hit_peak_wavenumber = 4.0
hit_integral_length = 1.0
hit_kolmogorov_length = 0.02
hit_dealias_fraction = 0.6666666666666667
```

ケースYAMLから生成する場合は次の形で指定する。

```yaml
flow:
  type: hit_spectral
  hit:
    spectrum: johnsen
    random_seed: 13579
    rms_velocity: 0.1
    peak_wavenumber: 4.0
    integral_length: 1.0
    kolmogorov_length: 0.02
    dealias_fraction: 0.6666666666666667
```

NSE共通雛形は`ScriptLibrary/RunEnvironment/case_templates/nse.yaml`にある。

| パラメータ | 意味 |
|---|---|
| `hit_spectrum` | `johnsen`または`pope` |
| `hit_seed` | 波数番号ベースの再現可能な乱数seed。MPI分割数に依存しない |
| `hit_rms_velocity` | 無次元の1速度成分RMS。m/sではない。0以下なら`mach / sqrt(3)`を使用 |
| `hit_peak_wavenumber` | Johnsen型スペクトルのピーク波数 |
| `hit_integral_length` | Pope型スペクトルの積分長さ |
| `hit_kolmogorov_length` | Pope型スペクトルのKolmogorov長さ |
| `hit_dealias_fraction` | 各軸のNyquist波数に対する保持率 |

動作確認用入力は`tests/data/input_hit_2decomp_small.dat`にある。

```powershell
mpiexec -n 4 `
  .\build\hit-2decomp\bin\solver.exe `
  .\tests\data\input_hit_2decomp_small.dat
```

正常なら、標準出力に次の診断値が表示される。

- 2DECOMP process grid
- max spectral divergence
- max inverse FFT imaginary residual
- unscaled component RMS
- target component RMS
- initial turbulent Mach number

## 現在の物理モデル

- 速度場は完全にsolenoidalで、初期のdilatational成分はゼロ。
- 平均速度はゼロ。
- 密度は`rho0`で一様。
- 圧力は`rho0 / gamma`で一様。この無次元化では初期音速が1になる。
- 1成分RMSを`u_rms`とすると、初期乱流Mach数は`sqrt(3) * u_rms`になる。
- `hit_rms_velocity <= 0`の場合は、指定した`mach`が初期乱流Mach数になるよう
  `hit_rms_velocity = mach / sqrt(3)`を使用する。
- 例えば`hit_rms_velocity = 100`は100 m/sではなく、初期乱流Mach数約173の
  超音速条件になるため、通常の低Mach HITには使用しない。

初期条件生成だけが擬スペクトル処理であり、その後の時間発展は現在の
有限体積NSEソルバーを使用する。
