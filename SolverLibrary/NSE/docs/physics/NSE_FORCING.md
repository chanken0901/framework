# 統計的定常乱流のForcing

負圧対策と外力を含む更新前の時間刻み検査は
[`NSE_POSITIVITY.md`](../numerics/NSE_POSITIVITY.md)を参照。

> Windows（PowerShell）とLinux（bash）のコマンド対応は[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

## 実装範囲

`forcing.type`は`none`、`petersen_livescu`、`watanabe2021`から選ぶ。
以下のCPU/GPU対応表と説明は従来の`petersen_livescu`についてのもの。
`watanabe2021`の追加範囲と制限は末尾の節を参照。

Petersen and Livescu (2010) の線形Forcingを実装している。密度重み付き速度
`w_i = sqrt(rho) u_i`をFFTし、波数空間のヘルムホルツ分解でソレノイダル成分と
ダイラテーショナル成分へ分ける。

Forcingは各SSPRK3段の運動量RHSへ加算される。論文の内部エネルギー補償
`f_e = -f_i u_i`に対応して、保存全エネルギーのRHSへ直接のForcingは加えない。

| 実行方式 | ヘルムホルツ分解 | 状態 |
|---|---|---|
| CPU + MPI/OpenMP | 2DECOMP&FFT + FFTW3 | 実装・動作確認済み |
| 単一GPU | cuFFT | 実装・動作確認済み |
| MPI + CUDA | cuFFTMp | 実装済み（LinuxマルチGPU環境で実機検証が必要） |

CPU版はroot集約を行わず、既存のx-pencil分割と2DECOMPの転置を使う。単一GPU版は
速度場、スペクトル、射影、逆FFT、RHS加算をGPU上で処理し、係数を求めるスカラー
だけをホストへ転送する。

MPI＋CUDA版もroot集約を行わない。局所状態とFFT記述子をGPU上に保持し、
Helmholtz射影後の分母と圧力膨張相関だけを`MPI_Allreduce`する。詳細は
[`NSE_CUFFTMP.md`](../parallel/NSE_CUFFTMP.md)を参照する。

## case.yaml

`case.yaml`へ次を追加する。

```yaml
forcing:
  type: petersen_livescu
  petersen_livescu:
    spectrum: low_wavenumber
    fft_backend: auto
    k_cutoff: 2.5
    target_dissipation: 0.1
    dilatational_ratio: 0.0
    denominator_floor: 1.0e-14
    max_coefficient: 0.0
    report_interval: 100
```

| 指定子 | 値・意味 |
|---|---|
| `type` | `none`または`petersen_livescu` |
| `petersen_livescu.spectrum` | `full_spectrum`または`low_wavenumber` |
| `petersen_livescu.fft_backend` | 通常は`auto`。明示時は`2decomp_fftw`または`cufft` |
| `petersen_livescu.k_cutoff` | `low_wavenumber`で保持する物理波数。条件は`abs(k) < k_cutoff` |
| `petersen_livescu.target_dissipation` | ソレノイダルとダイラテーショナルを合わせた目標注入率 |
| `petersen_livescu.dilatational_ratio` | 目標値`epsilon_d / epsilon_s` |
| `petersen_livescu.denominator_floor` | 係数計算でゼロ割を防ぐ下限 |
| `petersen_livescu.max_coefficient` | 係数の絶対値上限。`0`は制限なし |
| `petersen_livescu.report_interval` | Forcing評価回数ごとの診断間隔。`0`は非表示 |

`scheme`は既存caseとの後方互換用エイリアスとしてのみ残している。新しい
case設計では`type`を使用する。

目標注入率は`dilatational_ratio=r`から次のように分ける。

```text
epsilon_s = target_dissipation / (1 + r)
epsilon_d = target_dissipation - epsilon_s
```

圧力膨張相関は`PD=<p div(u)>`として6次中心差分で評価し、係数は次式で計算する。

```text
c_s = epsilon_s / <w_s . w>
c_d = (epsilon_d - PD) / <w_d . w>
```

`full_spectrum`では分母がそれぞれ`2 K_s`、`2 K_d`に一致する。

## プロファイル

environment設計書ではMPI/CUDAの使用有無だけを選ぶ。生成環境には互換profileが
同梱され、`forcing.type: petersen_livescu`を検出した`run_case.py`が、CPUでは
`cpu_mpi_2decomp_fftw`、単一GPUでは`cuda_single`を自動選択する。
`solver.profile`を重複指定する必要はない。

Windowsでは`python .\tools\run_case.py --prepare`、Linuxでは
`python3 ./tools/run_case.py --prepare`で`input.dat`へ展開される。MPIプロセス数と
OpenMPスレッド数は従来どおり`case.yaml`または実行時オプションで指定する。

## cuFFTMp実装

`cuda_mpi_cufftmp`ではGPU常駐、root集約禁止、3成分分散FFT、局所スペクトル
射影、係数の全ランク縮約、局所RHS加算を実装している。現在のWindows開発機には
cuFFTMp/NVSHMEMがないため、ここでは契約・構文と既存バックエンドの回帰を検証し、
分散FFTの実機検証はLinuxマルチGPU計算機で行う。

## 診断出力

`# forcing`行には次を出力する。

```text
評価回数 c_s c_d 分母_s 分母_d PD epsilon_s epsilon_d
```

CPUとGPUの同一小規模ケースでは、この値が丸め誤差の範囲で一致することを確認する。

## Watanabe・Tanaka・Nagata (2021)：ソレノイダル線形forcing

追加日：2026-10-23。Phys. Fluids 33, 095108, DOI: 10.1063/5.0062596、式(8)〜(14)。
既存方式とは独立して`forcing.type: watanabe2021`で選ぶ。

```yaml
forcing:
  type: watanabe2021
  watanabe2021:
    turbulent_mach_number: 0.3
    taylor_reynolds_number: 100.0
    feedback_exponent: 4.0
    fft_backend: auto
    report_interval: 100
    # mean_density: 1.0
    # mean_pressure: 0.7142857142857143
```

必須物理パラメータは目標乱流マッハ数Mtとテイラー尺度基準Re_lambda。
背景密度は解決済みrho0、背景圧力はrho0/gammaを既定とする。異なる平均温度では
対応するmean_density/mean_pressureを明示する。計算中の温度から目標を再計算しない。
feedback_exponentは論文のp、既定4。0なら分散フィードバックなし。負値は拒否する。

既存の無次元化・一定粘性mu=1/Reinputを維持し、次を入力生成時に計算する。

```text
a0² = gamma*p0/rho0
u0 = Mt*a0/sqrt(3)                    # 1成分の目標RMS
nu0 = 1/(Reinput*rho0)
Re_L0 = Re_lambda0²/15
L0 = Re_L0*nu0/u0                     # 論文の散逸ベース積分スケール
epsilon0 = u0³/L0
forcing_target_dissipation = rho0*epsilon0
forcing_target_velocity_variance = 3*u0²
```

Re_L0とL0はinput.datの`forcing_target_integral_reynolds`、
`forcing_target_integral_length`にも記録し、起動時に表示する（来歴用、力の計算には
変換済み散逸・分散を使用）。Reinputを上書きしない。HIT初期スペクトルの積分尺度や
その換算係数とは定義が違うため流用しない。Re_lambda=100ならRe_L0=666.6666667。
この換算は等方性に基づく目標設定であり、実測Re_lambdaの一致を保証しない。

各RK段でw=sqrt(rho)*uをFFTし、ゼロモードを除いて全波数でHelmholtz射影する。
体積平均を<>、D=<p div(u)>、V=<|u-<u>|²>、S=<|w_s|²>として

```text
alpha = (V/(3*u0²))^(-p)
A = alpha*(rho0*epsilon0-D)/S
f_i = A*sqrt(rho)*w_s,i
Q = D - <tau:grad(u)>
RHS(rho*u_i) += f_i
RHS(rho*E) += u_i*f_i + Q
```

論文のfeはrho*T方程式の源項なので、保存全エネルギーにはfe/(gamma-1)=Qを使う。
Qは空間一様で、平均内部エネルギーの上昇を抑える補償として常に適用する。
熱源項だけを無効にする設定は設けない。Qの符号を負にクリップする処理でもなく、
圧力膨張と粘性散逸の差をそのまま使う。
圧力膨張相関・応力散逸には6次中心の速度勾配と一定粘性を用いる。
速度平均・分散・散逸は全rankで縮約し、局所領域の統計で代用しない。
追加ログ`# forcing variance dissipation thermal_source`はV、応力散逸、Qを示す。
従来方式の運動量のみの加算経路は変更していない。

### 対応範囲・検証状況

- CPU/MPI + 2DECOMP&FFT/FFTW、単一CUDA/cuFFT、MPI+CUDA/cuFFTMpに実装済み。
  全方向周期、等間隔格子、単成分、CENTRAL6粘性を必要とする。
- profileはそれぞれ`cpu_mpi_2decomp_fftw`、`cuda_single`、`cuda_mpi_cufftmp`。
  CPU/OpenMPでは既存master領域で評価する。GPUでは統計・熱源・外力仕事をGPU上で計算し、
  MPI+CUDAでは統計を全rankで縮約する。
- FV2粘性、揺らぎ併用は明示的に拒否する。
- 低波数帯限定、膨張成分の加振、係数クリップはこの方式では許可しない。
- 論文の温度依存粘性・8次差分・フィルタ・RK方式を導入したものではない。
- 平均温度を保つ式の実装であって、有限dt・既存離散演算子のもとでの厳密な温度固定ではない。
  平均温度のドリフト、時間刻み・解像度依存を必ず検証する。
- 入力変換・係数の単体試験に加え、単一CUDA/cuFFTで短時間の温度・全エネルギー収支を確認。
  CPU/MPIでの新方式の温度維持・複数rank一致、cuFFTMp実機での新方式の時間発展、
  長時間の定常HIT統計は未検証。現状を本番検証済みとは扱わない。
- 線形forcingの実現積分尺度は領域長にも依存する。Mt/Re_lambdaを設定するだけで
  任意の領域に目標乱流を実現できるわけではなく、論文の領域長依存の議論も確認する。

新しく環境を生成し、Windowsでは`python tools/run_case.py --prepare --build`、
Linuxでは`python3 tools/run_case.py --prepare --build`で準備する。
既存環境での`--prepare`だけでは新しいFortranコードはコピーされない。
検証は小格子・短時間から始め、Mt、Re_lambda、平均温度、エネルギー収支を比較する。

### 温度・全エネルギー回帰試験（2026-10-27）

`tests/check_watanabe_temperature.py`は単一rankの保存変数SLFから
`T_v=<p/rho>`（体積平均）と`T_m=<p>/<rho>`（質量加重平均）を別々に評価する。
ここで温度は`p/rho`に比例する量であり、相対ドリフトを比較する。
周期領域の連続式でQが補償するのは平均内部エネルギー、したがって質量保存下の`T_m`である。
`T_v`や各点の温度を厳密に固定する条件ではない。

Taylor–Green、領域`[0,2*pi]^3`、gamma=1.4、mach=0.4、Reinput=10、Pr=0.72、
KEEP6/CENTRAL6/SSPRK3、目標散逸0.012、目標速度分散0.04、指数4、終了時刻0.1を使用。
初期・終端を含む11時刻すべてを検査し、負密度・負圧・非有限値・時刻の重複や欠落を拒否する。
各段のログから`Q=D-dissipation`と係数式を検査し、`A*S+Q`をSSPRK3の重み
`(1/6,1/6,2/3)`で積分して、各保存時刻の全エネルギー変化と照合する。
診断数が`3*nsteps`と異なる再試行・substepは、この固定dt試験では失敗として扱う。

| 格子 | dt | T_m最大相対ドリフト | T_v最大相対ドリフト |
|---|---|---|---|
| 16³ | 0.001 | 3.33633e-8 | 2.76354e-7 |
| 16³ | 0.0005 | 3.33633e-8 | 2.76354e-7 |
| 32³ | 0.001 | 5.33727e-10 | 3.00871e-7 |

現ソースを再ビルドした単一CUDA/cuFFTで上記を再確認した。平均密度の最大相対変化は
1.18e-14以下、全エネルギー収支の最大絶対残差は2.35e-13以下。
時間刻み半減ではT_mの誤差はほぼ変わらず、格子細分化で減少した。
これは短時間・指定条件での結果であり、時間精度次数や長時間定常状態の証明ではない。

合格基準はT_mが1e-6以下、T_vが1e-5以下、平均密度が1e-12以下。
細分化後のT_mドリフトは粗い条件+1e-10以下とする。
エネルギー収支の絶対許容誤差は`1e-10+5e-5*integral(|A*S|+|Q|)dt`。
係数・熱源照合の相対許容誤差は5e-5・3e-5（絶対1e-10）で、CPUログの表示桁数も考慮する。
入力、run.log、SLF、temperature_report.jsonは新規の検証出力ディレクトリに保持する。
既存ディレクトリは上書きしない。実行ファイルのDLL・ランタイムをPATHに設定する。

```powershell
python SolverLibrary/NSE/tests/check_watanabe_temperature.py --exe build/nse-positivity-cuda/bin/nse_cuda.exe --work-dir build/watanabe-temperature-new
python -m unittest discover -s SolverLibrary/NSE/tests -p test_watanabe_temperature.py
```

Linuxでは`python`を`python3`に、`--exe`をその環境の単一rank実行ファイルへ変更する。
NumPyが必要。検査器の単体試験は欠落出力、異常状態、誤った熱源・係数、外力仕事または
熱源を落とした全エネルギー履歴が合格しないことを確認する。
