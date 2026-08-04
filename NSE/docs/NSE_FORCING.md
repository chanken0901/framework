# 統計的定常乱流のForcing

## 実装範囲

Petersen and Livescu (2010) の線形Forcingを実装している。密度重み付き速度
`w_i = sqrt(rho) u_i`をFFTし、波数空間のヘルムホルツ分解でソレノイダル成分と
ダイラテーショナル成分へ分ける。

Forcingは各SSPRK3段の運動量RHSへ加算される。論文の内部エネルギー補償
`f_e = -f_i u_i`に対応して、保存全エネルギーのRHSへ直接のForcingは加えない。

| 実行方式 | ヘルムホルツ分解 | 状態 |
|---|---|---|
| CPU + MPI/OpenMP | 2DECOMP&FFT + FFTW3 | 実装・動作確認済み |
| 単一GPU | cuFFT | 実装・動作確認済み |
| MPI + CUDA | cuFFTMp | 契約テストのみ。実計算は無効 |

CPU版はroot集約を行わず、既存のx-pencil分割と2DECOMPの転置を使う。単一GPU版は
速度場、スペクトル、射影、逆FFT、RHS加算をGPU上で処理し、係数を求めるスカラー
だけをホストへ転送する。

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

CPUでは`cpu_mpi_2decomp_fftw`、単一GPUでは`cuda_single`を選ぶ。

```yaml
solver:
  profile: cpu_mpi_2decomp_fftw
```

または

```yaml
solver:
  profile: cuda_single
```

`python .\tools\run_case.py --prepare`で`input.dat`へ展開される。MPIプロセス数と
OpenMPスレッド数は従来どおり`case.yaml`または実行時オプションで指定する。

## cuFFTMp拡張境界

現在の環境にはcuFFTMpとCUDA-aware MPIがないため、MPI + CUDA実行ファイルは
生成しない。将来のバックエンドが満たす契約は
`src/forcing/cufftmp_forcing_contract.yaml`に固定し、CTestで検証する。
この契約はGPU常駐、root集約禁止、3成分バッチ分散FFT、局所スペクトル射影、
係数の全ランク縮約、局所RHS加算を要求する。

## 診断出力

`# forcing`行には次を出力する。

```text
評価回数 c_s c_d 分母_s 分母_d PD epsilon_s epsilon_d
```

CPUとGPUの同一小規模ケースでは、この値が丸め誤差の範囲で一致することを確認する。
