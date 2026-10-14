# SLFの空間FFT解析

利用者向け入口は実行環境直下の`tools/postprocess_case.py --task fft`。
内部の`SolverLibrary/NSE/tools/slf_fft.py`が単成分NSEのSLFを読み、指定した1変数の3次元FFTと
球殻積算スペクトルを計算する。時間方向のFFTではない。
元SLFは変更しない。NumPyが必要で、解析自体はCPU・単一プロセスで行う。

## Windows（PowerShell）

最新版から生成した実行環境のルートで実行する。以下のcase番号は実際の環境に変更する。
既存環境には自動反映されない。更新前に計算結果・設定を保全し、
計算結果がある環境への`prepare_environment.py --overwrite`は避ける。

```powershell
python -m pip install numpy
cd C:/Users/Owner/ResearchRuns/nse_case0022
python tools/postprocess_case.py --task fft --fft-field rho --steps latest --save-fft
```

## Linux（bash）

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy
# 生成した実行環境のルートから実行する。
cd /path/to/nse_case0022
python tools/postprocess_case.py --task fft --fft-field rho --steps latest --save-fft
```

新規生成される単成分NSE実行環境には本ツールとSLF読み込み依存ファイルが含まれる。
その実行環境のルートでは、Windows/Linux共通で次のように実行できる。

```text
python tools/postprocess_case.py --task fft --fft-field p --steps 0:1000:100
python tools/postprocess_case.py --task fft --fft-field u --steps all
python tools/postprocess_case.py --task all
```

入力は既定で`cases/<case_id>/output`、FFT出力は`cases/<case_id>/fft`。
`--input-dir`、`--meta`、`--fft-output`で変更できる。比熱比はcase.yamlから参照し、
必要時だけ`--gamma`で上書きする。`--dry-run`はコマンド表示のみ。
`--task all`はParaView・乱流統計・FFTの3処理。既定の対象はParaView/FFTがlatest、統計がall。
`--steps`指定時は3処理に共通で適用する。FFTは既定rho、平均除去、窓なし。
従来の`--output-dir`、`--fields`、`--stride`は可視化用で、FFTには適用しない。
圧力の導出は単成分・一定比熱比の理想気体に限る。
保存変数（`rho`, `rho_u`, `rho_v`, `rho_w`, `rho_E`など）は直接選択可能。
`u`, `v`, `w`, `p`は保存変数から導出する。1回につき1変数を指定する。
速度各成分のスペクトルは運動エネルギースペクトルではない。
体積平均の速度変動エネルギーを得るなら、同条件のu,v,wのshell_powerを足して1/2倍する。
これは圧縮性流れの密度重み付きエネルギーとは異なる。

## MPIと入力の注意

- `field_000100_rank00000.slf`等の通常の出力名を使用する。
- MPI出力は1rankのファイルではなく、全rankがある出力ディレクトリを渡す。
- 同じ計算の`meta.json`を自動読み込みする。別配置なら`--meta /path/to/meta.json`を指定。
- ペンシル／スラブいずれもrank_rangesで結合する。解析側のMPIプロセス数設定は不要。
- ghostを除去する。欠損、重複、重なり、時刻不一致、非有限値はエラーとする。
- `latest`は最大ステップを選ぶが、不完全な最新出力を以前のステップに自動差し替えない。
- globalとrankファイルが混在するときは`--layout global`または`--layout rank`で指定できる。
- 一様直交格子専用。一般座標・ノズル格子のデータには使用しない。

## 正規化と出力

既定では体積平均を引き、`F = fftn(f - mean(f))/N`を計算する。
`--fft-keep-mean`で平均を保持する。`--save-fft`がなければ大きなFFT配列は保存しない。

- `fft_000100_rho.csv`: `k, mode_count, shell_power, power_per_unit_k`
- `fft_000100_rho.json`: 入力、時刻、格子、平均値、正規化、Parseval確認値
- `fft_000100_rho.npz`: オプション。複素`fft`と角波数軸`kx,ky,kz`

`k`はSLFの長さ単位に対する角波数（2π/長さ）。NPZはfftshiftしていないFFT標準順序。
正負両方の波数を含み、片側スペクトルの2倍補正はしない。
球殻幅は`min(2π/Lx,2π/Ly,2π/Lz)`、最近接の殻中心へ積算する。
殻の端は中心±幅/2（最初の殻は0から）。
`sum(shell_power) = mean(processed_field**2)`であり、窓なし・平均除去時は分散になる。
`power_per_unit_k = shell_power / shell_width`。モード数平均ではない。
高波数殻は格子のNyquist境界で欠けるので、球殻全体が存在する範囲と区別する。
球殻積算は等方性を証明するものではない。

FFTは周期的延長を仮定する。衝撃波管や局所乱流では端の不連続による漏れがある。
必要なら`--fft-window hann`を指定する。各軸のHann窓を乗じ、窓の二乗平均による補正は行わない。
この場合Parsevalが比較するのは窓を掛けた場であり、元の場の分散ではない。

全領域をホストメモリに展開する。FFT・球殻作業配列だけでも目安100 byte/セル以上に加え、
SLF読み込み配列が必要。2048³級は1 TB規模になり得るため本ツールの通常対象外。
分散FFT／GPU解析は未対応。まず小さなデータで確認する。
再出力する場合だけ`--fft-overwrite`を指定する。
