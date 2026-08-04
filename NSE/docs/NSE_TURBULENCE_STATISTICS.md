# NSE乱流統計後処理

## 目的

`tools/nse_turbulence_statistics.py`は、NSEが保存したrank別または全領域SLFを
時刻ごとに統合し、乱流統計をCSV形式で出力する。ソルバーの時間発展には処理を
追加しないため、CPU、MPI/OpenMP、単一GPUのいずれの計算結果にも同じ後処理を
適用できる。

## 実行環境からの使用

全ての完全な保存時刻を処理する。

```powershell
python .\tools\postprocess_case.py --task statistics
```

保存ステップを限定する。

```powershell
python .\tools\postprocess_case.py --task statistics --steps 0:1000:100
```

ParaView変換も続けて行う。

```powershell
python .\tools\postprocess_case.py --task all
```

出力先は既定で次の2ファイルになる。

- `cases/<case_id>/statistics/turbulence_statistics.csv`
- `cases/<case_id>/statistics/turbulence_statistics_metadata.json`

`gamma`と基準Reynolds数は`case.yaml`の`physics.nse.gamma`および
`physics.nse.reynolds_number`から取得する。コマンドラインの`--gamma`、
`--reynolds`で一時的に上書きできる。

## 統計量

速度変動は各成分の体積平均を差し引いたReynolds変動として定義する。

```text
u'_rms = sqrt(<(u-<u>)^2>)
q_rms  = sqrt(<u_i' u_i'>)
K      = 0.5 <u_i' u_i'>
```

積分スケールは三次元FFTから得た離散エネルギースペクトルを用いる。

```text
L = (3*pi/4) * sum(E(k)/k) / sum(E(k)),  k > 0
```

散逸率はスペクトル微分を用いて次式から計算する。

```text
epsilon = nu <2 S_ij S_ij - (2/3) (div u)^2>
nu      = 1 / (Re * <rho>)
```

等方性を仮定する代表スケールとReynolds数は次式とする。

```text
u'        = sqrt(<u_i' u_i'>/3)
lambda    = sqrt(15 nu u'^2 / epsilon)
eta       = (nu^3 / epsilon)^(1/4)
Re_L      = u' L / nu
Re_lambda = u' lambda / nu
M_t       = q_rms / <sqrt(gamma p/rho)>
```

CSVには平均密度・圧力・音速、3方向の平均速度と変動RMS、Reynolds応力の
非対角成分、等方性誤差、散逸率、Parseval整合誤差も出力する。

## 注意事項

- 周期境界を前提にFFT微分と積分スケールを計算する。
- 統計計算は空間間引きを行わない。
- 非正の密度・圧力、NaN、rank欠落、領域重複を検出した場合は停止する。
- `lambda`、`eta`、`Re_lambda`は局所等方性を仮定した代表値であり、強い
  圧縮性・異方性がある場合は方向別Reynolds応力と併せて解釈する。
- 1時刻ずつ処理するが、三次元FFTのため全領域の速度3成分と複素スペクトルを
  保持できるメモリが必要になる。
