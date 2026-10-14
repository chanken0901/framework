# 衝撃波乱流干渉専用の後処理

## 入口と対応範囲

通常の可視化・統計・FFT用`postprocess_case.py`とは別に、
生成済み実行環境の`tools/analyze_shock_turbulence.py`から実行する。
解析本体は`SolverLibrary/NSE/tools/nse_shock_analysis.py`。
既存の後処理と同じPython/NumPyによるオフライン解析で、ソルバーやSLFを書き換えない。
NSEのmanifestと環境生成処理へ登録済み。新しく生成する環境には同梱される。
既存のResearchRunsへは自動で追加しない。

対象は単成分NSE、固定・一様な直交格子、x方向に伝播する衝撃波。
global SLFとrank分割SLFに対応し、rank_rangesとゴーストセルを考慮する。
rankデータの欠落・重複、時刻・格子不一致、非正の密度・圧力は拒否し、クリップしない。
rankレイアウトにはmeta.jsonの明示的なrank_rangesが必要。
SLFを2回読み、横断面の平均と中心化した二次モーメントを集計する。
解析自体のMPI/CUDA並列化は未実装。線ごとの波面検出用にglobal圧力配列を保持する。
追加メモリは少なくとも8*Nx*Ny*Nzバイト（2048³では64 GiB）と2D作業配列。
さらにSLFの読込み・変換用メモリが必要。global形式ではglobal SLF全体を読むため特に注意。

## 実行例

実行環境ルートで実行する。以下の範囲は**説明用**で、実際の領域長と波の位置に置き換える。
例は+x方向への衝撃波、探索範囲x=10〜30、衝撃波から2〜8の距離にある領域を統計対象とする。

Windows PowerShell：

```powershell
python .\tools\analyze_shock_turbulence.py --direction 1 --search 10 30 --upstream 2 8 --downstream 2 8 --steps all --dry-run
python .\tools\analyze_shock_turbulence.py --direction 1 --search 10 30 --upstream 2 8 --downstream 2 8 --steps all
```

Linux bash：

```bash
python3 tools/analyze_shock_turbulence.py --direction 1 --search 10 30 --upstream 2 8 --downstream 2 8 --steps all --dry-run
python3 tools/analyze_shock_turbulence.py --direction 1 --search 10 30 --upstream 2 8 --downstream 2 8 --steps all
```

`--dry-run`は呼出しコマンドを表示するだけで、SLF解析や出力作成はしない。
入力はenvironment.lock.jsonで指定されたケースのoutput、gammaはcase.yamlの値を使用する。
`--case-directory`、`--input-dir`、`--meta`、`--gamma`で明示指定も可能。

| 引数 | 意味 |
|---|---|
| --direction 1 / -1 | +x / -xへの伝播。自動推定しない |
| --search XMIN XMAX | 衝撃波を探す固定のx範囲。対象の全時刻の位置を含める |
| --upstream NEAR FAR | 衝撃波から進行方向側の未通過領域までの距離 |
| --downstream NEAR FAR | 衝撃波から逆側の通過後領域までの距離 |
| --max-shift D | 選択した隣接出力間の最大位置変化。指定時は追跡候補を制限 |
| --min-pressure-ratio R | 後方／前方の領域平均圧力比の下限。既定1.05 |
| --steps | all / latest / 10,20 / 0:1000:100。速度解析には複数時刻が必要 |
| --layout | auto / global / rank |
| --output-dir | 既定はcases/<case>/shock_analysis。既存ディレクトリへの上書きは禁止 |

座標・距離・時刻はSLFの値と同じ単位・無次元化で扱う。自動SI換算はしない。
NEAR>0、FAR>NEARが必要。衝撃波の数値的厚さ・波面の凹凸が統計窓に入らないようNEARを選ぶ。
窓の外端が計算領域外になる、または含まれる横断面が2未満の場合は停止する。
衝撃波が対象領域に入る前の時刻は`--steps`で除く。

再実行時は別の出力名を使う：

```bash
python tools/analyze_shock_turbulence.py --direction 1 --search 10 30 --upstream 2 8 --downstream 2 8 --output-dir cases/case0025/shock_analysis_v2
```

## 定義と出力

- `planes_<step>.csv`：x位置、横断面平均rho/p/u/v/w、Favre平均速度、Reynolds/Favre共分散6成分、乱流エネルギー。
- `shock_history.csv`：衝撃波位置・符号付き速度、前後領域統計、圧力比、乱流エネルギーと対角応力の増幅率。
- `analysis.json`：完了状態、入力ファイル一覧、解析パラメーター、各量の定義。
- `shock_surface_<step>.npz`：各(y,z)で検出した波面位置、検出マスク、圧力勾配・前後圧力比、座標と時刻。

完了判定は終了コード0とanalysis.jsonの`status=complete`。
途中停止では一部のCSVが残ることがあるが、完了メタデータは生成しない。
不完全なrank出力時刻は既存SLF選択処理に従って除外される。実際に使ったファイルはanalysis.jsonで確認する。

横断面平均を< >とすると、通常平均の変動はu'=u-<u>、
Favre平均はu_tilde=<rho*u>/<rho>、変動はu''=u-u_tilde。
R_ij=<u_i' u_j'>、F_ij=<rho*u_i''*u_j''>/<rho>。
k=0.5*(R_uu+R_vv+R_ww)、Favre kも同様。いずれも単位質量当たりの量。
前後領域の統計は**横断面ごとに平均流を除去した統計の平均**。
Favre量は各面平均密度で重み付けし、それ以外は面の算術平均。
x方向に平均速度が変化しても、その変化そのものを乱流エネルギーに加えない。

衝撃波位置は横断面平均圧力の差分から、進行方向に圧力が下がる最大勾配のセル面を採用する。
これは従来の代表位置shock_xの定義で、格子面単位の追跡。線ごとの波面検出は下記を参照。
速度は保存時刻を使う位置の数値微分。内部は非等間隔対応の中心差分、端は片側差分。
1時刻だけなら速度はNaN。位置の格子量子化と保存間隔により速度が揺れるため、
格子・出力時間間隔を変えた確認が必要。時間平滑化は自動適用しない。

増幅率は同時刻の後方／前方統計の比。分母<=1e-30の場合はNaNとし、0や無限大で代用しない。
同一流体塊の衝撃波通過前後を追跡した比ではない。
局所乱流では前方窓が乱流領域から外れると比が意味を失う。必ず横断面分布と併せて確認する。
圧縮性波動・熱揺らぎを乱流から分離するフィルターは未実装。

## 探索の注意と検証範囲

衝撃波管には接触面・膨張波・境界反射も存在する。圧力勾配だけで任意の流れの波を識別できるわけではない。
探索範囲とmax-shiftで対象波を隔離し、検出した位置を可視化で確認する。
折り返して単一値x_s(y,z)で表せない波面、複数衝撃波の同時追跡、y/z方向伝播、一般座標、反応流は対象外。
最低圧力比を満たさないときは、別の波を黙って採用せずエラーにする。

人工データで±x位置追跡、保存時刻からの速度、乱流増幅率、Favre重み付け、
rank分割＋ゴーストセルの一致、欠落・重複・不正状態・窓範囲外の拒否、
専用窓口からの実行と上書き防止を検証する。実計算データでの波の同定は別途必要。

## 各(y,z)の衝撃波位置 x_s(y,z,t)

従来と同じコマンドで、線ごとの波面検出も実行する。各(y_j,z_k)について、
横断面平均を取る前のp(x,y_j,z_k,t)から次の量を計算する。

```text
G[i+1/2,j,k] = -direction * (p[i+1,j,k]-p[i,j,k])/dx
x_s[j,k] = x[i_max+1/2], where i_max maximizes positive G within the search range
```

各線の検出位置を基準にupstream/downstream窓の圧力平均を取り、
圧力比がmin-pressure-ratio以上であることを検査する。
前後窓が領域外、2セル未満、正の勾配が見つからない、圧力比不足の場合、その線は無効。
平均圧力から得た位置で埋めたり、近隣線から補間したりはしない。
max-shiftは前の選択時刻の同じ(y,z)の有効位置に対して適用する。
前の位置が無効だった線は全探索範囲から再検出する。

NPZ配列は`shock_x[j,k]`、`valid[j,k]`、`pressure_gradient[j,k]`、`pressure_ratio[j,k]`。
`y[j]`、`z[k]`はセル中心座標、`time`と`step`は保存時刻とステップ。
無効線のshock_xおよびpressure_gradientはNaN。pressure_ratioには棄却判断用の値が残る場合がある。

shock_history.csvへ以下を追加する。

- surface_x_mean：有効な線の位置の算術平均。
- surface_x_std：有効な線の位置の母標準偏差（ddof=0）。
- surface_x_min / surface_x_max：有効な線の最小／最大位置。
- surface_valid_count / surface_total_count / surface_valid_fraction：有効検出数と割合。

有効な線の集合Vに対し、mean=sum_V(x_s)/|V|、std=sqrt(sum_V((x_s-mean)^2)/|V|)。
全て無効なら位置統計はNaN、有効率は0。完走しても有効率が低ければ波面統計をそのまま信用しない。
画面にも各時刻の有効数を表示する。

読み込み例（Windowsはpython、Linuxはpython3）：

```python
import numpy as np
with np.load('cases/case0025/shock_analysis/shock_surface_00001000.npz') as f:
    y, z, t = f['y'], f['z'], float(f['time'])
    xs, valid = f['shock_x'], f['valid']
    # xs[j,k] is x_s(y[j], z[k], t); invalid positions are NaN.
```

`shock_x`（平均圧力の最大勾配位置）と`surface_x_mean`（線ごとの検出位置の平均）は一般に一致しない。
既存のshock_speedと前後の乱流統計は引き続き従来のshock_xを基準とする。
波面に沿った前後領域統計、各線の速度、波面面積や傾斜は今回は追加していない。
位置は依然として格子面分解能。曲がった波面にも使えるが、x方向に単一値で検出できる必要がある。
各線の圧縮性揺らぎを別の衝撃波と誤認する可能性があるため、探索範囲・追跡距離・有効率を確認する。
