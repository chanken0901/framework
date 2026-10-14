# KEEP/WENOハイブリッド対流流束 設計書

## 1. 目的

滑らかな領域では低散逸なKEEPを使い、衝撃波を検出した面では
WENO5-Z/Roeへ連続的に切り替える。各セルが別々に離散化を選ぶのではなく、
共有面に一つの混合流束を定義することで有限体積法の保存性を維持する。

CPU、MPI/OpenMP、単一GPU CUDAで同じ離散式と設定を使用する。

## 2. 入力設定

運用上の正本はcaseの`case.yaml`であり、対流方式は
`numerics.convective_scheme`へ指定する。`run_case.py --prepare`が生成する
`input.dat`では、対応する値が`&nse`ナムリストへ展開される。
既定値は`keep6`のままであり、ハイブリッドは明示的に選択した場合だけ有効になる。

旧`numerics.flux`と未使用だった`numerics.reconstruction`は受け付けない。
ハイブリッド以外を選んだ場合、`hybrid`以下の設定は計算にも生成後の
`input.dat`にも使用されない。

### case.yamlでハイブリッドを使う場合

```yaml
numerics:
  convective_scheme: hybrid
  hybrid:
    smooth_scheme: keep6
    shock_scheme: weno5z_roe
    sensor: ducros_pressure
    sensor_onset: 0.01
    sensor_full: 0.10
```

### 生成後のinput.dat

```fortran
&nse
  convective_scheme = 'hybrid'
  hybrid_smooth_scheme = 'keep6'
  hybrid_shock_scheme = 'weno5z_roe'
  hybrid_sensor = 'ducros_pressure'
  hybrid_sensor_onset = 0.01
  hybrid_sensor_full = 0.10
/
```

### ハイブリッドを使わない場合

```fortran
&nse
  convective_scheme = 'keep6'
/
```

`convective_scheme`には次の値を指定できる。

| 値 | 計算方法 |
|---|---|
| `keep2` | 二次精度KEEP |
| `keep6` | 六次精度KEEP（既定値） |
| `weno5z_roe` | 五次精度WENO-Z再構築とRoe流束 |
| `hybrid` | センサーに基づく二つの面流束の連続混合 |

ハイブリッド専用項目は`convective_scheme`が`hybrid`のときだけ使用される。

| 項目 | 既定値 | 意味 |
|---|---:|---|
| `hybrid_smooth_scheme` | `keep6` | 滑らかな領域で使う流束 |
| `hybrid_shock_scheme` | `weno5z_roe` | 衝撃波領域で使う流束 |
| `hybrid_sensor` | `ducros_pressure` | 混合率を求めるセンサー |
| `hybrid_sensor_onset` | `0.01` | WENO混合を開始するセンサー値 |
| `hybrid_sensor_full` | `0.10` | WENOへ完全移行するセンサー値 |

現在、二つの構成流束には`keep2`、`keep6`、`weno5z_roe`を指定できる。
推奨構成は既定値の`keep6`と`weno5z_roe`である。センサー閾値は無次元で、
必ず`0 <= hybrid_sensor_onset < hybrid_sensor_full`を満たす必要がある。

## 3. 混合則

面`i+1/2`で滑らか側流束を`F_s`、衝撃波側流束を`F_w`、
混合率を`alpha`とすると、保存流束は次式である。

```text
F_h(i+1/2) = (1 - alpha(i+1/2)) F_s(i+1/2)
             + alpha(i+1/2) F_w(i+1/2)
```

`alpha=0`では`hybrid_smooth_scheme`、`alpha=1`では
`hybrid_shock_scheme`となる。既定構成ではそれぞれKEEP6とWENO5-Z/Roeである。
一つの共有面に一つの`F_h`を保存し、隣接セルは同じ流束を逆符号で使うため、
混合率が空間的に変化しても保存性は失われない。

## 4. Ducros-pressureセンサー

セルの圧力曲率と圧縮性Ducros係数を組み合わせる。

```text
kappa_p = abs(p(+1) - 2 p(0) + p(-1))
          / (p(+1) + 2 p(0) + p(-1) + epsilon)

D = min(div(u), 0)^2
    / (div(u)^2 + abs(curl(u))^2 + epsilon)

S_cell = clamp(kappa_p D, 0, 1)
S_face = max(S_left_cell, S_right_cell)
```

`min(div(u),0)`により圧縮領域だけを検出し、渦度が支配的な滑らかな乱流領域で
不要なWENO散逸が入りにくいようにする。面センサー値を二つの閾値間で正規化し、
三次smoothstepを適用する。

実装では発散と渦度をその最大絶対値で先に正規化してから二乗する。これにより
極めて小さい速度勾配でdenormalやunderflowフラグを不要に発生させない。

```text
x = clamp((S_face - onset) / (full - onset), 0, 1)
alpha = x^2 (3 - 2 x)
```

切替えを不連続な真偽値にせず、流束を連続的に変化させる。

## 5. 実装構造と拡張方法

CPU側は次の三層に分離している。

1. `mod_convective_keep`と`mod_convective_weno5z_roe`が面流束APIを提供する。
2. `mod_convective_leaf_registry`が方式名から面流束を選ぶ。
3. `mod_convective_hybrid`がセンサーと混合だけを担当する。

新しい構成流束を追加するときは、保存形の面流束APIを実装し、
`mod_convective_leaf_registry`へ一つの分岐を登録する。ハイブリッド本体の
ループや混合式を変更する必要はない。

新しいセンサーを追加するときは、`mod_convective_hybrid`のセンサー選択部と
入力検証へ追加する。構成流束の登録とは独立している。

CUDA側も`nse_cuda_hybrid.cuh`内で構成流束選択、センサー、混合カーネルを
分離しており、CPUと同じ設定値をC ABI経由で受け取る。

## 6. 制約

- 必要なghostセルは三層である。
- 保存変数は`[rho, rho*u, rho*v, rho*w, rho*E]`の5変数を前提とする。
- 現在の境界条件は周期境界である。
- WENO5-Z/Roeにpositivity-preserving limiterはまだ含まれない。
- 閾値を小さくすると頑健性と散逸が増え、大きくするとKEEP領域が増える。

## 7. 検証

自動テストでは次を確認する。

- 一定圧力の滑らかな場でハイブリッド流束がKEEP6と一致すること。
- 圧縮性圧力ジャンプで混合率が1となり、WENO5-Z/Roeと一致すること。
- 周期領域における面流束差分の総和が丸め誤差範囲で0となること。
- CPUとCUDAのSSPRK3更新結果が一致すること。
- ハイブリッド設定を入力ファイルから読み、実計算が正常終了すること。
