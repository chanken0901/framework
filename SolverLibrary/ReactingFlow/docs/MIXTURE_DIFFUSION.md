# R4k: 二成分係数を指定する混合平均拡散

独立ReactingFlowのCPU逐次・1Dソルバーに追加。計算・入力処理はFortran。
既存の `none`、`constant`、`species_constant` は変更せず、
`transport_model='mixture_averaged'` を選択した場合のみ有効になる。

## 入力

`&flow1d` に `binary_diffusivities(i,j)` を指定する。単位は m²/s、
添字は機構ファイルの種順。全要素が必要で、対角は0、非対角は正、対称成分には同じ値を指定する。
共通 `mass_diffusivity` は0、`species_diffusivities` は省略する。
例えば**3種の機構の場合のみ**、次の設定になる（数値は検証用で物性値ではない）。

```fortran
transport_model='mixture_averaged',
binary_diffusivities(:,1)=0,     1.e-4, 2.e-4,
binary_diffusivities(:,2)=1.e-4, 0,     3.e-4,
binary_diffusivities(:,3)=2.e-4, 3.e-4, 0,
```

10種の機構なら10×10を指定する。省略、非対称、負値、他の拡散モデルとの重複は入力エラー。
使用した種対ごとの係数は出力CSVのコメントに記録する。
粘性・熱伝導の既存モデル、化学反応、各境界条件と併用できる。

10種の実行例は `examples/flow1d_h2_mixture_diffusion.in`。
既存の手順でh2o2機構を `h2o2.rf` に変換した後、以下を実行する。
出力CSVは未作成のファイル名にする。例の係数は物理検証には使用しない。

```powershell
.\build\reactingflow-fortran-release\rf_flow1d.exe h2o2.rf SolverLibrary\ReactingFlow\examples\flow1d_h2_mixture_diffusion.in mixture.csv
```

```bash
./build/reactingflow-fortran-release/rf_flow1d h2o2.rf SolverLibrary/ReactingFlow/examples/flow1d_h2_mixture_diffusion.in mixture.csv
```

## R4l: 化学種名付き物性ファイル

行列の代わりに `binary_diffusion_file='h2_synthetic_binary.rf'` を指定できる。
`transport_model='mixture_averaged'` は必要で、`binary_diffusivities` との同時指定は禁止。
ファイル名は入力 `.in` ファイルのあるディレクトリを基準に解決する。絶対パスも可。
これは既存の粘性用 `transport_file` とは別の入力であり、両方を併用できる。

書式例（H2、O2、N2の3種だけを持つ機構の場合。係数は架空の検証値）：

```text
RF_BINARY_DIFFUSION_V1
3
H2 0.002016
O2 0.031998
N2 0.028014
H2 O2 1.e-4
H2 N2 2.e-4
O2 N2 3.e-4
```

1. バージョン識別子、化学種数をそれぞれ1行に書く。
2. 全化学種について `種名 モル質量[kg/mol]` を書く。機構との相対差は1e-8以内。
3. 異なる全種対について `種名 種名 Dij[m²/s]` を各1行書く。N種ならN(N−1)/2行。

種名は大文字小文字を区別する。種一覧と種対一覧の内部の行順は自由。
`H2 O2` と `O2 H2` は同じ種対で、両方書くと重複エラーになる。
自己対角は書かず、内部で0に設定する。空行と `#` / `!` コメントを使用できる。
欠落・重複・未知の種・質量不一致・非正値・余分な行は読込み時に拒否する。
出力CSVには解決されたファイルパスと、種名付きの全係数を記録する。

そのまま実行する例は `examples/flow1d_h2_binary_file.in`。
同じディレクトリの `h2_synthetic_binary.rf` も一緒に配置する。
上記のWindows/Linuxコマンドの `.in` ファイル名だけをこの例に置き換えればよい。
ファイル読込みも計算もFortranで、Pythonは必要ない（機構変換・検証は従来どおり別工程）。
この入力拡張は定数Dijの読込みであり、Cantera/CHEMKIN形式の直接読込みや詳細物性モデルではない。

## R4n: 温度テーブル入力

`binary_diffusion_model='tabulated'` と `binary_diffusion_file` を指定する。
テーブルは基準圧力でのDijを保持し、評価時には `pref/p` を掛ける。
定数モデル用ファイルとは識別子が異なり、誤ったモデルでの読込みは拒否する。

```text
RF_BINARY_TABLE_V1
3
2 101325
H2 0.002016
O2 0.031998
N2 0.028014
300
H2 O2 1.e-4
H2 N2 2.e-4
O2 N2 3.e-4
1200
H2 O2 8.e-4
H2 N2 1.e-3
O2 N2 2.e-3
```

これは3種だけの機構用の架空の例。3行目は「温度点数 基準圧力[Pa]」。
種一覧の後に「温度[K]の1行＋全種対の係数[m²/s]」を温度点数分並べる。
各温度で全種対が必要だが、種対の順序・向きは自由。温度は正・厳密な昇順で2点以上。
補間は log(T)–log(Dij) の区分線形補間で、種対ごとの温度依存を表現する。
評価面温度が表の範囲外なら、外挿や端値への丸めをせず明示的に停止する。
この停止は時間刻み縮小による再試行の対象ではない。初期・境界状態と予想される燃焼温度を
十分に含むテーブルを、物性モデルの有効範囲内で用意すること。
べき乗モデル用の3パラメータや全輸送共通の温度補正との併用は禁止。
粘性・熱伝導には影響しない。CSVコメントには入力パス、基準圧力、全温度と全係数を記録する。
既存の `binary_diffusivity_種名_種名` 行は最初の温度点の係数を表す。
時間刻みの上限評価には全温度点の最大Dijと圧力下限を使うため、温度依存が非単調でも過小評価を避ける。

今回追加したのは物性テーブルを利用する機構であり、衝突積分の計算や
Cantera/CHEMKINからのテーブル生成、実在混合気での妥当性検証は未実装。

## 流束の計算式

質量分率をY、モル分率をX、モル質量をMとすると、

\[
\bar M=(\sum_iY_i/M_i)^{-1},\quad X_i=Y_i\bar M/M_i,
\qquad D'_{im}=\frac{\sum_{j\ne i}Y_j}{\sum_{j\ne i}X_j/D_{ij}}.
\]

\[
J_i^*=-\rho\frac{M_i}{\bar M}D'_{im}\partial_xX_i,
\qquad J_i=J_i^*-Y_i\sum_jJ_j^*.
\]

これは[Canteraの混合平均拡散流束の定義](https://www.cantera.org/stable/reference/onedim/governing-equations.html#diffusive-fluxes)
と同じ形式。ただしCanteraの物性評価全体を実装したものではない。
内部面は質量分率の算術平均と中心差分を用い、モル分率勾配は
`grad(X_i)=Mbar/M_i*grad(Y_i)-X_i*Mbar*sum(grad(Y_j)/M_j)` で評価する。
Dirichlet面は指定組成と半セル幅を用いる。純物質面では支配種の未定義な未補正流束を0とし、
補正後の支配種流束を他種の合計の負値として得る。単一種は拡散流束0。
最後の種の流束を閉じて総質量流束を0にし、エネルギー流束には種エンタルピー輸送も含める。

時間刻みには `2*max(Dij)*r*(1+r)`（rは最大／最小モル質量比）の保守的な
拡散作用素ノルム上限を使用する。軽い種と重い種を併用すると非常に小さい刻みになる場合がある。
これは非線形計算全体の正値性保証ではなく、既存の段棄却・再試行も継続する。

## 範囲と残作業

### R4m: 温度・圧力補正（選択式）

`mixture_averaged` の行列入力・ファイル入力のどちらにも以下を追加できる。

```fortran
binary_diffusion_model='power_law',
binary_reference_temperature=300,
binary_reference_pressure=101325,
binary_temperature_exponent=1.75,
```

係数を `Dij(T,p)=Dij_ref*(T/Tref)**n*(pref/p)` として扱う。
TはK、pはPa、Dはm²/s。指数nは全種対共通・非負で、適用する物性に合わせて指定する。
上記1.75は設定例で、対象混合気への適合を保証する値ではない。
内部面の平均密度・温度・組成から理想気体式で圧力を求め、境界では指定状態を使用する。
種エンタルピー輸送も補正されるが、粘性と熱伝導は変更しない。
既存の全輸送共通 `transport_temperature_model='power_law'` とは併用禁止。
既定値 `binary_diffusion_model='constant'` は従来の結果を維持する。
基準温度・圧力は正値で、指数を含め3項目とも明示指定が必要。
時間刻みには最大温度と、最小密度・最小温度・最大モル質量から得る圧力下限を使う。
この保守的な制限は小さい刻みになる場合があり、非線形安定性の保証ではない。
これは経験的な補正モデルであり、衝突積分・Soret・圧力勾配による拡散は未実装。

- 定数Dijが既定。R4mでは上記の温度・圧力補正を選択できる。
- 衝突積分からの温度・圧力依存Dij、Soret、圧力拡散、完全なStefan–Maxwell多成分拡散は未実装。
- 二成分不等モル質量でのFick則一致、純物質面、三成分の独立勾配評価、入力拒否、周期計算の保存量を試験する。
- 実在混合気の火炎速度・デトネーションを検証済みという意味ではない。

R4全体の残作業は [R4_REMAINING.md](R4_REMAINING.md) を参照。
