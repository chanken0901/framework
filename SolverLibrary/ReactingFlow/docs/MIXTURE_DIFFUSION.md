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

## 計算式

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

- 今回のDijは指定定数。既存の共通温度べき乗を明示選択した場合は従来どおり全輸送流束を倍率補正する。
- 衝突積分からの温度・圧力依存Dij、Soret、圧力拡散、完全なStefan–Maxwell多成分拡散は未実装。
- 二成分不等モル質量でのFick則一致、純物質面、三成分の独立勾配評価、入力拒否、周期計算の保存量を試験する。
- 実在混合気の火炎速度・デトネーションを検証済みという意味ではない。

R4全体の残作業は [R4_REMAINING.md](R4_REMAINING.md) を参照。
