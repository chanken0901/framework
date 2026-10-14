# HITフォーシングの目標指定

## 方式ごとに設定をぶら下げて選択する

`target.type` の1行だけで切り替える。両方の設定を残してよく、選択した枝だけを読み込む。

```yaml
forcing:
  type: petersen_livescu
  petersen_livescu:
    target:
      type: direct  # direct または mach_reynolds
      direct:
        dissipation: 0.1
      mach_reynolds:
        turbulent_mach_number: 0.3
        taylor_reynolds_number: 100.0
        # mean_density: 1.0
        # mean_pressure: 0.7142857142857143
```

`direct` は散逸率直接指定、`mach_reynolds` は乱流パラメータからの換算。
選択された枝がない場合や、必須値が欠ける場合はエラー。
従来の平坦な `target_dissipation`、`target_mode` と `target_*` も互換入力として維持する。
旧形式で `target_mode` を省略した場合は直接指定。ただし新しい `target` と旧形式の併記はエラー。
背景密度の既定値は `physics.nse.rho0`、背景圧力はその密度/gamma。
異なる背景平均温度を使う場合は、その温度に対応した無次元平均圧力と密度を明示する。
これらは入力生成時の固定値であり、計算中の平均温度に追従しない。

1成分RMSを u、背景音速を a として、無次元量で

```text
a² = gamma*p/rho
u = Mt*a/sqrt(3)
nu = 1/(Reinput*rho)
epsilon_mass = 15*u^4/(nu*Re_lambda^2)
target_dissipation = rho*epsilon_mass
                   = (5/3)*Reinput*(gamma*p)^2*Mt^4/Re_lambda^2
```

等方性・低圧縮性を基準とした換算であり、強圧縮性での厳密な制御ではない。
膨張成分との配分は従来の `dilatational_ratio` を使用する。
実測の散逸率・Mt・Re_lambdaの定常平均が目標に一致する保証はない。

HIT初期化の粘性自動計算後のReinputを使い、換算自体は粘性を変更しない。
乱流読み込みでは読み込み先のReinputを使う。SLFから元の粘性を自動復元する機能ではない。
同じ物性で継続する場合、元の確定済みReinput、Prandtl数、比熱比、無次元化を維持する。
領域を長くしたという理由だけでReinputを再計算しない。

`run_case.py --prepare` による通常のinput.dat再生成で適用される。
Windows: `python .\tools\run_case.py --prepare`
Linux: `python3 ./tools/run_case.py --prepare`

生成されるFortran入力は従来と同じ `forcing_target_dissipation` の数値なので、
既存のCPU/MPI・CUDAフォーシング経路で共通に使用する。Fortranの時間発展処理は変更しない。
既存の実行環境はtoolsの古いコピーを持つ場合があるため、最新版から再生成または更新する。
直接input.datを編集する場合は従来どおり数値を指定する。
