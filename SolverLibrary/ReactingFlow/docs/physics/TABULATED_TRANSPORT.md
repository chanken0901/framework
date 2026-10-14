# 実在機構の輸送データ取り込みとFortran評価

Canteraの混合平均輸送モデルから、純種粘性・純種熱伝導率・二成分拡散係数を
事前に温度テーブルへ変換する経路を追加した。
**時間発展、補間、混合則、拡散流束、安定刻み評価はすべてFortran。**
Python/Canteraは既存の機構変換と同じく前処理・独立比較だけに使用する。
従来の定係数、Sutherland/Wilke、Eucken/WMS、二成分拡散の各設定は維持する。

## 選択方法

```fortran
transport_model='mixture_averaged',
viscosity_model='tabulated_wilke', conductivity_model='tabulated_mix',
transport_file='transport/properties.rf',
binary_diffusion_model='tabulated', binary_diffusion_file='transport/binary.rf',
```

粘性・熱伝導のテーブル選択は上記の組で指定する。
粘性・熱伝導の定数値、Sutherlandの種別値、共通温度べき乗とは併用できない。
二成分拡散は既存の `mixture_averaged`＋`tabulated` を利用する。
ファイルパスは作業ディレクトリではなく、入力 `.in` ファイルからの相対パス。

## 数値モデル

純種のmu_i(T)、lambda_i(T)とD_ij(T,p_ref)を、log(T)–log(係数)の区分線形補間で評価する。
範囲外では外挿や端値への丸めをせず停止する。表の範囲は初期値だけでなく、
予想される反応後の温度・境界状態を含むように指定する。

粘性は既存Wilke混合則、熱伝導率は次式を使う。

```
lambda_mix = 0.5 [sum(X_i lambda_i) + 1/sum(X_i/lambda_i)]
D_ij(T,p) = D_ij(T,p_ref) p_ref/p
```

これはCanteraの混合平均輸送と同じ混合則であり、既存のEucken/WMS近似とは区別する。
[Canteraのモデル定義](https://www.cantera.org/3.1/cxx/d9/d17/classCantera_1_1MixTransport.html)。
種拡散はモル分率勾配と質量補正速度を使う既存方式、エネルギー流束には
種エンタルピー輸送を含める。[式と制限](MIXTURE_DIFFUSION.md)

粘性と熱伝導の時間刻み上限は全テーブル点の最大値を使う。
非単調な表でも最大温度点だけを見て過小評価しない。保守的なため、小さい刻みになり得る。
輸送モデルの非線形安定性や正値性を保証するものではない。

## 生成・実行（FrameWorkルート）

[REACTORS.md](REACTORS.md) に従ってビルドと前処理用Cantera3.2を準備する。
以下は未使用の出力名で実行する。生成先ディレクトリ・機構・CSVは上書きしない。
入力例のコピー先に既存入力がある場合は、別名にするか内容を確認してからコピーする。

Windows PowerShell:

```powershell
cmake --build build/reactingflow --parallel 4
$mechanism = python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')"
# h2o2.rfをすでに作成済みなら次の行は省略
python SolverLibrary/ReactingFlow/tools/export_mechanism.py "$mechanism" build/reactingflow/h2o2.rf
python SolverLibrary/ReactingFlow/tools/export_transport.py "$mechanism" build/reactingflow/transport --tmin 300 --tmax 3500 --points 201
Copy-Item SolverLibrary/ReactingFlow/examples/flow_tabulated_transport.in build/reactingflow/flow_tabulated_transport.in -Confirm
.\build\reactingflow\rf_flow1d.exe build/reactingflow/h2o2.rf build/reactingflow/flow_tabulated_transport.in build/reactingflow/flow_transport.csv
```

Linux bash:

```bash
cmake --build build/reactingflow --parallel 4
mechanism=$(python -c "import cantera; from pathlib import Path; print(Path(cantera.__file__).parent/'data/h2o2.yaml')")
# h2o2.rfをすでに作成済みなら次の行は省略
python SolverLibrary/ReactingFlow/tools/export_mechanism.py "$mechanism" build/reactingflow/h2o2.rf
python SolverLibrary/ReactingFlow/tools/export_transport.py "$mechanism" build/reactingflow/transport --tmin 300 --tmax 3500 --points 201
cp -i SolverLibrary/ReactingFlow/examples/flow_tabulated_transport.in build/reactingflow/flow_tabulated_transport.in
./build/reactingflow/rf_flow1d build/reactingflow/h2o2.rf build/reactingflow/flow_tabulated_transport.in build/reactingflow/flow_transport.csv
```

`--phase` はCantera相の選択、`--pressure` は拡散表の基準圧力[Pa]（既定101325）。
温度点数は201が既定。粗い表と細かい表で解・物性の変化を確認すること。
物性未定義の機構・非理想気体・無効な温度区間・非正係数は拒否する。
元機構の全種に輸送データが必要。生成前に同じ機構を化学反応用エクスポータで検証する。

## ファイル仕様と再現性

`properties.rf` は次の形式。種レコード順は自由だが重複・不足・未知種を拒否する。

```
RF_TRANSPORT_TABLE_V1
種数 温度点数
種名 モル質量[kg/mol]            # 全種
温度[K]
種名 粘性[Pa.s] 熱伝導率[W/(m.K)] # 全種、各温度で繰返す
```

`binary.rf` は既存の `RF_BINARY_TABLE_V1`。[詳細](MIXTURE_DIFFUSION.md)
名前とモル質量を機構と照合する。コメントにはCanteraバージョン、元ファイルと機構の
ハッシュを記録するが、Fortran読込みはハッシュの自動照合までは行わない。
**必ず同じ元機構から両方を生成すること。** 実行CSVにも使用した全テーブル値を記録する。

## 検証と未対応

- H2/O2系10種の300/733/1500/3500K、50/202.65kPa、複数組成・純物質をCanteraと比較。
- 粘性・熱伝導・Dijは相対8e-5以内、種拡散・エンタルピー輸送は相対1e-4の許容値で比較。
- 21点→201点の温度格子細分化、非単調表の上限評価、入力拒否、反応あり周期CFDの保存量を確認。

これでCanteraの詳細輸送データをテーブルとして使えるが、**Fortranによる衝突積分の直接評価ではない**。
低密度理想気体の混合平均モデルに限定し、Soret、圧力拡散、完全な多成分Stefan–Maxwell、
高圧補正、火炎速度・消炎・反応境界の実用検証は未対応／未完了。
CFD例は接続確認用の短時間計算で、火炎やデトネーションの妥当性検証ではない。
