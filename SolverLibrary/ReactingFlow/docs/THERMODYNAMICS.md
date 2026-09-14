# R1 熱力学と入力

## 実装と単位

- NASA7、NASA9の複数温度区間、標準状態cp/h/s/g。
- 中性理想混合気体のcp/cv/h/e/s/g、密度、凍結音速、比熱比。
- エンタルピーには係数の生成エンタルピーを含める。反応熱の別加算は禁止。
- 混合エントロピーには組成と圧力の寄与を含む。ゼロ組成のlog(0)は評価しない。
- 質量分率の和は1（許容差1e-12）。自動正規化・負値クリップはしない。
- 温度外挿は拒否。エネルギーからの温度復元は共通有効温度区間内の二分法。
  係数の不連続等で収束できない場合はエラー。係数を自動で修正しない。
- 種別量：J/mol、J/(mol K)。混合量：J/kg、J/(kg K)。圧力Pa、温度K、密度kg/m³。
- `ReferenceScales`でrho_ref, U_ref, L_ref, T_refを明示。
  p_ref=rho_ref U_ref²、e_ref=U_ref²、t_ref=L_ref/U_ref。U_refを音速とは仮定しない。

## 外部機構

`import_cantera(path, phase=None)`はCantera 3.2.0で選択phaseを読み、
化学種・元素・NASA係数を独立した不変データに変換する。
未対応のthermoモデル、非理想気体、荷電種は拒否する。
反応物・生成物・可逆フラグをR0契約で検査し、反応形式を列挙する。
速度係数や輸送情報はCanteraが展開したcanonical YAMLに保持する。
R2以降は対応形式を独立した速度評価用データへコンパイルし、未対応形式を拒否する。
**速度評価は時間発展ではない。化学ODEの積分は未実装。** [R2仕様](KINETICS.md)参照。
選択phase以外の相をCFDへ取り込むものではない。

元ファイルのSHA-256と、展開済み機構のSHA-256を別々に記録する。
後者からは自動生成日時を除く。Cantera版等の変更ではハッシュが変わり得る。
元ファイルの外部参照先の内容は展開済み機構側に反映される。
CFD実行環境へのハッシュ登録は後段で実装する。

## Windows / Linux

FrameWorkルートで、Windowsは`python`、Linuxは`python3`（仮想環境推奨）：

```text
python -m pip install cantera==3.2.0 PyYAML
python -m unittest discover -s SolverLibrary/ReactingFlow/tests
python SolverLibrary/ReactingFlow/tools/inspect_thermo.py path/to/mechanism.yaml --temperature 1000 --pressure 101325 --mass-fractions '{"H2":0.1,"O2":0.9}'
```

PowerShellでもJSON引数を単一引用符で囲む。phase名が必要なら`--phase 名前`を指定。
出力は画面のJSONのみ。機構ファイルやケースを更新しない。
R2の瞬時反応速度も表示する場合は`--rates`を追加する。
`inspect_thermo.py`は開発段階の物性検証ツールで、計算のrun/postprocess入口ではない。

## 検証

Cantera 3.2.0同梱のh2o2.yaml（10種29反応）を入力検証・NASA7照合に使用。
NASA9はairNASA9.yamlから中性N2/O2/NOを選択して照合。
データファイルを本リポジトリへ複製せず、テスト時にCantera同梱資産を参照する。
Cantera未導入なら参照テストはskipされるため、全検証には上記依存が必要。
今回の検証ではskipなしで全テストを実施した。
cp/cv/h/e/s/g/密度/音速を350,800,1500,2800 Kで照合し、温度逆算も確認。
これは着火遅れ・燃焼速度・デトネーションの検証ではない。

式・標準状態の参照：[Cantera species thermodynamics](https://www.cantera.org/3.1/reference/thermo/species-thermo.html)。
計算時の入力アダプターは3.2.0で検証している。
