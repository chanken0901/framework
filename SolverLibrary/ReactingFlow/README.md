# ReactingFlow：独立反応流ライブラリ

状態：**R0基盤・R1熱力学と入力アダプター実装済み。詳細反応流ソルバーは未完成。**
従来NSEのStage 0〜10とは別の移行段階で管理する。

単成分NSE（熱揺らぎ拡張を含む）に依存しないライブラリとして再構築する。
既存の`NSE/src/extensions/multicomponent`は検証用・移行元として保持する。
現時点では既存ソルバーの物理的な移動、モデル名変更、実行環境の自動切替は行っていない。
NSE配下に同じ実装を複製して同期する方式にはしない。

## 構成

```text
ReactingFlow/
  src/reactingflow/    独立した機構データ契約（Python、NSE依存なし）
  tools/              データ検証の開発用入口
  examples/           実行可能な構造検証例（燃焼用機構ではない）
  tests/              回帰テスト
  docs/               移行計画・各段階の完了条件
```

R0では化学種・元素組成・反応物／生成物の係数を不変データに変換し、
元素保存を有理数で検査する。分子量は元素質量からkg/molで導出する。
未知のキー、誤った単位、重複、NaN、無効な係数は拒否する。
中性気体のみ。電離・表面反応は対象外。原子量は入力責任とし、自動補正しない。

この内部YAMLは**Cantera YAMLではない**。可逆フラグは保持するが逆反応速度はまだ計算しない。
反応速度・NASA物性・第三体等の未対応項目を黙って読み捨てることはしない。
R1の外部入力は別の`import_cantera`を使用する。対応範囲は下記のR1仕様を参照。

## 実行・検証

FrameWorkルートから、Windows（PowerShell）：

```powershell
python -m pip install PyYAML
python SolverLibrary/ReactingFlow/tools/validate_mechanism.py SolverLibrary/ReactingFlow/examples/topology_demo.yaml
python -m unittest discover -s SolverLibrary/ReactingFlow/tests
```

Linux（bash、利用中の仮想環境内）：

```bash
python3 -m pip install PyYAML
python3 SolverLibrary/ReactingFlow/tools/validate_mechanism.py SolverLibrary/ReactingFlow/examples/topology_demo.yaml
python3 -m unittest discover -s SolverLibrary/ReactingFlow/tests
```

`[OK] topology`は元素保存等の検査が通った意味であり、燃焼計算の検証済みを意味しない。
既存のケースは従来の`nse_multicomponent`で実行する。
新ライブラリの実行用manifestは、CFD入口を移植してから登録する。

詳細は[移行順序](docs/MIGRATION.md)を参照。

## R1：外部入力・熱力学

[R1仕様と実行手順](docs/THERMODYNAMICS.md)を追加した。
NASA-7/9物性、混合気体EOS、温度復元、基準量変換は独立Python実装。
入力アダプターと照合テストにのみCantera 3.2.0を使用する。
反応情報は保持するが、反応速度・化学時間発展・CFD・GPUはまだ新ライブラリに未実装。
