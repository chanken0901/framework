# ReactingFlow：独立反応流ライブラリ

状態：**R1熱力学・R2反応速度・R3の0次元BDF反応器の計算本体をFortran化済み。CFDは未実装。**
従来NSEのStage 0〜10とは別の移行段階で管理する。

単成分NSE（熱揺らぎ拡張を含む）に依存しないライブラリとして再構築する。
既存の`NSE/src/extensions/multicomponent`は検証用・移行元として保持する。
現時点では既存ソルバーの物理的な移動、モデル名変更、実行環境の自動切替は行っていない。
NSE配下に同じ実装を複製して同期する方式にはしない。

## 構成

```text
ReactingFlow/
  CMakeLists.txt      Fortranビルド
  src/fortran/        計算本体（NSE/Python依存なし）
  reference/python/   旧Python版：入力変換・比較検証用として保全
  tools/              機構入力変換・オフライン検証のみ
  examples/           構造検証例とFortran反応器入力例
  tests/              回帰テスト
  docs/               移行計画・各段階の完了条件
```

R0では化学種・元素組成・反応物／生成物の係数を不変データに変換し、
元素保存を有理数で検査する。分子量は元素質量からkg/molで導出する。
未知のキー、誤った単位、重複、NaN、無効な係数は拒否する。
中性気体のみ。電離・表面反応は対象外。原子量は入力責任とし、自動補正しない。

R0の内部YAMLは**Cantera YAMLではない**。構造だけを検証する形式で、速度評価には使わない。
反応速度・NASA物性・第三体等の未対応項目を黙って読み捨てることはしない。
R1の外部入力は別の`import_cantera`を使用する。対応範囲は下記のR1仕様を参照。

## 実行・検証

通常の計算は[Fortran版のビルド・実行手順](docs/REACTORS.md)を使用する。
以下は入力構造の検査のみ。以前のPython版は`reference/python/`へ移し、
Fortran計算時には呼び出さない。旧`run_reactor.py`は比較専用の
`run_reference_reactor.py`へ名称変更した。

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
NASA-7/9物性、混合気体EOS、温度復元、基準量変換はFortran実装。
以前のPython版はオフライン照合用に残す。
入力アダプターと照合テストにのみCantera 3.2.0を使用する。
R2の反応速度評価は[速度仕様](docs/KINETICS.md)を参照。
R3の定容・定圧断熱反応器、保存検査、着火比較は[反応器仕様と実行手順](docs/REACTORS.md)を参照。
CFD・MPI・GPUはまだ新ライブラリに未実装。
