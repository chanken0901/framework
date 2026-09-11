# ツールの入口

通常の作業で主に使用するのは次の3コマンドです。最初のコマンドはFrameWork直下、
後の2つは生成した実行環境直下で実行します（Windowsは`python`、Linuxは`python3`）。

```text
python ScriptLibrary/RunEnvironment/prepare_environment.py <環境設計書.yaml>
python tools/run_case.py --prepare
python tools/run_case.py --build
```

実行は`python tools/run_case.py --run`、後処理は生成環境の`tools/postprocess_case.py`を使います。

| 場所 | 役割 |
| --- | --- |
| [RunEnvironment](RunEnvironment/README.md) | 主入口。環境生成、入力作成、実行、後処理 |
| [BuildSolver](BuildSolver/README.md) | 共通ビルド処理。通常はrun_caseから呼び出す |
| SetupCase | ケース生成と互換GPE生成の内部部品。環境生成からも使用する |
| [WorkspaceLinks](WorkspaceLinks/README.md) | ワークスペースのリンク設定 |
| [Git](Git/README.md) | 現行Git手順とモノレポ移行記録 |
| [legacy](legacy/README.md) | 旧リポジトリ・旧研究プロジェクト生成ツールの保管場所 |

`case_input.py`、`case_configuration.py`、`profile_selection.py`等は独立した
利用者向け入口ではなく、上記コマンドを支えるモジュールです。
複数の`yaml_support.py`は単体配布先でも動作するための部品なので、今回統合・削除していません。
テスト、manifest、caseテンプレート、拡張設定も実行環境生成に必要なため保持しています。
