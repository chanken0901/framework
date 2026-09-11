# 旧運用ツールの保管場所

2026-09-11に配置だけ整理しました。ファイルは削除していません。
現在のFrameWorkでは[RunEnvironment](../RunEnvironment/README.md)と標準Gitを使ってください。

| 旧配置 | 現配置 | 用途 |
| --- | --- | --- |
| `ScriptLibrary/Git/manage_*.py`、`library_repositories.yaml` | `Git/` | 旧分離リポジトリ運用。現行モノレポへ適用しない |
| `ScriptLibrary/SetupResearch/` | `SetupResearch/` | 旧schema-driven研究プロジェクト生成 |
| `ScriptLibrary/SetupCase/setup_solverlibrary_github.py` | `SetupCase/` | 旧SolverLibrary単独リポジトリ作成 |

移動に伴うSetupCaseへの参照は補正しました。旧Git/GitHubツールの変更操作は
検証のためにも実行していません。履歴再現以外の目的で使わないでください。
モノレポ移行スクリプトと記録は引き続き`ScriptLibrary/Git/`にあります。
