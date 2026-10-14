# FrameWorkの構成と整理方針

更新：2026-10-14。利用者向け入口・開発履歴の分離に続き、NSE／ReactingFlowの
57文書を内容別フォルダへ移動した。物理モデル・実行コマンド・ソース配置・既存ケースは変更していない。

## ライブラリの区別

| 正本 | 役割 | 注意 |
|---|---|---|
| `SolverLibrary/NSE` | 単成分NSEと任意の熱揺らぎ拡張 | 不等間隔／境界の組合せはNSE資料一覧で確認 |
| `SolverLibrary/GPE/gp3d` | GPEソルバー | NSEとは独立 |
| `SolverLibrary/ReactingFlow` | 詳細反応機構の独立Fortranライブラリ | R4未完了。2Dも検証途上 |
| `SolverLibrary/NSE/src/extensions/multicomponent` | 既存Stage 0〜10系列 | 互換・比較用に保持。ReactingFlowと混同しない |
| `SolverLibrary/Shared` | マニフェストから参照する共通資産 | 単独で削除・移動しない |
| `ScriptLibrary/RunEnvironment` | 環境生成・実行・後処理の窓口 | スクリプト名と配置を維持 |
| `ScriptLibrary/BuildSolver` | 共通ビルド | 通常はrun_caseから呼ぶ |
| `ScriptLibrary/legacy` | 旧運用ツール | 現行の入口ではない |

## 保存場所

- `src/`：実装の正本。似た名前でも互換系列と新規系列を一括統合しない。
- 各ライブラリの`tests/`：回帰テストと入力。CMake・manifestから参照されるため本体と共に管理する。
- `examples/`・`case_templates/`：サンプルと環境生成用テンプレート。個別計算の結果を置かない。
- 各ライブラリの`docs/`：そのライブラリの仕様。NSE/ReactingFlowには用途別READMEを設けた。
- `docs/guides/`：共通運用手順。`docs/office/`は既存Word版の保管場所であり更新済みとは限らない。
- `docs/development/`：リポジトリ全体の運用・移行記録。
- `build/`・`.codex-build/`：既存のローカルビルド・検証成果物。今回移動・削除していない。
- 外部実行環境（ResearchRunsなど）：個別case・入力・計算結果。ライブラリ整理の対象外。

## 今回分離したもの

NSE/ReactingFlowのREADME冒頭の過去の追加記録を、それぞれの
NSEは`docs/history/DEVELOPMENT_HISTORY.md`、ReactingFlowは
`docs/development/DEVELOPMENT_HISTORY.md`へ移した。内容を捨てず、現状と履歴を区別した。
仕様書のファイル名は維持し、格納フォルダを変更した。リポジトリ内のリンクと
manifestのコピー元パス、生成READMEの案内を新配置へ更新した。
古い外部リンクは自動転送されない。各ライブラリの`docs/README.md`を入口にする。
既存の外部実行環境・Word本文・ブックマークは今回変更していない。
上位READMEの「不等間隔は未接続」「ReactingFlowはCFD未移植」という古い案内を修正した。

## 今回移動しないものと理由

ビルドディレクトリはCMakeキャッシュに絶対パスを持つため、そのまま移動すると
再ビルドできなくなる。ソース、テスト、Python内部モジュール、環境設計YAMLも
相対参照・コピー定義を持つため、見た目だけを理由に移動しない。
特に環境設計YAMLは配置を変えると参照先の解決も変わり得る。
今後物理配置を変更する場合は、manifest、CMake、import、テンプレート参照を同時更新し、
生成環境を含む回帰試験を行う。

ファイル削除・計算結果移動・Git設定変更・同期設定変更は行っていない。
したがって、この整理でファイル数やディスク使用量が減るわけではない。
