# 資料一覧

## 普段参照する手順書

- [学生向け導入・Git開発](guides/学生向け_FrameWork導入・Git開発手順書.md)
- [Git運用](guides/GIT運用マニュアル.md)
- [実行環境の生成・ビルド・実行](guides/外部実行環境_生成・ビルド・実行手順書.md)
- [ParaView後処理](guides/後処理_ParaView可視化手順書.md)
- [Windows／Linuxコマンド対応](WINDOWS_LINUX_COMMANDS.md)
- [NSE仕様書](../SolverLibrary/NSE/docs/NSE_SOLVER_SPECIFICATION.md)
- [GPE仕様書](guides/GPE_GP3Dコード全体仕様書.md)
- [実行環境の選択肢管理](../ScriptLibrary/RunEnvironment/ENVIRONMENT_OPTIONS.md)

## 保管資料

- `office/`：既存Word版7ファイル。内容は今回更新していないため、現行仕様はMarkdownとソースを確認する。
- `tables/`：ストレージ試算と実行環境選択肢台帳。内容は変更していない。
- `development/`：モノレポ移行記録。
- ソルバー固有の詳細仕様・検証手順は各ソルバーの`docs/`に保持する。

## 2026-09-11 配置整理

ファイル削除・機能削除は行っていない。ルートのMarkdown手順書は`guides/`、
Word版は`office/`、Excel資料は`tables/`へ移動した。Markdownの関連リンクを更新した。
Word・Excelの内容や埋め込みリンクは変更していない。古い本文中のルート配置は本一覧に読み替える。

旧運用ツールは`ScriptLibrary/legacy/`に保管し、現行の生成・ビルド・実行経路は保持した。
過去の`build-stage*-validation`6フォルダは`build/archived-validation/`へ移動した。
その中のCMakeキャッシュは移動前の絶対パスを含むため**検証記録として保管し、直接再ビルドしない**。
再利用が必要なら元のルート配置へ戻すか、新しいビルドフォルダでCMakeを構成し直す。
現在使用中の`build/nse-*`、`build/gpe-*`と`.codex-build/`は移動していない。

ソース、未コミットの開発成果、既存計算結果、参考文献、ショートカットは維持した。
配置のみの整理なのでディスク使用量は減らない。
