# FrameWork Git関連ツール

> 現在のFrameWorkは`ScriptLibrary`と`SolverLibrary`を統合したモノレポです。日常のGit操作はFrameWorkルートで標準Gitコマンドを使用してください。

現行の手順書:

- `../../学生向け_FrameWork導入・Git開発手順書.md`
- `../../GIT運用マニュアル.md`

## 現行モノレポで使用する操作

FrameWorkルートへ移動します。

```powershell
Set-Location C:\Research\FrameWork
git rev-parse --show-toplevel
git status
```

作業開始時:

```powershell
git switch main
git pull --ff-only origin main
git switch -c feature/github-id-topic
```

変更後:

```powershell
git status --short
git diff
git add -- path/to/changed-file
git diff --cached
git commit -m "具体的な変更内容"
$branch = git branch --show-current
git push -u origin $branch
```

push後はGitHubで`main`宛てのPull Requestを作成します。

## 移行用ファイル

次のファイルは2026年8月19日に実施したモノレポ移行の記録・再現用です。学生の日常操作では使用しません。

```text
migrate_to_monorepo.ps1
MIGRATE_TO_MONOREPO.md
```

移行は完了しています。既存のFrameWorkや学生のcloneに対して`migrate_to_monorepo.ps1`を再実行しないでください。

## 旧分離リポジトリ用ファイル

次のファイルは、旧`ScriptLibrary.git`と旧`SolverLibrary.git`を別々に管理していた時期の互換・履歴確認用です。

```text
manage_scriptlibrary.py
manage_solverlibrary.py
manage_library_repositories.py
library_repositories.yaml
```

現行モノレポでは使用しません。実行すると、`ScriptLibrary`や`SolverLibrary`を独立Gitリポジトリとして扱おうとするため、現在の構成と一致しません。

## 安全原則

- `main`へ直接pushしない
- force pushを使用しない
- 自動的な履歴書き換えを行わない
- 作業前後に`git status`と`git diff`を確認する
- 計算結果とビルド生成物をGitへ追加しない
- NAS上の共有コピーでGit操作しない

