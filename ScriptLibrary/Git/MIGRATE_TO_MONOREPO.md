# モノレポ移行記録（移行完了）

> この移行はWindows PowerShellスクリプトで実施済みです。日常のWindows／Linux Git運用には本スクリプトを再実行せず、[`../../docs/WINDOWS_LINUX_COMMANDS.md`](../../docs/WINDOWS_LINUX_COMMANDS.md)とルートのGit運用マニュアルを使用してください。

> 2026年8月19日に`framework.git`への移行は完了し、Pull Request #1として`main`へマージされました。本書と`migrate_to_monorepo.ps1`は履歴確認・再現用です。現行cloneや学生環境で再実行しないでください。日常運用は[Git運用マニュアル](../../docs/guides/GIT運用マニュアル.md)を参照してください。

`migrate_to_monorepo.ps1`は、次の3リポジトリを1つの`framework.git`へ統合するために使用したスクリプトです。

- `framework.git`
- 旧`ScriptLibrary.git`
- 旧`SolverLibrary.git`

既存の`FrameWork`を直接変更せず、新しい移行用cloneを作成します。ディレクトリ削除、reset、rebase、force pushは行いません。

取込処理にはGit標準の`merge`と`read-tree`を使用しました。`git subtree`や`git filter-repo`などの追加コマンドは使用していません。

## 実施済みの移行結果

- 移行日: 2026年8月19日
- 移行先: `https://github.com/chanken0901/framework.git`
- 移行Pull Request: `#1 Migrate ScriptLibrary and SolverLibrary into monorepo`
- 移行元FrameWork: `8acaa476aab60e9f7c2f24a48c4aa43284ce30e0`
- 移行元ScriptLibrary: `2960c9efb5ac1de41d97c1a01ff02835a474b880`
- 移行元SolverLibrary: `f525d15d8ecf1b34bfd470df589520bfbae90f74`
- 移行ブランチ最終コミット: `0bb1602179588e558d3e0f58a849cbe7480333c2`
- `main`の移行マージコミット: `54b1b494462d2fc02d7d0f8cd196ca4787d477d1`
- 移行前安全tag: `pre-monorepo-2026-08-19`
- 最初の学生配布tag: `v1.0.0-student`

ScriptLibraryとSolverLibraryは、それぞれの完全なコミットグラフをmerge parentとして取り込み、`read-tree`で同名プレフィックスの下へ配置しました。ルートにだけ`.git`があり、下位ディレクトリにネストした`.git`はありません。

## 現行運用

```text
FrameWork/
├─ .git/
├─ ScriptLibrary/
├─ SolverLibrary/
└─ docs/
```

現行開発では、ルートの`FrameWork`に対して1つのブランチ、コミット、Pull Requestを使用します。旧リポジトリ用の`../legacy/Git/manage_scriptlibrary.py`、`../legacy/Git/manage_solverlibrary.py`は使用しません。

## 旧リポジトリの扱い

旧`ScriptLibrary.git`と旧`SolverLibrary.git`は、移行前履歴の参照用として当面保持し、新規開発には使用しません。モノレポの動作確認と配布版確認後、READMEに移行先を記載してArchiveします。履歴確認が不要になるまで削除しません。

## 以下は移行当時の再現情報

この節は障害調査や移行設計の参照用です。通常は実行しません。

### DryRun

```powershell
Set-Location C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\Git
powershell -ExecutionPolicy Bypass -File .\migrate_to_monorepo.ps1
```

DryRunではファイル、tag、ブランチ、リモートを変更せず、当時の3リポジトリについて次を検査しました。

- `main`ブランチにいること
- 未コミット変更がないこと
- `origin`が想定したGitHub URLであること
- `HEAD`と`origin/main`が一致すること
- 移行先ディレクトリが存在しないこと

### 実行時に使用した主要オプション

```powershell
powershell -ExecutionPolicy Bypass -File .\migrate_to_monorepo.ps1 `
  -Apply `
  -CreateSafetyTags `
  -PushSafetyTags `
  -RunSmokeTests
```

| オプション | 当時の意味 |
|---|---|
| `-Apply` | 実際に移行する。省略時はDryRun |
| `-CreateSafetyTags` | 3リポジトリに移行前tagを作る |
| `-PushSafetyTags` | 安全tagをGitHubへpushする |
| `-RunSmokeTests` | PythonテストとGPE生成DryRunを実行する |
| `-PushBranch` | 移行ブランチをGitHubへpushする |
| `-DestinationRoot PATH` | 新しい移行先を指定する |

## 注意

- 現行`main`へ対して再移行しない
- 旧リポジトリを現行開発へ戻さない
- 移行履歴保持のため、移行コミットをsquashやrebaseで作り直さない
- `pre-monorepo-2026-08-19`と`v1.0.0-student`を付け替えない
- force pushで`main`の履歴を書き換えない

