# モノレポ移行スクリプト

`migrate_to_monorepo.ps1`は、次の3リポジトリを1つの`framework.git`へ統合します。

- `framework.git`
- `ScriptLibrary.git`
- `SolverLibrary.git`

既存の`FrameWork`を直接変更せず、新しい移行用cloneを作成します。ディレクトリ削除、
reset、rebase、force pushは行いません。

取込処理にはGit標準の`merge`と`read-tree`を使用します。`git subtree`や
`git filter-repo`などの追加コマンドは必要ありません。

## 1. 最初にDryRunする

PowerShellで次を実行します。

```powershell
Set-Location C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\Git

powershell -ExecutionPolicy Bypass -File .\migrate_to_monorepo.ps1
```

DryRunではファイル、タグ、ブランチ、リモートを変更しません。現在の各リポジトリについて、
次を検査します。

- `main`ブランチにいること
- 未コミット変更がないこと
- `origin`が想定したGitHub URLであること
- `HEAD`と`origin/main`が一致すること
- 移行先ディレクトリがまだ存在しないこと

親`FrameWork`で未追跡の`ScriptLibrary/`と`SolverLibrary/`だけは、現在の構成上必要なため
許可されます。それ以外の変更がある場合は停止します。

## 2. 表示された変更を整理する

各リポジトリで内容を確認してから、必要な変更だけをコミット・pushします。

```powershell
$frameworkRoot = "C:\Users\Owner\Documents\Codex\FrameWork"

git -C $frameworkRoot status
git -C "$frameworkRoot\ScriptLibrary" status
git -C "$frameworkRoot\SolverLibrary" status
```

この段階で親リポジトリから`git add ScriptLibrary`または`git add SolverLibrary`を
実行しないでください。

整理後にもう一度DryRunし、`Preflight passed.`と表示されることを確認します。

## 3. ローカル移行を実行する

安全タグを3リポジトリに作成し、GitHubへタグをpushしてから、ローカル移行を実行します。

```powershell
powershell -ExecutionPolicy Bypass -File .\migrate_to_monorepo.ps1 `
  -Apply `
  -CreateSafetyTags `
  -PushSafetyTags `
  -RunSmokeTests
```

既定の移行先は次です。

```text
C:\Users\Owner\Documents\Codex\FrameWork-monorepo-migration
```

移行先がすでに存在する場合、スクリプトは削除や上書きをせず停止します。再試行する場合は、
既存の移行先を手動で検査するか、別の移行先を指定します。

```powershell
.\migrate_to_monorepo.ps1 `
  -DestinationRoot C:\Users\Owner\Documents\Codex\FrameWork-monorepo-migration-2 `
  -Apply
```

## 4. 移行結果を確認する

```powershell
$monorepo = "C:\Users\Owner\Documents\Codex\FrameWork-monorepo-migration"

git -C $monorepo status
git -C $monorepo log --oneline --graph --decorate --all -30
git -C $monorepo diff main...chore/monorepo-migration --stat
```

次も確認します。

- `ScriptLibrary`と`SolverLibrary`が通常の追跡ディレクトリになっている
- 下位ディレクトリに`.git`がない
- Pythonテストが成功している
- GPEの実行環境生成DryRunが成功している
- 個人パス、ビルド成果物、計算結果が追加されていない
- `docs/development/MONOREPO_MIGRATION.md`のコミットIDが正しい

## 5. 移行ブランチをpushする

確認後、次を実行します。

```powershell
git -C $monorepo push -u origin chore/monorepo-migration
```

ローカル移行と同時にpushする場合だけ、最初の実行に`-PushBranch`を追加できます。

```powershell
.\migrate_to_monorepo.ps1 `
  -Apply `
  -CreateSafetyTags `
  -PushSafetyTags `
  -RunSmokeTests `
  -PushBranch
```

force pushは使用されません。

## 6. GitHubでPull Requestを作成する

次の向きでPull Requestを作成します。

```text
chore/monorepo-migration -> main
```

Pull Requestには次を記録します。

- 3リポジトリの移行元コミットID
- 安全タグ名
- 実行したテスト
- 移行後のディレクトリ構成
- 旧リポジトリをすぐ削除しないこと

## 7. マージ後の処理

Pull Requestをマージしてテストした後、初めてリリースタグを作ります。

```powershell
git clone https://github.com/chanken0901/framework.git FrameWork-release-check
Set-Location .\FrameWork-release-check

git tag -a v1.0.0-student -m "First monorepo release for student development"
git push origin v1.0.0-student
```

その後、旧`ScriptLibrary.git`と`SolverLibrary.git`のREADMEに移行先を書き、GitHub上で
読み取り専用としてArchiveします。モノレポの動作確認前に旧リポジトリを削除しないでください。

## オプション一覧

| オプション | 内容 |
|---|---|
| `-Apply` | 実際に移行する。省略時はDryRun |
| `-CreateSafetyTags` | 3リポジトリに移行前タグを作る |
| `-PushSafetyTags` | 安全タグをGitHubへpushする |
| `-RunSmokeTests` | PythonテストとGPE生成DryRunを実行する |
| `-PushBranch` | 移行ブランチをGitHubへpushする |
| `-DestinationRoot PATH` | 新しい移行先を指定する |
| `-SafetyTag NAME` | 安全タグ名を指定する |
| `-GitCommand PATH` | 使用するGit実行ファイルを指定する |
| `-PythonCommand PATH` | 使用するPython実行ファイルを指定する |

## 失敗時

スクリプトは既存環境を変更せず、新しい移行先にだけ処理します。途中で失敗した場合は、
エラーの対象を修正してから新しい`-DestinationRoot`で再実行します。既存移行先の自動削除、
`git reset --hard`、force pushは行わないでください。
