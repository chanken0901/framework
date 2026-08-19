# 学生向け FrameWork導入・Git開発手順書

更新日: 2026年8月19日  
対象リポジトリ: `https://github.com/chanken0901/framework.git`  
対象: FrameWorkを用いて計算する学生、およびソースコード・設定・文書の開発に参加する学生

> 重要: 現在のFrameWorkはモノレポです。`ScriptLibrary`と`SolverLibrary`は、同じ`framework`リポジトリの中にあります。旧`ScriptLibrary.git`と旧`SolverLibrary.git`はクローンしません。

## 1. この手順書でできること

本書を上から順に実行すると、次の作業ができるようになります。

- 安定版を取得して計算に使用する
- 開発用のリポジトリを自分のPCへクローンする
- 自分専用の作業ブランチを作る
- 変更をコミットしてGitHubへ送る
- Pull Requestを作ってレビューを依頼する
- 教員や他の学生が更新した最新版を安全に取り込む
- よくあるGitエラーから復旧する

計算だけを行う人は第5章、開発にも参加する人は第6章以降を使用してください。

## 2. 現在のリポジトリ構成

クローン後の基本構成は次のとおりです。

```text
FrameWork/
├─ .git/                  Git管理情報。編集・削除しない
├─ .gitignore             生成物や計算結果をGitから除外する規則
├─ README.md              最初に読む案内
├─ ScriptLibrary/         ケース生成、ビルド、実行環境、Git補助資料
├─ SolverLibrary/         NSE・GPEソルバー本体、設定、テスト
├─ docs/                  開発・移行記録
├─ 学生向け_FrameWork導入・Git開発手順書.md
├─ GIT運用マニュアル.md
└─ 各種仕様書・実行手順書
```

Gitリポジトリは最上位の`FrameWork/.git`だけです。`ScriptLibrary`と`SolverLibrary`を別々にクローン、pull、commit、pushしてはいけません。

## 3. 用語を最初に理解する

| 用語 | 意味 |
|---|---|
| リポジトリ | Gitで履歴管理されるプロジェクト一式 |
| クローン | GitHubのリポジトリをPCへ複製する操作 |
| `main` | 全員で共有する基準ブランチ。直接編集・直接pushしない |
| ブランチ | 個人の変更を安全に分離する作業線 |
| コミット | 変更内容と説明を履歴として記録する操作 |
| push | ローカルのコミットをGitHubへ送る操作 |
| pull | GitHubの更新をローカルへ取り込む操作 |
| Pull Request | 作業ブランチを`main`へ取り込むためのレビュー依頼 |
| fork | 公式リポジトリを自分のGitHubアカウントへ複製する操作 |
| `origin` | 通常、自分がクローンしたGitHubリポジトリの名前 |
| `upstream` | fork運用で使う公式`chanken0901/framework`の名前 |
| tag | 再現可能な配布版を示す固定ラベル |

## 4. 作業方法を選ぶ

### 4.1 計算だけを行う

安定版タグ`v1.0.0-student`を使用します。ソースを変更せず、入力ファイルや計算結果はリポジトリ外の実行環境へ置きます。第5章へ進んでください。

### 4.2 開発にも参加する（推奨）

自分のGitHubアカウントへforkし、作業ブランチから公式リポジトリへPull Requestを送ります。学生が多い場合に最も安全な方法です。第6章へ進んでください。

### 4.3 公式リポジトリへの書き込み権限を持つ

公式リポジトリを直接クローンできますが、`main`へ直接pushせず、必ず作業ブランチとPull Requestを使用します。第7章へ進んでください。

## 5. 計算利用者: 安定版をクローンする

### 5.1 必要なもの

WindowsではGit for WindowsとPowerShellを用意します。LinuxではGitと端末を用意します。ソルバーのコンパイラ、MPI、CUDA、FFTWなどは、使用する計算方式に応じて別途必要です。

Git for Windowsの公式入手先:

```text
https://git-scm.com/downloads/win
```

Gitが使えるか確認します。

```powershell
git --version
```

バージョン番号が表示されればGitは利用できます。`git: command not found`または「認識されていません」と表示された場合は、Gitをインストールして端末を開き直してください。

### 5.2 Windowsで安定版を取得する

ローカルSSD上の短いパスを使用します。OneDrive、ネットワークドライブ、NASへ直接クローンしないでください。

```powershell
New-Item -ItemType Directory -Path C:\Research -Force
Set-Location C:\Research
git clone https://github.com/chanken0901/framework.git FrameWork
Set-Location C:\Research\FrameWork
git switch --detach v1.0.0-student
```

### 5.3 Linuxで安定版を取得する

```bash
mkdir -p ~/research
cd ~/research
git clone https://github.com/chanken0901/framework.git FrameWork
cd FrameWork
git switch --detach v1.0.0-student
```

### 5.4 取得結果を確認する

```powershell
git status
git describe --tags --always
git remote -v
```

期待する状態は次のとおりです。

- `git describe`に`v1.0.0-student`が表示される
- `ScriptLibrary`と`SolverLibrary`の両方が存在する
- `git status`に意図しない変更が表示されない

> 安定版タグでは`detached HEAD`と表示されます。計算利用だけなら正常です。この状態でソースを編集・コミットしないでください。開発する場合は第6章または第7章の方法で`main`から作業ブランチを作ります。

### 5.5 計算を始める前の確認

```powershell
Test-Path .\ScriptLibrary\BuildSolver\build_model.py
Test-Path .\SolverLibrary\NSE\CMakeLists.txt
Test-Path .\SolverLibrary\GPE\gp3d\CMakeLists.txt
```

3行とも`True`なら基本ファイルは取得できています。具体的なビルド・実行方法は次を参照してください。

- `外部実行環境_生成・ビルド・実行手順書.md`
- `ScriptLibrary/BuildSolver/README.md`
- `SolverLibrary/NSE/docs/NSE_BUILD_AND_RUN.md`
- `SolverLibrary/GPE/gp3d/README.md`

## 6. 開発参加者: forkしてクローンする（推奨）

### 6.1 GitHubアカウントを準備する

GitHubアカウントを作成し、教員へGitHubユーザー名を連絡します。大学・研究室の規則に従って二要素認証を設定してください。パスワード、アクセストークン、秘密鍵をチャット、YAML、ソースコードへ貼り付けてはいけません。

### 6.2 公式リポジトリをforkする

ブラウザで次を開きます。

```text
https://github.com/chanken0901/framework
```

画面右上の`Fork`を押し、自分のアカウントに`framework`を作成します。`Copy the main branch only`は通常オンのままで構いません。

GitHub公式のfork説明:

```text
https://docs.github.com/en/pull-requests/how-tos/work-with-forks/fork-a-repo
```

### 6.3 自分のforkをクローンする

次の`YOUR_GITHUB_ID`を自分のGitHubユーザー名に置き換えます。

```powershell
New-Item -ItemType Directory -Path C:\Research -Force
Set-Location C:\Research
git clone https://github.com/YOUR_GITHUB_ID/framework.git FrameWork
Set-Location C:\Research\FrameWork
```

すでに`C:\Research\FrameWork`が存在する場合、上書きせず、まずそのフォルダの`git status`と`git remote -v`を確認してください。

### 6.4 公式リポジトリを`upstream`として登録する

```powershell
git remote add upstream https://github.com/chanken0901/framework.git
git remote -v
```

期待する表示は概ね次のとおりです。

```text
origin    https://github.com/YOUR_GITHUB_ID/framework.git
upstream  https://github.com/chanken0901/framework.git
```

`origin`は自分のfork、`upstream`は公式リポジトリです。この対応を入れ替えないでください。

### 6.5 Gitの名前とメールアドレスを設定する

初回だけ実行します。

```powershell
git config --global user.name "あなたの氏名またはGitHub表示名"
git config --global user.email "GitHubに登録したメールアドレス"
git config --global --get user.name
git config --global --get user.email
```

メールアドレスを公開したくない場合は、GitHubが提供する`noreply`アドレスを使用できます。

### 6.6 初回状態を確認する

```powershell
git status
git branch --show-current
git remote -v
git log -1 --oneline
```

この時点では`main`にいて、作業ツリーがcleanであることを確認します。

## 7. 書き込み権限を持つ学生: 公式リポジトリをクローンする

教員から公式`framework`リポジトリへの書き込み権限を付与された場合だけ、この方法を使用します。

```powershell
New-Item -ItemType Directory -Path C:\Research -Force
Set-Location C:\Research
git clone https://github.com/chanken0901/framework.git FrameWork
Set-Location C:\Research\FrameWork
git remote -v
```

この場合、`origin`が公式リポジトリです。権限があっても`main`へ直接pushしません。必ず次章の作業ブランチを使用します。

## 8. 日常の標準開発手順

以下はfork運用を基準に説明します。公式リポジトリを直接クローンした人は、`upstream`を`origin`に読み替えてください。

### 8.1 作業開始前に`main`を最新化する

```powershell
Set-Location C:\Research\FrameWork
git status
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git push origin main
```

`git status`に変更がある場合は、最新化を続けず、第13章の「変更が残っている」を確認してください。

### 8.2 作業ブランチを作る

ブランチ名は半角英数字とハイフンを使います。

```powershell
git switch -c feature/github-id-short-description
```

例:

```powershell
git switch -c feature/tanaka-add-gpe-test
git switch -c fix/sato-nse-boundary
git switch -c docs/yamada-update-build-guide
```

| 接頭辞 | 用途 |
|---|---|
| `feature/` | 新機能、機能追加 |
| `fix/` | バグ修正 |
| `docs/` | 文書だけの変更 |
| `test/` | テスト追加・修正 |
| `refactor/` | 動作を変えない構造整理 |

### 8.3 ファイルを編集する

編集対象は通常、`ScriptLibrary`、`SolverLibrary`、ルート文書のいずれかです。計算結果、ビルド生成物、個人用パス、資格情報をリポジトリへ置かないでください。

編集途中でも定期的に状態を確認します。

```powershell
git status --short
git diff --stat
git diff
```

`git diff`を終了するには`q`を押します。

### 8.4 必要なテストを実行する

変更範囲に対応するテストを実行します。少なくとも、変更した機能を直接確認できるテストを1つ実行してください。

ScriptLibraryのPythonテスト例:

```powershell
python -m unittest discover `
  -s .\ScriptLibrary\RunEnvironment\tests `
  -p "test_*.py"
```

NSE・GPEのビルドやテストは、各READMEと`外部実行環境_生成・ビルド・実行手順書.md`に従います。実行したコマンドと結果はPull Requestへ記載します。

### 8.5 変更ファイルを明示してステージする

最初は`git add -A`や`git add .`を使わず、必要なファイルを明示します。

```powershell
git add -- .\SolverLibrary\NSE\src\変更したファイル.f90
git add -- .\SolverLibrary\NSE\tests\追加したテスト.f90
git diff --cached --stat
git diff --cached
```

Wordを編集した場合はWordを閉じ、`~$`で始まる一時ファイルが含まれていないことを確認します。

### 8.6 コミットする

```powershell
git commit -m "NSE境界条件のコーナー処理を修正"
```

良いコミットメッセージは「何を変更したか」が分かります。

- `GPE虚時間発展の時間計測出力を追加`
- `KEEP6精度検証テストを追加`
- `学生向けGit導入手順を更新`

`update`、`fix`、`変更`だけの説明は避けてください。

### 8.7 自分のGitHubへpushする

```powershell
$branch = git branch --show-current
git push -u origin $branch
```

初回pushではブラウザ認証を求められる場合があります。GitHubのパスワードを端末へ直接入力する方式は使用できません。Git Credential Manager、SSH、またはGitHub CLIの認証を使用します。

### 8.8 Pull Requestを作る

push後に表示されるURL、またはGitHubの`Compare & pull request`からPull Requestを作ります。

- base repository: `chanken0901/framework`
- base branch: `main`
- head repository: 自分のfork
- compare branch: 自分の作業ブランチ

Pull Request本文には次を記載します。

- 何を変更したか
- なぜ必要か
- 結果や利用者への影響
- 実行したテストと結果
- 未確認事項や制約

開発途中ならDraft Pull Requestにします。レビュー可能になったら`Ready for review`へ変更します。

## 9. レビュー中の修正

レビューコメントに対応するときは、同じ作業ブランチで修正します。

```powershell
git switch feature/github-id-short-description
git status
```

修正、テスト、ステージ、コミット、pushを繰り返します。

```powershell
git add -- path/to/changed-file
git diff --cached
git commit -m "レビュー指摘に基づき入力検証を追加"
git push
```

同じPull Requestが自動的に更新されます。新しいPull Requestを作り直す必要はありません。

## 10. Pull Requestがマージされた後

### 10.1 fork運用

```powershell
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git push origin main
git branch -d feature/github-id-short-description
```

GitHub上の作業ブランチも不要なら削除できます。

```powershell
git push origin --delete feature/github-id-short-description
```

削除対象のブランチ名を必ず確認してください。`main`は削除しません。

### 10.2 公式リポジトリを直接クローンした場合

```powershell
git switch main
git pull --ff-only origin main
git branch -d feature/github-id-short-description
```

## 11. 計算データとソースコードを分ける

GitHubへ登録するもの:

- Fortran、CUDA、Pythonなどのソース
- CMake、YAML、JSONなどの設定
- 小さなテスト入力と期待値
- README、仕様書、手順書

GitHubへ登録しないもの:

- `build/`、CMake生成物、`*.o`、`*.mod`、`*.exe`
- `__pycache__/`、`*.pyc`
- SLF、VTI、PVD、checkpoint、計算結果
- 大容量ログ、動画、画像列
- 個人の絶対パスを含むローカル設定
- パスワード、トークン、秘密鍵

計算ケースと結果は、RunEnvironmentで作った外部実行環境または研究データ用ストレージへ保存します。NAS上の共有ミラーをGitの作業コピーとして使用しません。

## 12. 複数PC・計算サーバーで使う

各PCやサーバーに個別にクローンし、GitHubを介して同期します。PC AのフォルダをPC Bへ上書きコピーして履歴を同期しないでください。

作業開始時:

```powershell
git status
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
```

作業終了時:

```powershell
git status --short
git push
```

未pushの作業を別PCへ引き継ぐ場合は、作業ブランチをGitHubへpushしてから別PCで取得します。

```powershell
git fetch origin --prune
git switch --track origin/feature/github-id-short-description
```

## 13. よくある問題と対処

### 13.1 `destination path already exists`

クローン先がすでに存在します。既存フォルダを削除せず、まず確認します。

```powershell
Get-ChildItem C:\Research\FrameWork -Force
git -C C:\Research\FrameWork status
git -C C:\Research\FrameWork remote -v
```

既存フォルダが必要か分からない場合は、教員または管理者へ確認してください。

### 13.2 `not a git repository`

現在位置がクローンしたFrameWork内ではありません。

```powershell
Set-Location C:\Research\FrameWork
git rev-parse --show-toplevel
```

### 13.3 `detached HEAD`

安定版タグを使用中なら正常です。計算用cloneをそのまま開発用へ変更せず、第6章に従って自分のforkから開発用cloneを用意する方法が安全です。

すでにforkからクローンした開発環境で、一時的にtagを確認していただけなら、次で`main`へ戻します。

```powershell
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git switch -c feature/github-id-short-description
```

### 13.4 変更が残っていて更新できない

```powershell
git status --short
git diff
```

正式な変更なら先にコミットします。作業途中なら一時退避できます。

```powershell
git stash push -u -m "main更新前の一時退避"
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git stash pop
```

`stash pop`で競合した場合は、自動削除や強制復元を行わず、`git status`に表示されたファイルを確認します。

### 13.5 `non-fast-forward`でpushできない

同じブランチがGitHub側で更新されています。force pushせず、次を保存して管理者へ相談します。

```powershell
git status
git log --oneline --graph --decorate --all -20
git remote -v
```

### 13.6 マージ競合が発生した

```powershell
git status
```

競合ファイル内の`<<<<<<<`、`=======`、`>>>>>>>`を確認し、残す内容を手作業で決めます。解決後に次を実行します。

```powershell
git add -- path/to/resolved-file
git status
git commit
git push
```

両方の意味を理解できない場合は、その場で止めてレビュー担当者へ相談してください。

### 13.7 `detected dubious ownership`

NASや別ユーザー所有のコピーを操作していないか確認します。正しいローカルcloneであることを確認した場合だけ、エラーに表示された正確なパスを登録します。

```powershell
git config --global --add safe.directory C:/Research/FrameWork
```

`safe.directory "*"`は使用しません。

### 13.8 誤ってファイルを変更した

未コミット変更を破棄する操作は元に戻せない場合があります。まず差分を確認します。

```powershell
git diff -- path/to/file
```

内容を確実に破棄してよい場合だけ実行します。

```powershell
git restore -- path/to/file
```

`git reset --hard`や`git clean -fd`を通常作業で使用しないでください。

### 13.9 大容量ファイルが含まれた

コミット前ならステージだけ解除します。

```powershell
git restore --staged -- path/to/large-file
git status --short
```

計算結果をリポジトリ外へ移し、必要なら`.gitignore`追加をPull Requestで提案します。すでにpushした場合は自分で履歴を書き換えず、管理者へ連絡してください。

## 14. 禁止事項

- `main`への直接push
- `git push --force`および`git push -f`
- 理由を確認しない`git reset --hard`
- 理由を確認しない`git clean -fd`
- 競合時に片方を機械的に上書きすること
- 旧`ScriptLibrary.git`または旧`SolverLibrary.git`への新規開発push
- `migrate_to_monorepo.ps1`の再実行
- NASや共有フォルダ上での直接開発
- パスワード、トークン、秘密鍵、個人情報のコミット
- 計算結果やビルド生成物のコミット

## 15. 初回導入チェックリスト

- [ ] Gitのバージョンを確認した
- [ ] 正式な`framework`または自分のforkをクローンした
- [ ] `ScriptLibrary`と`SolverLibrary`が同じclone内にある
- [ ] fork利用者は`upstream`を登録した
- [ ] `git status`がcleanである
- [ ] Gitの名前とメールアドレスを設定した
- [ ] 作業ブランチを作成した
- [ ] 計算結果をリポジトリ外へ保存する場所を決めた
- [ ] ビルド・実行手順書を確認した
- [ ] Pull Requestの提出先が`chanken0901/framework`の`main`であることを確認した

## 16. 日常作業の最短手順

fork利用者の標準手順です。

```powershell
Set-Location C:\Research\FrameWork
git status
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git push origin main
git switch -c feature/github-id-short-description

# 編集とテスト

git status --short
git diff
git add -- path/to/changed-file
git diff --cached
git commit -m "具体的な変更内容"
$branch = git branch --show-current
git push -u origin $branch
```

push後、GitHubで公式`main`宛てのPull Requestを作成してください。

## 17. 困ったときに提出する情報

相談時はスクリーンショットだけでなく、次の出力をテキストで提出してください。パスワードやトークンが表示されていないことを確認します。

```powershell
git status
git branch -vv
git remote -v
git log --oneline --graph --decorate --all -20
git --version
```

あわせて、実行したコマンド、期待した結果、実際の結果、使用OS、対象ブランチを伝えてください。

## 18. 公式参考資料

- Git for Windows: `https://git-scm.com/downloads/win`
- GitHubでリポジトリをクローンする: `https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository`
- GitHubでforkを作成・同期する: `https://docs.github.com/en/pull-requests/how-tos/work-with-forks/fork-a-repo`
- GitHub認証の概要: `https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/about-authentication-to-github`

