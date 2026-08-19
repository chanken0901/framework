# FrameWork Git運用マニュアル

更新日: 2026年8月19日  
対象: `chanken0901/framework`の管理者、教員、共同開発者、学生開発者  
運用形態: ScriptLibrary・SolverLibrary統合モノレポ

> 本書は2026年8月19日のモノレポ移行後の運用を定義します。旧版に記載されていた`ScriptLibrary`と`SolverLibrary`の個別clone、個別sync、個別snapshotは現行運用では使用しません。

## 1. 目的

本書は、研究用FrameWorkのソース、設定、テスト、文書を1つのGitHubリポジトリで安全に共同開発するための標準運用を定めます。

目的は次のとおりです。

- ScriptLibraryとSolverLibraryを整合した1つの版として管理する
- 学生を含む複数開発者の変更をPull Requestでレビューする
- 再現可能な計算版をtagで固定する
- ソースコードと計算データを分離する
- 誤操作、履歴分岐、force push、大容量ファイル混入を防ぐ
- Windows、Linux、計算サーバー間で同じ履歴を共有する

学生の初回導入と日常操作は`学生向け_FrameWork導入・Git開発手順書.md`を参照してください。本書は全体設計、権限、レビュー、リリース、例外対応を扱います。

## 2. 現行の正本とリポジトリ境界

### 2.1 正式リポジトリ

```text
https://github.com/chanken0901/framework.git
```

既定ブランチは`main`です。正式なソース、設定、テスト、Markdown文書はこのリポジトリの`main`を正本とします。

### 2.2 モノレポ構成

```text
FrameWork/
├─ .git/
├─ .gitignore
├─ README.md
├─ ScriptLibrary/
├─ SolverLibrary/
├─ docs/
└─ 仕様書・手順書
```

`.git`はルートに1つだけあります。`ScriptLibrary`と`SolverLibrary`は独立リポジトリではなく、同一コミット、同一ブランチ、同一Pull Requestの対象です。

### 2.3 旧リポジトリ

旧`ScriptLibrary.git`と旧`SolverLibrary.git`は、移行前履歴の参照とロールバック確認のため当面保持します。

- 新規開発をpushしない
- 学生のclone先として案内しない
- モノレポ動作確認と配布版tag確認後にArchiveする
- 削除はせず、移行先URLをREADMEに明記する

### 2.4 移行用スクリプトと旧管理スクリプト

`ScriptLibrary/Git/migrate_to_monorepo.ps1`は移行記録用です。通常運用では再実行しません。

次のスクリプトは旧分離リポジトリ用で、現行モノレポの日常操作には使用しません。

```text
ScriptLibrary/Git/manage_scriptlibrary.py
ScriptLibrary/Git/manage_solverlibrary.py
ScriptLibrary/Git/manage_library_repositories.py
ScriptLibrary/Git/library_repositories.yaml
```

現行運用ではルートの`FrameWork`に対して標準Gitコマンドを使用します。

## 3. 文書の正本

Gitで差分レビューできるMarkdown版を文書内容の正本とします。Word版は配布・印刷用の同期成果物です。

| 文書 | 正本 | 配布版 |
|---|---|---|
| 学生向け導入・開発手順 | `学生向け_FrameWork導入・Git開発手順書.md` | 同名`.docx` |
| Git全体運用 | `GIT運用マニュアル.md` | 同名`.docx` |
| ビルド・実行 | `外部実行環境_生成・ビルド・実行手順書.md` | 同名`.docx` |

同じ変更内でMarkdown版とWord版を更新します。不一致が見つかった場合は、Pull RequestでレビューできるMarkdown版を基にWord版を再生成します。

## 4. 権限と役割

### 4.1 リポジトリ管理者

- GitHub設定、アクセス権、ブランチ保護を管理する
- `main`への取り込み方針とリリースtagを管理する
- 旧リポジトリのArchiveを管理する
- 秘密情報や大容量ファイル混入などの事故対応を行う

### 4.2 メンテナー・教員

- Pull Requestをレビューする
- 数値計算、テスト、文書の整合性を確認する
- 破壊的変更と互換性変更を承認する
- 学生の権限と担当範囲を調整する

### 4.3 学生開発者

- 原則としてforkと作業ブランチを使用する
- 1つのPull Requestを1つの目的に限定する
- 必要なテストを実行し、結果を記載する
- 指摘へ同じ作業ブランチで対応する

### 4.4 計算利用者

- 配布tagまたは管理者が指定したコミットを使用する
- ソースツリー内へ計算結果を保存しない
- 変更が必要になった場合は作業ブランチへ移行する

## 5. ブランチ運用

### 5.1 `main`

`main`は常に共有可能な基準状態とします。

- 直接pushしない
- Pull Requestを経由する
- force pushしない
- マージ前に対象範囲のテスト結果を確認する
- リリースtagは`main`上の確認済みコミットへ付ける

### 5.2 作業ブランチ

```text
feature/<github-id>-<topic>
fix/<github-id>-<topic>
docs/<github-id>-<topic>
test/<github-id>-<topic>
refactor/<github-id>-<topic>
hotfix/<topic>
```

例:

```text
feature/tanaka-add-gpe-observable
fix/sato-periodic-corner
docs/yamada-student-guide
```

ブランチ名は半角英数字、小文字、ハイフンを基本とし、空白、日本語、個人情報を避けます。

### 5.3 ブランチの寿命

作業ブランチは短期間で完了させます。長期ブランチを作る場合も、定期的に`main`との差分とテストを確認します。マージ後はローカル・リモートの作業ブランチを削除できます。

## 6. 開発方式

### 6.1 fork方式（学生の標準）

```text
公式 upstream/main
        ↑ Pull Request
学生 origin/feature/...
        ↑ push
学生ローカル作業ブランチ
```

初回設定:

```powershell
git clone https://github.com/YOUR_GITHUB_ID/framework.git FrameWork
Set-Location .\FrameWork
git remote add upstream https://github.com/chanken0901/framework.git
git remote -v
```

### 6.2 公式リポジトリ直接方式

管理者が書き込み権限を付与した共同開発者だけが使用します。

```powershell
git clone https://github.com/chanken0901/framework.git FrameWork
```

この場合も作業ブランチとPull Requestを必須とします。

## 7. 標準変更フロー

### 7.1 開始前確認

fork方式:

```powershell
Set-Location C:\Research\FrameWork
git status
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git push origin main
```

公式リポジトリ直接方式:

```powershell
Set-Location C:\Research\FrameWork
git status
git switch main
git pull --ff-only origin main
```

作業ツリーに変更がある場合、pullやmergeを続けません。変更をコミット、一時退避、または意図を確認して個別に整理します。

### 7.2 ブランチ作成

```powershell
git switch -c fix/github-id-topic
git branch --show-current
```

### 7.3 編集と差分確認

```powershell
git status --short
git diff --stat
git diff
```

1つの変更がScriptLibraryとSolverLibraryの両方へ及ぶ場合は、同じブランチとPull Requestに含めます。これがモノレポ化の主要な利点です。

### 7.4 テスト

変更に最も近いテストから実行し、必要に応じて上位のスモークテスト、ビルド、計算比較へ広げます。

ScriptLibraryの代表的なPythonテスト:

```powershell
python -m unittest discover `
  -s .\ScriptLibrary\RunEnvironment\tests `
  -p "test_*.py"
```

NSE、GPE、CUDA、MPI、OpenMPの変更では、使用したプロファイル、コンパイラ、並列数、入力ファイル、結果をPull Requestへ記録します。

### 7.5 ステージ

変更ファイルを明示します。

```powershell
git add -- path/to/source path/to/test
git diff --cached --stat
git diff --cached
```

混在した作業ツリーで`git add -A`や`git add .`を使用しません。ステージ後に生成物、個人設定、不要な文書一時ファイルが入っていないことを確認します。

### 7.6 コミット

```powershell
git commit -m "KEEP6の周期境界精度テストを追加"
```

コミットメッセージは動詞を含む具体的な説明にします。

- `NSEハイブリッド流束の切替判定を修正`
- `GPE虚時間発展の計測ログを追加`
- `学生向けclone手順をモノレポ構成へ更新`

### 7.7 pushとPull Request

```powershell
$branch = git branch --show-current
git push -u origin $branch
```

Pull Requestのbaseは`chanken0901/framework`の`main`です。本文に次を記載します。

- 変更概要
- 背景と目的
- 実装上の重要点
- 利用者・数値結果・互換性への影響
- 実行したテストと結果
- 未確認事項、既知の制約

開発途中はDraft、レビュー可能になった時点でReady for reviewへ変更します。

## 8. Pull Requestの範囲

### 8.1 1 PR 1目的

機能追加、無関係な整形、文書整理を同じPull Requestへ混在させません。レビュー不能な大規模変更は、機能を壊さない単位へ分割します。

### 8.2 ソースとテスト

挙動変更には原則としてテストを含めます。テスト追加が困難な場合は、理由と手動確認方法をPull Requestへ記載します。

### 8.3 文書同期

入力形式、設定値、ビルド方法、出力仕様を変更した場合、対応するREADME、設計書、手順書を同じPull Requestで更新します。

### 8.4 数値計算変更

離散化、境界条件、時間積分、並列化、GPU処理を変更した場合は、最低限次を記録します。

- 使用した入力条件
- 格子点数、時間刻み、ステップ数
- 並列構成
- 基準結果との比較
- 保存量、誤差、安定性、性能への影響
- CPU/GPUまたは並列数を比較した場合の測定条件

## 9. レビュー基準

レビュアーは次を確認します。

- 変更目的と実装が一致している
- 意図しないファイルが含まれていない
- 生成物、計算結果、秘密情報がない
- 境界条件、添字範囲、並列領域、入出力が破綻していない
- テストが変更内容を検証している
- 文書が実装と一致している
- 後方互換性または移行方法が説明されている
- `main`へ取り込める単位になっている

重大な数値変更は、作成者以外の計算確認を推奨します。

## 10. マージ方式

| 変更 | 推奨方式 | 理由 |
|---|---|---|
| 通常の小規模PR | Squash and merge | 1 PRを1つの論理コミットとして残す |
| 複数の意味あるコミットを保持するPR | Create a merge commit | コミット単位の履歴を保持する |
| 外部履歴や別リポジトリを統合する移行 | Create a merge commit必須 | 親履歴を到達可能な状態で保持する |
| Rebase and merge | 原則使用しない | コミットIDが変わり、履歴追跡が難しくなる |

2026年8月19日のモノレポ移行PRは、ScriptLibraryとSolverLibraryの履歴保持のためCreate a merge commitで統合されています。

## 11. リリースとtag

### 11.1 tagの目的

tagは論文、発表、学生配布、再現計算で使用したソース版を固定します。通常開発は`main`、再現計算は指定tagを使用します。

現在の最初の学生配布版:

```text
v1.0.0-student
```

### 11.2 命名

原則としてSemantic Versioning形式を使用します。

```text
vMAJOR.MINOR.PATCH
```

- MAJOR: 入力・API・結果に大きな非互換変更
- MINOR: 後方互換な機能追加
- PATCH: バグ修正、文書修正、小変更

用途を明示する必要がある場合は`-student`、`-benchmark`などを付けます。

### 11.3 tag作成

tag作成は管理者がcleanな新規cloneまたは確認済み`main`で行います。

```powershell
git switch main
git pull --ff-only origin main
git status
git tag -a v1.0.1-student -m "Student release v1.0.1"
git push origin v1.0.1-student
```

既存tagを移動・上書きしません。誤りがある場合は新しいPATCH版を作成します。

## 12. GitHub設定の推奨

`main`へbranch protection ruleまたはrulesetを設定します。

- Pull Requestを必須にする
- force pushを禁止する
- branch deletionを禁止する
- 会話の解決を必須にする
- 少なくとも1名の承認を必須にする
- 将来CIを導入したら必要チェックを必須にする
- 管理者にも原則として同じ規則を適用する

現在CIが未設定の場合、存在しないチェックを必須にしません。先にCIを追加して安定稼働を確認してからrulesetへ組み込みます。

## 13. `.gitignore`とファイル配置

### 13.1 Git管理するもの

- ソースコード
- CMake、YAML、JSON、スキーマ
- 小規模な検証データ
- テスト
- Markdown、Word、Excelの正式文書

### 13.2 Git管理しないもの

- `build/`、`cmake-build-*`
- `CMakeFiles/`、`CMakeCache.txt`
- `*.o`、`*.obj`、`*.mod`、`*.smod`
- `*.exe`、`*.dll`、`*.so`、`*.a`、`*.lib`
- `__pycache__/`、`*.pyc`、`.pytest_cache/`
- `output/`、SLF、VTI、PVD、checkpoint
- `fort.*`、大容量ログ、一時ファイル
- Officeの`~$`一時ファイル
- 個人用ショートカット、秘密情報、ローカル専用設定

新しい生成物を発見した場合、個人のグローバルignoreだけで隠さず、チーム全体に不要ならルート`.gitignore`更新をPull Requestで提案します。

### 13.3 大容量データ

計算結果はNAS、研究データストレージ、オブジェクトストレージなどへ保存します。ソースリポジトリへ直接追加しません。Git LFSを導入する場合は、対象形式、容量、保存料金、取得方法を管理者が決定します。

## 14. NAS・共有フォルダ・複数端末

Git操作は各端末のローカルSSD上のcloneで行います。NASは計算データまたは配布ミラーとして使用し、Git作業コピーの正本にしません。

標準同期:

```text
PC Aのclone → commit/push → GitHub → fetch/pull → PC Bのclone
```

フォルダ同期ソフトで`.git`を双方向同期しません。複数端末で同じブランチを使う場合、端末を移る前にcommitとpushを完了します。

## 15. 例外・障害対応

### 15.1 working treeがcleanでない

```powershell
git status --short
git diff
git diff --cached
```

正式変更、作業途中、生成物に分類します。判断できない状態で一括削除や一括復元を行いません。

### 15.2 pullがfast-forwardできない

```powershell
git status
git log --oneline --graph --decorate --all -30
git remote -v
```

自動merge、自動rebase、force pushを行わず、どの変更を残すかレビュー担当者と決めます。

### 15.3 一時退避

```powershell
git stash push -u -m "同期前の一時退避"
git stash list
```

同期後:

```powershell
git stash pop
git status
```

競合した場合、stashを消す前にすべての変更が戻ったことを確認します。

### 15.4 誤コミット

push前でも、共有予定の履歴を書き換える前に担当者へ確認します。push済みコミットは原則としてrevertで打ち消します。

```powershell
git revert COMMIT_SHA
git push
```

`git reset --hard`とforce pushで共有履歴を作り直しません。

### 15.5 秘密情報をコミットした

直ちに管理者へ連絡し、該当資格情報を失効・再発行します。ファイルを次のコミットで削除するだけでは、過去履歴から秘密情報は消えません。履歴除去が必要な場合は管理者が影響範囲と手順を決定します。

### 15.6 大容量ファイルをpushした

新たなpushを止め、管理者へ連絡します。共有履歴の削除はforce pushを伴う可能性があるため、個人判断で実行しません。

## 16. `detected dubious ownership`

最初に、対象がローカルcloneか、NAS・別ユーザー所有コピーかを確認します。

```powershell
git rev-parse --show-toplevel
```

正しいローカルcloneであることを確認した場合だけ、その正確なパスを登録します。

```powershell
git config --global --add safe.directory C:/Research/FrameWork
```

`safe.directory "*"`は使用しません。

## 17. 開発終了時チェックリスト

- [ ] `git status --short`を確認した
- [ ] 意図したファイルだけをステージした
- [ ] `git diff --cached`を確認した
- [ ] 対象範囲のテストが成功した
- [ ] コミットメッセージが具体的である
- [ ] 作業ブランチをpushした
- [ ] Pull Requestのbaseが公式`main`である
- [ ] 数値変更の条件と結果を記録した
- [ ] 関連文書を更新した
- [ ] 秘密情報と大容量生成物が含まれていない

## 18. 管理者のリリースチェックリスト

- [ ] Pull Requestが承認・マージ済みである
- [ ] `main`を新規cloneして確認した
- [ ] ScriptLibraryとSolverLibraryの必須テストが成功した
- [ ] 対象ソルバーを少なくとも1プロファイルでビルドした
- [ ] サンプル計算と出力確認が成功した
- [ ] 手順書と実装が一致している
- [ ] 作業ツリーがcleanである
- [ ] tag名とリリース説明を確認した
- [ ] tagをGitHubへpushした
- [ ] 学生へclone URL、tag、既知の制約を案内した

## 19. 禁止事項

- `main`への直接push
- `git push --force`または`git push -f`
- 通常作業での`git reset --hard`
- 通常作業での`git clean -fd`
- 旧ScriptLibrary・SolverLibraryへの新規開発push
- モノレポ移行スクリプトの再実行
- 理由を確認しない履歴書き換え
- NAS上の共有コピーでの直接開発
- 秘密情報、計算結果、ビルド生成物のコミット
- 既存tagの付け替え

## 20. 標準コマンド一覧

fork開発者:

```powershell
git status
git fetch upstream --prune
git switch main
git merge --ff-only upstream/main
git push origin main
git switch -c feature/github-id-topic

# 編集・テスト

git status --short
git diff
git add -- path/to/file
git diff --cached
git commit -m "具体的な変更内容"
$branch = git branch --show-current
git push -u origin $branch
```

公式リポジトリ直接開発者:

```powershell
git status
git switch main
git pull --ff-only origin main
git switch -c fix/github-id-topic

# 編集・テスト・コミット

$branch = git branch --show-current
git push -u origin $branch
```

マージ後:

```powershell
git switch main
git pull --ff-only origin main
git branch -d BRANCH_NAME
```

## 21. 相談時に必要な情報

```powershell
git status
git branch -vv
git remote -v
git log --oneline --graph --decorate --all -30
git --version
```

秘密情報を除いた上で、実行コマンド、エラー全文、期待結果、実際の結果、OS、ブランチ、対象コミットを添えて管理者へ連絡します。

## 22. 公式参考資料

- Git公式: `https://git-scm.com/`
- GitHub clone手順: `https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository`
- GitHub fork手順: `https://docs.github.com/en/pull-requests/how-tos/work-with-forks/fork-a-repo`
- GitHub認証: `https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/about-authentication-to-github`
- GitHub ruleset: `https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets`
