# ライブラリGit運用

ローカルFrameWork上の`ScriptLibrary`と`SolverLibrary`を、それぞれ独立して
安全に操作します。NASはGit操作後に同期するミラーであり、Gitの作業場所にはしません。
設定形式と安全検査だけを`manage_library_repositories.py`で共有します。

- `manage_scriptlibrary.py`: ScriptLibrary専用
- `manage_solverlibrary.py`: SolverLibrary専用

## 設定

`library_repositories.yaml`の`framework_root`をローカルFrameWorkの配置に合わせます。
現在の標準値は次です。

```yaml
framework_root: C:/Users/Owner/Documents/Codex/FrameWork
```

各リポジトリに
`origin`が設定済みなら`remote_url`は空で構いません。未設定の場合はGitHubのURLを指定します。

```yaml
remote_url: https://github.com/USER/ScriptLibrary.git
```

SSHを使用する場合:

```yaml
remote_url: git@github.com:USER/ScriptLibrary.git
```

## 初回確認

```powershell
python .\manage_scriptlibrary.py doctor
python .\manage_solverlibrary.py doctor
python .\manage_scriptlibrary.py status
python .\manage_solverlibrary.py status
```

Gitリポジトリまたは`origin`が未設定の場合、最初に接続します。

`doctor`は、ローカルの`.git`と`origin`だけでなく、GitHubへ到達できるか、指定した
ブランチがGitHub側に存在するかも`git ls-remote`で確認します。

```powershell
python .\manage_scriptlibrary.py connect
python .\manage_scriptlibrary.py connect --apply

python .\manage_solverlibrary.py connect
python .\manage_solverlibrary.py connect --apply
```

GitHub側に既存ファイルとコミットがあり、ローカル側が`not a Git repository`の場合は、
`connect`ではなく`adopt`を使用します。

```powershell
python .\manage_scriptlibrary.py adopt
python .\manage_scriptlibrary.py adopt --apply
```

`adopt`はGitHubの履歴をローカル側へ関連付けます。既存ローカルファイルのSHA-256を
前後で検証し、同名ファイルは上書きしません。GitHubにだけ存在するファイルは
ローカル側へ補完します。
処理後の差分はコミットせず、`status`で確認できる状態にします。

DryRunで表示されたURLが既存の`origin`と異なる場合は自動変更しません。内容を確認後、
`--update-remote --apply`を明示した場合だけ変更します。

## 日常操作

リモートの更新を取得し、fast-forwardだけを許可して同期します。

```powershell
python .\manage_scriptlibrary.py sync
python .\manage_scriptlibrary.py sync --apply

python .\manage_solverlibrary.py sync
python .\manage_solverlibrary.py sync --apply
```

各ライブラリの変更を別々のコミットとしてpushします。

```powershell
python .\manage_scriptlibrary.py snapshot --message "ケース生成スクリプトを更新"
python .\manage_scriptlibrary.py snapshot --message "ケース生成スクリプトを更新" --apply

python .\manage_solverlibrary.py snapshot --message "GPEモジュールを更新"
python .\manage_solverlibrary.py snapshot --message "GPEモジュールを更新" --apply
```

## 安全機構

- 書き込み操作は`--apply`を付けるまでDryRun
- pullは`--ff-only`のみ
- 自動merge、自動rebase、reset、force pushは実行しない
- ローカルとリモートが分岐していれば停止
- 未コミット変更がある状態でのpullを拒否
- 既定で50 MiBを超える変更ファイルのコミットを拒否
- 一方が失敗しても他方の結果を表示し、最後に失敗件数を返す

複数端末で作業する場合は、各端末にローカルFrameWorkを置き、GitHubを介して同期します。
作業前に`sync --apply`、作業後に`snapshot --apply`を行い、未コミット変更を残したまま
NAS同期を実行しないでください。

GitHub側に既存コミットがあるのにローカル側に`.git`履歴がない場合、`connect`後のsnapshotは
停止します。この場合は既存GitHubリポジトリを別フォルダへcloneし、現在のローカルファイルとの
差分を確認してからclone側へ変更を移してください。別履歴をforce pushしてはいけません。
