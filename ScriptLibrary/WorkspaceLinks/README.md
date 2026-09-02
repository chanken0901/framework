# ResearchRuns統合・リンク作成ツール

> このツールの既存リンク移行機能はWindows向けです。Linuxでの標準配置とbashコマンドは[`../../docs/WINDOWS_LINUX_COMMANDS.md`](../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。Linuxでは通常、`$HOME/ResearchRuns`と`$HOME/Research/FrameWork`を直接使用し、Windows互換ジャンクションは作成しません。

`C:\ResearchDesigns`に置いていた設計書を実行環境のルートへ統合し、
実行環境からローカルFrameWorkを参照できるようにします。

実行後の標準構成は次のとおりです。

```text
C:\Users\Owner\ResearchRuns\
  ├─ Designs\       # 設計書の実体
  ├─ FrameWork\     -> C:\Users\Owner\Documents\Codex\FrameWork
  ├─ gpe_case0001\
  ├─ gpe_case0002\
  └─ ...

C:\ResearchDesigns\ -> C:\Users\Owner\ResearchRuns\Designs
```

`C:\ResearchDesigns`は互換リンクとして残るため、以前のコマンドもそのまま
使用できます。FrameWork側にはResearchRunsへの逆リンクを作らないため、
再帰的なフォルダ走査は発生しません。

## 安全性

- `C:\ResearchDesigns`だけが存在する場合は、フォルダごと`Designs`へ移動します。
- 両方に内容がある場合は、同名項目がなければ統合します。
- 同名項目がある場合は、上書きせずエラーで停止します。
- 通常のフォルダがリンク予定位置を占有している場合も停止します。
- `--remove`はリンクだけを解除し、設計書、ソース、計算結果を削除しません。

## 設定

設定の正本は`workspace_links.yaml`です。

```yaml
schema_version: 1
execution_root: ${USERPROFILE}/ResearchRuns
link_type: auto

design_environment:
  legacy_path: C:/ResearchDesigns
  integrated_path: ${USERPROFILE}/ResearchRuns/Designs
  keep_legacy_link: true

links:
  framework:
    name: FrameWork
    target: ${USERPROFILE}/Documents/Codex/FrameWork
```

`link_type: auto`はWindowsでディレクトリ・ジャンクションを使用します。
Windowsのローカルフォルダでは通常、管理者権限は不要です。

## 操作

PowerShellでツールの場所を変数へ入れます。

```powershell
$tool = "C:\Users\Owner\Documents\Codex\FrameWork\ScriptLibrary\WorkspaceLinks"
```

まず現在の状態を確認します。

```powershell
python "$tool\setup_workspace_links.py" --status
```

移動・作成予定を確認します。この段階では変更しません。

```powershell
python "$tool\setup_workspace_links.py" --apply --dry-run
```

表示内容に問題がなければ、設計環境の統合とリンク作成を実行します。

```powershell
python "$tool\setup_workspace_links.py" --apply
```

統合後は次のパスを使用できます。

```powershell
Set-Location "$env:USERPROFILE\ResearchRuns"
Set-Location "$env:USERPROFILE\ResearchRuns\Designs"
Set-Location "$env:USERPROFILE\ResearchRuns\FrameWork"
```

接続先を設計書で変更した場合は、リンクだけを安全に張り替えます。

```powershell
python "$tool\setup_workspace_links.py" --apply --replace-links --dry-run
python "$tool\setup_workspace_links.py" --apply --replace-links
```

リンクを解除する場合は次を使います。`ResearchRuns\Designs`の実体は残ります。

```powershell
python "$tool\setup_workspace_links.py" --remove --dry-run
python "$tool\setup_workspace_links.py" --remove
```

## Linux（bash）の標準配置

`setup_workspace_links.py`の旧`C:\ResearchDesigns`移行処理はWindows専用なので、
Linuxでは実行しません。次のように実体ディレクトリを作成し、必要な場合だけ
FrameWorkへのシンボリックリンクを追加します。

```bash
mkdir -p "$HOME/ResearchRuns/Designs"
mkdir -p "$HOME/Research"

# FrameWorkは事前に$HOME/Research/FrameWorkへcloneしておく
if [ ! -e "$HOME/ResearchRuns/FrameWork" ]; then
  ln -s "$HOME/Research/FrameWork" "$HOME/ResearchRuns/FrameWork"
fi

ls -ld "$HOME/ResearchRuns/Designs" "$HOME/ResearchRuns/FrameWork"
readlink -f "$HOME/ResearchRuns/FrameWork"
```

リンクだけを解除する場合は、対象がシンボリックリンクであることを確認してから
次を実行します。リンク先のFrameWork本体は削除されません。

```bash
test -L "$HOME/ResearchRuns/FrameWork" && rm "$HOME/ResearchRuns/FrameWork"
```
