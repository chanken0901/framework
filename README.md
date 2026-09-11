# FrameWork

NSE（Navier-Stokes equations）およびGPE（Gross-Pitaevskii equation）のソルバー、ケース生成、ビルド、実行環境、後処理をまとめた研究用モノレポです。

## 最初に読む文書

文書全体は[資料一覧](docs/README.md)、実行するスクリプトは[ツール一覧](ScriptLibrary/README.md)から選べます。

### 学生・初回利用者

- [学生向け FrameWork導入・Git開発手順書](docs/guides/学生向け_FrameWork導入・Git開発手順書.md)
- [外部実行環境の生成・ビルド・実行](docs/guides/外部実行環境_生成・ビルド・実行手順書.md)
- [Windows／Linuxコマンド対応表](docs/WINDOWS_LINUX_COMMANDS.md)
- [case.yamlと拡張YAMLの使い方](ScriptLibrary/RunEnvironment/CASE_CONFIGURATION.md)

### 開発者・管理者

- [Git運用マニュアル](docs/guides/GIT運用マニュアル.md)
- [モノレポ移行記録](docs/development/MONOREPO_MIGRATION.md)

### ソルバー別

- [NSE](SolverLibrary/NSE/README.md)
- [NSEのビルドと実行](SolverLibrary/NSE/docs/NSE_BUILD_AND_RUN.md)
- [NSE多成分・反応流拡張](SolverLibrary/NSE/docs/NSE_MULTICOMPONENT_ROADMAP.md)
- [GPE](SolverLibrary/GPE/README.md)
- [GP3D](SolverLibrary/GPE/gp3d/README.md)
- [共通ビルドランナー](ScriptLibrary/BuildSolver/README.md)
- [実行環境生成](ScriptLibrary/RunEnvironment/README.md)

## リポジトリ構成

```text
FrameWork/
├─ ScriptLibrary/   ケース生成、ビルド、実行環境、補助ツール
├─ SolverLibrary/   NSE・GPEソルバー、設定、テスト
├─ docs/            guides（手順書）、office（Word版）、tables（台帳）、開発資料
└─ build/           ローカル生成物（Git管理外）
```

このリポジトリはモノレポです。`ScriptLibrary`と`SolverLibrary`は同じGit履歴で管理されています。旧`ScriptLibrary.git`、旧`SolverLibrary.git`は新規開発に使用しません。

## 安定版を取得する

### Windows（PowerShell）

```powershell
git clone https://github.com/chanken0901/framework.git FrameWork
Set-Location .\FrameWork
git switch --detach v1.0.0-student
```

### Linux（bash）

```bash
git clone https://github.com/chanken0901/framework.git FrameWork
cd FrameWork
git switch --detach v1.0.0-student
```

安定版タグは再現計算用です。開発する場合は`main`を最新化して作業ブランチを作成してください。

## 開発の原則

- `main`へ直接pushしない
- 作業ブランチとPull Requestを使用する
- 変更に対応するテストを実行する
- 計算結果、ビルド生成物、秘密情報をコミットしない
- ScriptLibraryとSolverLibraryを別々にclone・pushしない

詳細は[Git運用マニュアル](docs/guides/GIT運用マニュアル.md)を参照してください。

