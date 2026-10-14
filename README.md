# FrameWork

NSE、GPE、および開発中の独立反応流ライブラリReactingFlowと、ケース生成・ビルド・実行・後処理をまとめた研究用モノレポです。

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
- [NSEの用途別資料一覧](SolverLibrary/NSE/docs/README.md)
- [NSEのビルドと実行](SolverLibrary/NSE/docs/usage/NSE_BUILD_AND_RUN.md)
- [ReactingFlow：詳細反応機構の独立ライブラリ（開発中）](SolverLibrary/ReactingFlow/README.md)
- [旧NSE多成分系列：互換・比較用](SolverLibrary/NSE/docs/legacy-multicomponent/NSE_MULTICOMPONENT_ROADMAP.md)
- [GPE](SolverLibrary/GPE/README.md)
- [GP3D](SolverLibrary/GPE/gp3d/README.md)
- [共通ビルドランナー](ScriptLibrary/BuildSolver/README.md)
- [実行環境生成](ScriptLibrary/RunEnvironment/README.md)

## リポジトリ構成

```text
FrameWork/
├─ SolverLibrary/                 計算ライブラリの正本
│  ├─ NSE/                        単成分圧縮性NSE＋任意の熱揺らぎ
│  │  ├─ src/                    本体（旧多成分系列も互換用に保持）
│  │  ├─ docs/                   内容別の仕様・手順書
│  │  └─ tests/, examples/, tools/ 検証・入力例・後処理
│  ├─ GPE/gp3d/                  Gross–Pitaevskiiソルバー
│  ├─ ReactingFlow/              詳細反応機構の独立Fortranライブラリ
│  └─ Shared/                    共通資産
├─ ScriptLibrary/                 運用の正本
│  ├─ RunEnvironment/            設計YAML→環境生成→実行・後処理
│  ├─ BuildSolver/               共通ビルド処理
│  ├─ SetupCase/                 ケース生成の内部部品
│  └─ legacy/                    旧運用ツール（現行入口ではない）
├─ docs/                          共通資料
│  ├─ guides/                    導入・Git・実行・可視化
│  ├─ development/               構成・移行記録
│  └─ office/, tables/           既存Word版・台帳
└─ build/, .codex-build/          ローカル生成物（Git管理外）
```

### どのライブラリを使うか

| 対象 | 入口 | 位置付け |
|---|---|---|
| 単成分圧縮性流れ・衝撃波・乱流・熱揺らぎ | [NSE](SolverLibrary/NSE/README.md) | CPU/OpenMP・MPI・CUDA・MPI＋CUDA。組合せの制約は各仕様書を参照 |
| GPE・量子流体 | [GPE/GP3D](SolverLibrary/GPE/gp3d/README.md) | CPU/FFT・GPU/分散FFT。大規模実行前の実機確認は別途必要 |
| 詳細反応機構・反応流 | [ReactingFlow](SolverLibrary/ReactingFlow/README.md) | 独立ライブラリとして開発中。R4未完了、2Dも検証途上 |
| 従来の多成分Stage 0〜10 | [旧NSE多成分系列](SolverLibrary/NSE/docs/legacy-multicomponent/NSE_MULTICOMPONENT_ROADMAP.md) | 互換・比較用。ReactingFlowへ自動変換されるものではない |

熱揺らぎはNSEの任意拡張であり、ReactingFlowとは別です。
旧多成分Stage番号とReactingFlowのR番号を同じ進捗として扱わないでください。

### 設計・計算・解析の流れ

NSE/GPEと旧NSE多成分系列の標準窓口は次のとおりです。

```text
FrameWorkのソース＋環境設計YAML
  → prepare_environment.py：外部実行環境を生成
  → case.yaml／必要なconfig/*.yaml：計算条件を編集
  → tools/run_case.py：入力生成・ビルド・実行
  → tools/postprocess_case.py：共通後処理（NSEのFFTを含む）
  → tools/analyze_shock_turbulence.py：NSE衝撃波・乱流干渉の専用解析
```

計算条件・計算結果はResearchRuns等の外部実行環境に置き、ライブラリの正本と分けます。
ReactingFlowは現時点で専用CMake・実行入口を使用します。上の共通窓口へ自動切替しません。
手順の詳細は[環境生成](ScriptLibrary/RunEnvironment/README.md)と各ライブラリのREADMEを参照してください。

### 仕様書の探し方

- [NSE資料一覧](SolverLibrary/NSE/docs/README.md)：overview／usage／numerics／physics／parallel／postprocess／legacy-multicomponent／history。
- [ReactingFlow資料一覧](SolverLibrary/ReactingFlow/docs/README.md)：physics／numerics／flow／reference／validation／development。
- [共通資料一覧](docs/README.md)：導入、Git運用、OS別コマンド、Word版等。
- [GPE全体仕様](docs/guides/GPE_GP3Dコード全体仕様書.md)と[cuFFTMp実機テスト計画](SolverLibrary/GPE/gp3d/docs/cufftmp_test_plan.md)。

2026-10-14にNSE／ReactingFlowの57文書を内容別フォルダへ移しました。
本文の物理仕様はこの整理で変更していません。既存のWord版や外部生成環境は自動更新しません。
構成の管理方針は[整理方針](docs/development/REPOSITORY_LAYOUT.md)を参照してください。

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

