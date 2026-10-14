# SolverLibrary

FrameWorkモノレポ内でNSE・GPE・ReactingFlowを管理するディレクトリです。NSE/GPEのビルド条件は
各ソルバーの`solver_manifest.yaml`と`ScriptLibrary/BuildSolver`が管理し、計算ケースは
`ScriptLibrary/RunEnvironment`から`FrameWork`外の実行環境へ生成します。

## ソルバー

| ディレクトリ | 内容 | 入口 |
|---|---|---|
| `NSE` | 圧縮性Navier-Stokesソルバー | [`NSE/README.md`](NSE/README.md) |
| `GPE` | Gross-Pitaevskiiソルバー | [`GPE/gp3d/README.md`](GPE/gp3d/README.md) |
| `ReactingFlow` | 独立Fortran反応流ライブラリ。1D CFD等を実装、R4検証は未完了。2Dも開発中 | [状態・移行計画](ReactingFlow/README.md) |
| `Shared` | 複数ソルバーから利用する共通資産 | 各マニフェストから参照 |

## 標準運用

以下はNSE/GPEおよび旧NSE多成分系列の共通運用です。
ReactingFlowは専用CMake・実行入口を使用し、現段階で共通環境生成へ自動切替しません。
各入口と保存場所の区別は[リポジトリ整理方針](../docs/development/REPOSITORY_LAYOUT.md)を参照してください。

1. `SolverLibrary`と`ScriptLibrary`をローカル`FrameWork`で編集する。
2. `ScriptLibrary/RunEnvironment/environment.<model>.yaml`から外部実行環境を生成する。
3. 生成環境の`case.yaml`を編集し、`tools/run_case.py`で入力生成、検証、ビルド、実行を行う。
4. 確定した実装だけをローカル`FrameWork`とGitへ反映する。

NSEの最短手順、プロファイル、直接ビルドは
[`NSE/docs/usage/NSE_BUILD_AND_RUN.md`](NSE/docs/usage/NSE_BUILD_AND_RUN.md)を参照してください。

## 原則

- ソースやビルド設計の正本は`SolverLibrary`と`ScriptLibrary`に置く。
- 個別の計算条件と結果は外部実行環境で管理する。
- ビルド生成物はソースツリー外に出力し、Gitへ登録しない。
- NASは同期ミラーとして扱い、NAS上で直接編集・ビルド・実行しない。
