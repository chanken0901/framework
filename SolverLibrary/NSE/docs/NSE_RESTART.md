# 単成分NSEの再スタート

`restart.file`を指定すると、ghostなしの可搬SLFから5保存量、時刻、
ステップ番号を復元する。通常の初期条件生成（衝撃波・乱流の再配置も含む）は
実行しない。指定しない場合の従来動作は変更しない。

対象は等間隔格子の単成分NSE。CPU/OpenMP、MPI、CUDA、MPI+CUDAの開始処理に
接続している。MPI各rankは自分のy-z領域だけを読むので元のMPI並列数に依存しない。
不等間隔格子、反応流、多成分系のチェックポイントではない。

2026-10-08以降は、SLFに記録された実効gamma/Re/Pr/rho0/machも自動継承する。
読込み先YAMLのReや、HIT目標値から再計算したReより元データの値を優先する。
forcing・境界・数値方式は読込み先の選択を維持するため、必要な設定は引き続き記述する。
元の粘性あり／なしと異なる設定は拒否する。
旧SLFは実行時input.datを明示して再変換する必要がある。
詳細・Windows/Linuxの移行コマンドは
[物理パラメータの自動継承](NSE_IMPORTED_TURBULENCE.md)を参照する。

## 1. 保存出力を集約

新しい読込み先実行環境を用意し、そのルートで実行する。
元出力のmeta.jsonと、指定ステップの全rankのSLFが必要。
単一GPU出力も同じ変換を行う。変換は元の時刻・stepを保持する。

PowerShell（ケース番号を置換）:

```powershell
$sourceOutput = 'C:\Users\Owner\ResearchRuns\nse_case0025\cases\case0025\output'
New-Item -ItemType Directory -Force cases/case0026/initial_data | Out-Null
python ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py $sourceOutput --meta "$sourceOutput/meta.json" --step 1000 --output cases/case0026/initial_data/restart.slf
```

Linux/bash:

```bash
source_output="$HOME/ResearchRuns/nse_case0025/cases/case0025/output"
mkdir -p cases/case0026/initial_data
python3 ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py "$source_output" --meta "$source_output/meta.json" --step 1000 --output cases/case0026/initial_data/restart.slf
```

`--step latest`も使用可能。必ず変換成功を確認する。

## 2. case.yaml

元計算のphysics、flow、boundary、forcing、拡張機能の設定を維持し、追加する:

```yaml
restart:
  file: initial_data/restart.slf
```

相対パスはcase.yamlのディレクトリ基準。時刻・stepを手入力する項目はなく、
SLFヘッダーから取得する。元と同じ格子数・全方向の領域範囲が必須。
格子の補間、周期反復、速度オフセットは行わない。ghostは境界処理で再構築する。

`time.nsteps`と`time.t_max`は**通算の終了条件**。両方とも保存step・時刻より
大きくする。例: step=1000、t=0.1から再開してstep=2000まで進めるには
nsteps=2000とし、t_maxも十分大きくする。先に到達した条件で終了する。
固定dtは現在の設定、可変dtは再開した場から計算する。累積実行時間は0から計測。

出力先は新しいディレクトリを使う。既存meta.jsonがある出力先は拒否する。
write_initial=trueなら保存時刻・step番号で初期出力を作成する。

## 3. 入力生成と実行

PowerShell:

```powershell
python ./tools/run_case.py --prepare
python ./tools/run_case.py --validate-only
python ./tools/run_case.py --build
python ./tools/run_case.py --run
```

Linux/bashでは各行の`python`を`python3`に置換する。
`Restart loaded: step=..., time=...`を確認する。
既存の生成済み環境には変更が自動反映されないので、最新ライブラリから生成する。

## 再現性と検証範囲

- gamma、無次元化、輸送係数、境界条件、seed等はSLFに完全保存されない。
  同じ計算の継続には元のcase.yamlと拡張設定を保持する。物性の一致を自動保証しない。
- 揺らぎのカウンター型乱数は保存stepと現在のseedで継続する。CUDAにも復元stepを渡す。
  異なるバックエンド・並列数でのビット一致は保証しない。
- Petersen–Livescu forcingは読込み場から再評価する。評価回数による診断ログ間隔はリセットされる。
- ステップ境界の保存データが対象。RK段階途中の状態、FFT作業配列は復元しない。
- 非有限値、負の時刻/step、非物理的保存量、格子不一致、途中で切れたデータは拒否する。
- CPU/CUDA等のビルド確認と読込み単体テストは、マルチGPU継続計算の実機検証とは異なる。
