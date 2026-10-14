# 不等間隔SLFの後処理

## 実行方法

新しいFrameWorkから生成・更新した実行環境のルートで実行する。
保存形式は `output.format: slf` のままでよい。`NSEGRID1`格子情報を含む新しいSLFが対象。
等間隔格子の従来の操作・VTI出力は変更しない。

Windows PowerShell:

```powershell
python .\tools\postprocess_case.py --task paraview --steps all
python .\tools\postprocess_case.py --task statistics --steps all
```

Linux:

```bash
python3 ./tools/postprocess_case.py --task paraview --steps all
python3 ./tools/postprocess_case.py --task statistics --steps all
```

既定入力は `<case>/output`。別の出力先を使っている場合は、各コマンドに
`--input-dir cases/caseXXXX/output_run2` を追加する（実行環境ルート基準）。
可視化出力は `<case>/paraview/collection.pvd`、統計は
`<case>/statistics/turbulence_statistics.csv` と同名の `_metadata.json`。
`--task all` は不等間隔SLFでは可視化と統計を実行し、通常FFTをスキップした旨を表示する。
`--task fft` の明示指定は引き続き拒否する。補間やNUFFTは行わない。

## 可視化

ParaViewで `collection.pvd` を開く。不等間隔の場合、各時刻の実体は
座標付きRectilinearGrid (`.vtr`)。等間隔の場合は従来の `.vti`。
SLF本体は変更しない。MPIの全rankファイルを同じフォルダに置く。
ゴーストセルを除去し、SLF内の局所範囲で全領域を結合する。
欠損・重複・領域重複・格子/時刻/物理定数の不一致は拒否する。
rankファイル名・ステップ番号は保存時のまま維持する。

`--fields rho,u,v,w,p` 等の変数選択、`--steps latest`、時刻系列のPVDに対応。
不等間隔での `--stride` は1のみ。間引いた値を別のセル体積に置くことを避けるため、
2以上は停止する。格子情報付きSLFではgammaは保存値を使用し、CLI指定で上書きしない。
追加CellDataとして常に `cell_volume` と `integration_weight` を出力する。
KEEP6の状態は点値だが、VTRではソルバーの直接VTR出力と同じセル区画へ配置する。
これはセル平均への変換や高次積分を意味しない。

低レベルの変換・統計ツールはSLF内の情報だけでも動作する（meta.json不要）。
通常のpostprocess_case窓口では従来どおりcase.yamlとmeta.jsonを用意する。
大規模出力では1時刻の全領域配列をホストメモリへ結合するため、十分なRAMが必要。

## 統計の定義

物理セル体積を `V_i = Δx_i Δy_i Δz_i` とし、

```text
<f>_V = sum(V_i f_i) / sum(V_i)
u'_a = u_a - <u_a>_V
R_ab = <u'_a u'_b>_V
RMS_a = sqrt(R_aa)
k = (R_xx + R_yy + R_zz)/2
Mt = sqrt(2k) / <sqrt(gamma*p/rho)>_V
```

密度・圧力・音速・速度の平均、成分RMS、全速度RMS、k、相互応力を計算する。
これは全領域の**空間平均からの偏差**であり、時間平均・アンサンブル平均・Favre平均ではない。
噴流の平均せん断や非乱流域も含むので、kを局所の乱流強度と同一視しない。
TNTI条件付き統計や断面別の噴流統計は今回の対象外。

保存量Qについて、二つの積分を別の列に出力する。

- `volume_integral_* = sum(V_i Q_i)`：物理セル体積による積分近似。
- `solver_integral_* = sum(W_i Q_i)`：ソルバーと同じ重みでの積分。
  KEEP6系では `W_i = D6(x)_i D6(y)_i D6(z)_i`、それ以外は `W_i=V_i`。

`*` は mass / momentum_x / momentum_y / momentum_z / total_energy。
`physical_volume` と `solver_weight_sum` も出力する。
非周期境界では境界流束によって総量が変化するので、時間変化をただちに保存誤差と判断しない。

既存CSVとの共通列は残すが、不等間隔で未評価の粘性、散逸率、積分長、Taylor・Kolmogorov
スケール、Reスケール、スペクトル誤差は `nan` とする。既存のスペクトル計算を流用しない。
定義と未評価項目は付属metadata.jsonにも記録する。

## 配布・検証

新規スクリプトの追加配布は不要。既存の変換ツール・統計ツール・postprocess_caseを更新した。
生成済みResearchRuns環境は自動更新されないので、従来の手順で環境のソース/toolsを更新する。
SLF、計算結果、再スタート先は削除・変更しない。

`tests/test_mapped_slf_postprocess.py` でMPI結合、非等間隔座標、バイナリVTR、PVD、
物理体積重み、KEEP6重み、gamma継承、欠損・重複・切断ファイル、FFT拒否を検査する。

2026-10-13確認：後処理関連21テスト成功。周期・無反射それぞれについて、既存の
CPU/CUDA/4 MPI/4 MPI＋CUDA出力から求めた統計の一致を確認した。
MPI＋CUDAの元データは1GPU共有の試験であり、実複数GPUでの検証ではない。
環境生成の回帰テストには、今回変更していないforcingテンプレートの
target_dissipation記載を期待する既存テストが1件失敗する。後処理関連テストとは別項目。
