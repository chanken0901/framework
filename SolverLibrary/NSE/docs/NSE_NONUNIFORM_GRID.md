# 単成分NSEの不等間隔格子：実装状況と接続方針

2026-10-06。対象は単成分NSE。CPU、MPI、OpenMP、CUDA、MPI+CUDAの既存選択機能を
維持して直交不等間隔格子に対応させる作業。**現在は入力・格子生成まで。本計算は未対応。**
等間隔格子の生成・数値処理・環境選択は維持する。

## 入力・プレビュー（追加実装）

case.yamlで次を指定すると、環境の種類とは独立してinput.datへ変換される。
既存のnx/ny/nz、領域範囲などと同じ `grid` 内に記述する。

```yaml
grid:
  # nx/ny/nz、領域範囲などの既存項目も必要
  mapping:
    type: sinh
    strength: [0.0, 2.0, 2.0]
```

これは現時点では**幾何プレビュー専用**。本ソルバーに渡すと、
CPU/MPI/CUDA/MPI+CUDA共通の入力段階で「演算未接続」と明示して停止する。
solver側からこの停止を解除する入力スイッチはない。
`mapping` を省略すれば以前と同じinput.datを生成する。
`type: uniform` のstrengthは全て0に限る。
未知キー、3要素以外、非有限値、負値、20超、真偽値を拒否する。
GPEへsinhを転用する設定は拒否する。

伸長用に全域の軸情報から局所範囲のセル中心・面積・体積を生成する経路を追加した。
MPI通信を行う処理ではなく、各領域が同じ全域座標から切り出す構成。
現在の非周期ゴースト座標延長では物理端の幅を鏡映する。周期境界との接続は未完了。
sim%dx/dy/dzは伸長時には全域最小幅となるため、平均間隔とは区別する必要がある。

プレビューは各軸の左右境界・中心・幅をCSVへ出す（3D配列は作らない）。
ソルバー、MPIランチャー、GPU計算を起動せず、1プロセスで実行する。

Windows PowerShell、既存CPUビルドを使う例：

```powershell
cmake --build build/nse-positivity-cpu --target nse_grid_preview
.\build\nse-positivity-cpu\bin\nse_grid_preview.exe SolverLibrary/NSE/tests/input_grid_preview.dat build/grid_preview.csv
```

Linux、同じ名前でCPUビルドを構成済みの場合：

```bash
cmake --build build/nse-positivity-cpu --target nse_grid_preview
./build/nse-positivity-cpu/bin/nse_grid_preview SolverLibrary/NSE/tests/input_grid_preview.dat build/grid_preview.csv
```

生成環境ではサンプルの代わりに自分のinput.datを指定できる。
出力先の親フォルダが必要。既存CSVは上書きしない。
今回の確認はCPU構成のビルドのみ。実行例、Python回帰、MPI/CUDA実行は未検証。

## 今回追加したもの

`src/grid/mod_grid_axis.f90`（Fortran、MPI/CUDA依存なし）：

- `build_axis(edges,ng,periodic,axis)`：明示した昇順のセル境界座標から1方向の幾何を生成。
- `build_sinh_axis(n,ng,lower,upper,strength,periodic,axis)`：中心を細かくする対称伸長。
  強度0で等間隔。強度の範囲0〜20。極端な座標で境界が重なる場合は拒否。
- 物理境界0〜n、セル1〜nと、ng層のゴースト座標・セル幅・セル中心。
  周期は反対側の幅、非周期は端のセル幅列を鏡映して座標を延長する。
  これは幾何の延長のみで、物理境界条件を設定する処理ではない。
- 全域の最小セル幅。
- 7点の点値に対する1階・2階微分係数。CPU/GPUに共通の係数表を渡すためのデータ。
  係数は `(-3:3,global_cell)` の順に保持する。GPUへの転送・常駐処理自体は未実装。

伸長関数は、s=i/nとして
`x=lower+(upper-lower)*(1+sinh(strength*(2*s-1))/sinh(strength))/2`。
噴流中心が領域中央にない場合や片側だけの伸長は、明示座標で表す想定。
APIにはn>=ng>=3の制限がある。

微分係数は6次までの多項式点値の微分を再現するもの。
**有限体積のセル平均再構築係数でもKEEPの保存的流束係数でもない。**
任意の不等間隔格子で全演算の6次精度を保証するものではなく、
現行KEEP6/WENO5Z/CENTRAL6への単純置換は禁止する。
座標だけ変更して既存の等間隔用係数を流用する実装にはしない。

## 残作業（完了条件）

1. 共通入力：case.yaml→input.dat、テンプレートコメント、既定uniform、値検査は追加済み。
   実行環境生成を通した回帰確認は残る。
2. 格子と入出力：既存 `mod_grid_fvm` への接続、MPI局所範囲への切出し、
   SLF座標とメタデータ、初期化と境界の距離依存処理。伸長用の局所幾何生成までは追加済み。
3. 保存的演算：KEEP2/KEEP6、WENO5Z、ハイブリッド、粘性・熱流束・CFLを
   不等間隔用に整合させる。形式的な点微分の精度と保存性を別々に確認する。
4. GPU：同じ座標・係数をGPUへ常駐させ、CUDA/MPI+CUDAの演算を接続。
   既存CUDA-aware MPIとホスト経由の選択も保持する。
5. 拡張機能：揺らぎのセル体積・散逸離散化、FFTベースのHIT初期化・フォーシング、
   読込乱流の座標対応、後処理との整合を確認する。
   通常FFTは非等間隔の物理座標をそのまま扱えないため、格子写像や補間の定義が必要。
6. 検証：等間隔回帰、保存量、一様流保持、格子収束、境界、MPI分割数、CPU/GPU比較。

いずれかの環境を廃止したり、不等間隔指定を黙って等間隔へ戻したりしない。
上記が未完了である間は、不等間隔計算が可能になったとは案内しない。
入力のsinh指定はプレビューにのみ使用でき、未接続の計算へ入ることはない。

## 確認状況

- ソースを共通CMakeと環境生成用マニフェストに登録。
- 既存CMakeキャッシュでもソース一覧へ追加されるよう対応。
- CPU/OpenMP構成で係数生成の単体テスト実行ファイルまでビルド成功。
- テストコードは等間隔極限、中心伸長、周期ゴースト、多項式微分、重複座標拒否を含む。
- 実行テストは保留。MPI/CUDAのビルド・実機検証、数値スキームの検証は未実施。
