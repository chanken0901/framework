# 2D非一様初期条件

2026-10-02追加。Fortran `rf_flow2d` の `&flow2d` で `initial_mode` を選択する。
既定は従来の `uniform`。今回はビルド確認のみで、実行・物理検証は保留中。

| モード | 内容 | 追加項目 |
| --- | --- | --- |
| uniform | 全セル一様 | なし |
| split | 指定座標の小さい側だけ別状態 | split_axis=1または2、split_position |
| box | 指定矩形内だけ別状態 | region_box=xmin,xmax,ymin,ymax |
| profile | 同一格子のセルごとの任意分布を読む | initial_profile |

背景は `temperature,pressure,velocity(2),mass_fractions(ns)`。
split/boxの別状態は `region_temperature,region_pressure,region_velocity(2),region_y(ns)`。
温度・圧力・組成は明示する。組成の自動継承・正規化はしない。
座標は物理座標[m]、温度[K]、圧力[Pa]、速度[m/s]。

splitはセル重心座標が `split_position` **未満**のセルを別状態にする。
boxは `[xmin,xmax) × [ymin,ymax)` 内のセル重心を選ぶ。
セルを横切る界面の体積率平均は行わない。領域がセル重心を含まなければ背景のみとなる。
boxの高温部を指定する場合も、圧力・組成を含む初期状態はユーザーが決める。
自動的な定容加熱・圧力平衡・衝撃波跳躍条件の計算は行わない。

`examples/flow2d_split_h2air.in` は左右の圧力差を設定する粗格子の説明例。
検証済みの衝撃波管解やデトネーション条件ではない。
反応を有効にする場合は機構、状態、解像度、時間刻みの妥当性を別途確認する。

## 任意分布のファイル

`initial_mode='profile'` と `initial_profile='initial.rf'` を指定する。
uniform/split/box時にinitial_profileを併記したり、profile時に省略するとエラー。
profile時は背景の温度・圧力・組成ではなくファイル値を使用する。
相対パスは入力namelistファイルのフォルダ基準で、作業ディレクトリに依存しない。

```text
RF_INITIAL2D_V1
<mechanism canonical_hash>
<ns> <nx> <ny>
<x_node> <y_node>
... (nx+1)*(ny+1)行、iが先に進む順
1 <T> <p> <u> <v> <Y1> ... <Yns>
2 <T> <p> <u> <v> <Y1> ... <Yns>
... nx*ny行
```

機構ハッシュ、次元、全節点座標、セル番号と順序を照合する。
節点許容差は方向ごとに `1e-12*領域幅 + 64*機械epsilon*最大絶対座標`。
セル番号は `i+nx*(j-1)`。組成和は1、負値・非有限値は拒否する。
NASA温度範囲やEOSの検査は既存の状態変換と共通。
各レコードは指定数の数値だけを記載し、コメントや追加列は置かない。
格子間補間・SLF取込み・機構間変換は行わない。

## 初期分布の書出し

任意のモードで `initial_export='generated_initial.rf'` を指定すると、
計算開始前に同形式で初期場を書き出す。既存ファイルは上書きしない。
`end_time=0` と組み合わせれば時間発展をせず初期場を作成できる。
生成されたファイルのセル物理量を編集し、profileモードで読み込める。
これは基本変数形式の初期条件であり、時刻・保存変数を完全に維持する再始動ファイルではない。
出力CSVやSLFをそのままinitial_profileへ指定することはできない。

## 実行例

FrameWorkルート、機構ファイル `mechanism.rf` を別途用意した場合。
入力例専用の化学種順序はファイル冒頭を確認する。

Windows PowerShell：

```powershell
cmake --build build/reactingflow-fortran-release --target rf_flow2d
.\build\reactingflow-fortran-release\rf_flow2d.exe mechanism.rf SolverLibrary/ReactingFlow/examples/flow2d_split_h2air.in build/flow2d_split.csv
```

Linux：

```bash
cmake --build build/reactingflow-fortran-release --target rf_flow2d
./build/reactingflow-fortran-release/rf_flow2d mechanism.rf SolverLibrary/ReactingFlow/examples/flow2d_split_h2air.in build/flow2d_split.csv
```

実行は今回行っていない。split/boxの選択検査コードを既存の `rf_flow2d_unit` に追加し、
ビルドのみ確認。任意分布の読書きや衝撃波管の動作検証も未実施。
