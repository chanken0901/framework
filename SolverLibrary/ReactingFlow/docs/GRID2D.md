# 2次元構造格子・ノズル格子生成

追記：[2D流体・反応の実行入口](FLOW2D.md)へ接続した。格子専用 `rf_mesh2d` は引き続き時間発展を行わない。

2026-10-02。多次元有限体積基盤へ渡す幾何情報をFortranで生成する。
今回の確認はRelease/Debugビルドのみ。テスト実行と反応流の実用検証は保留中。

## 実装範囲

- 任意の構造格子節点 `nodes(2,0:nx,0:ny)` から面・セル情報を生成。
- 上下の壁面を指定した、2次元平面の収縮・拡大ノズル格子。
- 非一様なx方向節点間隔。上下壁の間はny等分。
- セルの面積・面積重心、面中心、向き付き面積ベクトル、面接続、境界区分。
- 可視化用のASCII VTK Structured Grid出力。

軸対称ではない。2次元の「体積」は単位奥行きの面積、面積ベクトルの大きさは辺長。
断面積の変化を流体方程式の準1Dソース項で模擬するものでもない。
壁は隣接節点間の直線として扱い、曲線補間や自動的な格子細分化は行わない。
既存の `rf_flow1d`、NSE、GPEの入力・実装は変更しない。

## 格子入力

例：`examples/nozzle_profile.rf`。座標単位はm。

```text
RF_NOZZLE_PROFILE_V1
3
0.00 -0.005 0.005
0.02 -0.002 0.002
0.04 -0.005 0.005
```

2行目はx方向の節点数。以後は `x y_lower y_upper` を節点数だけ記述する。
xは厳密な昇順、各断面の高さは正とする。nxは節点数−1。
各行は3個の数値のみとし、コメントや列追加は行わない。
同一高さの壁を指定すれば矩形格子、上下非対称の壁も指定可能。
入力例は格子生成の説明用で、反応帯を解像する計算格子ではない。

## ビルド・格子出力

FrameWorkルートで、既存のReactingFlowビルドディレクトリを使用する。
初回構成・DVODEの指定は [Fortranビルド手順](REACTORS.md) を参照。

Windows PowerShell（既存のReleaseビルドを使う例）：

```powershell
cmake --build build/reactingflow-fortran-release --target rf_mesh2d
.\build\reactingflow-fortran-release\rf_mesh2d.exe SolverLibrary/ReactingFlow/examples/nozzle_profile.rf 16 build/nozzle_grid.vtk
```

Linux（ビルドディレクトリ名は実際の構成に合わせる）：

```bash
cmake --build build/reactingflow-fortran-release --target rf_mesh2d
./build/reactingflow-fortran-release/rf_mesh2d SolverLibrary/ReactingFlow/examples/nozzle_profile.rf 16 build/nozzle_grid.vtk
```

第2引数は横断方向セル数ny。出力はParaView等で読める節点と `cell_area`。
出力親ディレクトリは事前に存在する必要がある。既存ファイルは上書きせずエラーにするため、
再出力時は別の出力ファイル名を選ぶ。今回はこの実行例を実行していない。
これは格子生成の実行入口であり、流体・化学の時間発展は行わない。

## 内部契約と次の接続

`mod_rf_grid2d` の `rf_grid2d` は、既存 `rf_face_mesh` を `grid%mesh` として所有する。
セル番号は `i + nx*(j-1)`。内部面は一度だけ生成する。
`grid%boundary` は0=内部、1=i最小、2=i最大、3=j最小、4=j最大。
公開定数 `rf_imin/rf_imax/rf_jmin/rf_jmax` を使う。
これは幾何学上の区分であり、流入／流出／鏡像等の物理境界条件はまだ割り当てない。
追加した `mod_rf_flow2d` がこの区分からghostを生成し、`finite_volume_rhs` に渡す。

セルは反時計回りの厳密な凸四角形であることを要求する。
逆転・凹形・退化、非有限座標、ゼロ高さ、面積ベクトル非閉鎖は拒否する。
極端に細長いセルも丸め誤差近傍の凸性検査で拒否される場合がある。
任意節点入力での離れたセル同士の重なりまで検出する検査ではない。
移動格子・3D押出し・周期面の接続・境界層格子自動生成は未実装。

`rf_grid2d_unit` に矩形・台形・せん断格子の幾何検査と異常入力拒否を追加した。
コンパイルのみ確認し、実行結果を合格とは記録していない。
