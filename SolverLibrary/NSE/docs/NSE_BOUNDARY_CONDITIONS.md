# NSE境界条件仕様

更新日: 2026-09-01  
状態: CPU/MPI/OpenMP、単一GPU CUDA、MPI＋CUDA版で実装済み

## 1. 目的

NSEソルバーで、x、y、zの各方向および各方向の下端・上端へ異なる物理境界条件を
指定できるようにする。本仕様の最初の実装対象は次の組合せである。

- `periodic`: 周期境界
- `non_reflecting`: 圧縮性流れの特性波に基づく無反射遠方境界
- `reflective`: 法線運動量を反転する自由滑り・断熱の鏡像境界
- `dirichlet`: 指定した密度・速度・圧力を全ghost層へ固定する流体リザーバ境界

代表的な用途は、x方向の両端だけを無反射境界とし、横方向のy、zを周期境界とする
孤立乱流・移流乱流計算である。

2026-08-28時点で、CPU/MPI/OpenMP、単一GPU CUDA、MPI＋CUDAの実行コード、
`case_input.py`、標準`case.yaml`テンプレートは本仕様に対応している。

## 2. 正規のcase.yaml形式

新しい正規形式は、6物理面を`boundary.faces`で明示する。

```yaml
boundary:
  faces:
    x_min:
      type: non_reflecting
      reference_state: far_field
    x_max:
      type: non_reflecting
      reference_state: far_field
    y_min:
      type: periodic
    y_max:
      type: periodic
    z_min:
      type: periodic
    z_max:
      type: periodic

  reference_states:
    far_field:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143

  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto
```

`x_min`、`x_max`などは物理座標の下端面・上端面を表す。MPI rankの局所端ではない。
新形式では6面をすべて指定し、省略時の暗黙境界を設けない。

### 2.1 面の設定

| キー | 型 | 必須条件 | 意味 |
|---|---|---|---|
| `type` | 文字列 | 全面で必須 | `periodic`、`non_reflecting`、`reflective`、`dirichlet` |
| `reference_state` | 文字列 | `non_reflecting`と`dirichlet`で必須 | `boundary.reference_states`の名前 |

周期境界は方向ごとの対で指定する。例えば`y_min: periodic`なら`y_max`も
`periodic`でなければならない。一方の面だけを周期境界にする入力は生成時および
実行開始時に拒否する。

`reference_state`は`non_reflecting`面と`dirichlet`面だけで使用する。`periodic`または`reflective`面へ
指定した場合は入力エラーとする。鏡像面は方向内で対にする必要はなく、反対面へ
無反射境界などを設定できる。

### 2.2 基準状態

`boundary.reference_states`には、無反射境界の外部一様状態またはDirichlet固定状態を名前付きで定義する。

| キー | 制約 | 意味 |
|---|---|---|
| `density` | 正値 | 無次元密度 `rho_ref` |
| `velocity` | 有限な3成分 | 無次元速度 `[u_ref,v_ref,w_ref]` |
| `pressure` | 正値 | 無次元圧力 `p_ref` |

値はNSEソルバー本体と同じ音響スケーリングによる無次元量とする。全エネルギーは
`gamma`と基準プリミティブ変数から計算し、入力では重複指定しない。

複数の名前付き状態を定義できるため、例えば流入側と流出側で異なる遠方状態を
割り当てられる。どの面からも参照されない状態は許可するが、未定義名の参照は
エラーとする。

### 2.3 無反射境界パラメータ

| キー | 既定値・制約 | 意味 |
|---|---|---|
| `formulation` | `characteristic_relaxation`固定 | 特性波・緩和型遠方境界 |
| `relaxation_strength` | `0.1`、0以上 | 基準状態へ戻す無次元強度 |
| `length_scale` | `auto`または正値 | 緩和率に使う代表長さ |

`length_scale: auto`は対象面の法線方向領域長を使う。基準音速を`c_ref`、代表長さを
`L_ref`、緩和強度を`sigma`として、緩和率の基準を

```text
kappa = sigma*c_ref/L_ref
```

とする。`relaxation_strength: 0`は基準状態への緩和を無効にするが、特性波の
流入・流出判定は有効なままとする。

### 2.4 鏡像境界の設定例

x両端を鏡像面、y-z両端を周期面にする例を次に示す。

```yaml
boundary:
  faces:
    x_min: {type: reflective}
    x_max: {type: reflective}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
  reference_states: {}
```

`reflective`面へ`reference_state`は指定しない。全6面を鏡像にする場合も、各面の
`type`を`reflective`へ変更するだけでよい。

## 3. 後方互換性

既存の

```yaml
numerics:
  boundary_condition: periodic
```

は互換入力として残し、`boundary`が存在しない場合に限り6面すべての
`type: periodic`へ展開する。

- `boundary`と`numerics.boundary_condition`の同時指定は、優先順位を設けずエラーにする。
- 旧キーで許可する値は`periodic`だけとする。
- 新規テンプレートは`boundary`形式を使用する。
- Fortran namelistの旧単一変数`boundary_condition`は移行期間だけ受理し、最終的には
  6面の正規設定へ内部展開する。

## 4. 無反射境界の数値仕様

### 4.1 特性波の判定

各面で領域外向き単位法線`n`を使い、境界内側セルのプリミティブ状態から

```text
u_n = velocity dot n
c   = sqrt(gamma*p/rho)
lambda = [u_n-c, u_n, u_n, u_n, u_n+c]
```

を評価する。`lambda < 0`を領域へ入る特性、`lambda >= 0`を領域から出る特性とする。
出る特性は境界隣接内点の値を保持し、入る特性だけを指定した基準状態へghost層の
距離に応じて緩和する。その後、特性量をプリミティブ変数へ戻し、保存変数
`[rho,rho*u,rho*v,rho*w,rho*E]`を再構築する。

これにより、亜音速流出、亜音速流入、超音速流出、超音速流入を同じ固有値判定で
切り分ける。単純なゼロ勾配コピーを`non_reflecting`として実装してはならない。

### 4.2 ghostセル

- KEEP6、WENO5-Z/Roe、ハイブリッドおよびcentral6粘性項に必要な3層を埋める。
- 無反射面では境界法線方向へ特性関係を外挿し、各ghost層の保存変数を構築する。
- 粘性流束は同じ特性ghost状態を使用する。壁面や規定温度条件を別に重ねない。
- `small_rho`と`small_p`は固有量評価の判定閾値に使用し、不正な保存状態を黙って
  正値へ置換しない。CPU版は面名を出して停止し、CUDA版は非有限状態として後段の
  backend検査へ伝播する。

### 4.3 面・辺・角の完成順序

一回の`apply_nse_boundary`呼出しが、対流・粘性演算に必要な面、辺、角ghostを
すべて完成させる。処理順を次のように固定する。

1. x、y、zの順に物理境界演算子を適用する。
2. 周期指定された物理端は、対応する方向の演算子内で対向面を接続する。
3. y、z方向のMPI内部隣接rank間haloを交換する。
4. `central6`では交換をもう一度行い、完成した物理面値を横方向haloへ伝播する。

後段の方向は、前段で完成した横方向ghostを入力として使用する。OpenMP、MPI分割数、
CPU/CUDAによってこの順序を変えてはならない。一様な基準状態は、周期面・無反射面が
交わる辺と角を含め機械精度で保存する。

## 5. 鏡像境界の数値仕様

`reflective`はセル中心格子上の鏡映でghostを構築する。法線方向の物理セル数を`N`、
境界から数えたghost層を`g = 1, 2, 3`とすると、下端と上端の保存状態はそれぞれ

```text
Q(1-g) = R_n Q(g)
Q(N+g) = R_n Q(N+1-g)
```

である。`Q = [rho,rho*u,rho*v,rho*w,rho*E]`に対し、`R_n`は面法線方向の運動量だけ
符号を反転し、密度、二つの接線運動量、全エネルギーを変えない。したがって境界面上の
法線速度は奇対称でゼロ、圧力、温度、接線速度は偶対称となる。

この条件は、直交格子に沿う不透過の自由滑り・断熱壁または対称面である。粘着壁
（no-slip）、移動壁、規定温度壁ではない。`central6`粘性項との組合せでは接線速度と
温度の法線勾配がゼロとなり、壁面せん断と熱流束を与えない。

辺と角ではx→y→zの順に鏡映を合成する。二つの鏡像面が交わる辺では二つの法線運動量、
三つが交わる角では三つの運動量成分を反転する。各方向は対応する内点深さから写すため、
3層すべてが同一セルの単純コピーになることはない。

## 6. 並列実行との関係

### 6.1 CPU、MPI、OpenMP

- x方向は分割されないため、全rankがxの物理境界を処理する。
- y、z方向は、領域内部のrank境界で従来どおりhalo交換する。
- `js=1`、`je=ny`、`ks=1`、`ke=nz`を持つrankだけが対応する物理境界を処理する。
- 周期方向だけが領域反対端rankと通信する。無反射／鏡像方向は反対端rankと通信しない。

### 6.2 CUDAとMPI＋CUDA

CPU版と同じ面設定、基準状態、特性関係、緩和係数を使用する。GPU常駐時間発展では
周期・無反射・鏡像境界を一つのCUDAカーネルで処理し、境界処理のために全場をhostへ
戻さない。面、辺、角ではx、y、zの順に状態変換を適用する。

MPI＋CUDAでは内部halo交換と物理端の境界カーネルを分離する。周期方向の物理端rankは
領域反対端rankと通信し、無反射／鏡像方向の物理端は`MPI_PROC_NULL`として通信せず、GPU上で
生成したghostを保持する。領域内部のrank同士は境界種別にかかわらず通常のhalo交換を行う。

段階導入の状況は次のとおりである。

1. CPU MPI/OpenMP版の混合周期・無反射・鏡像境界: 実装済み。
2. 単一GPU CUDA版のGPU常駐境界カーネル: 実装済み。
3. MPI＋CUDA版の物理端rank処理と内部halo交換: 実装済み。

未知の境界種別、方向内で対にならない周期面、無効な基準状態は入力生成時または
実行開始時に明示的に拒否し、周期境界へ暗黙に置き換えない。

## 7. 他機能との組合せ

| 機能 | 混合境界との関係 |
|---|---|
| Petersen–Livescu FFT Forcing | 全6面周期が必須。無反射面または鏡像面が1つでもあればエラー |
| HITスペクトル初期化 | 周期FFTで初期場を作る。非周期時間発展との併用時は警告を表示 |
| imported turbulence `embed` | x無反射、y-z周期を主用途とする。鏡像面との併用も入力上は可能 |
| imported turbulence `tile` | 境界設定とは独立。ただし物理的整合は利用者が確認 |
| finite-driver shock tube | `x_min=reflective`、`x_max=non_reflecting`。高圧状態をDirichlet供給しない |
| KEEP2／KEEP6 | 使用可能。必要ghostを境界処理が供給 |
| WENO5-Z/Roe／hybrid | 使用可能。必要ghostを境界処理が供給 |
| central6粘性項 | 使用可能。無反射面は遠方境界、鏡像面は自由滑り・断熱として閉じる |

## 8. 検証要件

CPU版とCUDA版では、少なくとも次を回帰確認する。音響反射率などの定量検証は継続課題とする。

- 全方向周期の既存面・辺・角テストが後方互換で成功する。
- x無反射・y-z周期で一様な基準状態を機械精度で保存する。
- x方向一次元微小音響パルスを流出させ、反射圧力を測定する。
- エントロピー波と渦度波が流出時に不自然な圧力波を生成しない。
- 亜音速・超音速の流入／流出で、特性の入出数が固有値符号と一致する。
- mixed boundaryの全辺・角ghostが有限かつ正の密度・圧力を持つ。
- MPI分割数を変更して同じ物理解と境界反射指標を得る。
- OpenMPスレッド数を変更して決定性を維持する。
- CPUとCUDAの境界状態を全6面、辺、角で許容誤差内に一致させる。
- 無反射面とPetersen–Livescu Forcingの併用が明示的に拒否される。
- 全6面鏡像で各ghost層が対応する内点深さを参照し、面法線運動量だけを反転する。
- 鏡像面が交わる辺・角で、交差する各面の法線運動量がそれぞれ反転する。
- CPU/MPI、単一GPU、MPI＋CUDAで鏡像境界を含む時間発展が正常終了する。

音響パルス試験では、境界到達前の入射圧力変動と、通過後に観測領域へ戻る圧力変動から

```text
R = norm(p_reflected) / norm(p_incident)
```

を記録する。受入れ閾値は格子数、時間刻み、対流方式とともにテストへ明記し、単なる
ゼロ勾配境界より反射が小さいことを最低条件とする。

## 9. 実装対象ファイル

実装箇所は次のとおりである。

- `src/common/mod_model_config.f90`: 6面種別、基準状態、無反射パラメータ
- `src/io/mod_input_reader.f90`: namelist入力と旧単一指定の展開
- `src/boundary/mod_boundary_runtime.f90`: 面別dispatcher、周期、無反射特性、鏡像演算子
- `src/gpu/nse_cuda_bridge.cu`: GPU常駐の面別周期／無反射／鏡像境界カーネル
- `src/gpu/mod_nse_gpu_cuda.f90`: 面設定と基準状態をCUDA contextへ渡すFortran API
- `src/gpu/mod_nse_gpu_mpi.f90`: MPI＋CUDAの内部隣接通信と物理端の切分け
- `ScriptLibrary/RunEnvironment/case_input.py`: YAML検証とnamelist生成
- `ScriptLibrary/RunEnvironment/case_templates/nse.yaml`: 新しい正規入力例
- `solver_manifest.yaml`、CMake、module catalog: runtime境界backend
- `tests/test_boundary_non_reflecting.f90`: MPI混合境界の一様場、非周期性、正値性テスト
- `tests/test_cuda_boundary_non_reflecting.f90`: CPU/CUDAの全面・辺・角境界値比較
- `tests/test_boundary_reflective.f90`: MPIの鏡像面・辺・角および内部haloテスト
- `tests/test_cuda_boundary_reflective.f90`: 鏡像ghostの定義確認とCPU/CUDA比較
- `tests/test_boundary_dirichlet.f90`: MPIのDirichlet面・辺・角と全ghost層の検査
- `tests/test_cuda_boundary_dirichlet.f90`: Dirichlet ghostのCPU/CUDA一致
- `ScriptLibrary/RunEnvironment/tests/test_case_input.py`: CPU/CUDA入力生成と不正構成の拒否テスト

## 10. 衝撃波流入用Dirichlet境界

平面衝撃波–乱流干渉では、衝撃波の背後状態を進行方向上流側のx面へ固定する。
`positive_x`では`x_min`、`negative_x`では`x_max`を`dirichlet`にする。YAMLでは
`source: planar_shock.downstream`を使い、初期条件と境界条件へ同じ値を重複入力しない。
詳細と完全な設定例は
[`NSE_SHOCK_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TURBULENCE_INTERACTION.md)を参照する。

有限高圧室の衝撃波管ケースでは`x_min`を`reflective`閉端、`x_max`を低圧状態参照の
`non_reflecting`とする。完全な設定は
[`NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)を参照する。
