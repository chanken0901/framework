# 単成分NSE：Landau–Lifshitz揺らぎ拡張（初版）

## 対応範囲

単成分理想気体・一定比熱比・一定粘性／熱伝導・体積粘性ゼロのLLNS拡張。
反応や化学種は不要。無効時は従来NSEの演算経路を維持する。

- 対応：直交一様格子、全周期境界、CPU MPI/OpenMP、CUDA／MPI+CUDA、固定時間刻み。
- 未対応：非周期境界、一般座標、多成分・反応流への揺らぎ追加。
- 有効時は粘性・熱伝導を専用の前進勾配／随伴発散演算子に置換する。
  **通常の六次精度粘性演算子へ単に乱数を加える仕様ではない**。
  交差微分を含む輸送全体は空間一次精度（軸方向単独の拡散は二次精度）。
- 対流はKEEP6。WENO/HYBRIDの数値散逸に対応する熱雑音は未実装なので拒否する。
- 決定論的SSPRK3の後に、ステップ開始状態で評価したItô増分を一回加える。
  確率系としては弱一次精度。有限dtで平衡分散が厳密になる方式ではない。
- 検証は一様・静止平衡で線形化した散逸／雑音の整合性が中心。
  非線形・強非平衡の不変分布や全波数の動的構造因子は未検証。
  初版は研究・検証用途とし、実ケースでは格子／dt収束と統計の確認が必要。

## Landau–Lifshitzの項

密度に直接雑音を加えず、運動量と全エネルギーへ追加する。

```text
momentum: div(Sigma)
energy:   div(Sigma . u + Xi)

<Sigma_ij(x,t) Sigma_kl(x',t')>
 = 2 k_B mu T (delta_ik delta_jl + delta_il delta_jk
              - (2/3) delta_ij delta_kl) delta(x-x') delta(t-t')
<Xi_i(x,t) Xi_j(x',t')>
 = 2 k_B kappa T^2 delta_ij delta(x-x') delta(t-t')
<Sigma_ij Xi_k> = 0
```

いずれも零平均のガウス雑音。応力は対称・トレースゼロ、熱流束とは独立である。
これは乱流forcingとは別の熱揺らぎである。
参考：[Donevほか、LLNSの離散化と揺動散逸関係](https://arxiv.org/html/0906.2425v2)。
同論文の高次RK方式をそのまま実装したものではない。

## 無次元化

本ソルバーで`theta = p/rho = R_specific*T_physical/U_ref^2`とすると、

```text
beta = k_B / (rho_ref * R_specific * L_ref^3)
mu_star = 1 / Re
kappa_star = gamma / ((gamma-1) * Re * Pr)
```

`rho_ref`、`L_ref`は有次元の基準密度・長さ、`R_specific`は気体の比気体定数。
**有次元のk_B = 1.380649e-23をそのままYAMLへ入れない**。
理想気体のセル内分子数の目安は`N_cell = rho_star * V_star / beta`。
十分多数の分子を含むセルと局所平衡を前提とし、極端な格子細分化では連続体近似が破綻する。

離散化ではデルタ関数を`1/(V_star*dt_star)`で置き換える。
応力の基本振幅は`sqrt(2*beta*mu_star*theta/(V_star*dt_star))`、
熱流束は`sqrt(2*beta*kappa_star*theta^2/(V_star*dt_star))`。
RHSの標準偏差は`dt^(-1/2)`、保存量増分は`dt^(1/2)`となる。

## 離散化・並列再現性

勾配を`D+`、発散を`D- = -(D+)^T`とし、同じ位置に対称応力を構成する。
静止平衡で線形化した運動量の散逸演算子Lと任意の速度プローブwに対し、
`Var(w^T noise_rhs) = -2*beta*theta/(V*dt) * w^T L w`となる。
熱伝導にも同じ随伴対を使う。既存CENTRAL6と異なるので有効時だけ決定論的輸送も置換する。
エネルギー流束は各方向の隣接速度平均と同じ応力を使用する。
周期領域では流束差が相殺され、全質量・全運動量・全エネルギーを丸め誤差範囲で保存する。

乱数はPhilox4x32-10＋Box–Muller。周期折返しした全体セル座標、ステップ、seed、
成分番号で決める。rank／thread番号には依存しない。
[Random123既知解](https://github.com/DEShawResearch/random123/blob/main/tests/kat_vectors)と照合する。
数学ライブラリが異なる機種間の浮動小数点ビット一致は保証しない。
厳密な再開には元の全体ステップも必要。初期場としてSLFを読み直す処理は確率過程の再開ではない。

固定dtの安定性を各ステップで確認する。負密度／負圧等は元状態へ戻して異常終了する。
受理・棄却による統計の偏りを避けるため、乱数の再抽選や雑音のクリップは行わない。
刻みを変更する場合は再実行する。

### CUDAとMPI＋CUDA

保存変数、応力・熱流束、ステップ開始時に評価する確率増分はGPU上で生成・保持する。
各RKステージでこれらの全配列をCPUへ戻さない。拡張有効時のみ、ghostを含む局所セル数Nに
対して流束12N要素と確率増分5N要素（float64、追加136N byte）を確保する。
既存の保存変数・RK作業配列等のメモリは別途必要。

MPI版は既存のy-z分割と保存変数のhalo交換を使用する。乱数は全体座標から再生成し、
応力や熱雑音をMPI通信しない。halo交換はCPUバッファ経由、または
[CUDA-aware MPI直接通信](NSE_CUDA_AWARE_MPI.md)を選択できる。実際のGPUDirect利用は実行環境に依存する。
安定性・異常判定のスカラー通信と出力用転送も残る。固定dtがいずれかのrankで
安定性制約を超えた場合、または更新後に負密度・負圧等が生じた場合は異常終了する。
通常運用は1 MPI rank／GPUとする。

## 設定

新規生成には`ScriptLibrary/RunEnvironment/environment.nse_fluctuating.yaml`を使う。
デモ初期場はTaylor–Green渦で、静止平衡ではない。次の参照と設定ファイルが生成される。
既存ケースでも同様に追加できる。

case.yaml（既存のschema_versionは2へ変更）：

```yaml
schema_version: 2
extensions:
  fluctuating_hydrodynamics: config/fluctuating_hydrodynamics.yaml
```

config/fluctuating_hydrodynamics.yaml：

```yaml
schema_version: 1
extension: fluctuating_hydrodynamics
config:
  enabled: true
  model: landau_lifshitz
  boltzmann_number: 1.0e-8  # デモ値。対象気体と基準量から再設定する
  seed: 13579
```

case.yamlにはKEEP6、CENTRAL6、SSPRK3、正のReとPr、全周期境界、
`time.use_fixed_dt: true`と正のdtを指定する。
無効化または拡張参照の省略で従来NSEへ戻る。
`enabled: true, boltzmann_number: 0`は専用輸送のまま雑音をゼロにする比較用であり、
拡張無効時のCENTRAL6計算とは異なる。

## WindowsとLinuxの実行手順

Windows / PowerShell（FrameWork直下）：

```powershell
python .\ScriptLibrary\RunEnvironment\prepare_environment.py .\ScriptLibrary\RunEnvironment\environment.nse_fluctuating.yaml --dry-run
python .\ScriptLibrary\RunEnvironment\prepare_environment.py .\ScriptLibrary\RunEnvironment\environment.nse_fluctuating.yaml
# 生成ログにある実行環境フォルダへ移動し、case.yamlとconfigの物理条件を確認後：
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

Linux / bash：設計書の`select.source`を`framework_relative`、`select.target`を
`linux_gnu_mpi`へ変更し、出力先もLinux用に変更する。設計書をRunEnvironment内に
置いたまま、例えば`select.destination: local_generated`と
`destination.root: generated/nse_fluctuating`を指定する。

```bash
python3 ./ScriptLibrary/RunEnvironment/prepare_environment.py ./ScriptLibrary/RunEnvironment/environment.nse_fluctuating.yaml --dry-run
python3 ./ScriptLibrary/RunEnvironment/prepare_environment.py ./ScriptLibrary/RunEnvironment/environment.nse_fluctuating.yaml
# 生成された実行環境フォルダへ移動し、条件確認後：
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --run
```

既存のResearchRunsには自動反映されない。更新版から実行環境を生成／更新し、
入力再生成とビルドが必要。実装変更だけで既存ケースの条件は変更しない。

### CUDA用の設計書

上記の生成コマンドの設計書名を`environment.nse_fluctuating.cuda.yaml`へ変更する。
既定は単一GPUであり、MPI＋CUDAには設計書の`parallel.use_mpi`を`true`に変更する。
`use_openmp: false`、`use_cuda: true`を維持する。profileはそれぞれ`cuda_single`、
`cuda_mpi`が自動選択される。LLNS自体はFFTライブラリを使用しない。
Linuxでは上記と同じtarget／出力先の変更が必要で、CUDA Toolkitと対応ホストコンパイラ、
MPI有効時はMPIライブラリも用意する。生成後の`run_case.py`の手順はCPU版と共通。
GPUの選択・MPI起動設定は[NSEビルド手順](NSE_BUILD_AND_RUN.md)も参照する。

## テスト

`nse_fluctuating`は乱数既知解、平均・共分散・応力対称性・トレース、周期保存、
分割再現性、dtスケーリング、無効時不変性、線形化した運動量／熱の離散FDTを検証。
`nse_fluctuating_smoke`は入力読込みから3ステップの正常終了を検証する。
MPIビルドでは`nse_fluctuating_mpi_smoke`を2rank×2threadで実施する。
`nse_cuda_fluctuating_compare`と`nse_cuda_fh_zero_compare`は非一様場で4ステップの
CPU／CUDA一致を非零雑音／零雑音で検証する。`nse_cuda_fh_mpi_2`／`_4`は
2／4rank分割と単一領域CUDAの一致・周期保存を検証する。
MPIテストの登録には`NSE_ENABLE_MULTI_GPU_TESTS=ON`を指定する。
開発機では1GPUを共有して2／4rankを検証済み。物理的に異なる複数GPU間の
通信・性能は利用先の実機で別途検証すること。

Windows/Linux共通（ビルドディレクトリは実環境に合わせる）：

```text
ctest --test-dir build/nse-positivity-cpu -R fluctuating --output-on-failure
ctest --test-dir build/nse-positivity-mpi -R fluctuating --output-on-failure
ctest --test-dir build/nse-positivity-cuda -R "fluctuating|fh_zero" --output-on-failure
ctest --test-dir build/nse-positivity-mpi-cuda -R nse_cuda_fh_mpi --output-on-failure
```
