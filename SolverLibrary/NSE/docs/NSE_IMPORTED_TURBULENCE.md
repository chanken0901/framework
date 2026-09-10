# 保存済み乱流場の読込みとx方向配置

> 本書ではWindows（PowerShell）とLinux（bash）のコマンドを併記します。共通の読み替えは[`../../../docs/WINDOWS_LINUX_COMMANDS.md`](../../../docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

## 概要

`flow.type: imported_turbulence`は、NSEのrank別SLF出力を一つの可搬SLFへ
変換して読み込み、より長いx方向計算領域へ初期乱流場を配置する。

次の三方式を選択できる。

- `embed`: 乱流場をx方向の一領域へ配置し、それ以外を指定した背景場にする。
- `tile`: 乱流場をx方向へ周期的に繰り返し、領域全体を埋める。
- `periodic_embed`: 指定長さの区間だけ乱流場を周期的に繰り返し、外側は背景場にする。

### 指定長さの局所乱流：periodic_embed

標準テンプレートには`x_length: null`を明示している。`periodic_embed`を選ぶ場合は
`null`を希望する正の物理長さへ変更する。古いケースに項目がなければ追加する。
長さを自動決定することはない。`embed`／`tile`では省略または`null`のままにする。
数値を指定できるのは`periodic_embed`のみである。
以下は配置部分の例であり、背景状態などの既存設定はそのまま併記する。

```yaml
flow:
  type: imported_turbulence
  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: periodic_embed
    x_start: 4.0
    x_length: 10.0
    blend_cells: 4
```

この例では`4.0 <= x < 14.0`だけに乱流を配置する。先頭から元データを読み、
末尾に達したら先頭へ戻る。`x_length`は元データの長さの整数倍でなくてもよく、
元データより短くてもよい。ただし開始位置はセル境界、長さは読込み先の`dx`の
正の整数倍とし、区間全体を計算領域内に収める。格子間隔・y/z領域などの互換条件は従来と同じ。

`blend_cells`は指定区間全体の両端だけに適用し、反復の継ぎ目には適用しない。
区間のセル数の半分以下を指定する。`0`では背景への接続が不連続になる場合がある。
同じ周期データの反復であり、独立した乱流を新たに生成する機能ではない。
元データはx方向に周期的な場を用意する（SLFから周期性を自動判定しない）。

`shock_turbulence_interaction`と`shock_tube_turbulence_interaction`でも同じ指定が可能。
既存の衝撃波／高圧室設定を維持し、乱流区間が衝撃波背後や高圧室に重ならないようにする。
`embed`と`tile`の動作は変更しない。`x_length`は`periodic_embed`専用である。
変更後は`input.dat`を再生成し、本機能を含むソルバーを再ビルドする。
既に生成済みの実行環境にはライブラリ更新が自動反映されないため、実行環境の更新も必要。

CPU MPI/OpenMP版、単一GPU CUDA版、MPI＋CUDA版で利用できる。MPI版では各rankが
自分のy-z局所領域だけをファイルから読み込む。初期データ生成時と本計算時の
MPI並列数は一致しなくてもよい。

## 1. 対応する計算の組合せ

元計算のバックエンドやMPI並列数は、読込み先へ引き継がない。元出力を一度
可搬SLFへ変換し、読込み先が現在のCPU MPI分割、単一GPU配置、またはMPI＋CUDA分割に合わせて
初期場を構築する。

| 元計算 | 読込み先 | 対応 | 読込み先で指定する並列設定 |
|---|---|---|---|
| CPU MPI | CPU MPI | 対応 | 読込み先で必要なMPIプロセス数を指定 |
| CPU MPI | CUDA | 対応 | 単一GPU、MPI/OpenMP無効 |
| CUDA | CPU MPI | 対応 | 読込み先で必要なMPIプロセス数を指定 |
| CUDA | CUDA | 対応 | 単一GPU、MPI/OpenMP無効 |
| CPU MPI／CUDA／MPI＋CUDA | MPI＋CUDA | 対応 | 読込み先のGPU数と同じMPIプロセス数 |

例えば、24 MPIプロセスのCPU計算から1 GPUへ移す場合や、1 GPUのCUDA計算から
16 MPIプロセスのCPU計算へ移す場合にも、同じ変換コマンドを使用する。

## 2. 最短実行手順

以下は、読込み先の実行環境がすでに生成されている場合のWindows／Linux手順である。
`nse_case0001`などの番号は実際の元ケースと読込み先ケースへ置き換える。

読込み先環境がまだない場合は、先に設計書をコピーしてCPU用またはCUDA用の
`parallel`設定へ編集し、実行環境を生成する。設定値は「2.5」「2.6」に示す。

```powershell
$framework = "$env:USERPROFILE\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$design = "$env:USERPROFILE\ResearchRuns\Designs\nse_imported_target.yaml"

New-Item -ItemType Directory -Force (Split-Path $design) | Out-Null
Copy-Item "$tool\environment.nse.yaml" $design
notepad $design

python "$tool\prepare_environment.py" $design --dry-run
python "$tool\prepare_environment.py" $design
```

Linux（bash）:

```bash
framework="$HOME/Research/FrameWork"
tool="$framework/ScriptLibrary/RunEnvironment"
design="$HOME/ResearchRuns/Designs/nse_imported_target.yaml"

mkdir -p "$(dirname "$design")"
cp "$tool/environment.nse.yaml" "$design"
${EDITOR:-vi} "$design"

# source/destination/targetもLinux用へ変更してから実行する
python3 "$tool/prepare_environment.py" "$design" --dry-run
python3 "$tool/prepare_environment.py" "$design"
```

生成ログに表示された`nse_caseNNNN`を、以下の`$targetRoot`へ指定する。

### 2.1 パスを設定する

```powershell
$sourceOutput = "$env:USERPROFILE\ResearchRuns\nse_case0001\cases\case0001\output"
$targetRoot = "$env:USERPROFILE\ResearchRuns\nse_case0002"
$targetCase = "case0002"
$portableSlf = Join-Path $targetRoot "cases\$targetCase\initial_data\turbulence.slf"

Set-Location $targetRoot
New-Item -ItemType Directory -Force (Split-Path $portableSlf) | Out-Null
```

Linux（bash）:

```bash
source_output="$HOME/ResearchRuns/nse_case0001/cases/case0001/output"
target_root="$HOME/ResearchRuns/nse_case0002"
target_case="case0002"
portable_slf="$target_root/cases/$target_case/initial_data/turbulence.slf"

cd "$target_root"
mkdir -p "$(dirname "$portable_slf")"
```

元計算がCPU MPIでもCUDAでも、`$sourceOutput`には元ケースの`output`ディレクトリを
指定する。

### 2.2 元出力を確認する

```powershell
if (-not (Test-Path "$sourceOutput\meta.json")) {
  throw "meta.jsonがありません: $sourceOutput"
}

$slfFiles = Get-ChildItem -LiteralPath $sourceOutput -Filter "*.slf"
if ($slfFiles.Count -eq 0) {
  throw "SLFファイルがありません: $sourceOutput"
}

$slfFiles | Select-Object Name, Length, LastWriteTime
```

Linux（bash）:

```bash
test -f "$source_output/meta.json" || {
  echo "meta.jsonがありません: $source_output" >&2
  exit 1
}

find "$source_output" -maxdepth 1 -type f -name '*.slf' -print -quit | grep -q . || {
  echo "SLFファイルがありません: $source_output" >&2
  exit 1
}

find "$source_output" -maxdepth 1 -type f -name '*.slf' -printf '%f %s bytes\n'
```

元ケースでは`output.write_meta: true`を使用する。CPU MPI／MPI＋CUDA出力では指定stepの
全rankファイルが必要であり、単一GPU出力では単一rankファイルと`meta.json`が必要である。

### 2.3 可搬SLFへ変換する

```powershell
python .\SolverLibrary\NSE\tools\nse_prepare_imported_turbulence.py `
  $sourceOutput `
  --meta "$sourceOutput\meta.json" `
  --step latest `
  --output $portableSlf

if ($LASTEXITCODE -ne 0) {
  throw "乱流SLFの変換に失敗しました"
}
```

Linux（bash）:

```bash
python3 ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py \
  "$source_output" \
  --meta "$source_output/meta.json" \
  --step latest \
  --output "$portable_slf" || {
    echo "乱流SLFの変換に失敗しました" >&2
    exit 1
  }
```

成功時には次の形式で表示される。

```text
[OK] Prepared imported turbulence: ...\turbulence.slf (step=<step>, grid=<nx>x<ny>x<nz>)
```

特定stepを使う場合は`--step latest`を`--step 1000`のように変更する。
CPU/CUDAの判別とCPU側の元MPI並列数は`meta.json`から自動取得されるため、通常は
`--layout`を指定しない。

### 2.4 読込み先のcase.yamlを編集する

Windowsでは`$targetRoot\cases\$targetCase\case.yaml`、Linuxでは
`$target_root/cases/$target_case/case.yaml`の`flow`を次のように設定する。

```yaml
flow:
  type: imported_turbulence
  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: embed
    x_start: 2.0
    blend_cells: 8
    velocity_offset: [0.5, 0.0, 0.0]
    background:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143
```

同じ`case.yaml`で`grid`、`physics.nse.gamma`、無次元化を元データと整合させる。
詳しい格子条件は「6. 格子と物理量の互換条件」を参照する。

### 2.5 読込み先がCPU MPIの場合

CPU用の実行環境は、生成前のenvironment設計書で次のように選ぶ。

```yaml
parallel:
  use_mpi: true
  use_openmp: true
  use_cuda: false
```

生成済みCPU環境では、`case.yaml`の実行時並列数を設定する。

```yaml
solver:
  type: nse_fvm
  use_openmp: true
  mpi_processes: 16
  omp_threads: 1
```

続けて実行環境ルートで次を実行する。

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

Linux（bash）:

```bash
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --run
```

元計算がCUDAまたは別のMPIプロセス数でも、`mpi_processes`は読込み先CPU計算で
使用したい値を指定する。ただし、NSE CPU版が要求するy-z領域分割条件を満たすこと。

### 2.6 読込み先がCUDAの場合

CUDA用の実行環境は、生成前のenvironment設計書で次のように選ぶ。

```yaml
parallel:
  use_mpi: false
  use_openmp: false
  use_cuda: true
```

生成済みCUDA環境の`case.yaml`では、実行時並列数を次のようにする。

```yaml
solver:
  type: nse_fvm
  use_openmp: false
  mpi_processes: 1
  omp_threads: 1
```

実行コマンドはCPU版と同じである。

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --run
```

Linux（bash）:

```bash
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --run
```

CUDA版はCPU上で可搬SLFを読み、初期場を完成させてから単一GPUへ転送する。
元計算のMPIプロセス数には依存しないが、初期場と計算用配列がGPUメモリへ収まる
必要がある。`case.yaml`だけでCPU環境をCUDA環境へ変更せず、必ずCUDA用として
生成した実行環境を使用する。

### 2.7 読込み先がMPI＋CUDAの場合

マルチGPU用のenvironment設計書は次のようにする。

```yaml
parallel:
  use_mpi: true
  use_openmp: false
  use_cuda: true
```

`case.yaml`の`solver.mpi_processes`を使用GPU数に設定する。各rankは可搬SLFから
自分のy-z局所範囲だけをCPU上へ読み、その局所場を担当GPUへ転送する。元計算の
バックエンドやMPIプロセス数と一致させる必要はない。

### 2.8 正常終了を確認する

標準出力の最後に次が表示され、最終stepのSLFがあることを確認する。

```text
NSE calculation completed successfully: step=<final step>, time=<final time>
```

MPI＋CUDA版では`NSE MPI+CUDA calculation completed successfully: ...`と表示される。

## 3. 初期乱流SLF変換の詳細

通常のNSE出力はrank別でghostセルを含むため、そのまま初期条件には使用しない。
出力ディレクトリの`meta.json`と全rankのSLFを次のツールへ渡す。

```powershell
python .\SolverLibrary\NSE\tools\nse_prepare_imported_turbulence.py `
  .\previous_case\output `
  --step latest `
  --output .\cases\caseNNNN\initial_data\turbulence.slf
```

Linux（bash）:

```bash
python3 ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py \
  ./previous_case/output \
  --step latest \
  --output ./cases/caseNNNN/initial_data/turbulence.slf
```

上記は生成された実行環境のルートで実行する例である。`SolverLibrary/NSE`を
カレントディレクトリにしている場合は、Windowsでは`python .\tools\...`、
Linuxでは`python3 ./tools/...`としてもよい。

特定stepを使用する場合は`--step 1000`のように指定する。ツールは次を実行する。

1. 指定stepの全rankファイルが揃っていることを確認する。
2. `meta.json`のrank範囲に従って全体場を復元する。
3. ghostセルを除去する。
4. `rho, rho_u, rho_v, rho_w, rho_E`を検査する。
5. MPI分割に依存しない単一のfloat64 SLFを書き出す。

変換先はメモリマップで作成し、元rankファイルを一つずつ配置するため、全体の
5変数場をRAM上へ同時展開しない。

変換時の`--gamma`は読込み先の`physics.nse.gamma`と同じ値にする。既定値は
`1.4`である。

## 4. embed配置

```yaml
flow:
  type: imported_turbulence
  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: embed
    x_start: 2.0
    blend_cells: 8
    velocity_offset: [0.5, 0.0, 0.0]
    background:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143
```

`x_start`は、乱流ブロック先頭セルの左側境界の物理x座標である。対象格子の
セル境界と一致しなければ実行を停止する。

`file`の相対パスは`case.yaml`があるケースディレクトリを基準に解決される。
ソルバーも同じケースディレクトリを作業ディレクトリとして実行されるため、上の例では
`initial_data/turbulence.slf`がFortran namelistへ渡される。ケースディレクトリ外の
ファイルを指定した場合は、実行場所に依存しない絶対パスへ変換される。

`blend_cells`は乱流ブロックの左右それぞれに設ける混合セル数である。混合係数
`w`にはraised-cosineを使用し、保存変数を次式で混合する。

```text
Q_initial = (1-w) Q_background + w Q_imported
```

`blend_cells: 0`では元データを厳密にコピーする。ただし背景場と乱流場が接続部で
一致しない場合、不連続が初期数値振動を発生させる可能性がある。

`velocity_offset`はSLFから復元した速度へ加える一定速度である。例えば平均速度が
ほぼゼロのHITへ`[U,0,0]`を加え、x方向平均流を持つ乱流塊を作れる。圧力を保持した
まま運動量と全エネルギーを再計算する。

## 5. tile配置

```yaml
flow:
  type: imported_turbulence
  imported_turbulence:
    file: initial_data/turbulence.slf
    mode: tile
    x_start: 0.0
    blend_cells: 0
    velocity_offset: [0.5, 0.0, 0.0]
```

`tile`では対象`nx`が元データ`nx`の整数倍でなければならない。この制約により、
x周期境界をまたぐセル列も元データの周期性と一致する。`x_start`は繰返し位相を
セル単位でずらすために使用できる。

## 6. 格子と物理量の互換条件

読込み時に次を検査し、不一致なら計算開始前に停止する。

- 元データと対象計算の`dx`, `dy`, `dz`が一致する。
- `ny`, `nz`とy-z物理領域が一致する。
- `embed`の乱流ブロック全体が対象x領域内に入る。
- `tile`の対象`nx`が元データ`nx`の整数倍である。
- 全セルで保存量が有限であり、密度が`small_rho`より大きい。
- 保存量から復元した圧力が`small_p`より大きい。

格子補間は行わない。補間やx方向の幾何学的引伸ばしは乱流スペクトルと渦スケールを
変えるため、必要な場合は別の前処理として明示的に実施する。

SLFは保存量を保持しているため、元計算と読込み先では同じ`gamma`と無次元化を使う。

## 7. 境界条件

CPU/MPI/OpenMP版では、`embed`とx方向両端の特性無反射境界を組み合わせられる。
背景状態と無反射境界の基準状態は同じ値にするのが基本である。

```yaml
boundary:
  faces:
    x_min: {type: non_reflecting, reference_state: background}
    x_max: {type: non_reflecting, reference_state: background}
    y_min: {type: periodic}
    y_max: {type: periodic}
    z_min: {type: periodic}
    z_max: {type: periodic}
  reference_states:
    background:
      density: 1.0
      velocity: [0.5, 0.0, 0.0]
      pressure: 0.7142857142857143
  non_reflecting:
    formulation: characteristic_relaxation
    relaxation_strength: 0.1
    length_scale: auto
```

これは初期配置した乱流塊を平均流で領域外へ通過させる用途であり、時間ごとに新しい
乱流を入口から供給する合成乱流流入条件ではない。上のx無反射・y-z周期設定は
CPU/MPI/OpenMP、単一GPU CUDA、MPI＋CUDAの全profileで使用できる。

## 8. 主なエラーメッセージ

| メッセージの要点 | 対処 |
|---|---|
| `prepare a ghost-free SLF` | 付属変換ツールでrank別SLFを可搬SLFへ変換する。 |
| `source and target ny/nz must match` | y-z格子数を元計算と同じにする。 |
| `source and target dx must match` | x格子間隔を元計算と同じにする。 |
| `x_start must lie on a target cell boundary` | `x_min + n*dx`となる値を指定する。 |
| `tile mode requires target nx ...` | 対象`nx`を元データ`nx`の整数倍にする。 |
| `density/pressure is below small_*` | 元データ、gamma、無次元化を確認する。 |

## 9. 保存乱流へ平面衝撃波を入射する場合

`flow.type: shock_turbulence_interaction`を選ぶと、`embed`配置した乱流の外側に
平面衝撃波を初期化できる。進行方向上流側のx面には衝撃波背後状態のDirichlet境界を
設定する。初期データを作ったCPU/CUDA方式やMPIプロセス数と、干渉計算側の方式・
プロセス数は一致させる必要がない。入力例、Rankine–Hugoniot関係、配置制約は
[`NSE_SHOCK_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TURBULENCE_INTERACTION.md)にまとめる。

有限長の高圧室から衝撃波と膨張波を発生させる場合は
`flow.type: shock_tube_turbulence_interaction`を選ぶ。配置と境界条件は
[`NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md`](NSE_SHOCK_TUBE_TURBULENCE_INTERACTION.md)にまとめる。
