# 後処理・ParaView可視化手順書

**版:** 2.2
**更新日:** 2026-09-02
**対象:** 研究フレームワークから生成したGPE/NSE実行環境  
**正本:** `C:\Users\Owner\Documents\Codex\FrameWork`

> コマンド表記: OS固有の操作はWindows（PowerShell）とLinux（bash）を併記します。共通の読み替えは[`docs/WINDOWS_LINUX_COMMANDS.md`](docs/WINDOWS_LINUX_COMMANDS.md)を参照してください。

---

## 1. 目的

本書は、ソルバーが出力した三次元SLFデータをParaView用のVTI/PVDへ変換し、
時系列として可視化する標準手順を説明します。逐次CPU、MPI、単一GPUの出力を
同じ入口から処理できます。

現行の後処理ツールはバージョン2.0.0です。主な機能は次のとおりです。

- GPE/NSE共通入口による変換
- global形式とMPI rank分割形式の自動判定
- `meta.json`に基づくrank結合とghost cell除去
- 変換する保存ステップの選択
- 出力変数の選択と空間間引き
- GPEの`density`、`phase`、`abs_psi`導出
- NSEの`u`、`v`、`w`、`p`導出
- NSEの`keep2`、`keep6`、`weno5z_roe`、`hybrid`に共通の変換手順
- 計算中に未完成のMPIステップを除外
- 書き込み前の`--inspect-only`検査
- VTI時系列をまとめる`collection.pvd`生成

ParaView本体のバージョンと変換器のバージョンは別物です。本書の「2.0.0」は
Python後処理ツールの版を指します。

NSEの対流流束を切り替えてもSLFの保存変数と後処理の呼び出し方は共通です。
`case.yaml`から生成した`input.dat`で計算し、完全に書き込まれた保存ステップを
本手順で変換します。

---

## 2. 入出力

入力はcaseの`output`ディレクトリにあるSLFと`meta.json`です。

```text
<run_root>\
├─ environment.lock.json
├─ tools\
│  └─ postprocess_case.py
└─ cases\
   └─ <case_id>\
      ├─ case.yaml
      ├─ output\
      │  ├─ meta.json
      │  └─ field_*.slf
      └─ paraview\
         ├─ collection.pvd
         └─ field_*.vti
```

出力の役割は次のとおりです。

| ファイル | 役割 |
|---|---|
| `field_XXXXXX.vti` | 1保存ステップの三次元場 |
| `collection.pvd` | 複数VTIを時系列として参照する一覧 |

ParaViewでは個々のVTIではなく、通常は`collection.pvd`を開きます。

---

## 3. 標準実行手順

### 3.1 実行環境ルートへ移動する

対象の実行環境ルートへ移動します。以下では`nse_case0012`を例にします。

Windows（PowerShell）:

```powershell
Set-Location "$env:USERPROFILE\ResearchRuns\nse_case0012"
```

現在位置を確認します。

```powershell
Get-Location
Test-Path .\environment.lock.json
Test-Path .\tools\postprocess_case.py
```

後ろの二つがどちらも`True`である必要があります。`environment.lock.json`がない
ディレクトリから実行すると、modelとcaseを正しく特定できません。

Linux（bash）:

```bash
cd "$HOME/ResearchRuns/nse_case0012"
pwd
test -f ./environment.lock.json && echo "environment.lock.json: OK"
test -f ./tools/postprocess_case.py && echo "postprocess_case.py: OK"
```

### 3.2 ツールの版を確認する

```powershell
python .\tools\postprocess_case.py --version
```

Linux（bash）:

```bash
python3 ./tools/postprocess_case.py --version
```

現行版では次のように表示されます。

```text
postprocess_case.py 2.0.0
```

`--version`が認識されない、または2.0.0未満の場合は、14章の手順で実行環境を
現行FrameWorkから再生成します。

### 3.3 実行予定コマンドだけを確認する

```powershell
python .\tools\postprocess_case.py --dry-run
```

Linux（bash）:

```bash
python3 ./tools/postprocess_case.py --dry-run
```

`--dry-run`は、内部で呼び出すモデル別変換器のコマンドを`[CMD]`行に表示します。
SLFは読み込まず、出力ファイルも作りません。case、入力、出力、モデルの解決結果を
確認したいときに使用します。

### 3.4 SLFと選択条件を検査する

```powershell
python .\tools\postprocess_case.py --inspect-only --steps latest
```

Linux（bash）:

```bash
python3 ./tools/postprocess_case.py --inspect-only --steps latest
```

`--inspect-only`はSLFヘッダーと`meta.json`を読み、次を表示します。

- 使用する変換器の版
- 入力ディレクトリと`meta.json`
- globalまたはrank-wiseの判定
- 検出したSLFファイル数
- 変換対象となる保存ステップ
- SLF内の一次変数名
- GPE/NSE導出モード
- 出力変数、全体格子数、空間間引き

VTI/PVDはまだ生成しません。`--dry-run`はコマンド確認、`--inspect-only`は
実データの検査、という違いがあります。

### 3.5 軽量プレビューを生成する

```powershell
python .\tools\postprocess_case.py
```

Linux（bash）:

```bash
python3 ./tools/postprocess_case.py
```

引数を省略した場合は次の条件です。

| 項目 | 既定値 |
|---|---|
| task | `paraview` |
| steps | `latest` |
| stride | `2` |
| layout | `auto` |
| NSE fields | `rho,u,v,w,p` |
| GPE fields | `density,phase` |
| output | `cases\<case_id>\paraview` |

最新の完全な保存ステップだけを各軸2点おきに変換するため、まず結果を確認する用途に
向いています。三次元格子点数は概ね`1/8`になります。

### 3.6 ParaViewで開く

```powershell
$lock = Get-Content .\environment.lock.json | ConvertFrom-Json
$pvd = Join-Path $lock.case_directory "paraview\collection.pvd"
Invoke-Item $pvd
```

ファイルダイアログから開く場合は、次を選択します。

```text
cases\<case_id>\paraview\collection.pvd
```

LinuxでParaViewを直接起動できる場合は次のように開きます。

```bash
case_id=case0012
paraview "./cases/$case_id/paraview/collection.pvd"
```

ParaViewのPropertiesで`Apply`を押し、上部の時刻操作ボタンまたはTime欄で
保存時刻を移動します。

---

## 4. 本番変換の例

### 4.1 全解像度で最新時刻を変換する

```powershell
python .\tools\postprocess_case.py --steps latest --stride 1
```

### 4.2 指定した保存ステップだけ変換する

```powershell
python .\tools\postprocess_case.py --steps 0,500,1000 --stride 1
```

### 4.3 範囲と間隔を指定する

```powershell
python .\tools\postprocess_case.py --steps 0:1000:100 --stride 1
```

`0:1000:100`は0から1000までを100刻みで選びます。終端1000も対象に含みます。
存在しない保存ステップを番号一覧で直接指定するとエラーになります。

### 4.4 全保存ステップを変換する

```powershell
python .\tools\postprocess_case.py --steps all --stride 1
```

大規模三次元計算ではVTIの容量と変換時間が大きくなります。先に`latest`と
`stride 2`または`stride 4`で確認してから、必要な保存ステップだけを全解像度で
変換してください。

### 4.5 出力先を分ける

```powershell
python .\tools\postprocess_case.py `
  --steps 0:1000:100 `
  --stride 1 `
  --output-dir .\cases\case0012\postprocess\paraview\full
```

Linux（bash）:

```bash
python3 ./tools/postprocess_case.py \
  --steps 0:1000:100 \
  --stride 1 \
  --output-dir ./cases/case0012/postprocess/paraview/full
```

既存プレビューを残したまま別条件で変換したい場合に使用します。

---

## 5. GPEの可視化

標準変換は次の二つです。

```powershell
python .\tools\postprocess_case.py `
  --fields density,phase `
  --steps latest `
  --stride 1
```

利用できる代表的な場は次のとおりです。

| 場 | 定義 | 用途 |
|---|---|---|
| `psi_real` | 波動関数の実部 | 波動関数の確認 |
| `psi_imag` | 波動関数の虚部 | 波動関数の確認 |
| `density` | `psi_real^2 + psi_imag^2` | 密度、渦芯、音波 |
| `phase` | `atan2(psi_imag, psi_real)` | 位相巻き込み |
| `abs_psi` | `sqrt(density)` | 波動関数振幅 |

量子渦を密度で見る場合は`Cell Data to Point Data`を適用し、`Contour`で背景密度の
0.1から0.3倍程度を試します。値はcaseによって異なるため、Histogramとデータ範囲を
確認して決めます。低密度領域では位相が不定になるため、`phase`だけで渦を判定せず、
同じ位置の`density`も確認してください。

---

## 6. NSEの可視化

標準変換は次の五つです。

```powershell
python .\tools\postprocess_case.py `
  --fields rho,u,v,w,p `
  --steps latest `
  --stride 1
```

| 場 | 内容 |
|---|---|
| `rho` | 密度 |
| `u`,`v`,`w` | 保存運動量を密度で割った速度成分 |
| `p` | 全エネルギーから導出した圧力 |

速度ベクトルが必要な場合はParaViewの`Calculator`を使用します。

```text
u*iHat + v*jHat + w*kHat
```

結果名を`velocity`とします。その後、`Gradient`のVorticity、`Stream Tracer`、
`Glyph`などへ渡せます。渦構造はVorticityだけでなく、必要に応じてQ criterionも
併用してください。

---

## 7. MPI出力の結合

逐次CPUと単一GPUは通常、一つの保存ステップにつき一つのglobal SLFを出力します。
MPIでは一つの保存ステップにつきrank数分のSLFを出力します。

| 形式 | 代表ファイル名 | layout |
|---|---|---|
| global | `field_000000.slf` | `global` |
| rank分割 | `field_000000_rank00000.slf` | `rank` |

通常は`--layout auto`を使用します。rank分割の結合には`meta.json`内の
`parallel.rank_ranges`が必要です。変換器は各rankの担当範囲を使って全領域を組み、
必要に応じてghost cellを中央から除去します。

計算中に変換した場合、ある保存ステップのrankファイルがまだ揃っていなければ、その
ステップを警告付きで除外します。`latest`は「見つかった最大番号」ではなく、
「全rankが揃った最大番号」から選ばれます。

`--layout rank`や`--layout global`は、同じ出力先に両形式が混在していて自動判定を
上書きしたい場合だけ使用します。

---

## 8. オプション一覧

```powershell
python .\tools\postprocess_case.py --help
```

| オプション | 内容 |
|---|---|
| `--version` | 共通後処理ツールの版を表示 |
| `--task paraview` | ParaView変換。既定値 |
| `--task statistics` | NSE乱流統計CSVを生成 |
| `--task all` | ParaView変換とNSE乱流統計を順に実行 |
| `--case-directory PATH` | lockファイルのcase位置を一時上書き |
| `--input-dir PATH` | SLF入力ディレクトリを上書き |
| `--output-dir PATH` | VTI/PVD出力ディレクトリを上書き |
| `--meta PATH` | 使用する`meta.json`を上書き |
| `--steps VALUE` | `all`、`latest`、一覧、範囲を指定 |
| `--fields LIST` | 出力する変数をカンマ区切りで指定 |
| `--stride N` | 各軸でN点おきに保持。1は全解像度 |
| `--derive MODE` | `auto`、`none`、`nse`、`gpe` |
| `--layout MODE` | `auto`、`global`、`rank` |
| `--gamma VALUE` | NSE圧力導出用gammaを一時上書き |
| `--pvd-name NAME` | PVDファイル名を変更 |
| `--dry-run` | 子コマンドだけ表示。SLFを読まない |
| `--inspect-only` | SLFと選択条件を検査。VTI/PVDを書かない |

PowerShellで複数行に分ける場合は、行末のバッククォートの後ろに空白を置かないで
ください。貼り付けが不安定な場合は一行で実行できます。

---

## 9. ParaViewでの基本操作

1. `collection.pvd`を開き、Propertiesの`Apply`を押します。
2. Coloringから`density`、`rho`、`p`などを選びます。
3. `Rescale to Data Range`で表示範囲をデータに合わせます。
4. 内部断面は`Slice`、等値面は`Contour`、領域抽出は`Threshold`を使います。
5. Cell Dataでは直接Contourできない場合があるため、`Cell Data to Point Data`を先に適用します。
6. 時系列では上部の再生ボタンとTime欄を使用します。
7. 時刻比較ではAutomatic Rescaleを無効にし、全時刻で同じ色範囲を固定します。
8. 表示設定を再利用する場合は`File > Save State`で`.pvsm`を保存します。

`Surface`表示は外表面しか見えません。三次元内部構造を確認したい場合は`Slice`、
`Clip`、`Contour`、`Volume`を使用してください。

---

## 10. 大規模データの扱い

VTIは後処理用の派生データです。SLFを正本として保持し、必要に応じてVTIを再生成します。

推奨順序は次のとおりです。

1. `--inspect-only --steps latest`で入力を検査する。
2. `--steps latest --stride 4`で軽量確認する。
3. 必要な変数だけ`--fields`で選ぶ。
4. 必要な保存ステップだけ`--steps`で選ぶ。
5. 最終画像に必要なデータだけ`--stride 1`で変換する。

空間間引きの目安です。

| stride | 格子点数の概算 | 用途 |
|---:|---:|---|
| 1 | 100% | 最終解析 |
| 2 | 12.5% | 通常プレビュー |
| 4 | 1.56% | 大規模caseの初期確認 |

VTIの容量は変数数にも比例します。GPEで密度だけ必要なら`--fields density`、NSEで
速度だけ必要なら`--fields u,v,w`とします。

---

## 11. 直接変換器を使う場合

通常は`tools\postprocess_case.py`を使用してください。モデル別変換器の直接実行は、
実行環境のlockを使わず入力と出力を明示したい場合や、障害調査向けです。

GPEの例です。

```powershell
python .\SolverLibrary\GPE\gp3d\tools\slf_to_paraview_merged_cropghost.py `
  .\cases\case0001\output `
  --output-dir .\cases\case0001\paraview `
  --meta .\cases\case0001\output\meta.json `
  --derive gpe `
  --layout auto `
  --fields density,phase `
  --steps latest `
  --stride 2
```

Linux（bash）:

```bash
python3 ./SolverLibrary/GPE/gp3d/tools/slf_to_paraview_merged_cropghost.py \
  ./cases/case0001/output \
  --output-dir ./cases/case0001/paraview \
  --meta ./cases/case0001/output/meta.json \
  --derive gpe \
  --layout auto \
  --fields density,phase \
  --steps latest \
  --stride 2
```

NSEの例です。

```powershell
python .\SolverLibrary\NSE\tools\slf_to_paraview_merged_cropghost.py `
  .\cases\case0012\output `
  --output-dir .\cases\case0012\paraview `
  --meta .\cases\case0012\output\meta.json `
  --derive nse `
  --layout auto `
  --fields rho,u,v,w,p `
  --steps latest `
  --stride 2
```

Linux（bash）:

```bash
python3 ./SolverLibrary/NSE/tools/slf_to_paraview_merged_cropghost.py \
  ./cases/case0012/output \
  --output-dir ./cases/case0012/paraview \
  --meta ./cases/case0012/output/meta.json \
  --derive nse \
  --layout auto \
  --fields rho,u,v,w,p \
  --steps latest \
  --stride 2
```

直接変換器も`--version`と`--inspect-only`に対応しています。

---

## 12. トラブルシューティング

### 12.1 `No .slf files found`

```powershell
$lock = Get-Content .\environment.lock.json | ConvertFrom-Json
$output = Join-Path $lock.case_directory "output"
Get-ChildItem $output -Filter *.slf | Select-Object -First 10 Name,Length
Test-Path (Join-Path $output "meta.json")
```

Linux（bash）:

```bash
output="$(python3 -c 'import json; p=json.load(open("environment.lock.json")); print(p["case_directory"] + "/output")')"
find "$output" -maxdepth 1 -type f -name '*.slf' -printf '%f %s bytes\n' | head
test -f "$output/meta.json" && echo "meta.json: OK"
```

SLFの実際の出力先と`environment.lock.json`の`case_directory`を確認します。

### 12.2 `unrecognized arguments: --steps`または`--inspect-only`

生成済み実行環境の変換器が古い状態です。`python .\tools\postprocess_case.py --version`
を確認し、14章に従って実行環境を再生成します。

### 12.3 `Rank-wise SLF conversion requires meta.json`

MPI rank分割SLFには正しい`output\meta.json`が必要です。`--meta`で別caseの
ファイルを指定していないか確認してください。

### 12.4 `No complete time steps are available`

全rankのファイルがまだ揃っていない、`meta.json`の`mpi_nprocs`が実行時rank数と
異なる、または出力が途中で中断しています。各保存ステップのrankファイル数を確認します。

```powershell
Get-ChildItem .\cases\<case_id>\output -Filter *rank*.slf |
  Group-Object { $_.BaseName -replace '_rank\d+$','' } |
  Select-Object Name,Count
```

Linux（bash）:

```bash
case_id=case0012
find "./cases/$case_id/output" -maxdepth 1 -type f -name '*rank*.slf' -printf '%f\n' |
  sed -E 's/_rank[0-9]+\.slf$//' |
  sort |
  uniq -c
```

### 12.5 `Requested steps are not available`

`--inspect-only --steps all`で利用可能な保存ステップを確認し、その番号を指定します。
ソルバーの計算stepとSLF保存stepは、出力間隔によって一致しない場合があります。

### 12.6 `cannot be cropped to meta rank range`

SLFと`meta.json`の格子、rank分割、ghost幅が一致していません。別実行の`meta.json`を
混ぜていないか、実行途中にMPIプロセス数を変えて同じ出力先へ追記していないかを
確認します。

### 12.7 ParaViewで何も見えない

- Propertiesで`Apply`を押す。
- `Reset Camera`を押す。
- `Rescale to Data Range`を実行する。
- SurfaceではなくSliceまたはContourを使う。
- strideと格子数を`--inspect-only`で確認する。

### 12.8 時系列の一部が欠ける

計算中の未完成MPIステップは意図的に除外されます。計算終了後に同じ変換を再実行して
ください。PVDとVTIを別々に移動すると参照が切れるため、出力ディレクトリ単位で扱います。

### 12.9 NumPyがない

Windows（PowerShell）:

```powershell
python -m pip install numpy
```

Linux（bash）:

```bash
python3 -m pip install --user numpy
```

スパコンではシステムPython、venv、module環境のいずれを使用するか運用方針に合わせます。

---

## 13. Linux・スパコンでの実行

実行環境ルートで次を実行します。

```bash
python3 tools/postprocess_case.py --version
python3 tools/postprocess_case.py --dry-run
python3 tools/postprocess_case.py --inspect-only --steps latest
python3 tools/postprocess_case.py --steps 0:1000:100 --stride 1
```

計算ノードよりログインノードまたは後処理ノードが適切な場合があります。大規模データを
ローカルPCへ転送する前に、スパコン側で必要な時刻と変数だけをVTIへ変換すると転送量を
抑えられます。ただし施設のログインノード利用規則に従ってください。

---

## 14. 古い実行環境の更新

実行環境内の`SolverLibrary`、`ScriptLibrary`、`tools`は、環境生成時にFrameWorkから
コピーされたスナップショットです。FrameWork側を更新しても、既存のResearchRuns配下へ
自動反映されません。

現在の版を確認します。

```powershell
python .\tools\postprocess_case.py --version
```

```bash
python3 ./tools/postprocess_case.py --version
```

古い場合は、編集元FrameWorkの`ScriptLibrary\RunEnvironment`から、使用している
`environment.nse.yaml`または`environment.gpe.yaml`に従って実行環境を再生成します。
caseの`case.yaml`と生データ`output`は、再生成前に退避または保持方法を確認してください。
同じ出力先を上書きするか、新しいcase番号へ生成するかは環境設計書の設定に従います。

再生成後、次をもう一度確認します。

```powershell
Test-Path .\environment.lock.json
python .\tools\postprocess_case.py --version
python .\tools\postprocess_case.py --inspect-only --steps latest
```

Linux（bash）:

```bash
test -f ./environment.lock.json && echo "environment.lock.json: OK"
python3 ./tools/postprocess_case.py --version
python3 ./tools/postprocess_case.py --inspect-only --steps latest
```

---

## 15. 最終確認チェックリスト

### 変換前

- [ ] 実行環境ルートにいる
- [ ] `environment.lock.json`と`case.yaml`がある
- [ ] 後処理ツールが2.0.0以上である
- [ ] `output`にSLFと`meta.json`がある
- [ ] `--dry-run`の入力caseと出力先が正しい
- [ ] `--inspect-only`の格子数、時刻、変数が正しい
- [ ] 必要な空き容量がある

### 変換後

- [ ] `field_*.vti`が生成された
- [ ] `collection.pvd`が生成された
- [ ] ParaViewで期待する保存時刻を選べる
- [ ] 座標、格子方向、物理変数が正しい
- [ ] 全時刻比較用の色範囲を固定した
- [ ] 解析に使ったコマンドとContour値を記録した

---

## 16. コマンド早見表

### Windows（PowerShell）

```powershell
# 版確認
python .\tools\postprocess_case.py --version

# 子コマンドだけ確認
python .\tools\postprocess_case.py --dry-run

# 入力と最新ステップを検査
python .\tools\postprocess_case.py --inspect-only --steps latest

# 最新の軽量プレビュー
python .\tools\postprocess_case.py

# 最新の全解像度
python .\tools\postprocess_case.py --steps latest --stride 1

# 0から1000まで100刻み、全解像度
python .\tools\postprocess_case.py --steps 0:1000:100 --stride 1

# GPEの密度と位相
python .\tools\postprocess_case.py --fields density,phase --stride 1

# NSEの密度、速度、圧力
python .\tools\postprocess_case.py --fields rho,u,v,w,p --stride 1
```

### Linux（bash）

```bash
# 版確認
python3 ./tools/postprocess_case.py --version

# 子コマンドだけ確認
python3 ./tools/postprocess_case.py --dry-run

# 入力と最新ステップを検査
python3 ./tools/postprocess_case.py --inspect-only --steps latest

# 最新の軽量プレビュー
python3 ./tools/postprocess_case.py

# 最新の全解像度
python3 ./tools/postprocess_case.py --steps latest --stride 1

# 0から1000まで100刻み、全解像度
python3 ./tools/postprocess_case.py --steps 0:1000:100 --stride 1

# GPEの密度と位相
python3 ./tools/postprocess_case.py --fields density,phase --stride 1

# NSEの密度、速度、圧力
python3 ./tools/postprocess_case.py --fields rho,u,v,w,p --stride 1
```
