# Windows／Linuxコマンド対応表

更新日: 2026年9月2日
対象: FrameWorkの取得、開発、実行環境生成、ビルド、実行、後処理

## 1. この文書の使い方

FrameWorkの手順書では、OS固有の操作を次の見出しで併記します。

- **Windows（PowerShell）**: PowerShell 7またはWindows PowerShellで実行する。
- **Linux（bash）**: bashで実行する。スパコンでは必要なmoduleを先に読み込む。
- **Windows／Linux共通**: GitやCMakeなど、引数とパスがOSに依存しないコマンド。

Linuxでは原則として`python3`を使用します。環境によって`python`がPython 3を指す
ことを確認できる場合だけ、`python3`を`python`へ読み替えて構いません。

## 2. 基本的な読み替え

| 操作 | Windows（PowerShell） | Linux（bash） |
|---|---|---|
| ホームディレクトリ | `$env:USERPROFILE` | `$HOME` |
| ディレクトリ移動 | `Set-Location C:\Research\FrameWork` | `cd "$HOME/Research/FrameWork"` |
| ディレクトリ作成 | `New-Item -ItemType Directory -Force C:\Research` | `mkdir -p "$HOME/Research"` |
| ファイルコピー | `Copy-Item source destination` | `cp source destination` |
| ディレクトリコピー | `Copy-Item source destination -Recurse` | `cp -a source destination` |
| 存在確認 | `Test-Path .\path` | `test -e ./path` |
| 一覧 | `Get-ChildItem -Force` | `ls -la` |
| 環境変数設定 | `$env:OMP_NUM_THREADS = "4"` | `export OMP_NUM_THREADS=4` |
| 複数行の継続 | 行末のバッククォート `` ` `` | 行末のバックスラッシュ `\` |
| 相対パス | `.\tools\run_case.py` | `./tools/run_case.py` |
| 実行ファイル | `.\solver.exe` | `./solver` |

パスに空白が含まれる可能性がある場合は、WindowsとLinuxのどちらでもパス全体を
二重引用符で囲みます。Linuxではファイル名の大文字と小文字が区別されます。

## 3. リポジトリを取得する

### Windows（PowerShell）

```powershell
New-Item -ItemType Directory -Path C:\Research -Force
Set-Location C:\Research
git clone https://github.com/chanken0901/framework.git FrameWork
Set-Location .\FrameWork
git status
```

### Linux（bash）

```bash
mkdir -p "$HOME/Research"
cd "$HOME/Research"
git clone https://github.com/chanken0901/framework.git FrameWork
cd FrameWork
git status
```

## 4. 作業開始時にmainを更新する

次のGitコマンド自体は両OS共通です。FrameWorkルートで実行します。

```console
git switch main
git pull --ff-only origin main
git status
```

fork運用で`upstream`から更新する場合も両OS共通です。

```console
git switch main
git fetch upstream
git merge --ff-only upstream/main
git push origin main
```

## 5. Python環境を準備する

### Windows（PowerShell）

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install pyyaml numpy
```

### Linux（bash）

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install pyyaml numpy
```

## 6. 外部実行環境を生成する

以下では、FrameWorkを`$HOME/Research/FrameWork`、計算環境を
`$HOME/ResearchRuns`へ置くLinux例を示します。設計書の`select.target`はLinuxでは
`linux_gnu_mpi`または利用機関用プロファイル、Slurm環境では`linux_hpc_slurm`を
選択します。

### Windows（PowerShell）

```powershell
$framework = "C:\Users\Owner\Documents\Codex\FrameWork"
$tool = "$framework\ScriptLibrary\RunEnvironment"
$designs = "$env:USERPROFILE\ResearchRuns\Designs"

New-Item -ItemType Directory -Path $designs -Force
Copy-Item "$tool\environment.nse.yaml" "$designs\nse_case0001.yaml"
python "$tool\prepare_environment.py" "$designs\nse_case0001.yaml" --dry-run
python "$tool\prepare_environment.py" "$designs\nse_case0001.yaml"
```

### Linux（bash）

```bash
framework="$HOME/Research/FrameWork"
tool="$framework/ScriptLibrary/RunEnvironment"
designs="$HOME/ResearchRuns/Designs"

mkdir -p "$designs"
cp "$tool/environment.nse.yaml" "$designs/nse_case0001.yaml"
```

コピーした設計書を開き、少なくとも次の項目をLinux用へ変更します。

```yaml
select:
  source: framework_relative
  destination: local_generated
  target: linux_gnu_mpi

source:
  framework_root: ${HOME}/Research/FrameWork

destination:
  root: ${HOME}/ResearchRuns
  case_index: ${HOME}/ResearchRuns/case_index.csv
  auto_case_number: true
```

変更後に生成します。

```bash
python3 "$tool/prepare_environment.py" "$designs/nse_case0001.yaml" --dry-run
python3 "$tool/prepare_environment.py" "$designs/nse_case0001.yaml"
```

生成先は設計書の`select.destination`と`destination`上書きで決まります。Linuxでは
Windows用の`local_framework`、`windows_research_runs`、`windows_gnu_msmpi`を残さないでください。

## 7. 入力生成、検証、ビルド、テスト、実行

生成済み実行環境のルートへ移動して実行します。

### Windows（PowerShell）

```powershell
python .\tools\run_case.py --prepare
python .\tools\run_case.py --validate-only
python .\tools\run_case.py --build
python .\tools\run_case.py --test
python .\tools\run_case.py --run
```

### Linux（bash）

```bash
python3 ./tools/run_case.py --prepare
python3 ./tools/run_case.py --validate-only
python3 ./tools/run_case.py --build
python3 ./tools/run_case.py --test
python3 ./tools/run_case.py --run
```

MPIプロセス数とOpenMPスレッド数を一時的に指定する場合は次のとおりです。

### Windows（PowerShell）

```powershell
python .\tools\run_case.py --run --processes 8 --omp-threads 2
```

### Linux（bash）

```bash
python3 ./tools/run_case.py --run --processes 8 --omp-threads 2
```

## 8. CMakeを直接使用する

通常は`run_case.py`または`build_model.py`を使用します。直接構成する場合、Windowsの
Visual Studioはmulti-config、一般的なLinuxのNinja/Makefilesはsingle-configである点に
注意してください。

### Windows（PowerShell、Ninja例）

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8
ctest --test-dir build --output-on-failure
```

### Linux（bash、Ninja例）

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel "$(nproc)"
ctest --test-dir build --output-on-failure
```

既存presetを使うNSE MPI版の例です。

### Windows（PowerShell）

```powershell
cmake --preset windows-msmpi-release
cmake --build --preset windows-msmpi-release
```

### Linux（bash）

```bash
cmake --preset linux-mpi-release -DCMAKE_Fortran_COMPILER=mpifort
cmake --build --preset linux-mpi-release
```

## 9. MPI、OpenMP、CUDA

手動起動する場合の代表例です。生成環境では`run_case.py --run`を優先してください。

### Windows（PowerShell）

```powershell
$env:OMP_NUM_THREADS = "2"
mpiexec -n 8 .\solver.exe .\input.dat
```

### Linux（bash）

```bash
export OMP_NUM_THREADS=2
mpirun -np 8 ./solver ./input.dat
```

LinuxのMPIランチャーは環境によって`mpirun`、`mpiexec`、またはSlurmの`srun`です。
コンパイルに用いたMPI実装と同じランチャーを使用します。MPI＋CUDAは原則として
1 MPI rankを1 GPUへ割り当て、利用機関の方法で`CUDA_VISIBLE_DEVICES`を設定します。

Slurmでは、生成された`submit.slurm`を確認してから投入します。

```bash
sbatch --nodes=1 --ntasks-per-node=8 --cpus-per-task=2 submit.slurm
squeue -u "$USER"
job_id=123456  # 実際のSlurm job IDへ変更する
tail -f "slurm-${job_id}.out"
```

## 10. 後処理・ParaView変換

### Windows（PowerShell）

```powershell
python .\tools\postprocess_case.py --inspect-only --steps latest
python .\tools\postprocess_case.py --steps latest --stride 1
```

### Linux（bash）

```bash
python3 ./tools/postprocess_case.py --inspect-only --steps latest
python3 ./tools/postprocess_case.py --steps latest --stride 1
```

GUIを使えない計算ノードでは変換のみ行い、生成されたVTI/PVDをParaViewが動作する
端末へ転送します。

## 11. 保存済み乱流データを準備する

### Windows（PowerShell）

```powershell
python .\SolverLibrary\NSE\tools\nse_prepare_imported_turbulence.py `
  .\previous_case\output `
  --step latest `
  --output .\cases\caseNNNN\initial_data\turbulence.slf
```

### Linux（bash）

```bash
python3 ./SolverLibrary/NSE/tools/nse_prepare_imported_turbulence.py \
  ./previous_case/output \
  --step latest \
  --output ./cases/caseNNNN/initial_data/turbulence.slf
```

## 12. トラブル調査用コマンド

### Windows（PowerShell）

```powershell
Get-Location
Get-ChildItem -Force
Test-Path .\environment.lock.json
git status --short --branch
python --version
cmake --version
```

### Linux（bash）

```bash
pwd
ls -la
test -e ./environment.lock.json && echo "environment.lock.json: OK"
git status --short --branch
python3 --version
cmake --version
```

Linuxで権限エラーが出た場合、まず所有者と権限を確認します。安易に`sudo`でビルドや
実行を行うと生成物の所有者がrootになり、その後の更新に失敗するため避けてください。

```bash
id
ls -ld . .git build 2>/dev/null
```
