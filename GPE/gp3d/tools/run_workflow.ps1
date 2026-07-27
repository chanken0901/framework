[CmdletBinding()]
<#
ケース生成器が作ったworkflow.jsonを読み、CMake構成、ビルド、CTest、計算実行を
順番に行う。CUDA版ではnvccのホストコンパイラとして必要なVisual Studio x64環境を
自動で取り込む。-DryRunを指定すると、外部コマンドを実行せず内容だけを表示する。
#>
param(
  [Parameter(Position = 0)]
  [string]$ConfigFile = "",

  [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Get-FullPath {
  param(
    [Parameter(Mandatory = $true)][string]$BaseDirectory,
    [Parameter(Mandatory = $true)][string]$Path
  )

  if ([System.IO.Path]::IsPathRooted($Path)) {
    return [System.IO.Path]::GetFullPath($Path)
  }
  return [System.IO.Path]::GetFullPath((Join-Path $BaseDirectory $Path))
}

function Convert-ToCMakeValue {
  param([Parameter(Mandatory = $true)]$Value)

  if ($Value -is [bool]) {
    if ($Value) { return "ON" }
    return "OFF"
  }
  return [string]$Value
}

function Format-CommandLine {
  param(
    [Parameter(Mandatory = $true)][string]$Command,
    [Parameter(Mandatory = $true)][string[]]$Arguments
  )

  $formatted = @($Command)
  foreach ($argument in $Arguments) {
    if ($argument -match '[\s"]') {
      $formatted += '"' + $argument.Replace('"', '\"') + '"'
    } else {
      $formatted += $argument
    }
  }
  return $formatted -join " "
}

function Invoke-CheckedCommand {
  param(
    [Parameter(Mandatory = $true)][string]$Command,
    [Parameter(Mandatory = $true)][object[]]$Arguments,
    [switch]$WhatIfOnly
  )

  [string[]]$stringArguments = @($Arguments | ForEach-Object { [string]$_ })
  Write-Host ("+ " + (Format-CommandLine -Command $Command -Arguments $stringArguments))
  if ($WhatIfOnly) { return }

  & $Command @stringArguments
  if ($LASTEXITCODE -ne 0) {
    throw "Command failed with exit code ${LASTEXITCODE}: $Command"
  }
}

function Get-DefinitionValue {
  param(
    [Parameter(Mandatory = $true)]$Definitions,
    [Parameter(Mandatory = $true)][string]$Name,
    [Parameter(Mandatory = $true)]$DefaultValue
  )

  $property = $Definitions.PSObject.Properties[$Name]
  if ($null -eq $property) { return $DefaultValue }
  return $property.Value
}

function Test-EnabledValue {
  param([Parameter(Mandatory = $true)]$Value)

  $text = ([string]$Value).Trim().ToUpperInvariant()
  return @("1", "ON", "TRUE", "YES") -contains $text
}

function Import-VisualStudioEnvironment {
  if ((Get-Command cl.exe -ErrorAction SilentlyContinue) -and
      (Get-Command lib.exe -ErrorAction SilentlyContinue)) {
    return
  }

  $vcvarsPath = ""
  $programFilesX86 = ${env:ProgramFiles(x86)}
  if (-not [string]::IsNullOrWhiteSpace($programFilesX86)) {
    $vswherePath = Join-Path $programFilesX86 `
      "Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path -LiteralPath $vswherePath -PathType Leaf) {
      $installationPath = (& $vswherePath -latest -products * `
        -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
        -property installationPath | Select-Object -First 1)
      if (-not [string]::IsNullOrWhiteSpace([string]$installationPath)) {
        $candidate = Join-Path ([string]$installationPath) `
          "VC\Auxiliary\Build\vcvars64.bat"
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
          $vcvarsPath = $candidate
        }
      }
    }
  }

  if ([string]::IsNullOrWhiteSpace($vcvarsPath)) {
    throw "CUDA build requires the Visual Studio x64 C++ build tools, but vcvars64.bat was not found"
  }

  $commandLine = "call `"$vcvarsPath`" >nul && set"
  $environmentLines = & $env:ComSpec /d /c $commandLine
  if ($LASTEXITCODE -ne 0) {
    throw "Failed to initialize the Visual Studio x64 build environment"
  }
  # Some Windows launchers expose both PATH and Path.  vcvars64.bat writes the
  # augmented value as PATH, followed by the original Path value.  Treat names
  # case-insensitively so the latter cannot erase the MSVC tool directories.
  $importedNames = @{}
  foreach ($line in $environmentLines) {
    if ($line -match '^([^=]+)=(.*)$') {
      $name = $matches[1]
      if (-not $importedNames.ContainsKey($name)) {
        [Environment]::SetEnvironmentVariable($name, $matches[2], "Process")
        $importedNames[$name] = $true
      }
    }
  }
  if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue) -or
      -not (Get-Command lib.exe -ErrorAction SilentlyContinue)) {
    throw "Visual Studio environment was loaded, but cl.exe or lib.exe is not on PATH"
  }
  Write-Host "  MSVC   : $vcvarsPath"
}

try {
  $scriptDirectory = Split-Path -Parent $MyInvocation.MyCommand.Path
  $defaultConfig = Join-Path (Split-Path -Parent $scriptDirectory) "workflow.json"
  if ([string]::IsNullOrWhiteSpace($ConfigFile)) {
    $configPath = $defaultConfig
  } else {
    $configPath = Get-FullPath -BaseDirectory (Get-Location).Path -Path $ConfigFile
  }

  if (-not (Test-Path -LiteralPath $configPath -PathType Leaf)) {
    throw "Workflow configuration was not found: $configPath"
  }

  $configPath = (Resolve-Path -LiteralPath $configPath).Path
  $configDirectory = Split-Path -Parent $configPath
  $config = Get-Content -LiteralPath $configPath -Raw | ConvertFrom-Json
  if ([int]$config.schema_version -ne 1) {
    throw "Unsupported workflow schema_version: $($config.schema_version)"
  }

  $projectRoot = Get-FullPath -BaseDirectory $configDirectory -Path ([string]$config.project_root)
  $sourceDirectory = Get-FullPath -BaseDirectory $projectRoot -Path ([string]$config.cmake.source_directory)
  $buildDirectory = Get-FullPath -BaseDirectory $projectRoot -Path ([string]$config.cmake.build_directory)
  $workingDirectory = Get-FullPath -BaseDirectory $projectRoot -Path ([string]$config.run.working_directory)

  if (-not (Test-Path -LiteralPath $sourceDirectory -PathType Container)) {
    throw "CMake source directory was not found: $sourceDirectory"
  }
  if (-not (Test-Path -LiteralPath $workingDirectory -PathType Container)) {
    throw "Run working directory was not found: $workingDirectory"
  }

  $useMpi = Test-EnabledValue (Get-DefinitionValue -Definitions $config.cmake.definitions `
    -Name "USE_MPI" -DefaultValue "OFF")
  $gpuBackend = [string](Get-DefinitionValue -Definitions $config.cmake.definitions `
    -Name "GPU_BACKEND" -DefaultValue "none")
  $useCuda = $gpuBackend.Trim().ToLowerInvariant() -eq "cuda"
  $useCufftmp = $gpuBackend.Trim().ToLowerInvariant() -eq "cufftmp"
  if ($useCuda -and $useMpi) {
    throw "GPU_BACKEND=cuda is currently single-GPU only; set USE_MPI=OFF"
  }
  if ($useCufftmp -and -not $useMpi) {
    throw "GPU_BACKEND=cufftmp requires USE_MPI=ON"
  }
  if ($useCufftmp -and $env:OS -eq "Windows_NT") {
    throw "GPU_BACKEND=cufftmp must be configured and run on Linux"
  }
  if ($useCuda -and $env:OS -eq "Windows_NT") {
    Import-VisualStudioEnvironment
  }
  if ([bool]$config.run.use_mpi_launcher -and -not $useMpi) {
    throw "run.use_mpi_launcher is true, but the CMake definition USE_MPI is not enabled"
  }
  if ([int]$config.run.processes -lt 1) {
    throw "run.processes must be at least one"
  }

  Write-Host "GP3D workflow"
  Write-Host "  config : $configPath"
  Write-Host "  source : $sourceDirectory"
  Write-Host "  build  : $buildDirectory"
  Write-Host "  MPI    : $useMpi"
  Write-Host "  GPU    : $gpuBackend"
  if ($DryRun) { Write-Host "  mode   : dry-run" }

  if ([bool]$config.stages.configure) {
    $configureArguments = @("-S", $sourceDirectory, "-B", $buildDirectory)
    if (-not [string]::IsNullOrWhiteSpace([string]$config.cmake.generator)) {
      $configureArguments += @("-G", [string]$config.cmake.generator)
    }
    foreach ($definition in $config.cmake.definitions.PSObject.Properties) {
      $value = Convert-ToCMakeValue $definition.Value
      $configureArguments += "-D$($definition.Name)=$value"
    }
    $configureArguments += @($config.cmake.configure_arguments)
    Invoke-CheckedCommand -Command ([string]$config.cmake.command) `
      -Arguments $configureArguments -WhatIfOnly:$DryRun
  }

  if ([bool]$config.stages.build) {
    $buildArguments = @("--build", $buildDirectory)
    if (-not [string]::IsNullOrWhiteSpace([string]$config.build.configuration)) {
      $buildArguments += @("--config", [string]$config.build.configuration)
    }
    if (-not [string]::IsNullOrWhiteSpace([string]$config.build.target)) {
      $buildArguments += @("--target", [string]$config.build.target)
    }
    if ([int]$config.build.parallel_jobs -gt 0) {
      $buildArguments += @("--parallel", [string]$config.build.parallel_jobs)
    }
    if ([bool]$config.build.clean_first) {
      $buildArguments += "--clean-first"
    }
    $buildArguments += @($config.build.arguments)
    Invoke-CheckedCommand -Command ([string]$config.cmake.command) `
      -Arguments $buildArguments -WhatIfOnly:$DryRun
  }

  if ([bool]$config.stages.test) {
    $testArguments = @("--test-dir", $buildDirectory)
    if (-not [string]::IsNullOrWhiteSpace([string]$config.test.configuration)) {
      $testArguments += @("--build-config", [string]$config.test.configuration)
    }
    $testArguments += @($config.test.arguments)
    Invoke-CheckedCommand -Command ([string]$config.test.command) `
      -Arguments $testArguments -WhatIfOnly:$DryRun
  }

  if ([bool]$config.stages.run) {
    $executableName = [string]$config.run.executable
    if ($executableName -eq "auto") {
      if ($useCufftmp) {
        $executableName = "gp3d_cufftmp"
      } elseif ($useCuda) {
        $executableName = "gp3d_cuda"
      } elseif ($useMpi) {
        $executableName = "gp3d_mpi"
      } else {
        $executableName = "gp3d_sequential"
      }
      if ($env:OS -eq "Windows_NT") { $executableName += ".exe" }
    }
    $executablePath = Get-FullPath -BaseDirectory $buildDirectory -Path $executableName

    $programArguments = @()
    if (-not [string]::IsNullOrWhiteSpace([string]$config.run.input_file)) {
      $inputPath = Get-FullPath -BaseDirectory $projectRoot -Path ([string]$config.run.input_file)
      if (-not (Test-Path -LiteralPath $inputPath -PathType Leaf)) {
        throw "Solver input file was not found: $inputPath"
      }
      $programArguments += $inputPath
    }
    $programArguments += @($config.run.program_arguments)

    foreach ($variable in $config.run.environment.PSObject.Properties) {
      [Environment]::SetEnvironmentVariable($variable.Name, [string]$variable.Value, "Process")
      Write-Host "  env    : $($variable.Name)=$($variable.Value)"
    }

    if (-not $DryRun -and -not (Test-Path -LiteralPath $executablePath -PathType Leaf)) {
      throw "Solver executable was not found: $executablePath"
    }

    Push-Location $workingDirectory
    try {
      if ([bool]$config.run.use_mpi_launcher) {
        $launcherArguments = @([string]$config.run.process_option, [string]$config.run.processes)
        $launcherArguments += @($config.run.launcher_arguments)
        $launcherArguments += $executablePath
        $launcherArguments += $programArguments
        Invoke-CheckedCommand -Command ([string]$config.run.launcher) `
          -Arguments $launcherArguments -WhatIfOnly:$DryRun
      } else {
        Invoke-CheckedCommand -Command $executablePath `
          -Arguments $programArguments -WhatIfOnly:$DryRun
      }
    } finally {
      Pop-Location
    }
  }

  Write-Host "Workflow completed successfully."
} catch {
  Write-Error $_.Exception.Message
  exit 1
}
