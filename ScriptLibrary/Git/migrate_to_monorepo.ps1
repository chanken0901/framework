#Requires -Version 5.1

<#
.SYNOPSIS
Safely imports ScriptLibrary and SolverLibrary into framework.git as a monorepo.

.DESCRIPTION
The script performs a read-only preflight by default.  With -Apply it:

1. Requires clean, synchronized main branches in all three source repositories.
2. Optionally creates and pushes annotated safety tags.
3. Clones framework.git into a new sibling directory.
4. Imports ScriptLibrary and SolverLibrary with history-preserving native Git commands.
5. Adds monorepo ignore rules and a migration report.
6. Optionally runs smoke tests and pushes the migration branch.

It never deletes directories, resets commits, rebases, force-pushes, or changes
the existing FrameWork checkout in place.

.EXAMPLE
.\migrate_to_monorepo.ps1

Runs preflight checks and prints the planned commands without changing anything.

.EXAMPLE
.\migrate_to_monorepo.ps1 -Apply -CreateSafetyTags -PushSafetyTags -RunSmokeTests

Creates the local migration checkout and safety tags, but does not push the
migration branch.

.EXAMPLE
.\migrate_to_monorepo.ps1 -Apply -CreateSafetyTags -PushSafetyTags `
  -RunSmokeTests -PushBranch

Also pushes chore/monorepo-migration to framework.git.  No force push is used.
#>

[CmdletBinding()]
param(
    [string]$FrameworkRoot,
    [string]$DestinationRoot,
    [string]$FrameworkUrl = "https://github.com/chanken0901/framework.git",
    [string]$ScriptLibraryUrl = "https://github.com/chanken0901/ScriptLibrary.git",
    [string]$SolverLibraryUrl = "https://github.com/chanken0901/SolverLibrary.git",
    [string]$BaseBranch = "main",
    [string]$MigrationBranch = "chore/monorepo-migration",
    [string]$SafetyTag = ("pre-monorepo-{0}" -f (Get-Date -Format "yyyy-MM-dd")),
    [string]$GitCommand = "git",
    [string]$PythonCommand = "python",
    [switch]$CreateSafetyTags,
    [switch]$PushSafetyTags,
    [switch]$RunSmokeTests,
    [switch]$PushBranch,
    [switch]$Apply
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"


function Write-Section {
    param([Parameter(Mandatory = $true)][string]$Message)
    Write-Host ""
    Write-Host "=== $Message ===" -ForegroundColor Cyan
}


function Format-Command {
    param(
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )

    $formatted = foreach ($argument in $Arguments) {
        if ($argument -match '[\s`"]') {
            '"{0}"' -f ($argument -replace '"', '\"')
        }
        else {
            $argument
        }
    }
    return "$Executable $($formatted -join ' ')"
}


function Invoke-Native {
    param(
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string[]]$Arguments,
        [switch]$AllowFailure,
        [switch]$ShowOutput
    )

    Write-Host ("[CMD] " + (Format-Command -Executable $Executable -Arguments $Arguments))
    # Windows PowerShell 5.1 converts native stderr records into non-terminating
    # PowerShell errors.  Git legitimately writes clone/fetch progress to stderr,
    # so temporarily allow those records and judge success only by the exit code.
    $savedErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        $output = @(& $Executable @Arguments 2>&1)
        $exitCode = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }

    if ($ShowOutput -and $output.Count -gt 0) {
        $output | ForEach-Object { Write-Host $_ }
    }
    if ($exitCode -ne 0 -and -not $AllowFailure) {
        $detail = ($output | Out-String).Trim()
        if (-not $detail) {
            $detail = "The command returned exit code $exitCode."
        }
        throw "Command failed:`n$detail"
    }

    return [pscustomobject]@{
        ExitCode = $exitCode
        Output = @($output | ForEach-Object { "$_" })
    }
}


function Invoke-Git {
    param(
        [string]$RepositoryPath,
        [Parameter(Mandatory = $true)][string[]]$Arguments,
        [switch]$AllowFailure,
        [switch]$ShowOutput
    )

    $gitArguments = @("-c", "core.quotepath=false")
    if ($RepositoryPath) {
        $gitArguments += @(
            "-c", "safe.directory=$RepositoryPath",
            "-C", $RepositoryPath
        )
    }
    $gitArguments += $Arguments
    return Invoke-Native -Executable $GitCommand -Arguments $gitArguments `
        -AllowFailure:$AllowFailure -ShowOutput:$ShowOutput
}


function Get-SingleLine {
    param([Parameter(Mandatory = $true)]$Result)
    return (($Result.Output -join "`n").Trim())
}


function Test-PathIsInside {
    param(
        [Parameter(Mandatory = $true)][string]$Candidate,
        [Parameter(Mandatory = $true)][string]$Parent
    )

    $separator = [IO.Path]::DirectorySeparatorChar
    $candidateFull = [IO.Path]::GetFullPath($Candidate).TrimEnd($separator)
    $parentFull = [IO.Path]::GetFullPath($Parent).TrimEnd($separator)
    return $candidateFull.Equals($parentFull, [StringComparison]::OrdinalIgnoreCase) -or
        $candidateFull.StartsWith($parentFull + $separator, [StringComparison]::OrdinalIgnoreCase)
}


function Show-PlannedCommands {
    param(
        [Parameter(Mandatory = $true)][string]$Destination,
        [Parameter(Mandatory = $true)][object[]]$Repositories
    )

    Write-Section "Planned mutating commands (not executed)"
    if ($CreateSafetyTags) {
        foreach ($repo in $Repositories) {
            Write-Host "git -C `"$($repo.Path)`" tag -a $SafetyTag -m `"State before monorepo migration`""
            if ($PushSafetyTags) {
                Write-Host "git -C `"$($repo.Path)`" push origin refs/tags/$SafetyTag"
            }
        }
    }
    Write-Host "git clone $FrameworkUrl `"$Destination`""
    Write-Host "git -C `"$Destination`" switch -c $MigrationBranch"
    Write-Host "git -C `"$Destination`" remote add scriptlibrary-import $ScriptLibraryUrl"
    Write-Host "git -C `"$Destination`" fetch --no-tags scriptlibrary-import $BaseBranch"
    Write-Host "git -C `"$Destination`" merge -s ours --no-commit --allow-unrelated-histories scriptlibrary-import/$BaseBranch"
    Write-Host "git -C `"$Destination`" read-tree --prefix=ScriptLibrary/ -u scriptlibrary-import/$BaseBranch"
    Write-Host "git -C `"$Destination`" commit -m `"Import ScriptLibrary history into monorepo`""
    Write-Host "git -C `"$Destination`" remote add solverlibrary-import $SolverLibraryUrl"
    Write-Host "git -C `"$Destination`" fetch --no-tags solverlibrary-import $BaseBranch"
    Write-Host "git -C `"$Destination`" merge -s ours --no-commit --allow-unrelated-histories solverlibrary-import/$BaseBranch"
    Write-Host "git -C `"$Destination`" read-tree --prefix=SolverLibrary/ -u solverlibrary-import/$BaseBranch"
    Write-Host "git -C `"$Destination`" commit -m `"Import SolverLibrary history into monorepo`""
    if ($PushBranch) {
        Write-Host "git -C `"$Destination`" push -u origin $MigrationBranch"
    }
}


if ([string]::IsNullOrWhiteSpace($FrameworkRoot)) {
    $FrameworkRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
}
else {
    $FrameworkRoot = [IO.Path]::GetFullPath($FrameworkRoot)
}

if ([string]::IsNullOrWhiteSpace($DestinationRoot)) {
    $frameworkParent = Split-Path -Parent $FrameworkRoot
    $frameworkName = Split-Path -Leaf $FrameworkRoot
    $DestinationRoot = Join-Path $frameworkParent ($frameworkName + "-monorepo-migration")
}
$DestinationRoot = [IO.Path]::GetFullPath($DestinationRoot)

if (Test-PathIsInside -Candidate $DestinationRoot -Parent $FrameworkRoot) {
    throw "DestinationRoot must be outside FrameworkRoot: $DestinationRoot"
}
if ($PushSafetyTags -and -not $CreateSafetyTags) {
    throw "-PushSafetyTags requires -CreateSafetyTags."
}
if ($MigrationBranch -notmatch '^[A-Za-z0-9._/-]+$') {
    throw "MigrationBranch contains unsupported characters: $MigrationBranch"
}
if ($SafetyTag -notmatch '^[A-Za-z0-9._/-]+$') {
    throw "SafetyTag contains unsupported characters: $SafetyTag"
}

$repositories = @(
    [pscustomobject]@{
        Name = "FrameWork"
        Path = $FrameworkRoot
        ExpectedUrl = $FrameworkUrl
        AllowedUntracked = @("?? ScriptLibrary/", "?? SolverLibrary/", "?? ScriptLibrary", "?? SolverLibrary")
    },
    [pscustomobject]@{
        Name = "ScriptLibrary"
        Path = (Join-Path $FrameworkRoot "ScriptLibrary")
        ExpectedUrl = $ScriptLibraryUrl
        AllowedUntracked = @()
    },
    [pscustomobject]@{
        Name = "SolverLibrary"
        Path = (Join-Path $FrameworkRoot "SolverLibrary")
        ExpectedUrl = $SolverLibraryUrl
        AllowedUntracked = @()
    }
)

Write-Section "Monorepo migration preflight"
Write-Host "Source      : $FrameworkRoot"
Write-Host "Destination : $DestinationRoot"
Write-Host "Branch      : $MigrationBranch"
Write-Host "Mode        : $(if ($Apply) { 'APPLY' } else { 'DRY-RUN' })"

$gitProbe = Invoke-Native -Executable $GitCommand -Arguments @("--version") -ShowOutput
if ($gitProbe.ExitCode -ne 0) {
    throw "Git is not available: $GitCommand"
}

$blockers = New-Object System.Collections.Generic.List[string]
$sourceCommits = @{}

foreach ($repo in $repositories) {
    Write-Section ("Inspect " + $repo.Name)
    if (-not (Test-Path -LiteralPath $repo.Path -PathType Container)) {
        $blockers.Add("$($repo.Name): directory does not exist: $($repo.Path)")
        continue
    }
    if (-not (Test-Path -LiteralPath (Join-Path $repo.Path ".git"))) {
        $blockers.Add("$($repo.Name): .git was not found: $($repo.Path)")
        continue
    }

    $topResult = Invoke-Git -RepositoryPath $repo.Path `
        -Arguments @("rev-parse", "--show-toplevel") -AllowFailure
    if ($topResult.ExitCode -ne 0) {
        $blockers.Add("$($repo.Name): Git cannot read the repository.")
        continue
    }

    $actualTop = [IO.Path]::GetFullPath((Get-SingleLine $topResult))
    if (-not $actualTop.Equals([IO.Path]::GetFullPath($repo.Path), [StringComparison]::OrdinalIgnoreCase)) {
        $blockers.Add("$($repo.Name): unexpected repository root: $actualTop")
    }

    $branch = Get-SingleLine (Invoke-Git -RepositoryPath $repo.Path `
        -Arguments @("branch", "--show-current"))
    Write-Host "branch: $branch"
    if ($branch -ne $BaseBranch) {
        $blockers.Add("$($repo.Name): switch to $BaseBranch before migration (current: $branch).")
    }

    $originResult = Invoke-Git -RepositoryPath $repo.Path `
        -Arguments @("remote", "get-url", "origin") -AllowFailure
    if ($originResult.ExitCode -ne 0) {
        $blockers.Add("$($repo.Name): origin is not configured.")
    }
    else {
        $originUrl = Get-SingleLine $originResult
        Write-Host "origin: $originUrl"
        if ($originUrl -ne $repo.ExpectedUrl) {
            $blockers.Add("$($repo.Name): origin differs from the expected URL: $($repo.ExpectedUrl)")
        }
    }

    $statusResult = Invoke-Git -RepositoryPath $repo.Path `
        -Arguments @("status", "--porcelain=v1", "--untracked-files=normal")
    $unexpectedStatus = @(
        $statusResult.Output | Where-Object {
            $_ -and ($repo.AllowedUntracked -notcontains $_)
        }
    )
    if ($unexpectedStatus.Count -gt 0) {
        Write-Host "Uncommitted entries:" -ForegroundColor Yellow
        $unexpectedStatus | ForEach-Object { Write-Host "  $_" -ForegroundColor Yellow }
        $blockers.Add("$($repo.Name): commit or intentionally remove all displayed changes.")
    }

    if ($Apply) {
        Invoke-Git -RepositoryPath $repo.Path `
            -Arguments @("fetch", "--prune", "origin", $BaseBranch) -ShowOutput | Out-Null
    }
    else {
        Write-Host "[DRY-RUN] git -C `"$($repo.Path)`" fetch --prune origin $BaseBranch"
    }

    $headResult = Invoke-Git -RepositoryPath $repo.Path `
        -Arguments @("rev-parse", "HEAD") -AllowFailure
    $remoteResult = Invoke-Git -RepositoryPath $repo.Path `
        -Arguments @("rev-parse", "origin/$BaseBranch") -AllowFailure
    if ($headResult.ExitCode -eq 0) {
        $headCommit = Get-SingleLine $headResult
        $sourceCommits[$repo.Name] = $headCommit
        Write-Host "HEAD: $headCommit"
    }
    if ($headResult.ExitCode -ne 0 -or $remoteResult.ExitCode -ne 0) {
        $blockers.Add("$($repo.Name): HEAD or origin/$BaseBranch cannot be resolved.")
    }
    elseif ((Get-SingleLine $headResult) -ne (Get-SingleLine $remoteResult)) {
        $blockers.Add("$($repo.Name): HEAD and origin/$BaseBranch differ; pull or push first.")
    }
}

if (Test-Path -LiteralPath $DestinationRoot) {
    $blockers.Add("Destination already exists; choose a new -DestinationRoot: $DestinationRoot")
}

if ($blockers.Count -gt 0) {
    Write-Section "Preflight blockers"
    $blockers | ForEach-Object { Write-Host "- $_" -ForegroundColor Yellow }
}
else {
    Write-Host "Preflight passed." -ForegroundColor Green
}

if (-not $Apply) {
    Show-PlannedCommands -Destination $DestinationRoot -Repositories $repositories
    Write-Host ""
    Write-Host "No files, tags, branches, or remotes were changed."
    if ($blockers.Count -gt 0) {
        Write-Host "Resolve the blockers, run this DryRun again, then add -Apply." -ForegroundColor Yellow
    }
    exit 0
}

if ($blockers.Count -gt 0) {
    throw "Migration was not started because preflight found $($blockers.Count) blocker(s)."
}

if ($CreateSafetyTags) {
    Write-Section "Create safety tags"
    foreach ($repo in $repositories) {
        $tagRef = "refs/tags/$SafetyTag"
        $existingTag = Invoke-Git -RepositoryPath $repo.Path `
            -Arguments @("rev-parse", "-q", "--verify", $tagRef) -AllowFailure
        if ($existingTag.ExitCode -eq 0) {
            $tagCommitResult = Invoke-Git -RepositoryPath $repo.Path `
                -Arguments @("rev-list", "-n", "1", $tagRef)
            if ((Get-SingleLine $tagCommitResult) -ne $sourceCommits[$repo.Name]) {
                throw "$($repo.Name): existing tag $SafetyTag points to a different commit."
            }
            Write-Host "$($repo.Name): existing safety tag is valid."
        }
        else {
            Invoke-Git -RepositoryPath $repo.Path `
                -Arguments @("tag", "-a", $SafetyTag, "-m", "State before monorepo migration") | Out-Null
        }
        if ($PushSafetyTags) {
            Invoke-Git -RepositoryPath $repo.Path `
                -Arguments @("push", "origin", $tagRef) -ShowOutput | Out-Null
        }
    }
}

Write-Section "Clone framework and create migration branch"
Invoke-Git -Arguments @("clone", $FrameworkUrl, $DestinationRoot) -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("switch", "-c", $MigrationBranch) -ShowOutput | Out-Null

foreach ($prefix in @("ScriptLibrary", "SolverLibrary")) {
    if (Test-Path -LiteralPath (Join-Path $DestinationRoot $prefix)) {
        throw "The fresh framework clone already contains $prefix. Nothing was deleted; inspect $DestinationRoot."
    }
}

Write-Section "Import ScriptLibrary"
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("remote", "add", "scriptlibrary-import", $ScriptLibraryUrl) | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("fetch", "--no-tags", "scriptlibrary-import", $BaseBranch) -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @(
        "merge", "-s", "ours", "--no-commit", "--allow-unrelated-histories",
        "scriptlibrary-import/$BaseBranch"
    ) -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("read-tree", "--prefix=ScriptLibrary/", "-u", "scriptlibrary-import/$BaseBranch") `
    -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("commit", "-m", "Import ScriptLibrary history into monorepo") `
    -ShowOutput | Out-Null

Write-Section "Import SolverLibrary"
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("remote", "add", "solverlibrary-import", $SolverLibraryUrl) | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("fetch", "--no-tags", "solverlibrary-import", $BaseBranch) -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @(
        "merge", "-s", "ours", "--no-commit", "--allow-unrelated-histories",
        "solverlibrary-import/$BaseBranch"
    ) -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("read-tree", "--prefix=SolverLibrary/", "-u", "solverlibrary-import/$BaseBranch") `
    -ShowOutput | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("commit", "-m", "Import SolverLibrary history into monorepo") `
    -ShowOutput | Out-Null

Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("remote", "remove", "scriptlibrary-import") | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("remote", "remove", "solverlibrary-import") | Out-Null

Write-Section "Add monorepo metadata"
$ignorePath = Join-Path $DestinationRoot ".gitignore"
$ignoreMarker = "# BEGIN MONOREPO GENERATED ARTIFACTS"
$ignoreBlock = @"

$ignoreMarker
__pycache__/
*.py[cod]
.pytest_cache/
build/
cmake-build-*/
CMakeFiles/
CMakeCache.txt
CTestTestfile.cmake
*.o
*.obj
*.mod
*.smod
*.exe
*.dll
*.so
*.a
*.lib
output/
*.slf
*.vti
*.pvd
fort.*
*.log
*.tmp
*.bak
~`$*
# END MONOREPO GENERATED ARTIFACTS
"@

$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
$existingIgnore = if (Test-Path -LiteralPath $ignorePath) {
    [IO.File]::ReadAllText($ignorePath)
}
else {
    ""
}
if (-not $existingIgnore.Contains($ignoreMarker)) {
    [IO.File]::WriteAllText($ignorePath, $existingIgnore.TrimEnd() + $ignoreBlock, $utf8NoBom)
}

$reportDirectory = Join-Path $DestinationRoot "docs\development"
$reportPath = Join-Path $reportDirectory "MONOREPO_MIGRATION.md"
if (Test-Path -LiteralPath $reportPath) {
    throw "Migration report already exists; nothing was overwritten: $reportPath"
}
New-Item -ItemType Directory -Path $reportDirectory -Force | Out-Null
$report = @"
# Monorepo migration record

- Migration date: $(Get-Date -Format "yyyy-MM-dd")
- Base repository: $FrameworkUrl
- Base branch: $BaseBranch
- Migration branch: $MigrationBranch
- FrameWork source commit: $($sourceCommits['FrameWork'])
- ScriptLibrary source commit: $($sourceCommits['ScriptLibrary'])
- SolverLibrary source commit: $($sourceCommits['SolverLibrary'])
- Safety tag: $(if ($CreateSafetyTags) { $SafetyTag } else { 'not created by script' })

`ScriptLibrary` and `SolverLibrary` were imported as merge parents with their
complete commit graphs, then placed below their prefixes with `git read-tree`.
The former standalone repositories must remain available until the monorepo
has passed review and a release tag has been created.
"@
[IO.File]::WriteAllText($reportPath, $report.TrimEnd() + "`n", $utf8NoBom)

Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("add", "--", ".gitignore", "docs/development/MONOREPO_MIGRATION.md") | Out-Null
Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("commit", "-m", "Document monorepo migration and ignore generated files") `
    -ShowOutput | Out-Null

Write-Section "Verify imported tree"
$requiredFiles = @(
    "ScriptLibrary/RunEnvironment/prepare_environment.py",
    "ScriptLibrary/Git/manage_library_repositories.py",
    "SolverLibrary/GPE/gp3d/CMakeLists.txt",
    "SolverLibrary/GPE/gp3d/solver_manifest.yaml",
    "SolverLibrary/NSE/CMakeLists.txt"
)
foreach ($requiredFile in $requiredFiles) {
    Invoke-Git -RepositoryPath $DestinationRoot `
        -Arguments @("ls-files", "--error-unmatch", "--", $requiredFile) | Out-Null
    Write-Host "tracked: $requiredFile"
}

foreach ($prefix in @("ScriptLibrary", "SolverLibrary")) {
    if (Test-Path -LiteralPath (Join-Path $DestinationRoot "$prefix\.git")) {
        throw "Nested Git metadata remains under $prefix."
    }
}

Invoke-Git -RepositoryPath $DestinationRoot -Arguments @("diff", "--check") | Out-Null
$finalStatus = Invoke-Git -RepositoryPath $DestinationRoot `
    -Arguments @("status", "--porcelain=v1", "--untracked-files=normal")
if ($finalStatus.Output.Count -gt 0 -and ($finalStatus.Output -join "").Trim()) {
    throw "The migration checkout is not clean:`n$($finalStatus.Output -join "`n")"
}

if ($RunSmokeTests) {
    Write-Section "Run smoke tests"
    Invoke-Native -Executable $PythonCommand -Arguments @(
        "-m", "unittest", "discover",
        "-s", (Join-Path $DestinationRoot "ScriptLibrary\RunEnvironment\tests"),
        "-p", "test_*.py"
    ) -ShowOutput | Out-Null
    Invoke-Native -Executable $PythonCommand -Arguments @(
        (Join-Path $DestinationRoot "ScriptLibrary\RunEnvironment\prepare_environment.py"),
        (Join-Path $DestinationRoot "ScriptLibrary\RunEnvironment\environment.gpe.yaml"),
        "--framework-root", $DestinationRoot,
        "--dry-run"
    ) -ShowOutput | Out-Null
}

if ($PushBranch) {
    Write-Section "Push migration branch"
    Invoke-Git -RepositoryPath $DestinationRoot `
        -Arguments @("push", "-u", "origin", $MigrationBranch) -ShowOutput | Out-Null
}

Write-Section "Migration completed"
Write-Host "Local monorepo : $DestinationRoot" -ForegroundColor Green
Write-Host "Branch         : $MigrationBranch"
if ($PushBranch) {
    Write-Host "The branch was pushed without force. Open a pull request into $BaseBranch."
}
else {
    Write-Host "Review the checkout, then push it with:"
    Write-Host "git -C `"$DestinationRoot`" push -u origin $MigrationBranch"
}
Write-Host "Do not archive the standalone repositories until the pull request is merged and tagged."
