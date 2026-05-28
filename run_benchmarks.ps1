# run_benchmarks.ps1 — FloodA5 thread-scaling benchmark sweep
#
# Runs benchmark_sim.jl once per thread count, appending all results to a
# single CSV so the file can be imported directly into Excel / Python for
# scaling analysis.
#
# Usage (from project root in PowerShell):
#   .\run_benchmarks.ps1
#   .\run_benchmarks.ps1 -Out test/my_results.csv -Warmup 5 -Steps 100
#   .\run_benchmarks.ps1 -MaxThreads 8          # override auto-detect
#
# The script detects the logical CPU count automatically.  Pass -MaxThreads N
# to override (useful when hyperthreading inflates the count beyond useful
# parallelism, or when you want a shorter run).

param(
    [string]  $Out        = "test/benchmark_results.csv",
    [int]     $Warmup     = 10,
    [int]     $Steps      = 200,
    [float]   $Dt         = 30.0,
    [string]  $Config     = "",           # path to benchmark_config.json; empty = use default
    [int]     $MaxThreads = 0,            # 0 = auto-detect logical CPU count
    [switch]  $DeleteExisting             # delete the CSV before starting (fresh run)
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ── Locate julia executable ───────────────────────────────────────────────────
$Julia = (Get-Command julia -ErrorAction SilentlyContinue)?.Source
if (-not $Julia) {
    Write-Error "julia not found on PATH. Add Julia\bin to your PATH and retry."
    exit 1
}

# ── Resolve thread count ──────────────────────────────────────────────────────
if ($MaxThreads -le 0) {
    $MaxThreads = [System.Environment]::ProcessorCount
    Write-Host "Detected $MaxThreads logical CPU cores."
}

# Build thread list: 1, 2, 4, 8, … up to MaxThreads, always including MaxThreads
$ThreadCounts = @(1)
$t = 2
while ($t -lt $MaxThreads) {
    $ThreadCounts += $t
    $t *= 2
}
if ($ThreadCounts[-1] -ne $MaxThreads) {
    $ThreadCounts += $MaxThreads
}

# ── Optionally delete existing results file ───────────────────────────────────
if ($DeleteExisting -and (Test-Path $Out)) {
    Remove-Item $Out
    Write-Host "Deleted existing $Out"
}

# ── Build common argument list ────────────────────────────────────────────────
$CommonArgs = @(
    "--out",    $Out,
    "--warmup", $Warmup,
    "--steps",  $Steps,
    "--dt",     $Dt
)
if ($Config -ne "") {
    $CommonArgs += @("--config", $Config)
}

# ── Run sweep ─────────────────────────────────────────────────────────────────
$TotalRuns = $ThreadCounts.Count
$RunIndex  = 0

Write-Host ""
Write-Host ("=" * 72)
Write-Host "FloodA5 benchmark sweep  —  $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Write-Host "  Thread counts : $($ThreadCounts -join ', ')"
Write-Host "  Output CSV    : $Out"
Write-Host "  Warmup steps  : $Warmup"
Write-Host "  Timed steps   : $Steps"
Write-Host "  Fixed dt      : ${Dt}s"
Write-Host ("=" * 72)

foreach ($N in $ThreadCounts) {
    $RunIndex++
    Write-Host ""
    Write-Host "── Run $RunIndex / $TotalRuns  (--threads $N) ──"

    $StartTime = Get-Date
    & $Julia --threads $N benchmark_sim.jl @CommonArgs
    $ExitCode = $LASTEXITCODE
    $Elapsed  = ((Get-Date) - $StartTime).TotalSeconds

    if ($ExitCode -ne 0) {
        Write-Warning "julia exited with code $ExitCode for --threads $N"
    } else {
        Write-Host ("  Completed in {0:F1} s" -f $Elapsed)
    }
}

Write-Host ""
Write-Host ("=" * 72)
Write-Host "Sweep complete.  Results appended to: $Out"
Write-Host ("=" * 72)
Write-Host ""
