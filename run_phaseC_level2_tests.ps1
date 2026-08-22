<#
.SYNOPSIS
    FloodA5 -- Phase C, Level 2 (diamond face-flux) test sequence.

.DESCRIPTION
    Runs the four-step test order from FloodA5_PhaseC_Level2_Summary.md:
      1. Existing unit/regression suite (face_flux_method untouched)
      2. Mesh-build smoke test -- diamond table construction only
      3. The actual hypothesis test -- planar-symmetry, diamond vs baseline
      4. Optional stretch tests -- point-spread, Carlisle regression

    Each step is a separate function so a failure in one step doesn't stop
    the others from running -- final summary lists pass/fail per step.
    Steps 2-4 are skipped automatically if their required mesh/output files
    are not found, with a clear message, rather than failing hard.

.PARAMETER RepoRoot
    Path to the FloodA5 repository root. Default: current directory.

.PARAMETER Threads
    Thread count for simulation runs (NOT mesh generation -- mesh
    generation must always use --threads 1 per Bugs 50/57; this script
    does not generate any mesh, so this applies to simulation runs only).
    Default: 'auto'.

.PARAMETER SkipRegression
    Skip step 1 (existing unit tests). Use if you've already confirmed
    these pass and just want to re-run the diamond-specific steps.

.PARAMETER RunStretchTests
    Run step 4 (point-spread + Carlisle regression). Off by default --
    these are longer runs and only meaningful once step 3 looks promising.

.PARAMETER DtMax
    dt_max (seconds) for the Step 3 planar-symmetry hypothesis run.
    Default: 2 -- matches the best dt-sensitivity result on record from
    the pentagon-chirality handover session.

.EXAMPLE
    .\run_phaseC_level2_tests.ps1

    Run steps 1-3 with defaults.

.EXAMPLE
    .\run_phaseC_level2_tests.ps1 -SkipRegression -RunStretchTests

    Skip the unit-test regression (already confirmed) and additionally
    run the point-spread and Carlisle stretch tests.
#>

[CmdletBinding()]
param(
    [string]$RepoRoot = ".",
    [string]$Threads = "auto",
    [switch]$SkipRegression,
    [switch]$RunStretchTests,
    [double]$DtMax = 2.0
)

$ErrorActionPreference = "Continue"   # keep going between steps; we track pass/fail ourselves
Set-Location $RepoRoot
$RepoRoot = (Get-Location).Path

$LogDir = Join-Path $RepoRoot "test\phaseC_level2_logs"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$Timestamp = Get-Date -Format "yyyyMMdd_HHmmss"

$Results = [System.Collections.Generic.List[object]]::new()

function Write-Section {
    param([string]$Title)
    Write-Host ""
    Write-Host ("=" * 78) -ForegroundColor Cyan
    Write-Host "  $Title" -ForegroundColor Cyan
    Write-Host ("=" * 78) -ForegroundColor Cyan
}

function Invoke-JuliaStep {
    <#
        Runs a julia command, tees output to a per-step log file, and
        records pass/fail (by exit code) into $Results.
    #>
    param(
        [string]$StepName,
        [string[]]$JuliaArgs,
        [string]$LogFileName
    )

    $logPath = Join-Path $LogDir "${Timestamp}_${LogFileName}"
    Write-Host "  -> julia $($JuliaArgs -join ' ')" -ForegroundColor DarkGray
    Write-Host "     log: $logPath" -ForegroundColor DarkGray

    & julia @JuliaArgs 2>&1 | Tee-Object -FilePath $logPath
    $exitCode = $LASTEXITCODE

    $status = if ($exitCode -eq 0) { "PASS" } else { "FAIL (exit $exitCode)" }
    $color  = if ($exitCode -eq 0) { "Green" } else { "Red" }
    Write-Host "  [$status] $StepName" -ForegroundColor $color

    $Results.Add([pscustomobject]@{
        Step   = $StepName
        Status = $status
        Log    = $logPath
    })
}

function Test-PathOrSkip {
    param([string]$Path, [string]$StepName)
    if (-not (Test-Path $Path)) {
        Write-Host "  [SKIP] $StepName -- required file not found: $Path" -ForegroundColor Yellow
        $Results.Add([pscustomobject]@{
            Step   = $StepName
            Status = "SKIP (missing: $Path)"
            Log    = ""
        })
        return $false
    }
    return $true
}

# =============================================================================
# Step 1 -- Existing unit / regression suite (face_flux_method untouched)
# =============================================================================
if (-not $SkipRegression) {
    Write-Section "Step 1: Regression suite (face_flux_method defaults everywhere)"

    $unitTests = @(
        "test\test_edge_geometry.jl",
        "test\test_noc_correction.jl",
        "test\test_sgs_unit.jl",
        "test\test_inflow_point.jl",
        "test\test_mirror_symmetry.jl",
        "test\test_cell_momentum.jl",
        "test\test_gradient_direction.jl"
    )

    foreach ($t in $unitTests) {
        $full = Join-Path $RepoRoot $t
        if (Test-PathOrSkip $full $t) {
            Invoke-JuliaStep -StepName $t `
                -JuliaArgs @("--project=.", $full) `
                -LogFileName ($t -replace '[\\\/]', '_')
        }
    }

    Write-Host ""
    Write-Host "  Step 1 done. All of these should match pre-patch pass counts" -ForegroundColor DarkGray
    Write-Host "  exactly -- none exercise face_flux_method=:diamond." -ForegroundColor DarkGray
}
else {
    Write-Host "Skipping Step 1 (regression suite) -- -SkipRegression set." -ForegroundColor Yellow
}

# =============================================================================
# Step 2 -- Mesh-build smoke test: diamond table construction only
# =============================================================================
Write-Section "Step 2: Diamond table construction smoke test (square res-18)"

$squareMesh = Join-Path $RepoRoot "test\square\square_mesh18_standard.parquet"
if (Test-PathOrSkip $squareMesh "Step 2: diamond table smoke test") {
    Invoke-JuliaStep -StepName "Step 2: diamond table smoke test" `
        -JuliaArgs @(
            "--project=.", "--threads", $Threads, "FloodModel.jl",
            "--meshload", $squareMesh,
            "--flow-model", "standard",
            "--gradient-correction", "on",
            "--face-flux-method", "diamond",
            "--sim-duration", "60",
            "--dt-max", "10"
        ) `
        -LogFileName "step2_diamond_smoke.log"

    Write-Host ""
    Write-Host "  Check the log above for the two new @info lines:" -ForegroundColor DarkGray
    Write-Host "    'Diamond flux: vertex table built -- N vertices (...)'" -ForegroundColor DarkGray
    Write-Host "    'Diamond flux: M/N_edges edges have a valid diamond record (...)'" -ForegroundColor DarkGray
    Write-Host "  Compare M and the fallback breakdown against your own" -ForegroundColor DarkGray
    Write-Host "  audit_diamond_gradient.jl run on this same mesh (square res-18:" -ForegroundColor DarkGray
    Write-Host "  227/11,523 k<3 vertex fallback, 0 degenerate diamonds)." -ForegroundColor DarkGray
}

# =============================================================================
# Step 3 -- The actual hypothesis test: planar-symmetry, diamond vs baseline
# =============================================================================
Write-Section "Step 3: Planar-symmetry hypothesis test (diamond, dt_max=$DtMax)"

$planarMesh    = Join-Path $RepoRoot "test\planar_embankment\planar_mesh18_std.parquet"
$baselineH5    = Join-Path $RepoRoot "test\planar_embankment\planar_res18_std_baseline.h5"
$diamondOutDir = Join-Path $RepoRoot "test\planar_embankment"
$diamondH5     = Join-Path $diamondOutDir "planar_res18_diamond_dtmax$($DtMax).h5"
$symmetryTest  = Join-Path $RepoRoot "test\test_planar_symmetry.jl"

if ((Test-PathOrSkip $planarMesh "Step 3: diamond run (planar res-18)") -and
    (Test-PathOrSkip $symmetryTest "Step 3: diamond run (planar res-18)")) {

    Invoke-JuliaStep -StepName "Step 3a: diamond simulation run" `
        -JuliaArgs @(
            "--project=.", "--threads", $Threads, "FloodModel.jl",
            "--meshload", $planarMesh,
            "--flow-model", "standard",
            "--injection-point", "51.0001,-0.0434,0.1",
            "--closed-boundaries",
            "--sim-duration", "72000",
            "--dt-max", "$DtMax",
            "--gradient-correction", "on",
            "--face-flux-method", "diamond",
            "--output", $diamondH5,
            "--output-interval", "300"
        ) `
        -LogFileName "step3a_diamond_run.log"

    if (Test-Path $baselineH5) {
        Invoke-JuliaStep -StepName "Step 3b: planar-symmetry comparison" `
            -JuliaArgs @(
                "--project=.", $symmetryTest,
                "--source-lat", "51.0001",
                "--source-lon", "-0.0434",
                "--baseline", $baselineH5,
                "--corrected", $diamondH5,
                "--sweep", "10"
            ) `
            -LogFileName "step3b_symmetry_compare.log"

        Write-Host ""
        Write-Host "  Compare the reported |asymmetry| against these prior results:" -ForegroundColor DarkGray
        Write-Host "    Baseline (no correction):            ~0.87-0.92, never flips" -ForegroundColor DarkGray
        Write-Host "    alpha=0 (orthogonal-only, legacy):    ~0.83-0.85, never flips" -ForegroundColor DarkGray
        Write-Host "    alpha=1 + cell-momentum (legacy):     flips to ~-0.58 to -0.94" -ForegroundColor DarkGray
        Write-Host "    Diamond (this run):                   <-- the number that matters" -ForegroundColor DarkGray
    }
    else {
        Write-Host "  [SKIP] Step 3b comparison -- baseline file not found: $baselineH5" -ForegroundColor Yellow
        Write-Host "         (Step 3a's diamond run still completed -- you can compare" -ForegroundColor Yellow
        Write-Host "          it manually once you have a baseline .h5 to point at.)" -ForegroundColor Yellow
        $Results.Add([pscustomobject]@{
            Step   = "Step 3b: planar-symmetry comparison"
            Status = "SKIP (missing baseline: $baselineH5)"
            Log    = ""
        })
    }
}

# =============================================================================
# Step 4 -- Optional stretch tests
# =============================================================================
if ($RunStretchTests) {
    Write-Section "Step 4: Stretch tests (point-spread, Carlisle regression)"

    # --- Point-spread, square mesh, diamond vs baseline -----------------
    $pointSpreadTest = Join-Path $RepoRoot "test\test_point_spread.jl"
    $squareBaselineH5 = Join-Path $RepoRoot "test\square\square_res18_baseline.h5"
    $squareDiamondH5  = Join-Path $RepoRoot "test\square\square_res18_diamond.h5"

    if ((Test-PathOrSkip $squareMesh "Step 4a: point-spread diamond run") -and
        (Test-PathOrSkip $pointSpreadTest "Step 4a: point-spread diamond run")) {

        Invoke-JuliaStep -StepName "Step 4a: point-spread diamond run" `
            -JuliaArgs @(
                "--project=.", "--threads", $Threads, "FloodModel.jl",
                "--meshload", $squareMesh,
                "--flow-model", "standard",
                "--injection-point", "0,0,50",
                "--closed-boundaries",
                "--sim-duration", "3600",
                "--dt-max", "10",
                "--gradient-correction", "on",
                "--face-flux-method", "diamond",
                "--output", $squareDiamondH5,
                "--output-interval", "60"
            ) `
            -LogFileName "step4a_pointspread_run.log"

        if (Test-Path $squareBaselineH5) {
            Invoke-JuliaStep -StepName "Step 4a: point-spread comparison" `
                -JuliaArgs @(
                    "--project=.", $pointSpreadTest,
                    "--baseline", $squareBaselineH5,
                    "--corrected", $squareDiamondH5,
                    "--frame", "50"
                ) `
                -LogFileName "step4a_pointspread_compare.log"
        }
        else {
            Write-Host "  [SKIP] point-spread comparison -- no baseline at $squareBaselineH5" -ForegroundColor Yellow
        }
    }

    # --- Carlisle regression, standard flow, diamond ---------------------
    $carlisleMesh = Join-Path $RepoRoot "test\carlisle\carlisle_mesh18_standard.parquet"
    $carlisleOut  = Join-Path $RepoRoot "test\carlisle\carlisle_diamond_res18.h5"

    if (Test-PathOrSkip $carlisleMesh "Step 4b: Carlisle regression (diamond)") {
        Invoke-JuliaStep -StepName "Step 4b: Carlisle regression (diamond)" `
            -JuliaArgs @(
                "--project=.", "--threads", $Threads, "FloodModel.jl",
                "--meshload", $carlisleMesh,
                "--flow-model", "standard",
                "--rainfall", "50",
                "--sim-duration", "36000",
                "--gradient-correction", "on",
                "--face-flux-method", "diamond",
                "--output", $carlisleOut,
                "--output-interval", "600"
            ) `
            -LogFileName "step4b_carlisle_diamond.log"

        Write-Host ""
        Write-Host "  Check the log for: mb_err = 0.0 throughout, dt stable" -ForegroundColor DarkGray
        Write-Host "  (declining or oscillating), no unexpected NaN_cells growth." -ForegroundColor DarkGray
    }
}
else {
    Write-Host ""
    Write-Host "Skipping Step 4 (stretch tests) -- pass -RunStretchTests to include." -ForegroundColor Yellow
}

# =============================================================================
# Summary
# =============================================================================
Write-Section "Summary"

$Results | Format-Table -AutoSize -Property Step, Status

$failCount = ($Results | Where-Object { $_.Status -like "FAIL*" }).Count
$skipCount = ($Results | Where-Object { $_.Status -like "SKIP*" }).Count
$passCount = ($Results | Where-Object { $_.Status -eq "PASS" }).Count

Write-Host ""
Write-Host "  $passCount passed, $failCount failed, $skipCount skipped." -ForegroundColor $(if ($failCount -gt 0) { "Red" } else { "Green" })
Write-Host "  Logs: $LogDir" -ForegroundColor DarkGray

if ($failCount -gt 0) {
    Write-Host ""
    Write-Host "  One or more steps failed -- do not proceed to Step 4 / production" -ForegroundColor Red
    Write-Host "  use of --face-flux-method diamond until these are resolved." -ForegroundColor Red
    exit 1
}
exit 0
