<#
.SYNOPSIS
    FloodA5 -- diamond face-flux method, dt-sensitivity sweep.

.DESCRIPTION
    Runs the planar-embankment hypothesis test at several --dt-max values
    with --face-flux-method diamond, compares each against the same
    baseline .h5 via test_planar_symmetry.jl, and parses the sweep table
    each run produces down to a single consolidated summary table:
    first-checkpoint asymmetry, first sign-flip frame (if any), and final
    asymmetry, side by side across dt values.

    This directly mirrors the {10, 4, 2, 1, 0.5} sweep already done for
    the legacy alpha=1 + cell-momentum configuration
    (FloodA5_GradientCorrection_PentagonChirality_Handover.md §4) --
    same mesh, same injection point, same comparison test -- so the two
    tables are directly comparable line for line.

    V-shaped response (error smallest at some intermediate dt, worse at
    both larger AND smaller dt) => explicit-stability-limit story.
    Flat or monotonically-worsening-as-dt-shrinks response => rules out
    ordinary stability error, points at a structural/geometric effect
    instead (or a bug -- see the companion diagnose_diamond_vs_legacy_
    gradient.jl script for a same-snapshot check that doesn't depend on
    running any new simulation).

.PARAMETER RepoRoot
    Path to the FloodA5 repository root. Default: current directory.

.PARAMETER Threads
    Thread count for simulation runs. Default: 'auto'.

.PARAMETER DtValues
    List of --dt-max values (seconds) to sweep. Default matches the
    legacy alpha=1 sweep exactly: 10, 4, 2, 1, 0.5.

.PARAMETER SimDuration
    Simulation duration in seconds. Default: 72000 (20h), matching all
    prior planar-symmetry runs on record.

.EXAMPLE
    .\run_diamond_dt_sweep.ps1

.EXAMPLE
    .\run_diamond_dt_sweep.ps1 -DtValues 2,1,0.5,0.25 -SimDuration 36000
#>

[CmdletBinding()]
param(
    [string]$RepoRoot = ".",
    [string]$Threads = "auto",
    [double[]]$DtValues = @(10, 4, 2, 1, 0.5),
    [double]$SimDuration = 72000
)

Set-Location $RepoRoot
$RepoRoot = (Get-Location).Path

$LogDir = Join-Path $RepoRoot "test\phaseC_dt_sweep_logs"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$Timestamp = Get-Date -Format "yyyyMMdd_HHmmss"

$PlanarMesh   = Join-Path $RepoRoot "test\planar_embankment\planar_mesh18_std.parquet"
$BaselineH5   = Join-Path $RepoRoot "test\planar_embankment\planar_res18_std_baseline.h5"
$SymmetryTest = Join-Path $RepoRoot "test\test_planar_symmetry.jl"
$InjectionPoint = "51.0001,-0.0434,0.1"
$InjectionLat = "51.0001"
$InjectionLon = "-0.0434"

foreach ($p in @($PlanarMesh, $BaselineH5, $SymmetryTest)) {
    if (-not (Test-Path $p)) {
        Write-Host "ERROR: required file not found: $p" -ForegroundColor Red
        Write-Host "Edit the path variables near the top of this script if your" -ForegroundColor Red
        Write-Host "layout differs, then re-run." -ForegroundColor Red
        exit 1
    }
}

# Regex matching the sweep table row format from test_planar_symmetry.jl:
#   frame      t(s) |    base Vn    base Vs base asym |    corr Vn    corr Vs corr asym | |asym| gain
$rowPattern = '^\s*(\d+)\s+(\d+)\s*\|\s*([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)\s*\|\s*([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)\s*\|\s*([\+\-][\d\.]+)'

function Parse-SweepRows {
    param([string[]]$Lines)
    $rows = @()
    foreach ($line in $Lines) {
        if ($line -match $rowPattern) {
            $rows += [pscustomobject]@{
                Frame    = [int]$Matches[1]
                T        = [int]$Matches[2]
                BaseVn   = [double]$Matches[3]
                BaseVs   = [double]$Matches[4]
                BaseAsym = [double]$Matches[5]
                CorrVn   = [double]$Matches[6]
                CorrVs   = [double]$Matches[7]
                CorrAsym = [double]$Matches[8]
                Gain     = [double]($Matches[9])
            }
        }
    }
    return $rows
}

function Find-FirstFlip {
    <#
        First frame (after the trivial t=0 row) where CorrAsym's sign
        differs from the sign of the first non-trivial CorrAsym value.
        Returns $null if no flip occurs across the sweep.
    #>
    param([array]$Rows)
    if ($Rows.Count -lt 2) { return $null }
    $refSign = [Math]::Sign($Rows[1].CorrAsym)
    if ($refSign -eq 0) { $refSign = 1 }   # treat an exact-zero first row as "positive reference"
    for ($i = 2; $i -lt $Rows.Count; $i++) {
        $s = [Math]::Sign($Rows[$i].CorrAsym)
        if ($s -ne 0 -and $s -ne $refSign) {
            return $Rows[$i]
        }
    }
    return $null
}

$Summary = [System.Collections.Generic.List[object]]::new()

foreach ($dt in $DtValues) {
    Write-Host ""
    Write-Host ("=" * 78) -ForegroundColor Cyan
    Write-Host "  dt_max = $dt" -ForegroundColor Cyan
    Write-Host ("=" * 78) -ForegroundColor Cyan

    $dtTag   = ($dt.ToString([System.Globalization.CultureInfo]::InvariantCulture)) -replace '\.', 'p'
    $outH5   = Join-Path $RepoRoot "test\planar_embankment\planar_res18_diamond_dtsweep_$dtTag.h5"
    $runLog  = Join-Path $LogDir "${Timestamp}_run_dt$dtTag.log"
    $cmpLog  = Join-Path $LogDir "${Timestamp}_compare_dt$dtTag.log"

    Write-Host "  Running simulation (dt_max=$dt) ..." -ForegroundColor DarkGray
    & julia --project=. --threads $Threads FloodModel.jl `
        --meshload $PlanarMesh `
        --flow-model standard `
        --injection-point $InjectionPoint `
        --closed-boundaries `
        --sim-duration $SimDuration `
        --dt-max $dt `
        --gradient-correction on `
        --face-flux-method diamond `
        --output $outH5 `
        --output-interval 300 2>&1 | Tee-Object -FilePath $runLog

    if ($LASTEXITCODE -ne 0) {
        Write-Host "  [FAIL] simulation run (dt_max=$dt), exit $LASTEXITCODE -- see $runLog" -ForegroundColor Red
        $Summary.Add([pscustomobject]@{
            DtMax = $dt; Status = "RUN FAILED"
            FirstT = $null; FirstAsym = $null
            FlipFrame = $null; FlipT = $null
            FinalT = $null; FinalAsym = $null; FinalGain = $null
        })
        continue
    }

    Write-Host "  Comparing against baseline ..." -ForegroundColor DarkGray
    $cmpOutput = & julia --project=. $SymmetryTest `
        --baseline $BaselineH5 `
        --corrected $outH5 `
        --source-lat $InjectionLat `
        --source-lon $InjectionLon `
        --sweep 10 2>&1
    $cmpOutput | Out-File -FilePath $cmpLog
    $cmpExit = $LASTEXITCODE

    if ($cmpExit -ne 0) {
        Write-Host "  [FAIL] comparison (dt_max=$dt), exit $cmpExit -- see $cmpLog" -ForegroundColor Red
        $Summary.Add([pscustomobject]@{
            DtMax = $dt; Status = "COMPARE FAILED"
            FirstT = $null; FirstAsym = $null
            FlipFrame = $null; FlipT = $null
            FinalT = $null; FinalAsym = $null; FinalGain = $null
        })
        continue
    }

    $rows = Parse-SweepRows -Lines $cmpOutput
    if ($rows.Count -lt 2) {
        Write-Host "  [WARN] could not parse sweep table from comparison output -- see $cmpLog" -ForegroundColor Yellow
        $Summary.Add([pscustomobject]@{
            DtMax = $dt; Status = "PARSE FAILED"
            FirstT = $null; FirstAsym = $null
            FlipFrame = $null; FlipT = $null
            FinalT = $null; FinalAsym = $null; FinalGain = $null
        })
        continue
    }

    $first = $rows[1]              # first non-trivial checkpoint (row 0 is t~=10s, ~zero)
    $final = $rows[$rows.Count - 1]
    $flip  = Find-FirstFlip -Rows $rows

    Write-Host "  [OK] dt_max=$dt -> first(t=$($first.T)) asym=$($first.CorrAsym)  " `
               "final(t=$($final.T)) asym=$($final.CorrAsym)  " `
               "flip=$(if ($flip) { "frame $($flip.Frame) (t=$($flip.T))" } else { "never" })" `
               -ForegroundColor Green

    $Summary.Add([pscustomobject]@{
        DtMax     = $dt
        Status    = "OK"
        FirstT    = $first.T
        FirstAsym = $first.CorrAsym
        FlipFrame = if ($flip) { $flip.Frame } else { $null }
        FlipT     = if ($flip) { $flip.T } else { $null }
        FinalT    = $final.T
        FinalAsym = $final.CorrAsym
        FinalGain = $final.Gain
    })
}

Write-Host ""
Write-Host ("=" * 78) -ForegroundColor Cyan
Write-Host "  Summary: diamond method, dt-sensitivity sweep" -ForegroundColor Cyan
Write-Host ("=" * 78) -ForegroundColor Cyan
$Summary | Format-Table -AutoSize -Property `
    @{Label="dt_max"; Expression={$_.DtMax}}, `
    @{Label="status"; Expression={$_.Status}}, `
    @{Label="first_t"; Expression={$_.FirstT}}, `
    @{Label="first_asym"; Expression={ if ($_.FirstAsym -ne $null) { "{0:+0.0000;-0.0000}" -f $_.FirstAsym } }}, `
    @{Label="flip_frame"; Expression={ if ($_.FlipFrame) { $_.FlipFrame } else { "never" } }}, `
    @{Label="flip_t"; Expression={$_.FlipT}}, `
    @{Label="final_t"; Expression={$_.FinalT}}, `
    @{Label="final_asym"; Expression={ if ($_.FinalAsym -ne $null) { "{0:+0.0000;-0.0000}" -f $_.FinalAsym } }}, `
    @{Label="final_gain"; Expression={ if ($_.FinalGain -ne $null) { "{0:+0.0000;-0.0000}" -f $_.FinalGain } }}

Write-Host ""
Write-Host "Compare this table directly against the legacy alpha=1 + cell-momentum" -ForegroundColor DarkGray
Write-Host "sweep already on record (pentagon-chirality handover, section 4):" -ForegroundColor DarkGray
Write-Host "  dt_max=10 -> asym@t~3000s = -0.206" -ForegroundColor DarkGray
Write-Host "  dt_max=4  -> asym@t~3000s = -0.073" -ForegroundColor DarkGray
Write-Host "  dt_max=2  -> asym@t~3000s = -0.010 (best)" -ForegroundColor DarkGray
Write-Host "  dt_max=1  -> asym@t~3000s = -0.145" -ForegroundColor DarkGray
Write-Host "  dt_max=0.5-> asym@t~3000s = -0.344" -ForegroundColor DarkGray
Write-Host "(V-shaped, minimum error at dt~2s -- NOT monotonic in dt)" -ForegroundColor DarkGray
Write-Host ""
Write-Host "Full per-dt logs: $LogDir" -ForegroundColor DarkGray

$Summary | Export-Csv -Path (Join-Path $LogDir "${Timestamp}_summary.csv") -NoTypeInformation
Write-Host "Summary also saved to: $(Join-Path $LogDir "${Timestamp}_summary.csv")" -ForegroundColor DarkGray
