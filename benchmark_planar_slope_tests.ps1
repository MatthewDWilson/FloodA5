<#
.SYNOPSIS
    Batch benchmark FloodA5 planar north/south symmetry against a baseline.

.DESCRIPTION
    Reads experiment definitions from experiment_manifest.csv, identifies
    the latest manifest record for each experiment ID, and runs
    test_planar_symmetry.jl against a user-selected baseline.

    If an experiment ID appears multiple times in the manifest, the final
    occurrence is treated as the authoritative/latest run. A warning is
    issued showing the duplicate count and the previous/latest statuses.

    Only experiments whose latest manifest record has Status=SUCCESS and
    a non-empty OutputFile are eligible for benchmarking. An older SUCCESS
    record is never used if a later record has a different status.

    The Julia benchmark produces machine-readable CSV output. This script
    combines those results with the experimental metadata from the manifest
    and calculates baseline-relative metrics.

.EXAMPLES

    # Benchmark all successful experiments against A01
    .\benchmark_planar_symmetry.ps1 -Baseline A01 -Frame 241

    # Benchmark a subset of experiments
    .\benchmark_planar_symmetry.ps1 `
        -Baseline A01 `
        -ExperimentId A02,B01,B02

    # Benchmark complete experimental sets
    .\benchmark_planar_symmetry.ps1 `
        -Baseline A01 `
        -Set B,C,D

    # Specify a different output CSV
    .\benchmark_planar_symmetry.ps1 `
        -Baseline A01 `
        -Frame 241 `
        -Output benchmark_symmetry.csv

    # Use a different symmetry-line tolerance
    .\benchmark_planar_symmetry.ps1 `
        -Baseline A01 `
        -Frame 241 `
        -LatToleranceM 25

    # Continue after an individual benchmark failure
    .\benchmark_planar_symmetry.ps1 `
        -Baseline A01 `
        -Frame 241 `
        -ContinueOnError
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory)]
    [string]$Baseline,

    [int]$Frame = 241,

    [double]$LatToleranceM = 50.0,

    [string]$Manifest = "test\planar_embankment\planar_slope_tests\experiment_manifest.csv",

    [string]$SymmetryScript = "test\test_planar_symmetry.jl",

    [string]$Output = "test\planar_embankment\planar_slope_tests\planar_symmetry_benchmark.csv",

    [string]$TempDirectory = "test\planar_embankment\planar_slope_tests\benchmark_tmp",

    [string[]]$Set,

    [string[]]$ExperimentId,

    [switch]$ContinueOnError
)

$ErrorActionPreference = "Stop"

# ============================================================================
# Paths
# ============================================================================

$ScriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path

function Resolve-RelativePath {
    param(
        [Parameter(Mandatory)]
        [string]$Path
    )

    if ([System.IO.Path]::IsPathRooted($Path)) {
        return $Path
    }

    return Join-Path $ScriptRoot $Path
}

$ManifestPath = Resolve-RelativePath $Manifest
$SymmetryScriptPath = Resolve-RelativePath $SymmetryScript
$OutputPath = Resolve-RelativePath $Output
$TempDirectoryPath = Resolve-RelativePath $TempDirectory

# ============================================================================
# Validation
# ============================================================================

if (-not (Test-Path $ManifestPath)) {
    throw "Manifest not found: $ManifestPath"
}

if (-not (Test-Path $SymmetryScriptPath)) {
    throw "Symmetry test script not found: $SymmetryScriptPath"
}

$ManifestRows = @(Import-Csv $ManifestPath)

if ($ManifestRows.Count -eq 0) {
    throw "Manifest contains no experiment records: $ManifestPath"
}

# ============================================================================
# Resolve latest manifest record for each experiment
#
# The manifest is treated as a run history. If an experiment ID occurs more
# than once, the LAST occurrence is authoritative.
# ============================================================================

$LatestExperimentRows = @{}
$ExperimentHistory = @{}

foreach ($Row in $ManifestRows) {

    $Id = [string]$Row.ID

    if ([string]::IsNullOrWhiteSpace($Id)) {
        Write-Warning "Manifest contains a row with no experiment ID; skipping that row."
        continue
    }

    if (-not $ExperimentHistory.ContainsKey($Id)) {
        $ExperimentHistory[$Id] = @()
    }

    $ExperimentHistory[$Id] += $Row

    # Deliberately overwrite earlier records: last row wins.
    $LatestExperimentRows[$Id] = $Row
}

# Warn once for every duplicate experiment ID.
foreach ($Id in ($ExperimentHistory.Keys | Sort-Object)) {

    $History = @($ExperimentHistory[$Id])

    if ($History.Count -le 1) {
        continue
    }

    $PreviousRows = @($History[0..($History.Count - 2)])
    $LatestRow = $History[$History.Count - 1]

    $PreviousStatuses = @(
        $PreviousRows |
            ForEach-Object {
                if ([string]::IsNullOrWhiteSpace($_.Status)) {
                    "<blank>"
                }
                else {
                    $_.Status
                }
            }
    )

    $LatestStatus = if ([string]::IsNullOrWhiteSpace($LatestRow.Status)) {
        "<blank>"
    }
    else {
        $LatestRow.Status
    }

    Write-Warning (
        "Experiment '{0}' appears {1} times in {2}; using the last row." -f `
        $Id,
        $History.Count,
        $ManifestPath
    )

    Write-Host (
        "         Previous status: {0}" -f `
        ($PreviousStatuses -join ", ")
    )

    Write-Host (
        "         Latest status:   {0}" -f `
        $LatestStatus
    )
}

# ============================================================================
# Helper functions
# ============================================================================

function Get-ManifestExperiment {
    param(
        [Parameter(Mandatory)]
        [string]$Id
    )

    if (-not $LatestExperimentRows.ContainsKey($Id)) {
        throw "Experiment '$Id' was not found in $ManifestPath."
    }

    return $LatestExperimentRows[$Id]
}


function Test-SuccessfulManifestRow {
    param(
        [Parameter(Mandatory)]
        $Row
    )

    # Only the latest manifest record counts.
    #
    # This deliberately does NOT search backwards for an older SUCCESS row.
    # If the latest run is RUNNING, FAILED, etc., that experiment is not
    # considered available for benchmarking.

    return (
        $Row.Status -eq "SUCCESS" -and
        -not [string]::IsNullOrWhiteSpace($Row.OutputFile)
    )
}


function Resolve-H5Path {
    param(
        [Parameter(Mandatory)]
        [string]$Path
    )

    # Paths in the manifest are normally relative to the FloodA5 project
    # root. Resolve them from the script directory.

    if ([System.IO.Path]::IsPathRooted($Path)) {
        return $Path
    }

    return Join-Path $ScriptRoot $Path
}


function Get-SelectedExperiments {

    # ------------------------------------------------------------------------
    # Explicit experiment IDs
    # ------------------------------------------------------------------------

    if ($ExperimentId) {

        $Selected = @()

        foreach ($Id in $ExperimentId) {

            $Row = Get-ManifestExperiment $Id

            if (-not (Test-SuccessfulManifestRow $Row)) {

                $Status = if ([string]::IsNullOrWhiteSpace($Row.Status)) {
                    "<blank>"
                }
                else {
                    $Row.Status
                }

                throw (
                    "Experiment '$Id' is not available for benchmarking: " +
                    "latest manifest status is '$Status'."
                )
            }

            $Selected += $Row
        }

        return $Selected
    }

    # ------------------------------------------------------------------------
    # Experimental sets
    # ------------------------------------------------------------------------

    if ($Set) {

        $Selected = @()

        foreach ($SetName in $Set) {

            $SetName = $SetName.ToUpper()

            $Rows = @(
                $LatestExperimentRows.Values |
                    Where-Object {
                        $_.Set -eq $SetName
                    }
            )

            if ($Rows.Count -eq 0) {
                throw "No experiments found in set '$SetName'."
            }

            $SuccessfulRows = @(
                $Rows |
                    Where-Object {
                        Test-SuccessfulManifestRow $_
                    }
            )

            foreach ($Row in $Rows) {

                if (-not (Test-SuccessfulManifestRow $Row)) {

                    $Status = if ([string]::IsNullOrWhiteSpace($Row.Status)) {
                        "<blank>"
                    }
                    else {
                        $Row.Status
                    }

                    Write-Host (
                        "Skipping {0} (latest status: {1})" -f `
                        $Row.ID,
                        $Status
                    ) -ForegroundColor DarkGray
                }
            }

            if ($SuccessfulRows.Count -eq 0) {
                throw "No successful experiments found in set '$SetName' (using latest manifest records)."
            }

            $Selected += $SuccessfulRows
        }

        return $Selected
    }

    # ------------------------------------------------------------------------
    # Default: all latest successful experiments
    # ------------------------------------------------------------------------

    $SuccessfulRows = @(
        $LatestExperimentRows.Values |
            Where-Object {
                Test-SuccessfulManifestRow $_
            }
    )

    if ($SuccessfulRows.Count -eq 0) {
        throw "No successful experiments were found in the latest manifest records."
    }

    # Report experiments that exist but whose latest run is not successful.
    foreach ($Row in ($LatestExperimentRows.Values | Sort-Object ID)) {

        if (-not (Test-SuccessfulManifestRow $Row)) {

            $Status = if ([string]::IsNullOrWhiteSpace($Row.Status)) {
                "<blank>"
            }
            else {
                $Row.Status
            }

            Write-Host (
                "Skipping {0} (latest status: {1})" -f `
                $Row.ID,
                $Status
            ) -ForegroundColor DarkGray
        }
    }

    return $SuccessfulRows
}


function Get-NumericValue {
    param(
        $Value
    )

    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) {
        return [double]::NaN
    }

    return [double]$Value
}


function Read-SymmetryCsv {
    param(
        [Parameter(Mandatory)]
        [string]$Path
    )

    if (-not (Test-Path $Path)) {
        throw "Julia benchmark did not create expected CSV: $Path"
    }

    $Rows = @(Import-Csv $Path)

    if ($Rows.Count -eq 0) {
        throw "Julia benchmark CSV is empty: $Path"
    }

    return $Rows
}


function Invoke-SymmetryBenchmark {
    param(
        [Parameter(Mandatory)]
        $BaselineRow,

        [Parameter(Mandatory)]
        $ExperimentRow,

        [Parameter(Mandatory)]
        [string]$TemporaryCsv
    )

    $BaselinePath = Resolve-H5Path $BaselineRow.OutputFile
    $ExperimentPath = Resolve-H5Path $ExperimentRow.OutputFile

    if (-not (Test-Path $BaselinePath)) {
        throw "Baseline HDF5 not found: $BaselinePath"
    }

    if (-not (Test-Path $ExperimentPath)) {
        throw "Experiment HDF5 not found: $ExperimentPath"
    }

    # The injection latitude used by all experiments is encoded in the
    # current test setup. Keep this in one place rather than duplicating it
    # across every manifest record.
    $SourceLat = "51.0001"

    $Arguments = @(
        "--project=."
        $SymmetryScriptPath
        "--baseline", $BaselinePath
        "--corrected", $ExperimentPath
        "--source-lat", $SourceLat
        "--source-lon", "-0.0434"
        "--frame", ([string]$Frame)
        "--lat-tolerance-m", ([string]$LatToleranceM)
        "--csv", $TemporaryCsv
    )

    Write-Host ""
    Write-Host "Running symmetry benchmark: $($ExperimentRow.ID)" -ForegroundColor Cyan
    Write-Host "  baseline : $($BaselineRow.ID) / $BaselinePath"
    Write-Host "  corrected: $($ExperimentRow.ID) / $ExperimentPath"
    Write-Host "  frame    : $Frame"
    Write-Host ""

    & julia @Arguments

    if ($LASTEXITCODE -ne 0) {
        throw "Julia benchmark failed with exit code $LASTEXITCODE."
    }

    return @(Read-SymmetryCsv $TemporaryCsv)
}


function New-ResultRow {
    param(
        [Parameter(Mandatory)]
        $ExperimentRow,

        [Parameter(Mandatory)]
        $SymmetryRow,

        [Parameter(Mandatory)]
        $BaselineSymmetry,

        [Parameter(Mandatory)]
        [string]$Status,

        [string]$ErrorMessage = ""
    )

    $Asymmetry = Get-NumericValue $SymmetryRow.asymmetry
    $AbsAsymmetry = [math]::Abs($Asymmetry)

    $BaselineAbsAsymmetry = [math]::Abs($BaselineSymmetry)

    $Gain = $BaselineAbsAsymmetry - $AbsAsymmetry

    if ($BaselineAbsAsymmetry -gt 1e-12) {
        $RelativeGain = $Gain / $BaselineAbsAsymmetry
    }
    else {
        $RelativeGain = [double]::NaN
    }

    if ($ExperimentRow.ID -eq $Baseline) {
        $Comparison = "BASELINE"
        $Improved = ""
        $Gain = 0.0
        $RelativeGain = 0.0
    }
    else {
        $Comparison = "COMPARED"

        if ($AbsAsymmetry -lt $BaselineAbsAsymmetry) {
            $Improved = "YES"
        }
        elseif ($AbsAsymmetry -gt $BaselineAbsAsymmetry) {
            $Improved = "NO"
        }
        else {
            $Improved = "UNCHANGED"
        }
    }

    return [PSCustomObject]@{
        ID                  = $ExperimentRow.ID
        Set                 = $ExperimentRow.Set
        Name                = $ExperimentRow.Name
        FaceFlux            = $ExperimentRow.FaceFlux
        Gradient            = $ExperimentRow.Gradient
        Alpha               = Get-NumericValue $ExperimentRow.Alpha
        QTheta              = Get-NumericValue $ExperimentRow.QTheta
        Momentum            = $ExperimentRow.Momentum

        Frame               = [int]$SymmetryRow.frame
        Time_s              = Get-NumericValue $SymmetryRow.t

        NorthCells          = [int]$SymmetryRow.n_north
        SouthCells          = [int]$SymmetryRow.n_south
        ExcludedCells       = [int]$SymmetryRow.n_excluded

        NorthVolume_m3      = Get-NumericValue $SymmetryRow.V_north
        SouthVolume_m3      = Get-NumericValue $SymmetryRow.V_south
        SymlineVolume_m3    = Get-NumericValue $SymmetryRow.V_symline
        NorthSouthVolume_m3 = Get-NumericValue $SymmetryRow.V_ns

        NorthPct             = Get-NumericValue $SymmetryRow.pct_north
        SouthPct             = Get-NumericValue $SymmetryRow.pct_south

        Asymmetry            = $Asymmetry
        AbsAsymmetry         = $AbsAsymmetry
        BaselineAbsAsymmetry = $BaselineAbsAsymmetry
        AsymmetryGain        = $Gain
        RelativeGain         = $RelativeGain
        Improved             = $Improved

        Comparison           = $Comparison
        Status               = $Status
        Error                = $ErrorMessage
    }
}


# ============================================================================
# Main
# ============================================================================

try {

    # ------------------------------------------------------------------------
    # Baseline
    # ------------------------------------------------------------------------

    $BaselineRow = Get-ManifestExperiment $Baseline

    if (-not (Test-SuccessfulManifestRow $BaselineRow)) {

        $BaselineStatus = if ([string]::IsNullOrWhiteSpace($BaselineRow.Status)) {
            "<blank>"
        }
        else {
            $BaselineRow.Status
        }

        throw (
            "Baseline '$Baseline' is not available for benchmarking: " +
            "latest manifest status is '$BaselineStatus'."
        )
    }

    $BaselinePath = Resolve-H5Path $BaselineRow.OutputFile

    if (-not (Test-Path $BaselinePath)) {
        throw "Baseline HDF5 not found: $BaselinePath"
    }

    # ------------------------------------------------------------------------
    # Select experiments
    # ------------------------------------------------------------------------

    $Selected = @(Get-SelectedExperiments)

    # Ensure the baseline is included in the final CSV even if the user
    # supplied -Set or -ExperimentId that does not include it.
    $SelectedIds = @($Selected | ForEach-Object { $_.ID })

    if ($Baseline -notin $SelectedIds) {
        $Selected = @($BaselineRow) + $Selected
    }

    # Remove duplicates while retaining selection order.
    $Seen = @{}
    $UniqueSelected = @()

    foreach ($Row in $Selected) {
        if (-not $Seen.ContainsKey($Row.ID)) {
            $Seen[$Row.ID] = $true
            $UniqueSelected += $Row
        }
    }

    $Selected = @($UniqueSelected)

    # ------------------------------------------------------------------------
    # Prepare directories
    # ------------------------------------------------------------------------

    $OutputParent = Split-Path -Parent $OutputPath

    if ($OutputParent -and -not (Test-Path $OutputParent)) {
        New-Item `
            -ItemType Directory `
            -Path $OutputParent `
            -Force | Out-Null
    }

    if (-not (Test-Path $TempDirectoryPath)) {
        New-Item `
            -ItemType Directory `
            -Path $TempDirectoryPath `
            -Force | Out-Null
    }

    # ------------------------------------------------------------------------
    # Summary of intended run
    # ------------------------------------------------------------------------

    Write-Host ""
    Write-Host "FloodA5 planar symmetry benchmark" -ForegroundColor Cyan
    Write-Host "==================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Baseline : $($BaselineRow.ID) / $($BaselineRow.Name)"
    Write-Host "Frame    : $Frame"
    Write-Host "Tolerance: $LatToleranceM m"
    Write-Host "Manifest : $ManifestPath"
    Write-Host "Output   : $OutputPath"
    Write-Host ""
    Write-Host "Experiments:"

    foreach ($Row in $Selected) {
        Write-Host (
            "  {0,-4} {1}" -f `
            $Row.ID,
            $Row.Name
        )
    }

    Write-Host ""

    # ------------------------------------------------------------------------
    # Establish baseline metrics directly from the baseline HDF5.
    # ------------------------------------------------------------------------

    $BaselineTemporaryCsv = Join-Path `
        $TempDirectoryPath `
        "baseline_$Baseline.csv"

    $BaselineArguments = @(
        "--project=."
        $SymmetryScriptPath
        "--baseline", $BaselinePath
        "--source-lat", "51.0001"
        "--source-lon", "-0.0434"
        "--frame", ([string]$Frame)
        "--lat-tolerance-m", ([string]$LatToleranceM)
        "--csv", $BaselineTemporaryCsv
    )

    Write-Host "Analysing baseline $Baseline..." -ForegroundColor Yellow

    & julia @BaselineArguments

    if ($LASTEXITCODE -ne 0) {
        throw "Julia baseline analysis failed with exit code $LASTEXITCODE."
    }

    $BaselineCsvRows = @(Read-SymmetryCsv $BaselineTemporaryCsv)

    $BaselineSymmetryRow = @(
        $BaselineCsvRows |
            Where-Object { $_.case -eq "baseline" }
    )[0]

    if ($null -eq $BaselineSymmetryRow) {
        throw "Could not find baseline row in Julia CSV output."
    }

    $BaselineAsymmetry = Get-NumericValue $BaselineSymmetryRow.asymmetry

    Write-Host ""
    Write-Host (
        "Baseline asymmetry: {0:+0.0000;-0.0000;0.0000}  |asymmetry| = {1:0.0000}" -f `
        $BaselineAsymmetry,
        [math]::Abs($BaselineAsymmetry)
    ) -ForegroundColor Yellow

    # ------------------------------------------------------------------------
    # Benchmark experiments
    # ------------------------------------------------------------------------

    $Results = @()

    foreach ($ExperimentRow in $Selected) {

        $TemporaryCsv = Join-Path `
            $TempDirectoryPath `
            "$($ExperimentRow.ID)_symmetry.csv"

        if (Test-Path $TemporaryCsv) {
            Remove-Item $TemporaryCsv -Force
        }

        # --------------------------------------------------------------------
        # Baseline row
        # --------------------------------------------------------------------

        if ($ExperimentRow.ID -eq $Baseline) {

            $Results += New-ResultRow `
                -ExperimentRow $ExperimentRow `
                -SymmetryRow $BaselineSymmetryRow `
                -BaselineSymmetry $BaselineAsymmetry `
                -Status "BASELINE"

            continue
        }

        # --------------------------------------------------------------------
        # Corrected experiment
        # --------------------------------------------------------------------

        try {

            $SymmetryRows = @(
                Invoke-SymmetryBenchmark `
                    -BaselineRow $BaselineRow `
                    -ExperimentRow $ExperimentRow `
                    -TemporaryCsv $TemporaryCsv
            )

            $CorrectedRow = @(
                $SymmetryRows |
                    Where-Object { $_.case -eq "corrected" }
            )[0]

            if ($null -eq $CorrectedRow) {
                throw "No corrected row was found in Julia benchmark output."
            }

            $Results += New-ResultRow `
                -ExperimentRow $ExperimentRow `
                -SymmetryRow $CorrectedRow `
                -BaselineSymmetry $BaselineAsymmetry `
                -Status "SUCCESS"
        }
        catch {

            $ErrorMessage = $_.Exception.Message

            Write-Host ""
            Write-Host (
                "Benchmark failed for {0}: {1}" -f `
                $ExperimentRow.ID,
                $ErrorMessage
            ) -ForegroundColor Red

            $Results += [PSCustomObject]@{
                ID                   = $ExperimentRow.ID
                Set                  = $ExperimentRow.Set
                Name                 = $ExperimentRow.Name
                FaceFlux             = $ExperimentRow.FaceFlux
                Gradient             = $ExperimentRow.Gradient
                Alpha                = Get-NumericValue $ExperimentRow.Alpha
                QTheta               = Get-NumericValue $ExperimentRow.QTheta
                Momentum             = $ExperimentRow.Momentum

                Frame                = $Frame
                Time_s               = [double]::NaN

                NorthCells           = 0
                SouthCells           = 0
                ExcludedCells        = 0

                NorthVolume_m3       = [double]::NaN
                SouthVolume_m3       = [double]::NaN
                SymlineVolume_m3     = [double]::NaN
                NorthSouthVolume_m3  = [double]::NaN

                NorthPct             = [double]::NaN
                SouthPct             = [double]::NaN

                Asymmetry            = [double]::NaN
                AbsAsymmetry         = [double]::NaN
                BaselineAbsAsymmetry = [math]::Abs($BaselineAsymmetry)
                AsymmetryGain        = [double]::NaN
                RelativeGain         = [double]::NaN
                Improved             = ""

                Comparison            = "COMPARED"
                Status                = "FAILED"
                Error                 = $ErrorMessage
            }

            if (-not $ContinueOnError) {
                throw
            }
        }
    }

    # =========================================================================
    # Write consolidated CSV
    # =========================================================================

    $Results |
        Export-Csv `
            -Path $OutputPath `
            -NoTypeInformation `
            -Encoding UTF8

    # =========================================================================
    # Summary
    # =========================================================================

    $Compared = @(
        $Results |
            Where-Object {
                $_.Comparison -eq "COMPARED" -and
                $_.Status -eq "SUCCESS"
            }
    )

    $Improved = @(
        $Compared |
            Where-Object {
                $_.Improved -eq "YES"
            }
    )

    $Worsened = @(
        $Compared |
            Where-Object {
                $_.Improved -eq "NO"
            }
    )

    Write-Host ""
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host "Planar symmetry benchmark summary" -ForegroundColor Cyan
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host ""

    Write-Host (
        "Baseline: {0} ({1})" -f `
        $BaselineRow.ID,
        $BaselineRow.Name
    )

    Write-Host (
        "Baseline asymmetry: {0:+0.0000;-0.0000;0.0000}" -f `
        $BaselineAsymmetry
    )

    Write-Host ""

    $SortedResults = @(
        $Results |
            Sort-Object `
                @{Expression={
                    if ($_.ID -eq $Baseline) {
                        [double]::PositiveInfinity
                    }
                    elseif ($_.Status -ne "SUCCESS") {
                        [double]::NegativeInfinity
                    }
                    else {
                        [double]$_.AsymmetryGain
                    }
                }; Descending=$true}
    )

    Write-Host (
        "{0,-5} {1,-30} {2,11} {3,11} {4,11} {5,12}" -f `
        "ID",
        "Name",
        "|asym|",
        "gain",
        "rel. gain",
        "result"
    )

    Write-Host ("-" * 86)

    foreach ($Row in $SortedResults) {

        if ($Row.Status -eq "FAILED") {

            Write-Host (
                "{0,-5} {1,-30} {2,11} {3,11} {4,11} {5,12}" -f `
                $Row.ID,
                $Row.Name,
                "FAILED",
                "",
                "",
                "FAILED"
            ) -ForegroundColor Red

            continue
        }

        $ResultText = switch ($Row.Improved) {
            "YES"       { "IMPROVED" }
            "NO"        { "WORSE" }
            "UNCHANGED" { "UNCHANGED" }
            default     { "BASELINE" }
        }

        $Colour = switch ($Row.Improved) {
            "YES" { "Green" }
            "NO"  { "Red" }
            default { "Gray" }
        }

        Write-Host (
            "{0,-5} {1,-30} {2,11:0.0000} {3,11:+0.0000;-0.0000;0.0000} {4,11:+0.1%;-0.1%;0.0%} {5,12}" -f `
            $Row.ID,
            $Row.Name,
            [double]$Row.AbsAsymmetry,
            [double]$Row.AsymmetryGain,
            [double]$Row.RelativeGain,
            $ResultText
        ) -ForegroundColor $Colour
    }

    Write-Host ""
    Write-Host "Compared successfully: $($Compared.Count)"
    Write-Host "Improved:              $($Improved.Count)" -ForegroundColor Green
    Write-Host "Worsened:              $($Worsened.Count)" -ForegroundColor Red
    Write-Host ""

    Write-Host "Results written to:"
    Write-Host "  $OutputPath" -ForegroundColor Cyan
    Write-Host ""

    Write-Host "Temporary Julia CSV files:"
    Write-Host "  $TempDirectoryPath" -ForegroundColor DarkGray
    Write-Host ""

    if ($Worsened.Count -gt 0) {
        Write-Host (
            "WARNING: {0} experiment(s) have larger |asymmetry| than the baseline." -f `
            $Worsened.Count
        ) -ForegroundColor Red
    }

    Write-Host ""
}
catch {

    Write-Host ""
    Write-Host "FATAL ERROR" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Write-Host ""

    exit 1
}
