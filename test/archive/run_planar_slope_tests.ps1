<#
.SYNOPSIS
    Run parameterised FloodA5 test experiments.

.DESCRIPTION
    Provides a single source of truth for the FloodA5 experimental matrix.

    The script:
      - Defines experiments by scientific parameters rather than raw commands.
      - Generates the corresponding Julia command automatically.
      - Supports running all experiments or selected experimental sets.
      - Provides a -List option to inspect the experiment matrix.
      - Generates systematic HDF5 output filenames.
      - Writes a CSV manifest recording every run and its status.
      - Prints each generated command before execution.

.EXAMPLES

    # List all available experiments
    .\run_tests.ps1 -List

    # Run all experiments
    .\run_tests.ps1 -All

    # Run one experimental set
    .\run_tests.ps1 -Set A

    # Run several experimental sets
    .\run_tests.ps1 -Set B,C

    # Run selected experiment IDs
    .\run_tests.ps1 -ExperimentId A01,A02,B01

    # Run all experiments without stopping if one fails
    .\run_tests.ps1 -All -ContinueOnError

.NOTES
    Assumes PowerShell 7+ and that Julia is available on PATH.

    All paths are relative to the directory containing this script.
#>

[CmdletBinding()]
param(
    [switch]$All,

    [switch]$List,

    [string[]]$Set,

    [string[]]$ExperimentId,

    [switch]$ContinueOnError
)

# ============================================================================
# Configuration
# ============================================================================

$ErrorActionPreference = "Stop"

$ScriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path

$Julia = "julia"
$Project = "."
$ModelScript = "FloodModel.jl"

$Mesh = "test\planar_embankment\planar_mesh18_std_smooth.parquet"

$OutputDirectory = "test\planar_embankment\planar_slope_tests"
$ManifestPath = Join-Path $OutputDirectory "experiment_manifest.csv"

$CommonArguments = @(
    "--meshload", $Mesh
    "--flow-model", "standard"
    "--injection-point", "51.0001,-0.0434,0.1"
    "--closed-boundaries"
    "--sim-duration", "72000"
    "--dt-max", "2"
    "--output-interval", "300"
)

# ============================================================================
# Experimental matrix
#
# Gradient values:
#   none       = gradient correction OFF
#   orthogonal = correction ON, alpha = 0
#   wlsq       = correction ON, alpha = 1
#
# QTheta:
#   0.9 = q-centred predictor
#   1.0 = pure Bates formulation
#
# FaceFlux:
#   legacy
#   diamond
#
# Momentum:
#   edge
#   cell
# ============================================================================

$Experiments = @(

    # ------------------------------------------------------------------------
    # Set A — Baseline / q-centering
    # ------------------------------------------------------------------------

    [ordered]@{
        ID          = "A01"
        Set         = "A"
        Name        = "baseline"
        Description = "Legacy face flux, no gradient correction, q-centred, edge momentum"
        FaceFlux    = "legacy"
        Gradient    = "none"
        QTheta      = 0.9
        Momentum    = "edge"
    },

    [ordered]@{
        ID          = "A02"
        Set         = "A"
        Name        = "pure_bates"
        Description = "Legacy face flux, no gradient correction, pure Bates, edge momentum"
        FaceFlux    = "legacy"
        Gradient    = "none"
        QTheta      = 1.0
        Momentum    = "edge"
    },


    # ------------------------------------------------------------------------
    # Set B — WLSQ / gradient correction
    # ------------------------------------------------------------------------

    [ordered]@{
        ID          = "B01"
        Set         = "B"
        Name        = "legacy_edge_wlsq_a0"
        Description = "Legacy face flux, orthogonal-only gradient correction, q-centred, edge momentum"
        FaceFlux    = "legacy"
        Gradient    = "orthogonal"
        QTheta      = 0.9
        Momentum    = "edge"
    },

    [ordered]@{
        ID          = "B02"
        Set         = "B"
        Name        = "legacy_edge_wlsq_a1"
        Description = "Legacy face flux, full WLSQ correction, q-centred, edge momentum"
        FaceFlux    = "legacy"
        Gradient    = "wlsq"
        QTheta      = 0.9
        Momentum    = "edge"
    },

    [ordered]@{
        ID          = "B03"
        Set         = "B"
        Name        = "legacy_edge_wlsq_a0_bates"
        Description = "Legacy face flux, orthogonal-only gradient correction, pure Bates, edge momentum"
        FaceFlux    = "legacy"
        Gradient    = "orthogonal"
        QTheta      = 1.0
        Momentum    = "edge"
    },

    [ordered]@{
        ID          = "B04"
        Set         = "B"
        Name        = "legacy_edge_wlsq_a1_bates"
        Description = "Legacy face flux, full WLSQ correction, pure Bates, edge momentum"
        FaceFlux    = "legacy"
        Gradient    = "wlsq"
        QTheta      = 1.0
        Momentum    = "edge"
    },


    # ------------------------------------------------------------------------
    # Set C — Face flux formulation
    # ------------------------------------------------------------------------

    [ordered]@{
        ID          = "C01"
        Set         = "C"
        Name        = "diamond_edge_qtheta09"
        Description = "Diamond face flux, orthogonal-only gradient treatment, q-centred, edge momentum"
        FaceFlux    = "diamond"
        Gradient    = "orthogonal"
        QTheta      = 0.9
        Momentum    = "edge"
    },

    [ordered]@{
        ID          = "C02"
        Set         = "C"
        Name        = "diamond_edge_qtheta10"
        Description = "Diamond face flux, orthogonal-only gradient treatment, pure Bates, edge momentum"
        FaceFlux    = "diamond"
        Gradient    = "orthogonal"
        QTheta      = 1.0
        Momentum    = "edge"
    },


    # ------------------------------------------------------------------------
    # Set D — Cell momentum
    # ------------------------------------------------------------------------

    [ordered]@{
        ID          = "D01"
        Set         = "D"
        Name        = "legacy_cellmom_baseline"
        Description = "Legacy face flux, no gradient correction, pure Bates, cell momentum"
        FaceFlux    = "legacy"
        Gradient    = "none"
        QTheta      = 1.0
        Momentum    = "cell"
    },

    [ordered]@{
        ID          = "D02"
        Set         = "D"
        Name        = "legacy_cellmom_wlsq_a0"
        Description = "Legacy face flux, orthogonal-only gradient correction, pure Bates, cell momentum"
        FaceFlux    = "legacy"
        Gradient    = "orthogonal"
        QTheta      = 1.0
        Momentum    = "cell"
    },

    [ordered]@{
        ID          = "D03"
        Set         = "D"
        Name        = "legacy_cellmom_wlsq_a1"
        Description = "Legacy face flux, full WLSQ correction, pure Bates, cell momentum"
        FaceFlux    = "legacy"
        Gradient    = "wlsq"
        QTheta      = 1.0
        Momentum    = "cell"
    },

    [ordered]@{
        ID          = "D04"
        Set         = "D"
        Name        = "diamond_cellmom"
        Description = "Diamond face flux, orthogonal-only gradient treatment, pure Bates, cell momentum"
        FaceFlux    = "diamond"
        Gradient    = "orthogonal"
        QTheta      = 1.0
        Momentum    = "cell"
    }
)

# ============================================================================
# Helper functions
# ============================================================================

function Get-ExperimentById {
    param(
        [Parameter(Mandatory)]
        [string]$Id
    )

    return $Experiments | Where-Object { $_.ID -eq $Id }
}


function Get-ExperimentOutputName {
    param(
        [Parameter(Mandatory)]
        [System.Collections.IDictionary]$Exp
    )

    return "flooda5_$($Exp.ID)_$($Exp.Name).h5"
}


function Get-GradientArguments {
    param(
        [Parameter(Mandatory)]
        [System.Collections.IDictionary]$Exp
    )

    switch ($Exp.Gradient) {

        "none" {
            return @(
                "--gradient-correction", "off"
                "--gradient-correction-alpha", "0.0"
            )
        }

        "orthogonal" {
            return @(
                "--gradient-correction", "on"
                "--gradient-correction-alpha", "0.0"
            )
        }

        "wlsq" {
            return @(
                "--gradient-correction", "on"
                "--gradient-correction-alpha", "1.0"
            )
        }

        default {
            throw "Unknown gradient mode '$($Exp.Gradient)' in experiment $($Exp.ID)."
        }
    }
}


function Get-ExperimentArguments {
    param(
        [Parameter(Mandatory)]
        [System.Collections.IDictionary]$Exp
    )

    $OutputName = Get-ExperimentOutputName $Exp
    $OutputPath = Join-Path $OutputDirectory $OutputName

    $Arguments = @()

    $Arguments += $CommonArguments

    $Arguments += Get-GradientArguments $Exp

    $Arguments += @(
        "--q-centre-theta", ([string]$Exp.QTheta)
        "--face-flux-method", $Exp.FaceFlux
        "--momentum-model", $Exp.Momentum
        "--output", $OutputPath
    )

    return $Arguments
}


function Format-Command {
    param(
        [Parameter(Mandatory)]
        [string[]]$Arguments
    )

    $Parts = @(
        "$Julia --project=$Project $ModelScript"
    )

    foreach ($Argument in $Arguments) {

        if ($Argument -match '[\s"]') {
            $Parts += '"' + $Argument.Replace('"', '\"') + '"'
        }
        else {
            $Parts += $Argument
        }
    }

    return ($Parts -join " ")
}


function Ensure-OutputDirectory {

    if (-not (Test-Path $OutputDirectory)) {
        New-Item `
            -ItemType Directory `
            -Path $OutputDirectory `
            -Force | Out-Null
    }
}


function Write-ManifestRow {
    param(
        [Parameter(Mandatory)]
        [System.Collections.IDictionary]$Row
    )

    $FileExists = Test-Path $ManifestPath

    $RowObject = [PSCustomObject]$Row

    if ($FileExists) {

        $RowObject | Export-Csv `
            -Path $ManifestPath `
            -NoTypeInformation `
            -Append `
            -Encoding UTF8
    }
    else {

        $RowObject | Export-Csv `
            -Path $ManifestPath `
            -NoTypeInformation `
            -Encoding UTF8
    }
}


function Get-GradientAlpha {
    param(
        [Parameter(Mandatory)]
        [System.Collections.IDictionary]$Exp
    )

    switch ($Exp.Gradient) {
        "none"       { return 0.0 }
        "orthogonal" { return 0.0 }
        "wlsq"       { return 1.0 }
    }
}


function Show-ExperimentList {

    Write-Host ""
    Write-Host "FloodA5 experiment matrix" -ForegroundColor Cyan
    Write-Host "========================="
    Write-Host ""

    $Groups = $Experiments | Group-Object Set

    foreach ($Group in $Groups) {

        switch ($Group.Name) {

            "A" {
                $SetDescription = "Baseline / q-centering"
            }

            "B" {
                $SetDescription = "WLSQ / gradient correction"
            }

            "C" {
                $SetDescription = "Face flux formulation"
            }

            "D" {
                $SetDescription = "Cell momentum"
            }

            default {
                $SetDescription = ""
            }
        }

        Write-Host "Set $($Group.Name) — $SetDescription" -ForegroundColor Yellow

        foreach ($Exp in $Group.Group) {

            $Alpha = Get-GradientAlpha $Exp

            Write-Host (
                "  {0,-4} {1,-32} face={2,-7} gradient={3,-11} alpha={4,-2} theta={5,-3} momentum={6}" -f `
                $Exp.ID,
                $Exp.Name,
                $Exp.FaceFlux,
                $Exp.Gradient,
                $Alpha,
                $Exp.QTheta,
                $Exp.Momentum
            )
        }

        Write-Host ""
    }

    Write-Host "Total experiments: $($Experiments.Count)" -ForegroundColor Green
    Write-Host ""
}


function Resolve-SelectedExperiments {

    # ------------------------------------------------------------------------
    # No selection supplied
    # ------------------------------------------------------------------------

    if (-not $All -and -not $Set -and -not $ExperimentId) {

        Write-Host ""
        Write-Host "No experiments selected." -ForegroundColor Yellow
        Write-Host ""
        Write-Host "Use one of:"
        Write-Host "  .\run_tests.ps1 -List"
        Write-Host "  .\run_tests.ps1 -All"
        Write-Host "  .\run_tests.ps1 -Set A"
        Write-Host "  .\run_tests.ps1 -Set B,C"
        Write-Host "  .\run_tests.ps1 -ExperimentId A01,A02"
        Write-Host ""

        exit 0
    }

    # ------------------------------------------------------------------------
    # -All
    # ------------------------------------------------------------------------

    if ($All) {
        return $Experiments
    }

    # ------------------------------------------------------------------------
    # Explicit experiment IDs
    # ------------------------------------------------------------------------

    if ($ExperimentId) {

        $Selected = @()

        foreach ($Id in $ExperimentId) {

            $Match = Get-ExperimentById $Id

            if (-not $Match) {
                throw "Unknown experiment ID '$Id'. Use -List to see available IDs."
            }

            $Selected += $Match
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

            $Matches = @(
                $Experiments |
                    Where-Object { $_.Set -eq $SetName }
            )

            if ($Matches.Count -eq 0) {
                throw "Unknown experiment set '$SetName'. Available sets: A, B, C, D."
            }

            $Selected += $Matches
        }

        return $Selected
    }
}


function Invoke-Experiment {
    param(
        [Parameter(Mandatory)]
        [System.Collections.IDictionary]$Exp
    )

    $OutputName = Get-ExperimentOutputName $Exp
    $OutputPath = Join-Path $OutputDirectory $OutputName

    $Arguments = Get-ExperimentArguments $Exp
    $Command = Format-Command $Arguments

    $StartTime = Get-Date

    Write-Host ""
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host "Running $($Exp.ID): $($Exp.Name)" -ForegroundColor Cyan
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host $Exp.Description -ForegroundColor DarkGray
    Write-Host ""
    Write-Host "Command:" -ForegroundColor Yellow
    Write-Host $Command
    Write-Host ""
    Write-Host "Output: $OutputPath"
    Write-Host ""

    # ------------------------------------------------------------------------
    # Record that the run has started.
    # ------------------------------------------------------------------------

    $StartRow = [ordered]@{
        ID         = $Exp.ID
        Set        = $Exp.Set
        Name       = $Exp.Name
        FaceFlux   = $Exp.FaceFlux
        Gradient   = $Exp.Gradient
        Alpha      = Get-GradientAlpha $Exp
        QTheta     = $Exp.QTheta
        Momentum   = $Exp.Momentum
        OutputFile = $OutputPath
        Status     = "RUNNING"
        StartTime  = $StartTime.ToString("o")
        EndTime    = ""
        ExitCode   = ""
        Command    = $Command
    }

    Write-ManifestRow $StartRow

    # ------------------------------------------------------------------------
    # Run Julia.
    # ------------------------------------------------------------------------

    try {

        & $Julia "--project=$Project" $ModelScript @Arguments

        $ExitCode = $LASTEXITCODE
        $EndTime = Get-Date

        if ($ExitCode -eq 0) {
            $Status = "SUCCESS"
            $StatusColour = "Green"
        }
        else {
            $Status = "FAILED"
            $StatusColour = "Red"
        }
    }
    catch {

        $ExitCode = -1
        $EndTime = Get-Date
        $Status = "ERROR"
        $StatusColour = "Red"

        Write-Host ""
        Write-Host "Exception while running $($Exp.ID):" -ForegroundColor Red
        Write-Host $_.Exception.Message -ForegroundColor Red
    }

    # ------------------------------------------------------------------------
    # Record final status.
    # ------------------------------------------------------------------------

    $FinalRow = [ordered]@{
        ID         = $Exp.ID
        Set        = $Exp.Set
        Name       = $Exp.Name
        FaceFlux   = $Exp.FaceFlux
        Gradient   = $Exp.Gradient
        Alpha      = Get-GradientAlpha $Exp
        QTheta     = $Exp.QTheta
        Momentum   = $Exp.Momentum
        OutputFile = $OutputPath
        Status     = $Status
        StartTime  = $StartTime.ToString("o")
        EndTime    = $EndTime.ToString("o")
        ExitCode   = $ExitCode
        Command    = $Command
    }

    Write-ManifestRow $FinalRow

    Write-Host ""
    Write-Host "Status: $Status" -ForegroundColor $StatusColour
    Write-Host "Exit code: $ExitCode" -ForegroundColor $StatusColour
    Write-Host "Elapsed: $($EndTime - $StartTime)"
    Write-Host ""

    return ($Status -eq "SUCCESS")
}


# ============================================================================
# Main
# ============================================================================

try {

    Ensure-OutputDirectory

    # ------------------------------------------------------------------------
    # -List takes precedence and does not run anything.
    # ------------------------------------------------------------------------

    if ($List) {
        Show-ExperimentList
        exit 0
    }

    # ------------------------------------------------------------------------
    # Resolve selected experiments.
    # ------------------------------------------------------------------------

    $SelectedExperiments = @(Resolve-SelectedExperiments)

    Write-Host ""
    Write-Host "FloodA5 experiment runner" -ForegroundColor Cyan
    Write-Host "========================="
    Write-Host ""
    Write-Host "Selected experiments: $($SelectedExperiments.Count)"
    Write-Host "Manifest: $ManifestPath"
    Write-Host ""

    foreach ($Exp in $SelectedExperiments) {

        Write-Host (
            "  {0,-4} {1}" -f `
            $Exp.ID,
            $Exp.Name
        )
    }

    Write-Host ""

    # ------------------------------------------------------------------------
    # Run experiments.
    # ------------------------------------------------------------------------

    $Results = @()

    foreach ($Exp in $SelectedExperiments) {

        $Success = Invoke-Experiment $Exp

        $Results += [PSCustomObject]@{
            ID      = $Exp.ID
            Name    = $Exp.Name
            Success = $Success
        }

        if (-not $Success -and -not $ContinueOnError) {

            Write-Host ""
            Write-Host "Experiment $($Exp.ID) failed." -ForegroundColor Red
            Write-Host "Stopping because -ContinueOnError was not specified." -ForegroundColor Red
            Write-Host ""

            exit 1
        }
    }

    # ------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------

    $Successful = @(
        $Results |
            Where-Object { $_.Success }
    ).Count

    $Failed = $Results.Count - $Successful

    Write-Host ""
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host "Experiment summary" -ForegroundColor Cyan
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host ""

    foreach ($Result in $Results) {

        if ($Result.Success) {
            $Symbol = "OK"
            $Colour = "Green"
        }
        else {
            $Symbol = "FAILED"
            $Colour = "Red"
        }

        Write-Host (
            "  {0,-8} {1,-32} {2}" -f `
            $Result.ID,
            $Result.Name,
            $Symbol
        ) -ForegroundColor $Colour
    }

    Write-Host ""
    Write-Host "Successful: $Successful" -ForegroundColor Green

    if ($Failed -gt 0) {
        Write-Host "Failed:     $Failed" -ForegroundColor Red
    }
    else {
        Write-Host "Failed:     $Failed" -ForegroundColor Green
    }

    Write-Host ""
    Write-Host "Manifest: $ManifestPath"
    Write-Host ""

    if ($Failed -gt 0) {
        exit 1
    }
    else {
        exit 0
    }
}
catch {

    Write-Host ""
    Write-Host "FATAL ERROR" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Write-Host ""

    exit 1
}
