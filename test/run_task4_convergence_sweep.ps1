# run_task4_convergence_sweep.ps1 — Phase A Task 4: CFL-held resolution
# convergence sweep (FloodA5_PhaseA_ImplementationScope.md).
#
# Determines whether the planar-slope north/south volume asymmetry shrinks
# with mesh refinement (ordinary truncation error) or stays flat (stronger
# evidence for genuine A5 topological chirality — consistent with Task 3's
# resolution-flat correlation finding and the "missing east-west edge
# bearing" structural result). Uses the CURRENT, UNCORRECTED solver plumbing
# throughout — this script does not touch the flux kernel, only runs it at
# different resolutions with a properly CFL-scaled --dt-max at each one.
#
# Runs, for each resolution in -Resolutions (default 14, 16, 18):
#   1. Generate the planar-slope DEM ONCE (embankment disabled via
#      --emb-height 0.0 — a pure 0.1% west-to-east slope, matching the
#      domain used throughout the Non-Orthogonal Correction Plan and both
#      Pentagon Chirality handovers), skipped if already present.
#   2. Generate the A5 mesh at this resolution from that SAME DEM
#      (--threads 1, required — mesh generation thread-safety, Bugs 50/57),
#      skipped if already present.
#   3. Compute a CFL-consistent --dt-max via compute_cfl_dt_max.jl, so every
#      resolution starts the sweep at the same Courant-number margin rather
#      than sharing one fixed dt_max that would be proportionally looser at
#      coarse resolution than at fine resolution during the dry startup
#      phase (Task 4's explicit requirement).
#   4. Run the simulation twice: --gradient-correction off (baseline) and
#      --gradient-correction on --gradient-correction-alpha 0.0 (the
#      documented interim default per the Pentagon Chirality handover —
#      NOT alpha=1.0, which that handover found overshoots into a mirrored
#      bias). Both skipped if their .h5 output already exists.
#   5. Run test_planar_symmetry.jl --sweep for a full qualitative view
#      (printed to console, not parsed) AND --frame at a small set of
#      matched checkpoints for a robust, parseable cross-resolution table.
#
# Final output: a table of |asymmetry| vs resolution at each checkpoint,
# for both configurations — the deliverable format Task 4 asks for.
#
# Usage (from project root in PowerShell):
#   .\test\run_task4_convergence_sweep.ps1
#   .\test\run_task4_convergence_sweep.ps1 -Resolutions 14,16,18 -SimDuration 72000
#   .\test\run_task4_convergence_sweep.ps1 -Force                # rerun everything
#   .\test\run_task4_convergence_sweep.ps1 -SkipSim              # only mesh-gen + dt_max
#
# Estimated wall-clock time: mesh generation is the slowest part at res 18
# (several minutes, single-threaded per Bugs 50/57); the six simulation runs
# themselves should be fast (comparable runs at res 18 standard flow have
# been benchmarked well under a minute each — see PROJECT_STATE.md's
# Stage 4 Carlisle regression table for reference figures). Re-running the
# script after a partial/interrupted run will skip anything already
# completed, unless -Force is passed.

param(
    [int[]]   $Resolutions      = @(14, 16, 18),
    [double]  $SimDuration      = 72000,      # seconds (20h) — matches prior sessions
    [double]  $OutputInterval   = 300,        # seconds — matches prior sessions
    [double]  $InjectionRateM3s = 0.1,        # matches Flow Direction Fixes handover
    [double]  $Courant          = 0.7,        # matches _cfl_dt's default
    [double]  $HRef             = 1.0,        # representative depth for dt_max scaling
    [double]  $GradAlpha        = 0.0,        # documented interim default (NOT 1.0)
    [double[]]$DomainKm         = @(4.0, 2.0),
    [double[]]$Centre           = @(-0.017, 51.0),
    [double]  $SlopePct         = 0.1,
    [int[]]   $CheckpointFrac   = @(25, 50, 100),   # % of sim duration, for the final table
    [int]     $SweepStep        = 20,               # frame step for the full --sweep printout
    [switch]  $Force,           # regenerate/rerun everything even if outputs exist
    [switch]  $SkipSim          # stop after mesh-gen + dt_max (useful for a dry run)
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ── Locate julia and python ────────────────────────────────────────────────
$Julia = (Get-Command julia -ErrorAction SilentlyContinue)?.Source
if (-not $Julia) { Write-Error "julia not found on PATH."; exit 1 }
$Python = (Get-Command python -ErrorAction SilentlyContinue)?.Source
if (-not $Python) { Write-Error "python not found on PATH (need the flooda5 conda env active)."; exit 1 }

# ── Paths ────────────────────────────────────────────────────────────────
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$TestDir     = Join-Path $ProjectRoot "test\planar_embankment"
$FloodModel  = Join-Path $ProjectRoot "FloodModel.jl"
$DemGenPy    = Join-Path $TestDir "generate_planar_embankment_dem.py"
$DemFile     = Join-Path $TestDir "planar_embankment_dem.tif"
$ParamsFile  = Join-Path $TestDir "planar_embankment_params.json"
$AoiFile     = Join-Path $TestDir "planar_embankment_aoi.geojson"
$CflScript   = Join-Path $PSScriptRoot "compute_cfl_dt_max.jl"
$SymScript   = Join-Path $PSScriptRoot "test_planar_symmetry.jl"

Write-Host ("=" * 78)
Write-Host "FloodA5 Phase A Task 4 — CFL-held resolution convergence sweep"
Write-Host "  Resolutions      : $($Resolutions -join ', ')"
Write-Host "  Sim duration     : $SimDuration s ($([math]::Round($SimDuration/3600,1)) h)"
Write-Host "  Output interval  : $OutputInterval s"
Write-Host "  Gradient alpha   : $GradAlpha (corrected config)"
Write-Host ("=" * 78)

# ── Stage 1 — DEM (once, embankment disabled: --emb-height 0.0) ───────────
if ($Force -or -not (Test-Path $DemFile) -or -not (Test-Path $ParamsFile) -or -not (Test-Path $AoiFile)) {
    Write-Host "`n[Stage 1] Generating planar-slope DEM (embankment disabled)..."
    & $Python $DemGenPy `
        --out $TestDir `
        --domain-km $DomainKm[0] $DomainKm[1] `
        --centre $Centre[0] $Centre[1] `
        --slope-pct $SlopePct `
        --emb-height 0.0
    if ($LASTEXITCODE -ne 0) { Write-Error "DEM generation failed."; exit 1 }
} else {
    Write-Host "`n[Stage 1] DEM already present — skipping (use -Force to regenerate)."
}

$Params       = Get-Content $ParamsFile | ConvertFrom-Json
$InjectionLat = $Params.injection_lat
$InjectionLon = $Params.injection_lon
Write-Host "  Injection point: lat=$InjectionLat  lon=$InjectionLon (read from $ParamsFile)"

# ── Per-resolution stages ──────────────────────────────────────────────────
$DtMaxByRes = @{}
$MeshByRes  = @{}

foreach ($Res in $Resolutions) {
    Write-Host ("`n" + ("─" * 78))
    Write-Host "Resolution $Res"
    Write-Host ("─" * 78)

    # ── Stage 2 — mesh generation (--threads 1 required) ───────────────────
    $MeshFile = Join-Path $TestDir "planar_mesh${Res}_task4.parquet"
    $MeshByRes[$Res] = $MeshFile
    if ($Force -or -not (Test-Path $MeshFile)) {
        Write-Host "[Stage 2] Generating res $Res mesh (--threads 1, may take a while at high res)..."
        & $Julia --threads 1 --project=$ProjectRoot $FloodModel `
            --meshgen $AoiFile --meshres $Res --dem $DemFile `
            --meshout $MeshFile --flow-model standard --mesh-only
        if ($LASTEXITCODE -ne 0) { Write-Error "Mesh generation failed for res $Res."; exit 1 }
    } else {
        Write-Host "[Stage 2] res $Res mesh already present — skipping."
    }

    # ── Stage 3 — CFL-consistent dt_max ─────────────────────────────────────
    Write-Host "[Stage 3] Computing CFL-consistent --dt-max for res $Res..."
    $DtMaxOut = & $Julia --project=$ProjectRoot $CflScript $MeshFile `
        --courant $Courant --h-ref $HRef
    Write-Host ($DtMaxOut -join "`n")
    $DtMaxQuiet = & $Julia --project=$ProjectRoot $CflScript $MeshFile `
        --courant $Courant --h-ref $HRef --quiet
    $DtMax = [double]($DtMaxQuiet | Select-Object -Last 1)
    $DtMaxByRes[$Res] = $DtMax
    Write-Host "  → using --dt-max $DtMax for res $Res"

    if ($SkipSim) { continue }

    # ── Stage 4 — simulation runs (baseline, corrected) ─────────────────────
    $H5Base = Join-Path $TestDir "planar_res${Res}_task4_baseline.h5"
    $H5Corr = Join-Path $TestDir "planar_res${Res}_task4_corrected.h5"

    $CommonArgs = @(
        "--meshload", $MeshFile,
        "--flow-model", "standard",
        "--injection-point", "$InjectionLat,$InjectionLon,$InjectionRateM3s",
        "--closed-boundaries",
        "--sim-duration", $SimDuration,
        "--dt-max", $DtMax,
        "--output-interval", $OutputInterval
    )

    if ($Force -or -not (Test-Path $H5Base)) {
        Write-Host "[Stage 4] Running res $Res BASELINE (--gradient-correction off)..."
        $StartTime = Get-Date
        & $Julia --threads auto --project=$ProjectRoot $FloodModel @CommonArgs `
            --output $H5Base --gradient-correction off
        if ($LASTEXITCODE -ne 0) { Write-Error "Baseline run failed for res $Res."; exit 1 }
        Write-Host ("  Completed in {0:F1} s" -f ((Get-Date) - $StartTime).TotalSeconds)
    } else {
        Write-Host "[Stage 4] res $Res baseline output already present — skipping."
    }

    if ($Force -or -not (Test-Path $H5Corr)) {
        Write-Host "[Stage 4] Running res $Res CORRECTED (--gradient-correction on --gradient-correction-alpha $GradAlpha)..."
        $StartTime = Get-Date
        & $Julia --threads auto --project=$ProjectRoot $FloodModel @CommonArgs `
            --output $H5Corr --gradient-correction on --gradient-correction-alpha $GradAlpha
        if ($LASTEXITCODE -ne 0) { Write-Error "Corrected run failed for res $Res."; exit 1 }
        Write-Host ("  Completed in {0:F1} s" -f ((Get-Date) - $StartTime).TotalSeconds)
    } else {
        Write-Host "[Stage 4] res $Res corrected output already present — skipping."
    }

    # ── Stage 5a — full qualitative sweep (printed, not parsed) ─────────────
    Write-Host "[Stage 5] Full asymmetry sweep for res $Res (for manual inspection):"
    & $Julia --project=$ProjectRoot $SymScript `
        --baseline $H5Base --corrected $H5Corr `
        --source-lat $InjectionLat --source-lon $InjectionLon `
        --sweep $SweepStep
}

if ($SkipSim) {
    Write-Host "`n-SkipSim was set — stopping before simulation/analysis stages."
    Write-Host "dt_max computed per resolution:"
    foreach ($Res in $Resolutions) { Write-Host "  res $Res : dt_max = $($DtMaxByRes[$Res])" }
    exit 0
}

# ── Stage 5b — targeted checkpoint frames for the final table ──────────────
# Uses --frame (not --sweep) at specific matched checkpoints so the output
# is a small number of clean "asymmetry = ±X.XXXX" lines to regex out,
# rather than parsing a full formatted sweep table.
#
# IMPORTANT: the number of frames actually written does NOT always equal
# floor(sim_duration / output_interval) — adaptive timestepping and
# end-of-run cadence handling can produce a few fewer (observed: 228 actual
# vs. 240 theoretical on a 72000s/300s run in initial testing). Checkpoint
# frame indices are therefore computed from each run's ACTUAL frame count
# (queried directly from the HDF5 file via count_h5_frames.jl), using the
# minimum of baseline's and corrected's counts per resolution as the safe
# upper bound, rather than assumed from sim_duration/output_interval.
$CountScript = Join-Path $PSScriptRoot "count_h5_frames.jl"
$Results = @{}   # key: "$Res|$Pct" -> @{base=...; corr=...; t=...}

Write-Host ("`n" + ("=" * 78))
Write-Host "Stage 6 — Targeted checkpoint comparison (for final table)"
Write-Host ("=" * 78)

foreach ($Res in $Resolutions) {
    $H5Base = Join-Path $TestDir "planar_res${Res}_task4_baseline.h5"
    $H5Corr = Join-Path $TestDir "planar_res${Res}_task4_corrected.h5"

    $NFramesBase = [int](& $Julia --project=$ProjectRoot $CountScript $H5Base | Select-Object -Last 1)
    $NFramesCorr = [int](& $Julia --project=$ProjectRoot $CountScript $H5Corr | Select-Object -Last 1)
    $NFrames = [math]::Min($NFramesBase, $NFramesCorr)
    Write-Host "res $Res : baseline has $NFramesBase frames, corrected has $NFramesCorr frames -> using $NFrames"

    foreach ($Pct in $CheckpointFrac) {
        $FrameIdx = [math]::Max(1, [math]::Min($NFrames, [math]::Round($NFrames * $Pct / 100.0)))
        $Out = & $Julia --project=$ProjectRoot $SymScript `
            --baseline $H5Base --corrected $H5Corr `
            --source-lat $InjectionLat --source-lon $InjectionLon `
            --frame $FrameIdx
        $OutText = $Out -join "`n"

        # Two "asymmetry = ±X.XXXX" lines appear in order: baseline, then corrected.
        $AsymMatches = [regex]::Matches($OutText, "asymmetry\s*=\s*([+\-]?[\d.]+)")
        if ($AsymMatches.Count -lt 2) {
            Write-Warning "Could not parse asymmetry values for res $Res, checkpoint $Pct% (frame $FrameIdx)."
            continue
        }
        $AsymBase = [double]$AsymMatches[0].Groups[1].Value
        $AsymCorr = [double]$AsymMatches[1].Groups[1].Value
        $Results["$Res|$Pct"] = @{ base = $AsymBase; corr = $AsymCorr; frame = $FrameIdx }
    }
}

# ── Final table ──────────────────────────────────────────────────────────
Write-Host ("`n" + ("=" * 78))
Write-Host "FINAL TABLE — |asymmetry| vs resolution (Task 4 deliverable)"
Write-Host ("=" * 78)
Write-Host "dt_max used per resolution (CFL-scaled, courant=$Courant, h_ref=$HRef m):"
foreach ($Res in $Resolutions) { Write-Host "  res $Res : dt_max = $([math]::Round($DtMaxByRes[$Res],2)) s" }
Write-Host ""

foreach ($Pct in $CheckpointFrac) {
    Write-Host "Checkpoint: $Pct% of sim duration"
    Write-Host ("{0,-8}{1,12}{2,12}{3,12}" -f "res", "|base|", "|corr|", "gain")
    foreach ($Res in $Resolutions) {
        $Key = "$Res|$Pct"
        if (-not $Results.ContainsKey($Key)) { continue }
        $R = $Results[$Key]
        $AbsBase = [math]::Abs($R.base)
        $AbsCorr = [math]::Abs($R.corr)
        $Gain    = $AbsBase - $AbsCorr
        $GainStr = if ($Gain -ge 0) { "+{0:F4}" -f $Gain } else { "{0:F4}" -f $Gain }
        Write-Host ("{0,-8}{1,12:F4}{2,12:F4}{3,12}" -f $Res, $AbsBase, $AbsCorr, $GainStr)
    }
    Write-Host ""
}

Write-Host ("=" * 78)
Write-Host "Interpretation (Phase A -> Phase B decision matrix):"
Write-Host "  - |base| SHRINKING clearly as resolution increases (14 -> 16 -> 18) at a"
Write-Host "    given checkpoint points toward ordinary truncation error."
Write-Host "  - |base| staying FLAT or not shrinking meaningfully, combined with Task 3's"
Write-Host "    resolution-flat correlation result, strengthens the genuine-chirality"
Write-Host "    reading and supports proceeding into Phase B with that expectation."
Write-Host "  - 'gain' > 0 at every resolution means the alpha=0.0 correction helps"
Write-Host "    uniformly across resolutions; a gain that shrinks toward zero at finer"
Write-Host "    resolution would suggest the correction is not resolution-robust."
Write-Host ("=" * 78)
