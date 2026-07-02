#!/usr/bin/env julia
# test_planar_symmetry.jl
# ------------------------
# Planar slope north/south VOLUME SYMMETRY benchmark — a more robust
# alternative to test_planar_slope.jl's per-cell velocity-direction metric
# (FloodA5_NonOrthogonal_Correction_Plan.md §10 Stage 3).
#
# Rationale: on a DEM that varies only with longitude (a pure east-west
# slope), the continuous physical problem is translation-invariant in
# latitude. A point source on such a slope therefore has NO physically
# preferred north/south direction at all — any imbalance between the total
# water volume north vs. south of the source's latitude is purely a
# numerical/mesh artifact (exactly the non-orthogonal directional bias this
# branch is trying to fix), not a feature of the physics.
#
# This sidesteps the problem found with the velocity-direction test
# (test_planar_slope.jl): that test averages instantaneous per-cell
# DIRECTION vectors, which turned out to be extremely noisy frame-to-frame
# (R never exceeded ~0.13 across a full sweep, with no trend). Volume is a
# conserved, INTEGRATED scalar quantity summed over hundreds of cells, so
# per-cell noise averages out — this should give a much cleaner signal.
#
# This test does NOT run the model. Generate the two comparison runs first
# (reuse the same baseline/corrected .h5 files as test_planar_slope.jl —
# no new simulation runs needed).
#
# Usage:
#   julia --project=. test/test_planar_symmetry.jl \
#       --baseline  test/planar_embankment/planar_baseline.h5 \
#       --corrected test/planar_embankment/planar_corrected.h5 \
#       --source-lat 51.00033 --source-lon -0.04329 \
#       --frame 120
#
#   julia --project=. test/test_planar_symmetry.jl \
#       --baseline  test/planar_embankment/planar_baseline.h5 \
#       --corrected test/planar_embankment/planar_corrected.h5 \
#       --source-lat 51.00033 --source-lon -0.04329 \
#       --sweep 10
#
#   julia --project=. test/test_planar_symmetry.jl \
#       --baseline test/planar_embankment/planar_baseline.h5 --list-frames
#
# Metric: cells are split into "north" (center_lat > source_lat + tol) and
# "south" (center_lat < source_lat − tol) using a latitude tolerance
# converted to metres, excluding cells whose centre falls within that
# tolerance of the source's latitude line (the user's suggestion — cells
# straddling the symmetry line shouldn't be forced into either bucket).
#
#   asymmetry = (V_north - V_south) / (V_north + V_south)
#
# asymmetry = 0 is the ideal (perfectly symmetric); |asymmetry| > 0
# indicates directional bias. Sign indicates which side is over-favoured,
# which is itself a useful diagnostic (a consistent sign across frames
# would point at a systematic bias direction in the mesh, e.g. linked to
# the dominant skew axis reported by the mesh-build diagnostics).
#
# --source-lon is accepted but not currently used by the metric itself
# (only latitude matters for a north/south split) — it's there so the
# same --injection-point values you already used for the simulation runs
# can be pasted in directly without re-deriving just the latitude.

using HDF5
using Statistics
using Printf

const M_PER_DEG_LAT = 111_320.0   # matches the convention used elsewhere
                                   # in the codebase (generate_planar_
                                   # embankment_dem.py, _build_wlsq_weights!)

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------

function _argval(args::Vector{String}, flag::String, default=nothing)
    idx = findfirst(==(flag), args)
    idx === nothing && return default
    idx == length(args) && error("$flag requires a value")
    return args[idx + 1]
end

const ARGS_V = String.(ARGS)
const DO_LIST = "--list-frames" in ARGS_V

# ---------------------------------------------------------------------------
# HDF5 access
# ---------------------------------------------------------------------------

function _frame_names(f)
    return sort(collect(keys(f["frames"])))
end

function _list_frames(h5_path::String)
    h5open(h5_path, "r") do f
        names = _frame_names(f)
        println("Frames in $h5_path:")
        @printf("%8s  %10s  %8s  %10s\n", "frame#", "t (s)", "n_wet", "max_depth")
        for nm in names
            grp = f["frames/$nm"]
            t = haskey(grp, "t") ? read(grp["t"]) : NaN
            sat = read(grp["saturation"])
            n_wet = count(>=(1.0), sat)
            depth = read(grp["water_depth"])
            md = isempty(depth) ? NaN : maximum(depth)
            @printf("%8s  %10.1f  %8d  %10.4f\n", nm, t, n_wet, md)
        end
    end
end

"""
    _read_frame_volume(h5_path, frame_idx) → (center_lats, volume, t)

Reads mesh latitudes (static — same for every frame) and the requested
frame's per-cell volume.
"""
function _read_frame_volume(h5_path::String, frame_idx::Int)
    h5open(h5_path, "r") do f
        center_lats = read(f["mesh/center_lats"])
        names = _frame_names(f)
        fname = lpad(string(frame_idx), 6, '0')
        haskey(f["frames"], fname) ||
            error("Frame $fname not found in $h5_path. Available: " *
                  "$(first(names, 3))…$(last(names, 3)) ($(length(names)) total). " *
                  "Run with --list-frames to inspect.")
        grp = f["frames/$fname"]
        volume = read(grp["volume"])
        t = haskey(grp, "t") ? read(grp["t"]) : NaN
        return center_lats, volume, t
    end
end

# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

"""
    analyse_symmetry(h5_path, frame_idx; source_lat, lat_tolerance_m=50.0)

Splits cells into north/south of `source_lat` (excluding a tolerance band)
and sums volume on each side.

NOTE on the denominator: `V_ns` is `V_north + V_south` — it deliberately
EXCLUDES the volume sitting in the tolerance-band ("symline") cells. This
is what makes `asymmetry` insensitive to the tolerance choice in the way
that matters: widening the tolerance band can only move volume between
"north/south" and "symline" (changing what's excluded), it can never sneak
extra always-included volume into the denominator that would shift the
ratio for a reason unrelated to north/south balance. `V_symline` is
reported separately, for transparency, but is never part of the ratio.
"""
function analyse_symmetry(h5_path::String, frame_idx::Int;
                           source_lat::Float64, lat_tolerance_m::Float64=50.0)
    center_lats, volume, t = _read_frame_volume(h5_path, frame_idx)
    n = length(center_lats)

    dlat_m = (center_lats .- source_lat) .* M_PER_DEG_LAT
    north_mask   = dlat_m .>  lat_tolerance_m
    south_mask   = dlat_m .< -lat_tolerance_m
    symline_mask = .!north_mask .& .!south_mask
    n_excluded   = count(symline_mask)

    V_north   = sum(volume[north_mask])
    V_south   = sum(volume[south_mask])
    V_symline = sum(volume[symline_mask])
    V_ns      = V_north + V_south   # denominator — excludes V_symline by construction

    asymmetry = V_ns > 1e-9 ? (V_north - V_south) / V_ns : 0.0

    return (t=t, n_north=count(north_mask), n_south=count(south_mask),
            n_excluded=n_excluded, V_north=V_north, V_south=V_south,
            V_symline=V_symline, V_ns=V_ns, asymmetry=asymmetry)
end

function _print_result(label::String, r)
    println(label, ":")
    @printf("  t = %.1f s\n", r.t)
    @printf("  cells:  north=%d  south=%d  excluded (within tolerance band)=%d\n",
            r.n_north, r.n_south, r.n_excluded)
    pct_n = r.V_ns > 1e-9 ? 100 * r.V_north / r.V_ns : NaN
    pct_s = r.V_ns > 1e-9 ? 100 * r.V_south / r.V_ns : NaN
    @printf("  volume: north=%.1f m³ (%.1f%%)   south=%.1f m³ (%.1f%%)   north+south=%.1f m³\n",
            r.V_north, pct_n, r.V_south, pct_s, r.V_ns)
    @printf("  (symline band volume, excluded from the ratio above: %.1f m³)\n", r.V_symline)
    @printf("  asymmetry = %+.4f   (0 = perfectly symmetric, i.e. 50%%/50%%; sign indicates which side is favoured)\n",
            r.asymmetry)
    println()
end

# ---------------------------------------------------------------------------
# Main — single frame
# ---------------------------------------------------------------------------

function main()
    if DO_LIST
        for flag in ("--baseline", "--corrected")
            path = _argval(ARGS_V, flag)
            path !== nothing && isfile(path) && _list_frames(path)
        end
        return
    end

    baseline_path   = _argval(ARGS_V, "--baseline")
    corrected_path  = _argval(ARGS_V, "--corrected")
    frame_str       = _argval(ARGS_V, "--frame")
    sweep_str       = _argval(ARGS_V, "--sweep")
    source_lat_str  = _argval(ARGS_V, "--source-lat")
    tol_str         = _argval(ARGS_V, "--lat-tolerance-m", "50.0")

    if baseline_path === nothing && corrected_path === nothing
        println("""
        Usage: julia test_planar_symmetry.jl --baseline FILE.h5 --corrected FILE.h5 \\
                   --source-lat LAT --frame N [--lat-tolerance-m T]

               julia test_planar_symmetry.jl --baseline FILE.h5 --corrected FILE.h5 \\
                   --source-lat LAT --sweep STEP [--lat-tolerance-m T]

               julia test_planar_symmetry.jl --baseline FILE.h5 --list-frames
        """)
        return
    end
    source_lat_str === nothing && error("--source-lat LAT is required (use the same latitude " *
                                          "you passed to --injection-point when running the simulation)")
    source_lat = parse(Float64, source_lat_str)
    lat_tolerance_m = parse(Float64, tol_str)

    if sweep_str !== nothing
        (baseline_path === nothing || corrected_path === nothing) &&
            error("--sweep requires both --baseline and --corrected")
        run_sweep(baseline_path, corrected_path, parse(Int, sweep_str);
                   source_lat=source_lat, lat_tolerance_m=lat_tolerance_m)
        return
    end

    frame_str === nothing && error("--frame N is required (use --list-frames to choose one, or use --sweep STEP)")
    frame_idx = parse(Int, frame_str)

    results = Dict{String,Any}()

    if baseline_path !== nothing
        r = analyse_symmetry(baseline_path, frame_idx;
                              source_lat=source_lat, lat_tolerance_m=lat_tolerance_m)
        results["baseline"] = r
        _print_result("Baseline (--gradient-correction off)", r)
    end

    if corrected_path !== nothing
        r = analyse_symmetry(corrected_path, frame_idx;
                              source_lat=source_lat, lat_tolerance_m=lat_tolerance_m)
        results["corrected"] = r
        _print_result("Corrected (--gradient-correction on)", r)
    end

    if haskey(results, "baseline") && haskey(results, "corrected")
        a, b = results["baseline"], results["corrected"]
        println("─" ^ 60)
        println("Comparison (corrected vs baseline):")
        @printf("  |asymmetry|:  %.4f → %.4f   (%s)\n", abs(a.asymmetry), abs(b.asymmetry),
                abs(b.asymmetry) < abs(a.asymmetry) ? "IMPROVED ✓" : "NOT improved ✗")
        println()
        println("Acceptance criterion (adapted for §10 Stage 3): the corrected run's " *
                 "|asymmetry| should be markedly smaller than baseline's, since a pure " *
                 "east-west slope has no physical north/south preference — any imbalance " *
                 "is a numerical artifact of mesh non-orthogonality. No mass-balance/" *
                 "stability regression (check mb_err separately from the model's own run logs).")
    end
end

# ---------------------------------------------------------------------------
# Sweep mode — track asymmetry across many frames
# ---------------------------------------------------------------------------

function run_sweep(baseline_path::String, corrected_path::String, step::Int;
                    source_lat::Float64, lat_tolerance_m::Float64)
    n_frames = h5open(baseline_path, "r") do f
        length(_frame_names(f))
    end
    println("Sweep — every $step frames, $n_frames total frames available")
    @printf("%6s  %8s | %10s %10s %9s | %10s %10s %9s | %9s\n",
            "frame", "t(s)", "base Vn", "base Vs", "base asym",
            "corr Vn", "corr Vs", "corr asym", "|asym| gain")
    for fidx in 1:step:(n_frames)
        local a, b
        try
            a = analyse_symmetry(baseline_path, fidx; source_lat=source_lat, lat_tolerance_m=lat_tolerance_m)
            b = analyse_symmetry(corrected_path, fidx; source_lat=source_lat, lat_tolerance_m=lat_tolerance_m)
        catch e
            @printf("%6d  (skipped — %s)\n", fidx, sprint(showerror, e))
            continue
        end
        gain = abs(a.asymmetry) - abs(b.asymmetry)   # positive = corrected is better
        @printf("%6d  %8.0f | %10.1f %10.1f %9.4f | %10.1f %10.1f %9.4f | %+9.4f\n",
                fidx, a.t, a.V_north, a.V_south, a.asymmetry,
                b.V_north, b.V_south, b.asymmetry, gain)
    end
    println()
    println("'|asym| gain' = |baseline asymmetry| − |corrected asymmetry|; positive means")
    println("the corrected run is closer to the physically-expected symmetric ideal at that")
    println("frame. Since volume is an integrated quantity, this should be far less noisy")
    println("frame-to-frame than the velocity-direction sweep — look for a consistently")
    println("positive gain (or at least a clear majority of frames) rather than scatter")
    println("around zero.")
end

main()
