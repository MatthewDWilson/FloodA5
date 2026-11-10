#!/usr/bin/env julia
# test_point_spread.jl
# ---------------------
# Point-spread directional-bias benchmark (FloodA5_NonOrthogonal_Correction_Plan.md
# §10 Stage 2 / §12 Step 12).
#
# Replicates the external R-based circularity assessment (convex hull of wet
# cell centres → Polsby–Popper circularity) as an in-repo, repeatable Julia
# test, plus an additional directional front-radius coefficient-of-variation
# metric to directly quantify the "preferential NW–SE axis" finding from the
# FOSS4G 2026 paper.
#
# This test is purely an HDF5-output analysis tool — it does NOT run the
# model itself. Generate the two comparison runs first (see Usage below),
# then point this script at both output files.
#
# Usage (from the FloodA5 project root):
#
#   Step 1 — generate the square mesh (once; reuse across both runs):
#     julia --threads 1 --project=. FloodModel.jl \
#         --meshgen test/square/square_domain.geojson --meshres 14 \
#         --meshout test/square/square_mesh_res14.parquet \
#         --flow-model standard --mesh-only
#
#   Step 2 — baseline run (uncorrected, legacy kernel):
#     julia --threads auto --project=. FloodModel.jl \
#         --meshload test/square/square_mesh_res14.parquet \
#         --flow-model standard \
#         --injection-point 0.0,0.0,50.0 \
#         --gradient-correction off \
#         --sim-duration 3600 --dt-max 30 \
#         --output test/square/square_baseline.h5 --output-interval 45
#
#   Step 3 — corrected run (WLSQ + skewness correction):
#     julia --threads auto --project=. FloodModel.jl \
#         --meshload test/square/square_mesh_res14.parquet \
#         --flow-model standard \
#         --injection-point 0.0,0.0,50.0 \
#         --gradient-correction on \
#         --sim-duration 3600 --dt-max 30 \
#         --output test/square/square_corrected.h5 --output-interval 45
#
#   Step 4 — compare:
#     julia --project=. test/test_point_spread.jl \
#         --baseline  test/square/square_baseline.h5 \
#         --corrected test/square/square_corrected.h5 \
#         --frame 80 \
#         --source-lon 0.0 --source-lat 0.0
#
# Acceptance criterion (§10 Stage 2): the corrected run should show HIGHER
# Polsby–Popper circularity and LOWER directional front-radius CV than the
# baseline, at matched output frames, with no mass-balance/stability
# regression (checked separately via the model's own mb_err logging).
#
# Pick a --frame index that is well into the spreading phase but before the
# wetted region reaches the domain boundary (boundary contact distorts the
# convex hull and is not informative about the model's intrinsic directional
# bias). Use --list-frames to inspect wet-cell-count growth and choose.
# --frame N is the frame number shown by --list-frames: the 1-based frame
# index (the same number as the old /frames/NNNNNN group name); output-format
# 1 and 2 files both work (io/FloodA5Output.jl).
#
# The metrics live in analysis/shape.jl (FloodA5Analysis, `FA.point_spread_shape`);
# this script is a CLI over them. Its conventions, passed explicitly: Earth
# radius R_EARTH = 6_371_000 m, coordinates about the NOMINAL source point
# (snap = false), bearing sectors [0°, 45°), [45°, 90°), … (centred = false; see
# test_octagonal_anisotropy.jl for sectors centred on the compass bearings).

using Printf
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis

# ---------------------------------------------------------------------------
# CLI parsing (minimal, no external deps)
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
# Core analysis
# ---------------------------------------------------------------------------

"""
    analyse_point_spread(h5_path, frame_idx; source_lon, source_lat,
                         saturation_threshold=1.0, n_dirs=8)
        → NamedTuple

Computes Polsby–Popper circularity and directional front-radius CV for one
run at one frame (`frame_idx` is the 1-based frame index --list-frames prints).
"""
function analyse_point_spread(h5_path::String, frame_idx::Int;
                               source_lon::Float64, source_lat::Float64,
                               saturation_threshold::Float64=1.0,
                               n_dirs::Int=8)
    r = FA.load_run(h5_path)
    FA.check_frame(r, frame_idx; hint = " Run with --list-frames to inspect.")
    # Project about the source location (not the wet-cell centroid), so the
    # two runs (baseline/corrected) and all frames share a common,
    # source-anchored coordinate system for the directional bearings.
    res = FA.point_spread_shape(r, frame_idx; lon0 = source_lon, lat0 = source_lat,
                                R = FA.R_EARTH, snap = false, centred = false,
                                saturation_threshold, n_dirs)
    res.n_wet < 3 && error("Only $(res.n_wet) wet cell(s) at frame $frame_idx in $h5_path " *
                           "— need ≥3 for a convex hull. Pick a later frame.")
    return res
end

function _print_result(label::String, r)
    println(label, ":")
    @printf("  t = %.1f s   n_wet = %d\n", r.t, r.n_wet)
    @printf("  Polsby–Popper circularity = %.4f   (1.0 = perfect circle)\n", r.circularity)
    if isnan(r.cv_radius)
        println("  cv_radius = NaN  (fewer than 2 non-empty bearing sectors — front too small/early)")
    else
        @printf("  directional front-radius CV = %.4f   (lower = more isotropic)\n", r.cv_radius)
    end
    if r.n_empty_sectors > 0
        @printf("  ⚠ %d/%d bearing sectors had no wet cells — front has not yet reached all directions; consider a later frame for a more reliable CV.\n",
                r.n_empty_sectors, length(r.bearings))
    end
    println("  per-sector radii (m):")
    for (b, rad) in zip(r.bearings, r.radii)
        @printf("    bearing %5.1f°: %8.1f m\n", b, rad)
    end
    println()
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    if DO_LIST
        for flag in ("--baseline", "--corrected")
            path = _argval(ARGS_V, flag)
            path !== nothing && isfile(path) && FA.print_frames(FA.load_run(path))
        end
        return
    end

    baseline_path  = _argval(ARGS_V, "--baseline")
    corrected_path = _argval(ARGS_V, "--corrected")
    frame_str      = _argval(ARGS_V, "--frame")
    src_lon_str    = _argval(ARGS_V, "--source-lon", "0.0")
    src_lat_str    = _argval(ARGS_V, "--source-lat", "0.0")
    sat_thresh_str = _argval(ARGS_V, "--saturation-threshold", "1.0")

    if baseline_path === nothing && corrected_path === nothing
        println("""
        Usage: julia test_point_spread.jl --baseline FILE.h5 --corrected FILE.h5 \\
                   --frame N [--source-lon LON] [--source-lat LAT] \\
                   [--saturation-threshold T]

               julia test_point_spread.jl --baseline FILE.h5 --list-frames
               (inspect frame indices / wet-cell growth before picking --frame)
        """)
        return
    end
    frame_str === nothing && error("--frame N is required (use --list-frames to choose one)")

    frame_idx = parse(Int, frame_str)
    source_lon = parse(Float64, src_lon_str)
    source_lat = parse(Float64, src_lat_str)
    sat_thresh = parse(Float64, sat_thresh_str)

    results = Dict{String,Any}()

    if baseline_path !== nothing
        r = analyse_point_spread(baseline_path, frame_idx;
                                  source_lon=source_lon, source_lat=source_lat,
                                  saturation_threshold=sat_thresh)
        results["baseline"] = r
        _print_result("Baseline (--gradient-correction off)", r)
    end

    if corrected_path !== nothing
        r = analyse_point_spread(corrected_path, frame_idx;
                                  source_lon=source_lon, source_lat=source_lat,
                                  saturation_threshold=sat_thresh)
        results["corrected"] = r
        _print_result("Corrected (--gradient-correction on)", r)
    end

    if haskey(results, "baseline") && haskey(results, "corrected")
        a, b = results["baseline"], results["corrected"]
        println("─" ^ 60)
        println("Comparison (corrected vs baseline):")
        @printf("  circularity:  %.4f → %.4f   (%s)\n", a.circularity, b.circularity,
                b.circularity > a.circularity ? "IMPROVED ✓" : "NOT improved ✗")
        if !isnan(a.cv_radius) && !isnan(b.cv_radius)
            @printf("  cv_radius:    %.4f → %.4f   (%s)\n", a.cv_radius, b.cv_radius,
                    b.cv_radius < a.cv_radius ? "IMPROVED ✓" : "NOT improved ✗")
        else
            println("  cv_radius:    not comparable (NaN in one or both runs — " *
                     "front too small/early; try a later --frame)")
        end
        println()
        println("Acceptance criterion (§10 Stage 2): corrected circularity > baseline " *
                 "AND corrected cv_radius < baseline, with no mass-balance/stability " *
                 "regression (check mb_err separately from the model's own run logs).")
    end
end

main()
