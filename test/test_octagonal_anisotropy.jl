#!/usr/bin/env julia
# test_octagonal_anisotropy.jl
# -----------------------------
# Cardinal-vs-diagonal directional bias check for 4-connected (Cartesian)
# grids — a companion to test_point_spread.jl, not a replacement for it.
#
# Motivation (from design discussion, cartesian-formulation branch):
# test_point_spread.jl's 8-sector front-radius metric bins wet cells into
# sectors CENTRED AT 22.5°, 67.5°, 112.5°, ... — i.e. each sector straddles
# the midpoint between a cardinal (N/E/S/W) and a diagonal (NE/SE/SW/NW)
# bearing, rather than being centred ON those bearings. That's a reasonable
# default for a bias with no particular preferred axis (A5's own pentagon
# chirality, say), but it can't cleanly test the specific hypothesis that
# matters for a 4-connected grid: does flow reach cardinal directions
# faster than diagonal ones, because a diagonal path through a 4-connected
# grid is a zig-zag ~1.41x longer than the straight-line distance, while a
# cardinal path is exact? (See FloodA5_PhaseA_Task3_Findings.md §6 for the
# equivalent analysis on the A5 side — same idea, different grid.)
#
# This script provides two things test_point_spread.jl doesn't:
#
#   1. Sectors CENTRED exactly on the 8 compass bearings
#      (0°/45°/90°/.../315° = N/NE/E/SE/S/SW/W/NW), so "cardinal" and
#      "diagonal" are literal sector labels, not an average across a
#      boundary. Reports the cardinal-mean-radius / diagonal-mean-radius
#      ratio directly — the number this whole exercise is about.
#
#   2. A fixed-radius, exact-bearing sample: rather than "furthest wet
#      cell in a sector" (which can be dominated by whichever one or two
#      discrete cells happen to be furthest — see the worked caution
#      below), pick ONE shared target radius and ask "at exactly this
#      distance, in exactly this direction, what does the mesh actually
#      look like" (wet/dry, depth) for each of the 8 bearings. This is
#      closer to the "sample points at equal distance in all 8 directions"
#      framing from the design discussion, and is less sensitive to a
#      single outlier cell than a per-sector maximum.
#
# CAUTION on interpreting either metric at low wet-cell counts: with only
# a few hundred–thousand wet cells split across 8 (or 16) bins, a sector's
# reported radius is often set by literally one or two specific cells, not
# a smooth front. Two adjacent sectors showing IDENTICAL radii (as seen in
# an early run of this branch's own planar test) is a tell — it usually
# means the same one or two near-diagonal-neighbour cells are the maximum
# for both bins, not a real, independent measurement in each direction.
# Prefer a later frame (more wet cells) and/or --n-dirs 16 before treating
# a difference between opposite sectors as a confirmed, systematic bias
# rather than discretisation noise. This script prints the raw wet-cell
# count per sector alongside the radius specifically so this can be judged
# rather than assumed either way.
#
# Usage:
#   julia --project=. test/test_octagonal_anisotropy.jl \
#       --input FILE.h5 --frame N --source-lon LON --source-lat LAT \
#       [--n-dirs 8|16] [--radius-m R] [--saturation-threshold T]
#
#   julia --project=. test/test_octagonal_anisotropy.jl --input FILE.h5 --list-frames
#
#   N is the 1-based frame index (the same number as the old /frames/NNNNNN
#   group name); output-format 1 and 2 files both work (io/FloodA5Output.jl).
#
# The metrics live in analysis/shape.jl (FloodA5Analysis:
# `FA.front_radius_by_sector(...; centred = true)`, `FA.sample_at_bearings`);
# this script is a CLI over them. Coordinates: equirectangular about the
# NOMINAL source point with R_EARTH = 6_371_000 m (as test_point_spread.jl).

using Statistics
using Printf
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis

function _argval(args::Vector{String}, flag::String, default=nothing)
    idx = findfirst(==(flag), args)
    idx === nothing && return default
    idx == length(args) && error("$flag requires a value")
    return args[idx + 1]
end

const ARGS_V = String.(ARGS)
const DO_LIST = "--list-frames" in ARGS_V

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    input_path = _argval(ARGS_V, "--input")
    input_path === nothing && error("--input FILE.h5 is required")

    if DO_LIST
        FA.print_frames(FA.load_run(input_path))
        return
    end

    frame_str   = _argval(ARGS_V, "--frame")
    frame_str === nothing && error("--frame N is required (use --list-frames to choose one)")
    frame_idx   = parse(Int, frame_str)
    source_lon  = parse(Float64, _argval(ARGS_V, "--source-lon", "0.0"))
    source_lat  = parse(Float64, _argval(ARGS_V, "--source-lat", "0.0"))
    n_dirs      = parse(Int, _argval(ARGS_V, "--n-dirs", "8"))
    sat_thresh  = parse(Float64, _argval(ARGS_V, "--saturation-threshold", "1.0"))
    radius_str  = _argval(ARGS_V, "--radius-m")

    n_dirs in (8, 16) || error("--n-dirs must be 8 or 16 (got $n_dirs)")

    run = FA.load_run(input_path)
    FA.check_frame(run, frame_idx; hint = " Run with --list-frames to inspect.")
    sat   = FA.field(run, "saturation", frame_idx)
    depth = FA.has_field(run, "water_depth") ? FA.depth(run, frame_idx) : fill(NaN, length(sat))
    t     = run.t[frame_idx]
    wet = sat .>= sat_thresh
    n_wet = count(wet)
    n_wet < 3 && error("Only $n_wet wet cell(s) at frame $frame_idx — pick a later frame.")

    all_xs, all_ys = FA.shape_xy(run.mesh, source_lon, source_lat; R = FA.R_EARTH, snap = false)
    wet_xs = all_xs[wet];    wet_ys = all_ys[wet]

    println("Input: $input_path   frame=$frame_idx   t=$(round(t,digits=1))s   n_wet=$n_wet")
    println()

    # ── Metric 1: centred-sector front radius ──────────────────────────────
    println("="^70)
    println("Metric 1: front radius, sectors centred on compass bearings")
    println("="^70)
    bearings, radii, counts = FA.front_radius_by_sector(wet_xs, wet_ys; n_dirs, centred = true)

    labels8  = ["N","NE","E","SE","S","SW","W","NW"]
    labels16 = ["N","NNE","NE","ENE","E","ESE","SE","SSE",
                "S","SSW","SW","WSW","W","WNW","NW","NNW"]
    labels = n_dirs == 8 ? labels8 : labels16
    # NOTE: at --n-dirs 16, "cardinal" still means exactly N/E/S/W, but
    # "diagonal" below then also sweeps in the 8 intermediate bearings
    # (NNE, ENE, ...) alongside the 4 true diagonals — the summary ratio at
    # n_dirs=16 is "cardinal vs everything else", not a clean 4-vs-4 split.
    # Use n_dirs=8 for the primary cardinal-vs-diagonal ratio; use n_dirs=16
    # only to eyeball the full per-bearing table for discretisation noise.
    is_cardinal = FA.cardinal_mask(bearings)

    @printf("%6s  %-4s  %10s  %8s\n", "bear°", "dir", "radius(m)", "n_cells")
    for (b, r, c, lab) in zip(bearings, radii, counts, labels)
        flag = c <= 2 ? "  ⚠ n<=2 — see caution in file header" : ""
        @printf("%6.1f  %-4s  %10.1f  %8d%s\n", b, lab, r, c, flag)
    end

    cdm = FA.cardinal_diagonal_means(bearings, radii)
    println()
    if cdm.n_cardinal >= 2 && cdm.n_diagonal >= 2
        @printf("Mean cardinal radius (N/E/S/W)   : %8.1f m  (n=%d)\n", cdm.cardinal, cdm.n_cardinal)
        @printf("Mean diagonal radius (NE/SE/SW/NW): %8.1f m  (n=%d)\n", cdm.diagonal, cdm.n_diagonal)
        @printf("Cardinal / diagonal ratio          : %6.4f\n", cdm.ratio)
        println("  (1.0 = no directional preference; >1 favours cardinal, matching " *
                 "the 4-connected zig-zag-path hypothesis)")
    else
        println("Not enough non-empty cardinal/diagonal sectors to compute a ratio " *
                "(front hasn't reached enough directions yet — try a later frame).")
    end

    # ── Metric 2: fixed-radius exact-bearing sample ────────────────────────
    println()
    println("="^70)
    println("Metric 2: fixed-radius sample at exact compass bearings")
    println("="^70)

    radius_m = if radius_str !== nothing
        parse(Float64, radius_str)
    else
        # Auto: 70% of the median non-zero centred-sector radius from Metric 1
        # — comfortably inside the front in most directions without needing
        # the user to guess a value.
        nz = radii[radii .> 0]
        isempty(nz) && error("No wet cells found in any sector — cannot auto-pick a radius; " *
                              "pass --radius-m explicitly.")
        0.7 * median(nz)
    end
    @printf("Sample radius: %.1f m %s\n", radius_m,
            radius_str === nothing ? "(auto: 70%% of median sector radius)" : "(user-specified)")
    println()

    samples = [(bearing = s.bearing, dist_to_target = s.dist,
                wet = sat[s.cell] >= 1.0, depth = depth[s.cell])
               for s in FA.sample_at_bearings(all_xs, all_ys, radius_m; n_dirs)]
    @printf("%6s  %-4s  %8s  %10s  %10s\n", "bear°", "dir", "wet?", "depth(m)", "cell-dist(m)")
    for (s, lab) in zip(samples, labels)
        @printf("%6.1f  %-4s  %8s  %10.4f  %10.1f\n",
                s.bearing, lab, s.wet ? "wet" : "dry",
                s.wet ? s.depth : NaN, s.dist_to_target)
    end

    card_depths = [s.depth for (s, c) in zip(samples, is_cardinal) if c && s.wet]
    diag_depths = [s.depth for (s, c) in zip(samples, is_cardinal) if !c && s.wet]
    println()
    if length(card_depths) >= 1 && length(diag_depths) >= 1
        @printf("Mean cardinal depth at r=%.0fm: %.4f m (n=%d)\n", radius_m, mean(card_depths), length(card_depths))
        @printf("Mean diagonal depth at r=%.0fm: %.4f m (n=%d)\n", radius_m, mean(diag_depths), length(diag_depths))
    else
        println("Not enough wet cardinal/diagonal samples at this radius to compare depths " *
                "(try a smaller --radius-m, well inside the current front).")
    end
end

main()
