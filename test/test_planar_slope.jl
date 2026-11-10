#!/usr/bin/env julia
# test_planar_slope.jl
# ---------------------
# Planar slope flow-direction benchmark (FloodA5_NonOrthogonal_Correction_Plan.md
# §10 Stage 3 / §12 Step 14).
#
# Measures the mean flow direction (relative to due east, the analytically
# correct downslope direction on this DEM) and its spread across wet cells,
# for a baseline (--gradient-correction off) and corrected (on) run on the
# same pure east-west planar slope DEM (embankment height 0 — see
# generate_planar_embankment_dem.py --emb-height 0).
#
# Unlike test_point_spread.jl (which infers direction from the SHAPE of the
# wetted region), this test reads the per-cell velocity components (vel_u,
# vel_v) directly — a much more direct measure of flow direction, and the
# reason vel_u/vel_v were added to _write_frame!'s HDF5 output.
#
# This test does NOT run the model. Generate the two comparison runs first
# (see the CLI commands accompanying this file), then point this script at
# both output files.
#
# Usage:
#   julia --project=. test/test_planar_slope.jl \
#       --baseline  test/planar_embankment/planar_baseline.h5 \
#       --corrected test/planar_embankment/planar_corrected.h5 \
#       --frame 30
#
#   julia --project=. test/test_planar_slope.jl \
#       --baseline test/planar_embankment/planar_baseline.h5 --list-frames
#
# --frame N is the 1-based frame index (the same number as the old
# /frames/NNNNNN group name); output-format 1 and 2 files both work
# (io/FloodA5Output.jl).
#
# Direction convention: angle = atand(vel_v, vel_u) — standard mathematical
# convention, NOT the compass-bearing convention used in test_point_spread.jl.
#   0°    = due east   (the analytically correct downslope direction here)
#   +90°  = due north
#   -90°  = due south
#   ±180° = due west   (upslope — should not occur in a stable run)
#
# Acceptance criterion (§10 Stage 3): the corrected run's circular mean
# direction should be within ±10° of due east (0°); the baseline is expected
# to reproduce the paper's reported 30–60° deviation. The corrected run's
# circular standard deviation should also be markedly lower (tighter
# clustering around the true downslope direction).
#
# The statistics live in analysis/shape.jl (FloodA5Analysis:
# `FA.flow_direction_metrics`, `FA.circular_mean_std`); this script is a CLI
# over them.

using Printf
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis

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
# Core analysis
# ---------------------------------------------------------------------------

"""
    analyse_planar_slope(run_or_path, frame_idx; depth_threshold, speed_threshold)

Circular statistics (unweighted and speed-weighted; `FA.direction_stats`) of
the velocity direction atand(vel_v, vel_u) over the cells that are both wet
(water_depth > depth_threshold) and moving (speed > speed_threshold) — a wet
but near-stationary cell has an ill-defined direction that adds noise, not
signal (`FA.flow_directions`).
"""
function analyse_planar_slope(r::FA.Run, frame_idx::Int;
                               depth_threshold::Float64=1e-3,
                               speed_threshold::Float64=1e-4)
    FA.check_frame(r, frame_idx; hint = " Run with --list-frames to inspect.")
    FA.has_field(r, "vel_u") || error(
        "Frame $(lpad(string(frame_idx), 6, '0')) in $(r.path) has no vel_u/vel_v dataset. This file " *
        "was written before vel_u/vel_v were added to _write_frame! — " *
        "re-run the simulation with the updated FloodModel.jl.")
    res = FA.flow_direction_metrics(r, frame_idx; depth_threshold, speed_threshold)
    res.n_moving < 3 && error("Only $(res.n_moving) moving cell(s) at frame $frame_idx in " *
                           "$(r.path) — need a meaningful sample. Pick a later frame " *
                           "(more time for the front to develop) or lower " *
                           "--speed-threshold.")
    return res
end
analyse_planar_slope(h5_path::String, frame_idx::Int; kw...) =
    analyse_planar_slope(FA.load_run(h5_path), frame_idx; kw...)

function _print_result(label::String, r)
    println(label, ":")
    @printf("  t = %.1f s   n_wet = %d   n_moving = %d\n", r.t, r.n_wet, r.n_moving)
    println("  ── unweighted (every moving cell counts equally) ──")
    @printf("    circular mean direction = %+.1f°  (0° = due east)\n", r.mean_dir)
    if isinf(r.circ_std)
        println("    circular std = Inf  (R≈0 — directions essentially uniform/bimodal around the circle; mean is not meaningful)")
    else
        @printf("    circular std            = %.1f°\n", r.circ_std)
    end
    @printf("    R (alignment strength)  = %.4f   (1.0 = perfectly aligned, 0.0 = no preferred direction)\n", r.R)
    println("  ── speed-weighted (cells with real momentum dominate) ──")
    @printf("    circular mean direction = %+.1f°  (0° = due east, the correct downslope direction)\n", r.mean_dir_w)
    if isinf(r.circ_std_w)
        println("    circular std = Inf  (R≈0 — even speed-weighted, no preferred direction)")
    else
        @printf("    circular std            = %.1f°   (lower = tighter clustering)\n", r.circ_std_w)
    end
    @printf("    R (alignment strength)  = %.4f\n", r.R_w)
    @printf("  mean |deviation| from east (naive, unweighted) = %.1f°   (≈90° is the value expected from purely random directions — a useful red flag, not a target)\n",
            r.mean_abs_dev)
    @printf("  fraction of moving cells within ±30° of east (unweighted) = %.1f%%\n", 100 * r.frac_within_30)
    if r.frac_within_30 < 0.5 && abs(r.mean_dir) < 15.0
        println("  ⚠ low unweighted circular std but most cells are NOT within ±30° of east — possible bimodal split rather than genuine eastward alignment. Trust the speed-weighted R above this unweighted mean.")
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
    sweep_str      = _argval(ARGS_V, "--sweep")
    depth_thresh_s = _argval(ARGS_V, "--depth-threshold", "1e-3")
    speed_thresh_s = _argval(ARGS_V, "--speed-threshold", "1e-4")
    depth_threshold = parse(Float64, depth_thresh_s)
    speed_threshold = parse(Float64, speed_thresh_s)

    if baseline_path === nothing && corrected_path === nothing
        println("""
        Usage: julia test_planar_slope.jl --baseline FILE.h5 --corrected FILE.h5 \\
                   --frame N [--depth-threshold T] [--speed-threshold T]

               julia test_planar_slope.jl --baseline FILE.h5 --corrected FILE.h5 \\
                   --sweep STEP [--depth-threshold T] [--speed-threshold T]
               (print the comparison every STEP frames across the whole run,
                instead of a single --frame snapshot)

               julia test_planar_slope.jl --baseline FILE.h5 --list-frames
               (inspect frame indices / wet-cell growth before picking --frame)
        """)
        return
    end

    if sweep_str !== nothing
        (baseline_path === nothing || corrected_path === nothing) &&
            error("--sweep requires both --baseline and --corrected")
        run_sweep(baseline_path, corrected_path, parse(Int, sweep_str);
                   depth_threshold=depth_threshold, speed_threshold=speed_threshold)
        return
    end

    frame_str === nothing && error("--frame N is required (use --list-frames to choose one, or use --sweep STEP)")

    frame_idx = parse(Int, frame_str)

    results = Dict{String,Any}()

    if baseline_path !== nothing
        r = analyse_planar_slope(baseline_path, frame_idx;
                                  depth_threshold=depth_threshold,
                                  speed_threshold=speed_threshold)
        results["baseline"] = r
        _print_result("Baseline (--gradient-correction off)", r)
    end

    if corrected_path !== nothing
        r = analyse_planar_slope(corrected_path, frame_idx;
                                  depth_threshold=depth_threshold,
                                  speed_threshold=speed_threshold)
        results["corrected"] = r
        _print_result("Corrected (--gradient-correction on)", r)
    end

    if haskey(results, "baseline") && haskey(results, "corrected")
        a, b = results["baseline"], results["corrected"]
        println("─" ^ 60)
        println("Comparison (corrected vs baseline) — using speed-weighted stats (more trustworthy than unweighted when R is low):")
        @printf("  |circular mean dev. from east| (weighted):  %.1f° → %.1f°   (%s)\n",
                abs(a.mean_dir_w), abs(b.mean_dir_w),
                abs(b.mean_dir_w) < abs(a.mean_dir_w) ? "IMPROVED ✓" : "NOT improved ✗")
        if !isinf(a.circ_std_w) && !isinf(b.circ_std_w)
            @printf("  circular std (weighted):                     %.1f° → %.1f°   (%s)\n",
                    a.circ_std_w, b.circ_std_w,
                    b.circ_std_w < a.circ_std_w ? "IMPROVED ✓" : "NOT improved ✗")
        else
            println("  circular std (weighted): not comparable (Inf in one or both runs)")
        end
        @printf("  R, alignment strength (weighted):            %.4f → %.4f\n", a.R_w, b.R_w)
        if max(a.R_w, b.R_w) < 0.2
            println("  ⚠ R is low (<0.2) in both runs even after speed-weighting — the directional signal at this single frame is weak. Treat the mean-direction comparison above with caution; consider --sweep to see whether R strengthens over time before drawing a conclusion.")
        end
        println()
        println("Acceptance criterion (§10 Stage 3): corrected circular mean direction " *
                 "should be within ±10° of due east; baseline is expected to reproduce " *
                 "the paper's reported 30–60° deviation. No mass-balance/stability " *
                 "regression (check mb_err separately from the model's own run logs).")
    end
end

# ---------------------------------------------------------------------------
# Sweep mode — track the comparison across many frames
# ---------------------------------------------------------------------------

"""
Runs the analysis every `step` frames from 0 up to the last available frame
in the baseline file (both files are expected to have the same frame
count/timing — generated from the same --sim-duration/--output-interval).
Prints a compact table so you can see whether the directional signal (R)
strengthens as the flood matures, rather than trusting a single snapshot
that might happen to fall at a low-R moment.
"""
function run_sweep(baseline_path::String, corrected_path::String, step::Int;
                    depth_threshold::Float64, speed_threshold::Float64)
    rb, rc = FA.load_run(baseline_path), FA.load_run(corrected_path)
    n_frames = FA.nframes(rb)
    println("Sweep — every $step frames, $n_frames total frames available")
    @printf("%6s  %8s | %8s %8s %7s | %8s %8s %7s | %8s\n",
            "frame", "t(s)", "base θw", "base σw", "base Rw",
            "corr θw", "corr σw", "corr Rw", "Rw gain")
    for fidx in 0:step:(n_frames - 1)
        local a, b
        try
            a = analyse_planar_slope(rb, fidx;
                                      depth_threshold=depth_threshold, speed_threshold=speed_threshold)
            b = analyse_planar_slope(rc, fidx;
                                      depth_threshold=depth_threshold, speed_threshold=speed_threshold)
        catch e
            @printf("%6d  (skipped — %s)\n", fidx, sprint(showerror, e))
            continue
        end
        rw_gain = b.R_w - a.R_w
        @printf("%6d  %8.0f | %8.1f %8.1f %7.4f | %8.1f %8.1f %7.4f | %+8.4f\n",
                fidx, a.t, a.mean_dir_w, a.circ_std_w, a.R_w,
                b.mean_dir_w, b.circ_std_w, b.R_w, rw_gain)
    end
    println()
    println("θw/σw/Rw = speed-weighted circular mean direction / std / alignment strength.")
    println("Watch Rw over time: if it trends upward in both runs as the flood matures, ")
    println("low-R readings early on reflect a genuinely weak/noisy signal at that time, ")
    println("not a measurement problem. If corrected consistently shows higher Rw and ")
    println("lower |θw| than baseline across most frames (not just one cherry-picked ")
    println("frame), that is much stronger evidence of a real improvement.")
end

main()
