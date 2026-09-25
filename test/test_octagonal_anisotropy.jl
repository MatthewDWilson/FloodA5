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

using HDF5
using Statistics
using Printf

const EARTH_R = 6_371_000.0

function _argval(args::Vector{String}, flag::String, default=nothing)
    idx = findfirst(==(flag), args)
    idx === nothing && return default
    idx == length(args) && error("$flag requires a value")
    return args[idx + 1]
end

const ARGS_V = String.(ARGS)
const DO_LIST = "--list-frames" in ARGS_V

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
            depth = haskey(grp, "water_depth") ? read(grp["water_depth"]) : Float64[]
            md = isempty(depth) ? NaN : maximum(depth)
            @printf("%8s  %10.1f  %8d  %10.4f\n", nm, t, n_wet, md)
        end
    end
end

"""
    _read_frame(h5_path, frame_idx) → (cell_lons, cell_lats, sat, depth, t)
"""
function _read_frame(h5_path::String, frame_idx::Int)
    h5open(h5_path, "r") do f
        cell_lons = read(f["mesh/center_lons"])
        cell_lats = read(f["mesh/center_lats"])
        names = _frame_names(f)
        fname = lpad(string(frame_idx), 6, '0')
        haskey(f["frames"], fname) ||
            error("Frame $fname not found in $h5_path. Available: " *
                  "$(first(names, 3))…$(last(names, 3)) ($(length(names)) total). " *
                  "Run with --list-frames to inspect.")
        grp = f["frames/$fname"]
        sat   = read(grp["saturation"])
        depth = haskey(grp, "water_depth") ? read(grp["water_depth"]) : fill(NaN, length(sat))
        t     = haskey(grp, "t") ? read(grp["t"]) : NaN
        return cell_lons, cell_lats, sat, depth, t
    end
end

# ---------------------------------------------------------------------------
# Metric 1 — sectors centred exactly on the 8 (or 16) compass bearings
# ---------------------------------------------------------------------------

"""
    _front_radius_centred(xs, ys, n_dirs) → (bearings, radii, counts)

Like test_point_spread.jl's _front_radius_by_direction, but sectors are
centred ON each compass bearing (bearing_k ± half-width), not offset by
half a sector. `counts[k]` is the number of wet cells that fell in sector
k — print this alongside the radius; a sector controlled by 1-2 cells is
not a reliable per-direction measurement (see module docstring).
"""
function _front_radius_centred(xs::Vector{Float64}, ys::Vector{Float64}; n_dirs::Int=8)
    radii  = zeros(n_dirs)
    counts = zeros(Int, n_dirs)
    sector_width = 360.0 / n_dirs
    half = sector_width / 2
    for (x, y) in zip(xs, ys)
        r = sqrt(x^2 + y^2)
        r < 1e-6 && continue
        bearing = mod(atand(x, y), 360.0)               # 0°=N, clockwise
        shifted = mod(bearing + half, 360.0)
        sector  = clamp(Int(floor(shifted / sector_width)) + 1, 1, n_dirs)
        counts[sector] += 1
        radii[sector] = max(radii[sector], r)
    end
    bearings = [(k - 1) * sector_width for k in 1:n_dirs]   # 0, 45, 90, ...
    return bearings, radii, counts
end

# ---------------------------------------------------------------------------
# Metric 2 — fixed-radius, exact-bearing point sample
# ---------------------------------------------------------------------------

"""
    _sample_at_bearings(cell_xs, cell_ys, sat, depth, radius_m, n_dirs)
        → per-bearing (bearing, nearest_dist_to_target, wet, depth)

For each of n_dirs exact compass bearings, computes the target point at
(radius_m, bearing) from the source (origin, by construction of the
caller's projection) and finds the NEAREST actual cell centre to that
point — over ALL cells, not just wet ones, since we need to know whether
the mesh is wet or dry there. `nearest_dist_to_target` is reported so you
can sanity-check resolution (should be well under one cell diagonal;
markedly larger means the target point fell outside the loaded cell set,
e.g. radius_m too large for this mesh).
"""
function _sample_at_bearings(cell_xs::Vector{Float64}, cell_ys::Vector{Float64},
                              sat::Vector{Float64}, depth::Vector{Float64},
                              radius_m::Float64, n_dirs::Int)
    bearings = [(k - 1) * (360.0 / n_dirs) for k in 1:n_dirs]
    results = NamedTuple[]
    for b in bearings
        tx = radius_m * sind(b)
        ty = radius_m * cosd(b)
        best_i = 0
        best_d2 = Inf
        for i in eachindex(cell_xs)
            d2 = (cell_xs[i] - tx)^2 + (cell_ys[i] - ty)^2
            if d2 < best_d2
                best_d2 = d2
                best_i = i
            end
        end
        push!(results, (bearing=b, dist_to_target=sqrt(best_d2),
                         wet=sat[best_i] >= 1.0, depth=depth[best_i]))
    end
    return results
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    input_path = _argval(ARGS_V, "--input")
    input_path === nothing && error("--input FILE.h5 is required")

    if DO_LIST
        _list_frames(input_path)
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

    cell_lons, cell_lats, sat, depth, t = _read_frame(input_path, frame_idx)
    wet = sat .>= sat_thresh
    n_wet = count(wet)
    n_wet < 3 && error("Only $n_wet wet cell(s) at frame $frame_idx — pick a later frame.")

    cos_lat0 = cosd(source_lat)
    to_xy(lon, lat) = (deg2rad(lon - source_lon) * EARTH_R * cos_lat0,
                        deg2rad(lat - source_lat) * EARTH_R)
    all_xy = [to_xy(lo, la) for (lo, la) in zip(cell_lons, cell_lats)]
    all_xs = first.(all_xy); all_ys = last.(all_xy)
    wet_xs = all_xs[wet];    wet_ys = all_ys[wet]

    println("Input: $input_path   frame=$frame_idx   t=$(round(t,digits=1))s   n_wet=$n_wet")
    println()

    # ── Metric 1: centred-sector front radius ──────────────────────────────
    println("="^70)
    println("Metric 1: front radius, sectors centred on compass bearings")
    println("="^70)
    bearings, radii, counts = _front_radius_centred(wet_xs, wet_ys; n_dirs=n_dirs)

    labels8  = ["N","NE","E","SE","S","SW","W","NW"]
    labels16 = ["N","NNE","NE","ENE","E","ESE","SE","SSE",
                "S","SSW","SW","WSW","W","WNW","NW","NNW"]
    labels = n_dirs == 8 ? labels8 : labels16

    is_cardinal = n_dirs == 8 ?
        [true,false,true,false,true,false,true,false] :        # N,E,S,W
        [b % 90 == 0 for b in bearings]
    # NOTE: at --n-dirs 16, "cardinal" still means exactly N/E/S/W, but
    # "diagonal" below then also sweeps in the 8 intermediate bearings
    # (NNE, ENE, ...) alongside the 4 true diagonals — the summary ratio at
    # n_dirs=16 is "cardinal vs everything else", not a clean 4-vs-4 split.
    # Use n_dirs=8 for the primary cardinal-vs-diagonal ratio; use n_dirs=16
    # only to eyeball the full per-bearing table for discretisation noise.

    @printf("%6s  %-4s  %10s  %8s\n", "bear°", "dir", "radius(m)", "n_cells")
    for (b, r, c, lab) in zip(bearings, radii, counts, labels)
        flag = c <= 2 ? "  ⚠ n<=2 — see caution in file header" : ""
        @printf("%6.1f  %-4s  %10.1f  %8d%s\n", b, lab, r, c, flag)
    end

    cardinal_r = radii[is_cardinal]
    diagonal_r = radii[.!is_cardinal]
    cardinal_valid = cardinal_r[cardinal_r .> 0]
    diagonal_valid = diagonal_r[diagonal_r .> 0]

    println()
    if length(cardinal_valid) >= 2 && length(diagonal_valid) >= 2
        mc, md = mean(cardinal_valid), mean(diagonal_valid)
        @printf("Mean cardinal radius (N/E/S/W)   : %8.1f m  (n=%d)\n", mc, length(cardinal_valid))
        @printf("Mean diagonal radius (NE/SE/SW/NW): %8.1f m  (n=%d)\n", md, length(diagonal_valid))
        @printf("Cardinal / diagonal ratio          : %6.4f\n", mc / md)
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

    samples = _sample_at_bearings(all_xs, all_ys, sat, depth, radius_m, n_dirs)
    @printf("%6s  %-4s  %8s  %10s  %10s\n", "bear°", "dir", "wet?", "depth(m)", "cell-dist(m)")
    for (s, lab) in zip(samples, labels)
        @printf("%6.1f  %-4s  %8s  %10.4f  %10.1f\n",
                s.bearing, lab, s.wet ? "wet" : "dry",
                s.wet ? s.depth : NaN, s.dist_to_target)
    end

    is_cardinal2 = n_dirs == 8 ?
        [true,false,true,false,true,false,true,false] :
        [b % 90 == 0 for b in bearings]
    card_depths = [s.depth for (s, c) in zip(samples, is_cardinal2) if c && s.wet]
    diag_depths = [s.depth for (s, c) in zip(samples, is_cardinal2) if !c && s.wet]
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
