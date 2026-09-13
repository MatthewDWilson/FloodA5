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

using HDF5
using Statistics
using Printf

const EARTH_R = 6_371_000.0   # mean Earth radius (m) — matches A5Grid._EARTH_R

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
# Geometry helpers
# ---------------------------------------------------------------------------

"""
    _project_to_local_xy(lons, lats) → (xs, ys)

Local equirectangular projection (metres) centred on the mean of the
supplied points. Consistent with the projection used throughout
A5Grid.jl/FloodModel.jl for local geometry work (_edge_geometry,
_build_wlsq_weights!).
"""
function _project_to_local_xy(lons::Vector{Float64}, lats::Vector{Float64})
    lat0 = mean(lats)
    lon0 = mean(lons)
    cos_lat0 = cosd(lat0)
    xs = deg2rad.(lons .- lon0) .* EARTH_R .* cos_lat0
    ys = deg2rad.(lats .- lat0) .* EARTH_R
    return xs, ys
end

"""
    _convex_hull_2d(pts) → Vector{Tuple{Float64,Float64}}

Andrew's monotone chain convex hull. Returns hull vertices in
counter-clockwise order, open (first point not repeated at the end).
"""
function _convex_hull_2d(pts::Vector{Tuple{Float64,Float64}})
    pts = sort(unique(pts))
    n = length(pts)
    n <= 2 && return pts

    cross(o, a, b) = (a[1]-o[1])*(b[2]-o[2]) - (a[2]-o[2])*(b[1]-o[1])

    lower = Tuple{Float64,Float64}[]
    for p in pts
        while length(lower) >= 2 && cross(lower[end-1], lower[end], p) <= 0
            pop!(lower)
        end
        push!(lower, p)
    end
    upper = Tuple{Float64,Float64}[]
    for p in reverse(pts)
        while length(upper) >= 2 && cross(upper[end-1], upper[end], p) <= 0
            pop!(upper)
        end
        push!(upper, p)
    end
    return vcat(lower[1:end-1], upper[1:end-1])
end

"""
    _polygon_area_perimeter(hull) → (area, perimeter)

Shoelace area + summed edge length on the closed hull polygon.
"""
function _polygon_area_perimeter(hull::Vector{Tuple{Float64,Float64}})
    n = length(hull)
    n < 3 && return (0.0, 0.0)
    area = 0.0
    perim = 0.0
    for k in 1:n
        x1, y1 = hull[k]
        x2, y2 = hull[mod1(k+1, n)]
        area  += x1 * y2 - x2 * y1
        perim += sqrt((x2-x1)^2 + (y2-y1)^2)
    end
    return (abs(area) / 2.0, perim)
end

polsby_popper(area::Float64, perimeter::Float64) =
    perimeter <= 0.0 ? 0.0 : (4 * π * area) / perimeter^2

"""
    _front_radius_by_direction(wet_xs, wet_ys, source_x, source_y; n_dirs=8)
        → (bearings_deg, radii_m)

For each of n_dirs equally-spaced bearing sectors (0° = north, clockwise),
the maximum distance from the source to a wet cell whose bearing falls in
that sector. Sectors with no wet cells get radius 0.0 (excluded from the
CV calculation, with a warning, since they are uninformative early in a
run before the front has reached that direction).
"""
function _front_radius_by_direction(wet_xs::Vector{Float64}, wet_ys::Vector{Float64},
                                     source_x::Float64, source_y::Float64;
                                     n_dirs::Int=8)
    radii = zeros(n_dirs)
    sector_width = 360.0 / n_dirs
    for (x, y) in zip(wet_xs, wet_ys)
        dx, dy = x - source_x, y - source_y
        r = sqrt(dx^2 + dy^2)
        r < 1e-6 && continue
        bearing = mod(atand(dx, dy), 360.0)         # 0° = north, clockwise
        sector  = Int(floor(bearing / sector_width)) + 1
        sector  = clamp(sector, 1, n_dirs)
        radii[sector] = max(radii[sector], r)
    end
    bearings = [(k - 0.5) * sector_width for k in 1:n_dirs]
    return bearings, radii
end

# ---------------------------------------------------------------------------
# HDF5 access
# ---------------------------------------------------------------------------

function _frame_names(f)
    return sort(collect(keys(f["frames"])))
end

"""
    _read_wet_cells(h5_path, frame_idx; saturation_threshold=1.0)
        → (cell_lons, cell_lats, wet_mask, t)

Reads mesh centres and the requested frame's saturation field. `frame_idx`
is the zero-based frame index as written by _write_frame! (lpad(idx,6,'0')),
NOT a 1-based Julia array index into the sorted frame-name list.
"""
function _read_wet_cells(h5_path::String, frame_idx::Int; saturation_threshold::Float64=1.0)
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
        sat = read(grp["saturation"])
        t   = haskey(grp, "t") ? read(grp["t"]) : NaN
        wet = sat .>= saturation_threshold
        return cell_lons, cell_lats, wet, t
    end
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

# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

"""
    analyse_point_spread(h5_path, frame_idx; source_lon, source_lat,
                         saturation_threshold=1.0, n_dirs=8)
        → NamedTuple

Computes Polsby–Popper circularity and directional front-radius CV for one
run at one frame.
"""
function analyse_point_spread(h5_path::String, frame_idx::Int;
                               source_lon::Float64, source_lat::Float64,
                               saturation_threshold::Float64=1.0,
                               n_dirs::Int=8)
    cell_lons, cell_lats, wet, t = _read_wet_cells(h5_path, frame_idx;
                                                    saturation_threshold=saturation_threshold)
    n_wet = count(wet)
    n_wet < 3 && error("Only $n_wet wet cell(s) at frame $frame_idx in $h5_path " *
                        "— need ≥3 for a convex hull. Pick a later frame.")

    wet_lons, wet_lats = cell_lons[wet], cell_lats[wet]
    # Project using a fixed origin (the source location) rather than the
    # wet-cell centroid, so the two runs (baseline/corrected) and all frames
    # share a common, source-anchored coordinate system — required for the
    # directional bearing calculation to mean the same thing across runs.
    lat0 = source_lat
    cos_lat0 = cosd(lat0)
    to_xy(lon, lat) = (deg2rad(lon - source_lon) * EARTH_R * cos_lat0,
                        deg2rad(lat - source_lat) * EARTH_R)

    xs = [to_xy(lo, la)[1] for (lo, la) in zip(wet_lons, wet_lats)]
    ys = [to_xy(lo, la)[2] for (lo, la) in zip(wet_lons, wet_lats)]
    src_x, src_y = 0.0, 0.0   # source is the projection origin by construction

    hull = _convex_hull_2d(collect(zip(xs, ys)))
    area, perim = _polygon_area_perimeter(hull)
    circularity = polsby_popper(area, perim)

    bearings, radii = _front_radius_by_direction(xs, ys, src_x, src_y; n_dirs=n_dirs)
    nonzero = radii .> 0.0
    n_empty_sectors = count(!, nonzero)
    cv_radius = count(nonzero) >= 2 ? std(radii[nonzero]) / mean(radii[nonzero]) : NaN

    return (t=t, n_wet=n_wet, circularity=circularity, area=area, perimeter=perim,
            cv_radius=cv_radius, n_empty_sectors=n_empty_sectors,
            bearings=bearings, radii=radii)
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
            path !== nothing && isfile(path) && _list_frames(path)
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
