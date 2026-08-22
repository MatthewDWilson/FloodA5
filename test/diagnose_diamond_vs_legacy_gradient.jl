#!/usr/bin/env julia
# test/diagnose_diamond_vs_legacy_gradient.jl
#
# Same-snapshot comparison of the legacy WLSQ+skew dWSE_n construction vs
# the new B2/B3 diamond dWSE_n construction, evaluated on a REAL WSE field
# taken from an existing simulation output — no new simulation is run.
#
# Purpose: test the "diamond is less-smoothed / less-damped than the
# cell-averaged legacy construction, and that's why it flips sign faster
# and harder near the front" hypothesis directly, by comparing the two
# constructions' actual numerical output edge-by-edge on the same field,
# rather than inferring it indirectly from two different simulations'
# outcomes.
#
# Usage:
#   julia --project=. test/diagnose_diamond_vs_legacy_gradient.jl \
#       <mesh.parquet> <output.h5> [frame|last] [source_lat source_lon]
#
# `mesh.parquet` MUST be the exact mesh file used to produce `output.h5`
# (the script asserts /mesh/cell_ids matches the loaded mesh's cell order
# and errors out clearly if not, rather than silently misaligning arrays).
#
# `source_lat`/`source_lon` are optional; if omitted, the cell with
# maximum water depth at the chosen frame is used as a proxy source
# location for the distance-from-source binning.

using Statistics, HDF5, LinearAlgebra, Printf

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

length(ARGS) >= 2 || error("Usage: julia diagnose_diamond_vs_legacy_gradient.jl " *
                            "<mesh.parquet> <output.h5> [frame|last] [source_lat source_lon]")

mesh_path  = ARGS[1]
h5_path    = ARGS[2]
frame_arg  = length(ARGS) >= 3 ? ARGS[3] : "last"
have_source = length(ARGS) >= 5
source_lat_arg = have_source ? parse(Float64, ARGS[4]) : nothing
source_lon_arg = have_source ? parse(Float64, ARGS[5]) : nothing

mesh = load_mesh_geoparquet(mesh_path)
println("Loaded $(length(mesh.cells)) cells from $mesh_path")

# Build a FlowState with BOTH the legacy WLSQ infrastructure (built
# unconditionally by initialise_flow_model whenever gradient_correction
# is requested) and the diamond table (built because face_flux_method is
# :diamond here) in a single call — one state gives us everything needed
# to evaluate both constructions on the same field.
state = initialise_flow_model(mesh, StandardFlow();
                               gradient_correction       = true,
                               gradient_correction_alpha = 1.0,
                               face_flux_method           = :diamond)

n     = length(state.cell_ids)
edges = state.edges
ne    = edges.n_edges

# ── read WSE snapshot from the existing HDF5 output (no do-block: avoid
#    the top-level soft-scope pitfall this project's own audit scripts
#    have hit before) ──────────────────────────────────────────────────
fid = HDF5.h5open(h5_path, "r")
h5_ids = read(fid["mesh/cell_ids"])
if h5_ids != state.cell_ids
    error("mesh/cell_ids in $h5_path does not match the cell order of the " *
          "loaded mesh $mesh_path — point this script at the exact mesh " *
          "file that produced this .h5 output, or the per-cell arrays " *
          "below will be silently misaligned.")
end

frame_names = sort(collect(keys(fid["frames"])))
frame_name  = frame_arg == "last" ? frame_names[end] : lpad(frame_arg, 6, '0')
if !(frame_name in frame_names)
    error("frame '$frame_name' not found. Available: $(frame_names[1]) .. " *
          "$(frame_names[end]) ($(length(frame_names)) frames total).")
end

t     = read(fid["frames/$frame_name/t"])
depth = read(fid["frames/$frame_name/water_depth"])
close(fid)

println("Using frame $frame_name, t = $(round(t, digits=1))s, " *
        "wet cells (depth>1e-4) = $(count(>(1e-4), depth))/$n")

wse_all = state.elevation .+ depth

# ── legacy gradient pre-pass (exactly what step_standard! Phase A does
#    before its edge loop when gradient_correction is true) ────────────
_compute_wse_gradients!(state, wse_all)
alpha = state.gradient_correction_alpha

# ── source location for distance binning ────────────────────────────
if have_source
    src_lon0, src_lat0 = source_lon_arg, source_lat_arg
    println("Source (from CLI): lon=$src_lon0, lat=$src_lat0")
else
    src_idx = argmax(depth)
    src_lon0, src_lat0 = state.cell_lons[src_idx], state.cell_lats[src_idx]
    println("No source given — using max-depth cell as proxy: " *
            "lon=$(round(src_lon0, digits=5)), lat=$(round(src_lat0, digits=5))")
end

R = A5Grid._EARTH_R

legacy_dwse  = fill(NaN, ne)
diamond_dwse = fill(NaN, ne)
bearing_deg  = fill(NaN, ne)
dist_src_m   = fill(NaN, ne)
is_valid     = falses(ne)

for e in 1:ne
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]

    isfinite(wse_ci) && isfinite(wse_cj) || continue

    # ── legacy dWSE_n — must match step_standard! Phase A exactly ──
    gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
    gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
    vhat_dot = gx_f * edges.skew_x[e] + gy_f * edges.skew_y[e]
    legacy_dwse[e] = edges.cos_theta[e] * (wse_ci - wse_cj) - alpha * edges.L[e] * vhat_dot

    # ── diamond dWSE_n ──
    dn = diamond_dWSE_n(state.diamond_table, e, wse_ci, wse_cj, wse_all, edges.L[e])
    isnan(dn) && continue   # this edge fell back / had no valid diamond record

    diamond_dwse[e] = dn
    is_valid[e] = true

    # ── edge bearing (0=N, clockwise), i -> j ──
    lat_mid = 0.5 * (state.cell_lats[ci] + state.cell_lats[cj])
    dlat = state.cell_lats[cj] - state.cell_lats[ci]
    dlon = (state.cell_lons[cj] - state.cell_lons[ci]) * cosd(lat_mid)
    bearing_deg[e] = mod(atand(dlon, dlat), 360.0)

    # ── crude great-circle-ish distance from source to edge's cell i ──
    dlat_s = state.cell_lats[ci] - src_lat0
    dlon_s = (state.cell_lons[ci] - src_lon0) * cosd(state.cell_lats[ci])
    dist_src_m[e] = sqrt(dlat_s^2 + dlon_s^2) * deg2rad(1.0) * R
end

valid_idx = findall(is_valid)
n_valid   = length(valid_idx)
n_skipped = ne - n_valid
println("
$n_valid/$ne edges compared (diamond had a valid record); " *
        "$n_skipped skipped (invalid diamond record or non-finite WSE).")

diff = diamond_dwse[valid_idx] .- legacy_dwse[valid_idx]

println("
--- dWSE_n magnitude: diamond vs legacy (valid edges only) ---")
@printf("  legacy:  mean=%.6f  std=%.6f  max|.|=%.6f\n",
        mean(legacy_dwse[valid_idx]), std(legacy_dwse[valid_idx]),
        maximum(abs.(legacy_dwse[valid_idx])))
@printf("  diamond: mean=%.6f  std=%.6f  max|.|=%.6f\n",
        mean(diamond_dwse[valid_idx]), std(diamond_dwse[valid_idx]),
        maximum(abs.(diamond_dwse[valid_idx])))
@printf("  diff (diamond - legacy): mean=%.6f  std=%.6f  max|diff|=%.6f\n",
        mean(diff), std(diff), maximum(abs.(diff)))
println("  If diamond is systematically LESS damped than legacy near a " *
        "non-linear (front) region, expect std(diamond) > std(legacy) " *
        "and/or a heavy tail in |diff| concentrated at short source " *
        "distance — see the distance-binned table below.")

# ── bearing correlation of the divergence (mirrors Task 3's original
#    cos(2*bearing) correlation check on skew_x, applied here to the
#    diamond-vs-legacy DIFFERENCE rather than to skew_x directly) ──────
bearings_valid = bearing_deg[valid_idx]
cos2b = cosd.(2 .* bearings_valid)
corr_diff_bearing = cor(diff, cos2b)
println("
--- Bearing correlation of (diamond - legacy) ---")
@printf("  corr(diff, cos(2*bearing)) = %.4f\n", corr_diff_bearing)
println("  Task 3's original skew_x/heading correlation was ~0.30 and " *
        "resolution-flat (noise floor 0.005-0.02) — treat that as a rough " *
        "scale for 'a real signal' vs 'noise' here.")

# ── explicit north/south split, since the planar-slope domain's known
#    bias is specifically north/south, not a generic 2*bearing pattern ─
is_north = (bearings_valid .< 90.0) .| (bearings_valid .> 270.0)
is_south = .!is_north
println("
--- North/South split of (diamond - legacy) ---")
@printf("  mean diff | north-pointing edges (n=%d): %.6f\n",
        count(is_north), mean(diff[is_north]))
@printf("  mean diff | south-pointing edges (n=%d): %.6f\n",
        count(is_south), mean(diff[is_south]))
println("  A clean sign flip north vs south (mirroring Task 8's skew_x " *
        "finding: north mean -0.055, south mean +0.057) would directly " *
        "implicate the SAME underlying mesh chirality in both the legacy " *
        "V̂ term and the diamond construction, just expressed with less " *
        "damping in the diamond case.")

# ── distance-from-source binning: tests whether the divergence is a
#    general property of the construction, or concentrated specifically
#    near the front (where the WSE field is least linear, and where the
#    diamond's tighter stencil vs. legacy's 5-neighbour-averaged stencil
#    would be expected to diverge most) ─────────────────────────────────
dist_valid = dist_src_m[valid_idx]
n_bins = 8
lo_d, hi_d = extrema(dist_valid)
bin_edges = range(lo_d, hi_d, length = n_bins + 1)
println("
--- (diamond - legacy) by distance from source (proxy: max-depth cell " *
        "unless --source given) ---")
for b in 1:n_bins
    lo, hi = bin_edges[b], bin_edges[b+1]
    in_bin = (dist_valid .>= lo) .& (dist_valid .< (b == n_bins ? hi + 1.0 : hi))
    nb = count(in_bin)
    nb == 0 && continue
    @printf("  [%7.0f - %7.0f m]  n=%5d  mean=%+.6f  std=%.6f\n",
            lo, hi, nb, mean(diff[in_bin]), std(diff[in_bin]))
end

println("
Done.")
