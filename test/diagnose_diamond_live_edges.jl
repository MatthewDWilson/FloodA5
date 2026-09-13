#!/usr/bin/env julia
# test/diagnose_diamond_live_edges.jl
#
# Corrected version of diagnose_diamond_vs_legacy_gradient.jl /
# diagnose_diamond_front_localisation.jl. Those scripts compared raw
# dWSE_n magnitudes across ALL edges, including dry-dry pairs where, for
# standard flow, h_flow = max(WSE_i,WSE_j) - z_sill = 0 EXACTLY (since
# z_sill = max(elev_i,elev_j) and WSE=elevation when there's no water) --
# meaning the flux kernel returns Q=0 on those edges regardless of
# dWSE_n's value. The earlier top-20-worst-edges list was dominated by
# exactly this class (18/20 'both_dry'), which is numerically inert noise
# on flat/gently-sloped dry terrain, not a real effect.
#
# This script restricts EVERY comparison (magnitude, bearing correlation,
# north/south split, distance-from-source binning, wet/dry-class
# breakdown) to LIVE edges only: h_flow > HFLOW_THRESHOLD, computed
# identically to the production kernel (same edges.sill[e], same
# threshold constant) -- i.e. edges that can actually produce non-zero
# flux and therefore actually matter to the simulation's output.
#
# It reports both the gated and ungated picture side by side, so the
# scale of the "irrelevant noise" problem in the earlier scripts is
# visible directly, not just asserted.
#
# Usage:
#   julia --project=. test/diagnose_diamond_live_edges.jl \
#       <mesh.parquet> <output.h5> [frame|last] [wet_threshold_m] [source_lat source_lon]

using Statistics, HDF5, Printf

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

length(ARGS) >= 2 || error("Usage: julia diagnose_diamond_live_edges.jl " *
                            "<mesh.parquet> <output.h5> [frame|last] [wet_threshold_m] [source_lat source_lon]")

mesh_path     = ARGS[1]
h5_path       = ARGS[2]
frame_arg     = length(ARGS) >= 3 ? ARGS[3] : "last"
wet_threshold = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1e-4
have_source   = length(ARGS) >= 6
source_lat_arg = have_source ? parse(Float64, ARGS[5]) : nothing
source_lon_arg = have_source ? parse(Float64, ARGS[6]) : nothing

mesh = load_mesh_geoparquet(mesh_path)
println("Loaded $(length(mesh.cells)) cells from $mesh_path")

state = initialise_flow_model(mesh, StandardFlow();
                               gradient_correction       = true,
                               gradient_correction_alpha = 1.0,
                               face_flux_method           = :diamond)

n     = length(state.cell_ids)
edges = state.edges
ne    = edges.n_edges

fid = HDF5.h5open(h5_path, "r")
h5_ids = read(fid["mesh/cell_ids"])
if h5_ids != state.cell_ids
    error("mesh/cell_ids in $h5_path does not match the loaded mesh's cell " *
          "order -- point this script at the exact mesh file that produced " *
          "this .h5 output.")
end

frame_names = sort(collect(keys(fid["frames"])))
frame_name  = frame_arg == "last" ? frame_names[end] : lpad(frame_arg, 6, '0')
if !(frame_name in frame_names)
    error("frame '$frame_name' not found. Available: $(frame_names[1]) .. " *
          "$(frame_names[end]).")
end

t     = read(fid["frames/$frame_name/t"])
depth = read(fid["frames/$frame_name/water_depth"])
close(fid)

println("Using frame $frame_name, t = $(round(t, digits=1))s, " *
        "wet threshold = $(wet_threshold) m, " *
        "wet cells = $(count(>(wet_threshold), depth))/$n")

wse_all = state.elevation .+ depth
is_wet  = depth .> wet_threshold

_compute_wse_gradients!(state, wse_all)
alpha = state.gradient_correction_alpha

# ── source location for distance binning ─────────────────────────────
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
h_flow_arr   = fill(NaN, ne)
edge_class   = fill(:invalid, ne)   # :both_wet, :both_dry, :mixed
is_diamond_valid = falses(ne)

for e in 1:ne
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]
    isfinite(wse_ci) && isfinite(wse_cj) || continue

    # ── h_flow, computed EXACTLY as the production kernel does ──
    h_flow = max(wse_ci, wse_cj) - edges.sill[e]
    h_flow_arr[e] = h_flow

    gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
    gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
    vhat_dot = gx_f * edges.skew_x[e] + gy_f * edges.skew_y[e]
    legacy_dwse[e] = edges.cos_theta[e] * (wse_ci - wse_cj) - alpha * edges.L[e] * vhat_dot

    dn = diamond_dWSE_n(state.diamond_table, e, wse_ci, wse_cj, wse_all, edges.L[e])
    isnan(dn) && continue
    diamond_dwse[e] = dn
    is_diamond_valid[e] = true

    lat_mid = 0.5 * (state.cell_lats[ci] + state.cell_lats[cj])
    dlat = state.cell_lats[cj] - state.cell_lats[ci]
    dlon = (state.cell_lons[cj] - state.cell_lons[ci]) * cosd(lat_mid)
    bearing_deg[e] = mod(atand(dlon, dlat), 360.0)

    dlat_s = state.cell_lats[ci] - src_lat0
    dlon_s = (state.cell_lons[ci] - src_lon0) * cosd(state.cell_lats[ci])
    dist_src_m[e] = sqrt(dlat_s^2 + dlon_s^2) * deg2rad(1.0) * R

    if is_wet[ci] && is_wet[cj]
        edge_class[e] = :both_wet
    elseif !is_wet[ci] && !is_wet[cj]
        edge_class[e] = :both_dry
    else
        edge_class[e] = :mixed
    end
end

all_valid = findall(is_diamond_valid)
live_mask = h_flow_arr[all_valid] .> HFLOW_THRESHOLD
live_idx  = all_valid[live_mask]

println("
$(length(all_valid))/$ne edges have a valid diamond record.")
println("Of those, $(length(live_idx)) ($(round(100*length(live_idx)/length(all_valid),digits=1))%) " *
         "are LIVE (h_flow > $HFLOW_THRESHOLD m) -- i.e. can actually produce non-zero flux.")

function report_block(idx::Vector{Int}, label::String)
    isempty(idx) && (println("  [$label] no edges in this set"); return)
    diff = diamond_dwse[idx] .- legacy_dwse[idx]
    println("
=== $label (n=$(length(idx))) ===")
    @printf("  legacy:  mean=%.6f  std=%.6f  max|.|=%.6f\n",
            mean(legacy_dwse[idx]), std(legacy_dwse[idx]), maximum(abs.(legacy_dwse[idx])))
    @printf("  diamond: mean=%.6f  std=%.6f  max|.|=%.6f\n",
            mean(diamond_dwse[idx]), std(diamond_dwse[idx]), maximum(abs.(diamond_dwse[idx])))
    @printf("  diff (diamond-legacy): mean=%.6f  std=%.6f  max|diff|=%.6f\n",
            mean(diff), std(diff), maximum(abs.(diff)))

    bearings = bearing_deg[idx]
    cos2b = cosd.(2 .* bearings)
    @printf("  corr(diff, cos(2*bearing)) = %.4f\n", cor(diff, cos2b))

    is_n = (bearings .< 90.0) .| (bearings .> 270.0)
    is_s = .!is_n
    if count(is_n) > 0 && count(is_s) > 0
        @printf("  mean diff | north (n=%d): %+.6f   south (n=%d): %+.6f\n",
                count(is_n), mean(diff[is_n]), count(is_s), mean(diff[is_s]))
    end

    cls = edge_class[idx]
    for c in (:both_wet, :both_dry, :mixed)
        in_c = cls .== c
        nb = count(in_c)
        nb == 0 && continue
        @printf("  class %-10s n=%6d  mean|diff|=%.6f  (%.1f%% of set, %.1f%% of total|diff|)\n",
                String(c), nb, mean(abs.(diff[in_c])), 100*nb/length(idx),
                100*sum(abs.(diff[in_c]))/sum(abs.(diff)))
    end
end

report_block(all_valid, "UNGATED (all valid edges, including inert dry-dry pairs)")
report_block(live_idx,  "LIVE ONLY (h_flow > HFLOW_THRESHOLD -- what actually matters)")

# ── distance-from-source binning, live edges only ──
if !isempty(live_idx)
    diff_live = diamond_dwse[live_idx] .- legacy_dwse[live_idx]
    dist_live = dist_src_m[live_idx]
    n_bins = 8
    lo_d, hi_d = extrema(dist_live)
    if hi_d > lo_d
        bin_edges = range(lo_d, hi_d, length = n_bins + 1)
        println("
--- LIVE edges: (diamond - legacy) by distance from source ---")
        for b in 1:n_bins
            lo, hi = bin_edges[b], bin_edges[b+1]
            in_bin = (dist_live .>= lo) .& (dist_live .< (b == n_bins ? hi + 1.0 : hi))
            nb = count(in_bin)
            nb == 0 && continue
            @printf("  [%7.0f - %7.0f m]  n=%5d  mean=%+.6f  std=%.6f\n",
                    lo, hi, nb, mean(diff_live[in_bin]), std(diff_live[in_bin]))
        end
    end

    # top-20 worst LIVE edges
    order = sortperm(abs.(diff_live), rev=true)
    top_n = min(20, length(order))
    println("
--- Top $top_n LIVE edges by |diamond - legacy|, with class and h_flow ---")
    @printf("  %6s  %10s  %10s  %s\n", "edge#", "|diff|", "h_flow", "class")
    for k in 1:top_n
        idx = live_idx[order[k]]
        @printf("  %6d  %10.6f  %10.6f  %s\n",
                idx, abs(diff_live[order[k]]), h_flow_arr[idx], String(edge_class[idx]))
    end
end

println("
Done.")
