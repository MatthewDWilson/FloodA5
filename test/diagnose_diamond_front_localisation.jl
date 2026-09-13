#!/usr/bin/env julia
# test/diagnose_diamond_front_localisation.jl
#
# Follow-up to diagnose_diamond_vs_legacy_gradient.jl. That script found
# the diamond-vs-legacy dWSE_n divergence is NOT bearing-correlated and
# NOT north/south split (ruling out the Task 3/8 mesh-chirality signature
# as the driver of THIS effect), but IS sharply localized to a narrow
# distance-from-source band -- consistent with a front-proximity effect
# rather than a domain-wide geometric one.
#
# This script tests that directly: does |diamond - legacy| correlate with
# an edge connecting a WET cell to a DRY (or near-dry) cell, rather than
# just happening to sit in a particular distance range? This is the
# precise signature predicted by the "unlimited linear reconstruction
# behaves like an unlimited higher-order FVM scheme near a discontinuity"
# hypothesis (Gibbs-type overshoot at a sharp wet/dry front) -- distance
# from source was only ever a coarse proxy for "near the current front."
#
# Usage:
#   julia --project=. test/diagnose_diamond_front_localisation.jl \
#       <mesh.parquet> <output.h5> [frame|last] [wet_threshold_m]
#
# wet_threshold_m: depth threshold (m) below which a cell counts as "dry"
# for the purposes of edge classification. Default 1e-4 m, matching the
# convention already used elsewhere in this codebase (e.g. water_depth
# wet/dry checks in test_planar_symmetry.jl).

using Statistics, HDF5, Printf

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

length(ARGS) >= 2 || error("Usage: julia diagnose_diamond_front_localisation.jl " *
                            "<mesh.parquet> <output.h5> [frame|last] [wet_threshold_m]")

mesh_path      = ARGS[1]
h5_path        = ARGS[2]
frame_arg      = length(ARGS) >= 3 ? ARGS[3] : "last"
wet_threshold  = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1e-4

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

wse_all  = state.elevation .+ depth
is_wet   = depth .> wet_threshold

_compute_wse_gradients!(state, wse_all)
alpha = state.gradient_correction_alpha

legacy_dwse  = fill(NaN, ne)
diamond_dwse = fill(NaN, ne)
edge_class   = fill(:invalid, ne)   # :both_wet, :both_dry, :mixed, :invalid

for e in 1:ne
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]
    isfinite(wse_ci) && isfinite(wse_cj) || continue

    gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
    gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
    vhat_dot = gx_f * edges.skew_x[e] + gy_f * edges.skew_y[e]
    legacy_dwse[e] = edges.cos_theta[e] * (wse_ci - wse_cj) - alpha * edges.L[e] * vhat_dot

    dn = diamond_dWSE_n(state.diamond_table, e, wse_ci, wse_cj, wse_all, edges.L[e])
    isnan(dn) && continue
    diamond_dwse[e] = dn

    if is_wet[ci] && is_wet[cj]
        edge_class[e] = :both_wet
    elseif !is_wet[ci] && !is_wet[cj]
        edge_class[e] = :both_dry
    else
        edge_class[e] = :mixed
    end
end

valid_idx = findall(edge_class .!= :invalid)
diff_abs  = abs.(diamond_dwse[valid_idx] .- legacy_dwse[valid_idx])
classes   = edge_class[valid_idx]

println("
$(length(valid_idx))/$ne edges classified.")

println("
--- |diamond - legacy| by wet/dry edge class ---")
for cls in (:both_wet, :both_dry, :mixed)
    in_cls = classes .== cls
    nb = count(in_cls)
    nb == 0 && continue
    d = diff_abs[in_cls]
    @printf("  %-10s  n=%6d  mean|diff|=%.6f  median|diff|=%.6f  max|diff|=%.6f  (%.1f%% of edges, %.1f%% of total |diff|)\n",
            String(cls), nb, mean(d), median(d), maximum(d),
            100*nb/length(valid_idx), 100*sum(d)/sum(diff_abs))
end

println("
If 'mixed' edges show a mean|diff| an order of magnitude (or more) above " *
        "'both_wet'/'both_dry', and account for a share of total |diff| far " *
        "exceeding their share of edge count, that directly confirms the " *
        "front-discontinuity hypothesis: the diamond-legacy divergence is " *
        "driven by cells straddling the wet/dry boundary, not by a domain-" *
        "wide geometric/chirality effect.")

# ── zoom in: rank the top-20 worst edges by |diff| and report their class ──
order = sortperm(diff_abs, rev=true)
top_n = min(20, length(order))
println("
--- Top $top_n edges by |diamond - legacy|, with class ---")
@printf("  %6s  %10s  %s\n", "edge#", "|diff|", "class")
for k in 1:top_n
    idx = valid_idx[order[k]]
    @printf("  %6d  %10.6f  %s\n", idx, diff_abs[order[k]], String(edge_class[idx]))
end
n_mixed_in_top = count(==(:mixed), edge_class[valid_idx[order[1:top_n]]])
println("
$n_mixed_in_top / $top_n of the worst edges are 'mixed' (wet/dry) edges.")

println("
Done.")
