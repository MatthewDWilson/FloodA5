#!/usr/bin/env julia
# test/diagnose_diamond_boundary_proximity.jl
#
# Follow-up to diagnose_diamond_live_edges.jl. That script found live-edge
# (h_flow > HFLOW_THRESHOLD) divergence between diamond and legacy dWSE_n
# is near-zero at t=35702s (well before the front reaches the closed
# boundary) but ~140x larger in std by t=72000s (at/after wall contact).
#
# This tests the specific hypothesis that follows from that: boundary
# outflow (step_standard! Phase D, _bates_ghost_flux) is NOT wired to
# face_flux_method -- it always uses the plain, uncorrected legacy Bates
# kernel regardless of --gradient-correction / --face-flux-method. If
# that inconsistency (direction-aware diamond flux on interior edges vs.
# isotropic uncorrected flux at the boundary) is what's driving the
# late-time divergence, live-edge |diamond - legacy| should concentrate
# specifically near boundary cells (state.boundary_mask), not be uniform
# across the domain.
#
# Method: BFS hop-distance from every boundary cell (state.boundary_mask
# == true) across the interior adjacency graph, then bin LIVE edges'
# |diamond - legacy| by min(hop_dist[ci], hop_dist[cj]).
#
# Usage:
#   julia --project=. test/diagnose_diamond_boundary_proximity.jl \
#       <mesh.parquet> <output.h5> [frame|last] [max_hops]

using Statistics, HDF5, Printf

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

length(ARGS) >= 2 || error("Usage: julia diagnose_diamond_boundary_proximity.jl " *
                            "<mesh.parquet> <output.h5> [frame|last] [max_hops]")

mesh_path = ARGS[1]
h5_path   = ARGS[2]
frame_arg = length(ARGS) >= 3 ? ARGS[3] : "last"
max_hops  = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 8

mesh = load_mesh_geoparquet(mesh_path)
println("Loaded $(length(mesh.cells)) cells from $mesh_path")

state = initialise_flow_model(mesh, StandardFlow();
                               gradient_correction       = true,
                               gradient_correction_alpha = 1.0,
                               face_flux_method           = :diamond)

n     = length(state.cell_ids)
edges = state.edges
ne    = edges.n_edges

n_boundary = count(state.boundary_mask)
println("Boundary cells: $n_boundary / $n ($(round(100*n_boundary/n, digits=1))%)")

# ── BFS hop-distance from any boundary cell, over the interior adjacency
#    graph (adj_matrix -- same graph the interior flux loop uses; ghost
#    edges themselves are a separate structure and not traversed here,
#    we only need "how many interior hops to the nearest boundary cell") ──
hop_dist = fill(typemax(Int), n)
queue = Int[]
for i in 1:n
    if state.boundary_mask[i]
        hop_dist[i] = 0
        push!(queue, i)
    end
end
qhead = 1
while qhead <= length(queue)
    global qhead
    ci = queue[qhead]; qhead += 1
    d = hop_dist[ci]
    d >= max_hops && continue
    for s in 1:size(state.adj_matrix, 1)
        cj = state.adj_matrix[s, ci]
        cj == 0 && continue
        if hop_dist[cj] > d + 1
            hop_dist[cj] = d + 1
            push!(queue, cj)
        end
    end
end
n_unreached = count(==(typemax(Int)), hop_dist)
println("Cells with hop_dist > $max_hops (or unreached): $n_unreached / $n")

# ── read WSE snapshot ──────────────────────────────────────────────────
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
    error("frame '$frame_name' not found. Available: $(frame_names[1]) .. $(frame_names[end]).")
end
t     = read(fid["frames/$frame_name/t"])
depth = read(fid["frames/$frame_name/water_depth"])
close(fid)

println("Using frame $frame_name, t = $(round(t, digits=1))s, " *
        "wet cells = $(count(>(1e-4), depth))/$n")

wse_all = state.elevation .+ depth
_compute_wse_gradients!(state, wse_all)
alpha = state.gradient_correction_alpha

diff_live      = Float64[]
hopdist_live   = Int[]
n_live = 0

for e in 1:ne
    global n_live
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]
    isfinite(wse_ci) && isfinite(wse_cj) || continue

    h_flow = max(wse_ci, wse_cj) - edges.sill[e]
    h_flow > HFLOW_THRESHOLD || continue

    gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
    gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
    vhat_dot = gx_f * edges.skew_x[e] + gy_f * edges.skew_y[e]
    legacy = edges.cos_theta[e] * (wse_ci - wse_cj) - alpha * edges.L[e] * vhat_dot

    dn = diamond_dWSE_n(state.diamond_table, e, wse_ci, wse_cj, wse_all, edges.L[e])
    isnan(dn) && continue

    n_live += 1
    push!(diff_live, abs(dn - legacy))
    push!(hopdist_live, min(hop_dist[ci], hop_dist[cj]))
end

println("
$n_live live edges compared.")

println("
--- LIVE edges: |diamond - legacy| by hop-distance from nearest boundary cell ---")
@printf("  %8s  %6s  %10s  %10s  %10s\n", "hop_dist", "n", "mean|diff|", "median", "max|diff|")
for h in 0:max_hops
    in_h = hopdist_live .== h
    nb = count(in_h)
    nb == 0 && continue
    d = diff_live[in_h]
    @printf("  %8d  %6d  %10.6f  %10.6f  %10.6f\n", h, nb, mean(d), median(d), maximum(d))
end
in_far = hopdist_live .> max_hops
nb_far = count(in_far)
if nb_far > 0
    d = diff_live[in_far]
    @printf("  %8s  %6d  %10.6f  %10.6f  %10.6f\n", ">$max_hops", nb_far, mean(d), median(d), maximum(d))
end

println("
If mean|diff| is clearly elevated at hop_dist=0-2 and decays toward the " *
        "interior, that directly confirms the boundary/ghost-edge method-" *
        "mismatch hypothesis: interior edges use direction-aware diamond " *
        "flux, ghost edges use isotropic uncorrected flux, and that seam " *
        "is producing the divergence -- not a domain-wide geometric effect.")

# ── overall correlation, as a single summary number ──
if length(diff_live) > 2
    c = cor(Float64.(hopdist_live), diff_live)
    println("
Overall corr(hop_dist, |diff|) = $(round(c, digits=4)) " *
            "(negative = divergence decays away from the boundary, " *
            "supporting the hypothesis)")
end

println("
Done.")
