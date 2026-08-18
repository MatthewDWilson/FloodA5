#!/usr/bin/env julia
# ============================================================================
# diagnose_skew_bias.jl
#
# Loads a real mesh and builds the real EdgeList (via _build_edge_list ->
# _edge_geometry, exactly as initialise_flow_model does) -- NO simulation is
# run. Checks whether V̂ = (skew_x, skew_y), the WLSQ non-orthogonal
# correction vector, has a systematic north/south bias on the real,
# irregular A5 mesh.
#
# Motivation (2026-08-18 session): test_mirror_symmetry.jl hand-built a
# synthetic EdgeList with V̂ set directly to be an exact mirror pair -- it
# proved the Phase A FORMULA is symmetric given a symmetric V̂ input, but it
# never exercised _edge_geometry's actual per-edge orientation-flip logic
# (`if c_raw < 0.0; flip n̂; end`) on real, non-mirror-symmetric pentagon
# geometry. Real-mesh runs with --gradient-correction-alpha 1.0 show a
# CONSISTENT, non-saturating southward drift (not a bounded oscillation) --
# a pattern much more consistent with a small per-edge systematic bias
# compounding over thousands of edges/steps than with a resonant instability
# that should eventually settle. This script checks for that bias directly,
# at mesh-build time, with no simulation confound.
#
# What "unbiased" looks like on THIS mesh: the planar-embankment domain is
# translation-invariant in latitude (uniform east-west slope, ignoring the
# thin embankment strip) -- so if the A5 tiling and _edge_geometry are both
# free of systematic y-direction bias, skew_y should have mean ~0 across all
# edges, and critically, mean(skew_y | dy_m > 0) should be the negative of
# mean(skew_y | dy_m < 0) (i.e. "north-pointing" and "south-pointing" edges
# should be statistical mirror images of each other in aggregate, even
# though individual pentagon cells are irregular).
#
# Usage: julia --project=. test\diagnose_skew_bias.jl <mesh.parquet>
#   e.g. julia --project=. test\diagnose_skew_bias.jl test/planar_embankment/planar_mesh18_std.parquet
# ============================================================================

using Statistics

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_path = length(ARGS) >= 1 ? ARGS[1] :
            joinpath("test", "planar_embankment", "planar_mesh18_std.parquet")

println("=" ^ 76)
println("diagnose_skew_bias.jl -- real-mesh V̂ north/south bias check")
println("Mesh: $mesh_path")
println("=" ^ 76)

mesh = load_mesh_geoparquet(mesh_path)
# Build the flow model purely to get the real EdgeList out of
# _build_edge_list/_edge_geometry -- gradient_correction doesn't affect edge
# geometry construction (that happens regardless of the flag), so this is
# safe and touches no simulation state.
state = initialise_flow_model(mesh, StandardFlow(); gradient_correction=true)

edges = state.edges
n_edges = edges.n_edges
skew_x, skew_y = edges.skew_x, edges.skew_y
dy_m, dx_m     = edges.dy_m, edges.dx_m
width          = edges.width

println("\nTotal edges: $n_edges")

# ── 1. Global mean/weighted-mean of skew_y ──────────────────────────────────
mean_skew_y   = mean(skew_y)
wmean_skew_y  = sum(skew_y .* width) / sum(width)
println("\n--- Global skew_y statistics ---")
println("  mean(skew_y)                    = $mean_skew_y")
println("  width-weighted mean(skew_y)     = $wmean_skew_y")
println("  std(skew_y)                     = $(std(skew_y))")
println("  (an unbiased mesh/formula should give ~0 here, well within the")
println("   noise floor set by std/sqrt(n_edges) = $(std(skew_y)/sqrt(n_edges)))")

# ── 2. North-pointing vs south-pointing edges ───────────────────────────────
# Every edge appears once (cell_i < cell_j convention), so "north-pointing"
# here means cell_j is north of cell_i for that specific edge -- there is no
# double-counting to worry about, but also no guarantee of an exact 50/50
# split (mesh isn't perfectly regular), so we look at MEANS, not just counts.
north_mask = dy_m .> 0
south_mask = dy_m .< 0
n_north, n_south = count(north_mask), count(south_mask)

mean_sy_north = mean(skew_y[north_mask])
mean_sy_south = mean(skew_y[south_mask])

println("\n--- North-pointing vs south-pointing edges (by sign of dy_m) ---")
println("  n_north (dy_m>0) = $n_north   mean(skew_y) = $mean_sy_north")
println("  n_south (dy_m<0) = $n_south   mean(skew_y) = $mean_sy_south")
println("  sum                              = $(mean_sy_north + mean_sy_south)")
println("  (unbiased: north and south means should be near-exact negatives")
println("   of each other, i.e. the sum above should be ~0. A nonzero sum")
println("   means north-pointing and south-pointing edges are NOT behaving")
println("   as statistical mirror images -- direct evidence of a systematic")
println("   bias in _edge_geometry's real orientation-flip logic, not just")
println("   random irregular-mesh noise.)")

# ── 3. Correlation between sign(dy_m) and sign(skew_y) ──────────────────────
# If the orientation-flip logic is unbiased, skew_y's sign should be roughly
# independent of dy_m's sign in aggregate terms of *how often it agrees*,
# even though individual pentagons vary -- check the fraction where
# sign(skew_y) == sign(dy_m) among edges with |dy_m| large enough to have a
# clear north/south character (filters out near-east-west edges where sign
# is noise-dominated).
clear_ns = abs.(dy_m) .> 0.3 .* sqrt.(dx_m.^2 .+ dy_m.^2)   # |dy_m| > 30% of L
same_sign = count(sign.(skew_y[clear_ns]) .== sign.(dy_m[clear_ns]))
total_clear = count(clear_ns)
println("\n--- Sign agreement (edges with a clear N/S character) ---")
println("  edges with |dy_m| > 30% of L : $total_clear / $n_edges")
println("  sign(skew_y) == sign(dy_m)  : $same_sign / $total_clear " *
        "($(round(100*same_sign/max(total_clear,1), digits=1))%)")
println("  (unbiased: should be close to 50%. Far from 50% means the")
println("   orientation-flip has a systematic preference correlated with")
println("   compass direction, not just per-edge geometric noise.)")

# ── 4. Per-latitude-band breakdown (does bias vary with distance from the
#      injection source, or is it uniform across the whole domain?) ────────
lat0 = mean(state.cell_lats)
edge_lat = [0.5 * (state.cell_lats[edges.cell_i[e]] + state.cell_lats[edges.cell_j[e]])
            for e in 1:n_edges]
bands = range(minimum(edge_lat), maximum(edge_lat), length=6)
println("\n--- skew_y mean by latitude band (checks spatial uniformity) ---")
for b in 1:5
    m = (edge_lat .>= bands[b]) .& (edge_lat .< bands[b+1])
    count(m) == 0 && continue
    println("  lat [$(round(bands[b],digits=4)), $(round(bands[b+1],digits=4))): " *
            "n=$(count(m))  mean(skew_y)=$(round(mean(skew_y[m]), digits=6))")
end

println("\n" * "=" ^ 76)
println("If the global/north-south means are near-zero and sign agreement is")
println("near 50%: no systematic geometric bias found here -- the alpha=1")
println("drift likely comes from somewhere else (e.g. compounding through")
println("the SIMULATION's volume/WSE feedback rather than the static")
println("geometry itself), and the fix is more likely a per-edge adaptive")
println("or capped alpha than a geometry bug.")
println()
println("If north/south means clearly don't cancel, or sign agreement is far")
println("from 50%: this is a real, fixable geometry-construction bias in")
println("_edge_geometry's orientation-flip logic on irregular pentagons --")
println("worth a targeted look at that function specifically.")
println("=" ^ 76)
