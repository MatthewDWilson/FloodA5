#!/usr/bin/env julia
# ============================================================================
# audit_edge_correlation.jl
#
# Phase A, Task 3 of FloodA5_PhaseA_ImplementationScope.md.
#
# Extends test/diagnose_skew_bias.jl (which already computes the skew_x/
# skew_y north-south mirror-statistics check) rather than duplicating it —
# this script does NOT repeat that check. Instead it adds the three things
# Task 3 asks for that diagnose_skew_bias.jl does not currently produce:
#
#   1. The FULL per-edge non-orthogonality angle (δ_e) distribution — not
#      just mean/p95/max summary statistics (already logged at mesh-build
#      time and recorded in prior sessions' handover docs).
#   2. Correlation of δ_e, and of |V̂| = |(skew_x, skew_y)|, against the
#      COMPASS BEARING of each edge — is non-orthogonality itself
#      heading-dependent, independent of the skew_x sign-flip already found?
#   3. A resolution comparison — run against two (or more) real meshes at
#      different resolutions in one invocation, to see whether any
#      correlation strengthens, weakens, or stays flat with refinement.
#
# It also numerically confirms the documented |V̂| = sin(δ_e) relationship
# holds edge-by-edge (not just as a mean statistic), since that identity is
# asserted in several prior handover docs but has not been checked at the
# per-edge level before.
#
# Usage:
#   julia --project=. test/audit_edge_correlation.jl <mesh1.parquet> [mesh2.parquet ...]
#
#   e.g. julia --project=. test/audit_edge_correlation.jl \
#            test/square/square_mesh_res14.parquet \
#            test/planar_embankment/planar_mesh16_std.parquet \
#            test/planar_embankment/planar_mesh18_std.parquet
#
# Each mesh is reported independently, followed by a cross-mesh resolution
# summary table at the end. Per-edge data is written to a small CSV per
# mesh (<meshname>_edge_correlation.csv) alongside stdout output, so it can
# be plotted/histogrammed separately if useful — this script itself only
# prints numeric summaries, no plotting dependency required.
# ============================================================================

using Statistics

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_paths = length(ARGS) >= 1 ? ARGS :
             [joinpath("test", "planar_embankment", "planar_mesh18_std.parquet")]

"""
    _analyse_mesh(path) → NamedTuple

Loads one mesh, builds the real EdgeList via initialise_flow_model (exactly
as production does — gradient_correction doesn't affect edge geometry
construction, so this is safe and touches no simulation state), and returns
the per-edge arrays plus summary statistics needed for the cross-mesh table.
"""
function _analyse_mesh(path::String)
    println("\n" * "=" ^ 78)
    println("Mesh: $path")
    println("=" ^ 78)

    mesh  = load_mesh_geoparquet(path)
    state = initialise_flow_model(mesh, StandardFlow(); gradient_correction=true)
    edges = state.edges
    n_edges = edges.n_edges

    skew_x, skew_y = edges.skew_x, edges.skew_y
    dx_m, dy_m     = edges.dx_m, edges.dy_m
    cos_theta      = edges.cos_theta
    width          = edges.width

    valid = isfinite.(cos_theta) .& (cos_theta .>= -1.0) .& (cos_theta .<= 1.0)
    n_valid = count(valid)
    println("Total edges: $n_edges   valid geometry: $n_valid")

    # ── 1. Full δ_e (non-orthogonality angle) distribution ─────────────────
    delta_deg = acosd.(clamp.(cos_theta[valid], -1.0, 1.0))
    println("\n--- Non-orthogonality angle δ_e, FULL distribution (degrees) ---")
    println("  min    = $(round(minimum(delta_deg), digits=2))")
    println("  p5     = $(round(quantile(delta_deg, 0.05), digits=2))")
    println("  p25    = $(round(quantile(delta_deg, 0.25), digits=2))")
    println("  mean   = $(round(mean(delta_deg), digits=2))")
    println("  median = $(round(median(delta_deg), digits=2))")
    println("  p75    = $(round(quantile(delta_deg, 0.75), digits=2))")
    println("  p95    = $(round(quantile(delta_deg, 0.95), digits=2))")
    println("  max    = $(round(maximum(delta_deg), digits=2))")
    println("  std    = $(round(std(delta_deg), digits=2))")
    # Coarse histogram (10 bins, 0-90°) so the shape (not just quantiles) is visible.
    println("  Histogram (10° bins, 0-90°):")
    for lo in 0:10:80
        cnt = count(d -> lo <= d < lo+10, delta_deg)
        bar = "#" ^ min(60, round(Int, 60 * cnt / max(n_valid,1)))
        println("    [$lo,$(lo+10))° : $(lpad(cnt, 6))  $bar")
    end

    # ── 2a. δ_e vs bearing correlation ──────────────────────────────────────
    bearing_deg = mod.(atand.(dx_m[valid], dy_m[valid]), 360.0)   # 0°=N, clockwise
    println("\n--- δ_e (non-orthogonality) by bearing octant ---")
    octant_names = ["N","NE","E","SE","S","SW","W","NW"]
    octant_delta_means = Float64[]
    for k in 0:7
        lo_b, hi_b = k*45.0, (k+1)*45.0
        m = (bearing_deg .>= lo_b) .& (bearing_deg .< hi_b)
        cnt = count(m)
        md = cnt > 0 ? mean(delta_deg[m]) : NaN
        push!(octant_delta_means, md)
        cnt > 0 && println("  $(octant_names[k+1]) [$(lo_b)°,$(hi_b)°): n=$cnt   mean δ_e = $(round(md, digits=2))°")
    end
    delta_octant_range = maximum(filter(isfinite, octant_delta_means)) -
                          minimum(filter(isfinite, octant_delta_means))
    println("  Range across octants: $(round(delta_octant_range, digits=2))° " *
            "(large relative to overall std=$(round(std(delta_deg),digits=2))° " *
            "⇒ heading-dependent; small ⇒ heading-independent)")

    # Pearson correlation between bearing (as a 2D unit vector, since bearing
    # itself is circular/non-linear) and δ_e: use cos(2·bearing) and
    # sin(2·bearing) as linear proxies for "axis alignment" (period-180°,
    # appropriate since an edge and its reverse have the same δ_e).
    cos2b = cosd.(2 .* bearing_deg)
    sin2b = sind.(2 .* bearing_deg)
    corr_cos2b = cor(delta_deg, cos2b)
    corr_sin2b = cor(delta_deg, sin2b)
    println("  Correlation(δ_e, cos(2·bearing)) = $(round(corr_cos2b, digits=4))")
    println("  Correlation(δ_e, sin(2·bearing)) = $(round(corr_sin2b, digits=4))")
    println("  (both near 0 ⇒ no linear axis-alignment relationship; a value with")
    println("   |r| notably above the sampling noise floor ~1/sqrt(n) = " *
            "$(round(1/sqrt(n_valid), digits=4)) is worth a closer look)")

    # ── 2b. |V̂| vs δ_e — confirm |V̂| = sin(δ_e) edge-by-edge ────────────────
    vhat_mag = sqrt.(skew_x[valid].^2 .+ skew_y[valid].^2)
    sin_delta = sind.(delta_deg)
    resid = vhat_mag .- sin_delta
    println("\n--- |V̂| vs sin(δ_e): per-edge identity check ---")
    println("  mean(|V̂|)      = $(round(mean(vhat_mag), digits=4))")
    println("  mean(sin δ_e)  = $(round(mean(sin_delta), digits=4))")
    println("  mean(|V̂| - sin δ_e) = $(round(mean(resid), digits=6))")
    println("  max |residual|      = $(round(maximum(abs.(resid)), digits=6))")
    println("  (should be ~0 at every edge, not just on average, if the")
    println("   documented |V̂|=sinθ identity holds exactly — a nonzero MAX")
    println("   residual with a near-zero MEAN would indicate the identity holds")
    println("   only in aggregate, hiding per-edge cancellation.)")

    # ── 2c. |V̂| magnitude vs bearing (does the MAGNITUDE of the correction,")
    #      not just skew_x's sign, vary by heading?) ────────────────────────
    println("\n--- |V̂| by bearing octant ---")
    for k in 0:7
        lo_b, hi_b = k*45.0, (k+1)*45.0
        m = (bearing_deg .>= lo_b) .& (bearing_deg .< hi_b)
        cnt = count(m)
        cnt == 0 && continue
        println("  $(octant_names[k+1]) [$(lo_b)°,$(hi_b)°): n=$cnt   mean |V̂| = " *
                "$(round(mean(vhat_mag[m]), digits=4))")
    end

    # ── 3. Save per-edge CSV for external plotting if desired ──────────────
    csv_path = replace(basename(path), ".parquet" => "") * "_edge_correlation.csv"
    open(csv_path, "w") do io
        println(io, "delta_deg,bearing_deg,skew_x,skew_y,vhat_mag,sin_delta,width")
        for (d, b, sx, sy, vm, sd, w) in zip(delta_deg, bearing_deg,
                                               skew_x[valid], skew_y[valid],
                                               vhat_mag, sin_delta, width[valid])
            println(io, "$d,$b,$sx,$sy,$vm,$sd,$w")
        end
    end
    println("\nPer-edge data written to: $csv_path ($n_valid rows)")

    return (path=path, n_cells=length(mesh.cells), n_edges=n_edges,
            mean_delta=mean(delta_deg), p95_delta=quantile(delta_deg, 0.95),
            max_delta=maximum(delta_deg),
            mean_vhat=mean(vhat_mag), max_vhat=maximum(vhat_mag),
            corr_cos2b=corr_cos2b, corr_sin2b=corr_sin2b,
            delta_octant_range=delta_octant_range)
end

results = [_analyse_mesh(p) for p in mesh_paths]

if length(results) > 1
    println("\n\n" * "=" ^ 78)
    println("CROSS-MESH RESOLUTION COMPARISON")
    println("=" ^ 78)
    println(rpad("mesh", 40), rpad("n_cells", 10), rpad("n_edges", 10),
            rpad("mean δ°", 9), rpad("p95 δ°", 8), rpad("max δ°", 8),
            rpad("mean|V̂|", 9), rpad("corr(cos2b)", 12), "octant_range°")
    for r in results
        println(rpad(basename(r.path), 40), rpad(r.n_cells, 10), rpad(r.n_edges, 10),
                rpad(round(r.mean_delta, digits=1), 9), rpad(round(r.p95_delta, digits=1), 8),
                rpad(round(r.max_delta, digits=1), 8), rpad(round(r.mean_vhat, digits=3), 9),
                rpad(round(r.corr_cos2b, digits=4), 12), round(r.delta_octant_range, digits=2))
    end
    println()
    println("Interpretation guide (per FloodA5_PhaseA_ImplementationScope.md Task 3):")
    println("  - If mean δ°, p95 δ°, mean|V̂| stay roughly CONSTANT across resolutions:")
    println("    non-orthogonality is a fixed property of the A5 pentagon tiling itself,")
    println("    not a refinement artefact — consistent with genuine topological")
    println("    chirality rather than ordinary truncation error.")
    println("  - If octant_range or |corr(cos2b)| GROWS or SHRINKS systematically with")
    println("    resolution: combine with Task 4's convergence-sweep result before")
    println("    drawing a conclusion — a shrinking bias with refinement points to")
    println("    truncation error; a flat or growing one strengthens the chirality")
    println("    hypothesis (see the Phase A → Phase B decision matrix).")
end

println("\nDone.")
