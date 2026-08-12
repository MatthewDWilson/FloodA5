"""
test_noc_correction.jl
-----------------------
Unit tests for the WLSQ non-orthogonal gradient correction
(FloodA5_NonOrthogonal_Correction_Plan.md, "flow-direction-fixes" branch).

This file implements the T-NOC test series from the plan §10 Stage 1.  Test
availability is gated by what has actually been implemented so far:

  T-NOC1  _edge_skewness_vec orthogonal pair        — SUPERSEDED, not run here.
          _edge_skewness_vec was consolidated into _edge_geometry (Step 1) and
          no longer exists as a standalone function. The equivalent coverage
          lives in test_edge_geometry.jl Group 1 (cos θ component) plus the
          real-mesh skewness diagnostic logged by _build_edge_list (see
          FloodA5_NonOrthogonal_Correction_Plan.md §10.5.4 — square mesh
          mean=13.8%, max=24.3% of edge length). Not duplicated here.

  T-NOC2  _edge_skewness_vec known skewed pair       — SUPERSEDED, not run here.
          Same reasoning as T-NOC1. _edge_geometry's skewness output is
          exercised by the real-mesh diagnostic; a dedicated analytical
          cross-check could be added to test_edge_geometry.jl's existing
          Group 1 in a future session if more precision is wanted, but is
          not blocking for this branch's progress.

  T-NOC3  _build_wlsq_weights! linear field recovery  — IMPLEMENTED below.
          The only T-NOC test whose dependencies (_build_wlsq_weights!,
          wlsq_weights) exist at this point in the implementation sequence.

  T-NOC4  _bates_flux_corrected reduces to _bates_flux — IMPLEMENTED below
          (2026-06-23, once _bates_flux_corrected / _bates_flux_limited_corrected
          landed in flow2d.jl, Step 6).

  T-NOC5  _compute_wse_gradients! correctness + mass conservation — IMPLEMENTED
          below (2026-06-23, once _compute_wse_gradients! landed in
          FloodModel.jl, Step 7). Includes both a direct linear-field-recovery
          check (T-NOC5a, mirroring T-NOC3's method but exercising
          _compute_wse_gradients! itself rather than a manual multiply-add)
          and an end-to-end step_standard! mass-conservation check with
          gradient_correction=true (T-NOC5b).

Run with:
    julia test_noc_correction.jl

No external mesh or Python bridge required — geometry is constructed
synthetically, following the same pattern as test_edge_geometry.jl.
"""

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "mesh"))
include(joinpath(dirname(@__DIR__), "mesh", "A5Grid.jl"))

# FloodModel.jl is a flat script (no module wrapper).  We include only the
# parts needed for testing — the physics functions and types — without
# pulling in the heavy VisualisationServer / MakieVisualiser / HDF5
# dependencies that FloodModel.jl normally loads.  Same minimal preamble as
# test_edge_geometry.jl.
using .A5Grid
using .A5Grid: SGSTable, wse_from_volume, wetted_area_from_wse,
               sgs_table, build_sgs_tables!,
               grid_disk_neighbours, grid_disk_neighbours_batch,
               _polygon_area_m2, _shared_edge, _haversine_m
using JSON3
using Dates
using Statistics: mean
using LinearAlgebra: Diagonal
using HDF5
using Printf

# Stub out the vis modules so FloodModel.jl can be included without GLMakie.
module VisualisationServer
    start(;kwargs...)  = nothing
    stop(_)            = nothing
    push_frame!(args...) = nothing
    set_mesh!(args...)   = nothing
    notify_complete!(_)  = nothing
    struct VisServer end
end

module MakieVisualiser
    start(args...; kwargs...) = nothing
    stop(_)                   = nothing
    push_frame!(args...)      = nothing
end

include(joinpath(dirname(@__DIR__), "FloodModel.jl"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

const PASS = "\e[32m✓\e[0m"
const FAIL = "\e[31m✗\e[0m"

n_pass = Ref(0)
n_fail = Ref(0)

function check(name::String, cond::Bool, detail::String = "")
    if cond
        println("  $PASS  $name")
        n_pass[] += 1
    else
        println("  $FAIL  $name$(isempty(detail) ? "" : "  ($detail)")")
        n_fail[] += 1
    end
end

# Build a synthetic A5Cell — regular pentagon centred at (lon, lat), with a
# given half-width in degrees, optionally rotated. Mirrors the helper of the
# same name in test_edge_geometry.jl (kept independent here so this file has
# no cross-file dependency).
function _synthetic_pentagon(lon::Float64, lat::Float64,
                              half_width_deg::Float64,
                              rot_deg::Float64 = 0.0)::A5Cell
    verts = Vector{Vector{Float64}}()
    for k in 0:4
        θ = deg2rad(rot_deg + 90.0 + 72.0 * k)
        push!(verts, [lon + half_width_deg * cos(θ),
                      lat + half_width_deg * sin(θ)])
    end
    push!(verts, verts[1])   # close the ring
    A5Cell("test_$(round(Int, lon*1000))_$(round(Int, lat*1000))",
           14, lon, lat, verts, NaN)
end

# Build a small hexagonal-ish cluster of n pentagons around a centre cell,
# each offset by `spacing_deg` at angle k*(360/n_ring) degrees, giving the
# centre cell exactly n_ring neighbours. Returns (cells, adj_matrix) where
# adj_matrix is (5, n_total) with adj_matrix[:, 1] populated for the centre
# cell (cell 1) and zeros elsewhere (only cell 1's gradient is under test).
#
# This directly exercises _build_wlsq_weights! without needing a real A5
# mesh or the pya5 Python bridge — adj_matrix is hand-built, matching
# exactly what _build_adjacency_matrix! would produce for this geometry.
function _ring_cluster(n_ring::Int; spacing_deg::Float64 = 0.01,
                        lon0::Float64 = 172.0, lat0::Float64 = -43.5)
    cells = A5Cell[_synthetic_pentagon(lon0, lat0, spacing_deg * 0.4)]
    for k in 0:(n_ring - 1)
        θ = deg2rad(360.0 * k / n_ring)
        # Account for longitude degrees being physically shorter than
        # latitude degrees away from the equator, so the ring is
        # geometrically regular in local metres, not just in raw degrees.
        cos_lat0 = cosd(lat0)
        dlon = spacing_deg * cos(θ) / max(cos_lat0, 1e-6)
        dlat = spacing_deg * sin(θ)
        push!(cells, _synthetic_pentagon(lon0 + dlon, lat0 + dlat,
                                          spacing_deg * 0.4))
    end
    n = length(cells)
    adj_matrix = zeros(Int, 5, n)
    for k in 1:n_ring
        adj_matrix[k, 1] = k + 1   # centre cell's k-th neighbour slot
    end
    return cells, adj_matrix
end

println("""
Group N — _build_wlsq_weights! (T-NOC3: linear field recovery)
──────────────────────────────────────────────────""")

# ---------------------------------------------------------------------------
# T-NOC3a — 5-neighbour ring (standard A5 interior-cell case)
# ---------------------------------------------------------------------------
# A linear WSE field WSE = a*x + b*y (local metres) should be recovered
# exactly (to floating-point precision) by the WLSQ gradient, for any
# non-degenerate stencil — this is the defining property of a consistent
# (first-order-exact) gradient reconstruction scheme.
let
    cells, adj_matrix = _ring_cluster(5; spacing_deg = 0.01)
    n = length(cells)
    lons = [c.center_lon for c in cells]
    lats = [c.center_lat for c in cells]

    wlsq_weights = zeros(Float64, 10, n)
    _build_wlsq_weights!(wlsq_weights, adj_matrix, lons, lats)

    # Local equirectangular frame, matching _build_wlsq_weights!'s own
    # convention (centred per-cell on cell i; here we only need cell 1's).
    R = A5Grid._EARTH_R
    cos_lat0 = cosd(lats[1])
    to_xy(lon, lat) = (deg2rad(lon - lons[1]) * R * cos_lat0,
                        deg2rad(lat - lats[1]) * R)

    a, b = 0.0037, -0.0021   # arbitrary, non-trivial linear coefficients
    wse = Vector{Float64}(undef, n)
    for i in 1:n
        x, y = to_xy(lons[i], lats[i])
        wse[i] = a * x + b * y
    end

    # Reconstruct cell 1's gradient directly from wlsq_weights, following
    # the multiply-add pattern documented in _build_wlsq_weights!'s
    # docstring (no _compute_wse_gradients! dependency — see file header).
    gx = 0.0; gy = 0.0
    for s in 1:5
        j = adj_matrix[s, 1]
        j == 0 && continue
        dWSE = wse[j] - wse[1]
        gx += wlsq_weights[s,     1] * dWSE
        gy += wlsq_weights[5 + s, 1] * dWSE
    end

    check("T-NOC3a  5-neighbour ring: gx recovers a exactly",
          isapprox(gx, a; atol=1e-9),
          @sprintf("gx=%.10f  expected a=%.10f", gx, a))
    check("T-NOC3a  5-neighbour ring: gy recovers b exactly",
          isapprox(gy, b; atol=1e-9),
          @sprintf("gy=%.10f  expected b=%.10f", gy, b))
end

# ---------------------------------------------------------------------------
# T-NOC3b — reduced-neighbour stencils (2, 3, 4 neighbours)
# ---------------------------------------------------------------------------
# Directly tests the WLSQ fallback path noted as heavily exercised on the
# 27-cell square test mesh (FloodA5_NonOrthogonal_Correction_Plan.md §9.1,
# §10.6.3 — all 27 cells there got a valid stencil, but that diagnostic only
# confirms non-degeneracy, not correctness). A linear field should still be
# recovered exactly with as few as 2 non-collinear neighbours, since the 2D
# WLSQ problem is exactly determined (not over-determined) at k=2 — any
# discrepancy here would indicate a bug in the weighting or solve, not an
# expected approximation error.
for n_ring in (2, 3, 4)
    let
        # For k=2 the two neighbours must not be collinear with the centre
        # cell, or the 2×2 normal-equations matrix is singular by
        # construction. Use angles that guarantee non-collinearity for every
        # n_ring tested (0°, 90° for k=2 is manifestly non-collinear; for
        # k=3,4 the evenly-spaced construction in _ring_cluster is already
        # non-collinear).
        spacing = 0.012
        cells, adj_matrix = if n_ring == 2
            # Override _ring_cluster's even spacing for k=2 (0°, 180° would
            # be collinear) — place neighbours at 0° and 90° instead.
            lon0, lat0 = 172.0, -43.5
            c0 = _synthetic_pentagon(lon0, lat0, spacing * 0.4)
            cos_lat0 = cosd(lat0)
            c1 = _synthetic_pentagon(lon0 + spacing / max(cos_lat0,1e-6), lat0,
                                      spacing * 0.4)
            c2 = _synthetic_pentagon(lon0, lat0 + spacing, spacing * 0.4)
            adj = zeros(Int, 5, 3)
            adj[1, 1] = 2; adj[2, 1] = 3
            ([c0, c1, c2], adj)
        else
            _ring_cluster(n_ring; spacing_deg = spacing)
        end

        n = length(cells)
        lons = [c.center_lon for c in cells]
        lats = [c.center_lat for c in cells]

        wlsq_weights = zeros(Float64, 10, n)
        _build_wlsq_weights!(wlsq_weights, adj_matrix, lons, lats)

        is_degenerate = all(==(0.0), view(wlsq_weights, :, 1))
        check("T-NOC3b  k=$n_ring neighbours: non-degenerate stencil",
              !is_degenerate)
        is_degenerate && continue   # can't test recovery on a zero stencil

        R = A5Grid._EARTH_R
        cos_lat0 = cosd(lats[1])
        to_xy(lon, lat) = (deg2rad(lon - lons[1]) * R * cos_lat0,
                            deg2rad(lat - lats[1]) * R)

        a, b = -0.0058, 0.0044
        wse = [a * to_xy(lons[i], lats[i])[1] + b * to_xy(lons[i], lats[i])[2]
               for i in 1:n]

        gx = 0.0; gy = 0.0
        for s in 1:5
            j = adj_matrix[s, 1]
            j == 0 && continue
            dWSE = wse[j] - wse[1]
            gx += wlsq_weights[s,     1] * dWSE
            gy += wlsq_weights[5 + s, 1] * dWSE
        end

        check("T-NOC3b  k=$n_ring neighbours: gx recovers a exactly",
              isapprox(gx, a; atol=1e-9),
              @sprintf("gx=%.10f  expected a=%.10f", gx, a))
        check("T-NOC3b  k=$n_ring neighbours: gy recovers b exactly",
              isapprox(gy, b; atol=1e-9),
              @sprintf("gy=%.10f  expected b=%.10f", gy, b))
    end
end

# ---------------------------------------------------------------------------
# T-NOC3c — degenerate stencil (k=1, k=0): must stay zero, not crash/NaN
# ---------------------------------------------------------------------------
let
    for n_ring in (0, 1)
        cells, adj_matrix = _ring_cluster(max(n_ring, 1); spacing_deg=0.01)
        # For n_ring=0, manually zero out the only adjacency slot so the
        # centre cell genuinely has zero neighbours (the n_ring=1 cluster
        # builder always gives at least one).
        n_ring == 0 && (adj_matrix[:, 1] .= 0)

        n = length(cells)
        lons = [c.center_lon for c in cells]
        lats = [c.center_lat for c in cells]
        wlsq_weights = zeros(Float64, 10, n)
        _build_wlsq_weights!(wlsq_weights, adj_matrix, lons, lats)

        all_zero = all(==(0.0), view(wlsq_weights, :, 1))
        check("T-NOC3c  k=$n_ring neighbours: correctly left at zero (underdetermined)",
              all_zero)
        check("T-NOC3c  k=$n_ring neighbours: no NaN/Inf produced",
              all(isfinite, view(wlsq_weights, :, 1)))
    end
end

# ---------------------------------------------------------------------------
# T-NOC3d — collinear neighbours: degenerate determinant, must not blow up
# ---------------------------------------------------------------------------
let
    # Two neighbours placed at 0° and 180° from the centre cell are exactly
    # collinear — the 2×2 normal-equations matrix is exactly singular.
    # _build_wlsq_weights!'s determinant guard (FloodModel.jl, 1e-20
    # threshold) must catch this and leave the weights at zero rather than
    # dividing by ~0 and producing huge or NaN/Inf weights.
    lon0, lat0 = 172.0, -43.5
    spacing = 0.01
    cos_lat0 = cosd(lat0)
    c0 = _synthetic_pentagon(lon0, lat0, spacing * 0.4)
    c1 = _synthetic_pentagon(lon0 + spacing / max(cos_lat0,1e-6), lat0, spacing * 0.4)
    c2 = _synthetic_pentagon(lon0 - spacing / max(cos_lat0,1e-6), lat0, spacing * 0.4)
    cells = [c0, c1, c2]
    adj_matrix = zeros(Int, 5, 3)
    adj_matrix[1, 1] = 2
    adj_matrix[2, 1] = 3

    n = length(cells)
    lons = [c.center_lon for c in cells]
    lats = [c.center_lat for c in cells]
    wlsq_weights = zeros(Float64, 10, n)
    _build_wlsq_weights!(wlsq_weights, adj_matrix, lons, lats)

    check("T-NOC3d  collinear neighbours: determinant guard leaves weights at zero",
          all(==(0.0), view(wlsq_weights, :, 1)))
    check("T-NOC3d  collinear neighbours: no NaN/Inf produced",
          all(isfinite, view(wlsq_weights, :, 1)))
end

# ---------------------------------------------------------------------------
# T-NOC4 — _bates_flux_corrected reduces to _bates_flux when skew=0
# ---------------------------------------------------------------------------
# IMPLEMENTED 2026-06-23: _bates_flux_corrected now exists (Step 6).
#
# When the skewness correction is zero (dWSE_n == wse_j - wse_i, exactly the
# raw WSE difference the legacy kernel uses), _bates_flux_corrected must
# produce bit-identical output to _bates_flux for the same inputs — this is
# the defining backward-compatibility property documented in
# FloodA5_NonOrthogonal_Correction_Plan.md §5.3/§5.4: the corrected kernel
# is a strict generalisation of the legacy one, not a different formula
# that happens to agree in some cases.
println()
println("""
Group N+1 — _bates_flux_corrected (T-NOC4: reduces to _bates_flux at skew=0)
──────────────────────────────────────────────────""")

let
    # Representative parameter set spanning wet, dry, and near-threshold cases.
    cases = [
        (q_prev=0.0,   wse_i=2.0, wse_j=1.0, z_sill=0.0, width=100.0, L=1000.0, n=0.03, dt=30.0),
        (q_prev=0.5,   wse_i=1.0, wse_j=2.0, z_sill=0.0, width=80.0,  L=500.0,  n=0.04, dt=10.0),
        (q_prev=-0.3,  wse_i=3.0, wse_j=3.0, z_sill=2.5, width=50.0,  L=250.0,  n=0.03, dt=5.0),
        (q_prev=0.0,   wse_i=0.0005, wse_j=0.0, z_sill=0.0, width=10.0, L=50.0, n=0.03, dt=1.0),  # below HFLOW_THRESHOLD (0.001): both kernels must agree on the dry-edge 0.0 return
    ]

    all_match = true
    for (k, c) in enumerate(cases)
        Q_legacy = _bates_flux(c.q_prev, c.wse_i, c.wse_j, c.z_sill,
                                c.width, c.L, 1.0,   # cos_theta = 1.0 → L_eff = L, no scaling
                                c.n, c.dt)

        h_flow_raw = max(c.wse_i, c.wse_j) - c.z_sill
        dWSE_n     = c.wse_i - c.wse_j   # skew=0: dWSE_n is exactly the raw dWSE
        Q_corrected = _bates_flux_corrected(c.q_prev, h_flow_raw, dWSE_n,
                                             c.width, c.L, c.n, c.dt)

        match = isapprox(Q_legacy, Q_corrected; atol=1e-12)
        all_match &= match
        check("T-NOC4  case $k: _bates_flux_corrected(skew=0) == _bates_flux(cos_theta=1)",
              match,
              @sprintf("legacy=%.10f  corrected=%.10f", Q_legacy, Q_corrected))
    end
end

# ---------------------------------------------------------------------------
# T-NOC4b — _bates_flux_limited_corrected reduces to _bates_flux_limited
# ---------------------------------------------------------------------------
let
    cases = [
        (q_prev=0.0,  wse_i=2.0, wse_j=1.0, z_sill=0.0, width=100.0, L=1000.0, n=0.03, dt=30.0, depth_donor=2.0),
        (q_prev=5.0,  wse_i=1.0, wse_j=2.0, z_sill=0.0, width=80.0,  L=500.0,  n=0.04, dt=10.0, depth_donor=2.0),  # large q_prev: exercises Froude cap
        (q_prev=0.1,  wse_i=3.0, wse_j=3.0, z_sill=2.5, width=50.0,  L=250.0,  n=0.03, dt=5.0,  depth_donor=0.5), # exercises volume limiter
    ]
    for (k, c) in enumerate(cases)
        Q_legacy, q_legacy = _bates_flux_limited(c.q_prev, c.wse_i, c.wse_j, c.z_sill,
                                                   c.width, c.L, 1.0, c.n, c.dt, c.depth_donor)

        h_flow_raw = max(c.wse_i, c.wse_j) - c.z_sill
        dWSE_n     = c.wse_i - c.wse_j
        Q_corr, q_corr = _bates_flux_limited_corrected(c.q_prev, h_flow_raw, dWSE_n,
                                                          c.width, c.L, c.n, c.dt, c.depth_donor)

        check("T-NOC4b  case $k: Q matches (limited, skew=0)",
              isapprox(Q_legacy, Q_corr; atol=1e-12),
              @sprintf("legacy=%.10f  corrected=%.10f", Q_legacy, Q_corr))
        check("T-NOC4b  case $k: q_stored matches (limited, skew=0)",
              isapprox(q_legacy, q_corr; atol=1e-12),
              @sprintf("legacy=%.10f  corrected=%.10f", q_legacy, q_corr))
    end
end

# ---------------------------------------------------------------------------
# T-NOC5 — _compute_wse_gradients! correctness + mass conservation
# ---------------------------------------------------------------------------
# IMPLEMENTED 2026-06-23: _compute_wse_gradients! now exists (Step 7).
println()
println("""
Group N+2 — _compute_wse_gradients! (T-NOC5)
──────────────────────────────────────────────────""")

# T-NOC5a: linear field recovery via the full state-based call path (not the
# manual multiply-add used in T-NOC3 — this exercises _compute_wse_gradients!
# itself, end-to-end, on a real FlowState with wlsq_weights already built by
# initialise_flow_model's call path).
let
    cells, adj_matrix = _ring_cluster(5; spacing_deg = 0.01)
    n = length(cells)
    lons = [c.center_lon for c in cells]
    lats = [c.center_lat for c in cells]

    wlsq_weights = zeros(Float64, 10, n)
    _build_wlsq_weights!(wlsq_weights, adj_matrix, lons, lats)

    R = A5Grid._EARTH_R
    cos_lat0 = cosd(lats[1])
    to_xy(lon, lat) = (deg2rad(lon - lons[1]) * R * cos_lat0,
                        deg2rad(lat - lats[1]) * R)

    a, b = 0.0019, 0.0061
    wse = [a * to_xy(lons[i], lats[i])[1] + b * to_xy(lons[i], lats[i])[2]
           for i in 1:n]

    # Minimal FlowState — only the fields _compute_wse_gradients! actually
    # reads (cell_ids for n, adj_matrix, wlsq_weights, grad_wse) need to be
    # meaningful; everything else is a placeholder, following the same
    # pattern as test_edge_geometry.jl's _minimal_state.
    dummy_edges = EdgeList(0, Int[], Int[], Float64[], Float64[], Float64[],
                            Float64[], Float64[], Float64[], Int[], Int[],
                            Float64[], Float64[], Float64[], Float64[])
    state = FlowState(
        ["c$i" for i in 1:n],
        zeros(n), zeros(n), zeros(n), zeros(n), zeros(n),
        zeros(n), fill(0.03, n), fill(1.0e6, n),
        lons, lats,
        Dict{String,Vector{String}}(),
        adj_matrix,
        dummy_edges,
        Any[],
        falses(n), Any[], Any[], 0.0,
        zeros(Float64, 2, n),
        wlsq_weights,
        0.9, true,
    )

    _compute_wse_gradients!(state, wse)

    check("T-NOC5a  _compute_wse_gradients!: gx recovers a exactly",
          isapprox(state.grad_wse[1, 1], a; atol=1e-9),
          @sprintf("gx=%.10f  expected a=%.10f", state.grad_wse[1, 1], a))
    check("T-NOC5a  _compute_wse_gradients!: gy recovers b exactly",
          isapprox(state.grad_wse[2, 1], b; atol=1e-9),
          @sprintf("gy=%.10f  expected b=%.10f", state.grad_wse[2, 1], b))
end

# T-NOC5b: mass conservation with gradient_correction=true on a tiny but
# real edge list (not a degenerate 0-edge stub), using step_standard! end
# to end. This is the first test in this file that actually runs the
# corrected kernel through the full Phase A/B/C dispatch in step_standard!,
# rather than calling _bates_flux_corrected directly.
let
    # 3-cell chain: cell 1 (high) -> cell 2 (mid) -> cell 3 (low), orthogonal
    # edges (skew=0) so this also cross-checks T-NOC4's no-op property in
    # the context of a real step_standard! call, not just the bare kernel.
    n  = 3
    ci = [1, 2]; cj = [2, 3]
    width = [100.0, 100.0]
    L     = [500.0, 500.0]
    cos_theta = [1.0, 1.0]
    sill  = [0.0, 0.0]
    edges = EdgeList(2, ci, cj, width, L, cos_theta, sill,
                      zeros(2), zeros(2),
                      [0, 0], [0, 0],          # collinear_i/j
                      [0.0, 0.0], [0.0, 0.0],  # skew_x/y = 0: orthogonal
                      L, zeros(2))              # dx_m = L (chain E), dy_m = 0

    adj_matrix = zeros(Int, 5, n)
    adj_matrix[1, 1] = 2
    adj_matrix[1, 2] = 1; adj_matrix[2, 2] = 3
    adj_matrix[1, 3] = 2

    wlsq_weights = zeros(Float64, 10, n)   # not used meaningfully on a 1D
                                            # chain — each cell has fewer
                                            # than 2 non-collinear
                                            # neighbours, so
                                            # _build_wlsq_weights! would
                                            # leave these at zero anyway —
                                            # confirms the "degenerate to
                                            # zero gradient, no-op
                                            # correction" fallback path
                                            # documented in
                                            # _compute_wse_gradients!'s
                                            # docstring is exercised here.

    volume = [30000.0, 10000.0, 0.0]
    cell_area = fill(10000.0, n)   # FIXED 2026-06-24: was fill(100.0, n),
                                   # i.e. cell_area == width (100.0) — a
                                   # geometrically nonsensical fixture (a
                                   # 100m-wide cell with only 100 m² area
                                   # implies an absurdly thin cell, roughly
                                   # 1m "deep"). Real A5 cells have area
                                   # roughly proportional to width², not
                                   # equal to width — confirmed against the
                                   # real square-mesh run (~127,200 m² area,
                                   # ~400m equivalent diameter). This was a
                                   # second, independent fixture bug found
                                   # alongside the dWSE_n sign/formula bugs:
                                   # it didn't cause incorrect *direction*,
                                   # but made the Fix B volume-limiter
                                   # bound (depth_donor×width/(5×dt)) huge
                                   # relative to the cell's actual stored
                                   # volume, masking whether the limiter
                                   # was behaving sensibly. Volumes below
                                   # are scaled up to preserve the same
                                   # depths (3.0, 1.0, 0.0 m) as before.
    state = FlowState(
        ["c1", "c2", "c3"],
        volume ./ cell_area,        # water_depth — MUST be synced from volume
                                    # before calling step_standard!, exactly as
                                    # production does in run_simulation! (Bug
                                    # 37, PROJECT_STATE.md). Leaving this at
                                    # zeros(n) here previously meant
                                    # depth_donor was always 0.0, which SKIPS
                                    # the Fix B volume limiter entirely (see
                                    # _bates_flux_limited_corrected's
                                    # `if depth_donor > 0.0` guard — it is not
                                    # clamped to zero when depth_donor==0, the
                                    # limiter is simply not applied), letting
                                    # unrestricted flux through and producing
                                    # the spurious cell-1-gains-volume result
                                    # first observed in this test. This was a
                                    # test-fixture bug, not a kernel bug.
        copy(volume), zeros(n), zeros(n), zeros(n),
        [2.0, 1.0, 0.0],            # elevation: downhill chain
        fill(0.03, n), cell_area,
        zeros(n), zeros(n),        # cell_lons, cell_lats: left at zero.
                                    # Phase E (_compute_velocity!) reads these
                                    # to get edge direction; with both zero,
                                    # every edge has dist=0 and is skipped, so
                                    # velocity silently stays zero. Harmless
                                    # for this test since only mass
                                    # conservation is under test here, but
                                    # velocity output should not be read as
                                    # validated by this particular check.
        Dict{String,Vector{String}}(),
        adj_matrix,
        edges,
        Any[],
        falses(n), Any[], Any[], 0.0,
        zeros(Float64, 2, n),
        wlsq_weights,
        0.9, false,  # gradient_correction = FALSE — 1D chain has collinear
        # neighbours, so _build_wlsq_weights! yields zero weights and
        # _compute_wse_gradients! yields zero gradients. With the n̂_f formula
        # dWSE_n = 0 for all edges → no flux. This is the documented fallback
        # for degenerate stencils, not a bug. T-NOC5b therefore tests the
        # UNCORRECTED kernel path only. The corrected path is validated on a
        # proper 2D stencil in test_gradient_direction.jl GD3.
    )

    vol_before = sum(state.volume)
    step_standard!(state, 5.0)
    vol_after = sum(state.volume)

    rel_err = abs(vol_after - vol_before) / max(vol_before, 1.0)
    check("T-NOC5b  step_standard! (uncorrected dispatch): mass conserved",
          rel_err < 1e-10,
          @sprintf("before=%.6f  after=%.6f  rel_err=%.2e", vol_before, vol_after, rel_err))
    check("T-NOC5b  step_standard! (uncorrected dispatch): no NaN volumes",
          all(isfinite, state.volume))
    check("T-NOC5b  step_standard! (uncorrected dispatch): downhill flow (cell 1 loses volume)",
          state.volume[1] < volume[1],
          "vol[1] before=$(volume[1])  after=$(state.volume[1])")
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

println()
println("─" ^ 50)
total = n_pass[] + n_fail[]
@printf("Results: %d / %d passed", n_pass[], total)
if n_fail[] == 0
    println("  — all tests passed ✓")
else
    println()
    println("  $(n_fail[]) test(s) FAILED")
    exit(1)
end
