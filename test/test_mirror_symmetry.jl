#!/usr/bin/env julia
# ============================================================================
# test_mirror_symmetry.jl
#
# Regression test for the north/south symmetry of the corrected flux
# formula. Originally written (2026-08-17) as a diagnostic against the
# "purely gradient-driven" dWSE_n formula, which passed this test but still
# produced a severe real-mesh directional instability when combined with
# cell-vector momentum -- the root cause turned out to be structural (no
# locally-anchored term in the driving force at all), not a per-step
# symmetry bug, so this script alone could not have caught it.
#
# Updated (2026-08-18) to validate the REPLACEMENT formula: a standard
# OpenFOAM-style local-orthogonal + non-orthogonal-correction decomposition,
#   dWSE_n = c·(wse_ci − wse_cj) − alpha·L·(∇WSE_f·V̂)
# where the dominant term is the direct, locally-measured WSE difference
# between the two adjacent cells, and only the smaller V̂ (tangential)
# remainder comes from the WLSQ-reconstructed gradient. See FloodModel.jl
# step_standard! Phase A for the full rationale.
#
# This script builds a hand-crafted 6-cell "star" mesh: one centre cell C
# with 5 neighbours (E, NE, SE, NW, SW) placed so that the NE/SE pair and
# NW/SW pair are EXACT mirror images across the east-west (longitude) axis
# -- same geometry magnitudes, y-components negated. The WSE field assigned
# to every cell is a pure function of longitude (WSE = WSE0 - slope*x), i.e.
# by construction it has NO north-south dependence at all.
#
# Given perfectly mirror-symmetric geometry and a perfectly y-independent
# input field, the physically correct answer is UNAMBIGUOUS: the NE and SE
# edges must produce identical dWSE_n and identical flux magnitude (and same
# for NW/SW). Any measured difference is not physics, not accumulation, not
# resolution -- it is a bug in the formula or its implementation.
#
# Kept in the suite going forward as a cheap regression guard: unlike the
# old formula, this one's dominant term (wse_ci-wse_cj) is trivially
# mirror-symmetric by construction, so MS2/MS3 passing is close to a
# tautology for THAT term -- what remains genuinely informative is MS1 (is
# the WLSQ gradient reconstruction itself symmetric on an irregular
# stencil?) and MS4 (does the Phase F momentum reconstruction, a separate
# WLSQ operator, stay symmetric?). A future change to either reconstruction
# should keep passing both.
#
# This script does NOT call step_standard! or run a real simulation. It
# calls the individual building blocks (_build_wlsq_weights!,
# _compute_wse_gradients!, _bates_flux_limited_corrected, the Phase F
# reconstruction formula) directly, using the exact same code paths and
# formulas as the production Phase A/F implementation, replicated inline
# where step_standard! only exists as a threaded loop body.
#
# Run: julia --project=. test\test_mirror_symmetry.jl
# ============================================================================

using Test
using LinearAlgebra

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

println("=" ^ 76)
println("test_mirror_symmetry.jl")
println("Diagnostic: does the corrected kernel preserve north/south symmetry")
println("on a hand-built, exactly-mirror-symmetric synthetic mesh?")
println("=" ^ 76)

# ── 1. Build the synthetic 6-cell star ─────────────────────────────────────
# Cells: 1=C (centre), 2=E, 3=NE, 4=SE, 5=NW, 6=SW
# Placed via metre offsets from C, converted to lon/lat with the SAME local
# equirectangular projection _build_wlsq_weights! uses internally -- this
# guarantees exact symmetry regardless of floating-point rounding, since the
# NE/SE pair share an identical dx and exactly negated dy by construction.

const R_EARTH = 6_371_000.0
lat0, lon0 = 51.0, -0.02
coslat0 = cosd(lat0)
to_lonlat(dx, dy) = (lon0 + rad2deg(dx / (R_EARTH * coslat0)),
                      lat0 + rad2deg(dy / R_EARTH))

# (dx, dy) offsets in metres -- an irregular-but-mirror-symmetric pentagon
# star, loosely modelled on real A5 cell proportions.
offsets = Dict(
    :E  => ( 400.0,    0.0),
    :NE => ( 150.0,  380.0),
    :SE => ( 150.0, -380.0),
    :NW => (-350.0,  300.0),
    :SW => (-350.0, -300.0),
)

idx = Dict(:C => 1, :E => 2, :NE => 3, :SE => 4, :NW => 5, :SW => 6)
n = 6

cell_lons = zeros(Float64, n)
cell_lats = zeros(Float64, n)
cell_lons[idx[:C]], cell_lats[idx[:C]] = lon0, lat0
for k in (:E, :NE, :SE, :NW, :SW)
    dx, dy = offsets[k]
    cell_lons[idx[k]], cell_lats[idx[k]] = to_lonlat(dx, dy)
end

# ── 2. WSE field: pure function of longitude (no y-dependence at all) ──────
# WSE(x) = WSE0 - slope*x, matching the real 0.1%-slope planar test's
# west-to-east downhill convention. Implemented via elevation = target WSE,
# volume = 0, so wse_all[i] = elevation[i] + volume[i]/area[i] = elevation[i].
WSE0, slope = 5.0, 0.001
target_wse = Dict(k => WSE0 - slope * offsets[k][1] for k in (:E, :NE, :SE, :NW, :SW))
target_wse[:C] = WSE0

elevation = zeros(Float64, n)
for (k, i) in idx
    elevation[i] = target_wse[k]
end
volume     = zeros(Float64, n)
cell_area  = fill(1.0, n)
water_depth = zeros(Float64, n)   # dry bed -> depth_donor=0 -> Fix B limiter inert,
                                   # isolating the driving-term/Froude-limiter path
manning_n  = fill(0.03, n)

println("\nInput WSE (function of longitude only, by construction):")
for k in (:C, :E, :NE, :SE, :NW, :SW)
    println("  $k: wse=$(elevation[idx[k]])")
end
@assert elevation[idx[:NE]] == elevation[idx[:SE]] "fixture bug: NE/SE WSE not equal"
@assert elevation[idx[:NW]] == elevation[idx[:SW]] "fixture bug: NW/SW WSE not equal"

# ── 3. adj_matrix: only C has neighbours (slots 1..5 = E,NE,SE,NW,SW) ──────
adj_matrix = zeros(Int, N_SIDES, n)
adj_matrix[1, idx[:C]] = idx[:E]
adj_matrix[2, idx[:C]] = idx[:NE]
adj_matrix[3, idx[:C]] = idx[:SE]
adj_matrix[4, idx[:C]] = idx[:NW]
adj_matrix[5, idx[:C]] = idx[:SW]
# Neighbour cells deliberately have NO adjacency of their own (isolated) --
# this gives them grad_wse=(0,0) trivially (underdetermined stencil, the
# same safe fallback _build_wlsq_weights! already uses), which means the
# Phase A face-averaged gradient 0.5*(grad[ci]+grad[cj]) reduces to
# 0.5*grad[C] for EVERY one of C's five edges -- identical for all five,
# so this cannot itself be a source of NE/SE asymmetry. Any asymmetry we
# see must come from (a) grad_wse[:,C] itself having a spurious non-zero
# y-component, or (b) the per-edge formula/kernel treating the mirrored
# geometry inconsistently.

# ── 4. Hand-built EdgeList: 5 star edges, ci=C for all ─────────────────────
# Field order (17 total): n_edges, cell_i, cell_j, width, L, cos_theta, sill,
# flux, flux_Q, collinear_i, collinear_j, skew_x, skew_y, dx_m, dy_m, nf_x, nf_y
edge_order = [:E, :NE, :SE, :NW, :SW]
ne = 5
cell_i = fill(idx[:C], ne)
cell_j = [idx[k] for k in edge_order]

width     = Dict(:E=>35.0, :NE=>30.0, :SE=>30.0, :NW=>28.0, :SW=>28.0)
cos_theta = Dict(:E=>1.00, :NE=>0.85, :SE=>0.85, :NW=>0.80, :SW=>0.80)
# V̂ (skew_x, skew_y) -- the WLSQ correction DIRECTION vector (not a
# positional offset -- see FloodA5_NonOrthogonal_Correction_Plan.md §10.6.7).
# Constructed exactly mirror-symmetric: same x-component, negated y-component.
skewv = Dict(:E=>(0.00,0.00), :NE=>(0.15,0.30), :SE=>(0.15,-0.30),
             :NW=>(-0.10,0.25), :SW=>(-0.10,-0.25))

widths_v = [width[k] for k in edge_order]
ct_v     = [cos_theta[k] for k in edge_order]
sill_v   = [min(elevation[idx[:C]], elevation[idx[k]]) - 1.0 for k in edge_order]
dx_v     = [offsets[k][1] for k in edge_order]
dy_v     = [offsets[k][2] for k in edge_order]
L_v      = [sqrt(dx_v[m]^2 + dy_v[m]^2) for m in 1:ne]
skx_v    = [skewv[k][1] for k in edge_order]
sky_v    = [skewv[k][2] for k in edge_order]
# nf = cos_theta * d̂ + V̂ (matches the real _build_edge_list derivation
# described in the cell-momentum handover: "nf_x = cos_theta × dx_m/L + skew_x")
nfx_v    = [ct_v[m] * dx_v[m] / L_v[m] + skx_v[m] for m in 1:ne]
nfy_v    = [ct_v[m] * dy_v[m] / L_v[m] + sky_v[m] for m in 1:ne]

edges = EdgeList(
    ne, cell_i, cell_j,
    widths_v, L_v, ct_v, sill_v,
    zeros(Float64, ne),           # flux (q_prev = 0 -- isolates the driving term)
    zeros(Float64, ne),           # flux_Q
    zeros(Int, ne), zeros(Int, ne),  # collinear_i, collinear_j (unused here)
    skx_v, sky_v,
    dx_v, dy_v,
    nfx_v, nfy_v,
)

println("\nEdge geometry (by construction, NE/SE and NW/SW are exact mirrors):")
for (m, k) in enumerate(edge_order)
    println("  $k: dx=$(dx_v[m]) dy=$(dy_v[m]) L=$(round(L_v[m],digits=3)) " *
            "cos_theta=$(ct_v[m]) skew=($(skx_v[m]),$(sky_v[m])) " *
            "nf=($(round(nfx_v[m],digits=4)),$(round(nfy_v[m],digits=4))) sill=$(sill_v[m])")
end

# ── 5. Minimal FlowState (only fields _compute_wse_gradients! touches) ─────
cell_edge_index = zeros(Int, N_SIDES, n)
_build_cell_edge_index!(cell_edge_index, adj_matrix, edges, n)
mom_weights = zeros(Float64, 10, n)
_build_mom_weights!(mom_weights, cell_edge_index, edges, n)

wlsq_weights = zeros(Float64, 10, n)
_build_wlsq_weights!(wlsq_weights, adj_matrix, cell_lons, cell_lats)

state = FlowState(
    ["cell$i" for i in 1:n],
    copy(water_depth), copy(volume),
    zeros(n), zeros(n), zeros(n),          # velocity, vel_u, vel_v
    copy(elevation), copy(manning_n), copy(cell_area),
    copy(cell_lons), copy(cell_lats),
    Dict{String,Vector{String}}(),
    adj_matrix, edges, Any[],
    falses(n), Any[], Any[], 0.0,
    zeros(Float64, 2, n), wlsq_weights,
    0.9, true, 1.0,                        # q_centre_theta, gradient_correction, alpha
    zeros(Float64, n), zeros(Float64, n),  # qvec_u, qvec_v
    cell_edge_index, mom_weights,
    :cell,
    :legacy, nothing,   # face_flux_method, diamond_table (Phase C)
)

# ── 6. Run the real _compute_wse_gradients! and check grad_wse[2,C] ────────
wse_all = copy(elevation)   # volume=0, so wse_all == elevation exactly
_compute_wse_gradients!(state, wse_all)

gx_C, gy_C = state.grad_wse[1, idx[:C]], state.grad_wse[2, idx[:C]]
println("\n--- Step A: _build_wlsq_weights! / _compute_wse_gradients! ---")
println("  grad_wse[:,C] = ($gx_C, $gy_C)")
println("  (gy should be ~0: the input WSE field has zero y-dependence and")
println("   the stencil geometry is exactly mirror-symmetric about y=0)")

@testset "MS1 -- WLSQ gradient has zero y-component on symmetric stencil" begin
    @test isapprox(gy_C, 0.0, atol=1e-9)
end

# ── 7. Replicate the exact Phase A dWSE_n formula for each edge ────────────
# (2026-08-18 formula: local-orthogonal term uses the direct wse_ci-wse_cj
# difference; only the V̂ non-orthogonal correction comes from the
# WLSQ-reconstructed gradient. See FloodModel.jl step_standard! Phase A.)
alpha = state.gradient_correction_alpha
dWSE_n = zeros(Float64, ne)
for m in 1:ne
    j = cell_j[m]
    gx_f = 0.5 * (state.grad_wse[1, idx[:C]] + state.grad_wse[1, j])
    gy_f = 0.5 * (state.grad_wse[2, idx[:C]] + state.grad_wse[2, j])
    Vhat_dot = gx_f * skx_v[m] + gy_f * sky_v[m]
    dWSE_n[m] = ct_v[m] * (elevation[idx[:C]] - elevation[j]) -
                alpha * L_v[m] * Vhat_dot
end

println("\n--- Step B: dWSE_n per edge (Phase A formula, replicated) ---")
for (m, k) in enumerate(edge_order)
    println("  $k: dWSE_n = $(dWSE_n[m])")
end

i_NE, i_SE = findfirst(==(:NE), edge_order), findfirst(==(:SE), edge_order)
i_NW, i_SW = findfirst(==(:NW), edge_order), findfirst(==(:SW), edge_order)

@testset "MS2 -- dWSE_n identical for mirror-image edge pairs" begin
    @test isapprox(dWSE_n[i_NE], dWSE_n[i_SE], atol=1e-9)
    @test isapprox(dWSE_n[i_NW], dWSE_n[i_SW], atol=1e-9)
end

# ── 8. Call the real flux kernel per edge, compare NE/SE and NW/SW ─────────
# Sweep two very different dt values, since the real-mesh finding was that
# asymmetry magnitude is dt-DEPENDENT (V-shaped, not monotonic) -- if a bug
# exists here it may only be visible (or may change magnitude) at certain dt.
println("\n--- Step C: _bates_flux_limited_corrected per edge ---")
for dt in (10.0, 0.5)
    println("\n  dt = $dt s:")
    Q = zeros(Float64, ne)
    qs = zeros(Float64, ne)
    for m in 1:ne
        h_flow = max(elevation[idx[:C]], elevation[cell_j[m]]) - sill_v[m]
        depth_donor = elevation[idx[:C]] >= elevation[cell_j[m]] ?
                      water_depth[idx[:C]] : water_depth[cell_j[m]]
        Q[m], qs[m] = _bates_flux_limited_corrected(
            0.0, h_flow, dWSE_n[m], widths_v[m], L_v[m], 0.03, dt, depth_donor)
    end
    for (m, k) in enumerate(edge_order)
        println("    $k: Q=$(round(Q[m],digits=6))  q_stored=$(round(qs[m],digits=6))")
    end

    @testset "MS3 (dt=$dt) -- flux magnitude identical for mirror-image pairs" begin
        @test isapprox(abs(Q[i_NE]), abs(Q[i_SE]), rtol=1e-9)
        @test isapprox(abs(Q[i_NW]), abs(Q[i_SW]), rtol=1e-9)
        # Signs should also be OPPOSITE about the axis in the sense that the
        # underlying q_prev_eff-driven momentum term (tested separately,
        # Step D) sees this as a symmetric input -- here we only assert
        # magnitude, since q_prev=0 removes the momentum term from play.
    end

    # ── 9. Phase F reconstruction: qvec_v[C] should come out ~0 ────────────
    edges.flux .= qs   # write q_stored back, as Phase A would before Phase F
    qu, qv = 0.0, 0.0
    for s in 1:N_SIDES
        e = state.cell_edge_index[s, idx[:C]]
        e == 0 && continue
        obs = -edges.flux[e]
        qu += state.mom_weights[s,     idx[:C]] * obs
        qv += state.mom_weights[5 + s, idx[:C]] * obs
    end
    println("    Phase F reconstruction at C: qvec_u=$(round(qu,digits=6))  " *
            "qvec_v=$(round(qv,digits=6))")

    @testset "MS4 (dt=$dt) -- Phase F qvec_v ~0 at centre cell" begin
        @test isapprox(qv, 0.0, atol=1e-9)
    end
end

println("\n" * "=" ^ 76)
println("If MS1/MS2/MS3/MS4 all PASS at both dt values: the new local-anchor")
println("formula is symmetric in isolation on this synthetic mesh, as")
println("expected -- this is a necessary but not sufficient condition. The")
println("real test is the res-18/res-16 planar-symmetry dt sweep on the")
println("actual simulation: does removing the purely-gradient-driven design")
println("(the missing local anchor identified 2026-08-18) fix the real-mesh")
println("north/south instability, or at least remove its V-shaped, worse-at-")
println("small-dt tail?")
println()
println("If any assertion FAILS: something is still wrong with the new")
println("formula's implementation even in this trivial symmetric case -- the")
println("printed values above pinpoint exactly which term diverges.")
println("=" ^ 76)
