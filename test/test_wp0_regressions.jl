# test/test_wp0_regressions.jl
# ----------------------------
# Regression tests for the WP0 correctness fixes (branch wp0-correctness,
# 2026-10 — see FloodA5_Review_and_ForwardPlan_2026-10.md §4.4) plus the
# lake-at-rest (well-balancing) check the review asked for.
#
#   W1  _cfl_dt applies the SGS terrain-range depth cap (was never applied:
#       `method == SGSFlow` compared an instance with a type)
#   W2  _adjacency_slot searches every slot of the mesh (was capped at 5)
#   W3  eta_at_vertex returns NaN for a :degenerate vertex record (was 0.0)
#   W4  step_standard! Phase A is a Jacobi update: every edge's new q is
#       computed from the START-OF-STEP momentum (was a data race)
#   W5  _compute_velocity! recovers a uniform velocity field exactly (m/s)
#   W6  Lake at rest: a still, non-flat water surface stays still under every
#       flux/momentum option (legacy, WLSQ-corrected, diamond, cell momentum)
#
# Uses a small Cartesian mesh generated in-process — no Python mesh bridge,
# no data files.  (pya5 must still be importable, because A5Grid's __init__
# imports it.)
#
# Run standalone (headless — no GLMakie needed):
#   FLOODA5_NO_VIS=1 julia --project=. test/test_wp0_regressions.jl
# PowerShell:
#   $env:FLOODA5_NO_VIS="1"; julia --project=. test\test_wp0_regressions.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test

# ---------------------------------------------------------------------------
# Fixture: an nr × nc Cartesian mesh around (0°, 0°) with ~dx m cells.
# ---------------------------------------------------------------------------
function _cart_mesh(; nr::Int = 9, nc::Int = 9, dx::Float64 = 20.0)
    R    = A5Grid._EARTH_R
    dlat = rad2deg(dx / R)
    dlon = dlat                      # lat_ref = 0 → cos = 1
    # Bounds are set a hair inside nr×nc cells so ceil() gives exactly nr, nc.
    return generate_cartesian_mesh_bounds(0.0, 0.0,
                                          (nc - 0.5) * dlon, (nr - 0.5) * dlat,
                                          dx; lat_ref = 0.0)
end

function _closed!(state::FlowState)
    for k in eachindex(state.ghost_cell_bc)
        state.ghost_cell_bc[k] = Closed
    end
    return state
end

_wse(state, i) = state.elevation[i] + state.volume[i] / max(state.cell_area[i], 1.0)

function _set_depths!(state::FlowState, h::Vector{Float64})
    for i in eachindex(h)
        state.volume[i]      = h[i] * state.cell_area[i]
        state.water_depth[i] = h[i]
    end
end

# ---------------------------------------------------------------------------
@testset "W1 — _cfl_dt applies the SGS depth cap" begin
    mesh  = _cart_mesh()
    state = initialise_flow_model(mesh, StandardFlow())
    n     = length(state.cell_ids)
    nb    = 10
    # Minimal per-cell SGS tables with a 0.5 m terrain range (z_min 0, z_max 0.5).
    tbls = Any[]
    for i in 1:n
        elev = collect(range(0.0, 0.5; length = nb))
        A    = state.cell_area[i]
        push!(tbls, SGSTable(elev, A .* elev, fill(A, nb), A, 0.0, 0.5,
                             zeros(nb, 5), zeros(nb, 5)))
    end
    state.sgs_tables = tbls
    fill!(state.water_depth, 4.0)          # every cell "overfull": 4 m ≫ 0.5 m range

    dt_std = _cfl_dt(state, StandardFlow())
    dt_sgs = _cfl_dt(state, SGSFlow())
    # Standard uses h = 4 m; SGS must cap at the 0.5 m terrain range → dt
    # larger by √(4/0.5) = √8.
    @test dt_sgs > dt_std
    @test isapprox(dt_sgs / dt_std, sqrt(8.0); rtol = 1e-10)
end

# ---------------------------------------------------------------------------
@testset "W2 — _adjacency_slot searches every slot" begin
    adj = zeros(Int, 6, 2)
    adj[6, 1] = 2                          # neighbour in slot 6 of a 6-slot mesh
    adj[1, 2] = 1
    @test _adjacency_slot(adj, 1, 2) == 6  # default max_nb = size(adj, 1)
    @test _adjacency_slot(adj, 1, 2, 6) == 6
    @test _adjacency_slot(adj, 2, 1) == 1
end

# ---------------------------------------------------------------------------
@testset "W3 — degenerate diamond vertex gives NaN" begin
    rec_deg = VertexRecord([1, 2, 3], zeros(3), 1e12, :degenerate)
    rec_ok  = VertexRecord([1, 2, 3], [0.2, 0.3, 0.5], 10.0, :none)
    eta     = [10.0, 20.0, 30.0]
    @test isnan(eta_at_vertex(rec_deg, eta))
    @test eta_at_vertex(rec_ok, eta) ≈ 0.2*10 + 0.3*20 + 0.5*30
end

# ---------------------------------------------------------------------------
@testset "W4 — Phase A momentum update is Jacobi (race fix)" begin
    mesh  = _cart_mesh(nr = 11, nc = 11)
    state = _closed!(initialise_flow_model(mesh, StandardFlow();
                                           q_centre_theta = 0.9))
    n     = length(state.cell_ids)
    edges = state.edges
    # A non-trivial depth field (a smooth mound) and non-trivial prior momentum.
    xs = [c.center_lon for c in mesh.cells]; ys = [c.center_lat for c in mesh.cells]
    dd = rad2deg(20.0 / A5Grid._EARTH_R)      # one cell, in degrees
    cx = sum(xs) / n; cy = sum(ys) / n
    h  = [0.05 + 0.5 * exp(-((xs[i]-cx)^2 + (ys[i]-cy)^2) / (2 * (2dd)^2)) for i in 1:n]
    _set_depths!(state, h)
    for e in 1:edges.n_edges
        edges.flux[e] = 0.01 * sin(0.7 * e)    # arbitrary, deterministic, sign-varying
    end
    flux0 = copy(edges.flux)
    wse0  = [_wse(state, i) for i in 1:n]
    dep0  = copy(state.water_depth)
    dt    = 0.5

    step_standard!(state, dt)

    max_err = 0.0
    for e in 1:edges.n_edges
        ci = edges.cell_i[e]; cj = edges.cell_j[e]
        q_prev = _q_centred(flux0, e, edges.collinear_i[e], edges.collinear_j[e],
                            state.q_centre_theta)
        depth_donor = wse0[ci] >= wse0[cj] ? dep0[ci] : dep0[cj]
        _, q_exp = _bates_flux_limited(q_prev, wse0[ci], wse0[cj], edges.sill[e],
                                       edges.width[e], edges.L[e], edges.cos_theta[e],
                                       min(state.manning_n[ci], state.manning_n[cj]),
                                       dt, depth_donor;
                                       vol_limit_divisor = Float64(state.max_neighbours))
        max_err = max(max_err, abs(edges.flux[e] - q_exp))
    end
    # Same kernel, same inputs → identical to rounding. Under the old race the
    # error was O(1e-3) relative even single-threaded (in-order Gauss–Seidel).
    @test max_err < 1e-13
end

# ---------------------------------------------------------------------------
@testset "W5 — _compute_velocity! recovers a uniform field (m/s)" begin
    mesh  = _cart_mesh()
    state = initialise_flow_model(mesh, StandardFlow())
    n     = length(state.cell_ids)
    edges = state.edges
    qx, qy, h = 0.30, -0.12, 0.5              # m²/s, m²/s, m
    fill!(state.water_depth, h)
    for e in 1:edges.n_edges
        # Observation convention (Phase F / _build_mom_weights!): q·n̂ = −flux[e]
        edges.flux[e] = -(qx * edges.nf_x[e] + qy * edges.nf_y[e])
    end
    fill!(edges.flux_Q, 0.0)
    _compute_velocity!(state)
    interior = [i for i in 1:n if count(!iszero, view(state.adj_matrix, :, i)) == 4]
    @test !isempty(interior)
    @test maximum(abs(state.vel_u[i] - qx / h) for i in interior) < 1e-9
    @test maximum(abs(state.vel_v[i] - qy / h) for i in interior) < 1e-9
    @test maximum(abs(state.velocity[i] - hypot(qx, qy) / h) for i in interior) < 1e-9
end

# ---------------------------------------------------------------------------
# W6 — lake at rest. Bed varies cell to cell (smooth + checker component);
# every cell is submerged under a flat WSE of 2.0 m. Any spurious flux means
# the scheme is not well-balanced for that option.
function _lake_case(; kwargs...)
    mesh = _cart_mesh(nr = 9, nc = 9)
    n    = length(mesh.cells)
    mesh.static_vars["elevation"] =
        [0.6 + 0.3 * sin(0.9 * i) + 0.2 * isodd(i) for i in 1:n]   # 0.3 … 1.1 m
    state = _closed!(initialise_flow_model(mesh, StandardFlow(); kwargs...))
    wse_target = 2.0
    _set_depths!(state, [wse_target - state.elevation[i] for i in 1:n])
    for _ in 1:200
        step_standard!(state, 1.0)
    end
    max_dev  = maximum(abs(_wse(state, i) - wse_target) for i in 1:n)
    max_flux = maximum(abs, state.edges.flux)
    return max_dev, max_flux
end

@testset "W6 — lake at rest ($label)" for (label, kw) in [
        ("legacy, edge momentum",           (;)),
        ("legacy, θ = 1",                   (; q_centre_theta = 1.0)),
        ("WLSQ corrected",                  (; gradient_correction = true)),
        ("diamond",                         (; gradient_correction = true,
                                               face_flux_method = :diamond)),
        ("cell momentum",                   (; momentum_model = :cell)),
        ("diamond + cell momentum",         (; gradient_correction = true,
                                               face_flux_method = :diamond,
                                               momentum_model = :cell)),
    ]
    max_dev, max_flux = _lake_case(; kw...)
    @test max_dev  < 1e-10
    @test max_flux < 1e-10
end

println("\nWP0 regression tests complete.")
