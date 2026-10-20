# test/test_wp2_momentum_limiters.jl
# ----------------------------------
# Tests for WP2 (FloodA5_Review_and_ForwardPlan_2026-10.md §4.2, §9 WP2):
#
#   B1  --momentum-model blend: θ = 0 reproduces :cell exactly and θ = 1
#       reproduces :edge with Q-centring off (θ = 1) exactly
#   B2  Frozen-coefficient isotropy with blend at the default θ = 0.9
#       (vector friction, diamond): no deflection
#   P1  _positivity_scale!: scales every outflow of an overdrawn donor by
#       V/Σout, leaves other edges alone, writes back the stored momentum
#   P2  --limiter positivity|hybrid in a point-source run: exact mass balance,
#       no negative volume
#   C1  --cfl-mode local: ℓ = 4A/P equals dx on every Cartesian cell
#       (interior and edge cells), dt = 0.7·dx/√(g·h_max)
#   D1  diamond_stencil_wet, and diamond-wetdry == diamond when all cells are wet
#   G1  Cartesian ghost edges: one per missing side (phantom-edge fix)
#   L1  Lake at rest with the WP2 options on
#
# Uses an in-process Cartesian mesh — no data files.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_wp2_momentum_limiters.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
include(joinpath(@__DIR__, "frozen_isotropy_lib.jl"))

using Test

function _cart_mesh(; nr::Int = 9, nc::Int = 9, dx::Float64 = 20.0)
    dlat = rad2deg(dx / A5Grid._EARTH_R)
    return generate_cartesian_mesh_bounds(0.0, 0.0,
                                          (nc - 0.5) * dlat, (nr - 0.5) * dlat,
                                          dx; lat_ref = 0.0)
end

function _closed!(state)
    for k in eachindex(state.ghost_cell_bc); state.ghost_cell_bc[k] = Closed; end
    return state
end

# A short point-source run on a gently sloping bed; returns the final state.
function _spread(; nsteps = 150, dt = 0.5, kw...)
    mesh = _cart_mesh(nr = 15, nc = 15, dx = 10.0)
    n = length(mesh.cells)
    mesh.static_vars["elevation"] = [0.002 * (i % 15) + 0.01 * sin(0.7 * i) for i in 1:n]
    state = _closed!(initialise_flow_model(mesh, StandardFlow(); kw...))
    src = div(n, 2) + 1
    for _ in 1:nsteps
        state.volume[src] += 1.5 * dt
        for i in 1:n
            state.water_depth[i] = state.volume[i] / state.cell_area[i]
        end
        step_standard!(state, dt)
    end
    return state, 1.5 * dt * nsteps
end

# ---------------------------------------------------------------------------
@testset "B1 — blend reproduces cell (θ = 0) and edge (θ = 1)" begin
    s_cell,  _ = _spread(; momentum_model = :cell)
    s_b0,    _ = _spread(; momentum_model = :blend, q_centre_theta = 0.0)
    @test s_b0.volume == s_cell.volume
    @test s_b0.edges.flux == s_cell.edges.flux
    s_edge1, _ = _spread(; momentum_model = :edge,  q_centre_theta = 1.0)
    s_b1,    _ = _spread(; momentum_model = :blend, q_centre_theta = 1.0)
    @test s_b1.volume == s_edge1.volume
    @test s_b1.edges.flux == s_edge1.edges.flux
    # The default θ = 0.9 blend differs from Q-centring at θ = 0.9
    s_edge9, _ = _spread(; momentum_model = :edge)
    s_b9,    _ = _spread(; momentum_model = :blend)
    @test s_b9.volume != s_edge9.volume
end

@testset "B2 — frozen isotropy, blend θ = 0.9 + vector + diamond" begin
    mesh = _cart_mesh(nr = 15, nc = 15, dx = 100.0)
    r = frozen_isotropy(mesh; angles_deg = 0:5:45, friction_model = :vector,
                        momentum_model = :blend, gradient_correction = true,
                        face_flux_method = :diamond)
    @test r.max_abs_mean_defl < 0.1
    @test r.cond_ratio < 1.01
end

# ---------------------------------------------------------------------------
@testset "P1 — _positivity_scale! unit behaviour" begin
    # cells 1..3; edges (1,2), (1,3), (2,3). Cell 1 holds 1.0 m³ and tries to
    # send 0.8 to each of 2 and 3 (edge_vol < 0: cell_i donates).
    ci = [1, 1, 2]; cj = [2, 3, 3]
    vol = [1.0, 5.0, 0.0]
    ev  = [-0.8, -0.8, 0.3]          # edge 3: cell 3 (empty) "donates" 0.3 to cell 2
    q   = [-4.0, -4.0, 1.5]
    _positivity_scale!(ev, q, vol, ci, cj)
    @test ev[1] ≈ -0.5 && ev[2] ≈ -0.5          # scaled by 1.0/1.6
    @test q[1]  ≈ -2.5 && q[2]  ≈ -2.5          # momentum written back
    @test ev[3] == 0.0 && q[3] == 0.0           # empty donor: nothing leaves
    # A donor with enough water is untouched
    ev2 = [-0.2, -0.3, 0.0]; q2 = [-1.0, -1.0, 0.0]
    _positivity_scale!(ev2, q2, vol, ci, cj)
    @test ev2 == [-0.2, -0.3, 0.0] && q2 == [-1.0, -1.0, 0.0]
end

@testset "P2 — positivity limiter: mass balance and positivity ($label)" for (label, kw) in [
        ("legacy friction/momentum", (; limiter = :positivity)),
        ("vector + blend + diamond", (; limiter = :positivity, friction_model = :vector,
                                        momentum_model = :blend, gradient_correction = true,
                                        face_flux_method = :diamond)),
        ("hybrid limiter",           (; limiter = :hybrid)),
        ("hybrid + diamond-wetdry + vector + local CFL",
                                     (; limiter = :hybrid, friction_model = :vector,
                                        gradient_correction = true,
                                        face_flux_method = :diamond_wetdry, cfl_mode = :local)),
    ]
    state, added = _spread(; kw...)
    @test minimum(state.volume) >= 0.0
    @test isapprox(sum(state.volume), added; rtol = 1e-12)
    @test count(>(1e-4), state.water_depth) > 5             # it actually spread
end

# ---------------------------------------------------------------------------
@testset "C1 — local CFL length and timestep" begin
    dx = 20.0
    state = _closed!(initialise_flow_model(_cart_mesh(dx = dx), StandardFlow(); cfl_mode = :local))
    @test all(isapprox.(state.cell_length, dx; rtol = 1e-6))
    fill!(state.water_depth, 0.0)
    state.water_depth[5] = 0.8; state.water_depth[40] = 0.2
    @test _cfl_dt(state, StandardFlow()) ≈ 0.7 * dx / sqrt(_G * 0.8) rtol = 1e-6
    # Default (legacy) mode is unchanged
    legacy = _closed!(initialise_flow_model(_cart_mesh(dx = dx), StandardFlow()))
    @test legacy.cfl_mode === :legacy && legacy.limiter === :legacy
    fill!(legacy.water_depth, 0.0); legacy.water_depth[5] = 0.8
    @test _cfl_dt(legacy, StandardFlow()) ≈ 0.7 * sqrt(minimum(legacy.cell_area)) / sqrt(_G * 0.8)
    # Option validation
    @test_throws ArgumentError initialise_flow_model(_cart_mesh(), StandardFlow(); limiter = :none)
    @test_throws ArgumentError initialise_flow_model(_cart_mesh(), StandardFlow(); cfl_mode = :global)
    @test_throws ArgumentError initialise_flow_model(_cart_mesh(), StandardFlow(); momentum_model = :vertex)
end

@testset "D1 — diamond-wetdry: stencil check and all-wet equivalence" begin
    mesh = _cart_mesh(nr = 9, nc = 9, dx = 20.0)
    n = length(mesh.cells)
    mesh.static_vars["elevation"] = [0.01 * (i % 9) for i in 1:n]
    st = initialise_flow_model(mesh, StandardFlow(); gradient_correction = true,
                               face_flux_method = :diamond_wetdry)
    tbl = st.diamond_table
    depth = fill(0.5, n)
    valid = findall(tbl.valid)
    @test !isempty(valid)
    @test all(e -> diamond_stencil_wet(tbl, e, depth, HFLOW_THRESHOLD), valid)
    dry = div(n, 2) + 1; depth[dry] = 0.0
    touching = [e for e in valid
                if any(k -> dry in tbl.vertex_table[k].cell_indices, (tbl.va_key[e], tbl.vb_key[e]))]
    @test !isempty(touching)
    @test all(e -> !diamond_stencil_wet(tbl, e, depth, HFLOW_THRESHOLD), touching)
    # Fully wet: diamond-wetdry is bit-identical to diamond
    r_d  = frozen_isotropy(_cart_mesh(nr = 15, nc = 15, dx = 100.0); angles_deg = [10.0, 30.0],
                           friction_model = :vector, gradient_correction = true,
                           face_flux_method = :diamond)
    r_wd = frozen_isotropy(_cart_mesh(nr = 15, nc = 15, dx = 100.0); angles_deg = [10.0, 30.0],
                           friction_model = :vector, gradient_correction = true,
                           face_flux_method = :diamond_wetdry)
    @test r_wd.mean_defl == r_d.mean_defl && r_wd.cond == r_d.cond
    @test_throws ArgumentError initialise_flow_model(_cart_mesh(), StandardFlow();
                                                     face_flux_method = :exact)
end

@testset "G1 — Cartesian ghost edges: one per missing side" begin
    # Regression for the phantom-ghost-edge bug (shared north–south sides were
    # matched by a longitude-only key and could be missed): a 9×9 grid has
    # 4·9 = 36 missing sides, two on each corner cell and one on each other
    # edge cell, each 20 m wide.
    state = initialise_flow_model(_cart_mesh(nr = 9, nc = 9, dx = 20.0), StandardFlow())
    @test length(state.ghost_edges) == 36
    @test all(ge -> isapprox(ge.width, 20.0; rtol = 1e-6), state.ghost_edges)
    @test count(state.boundary_mask) == 32
end

@testset "L1 — lake at rest, WP2 options ($label)" for (label, kw) in [
        ("blend + positivity + diamond", (; momentum_model = :blend, limiter = :positivity,
                                             face_flux_method = :diamond)),
        ("hybrid + diamond-wetdry",      (; limiter = :hybrid, face_flux_method = :diamond_wetdry)),
    ]
    mesh = _cart_mesh()
    n    = length(mesh.cells)
    mesh.static_vars["elevation"] = [0.6 + 0.3 * sin(0.9 * i) + 0.2 * isodd(i) for i in 1:n]
    state = _closed!(initialise_flow_model(mesh, StandardFlow();
                     friction_model = :vector, cfl_mode = :local, gradient_correction = true,
                     kw...))
    for i in 1:n
        state.water_depth[i] = 2.0 - state.elevation[i]
        state.volume[i]      = state.water_depth[i] * state.cell_area[i]
    end
    for _ in 1:200; step_standard!(state, 1.0); end
    wse = [state.elevation[i] + state.volume[i] / state.cell_area[i] for i in 1:n]
    @test maximum(abs, wse .- 2.0) < 1e-10
    @test maximum(abs, state.edges.flux) < 1e-10
end

println("\nWP2 momentum/limiter tests complete.")
