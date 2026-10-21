# test/test_wp1_vector_friction.jl
# --------------------------------
# Tests for WP1 vector friction (--friction-model vector, 2026-10 — see
# FloodA5_Review_and_ForwardPlan_2026-10.md §3, §9 WP1).
#
#   V1  Kernels: defaults are bit-identical to component friction; q_t is
#       ignored unless vector_friction = true
#   V2  Kernels: vector friction uses √(q_n² + q_t²) in the friction term
#   V3  Kernels: the vector Froude cap scales q_n by q_max/|q|
#   V4  R-A kernels (SGS): vector friction uses √(Q_n² + Q_t²)
#   V5  Frozen-coefficient friction-limit isotropy on a Cartesian grid,
#       through step_standard!: component friction deflects the flow by up to
#       ≈12.8° (theory) and has a 2^¼ conductance anisotropy; vector friction
#       removes both, for edge and cell momentum and for the diamond gradient
#   V6  Lake at rest under vector friction
#   V7  initialise_flow_model rejects an unknown friction model
#
# Uses an in-process Cartesian mesh — no data files.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_wp1_vector_friction.jl

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

# ---------------------------------------------------------------------------
@testset "V1 — defaults are component friction, bit-identical" begin
    args = (0.05, 1.20, 1.00, 0.70, 10.0, 12.0, 0.9, 0.03, 2.0, 0.5)
    @test _bates_flux_limited(args...) ===
          _bates_flux_limited(args...; q_t = 0.3)                  # q_t ignored
    @test _bates_flux_limited(args...) ===
          _bates_flux_limited(args...; q_t = 0.0, vector_friction = true)
    @test _bates_flux(args[1:9]...) === _bates_flux(args[1:9]...; q_t = 0.3)
    cargs = (0.05, 0.4, 0.15, 10.0, 12.0, 0.03, 2.0, 0.5)
    @test _bates_flux_limited_corrected(cargs...) ===
          _bates_flux_limited_corrected(cargs...; q_t = 0.3)
    @test _bates_flux_corrected(cargs[1:7]...) ===
          _bates_flux_corrected(cargs[1:7]...; q_t = 0.3)
    @test _friction_q(-0.2, 0.7, false) == 0.2
    @test _friction_q(-0.3, 0.4, true)  ≈ 0.5
end

@testset "V2 — vector friction denominator uses |q|" begin
    q_prev, q_t, h, dWSE_n, L, n, dt = 0.04, 0.03, 0.4, 0.002, 10.0, 0.03, 1.0
    num  = q_prev - _G * h * dt * dWSE_n / L
    den(qm) = 1.0 + _G * h * dt * n^2 * qm / h^(10/3)
    q_c  = _bates_flux_corrected(q_prev, h, dWSE_n, 1.0, L, n, dt)
    q_v  = _bates_flux_corrected(q_prev, h, dWSE_n, 1.0, L, n, dt;
                                 q_t = q_t, vector_friction = true)
    @test q_c ≈ num / den(abs(q_prev))
    @test q_v ≈ num / den(hypot(q_prev, q_t))
    @test abs(q_v) < abs(q_c)                  # more friction with a tangential flow
    # Legacy (cos θ) kernel, unlimited and limited, agree with the corrected one
    # at orthogonality (dWSE_n = wse_i − wse_j, c = 1).
    wi, wj, zs = 1.402, 1.400, 1.002
    Q1 = _bates_flux(q_prev, wi, wj, zs, 1.0, L, 1.0, n, dt; q_t = q_t, vector_friction = true)
    Q2 = _bates_flux_corrected(q_prev, wi - zs, wi - wj, 1.0, L, n, dt; q_t = q_t, vector_friction = true)
    @test Q1 ≈ Q2
end

@testset "V3 — vector Froude cap acts on |q|" begin
    q_max = 1.0
    @test _froude_cap(1.5, 5.0, q_max, false) == 1.0          # component: clamp q_n only
    @test _froude_cap(0.6, 0.8, q_max, true)  == 0.6          # |q| = 1.0: not capped
    q = _froude_cap(0.9, 1.2, q_max, true)                     # |q| = 1.5 → scale by 2/3
    @test q ≈ 0.6
    @test _froude_cap(-0.9, 1.2, q_max, true) ≈ -0.6           # sign preserved
    # Through the limited kernel: a deep, fast edge with a large tangential flow
    h = 0.2; L = 10.0
    Q, qs = _bates_flux_limited_corrected(0.0, h, 0.5, 1.0, L, 0.01, 1.0, 0.0;
                                          q_t = 2.0, vector_friction = true)
    qmax = h * sqrt(_G * h) * FROUDE_LIMIT
    @test hypot(qs, 2.0) ≥ qmax - 1e-12                       # scaled onto the cap...
    @test abs(qs) < qmax                                      # ...so |q_n| is below the component cap
end

@testset "V4 — R-A kernels use |Q| for vector friction" begin
    Qp, Qt, A, R, L, n, dt = 2.0, 1.5, 3.0, 0.5, 20.0, 0.03, 1.0
    wi, wj, zs = 2.01, 2.00, 1.0
    Qc = _manning_flux_ra(Qp, wi, wj, zs, A, R, L, 1.0, n, dt)
    @test Qc === _manning_flux_ra(Qp, wi, wj, zs, A, R, L, 1.0, n, dt; Q_t = Qt)
    Qv = _manning_flux_ra(Qp, wi, wj, zs, A, R, L, 1.0, n, dt; Q_t = Qt, vector_friction = true)
    num = Qp - _G * A * dt * (wi - wj) / L
    @test Qv ≈ num / (1.0 + _G * dt * n^2 * hypot(Qp, Qt) / (R^(4/3) * A))
    Qvc = _manning_flux_ra_corrected(Qp, wi - zs, A, R, wi - wj, L, n, dt;
                                     Q_t = Qt, vector_friction = true)
    @test Qvc ≈ Qv
end

# ---------------------------------------------------------------------------
# V5 — frozen-coefficient isotropy. 15×15 cells of 100 m; bearings 0°–45° in
# 5° steps (the Cartesian pattern repeats every 45°). Theory for component
# friction: ψ = atan(√tan φ), max deflection 12.8° at φ ≈ 10.7°, conductance
# ratio 2^¼ ≈ 1.189 between diagonal and cardinal bearings.
const _V5_MESH = _cart_mesh(nr = 15, nc = 15, dx = 100.0)
_v5(; kw...) = frozen_isotropy(_V5_MESH; angles_deg = 0:5:45, kw...)

@testset "V5 — frozen-coefficient isotropy ($label)" for (label, kw, biased) in [
        ("component, edge",            (; friction_model = :component), true),
        ("vector, edge",               (; friction_model = :vector), false),
        ("component, cell",            (; friction_model = :component, momentum_model = :cell), true),
        ("vector, cell",               (; friction_model = :vector, momentum_model = :cell), false),
        ("vector, diamond + cell",     (; friction_model = :vector, momentum_model = :cell,
                                          gradient_correction = true, face_flux_method = :diamond), false),
    ]
    r = _v5(; kw...)
    @test r.n_interior > 50
    if biased
        theory = maximum(abs(atand(sqrt(tand(φ))) - φ) for φ in 0:5:45)   # 12.6° on this grid
        @test r.max_abs_mean_defl > 0.8 * theory
        @test r.cond_ratio > 1.15
    else
        @test r.max_abs_mean_defl < 0.1
        @test r.cond_ratio < 1.01
    end
end

# ---------------------------------------------------------------------------
@testset "V6 — lake at rest, vector friction ($label)" for (label, kw) in [
        ("edge momentum",           (; friction_model = :vector)),
        ("diamond + cell momentum", (; friction_model = :vector, gradient_correction = true,
                                       face_flux_method = :diamond, momentum_model = :cell)),
    ]
    mesh = _cart_mesh()
    n    = length(mesh.cells)
    mesh.static_vars["elevation"] = [0.6 + 0.3 * sin(0.9 * i) + 0.2 * isodd(i) for i in 1:n]
    state = initialise_flow_model(mesh, StandardFlow(); kw...)
    for k in eachindex(state.ghost_cell_bc); state.ghost_cell_bc[k] = Closed; end
    for i in 1:n
        state.water_depth[i] = 2.0 - state.elevation[i]
        state.volume[i]      = state.water_depth[i] * state.cell_area[i]
    end
    for _ in 1:200; step_standard!(state, 1.0); end
    wse = [state.elevation[i] + state.volume[i] / state.cell_area[i] for i in 1:n]
    @test maximum(abs, wse .- 2.0) < 1e-10
    @test maximum(abs, state.edges.flux) < 1e-10
    @test maximum(abs, state.qvec_u) < 1e-10 && maximum(abs, state.qvec_v) < 1e-10
end

@testset "V7 — unknown friction model is rejected" begin
    @test_throws ArgumentError initialise_flow_model(_cart_mesh(), StandardFlow();
                                                     friction_model = :decoupled)
    @test initialise_flow_model(_cart_mesh(), StandardFlow()).friction_model === :component
end

println("\nWP1 vector-friction tests complete.")
