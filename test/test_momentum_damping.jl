# test/test_momentum_damping.jl
# -----------------------------
# D6′ (2026-10): divergence-free damping of circulating face discharges
# (`--momentum-damping C_m`, surfacewater/momentum_damping.jl).
#
#   M1  off by default (initialise_flow_model, the 39-argument constructor);
#       argument checks (SGS, no diamond gradient, component friction)
#   M2  the change is exactly divergence-free: no cell's net discharge
#       changes, for a random discharge field (A5 and Cartesian)
#   M3  a uniform velocity field is untouched, also on the skewed A5 mesh
#   M4  a purely circulating field decays, monotonically, at a rate per unit
#       time independent of Δt; :flow and :wave both damp
#   M5  in the full step: a lake at rest stays at rest; volume is conserved
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_momentum_damping.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging, Random

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))

const _DX = 25.0
function _cart(nrc = 21)
    d = rad2deg(_DX / A5Grid._EARTH_R)
    return generate_cartesian_mesh_bounds(0.0, 0.0, (nrc - 0.5) * d, (nrc - 0.5) * d, _DX; lat_ref = 0.0)
end
const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _A5 = _quiet(() -> A5Grid.mesh_for_aoi(_AOI, 16))
const _MESHES = [("A5", _A5), ("Cartesian", _cart())]

function _state(mesh; cm = 0.1, scale = :wave, z = nothing, kw...)
    m = deepcopy(mesh)
    m.static_vars["elevation"] = z === nothing ? zeros(length(m.cells)) : z
    s = _quiet(() -> initialise_flow_model(m, StandardFlow();
                                          merge(scheme_options(:recommended, :standard),
                                                (; limiter = :positivity, q_centre_theta = 1.0,
                                                   momentum_damping = cm, momentum_damping_scale = scale),
                                                NamedTuple(kw))...))
    fill!(s.ghost_cell_bc, Closed)
    return s
end
_fill!(s, H) = (s.volume .= H .* s.cell_area; s.water_depth .= H; s)
_damp!(s, dt) = (_reconstruct_qvec!(s, false); _momentum_damping!(s, _workspace!(s), dt))
function _netQ(s)                                   # per-cell net outflow of the face discharges
    E = s.edges; d = zeros(length(s.cell_ids))
    for e in 1:E.n_edges
        Q = E.flux[e] * E.width[e]; d[E.cell_i[e]] += Q; d[E.cell_j[e]] -= Q
    end
    return d
end
# a purely circulating field from a random streamfunction at the interior vertices
function _vortex!(s, rng; amp = 100.0)          # amp: ~0.5 m/s face velocities
    ws = _mom_damp_ws!(s, _workspace!(s)); E = s.edges
    ψ = [ws.m_interior[v] ? amp * randn(rng) : 0.0 for v in eachindex(ws.m_interior)]
    for e in 1:E.n_edges
        E.flux[e] = ws.m_tail[e] == 0 ? 0.0 : (ψ[ws.m_head[e]] - ψ[ws.m_tail[e]]) / E.width[e]
    end
    return s
end
_enstrophy(s) = sum(abs2, s.edges.flux)

@testset "M1 — off by default, argument checks" begin
    s = _quiet(() -> initialise_flow_model(deepcopy(_A5), StandardFlow()))
    @test s.mom_damping_cs == 0.0 && s.mom_damping_scale === :wave
    f = fieldnames(FlowState)
    s39 = FlowState((getfield(s, k) for k in f[1:39])...)
    @test s39.mom_damping_cs == 0.0
    init(m; kw...) = _quiet(() -> initialise_flow_model(deepcopy(_A5), m; kw...))
    @test_throws ArgumentError init(StandardFlow(); momentum_damping = -1.0)
    @test_throws ArgumentError init(SGSFlow(); momentum_damping = 0.1)
    @test_throws ArgumentError init(StandardFlow(); momentum_damping = 0.1)           # no diamond gradient
    @test_throws ArgumentError init(StandardFlow(); merge(scheme_options(:recommended, :standard),
                                    (; friction_model = :component, momentum_damping = 0.1))...)
    @test_throws ArgumentError init(StandardFlow(); merge(scheme_options(:recommended, :standard),
                                    (; momentum_damping = 0.1, momentum_damping_scale = :other))...)
end

@testset "M2 — divergence-free ($lab)" for (lab, mesh) in _MESHES
    s = _fill!(_state(mesh), 2.0)
    rng = MersenneTwister(1)
    s.edges.flux .= randn(rng, s.edges.n_edges)
    d0 = _netQ(s); q0 = copy(s.edges.flux)
    _damp!(s, 5.0)
    Δ = maximum(abs.(s.edges.flux .- q0) .* s.edges.width)
    @test Δ > 1e-3                                              # it acts
    @test maximum(abs.(_netQ(s) .- d0)) < 1e-11 * max(Δ, 1.0)  # but moves no water
end

@testset "M3 — uniform velocity untouched ($lab)" for (lab, mesh) in _MESHES
    s = _fill!(_state(mesh), 2.0)
    U = (0.6, -0.35); E = s.edges
    for e in 1:E.n_edges
        E.flux[e] = -2.0 * (U[1] * E.nf_x[e] + U[2] * E.nf_y[e])   # flux > 0: cell_j → cell_i
    end
    q0 = copy(E.flux)
    _damp!(s, 5.0)
    @test maximum(abs.(E.flux .- q0)) < 1e-5 * 1.4          # local-frame closure of the d_e loops
end

@testset "M4 — vortex decay, Δt-independent ($lab)" for (lab, mesh) in _MESHES
    for scale in (:wave, :flow)
        s = _vortex!(_fill!(_state(mesh; scale = scale), 2.0), MersenneTwister(2))
        e0 = _enstrophy(s); es = Float64[]
        for _ in 1:10; _damp!(s, 2.0); push!(es, _enstrophy(s)); end
        @test es[end] < 0.95 * e0 && all(diff([e0; es]) .< 0)
    end
    function decay(dt, nsteps)
        s = _vortex!(_fill!(_state(mesh; cm = 0.05), 2.0), MersenneTwister(3))
        e0 = _enstrophy(s)
        for _ in 1:nsteps; _damp!(s, dt); end
        return _enstrophy(s) / e0
    end
    a = decay(0.5, 40); b = decay(1.0, 20); c = decay(2.0, 10)
    @test a < 0.9
    @test isapprox(a, b; rtol = 0.03) && isapprox(a, c; rtol = 0.08)
end

@testset "M5 — full step: rest and mass ($lab)" for (lab, mesh) in _MESHES
    n = length(mesh.cells)
    z = [0.5 * sin(0.37 * i) for i in 1:n]
    s = _state(mesh; z = z, cm = 0.1, eta_filter = 0.05)
    s.volume .= max.(2.0 .- s.elevation, 0.0) .* s.cell_area; s.water_depth .= s.volume ./ s.cell_area
    V0 = copy(s.volume)
    for _ in 1:50; step_standard!(s, 3.0); end
    @test maximum(abs.(s.volume .- V0) ./ s.cell_area) < 1e-12
    s.volume .+= [1e-2 * (-1)^i for i in 1:n] .* s.cell_area; M0 = sum(s.volume)
    for _ in 1:200
        s.water_depth .= s.volume ./ s.cell_area
        step_standard!(s, 3.0)
    end
    @test isapprox(sum(s.volume), M0; rtol = 1e-12)
    @test all(s.volume .>= 0.0)
end

println("\nMomentum-damping tests complete.")
