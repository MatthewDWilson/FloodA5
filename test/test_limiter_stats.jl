# test/test_limiter_stats.jl
# --------------------------
# D1b (2026-10): `--limiter-stats` must not change results, and must record
# sensible flags.
#
#   L1  steps with stats on are bit-identical to steps with stats off
#       (standard: recommended, legacy limiter, positivity; SGS) on a partly wet
#       A5 domain draining down a slope
#   L2  the flags: wet edges are counted, a Froude-capped or volume-limited
#       edge is detected, the report runs
#   L3  `--boundary-limiter positivity` (D1b experiment): off by default; when
#       on it changes boundary outflow near CFL, keeps volumes non-negative and
#       conserves mass
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_limiter_stats.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _A5 = A5Grid.mesh_for_aoi(_AOI, 16)
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

function _state(method; kw...)
    mesh = deepcopy(_A5)
    x = [c.center_lon for c in mesh.cells]
    z = 5e-3 .* 111_000 .* (maximum(x) .- x)                 # steep: fast, Froude-capped flow
    mesh.static_vars["elevation"] = z
    s = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); kw...))
    fill!(s.ghost_cell_bc, ZeroGradient)
    top = sortperm(x)[1:40]
    fill!(s.volume, 0.0); s.volume[top] .= 1.5 .* s.cell_area[top]
    s.water_depth .= s.volume ./ s.cell_area
    if method isa SGSFlow
        nb = 8; nsl = s.max_neighbours
        s.sgs_tables = [begin
            zi = z[i]; A = s.cell_area[i]; eb = collect(range(zi, zi + 2.0; length = nb))
            SGSTable(eb, A .* (eb .- zi), fill(A, nb), A, zi, zi + 2.0,
                     repeat(10.0 .* (eb .- zi) .+ 0.5, 1, nsl), repeat(10.0 .+ 2.0 .* (eb .- zi), 1, nsl))
        end for i in eachindex(s.cell_ids)]
    end
    return s
end

function _run!(s, method, stats::Bool; steps = 150, dt = 2.0)
    LIMITER_STATS[] = stats
    stats && _lstats_reset!(s.edges.n_edges)
    for k in 1:steps
        method isa SGSFlow || (s.water_depth .= s.volume ./ s.cell_area)
        method isa SGSFlow ? step_sgs!(s, dt) : step_standard!(s, dt)
        stats && _lstats_step!(dt * k)
    end
    LIMITER_STATS[] = false
    return s
end

@testset "L1 — stats on == stats off ($lab)" for (lab, method, kw) in [
        ("std recommended",    StandardFlow(), scheme_options(:recommended, :standard)),
        ("std legacy limiter", StandardFlow(), merge(scheme_options(:recommended, :standard), (; limiter = :legacy))),
        ("std positivity",     StandardFlow(), merge(scheme_options(:recommended, :standard), (; limiter = :positivity))),
        ("sgs recommended",    SGSFlow(),      scheme_options(:recommended, :sgs)),
    ]
    a = _run!(_state(method; kw...), method, false)
    b = _run!(_state(method; kw...), method, true)
    @test a.volume == b.volume && a.edges.flux == b.edges.flux && a.edges.flux_Q == b.edges.flux_Q
    @test a.vol_removed == b.vol_removed
end

@testset "L2 — flags recorded" begin
    s = _run!(_state(StandardFlow(); merge(scheme_options(:recommended, :standard),
                                            (; limiter = :legacy))...), StandardFlow(), true)
    st = _LSTATS[]
    @test st.steps == 150
    @test sum(st.wet) > 0
    @test sum(st.count[1, :]) + sum(st.count[2, :]) > 0        # Froude cap or volume limiter bound
    @test all(st.count .<= reshape(st.wet, 1, 4))              # no flag on more edges than are wet
    @test (_quiet(_lstats_report); true)
end

@testset "L3 — --boundary-limiter positivity ($lab)" for (lab, method, kw) in [
        ("std positivity", StandardFlow(), merge(scheme_options(:recommended, :standard), (; limiter = :positivity))),
        ("sgs recommended", SGSFlow(),     scheme_options(:recommended, :sgs)),
    ]
    @test BOUNDARY_POSITIVITY[] == false                       # off by default
    s0 = _state(method; kw...); M0 = sum(s0.volume)
    a = _run!(deepcopy(s0), method, false; steps = 100, dt = 10.0)   # near CFL: the ghost limits bind
    BOUNDARY_POSITIVITY[] = true
    b = try _run!(deepcopy(s0), method, false; steps = 100, dt = 10.0) finally BOUNDARY_POSITIVITY[] = false end
    @test a.vol_removed > 0 && b.vol_removed != a.vol_removed  # the option acts on boundary outflow
    @test all(b.volume .>= 0.0)
    @test isapprox(sum(b.volume) + b.vol_removed, M0; rtol = 1e-12)
end

println("\nLimiter-stats tests complete.")
