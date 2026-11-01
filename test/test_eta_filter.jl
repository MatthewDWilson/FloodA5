# test/test_eta_filter.jl
# -----------------------
# D4 (2026-10): the conservative η high-pass filter (`--eta-filter C_s`,
# surfacewater/eta_filter.jl).
#
#   E1  off by default (initialise_flow_model, the 38-argument constructor);
#       rejected for SGS and for C_s < 0
#   E2  a lake at rest over an uneven bed and a planar water surface are
#       untouched (residual zero; changes at round-off only)
#   E3  checkerboard on a square grid: r = ±1.6δ in the interior, the face
#       exchange assembled exactly from r, and the decay δ' = (1 − 12.8·C_s·c·Δt/Δx)·δ
#       away from the boundary; volume conserved
#   E4  Δt-consistency: the filter alone damps a field by the same amount
#       over the same time, whether in N steps of Δt or N/2 steps of 2Δt
#   E5  wet/dry: no water is pushed into a dry bank cell, a thin cell cannot go
#       negative, volume is conserved
#   E6  in the full step: a lake at rest stays exactly at rest, and the
#       filtered step conserves volume (closed domain)
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_eta_filter.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))

const _DX = 25.0
function _cart(nrc = 21)
    d = rad2deg(_DX / A5Grid._EARTH_R)
    return generate_cartesian_mesh_bounds(0.0, 0.0, (nrc - 0.5) * d, (nrc - 0.5) * d, _DX; lat_ref = 0.0)
end
const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _A5 = _quiet(() -> A5Grid.mesh_for_aoi(_AOI, 16))

function _state(mesh, z; cs = 0.05, kw...)
    m = deepcopy(mesh)
    m.static_vars["elevation"] = z
    s = _quiet(() -> initialise_flow_model(m, StandardFlow();
                                          merge(scheme_options(:v02, :standard),
                                                (; limiter = :positivity, q_centre_theta = 1.0,
                                                   eta_filter = cs), NamedTuple(kw))...))
    fill!(s.ghost_cell_bc, Closed)
    return s
end
_setwse!(s, η) = (s.volume .= max.(η .- s.elevation, 0.0) .* s.cell_area;
                  s.water_depth .= s.volume ./ s.cell_area; s)
_filter!(s, dt) = _eta_filter!(s, _workspace!(s), dt)
# local metric coordinates (m) of the cell centres
_xy(s) = (deg2rad.(s.cell_lons .- s.cell_lons[1]) .* A5Grid._EARTH_R .* cosd.(s.cell_lats),
          deg2rad.(s.cell_lats .- s.cell_lats[1]) .* A5Grid._EARTH_R)

@testset "E1 — off by default, argument checks" begin
    s = _quiet(() -> initialise_flow_model(deepcopy(_A5), StandardFlow()))
    @test s.eta_filter_cs == 0.0
    f = fieldnames(FlowState)
    s38 = FlowState((getfield(s, k) for k in f[1:38])...)        # pre-D4 positional form
    @test s38.eta_filter_cs == 0.0
    @test_throws ArgumentError _quiet(() -> initialise_flow_model(deepcopy(_A5), StandardFlow(); eta_filter = -0.1))
    @test_throws ArgumentError _quiet(() -> initialise_flow_model(deepcopy(_A5), SGSFlow(); eta_filter = 0.05))
end

@testset "E2 — rest and planar surfaces untouched ($lab)" for (lab, mesh) in [("A5", _A5), ("Cartesian", _cart())]
    n = length(mesh.cells)
    z = [0.3 * sin(0.7 * i) for i in 1:n]                         # uneven bed, all below the water
    s = _state(mesh, z)
    _setwse!(s, fill(2.0, n)); V0 = copy(s.volume)
    _filter!(s, 5.0)
    @test maximum(abs.(s.volume .- V0) ./ s.cell_area) < 1e-12
    # tilted plane, linear in (lon, lat) — exactly planar in each cell's local
    # equirectangular frame, which the fit uses (~1e-3 / −4e-4 m per m)
    R = A5Grid._EARTH_R; c0 = cosd(s.cell_lats[1])
    _setwse!(s, 2.0 .+ 1e-3 .* deg2rad.(s.cell_lons .- s.cell_lons[1]) .* R .* c0 .-
                4e-4 .* deg2rad.(s.cell_lats .- s.cell_lats[1]) .* R); V0 = copy(s.volume)
    _filter!(s, 5.0)
    @test maximum(abs.(s.volume .- V0) ./ s.cell_area) < 1e-11
end

@testset "E3 — checkerboard decay on a square grid" begin
    mesh = _cart(); n = length(mesh.cells)
    cs = 0.05; dt = 2.0; H = 3.0; δ = 0.01
    s = _state(mesh, zeros(n); cs = cs)
    x, y = _xy(s)
    sgn = [iseven(round(Int, x[i] / _DX) + round(Int, y[i] / _DX)) ? 1.0 : -1.0 for i in 1:n]
    _setwse!(s, H .+ δ .* sgn); V0 = sum(s.volume); η0 = H .+ δ .* sgn
    _filter!(s, dt)
    @test isapprox(sum(s.volume), V0; rtol = 1e-14)
    ws = s.ws
    nnb(i) = count(!iszero, view(s.adj_matrix, :, i))
    inner = [i for i in 1:n if nnb(i) == 4]
    @test maximum(abs.(ws.f_r[inner] .- 1.6δ .* sgn[inner])) < 1e-12   # r = ±1.6δ
    corner = [i for i in 1:n if nnb(i) == 2]
    @test length(corner) == 4 && maximum(abs.(ws.f_r[corner])) < 1e-8 * δ  # 3 points: exact fit
    E = s.edges; c = sqrt(_G * (H - δ)); r = ws.f_r
    expect = copy(η0)
    for e in 1:E.n_edges                                           # F_f = C_s c w (r_j − r_i) dt
        i, j = E.cell_i[e], E.cell_j[e]
        F = cs * c * E.width[e] * (r[j] - r[i]) * dt
        expect[i] += F / s.cell_area[i]; expect[j] -= F / s.cell_area[j]
    end
    η1 = s.elevation .+ s.volume ./ s.cell_area
    @test maximum(abs.(η1 .- expect)) < 1e-12
    deep = [i for i in inner if all(j -> j == 0 || nnb(j) == 4, view(s.adj_matrix, :, i))]
    β = 12.8 * cs * c * dt / _DX                                   # δ' = (1 − β)δ
    @test maximum(abs.((η1[deep] .- H) .- (1 - β) .* δ .* sgn[deep])) < 1e-6 * δ
end

@testset "E4 — damping per unit time independent of Δt ($lab)" for (lab, mesh) in [("A5", _A5), ("Cartesian", _cart())]
    n = length(mesh.cells)
    pert = [1e-3 * sin(2.1 * i) for i in 1:n]                      # rough field
    function damp(dt, nsteps)
        s = _state(mesh, zeros(n); cs = 0.02)
        _setwse!(s, 3.0 .+ pert)
        for _ in 1:nsteps; _filter!(s, dt); end
        η = s.volume ./ s.cell_area
        return sqrt(sum(abs2, η .- sum(η .* s.cell_area) / sum(s.cell_area)))
    end
    a = damp(0.5, 80); b = damp(1.0, 40); c = damp(2.0, 20)
    r0 = sqrt(sum(abs2, pert .- sum(pert) / n))
    @test a < 0.9 * r0                                             # it does damp
    @test isapprox(a, b; rtol = 0.02) && isapprox(a, c; rtol = 0.05)
end

@testset "E5 — wet/dry and positivity" begin
    mesh = _cart(); n = length(mesh.cells)
    s = _state(mesh, zeros(n); cs = 0.05)
    x, _ = _xy(s)
    bank = x .> 10 * _DX                                           # dry bank above the water
    s.elevation[bank] .= 5.0
    E = s.edges
    for e in 1:E.n_edges; E.sill[e] = max(s.elevation[E.cell_i[e]], s.elevation[E.cell_j[e]]); end
    η = [bank[i] ? 5.0 : 2.0 + 0.2 * (-1)^i for i in 1:n]
    _setwse!(s, η)
    thin = findfirst(i -> !bank[i] && η[i] < 2.0, 1:n)            # thin-film cell with neighbours above it
    s.volume[thin] = 1e-6 * s.cell_area[thin]; s.water_depth[thin] = 1e-6
    V0 = sum(s.volume)
    for _ in 1:20; _filter!(s, 10.0); end
    @test all(s.volume[bank] .== 0.0)
    @test all(s.volume .>= 0.0)
    @test isapprox(sum(s.volume), V0; rtol = 1e-13)
end

@testset "E6 — full step: rest preserved, volume conserved ($lab)" for (lab, mesh) in [("A5", _A5), ("Cartesian", _cart())]
    n = length(mesh.cells)
    z = [0.5 * sin(0.37 * i) for i in 1:n]
    s = _state(mesh, z; cs = 0.05)
    _setwse!(s, fill(2.0, n)); V0 = copy(s.volume)
    for _ in 1:50; step_standard!(s, 3.0); end
    @test maximum(abs.(s.volume .- V0) ./ s.cell_area) < 1e-12
    _setwse!(s, 2.0 .+ [1e-2 * (-1)^i for i in 1:n]); M0 = sum(s.volume)
    for _ in 1:200
        s.water_depth .= s.volume ./ s.cell_area
        step_standard!(s, 3.0)
    end
    @test isapprox(sum(s.volume), M0; rtol = 1e-12)
    @test all(s.volume .>= 0.0)
end

println("\nη-filter tests complete.")
