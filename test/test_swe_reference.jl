# test/test_swe_reference.jl
# --------------------------
# R1 (2026-10): the full shallow-water reference solver (`--flow-model swe`,
# surfacewater/swe_reference.jl).
#
#   S1  well-balanced: a lake at rest over an uneven bed with an emerged,
#       dry island stays at rest to round-off (A5 and Cartesian, orders 1, 2)
#   S2  closed domain sloshing over the island: volume conserved to round-off,
#       volumes non-negative, nothing clamped
#   S3  Stoker dam break (wet bed, frictionless) against the exact solution:
#       the error falls with resolution, order 2 beats order 1
#   S4  Ritter dam break (dry bed): front position and profile
#   S5  open boundaries (outflow-only, as in the standard solver): storage +
#       outflow = initial volume
#   S6  uniform Manning flow on a slope (fully transmissive edges) holds
#       normal depth and velocity away from the boundaries
#   S7  local-inertial states carry no SWE data; the 41-argument constructor
#   S8  the eddy viscosity (--swe-viscosity) removes shear energy, keeps mass
#   S9  HDF5 frames (including the last, written after the loop) carry the
#       SWE cell velocity
#   S10 thin sheet flow (2 cm) at normal depth: second-order face bed and
#       stage-start friction
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_swe_reference.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))
const _R = A5Grid._EARTH_R
_deg(m) = rad2deg(m / _R)
const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _A5 = _quiet(() -> A5Grid.mesh_for_aoi(_AOI, 16))
_cart(nx, ny, dx) = generate_cartesian_mesh_bounds(0.0, 0.0, (nx - 0.5) * _deg(dx), (ny - 0.5) * _deg(dx), dx; lat_ref = 0.0)
const _MESHES = [("A5", _A5), ("Cartesian", _cart(31, 31, 25.0))]

function _swe(mesh, z; order = 2, n = 0.03, bc = Closed)
    m = deepcopy(mesh); m.static_vars["elevation"] = z
    s = _quiet(() -> initialise_flow_model(m, SWEFlow(); swe_order = order, manning_n = n))
    fill!(s.ghost_cell_bc, bc)
    return s
end
_setwse!(s, η) = (s.volume .= max.(η .- s.elevation, 0.0) .* s.cell_area;
                  s.water_depth .= s.volume ./ s.cell_area; s)
function _run!(s, T; courant = 0.7)
    t = 0.0
    while t < T
        dt = min(_cfl_dt(s, SWEFlow(); courant = courant), T - t)
        step_swe!(s, dt); t += dt
    end
    return s
end
_xn(mesh) = (x = [c.center_lon for c in mesh.cells]; y = [c.center_lat for c in mesh.cells];
             ((x .- minimum(x)) ./ (maximum(x) - minimum(x)), (y .- minimum(y)) ./ (maximum(y) - minimum(y))))
_island(mesh) = ((xm, ym) = _xn(mesh);
                 1.5 .* exp.(-((xm .- 0.5) .^ 2 .+ (ym .- 0.5) .^ 2) ./ 0.02) .+ 0.2 .* sin.(7 .* xm))

@testset "S1 — lake at rest ($lab, order $o)" for (lab, mesh) in _MESHES, o in (1, 2)
    s = _setwse!(_swe(mesh, _island(mesh); order = o), fill(1.0, length(mesh.cells)))
    V0 = copy(s.volume)
    @test count(==(0.0), V0) > 0                                   # the island is dry
    for _ in 1:200; step_swe!(s, _cfl_dt(s, SWEFlow())); end
    @test maximum(abs.(s.volume .- V0) ./ s.cell_area) < 1e-12
    @test maximum(hypot.(s.qvec_u, s.qvec_v)) < 1e-12
end

@testset "S2 — mass and positivity, closed sloshing ($lab, order $o)" for (lab, mesh) in _MESHES, o in (1, 2)
    xm, _ = _xn(mesh)
    s = _setwse!(_swe(mesh, _island(mesh); order = o), 1.0 .+ 0.4 .* (xm .- 0.5))
    M0 = sum(s.volume)
    _run!(s, 600.0)
    @test isapprox(sum(s.volume), M0; rtol = 1e-12)
    @test all(s.volume .>= 0.0) && sum(s.swe.clamp_vol) == 0.0
    @test s.vol_removed == 0.0
end

# exact dam-break solutions (frictionless, flat bed, dam at x = 0)
function _stoker(x, t, hL, hR)
    g = _G; cL = sqrt(g * hL)
    hR == 0.0 && return x < -cL * t ? hL : x <= 2cL * t ? (2cL - x / t)^2 / (9g) : 0.0
    f(h) = 2 * (sqrt(g * h) - cL) + (h - hR) * sqrt(g / 2 * (1 / h + 1 / hR))
    a, b = hR, hL
    for _ in 1:200; m = 0.5(a + b); f(m) > 0 ? (b = m) : (a = m); end
    hs = 0.5(a + b); cs = sqrt(g * hs); us = 2(cL - cs); S = hs * us / (hs - hR)
    x < -cL * t && return hL
    x < (us - cs) * t && return (2cL - x / t)^2 / (9g)
    x < S * t && return hs
    return hR
end
function _dambreak(dx, order, hL, hR; T = 20.0)
    nx = round(Int, 1000.0 / dx)
    mesh = _cart(nx, 3, dx)
    s = _swe(mesh, zeros(length(mesh.cells)); order = order, n = 0.0)
    x = deg2rad.(s.cell_lons .- s.cell_lons[1]) .* _R .- 500.0
    _setwse!(s, [xi < 0 ? hL : hR for xi in x])
    _run!(s, T; courant = 0.5)
    h = s.volume ./ s.cell_area
    ex = [_stoker(xi, T, hL, hR) for xi in x]
    return sum(abs.(h .- ex)) / length(h), h, x, s
end

@testset "S3 — Stoker dam break" begin
    e1c, = _dambreak(10.0, 1, 2.0, 0.5); e1f, = _dambreak(5.0, 1, 2.0, 0.5)
    e2c, = _dambreak(10.0, 2, 2.0, 0.5); e2f, = _dambreak(5.0, 2, 2.0, 0.5)
    @info "Stoker L1(h): order 1 $(e1c) → $(e1f), order 2 $(e2c) → $(e2f)"
    @test e1f < 0.75 * e1c && e2f < 0.75 * e2c                      # converges
    @test e2c < e1c && e2f < e1f                                    # order 2 more accurate
    @test e2f < 0.02                                                # m, mean abs error
end

@testset "S4 — Ritter dam break (dry bed)" begin
    # front position where the depth first exceeds th: exact (2c₀ − √(9g·th))·t
    front(h, x, th) = maximum(x[i] for i in eachindex(h) if h[i] > th)
    exact(th) = (2 * sqrt(_G * 2.0) - sqrt(9 * _G * th)) * 20.0
    ec, hc, xc, sc = _dambreak(5.0, 2, 2.0, 0.0)
    ef, hf, xf, sf = _dambreak(2.5, 2, 2.0, 0.0)
    @test ef < 0.6 * ec && ef < 0.005                               # converges (m, mean abs)
    @test abs(front(hf, xf, 0.1) - exact(0.1)) / exact(0.1) < 0.05  # the bulk of the wave
    # the thin tip lags (typical of finite volumes on a dry bed) but closes with resolution
    @test front(hc, xc, 1e-2) < front(hf, xf, 1e-2) <= exact(1e-2) + 2.5
    @test abs(front(hf, xf, 1e-2) - exact(1e-2)) / exact(1e-2) < 0.15
    @test all(sc.volume .>= 0.0) && all(sf.volume .>= 0.0)
end

@testset "S5 — open boundaries: storage + outflow ($lab)" for (lab, mesh) in _MESHES
    xm, _ = _xn(mesh)
    s = _swe(mesh, 0.5 .* (1 .- xm); bc = ZeroGradient)
    _setwse!(s, s.elevation .+ 0.5)                                 # uniform depth on a slope: drains
    M0 = sum(s.volume)
    _run!(s, 900.0)
    @test s.vol_removed > 0.01 * M0                                 # water left
    @test isapprox(sum(s.volume) + s.vol_removed, M0; rtol = 1e-12)
    @test all(s.volume .>= 0.0)
end

@testset "S6 — uniform Manning flow at normal depth ($lab)" for (lab, mesh) in _MESHES
    # Fully transmissive edges (swe_open_inflow) so that the flowing state is
    # not disturbed at the inflow side; cells that boundary signals (speed
    # |u| + c, with margin) cannot have reached in T must keep the normal
    # depth and velocity: gravity and friction stay in balance on the stepped
    # (cell-constant) bed. Measured (order 2, T = 50 s): A5 0.4% depth, 1.1%
    # velocity, 0.3% cross-flow; Cartesian 0.07%, 0.3%, 0 — the stepped bed's
    # discretisation error, larger on the irregular mesh.
    S = 1e-3; n = 0.03; h0 = 1.0; T = 50.0
    st = _swe(mesh, zeros(length(mesh.cells)); n = n, bc = ZeroGradient)
    x = deg2rad.(st.cell_lons .- minimum(st.cell_lons)) .* _R .* cosd.(st.cell_lats)
    y = deg2rad.(st.cell_lats .- minimum(st.cell_lats)) .* _R
    m = deepcopy(mesh); m.static_vars["elevation"] = S .* (maximum(x) .- x)
    s = _quiet(() -> initialise_flow_model(m, SWEFlow(); manning_n = n, swe_open_inflow = true))
    fill!(s.ghost_cell_bc, ZeroGradient)
    _setwse!(s, s.elevation .+ h0)
    U = h0^(2 / 3) * sqrt(S) / n                                    # Manning, towards +x
    s.qvec_u .= h0 * U; s.qvec_v .= 0.0
    _run!(s, T)
    bnd = findall(s.boundary_mask)
    reach = 1.3 * (U + sqrt(_G * h0)) * T
    inner = [i for i in eachindex(x) if minimum(hypot(x[i] - x[b], y[i] - y[b]) for b in bnd) > reach]
    @test length(inner) >= 10
    h = s.volume[inner] ./ s.cell_area[inner]; u = s.qvec_u[inner] ./ h
    @info "$lab uniform flow ($(length(inner)) cells): depth $(extrema(h)), u $(extrema(u)) (normal $h0, $(round(U, digits = 4)))"
    @test maximum(abs.(h .- h0)) < 0.01 * h0
    @test maximum(abs.(u .- U)) < 0.02 * U
    @test maximum(abs.(s.qvec_v[inner])) < 0.005 * h0 * U
end

@testset "S10 — thin sheet flow at normal depth ($lab)" for (lab, mesh) in _MESHES
    # S6 at h0 = 2 cm, where the bed step across a face (2.5 cm on the 25 m
    # Cartesian grid, ~9 cm on A5 res-16) exceeds the depth. Needs the
    # second-order face bed (z_L = η_L − h_L) and the semi-implicit friction
    # with |q| from the stage start; before both, the sheet ran 34–80% slow.
    # Measured (order 2, T = 300 s): Cartesian 0.04% depth, 0.1% velocity;
    # A5 up to 7.9% / 7.5% (the irregular stepped bed's WLSQ face bed).
    S = 1e-3; n = 0.03; h0 = 0.02; T = 300.0
    st = _swe(mesh, zeros(length(mesh.cells)); n = n, bc = ZeroGradient)
    x = deg2rad.(st.cell_lons .- minimum(st.cell_lons)) .* _R .* cosd.(st.cell_lats)
    y = deg2rad.(st.cell_lats .- minimum(st.cell_lats)) .* _R
    m = deepcopy(mesh); m.static_vars["elevation"] = S .* (maximum(x) .- x)
    s = _quiet(() -> initialise_flow_model(m, SWEFlow(); manning_n = n, swe_open_inflow = true))
    fill!(s.ghost_cell_bc, ZeroGradient)
    _setwse!(s, s.elevation .+ h0)
    U = h0^(2 / 3) * sqrt(S) / n
    s.qvec_u .= h0 * U; s.qvec_v .= 0.0
    _run!(s, T)
    bnd = findall(s.boundary_mask)
    reach = 1.3 * (U + sqrt(_G * h0)) * T
    inner = [i for i in eachindex(x) if minimum(hypot(x[i] - x[b], y[i] - y[b]) for b in bnd) > reach]
    @test length(inner) >= 10
    h = s.volume[inner] ./ s.cell_area[inner]; u = s.qvec_u[inner] ./ h
    @info "$lab thin sheet ($(length(inner)) cells): depth $(extrema(h)), u $(extrema(u)) (normal $h0, $(round(U, digits = 5)))"
    tol = lab == "Cartesian" ? 0.005 : 0.10
    @test maximum(abs.(h .- h0)) < tol * h0
    @test maximum(abs.(u .- U)) < tol * U
end

@testset "S8 — eddy viscosity damps a shear layer ($lab)" for (lab, mesh) in _MESHES
    _, ym = _xn(mesh)
    function ke(ce)
        m = deepcopy(mesh); m.static_vars["elevation"] = zeros(length(m.cells))
        s = _quiet(() -> initialise_flow_model(m, SWEFlow(); swe_viscosity = ce, manning_n = 0.03))
        fill!(s.ghost_cell_bc, Closed); _setwse!(s, fill(1.0, length(m.cells)))
        s.qvec_u .= [y > 0.5 ? 0.5 : -0.5 for y in ym]; s.qvec_v .= 0.0
        M0 = sum(s.volume); _run!(s, 600.0)
        return sum(s.cell_area .* (s.qvec_u .^ 2 .+ s.qvec_v .^ 2) ./ max.(s.water_depth, 1e-6)),
               isapprox(sum(s.volume), M0; rtol = 1e-12)
    end
    k0, m0 = ke(0.0); k1, m1 = ke(50.0)          # ν ≈ 2.3 m²/s; at C_E ~ 0.1–1 it is ~100× weaker
    @test m0 && m1
    @test k1 < 0.99 * k0                                            # the viscosity removes shear energy
end

@testset "S7 — local-inertial states unaffected" begin
    s = _quiet(() -> initialise_flow_model(deepcopy(_A5), StandardFlow()))
    @test s.swe === nothing
    f = fieldnames(FlowState)
    s41 = FlowState((getfield(s, k) for k in f[1:41])...)
    @test s41.swe === nothing
    @test_throws ArgumentError step_swe!(s, 1.0)
    @test_throws ArgumentError _quiet(() -> initialise_flow_model(deepcopy(_A5), SWEFlow(); swe_order = 3))
end

@testset "S9 — HDF5 frames carry the SWE velocity, the last frame too" begin
    using HDF5
    m = deepcopy(_MESHES[2][2]); s = _swe(m, _island(m))
    _setwse!(s, 1.0 .+ 0.3 .* _xn(m)[1])                    # sloshing: velocities non-zero
    p = joinpath(mktempdir(), "s9.h5")
    out = SimOutput(path = p, interval = 60.0, enabled = true)
    _write_mesh_metadata!(out, m)
    _quiet(() -> run_simulation!(s, m, 100.0, 10.0, nothing, :none, out; method = SWEFlow()))
    _swe_velocity!(s)
    h5open(p) do f
        ks = sort(collect(keys(f["frames"])))
        @test length(ks) >= 2
        last = f["frames"][ks[end]]
        @test read(last["t"]) == 100.0                       # written after the loop
        @test maximum(abs, read(last["vel_u"])) > 1e-3
        @test read(last["vel_u"]) == s.vel_u && read(last["vel_v"]) == s.vel_v
        @test maximum(abs, read(f["frames"][ks[end - 1]]["vel_u"])) > 1e-3
    end
end

println("\nSWE reference tests complete.")
