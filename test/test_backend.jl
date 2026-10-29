# test/test_backend.jl
# --------------------
# WP5 (2026-10): compute backends (surfacewater/backend.jl, device_run.jl,
# std_device.jl).
#
#   K1  SWE: run_simulation! on KADevice(KernelAbstractions.CPU()) — the
#       device code path, run on the CPU — is bit-identical to CPUThreads
#       (rain, a point source, open boundaries; A5 and Cartesian)
#   K2  standard flow, likewise, for --scheme recommended (η filter,
#       momentum damping, positivity at boundaries) and v0.2 (consistent
#       Q-centring, hybrid limiter)
#   K3  unsupported device options are refused (cell momentum, legacy CFL,
#       SGS)
#   K4  (FLOODA5_TEST_GPU=1 only) the GPU agrees with the CPU after a short
#       wet run: SWE to 1e-12, standard flow to 1e-9 (relative volume)
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_backend.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))
const _R = A5Grid._EARTH_R
_deg(m) = rad2deg(m / _R)
const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _MESHES = [("A5", _quiet(() -> A5Grid.mesh_for_aoi(_AOI, 16))),
                 ("Cartesian", _quiet(() -> generate_cartesian_mesh_bounds(0.0, 0.0, 30.5 * _deg(25.0),
                                                                         30.5 * _deg(25.0), 25.0; lat_ref = 0.0)))]

_xn(m) = (x = [c.center_lon for c in m.cells]; y = [c.center_lat for c in m.cells];
          ((x .- minimum(x)) ./ (maximum(x) - minimum(x)), (y .- minimum(y)) ./ (maximum(y) - minimum(y))))

# a sloping bed with a bump, water in the lower half, a point source upslope
function _state(mesh, method; kw...)
    m = deepcopy(mesh); xm, ym = _xn(m)
    m.static_vars["elevation"] = 2.0 .* (1 .- xm) .+ 0.5 .* exp.(-((xm .- 0.4) .^ 2 .+ (ym .- 0.5) .^ 2) ./ 0.01)
    s = _quiet(() -> initialise_flow_model(m, method; manning_n = 0.03, kw...))
    fill!(s.ghost_cell_bc, ZeroGradient)
    s.volume .= max.(1.0 .- s.elevation, 0.0) .* s.cell_area
    s.water_depth .= s.volume ./ s.cell_area
    k = argmin(abs.(xm .- 0.8) .+ abs.(ym .- 0.5))
    src = InjectionPoint(k, s.cell_ids[k], 5.0, m.cells[k].center_lon, m.cells[k].center_lat)
    return m, s, src
end

function _run!(m, s, src, method, backend; T = 300.0)
    _quiet(() -> run_simulation!(s, m, T, 10.0, nothing, :none; method = method, rainfall_rate = 20e-3 / 3600,
                                 all_sources = AbstractSource[src], backend = backend))
    return s
end

const KACPU = KADevice(KernelAbstractions.CPU())

@testset "K1 — SWE: device path on the CPU = CPU threads ($lab)" for (lab, mesh) in _MESHES
    m, s1, src = _state(mesh, SWEFlow()); s2 = deepcopy(s1)
    _run!(m, s1, src, SWEFlow(), CPUThreads()); _run!(m, s2, src, SWEFlow(), KACPU)
    @test s1.vol_removed > 0.0
    @test s1.volume == s2.volume && s1.qvec_u == s2.qvec_u && s1.qvec_v == s2.qvec_v
    @test s1.water_depth == s2.water_depth && s1.vol_removed == s2.vol_removed
end

@testset "K2 — standard: device path on the CPU = CPU threads ($lab, $sch)" for (lab, mesh) in _MESHES,
                                                                          sch in (:recommended, :v02)
    opts = scheme_options(sch, :standard)
    BOUNDARY_POSITIVITY[] = scheme_boundary_limiter(sch, :standard) === :positivity
    try
        m, s1, src = _state(mesh, StandardFlow(); opts...); s2 = deepcopy(s1)
        _run!(m, s1, src, StandardFlow(), CPUThreads()); _run!(m, s2, src, StandardFlow(), KACPU)
        @test s1.vol_removed > 0.0
        @test s1.volume == s2.volume && s1.edges.flux == s2.edges.flux
        @test s1.qvec_u == s2.qvec_u && s1.qvec_v == s2.qvec_v
        @test s1.water_depth == s2.water_depth && s1.vol_removed == s2.vol_removed
        @test [g.flux_prev for g in s1.ghost_edges] == [g.flux_prev for g in s2.ghost_edges]
    finally
        BOUNDARY_POSITIVITY[] = false
    end
end

@testset "K3 — unsupported device options are refused" begin
    m = _MESHES[2][2]
    for kw in ((momentum_model = :cell,), (cfl_mode = :legacy,))
        _, s, src = _state(m, StandardFlow(); merge(scheme_options(:v02, :standard), kw)...)
        @test_throws ArgumentError DeviceRun(s, StandardFlow(), KACPU, AbstractSource[src])
    end
    _, s, src = _state(m, StandardFlow(); scheme_options(:v02, :standard)...)
    @test_throws ArgumentError DeviceRun(s, SGSFlow(), KACPU, AbstractSource[src])
end

if get(ENV, "FLOODA5_TEST_GPU", "") == "1"
    gpu = _cuda_backend()
    @testset "K4 — GPU agrees with the CPU ($lab, $(typeof(method)))" for (lab, mesh) in _MESHES,
                                                                         method in (SWEFlow(), StandardFlow())
        kw = method isa StandardFlow ? scheme_options(:recommended, :standard) : (;)
        m, s1, src = _state(mesh, method; kw...); s2 = deepcopy(s1)
        _run!(m, s1, src, method, CPUThreads(); T = 60.0)
        Base.invokelatest(_run!, m, s2, src, method, gpu; T = 60.0)
        rel = maximum(abs.(s1.volume .- s2.volume) ./ max.(s1.volume, 1.0))
        @info "$lab $(typeof(method)): max relative volume difference GPU vs CPU $rel"
        @test rel < (method isa SWEFlow ? 1e-12 : 1e-9)
    end
end

println("\nBackend tests complete.")
