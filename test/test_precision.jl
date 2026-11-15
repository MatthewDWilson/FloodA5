# test/test_precision.jl
# ----------------------
# 0.6.0 (2026-11): `--precision float32` on the device backends
# (surfacewater/precision.jl): Float32 fluxes and momentum, Float64 volume,
# water levels relative to a local datum (RelLevels). Run on the device code
# path on the CPU (KADevice(KernelAbstractions.CPU())), so CI needs no GPU.
#
#   P1  `@fp` and `_to_precision`: literals and named constants take the
#       kernel's float type; Float64 arguments are passed through untouched
#       and the volume arrays stay Float64 in a Float32 tuple
#   P2  RelLevels in Float64 = the original arithmetic (AbsLevels) to
#       round-off: SWE, and standard flow on the Cartesian mesh. On A5 the
#       standard solver amplifies round-off (HYDRAULICS.md §15.5), so the
#       two agree there only as two Float64 runs with different round-off do
#   P3  Float32 against the Float64 CPU run, plan of record §7.3: max |Δh| <
#       1e-4 m, relative mass error < 1e-6, wet-extent CSI > 0.99 (SWE both
#       meshes and orders; standard flow, recommended and v0.2, Cartesian);
#       on A5 standard flow: mass error and CSI, and depth within the bound
#       of two Float64 runs
#   P4  the result does not depend on the datum: the same case 1000 m higher
#       gives the same Float32 depths (it would not with levels z + h, which
#       Float32 resolves to 6e-5 m there)
#   P5  a lake at rest over an uneven bed at z ≈ 1000 m stays at rest in
#       Float32 (SWE, and standard flow with the diamond gradient and the η
#       filter)
#   P6  refusals: Float32 with AbsLevels, with SGS, on the CPU-thread backend
#   P7  the command line: --precision float32 --backend ka-cpu is recorded in
#       the output and --rerun reproduces it bit for bit; float32 without a
#       device backend is refused
#   P8  (FLOODA5_TEST_GPU=1 only) the GPU in Float32 against the Float64 CPU
#       run, as P3
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_precision.jl

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))

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

# a sloping bed with a bump at datum z0, water in the lower half, a point source upslope
function _state(mesh, method; z0 = 0.0, kw...)
    m = deepcopy(mesh); xm, ym = _xn(m)
    m.static_vars["elevation"] = z0 .+ 2.0 .* (1 .- xm) .+ 0.5 .* exp.(-((xm .- 0.4) .^ 2 .+ (ym .- 0.5) .^ 2) ./ 0.01)
    s = _quiet(() -> initialise_flow_model(m, method; manning_n = 0.03, kw...))
    fill!(s.ghost_cell_bc, ZeroGradient)
    s.volume .= max.(z0 + 1.0 .- s.elevation, 0.0) .* s.cell_area
    s.water_depth .= s.volume ./ s.cell_area
    k = argmin(abs.(xm .- 0.8) .+ abs.(ym .- 0.5))
    src = InjectionPoint(k, s.cell_ids[k], 5.0, m.cells[k].center_lon, m.cells[k].center_lat)
    return m, s, src
end

const RAIN = 20e-3 / 3600
function _run!(m, s, src, method, backend; T = 300.0)
    _quiet(() -> run_simulation!(s, m, T, 10.0, nothing, :none; method = method, rainfall_rate = RAIN,
                                 all_sources = AbstractSource[src], backend = backend))
    return s
end
# conservation: what was there, the source and the rain, against what is stored and what left
_mass_err(s, v0, T) = abs(sum(s.volume) + s.vol_removed - v0 - 5.0 * T - RAIN * T * sum(s.cell_area)) /
                      (v0 + 5.0 * T + RAIN * T * sum(s.cell_area))
_dh(a, b) = maximum(abs.(a.water_depth .- b.water_depth))
function _csi(a, b; h = 0.01)
    wa = a.water_depth .> h; wb = b.water_depth .> h
    return count(wa .& wb) / max(count(wa .| wb), 1)
end

const KA = KernelAbstractions.CPU()
const REL64 = KADevice(KA; levels = RelLevels())
const F32 = KADevice(KA; precision = Float32)

_std(sch) = (StandardFlow(), scheme_options(sch, :standard), scheme_boundary_limiter(sch, :standard) === :positivity)
const _CONFIGS = [("SWE order 2", SWEFlow(), (;), false), ("SWE order 1", SWEFlow(), (swe_order = 1,), false),
                  ("standard recommended", _std(:recommended)...), ("standard v0.2", _std(:v02)...)]

# CPU-thread reference and a device run of the same case
function _pair(mesh, method, opts, bpos, backend; z0 = 0.0, T = 300.0)
    BOUNDARY_POSITIVITY[] = bpos
    try
        m, s0, src = _state(mesh, method; z0 = z0, opts...)
        v0 = sum(s0.volume)
        ref = _run!(m, deepcopy(s0), src, method, CPUThreads(); T = T)
        dev = _run!(m, deepcopy(s0), src, method, backend; T = T)
        return ref, dev, v0
    finally
        BOUNDARY_POSITIVITY[] = false
    end
end

@testset "P1 — @fp and _to_precision" begin
    f(x::T) where {T} = @fp T (0.5 * x + _G * 1e-3)
    @test f(2.0f0) isa Float32 && f(2.0) isa Float64
    @test f(2.0) === 0.5 * 2.0 + _G * 1e-3
    @test f(2.0f0) === 0.5f0 * 2.0f0 + Float32(_G) * 1f-3
    g(x::T) where {T} = @fp T (isnan(x) ? NaN : x > HFLOW_THRESHOLD ? Inf : 0.0)
    @test g(1.0f0) === Inf32 && g(0.0f0) === 0.0f0 && g(NaN32) === NaN32 && g(1.0) === Inf
    # the rewrite touches Float64 literals and the listed constants only
    ex = macroexpand(Main, :(@fp T (a.NaN + 2 * x^2 + 1.5 + f(y; k = 0.25))))
    @test ex == :(a.NaN + 2 * x^2 + convert(T, 1.5) + f(y; k = convert(T, 0.25)))

    V = [1.0, 2.0]; z = [0.5, 0.25]
    a = (volume = V, z = z, n = 2, theta = 0.9, flags = Bool[true], idx = [1, 2],
         nested = (w = [1.0 2.0; 3.0 4.0],), lev = RelLevels())
    @test _to_precision(Float64, a) === a
    b = _to_precision(Float32, a)
    @test b.volume === V                                     # the volume stays Float64, and shared
    @test b.z isa Vector{Float32} && b.z == Float32.(z)
    @test b.theta === 0.9f0 && b.n === 2 && b.idx === a.idx && b.flags === a.flags
    @test b.nested.w isa Matrix{Float32} && b.lev === RelLevels()
    @test _precision(CPUThreads()) === Float64 && _precision(F32) === Float32
    @test _levels(KADevice(KA)) === AbsLevels() && _levels(F32) === RelLevels() && _levels(REL64) === RelLevels()
end

@testset "P2 — RelLevels in Float64 = AbsLevels to round-off ($lab, $name)" for (lab, mesh) in _MESHES,
        (name, method, opts, bpos) in _CONFIGS
    ref, dev, v0 = _pair(mesh, method, opts, bpos, REL64)
    amplifies = lab == "A5" && method isa StandardFlow      # HYDRAULICS.md §15.5
    @test _dh(ref, dev) < (amplifies ? 1e-3 : 1e-10)
    @test _mass_err(dev, v0, 300.0) < 1e-12
    @test abs(dev.vol_removed - ref.vol_removed) <= (amplifies ? 1e-3 : 1e-9) * ref.vol_removed
    @test ref.vol_removed > 0.0
end

@testset "P3 — Float32 against Float64, §7.3 tolerances ($lab, $name)" for (lab, mesh) in _MESHES,
        (name, method, opts, bpos) in _CONFIGS
    ref, dev, v0 = _pair(mesh, method, opts, bpos, F32)
    amplifies = lab == "A5" && method isa StandardFlow
    d = _dh(ref, dev)
    @info "$lab, $name: Float32 vs Float64 max |Δh| = $d m, mass error $(_mass_err(dev, v0, 300.0)), CSI $(_csi(ref, dev))"
    @test d < (amplifies ? 1e-3 : 1e-4)
    @test _mass_err(dev, v0, 300.0) < 1e-6
    @test _csi(ref, dev) > 0.99
    @test all(isfinite, dev.volume) && all(>=(0.0), dev.volume)
end

@testset "P4 — Float32 does not depend on the datum ($lab, $name)" for (lab, mesh) in _MESHES,
        (name, method, opts, bpos) in _CONFIGS[[1, 3]]
    _, lo, _ = _pair(mesh, method, opts, bpos, F32; z0 = 0.0, T = 120.0)
    _, hi, _ = _pair(mesh, method, opts, bpos, F32; z0 = 1000.0, T = 120.0)
    # the bed offsets are differenced in Float64 before they are rounded, so they carry
    # the datum only through Float64 round-off (1e-13 m at 1000 m)
    @test _dh(lo, hi) < (lab == "A5" && method isa StandardFlow ? 1e-4 : 1e-6)
end

@testset "P5 — Float32 lake at rest at z ≈ 1000 m ($lab, $name)" for (lab, mesh) in _MESHES,
        (name, method, opts, bpos) in _CONFIGS[[1, 3]]
    m = deepcopy(mesh); xm, ym = _xn(m)
    m.static_vars["elevation"] = 1000.0 .+ 0.4 .* sin.(7 .* xm) .* cos.(5 .* ym) .+ 0.3 .* xm
    s = _quiet(() -> initialise_flow_model(m, method; manning_n = 0.03, opts...))
    fill!(s.ghost_cell_bc, Closed)
    s.volume .= (1001.5 .- s.elevation) .* s.cell_area       # level surface, 0.8–1.9 m deep
    s.water_depth .= s.volume ./ s.cell_area
    h0 = copy(s.water_depth); v0 = sum(s.volume)
    _quiet(() -> run_simulation!(s, m, 300.0, 10.0, nothing, :none; method = method, backend = F32))
    @test maximum(abs.(s.water_depth .- h0)) < 2e-6           # round-off of a 1–2 m depth in Float32
    @test abs(sum(s.volume) - v0) < 1e-9 * v0
    q = method isa SWEFlow ? max(maximum(abs, s.qvec_u), maximum(abs, s.qvec_v)) : maximum(abs, s.edges.flux)
    # Float32 round-off of the heads (1e-7 m on a 1–2 m depth) drives discharges of
    # order 1e-5 m²/s on A5 — velocities of 1e-5 m/s; Float64 gives exactly none
    @test q < 1e-4
end

@testset "P6 — refusals" begin
    @test_throws ArgumentError KADevice(KA; precision = Float32, levels = AbsLevels())
    @test_throws ArgumentError KADevice(KA; precision = Float16)
    m, s, src = _state(_MESHES[2][2], StandardFlow(); scheme_options(:recommended, :sgs)...)
    @test_throws ArgumentError DeviceRun(s, SGSFlow(), F32, AbstractSource[src])
    @test_throws ArgumentError DeviceRun(s, SGSFlow(), REL64, AbstractSource[src])
    @test _check_precision(:float64, :cpu, :sgs) === Float64
    @test _check_precision(:float32, :gpu, :standard) === Float32
    @test _check_precision(:float32, :ka_cpu, :swe) === Float32
    @test_throws ArgumentError _check_precision(:float32, :cpu, :standard)
    @test_throws ArgumentError _check_precision(:float32, :gpu, :sgs)
    @test_throws ArgumentError _check_precision(:float16, :gpu, :standard)
end

@testset "P7 — command line: recorded, and --rerun reproduces it" begin
    dir = mktempdir()
    aoi = joinpath(dir, "aoi.geojson")
    write(aoi, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
      "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.0027,0.0],[0.0027,0.0027],[0.0,0.0027],[0.0,0.0]]]}}]}""")
    run_cli(args) = with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            main(copy(args))
        end
    end
    same(a, b) = FO.frame_times(a) == FO.frame_times(b) &&
                 all(v -> isequal(FO.read_var(a, v), FO.read_var(b, v)), FO.frame_vars(a))
    base = ["--meshgen-cartesian", aoi, "--dx", "25", "--meshout", joinpath(dir, "m.parquet"),
            "--flow-model", "standard", "--rainfall", "40", "--sim-duration", "300", "--output-interval", "60"]
    o32 = joinpath(dir, "f32.h5"); o64 = joinpath(dir, "f64.h5"); r32 = joinpath(dir, "r32.h5")
    run_cli(vcat(base, ["--backend", "ka-cpu", "--precision", "float32", "--output", o32]))
    run_cli(vcat(base, ["--backend", "ka-cpu", "--output", o64]))
    ri = FO.run_info(o32)
    @test ri["precision"] == "float32" && ri["options"]["precision"] == "float32"
    @test FO.run_info(o64)["precision"] == "float64"
    @test ri["summary"]["completed"] == true && ri["summary"]["wall_first5_steps_s"] > 0
    @test abs(ri["summary"]["conservation_error_m3"]) < 1e-6 * ri["summary"]["vol_applied_m3"]
    h32 = FO.read_frame(o32, "water_depth"); h64 = FO.read_frame(o64, "water_depth")
    @test eltype(h32) == Float64                             # output files are Float64
    @test h32 != h64 && maximum(abs.(h32 .- h64)) < 1e-4
    run_cli(["--rerun", o32, "--output", r32])
    @test same(o32, r32) && FO.run_info(r32)["precision"] == "float32"
    # float32 without a device backend is refused before any work is done
    @test_throws ArgumentError run_flood_model(; mesh_source = (:load, joinpath(dir, "m.parquet"), nothing),
                                               flow_method = :standard, precision = :float32,
                                               sim_duration = 10.0)
end

if get(ENV, "FLOODA5_TEST_GPU", "") == "1"
    gpu32 = _cuda_backend(; precision = Float32)
    @testset "P8 — GPU Float32 against Float64 CPU ($lab, $name)" for (lab, mesh) in _MESHES,
            (name, method, opts, bpos) in _CONFIGS[[1, 3]]
        ref, dev, v0 = Base.invokelatest(_pair, mesh, method, opts, bpos, gpu32)
        amplifies = lab == "A5" && method isa StandardFlow
        @info "$lab, $name: GPU Float32 vs CPU Float64 max |Δh| = $(_dh(ref, dev)) m"
        @test _dh(ref, dev) < (amplifies ? 1e-3 : 1e-4)
        @test _mass_err(dev, v0, 300.0) < 1e-6
        @test _csi(ref, dev) > 0.99
    end
end

println("\nPrecision tests complete.")
