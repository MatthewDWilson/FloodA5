# test/test_default_scheme.jl
# ---------------------------
# --scheme recommended|v0.2|legacy (WP2b, 2026-10; 0.3.0): the default
# numerical options.
#
#   D1  scheme_options and scheme_boundary_limiter for every flow model and
#       scheme; every scheme's options are accepted by initialise_flow_model
#   D2  CLI end to end on a small Cartesian mesh: the recommended (0.3.0)
#       options by default, the 0.2.x options with --scheme v0.2, the
#       pre-2026-10 options with --scheme legacy, explicit flags overriding
#       the scheme, and the scheme's momentum damping dropped when explicit
#       flags remove what it needs
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_default_scheme.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test

const _LEGACY = (gradient_correction = false, face_flux_method = :legacy,
                 friction_model = :component, limiter = :legacy, cfl_mode = :legacy,
                 momentum_model = :edge, q_centre_mode = :legacy, q_centre_theta = 0.9,
                 eta_filter = 0.0, momentum_damping = 0.0, momentum_damping_scale = :wave)

@testset "D1 — scheme_options" begin
    r = scheme_options(:recommended, :standard)
    @test r == merge(_LEGACY, (gradient_correction = true, face_flux_method = :diamond_wetdry,
                               friction_model = :vector, limiter = :positivity, cfl_mode = :local,
                               q_centre_mode = :consistent, q_centre_theta = 1.0,
                               eta_filter = 0.05, momentum_damping = 0.05))
    v = scheme_options(:v02, :standard)
    @test v == merge(_LEGACY, (gradient_correction = true, face_flux_method = :diamond_wetdry,
                               friction_model = :vector, limiter = :hybrid, cfl_mode = :local,
                               q_centre_mode = :consistent))
    for s in (:recommended, :v02)                     # SGS: unchanged in 0.3.0
        g = scheme_options(s, :sgs)
        @test g == merge(_LEGACY, (friction_model = :vector, cfl_mode = :local, q_centre_theta = 1.0))
    end
    for f in (:standard, :sgs, :swe), s in (:legacy, :recommended, :v02)
        (s === :legacy || f === :swe) && @test scheme_options(s, f) == _LEGACY
    end
    @test scheme_boundary_limiter(:recommended, :standard) === :positivity
    for (s, f) in ((:recommended, :sgs), (:recommended, :swe), (:v02, :standard), (:legacy, :standard))
        @test scheme_boundary_limiter(s, f) === :legacy
    end
    @test_throws ArgumentError scheme_options(:fast, :standard)
    # Every scheme's options are accepted by initialise_flow_model
    dlat = rad2deg(20.0 / A5Grid._EARTH_R)
    mesh = generate_cartesian_mesh_bounds(0.0, 0.0, 4.5dlat, 4.5dlat, 20.0; lat_ref = 0.0)
    st = initialise_flow_model(mesh, StandardFlow(); r...)
    @test st.friction_model === :vector && st.limiter === :positivity && st.cfl_mode === :local
    @test st.q_centre_theta == 1.0 && st.eta_filter_cs == 0.05 && st.mom_damping_cs == 0.05
    st = initialise_flow_model(mesh, StandardFlow(); v...)
    @test st.limiter === :hybrid && st.q_centre_mode === :consistent && st.eta_filter_cs == 0.0
end

# D2 — the CLI resolves unset options from the scheme
const AOI = joinpath(mktempdir(), "aoi.geojson")
write(AOI, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
  "geometry":{"type":"Polygon","coordinates":[[[0,0],[0.0015,0],[0.0015,0.0015],[0,0.0015],[0,0]]]}}]}""")

function _cli_log(extra...)
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(joinpath(@__DIR__, "..", "FloodModel.jl"))
           --meshgen-cartesian $AOI --dx 25 --flow-model standard --injection-point 0.0007,0.0007,1.0
           --sim-duration 20 --dt-max 2 $extra`
    out = IOBuffer()
    run(pipeline(addenv(cmd, "FLOODA5_NO_VIS" => "1"); stdout = out, stderr = out))
    log = String(take!(out))
    i = findfirst("Flow model ready", log)
    @assert i !== nothing "no 'Flow model ready' line in the CLI output:\n$log"
    return log[i.start:min(end, i.start + 900)], log
end

@testset "D2 — CLI uses the scheme for unset options" begin
    rec, _ = _cli_log()
    @test occursin("gradient_correction = true", rec)
    @test occursin("face_flux_method = diamond_wetdry", rec)
    @test occursin("friction_model = vector", rec)
    @test occursin("limiter = positivity, cfl_mode = local", rec)
    @test occursin("scheme = recommended", rec)
    @test occursin("q_centre_theta = 1.0", rec)
    @test occursin("eta_filter = 0.05", rec) && occursin("momentum_damping = 0.05 (wave)", rec)
    @test occursin("boundary_limiter = positivity", rec)
    # The 0.2.x default
    v02, _ = _cli_log("--scheme", "v0.2")
    @test occursin("limiter = hybrid, cfl_mode = local", v02) && occursin("scheme = v02", v02)
    @test occursin("q_centre_mode = consistent", v02) && occursin("q_centre_theta = 0.9", v02)
    @test !occursin("eta_filter", v02) && !occursin("momentum_damping", v02)
    @test !occursin("boundary_limiter", v02)
    # Legacy scheme with one explicit override: the override wins, every other
    # unset option takes the legacy value.
    mix, _ = _cli_log("--scheme", "legacy", "--friction-model", "vector")
    @test occursin("gradient_correction = false", mix) && occursin("face_flux_method = legacy", mix)
    @test occursin("limiter = legacy, cfl_mode = legacy", mix)
    @test occursin("friction_model = vector", mix)
    @test occursin("scheme = legacy", mix)
    @test occursin("q_centre_mode = legacy", mix)
    @test !occursin("eta_filter", mix) && !occursin("momentum_damping", mix)
    # An explicit flag that removes what the scheme's momentum damping needs
    # switches the damping off (with a notice); the η filter stays.
    off, full = _cli_log("--face-flux-method", "legacy")
    @test occursin("face_flux_method = legacy", off) && !occursin("momentum_damping", off)
    @test occursin("eta_filter = 0.05", off)
    @test occursin("Momentum damping of --scheme recommended switched off", full)
    # ... but an explicit --momentum-damping is honoured
    exp, _ = _cli_log("--momentum-damping", "0.1", "--eta-filter", "0")
    @test occursin("momentum_damping = 0.1 (wave)", exp) && !occursin("eta_filter", exp)
end

println("\nDefault-scheme tests complete.")
