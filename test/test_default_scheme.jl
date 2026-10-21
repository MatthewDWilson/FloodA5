# test/test_default_scheme.jl
# ---------------------------
# --scheme recommended|legacy (WP2b, 2026-10): the default numerical options.
#
#   D1  scheme_options values for both flow models and both schemes
#   D2  CLI end to end on a small Cartesian mesh: the recommended options by
#       default, the pre-2026-10 options with --scheme legacy, and explicit
#       flags overriding the scheme
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_default_scheme.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test

@testset "D1 — scheme_options" begin
    r = scheme_options(:recommended, :standard)
    @test r == (gradient_correction = true, face_flux_method = :diamond_wetdry,
                friction_model = :vector, limiter = :hybrid, cfl_mode = :local,
                momentum_model = :edge)
    g = scheme_options(:recommended, :sgs)
    @test g.friction_model === :vector && g.cfl_mode === :local
    @test g.limiter === :legacy && !g.gradient_correction && g.face_flux_method === :legacy
    for f in (:standard, :sgs)
        l = scheme_options(:legacy, f)
        @test l == (gradient_correction = false, face_flux_method = :legacy,
                    friction_model = :component, limiter = :legacy, cfl_mode = :legacy,
                    momentum_model = :edge)
    end
    @test_throws ArgumentError scheme_options(:fast, :standard)
    # Every scheme's options are accepted by initialise_flow_model
    dlat = rad2deg(20.0 / A5Grid._EARTH_R)
    mesh = generate_cartesian_mesh_bounds(0.0, 0.0, 4.5dlat, 4.5dlat, 20.0; lat_ref = 0.0)
    st = initialise_flow_model(mesh, StandardFlow(); r...)
    @test st.friction_model === :vector && st.limiter === :hybrid && st.cfl_mode === :local
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
    return log[i.start:min(end, i.start + 600)]
end

@testset "D2 — CLI uses the scheme for unset options" begin
    rec = _cli_log()
    @test occursin("gradient_correction = true", rec)
    @test occursin("face_flux_method = diamond_wetdry", rec)
    @test occursin("friction_model = vector", rec)
    @test occursin("limiter = hybrid", rec) && occursin("cfl_mode = local", rec)
    @test occursin("scheme = recommended", rec)
    # Legacy scheme with one explicit override: the override wins, every other
    # unset option takes the legacy value.
    mix = _cli_log("--scheme", "legacy", "--friction-model", "vector")
    @test occursin("gradient_correction = false", mix) && occursin("face_flux_method = legacy", mix)
    @test occursin("limiter = legacy", mix) && occursin("cfl_mode = legacy", mix)
    @test occursin("friction_model = vector", mix)
    @test occursin("scheme = legacy", mix)
end

println("\nDefault-scheme tests complete.")
