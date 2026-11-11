# test/test_plots_render.jl — render every FloodA5Plots (FP) figure type (local only)
#
#   julia --project=. test/test_plots_render.jl [OUTDIR]
#
# Needs CairoMakie (the full project environment); skipped (exit 0) without
# it, so it is not in ci/run_tests.sh. Builds two small format-2 fixture runs
# (Cartesian 25 m, a point inflow on a gentle slope, local-inertial and SWE),
# then renders fieldmap, fieldmaps, difference_map, extent_map, timeseries,
# storage_plot, mass_balance_plot and a composed figure with mark_points! to
# PNG in OUTDIR (default: a temp dir) and checks each file is a PNG.
# The Makie-free logic is tested in CI by test/test_plots.jl.

if Base.find_package("CairoMakie") === nothing || get(ENV, "FLOODA5_PLOTS_NO_MAKIE", "") == "1"
    println("test_plots_render.jl: CairoMakie not available — skipped")
    exit(0)
end

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
include(joinpath(@__DIR__, "..", "visualisation", "FloodA5Plots.jl"))
const FP = FloodA5Plots
const FA = FP.FA
using Test, Logging, CairoMakie

const OUT = isempty(ARGS) ? mktempdir() : mkpath(ARGS[1])
const FIX = mktempdir()
const SRC = (0.0012, 0.0014)                               # (lon, lat)

function fixture_run(model::String)
    aoi, hyd = joinpath(FIX, "aoi.geojson"), joinpath(FIX, "hydro.csv")
    mesh, h5 = joinpath(FIX, "m.parquet"), joinpath(FIX, "$model.h5")
    write(aoi, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
      "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.0036,0.0],[0.0036,0.0027],[0.0,0.0027],[0.0,0.0]]]}}]}""")
    write(hyd, "t_s,Q\n0,3\n1200,3\n")
    args = isfile(mesh) ? ["--meshload", mesh] :
           ["--meshgen-cartesian", aoi, "--dx", "25", "--meshout", mesh]
    append!(args, ["--flow-model", model, "--closed-boundaries",
                   "--inflow-point", "$(SRC[2]),$(SRC[1]),$hyd",
                   "--sim-duration", "1200", "--output", h5, "--output-interval", "120"])
    with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            main(copy(args))
        end
    end
    return FA.load_run(h5; label = model == "standard" ? "cart_li" : "cart_swe")
end

li = fixture_run("standard")
swe = fixture_run("swe")

function check_png(name, fig)
    p = joinpath(OUT, name * ".png")
    save(p, fig; px_per_unit = 1)
    @test isfile(p) && filesize(p) > 1000
    @test read(p, 8) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
    return p
end

@testset "R — every plot type renders" begin
    @test FP.HAS_MAKIE
    @test FP.fieldmap(li) isa Figure
    check_png("fieldmap_depth", FP.fieldmap(li; t = 600))
    check_png("fieldmap_speed_local", FP.fieldmap(swe, 5; var = :speed, coords = :local, origin = :source))
    check_png("fieldmap_wse_grey", FP.fieldmap(li; var = :wse, basemap = :grey, colorrange = (0, 0.5)))
    check_png("fieldmaps", FP.fieldmaps([li, swe], [300.0, 600.0, 1200.0]; coords = :local, origin = :source))
    check_png("fieldmaps_one_run", FP.fieldmaps(li, [600.0, 1200.0]))
    check_png("difference_map", FP.difference_map(swe, li; t = 1200.0))
    check_png("extent_map", FP.extent_map(swe, li; t = 1200.0, threshold = 0.01))
    check_png("timeseries_runs", FP.timeseries([li, swe], [(lon = SRC[1], lat = SRC[2], label = "source"), 1]))
    check_png("timeseries_one_run", FP.timeseries(li, [(SRC[1], SRC[2]), (0.0025, 0.0014)]))
    check_png("storage_plot", FP.storage_plot([li, swe]))
    check_png("mass_balance", FP.mass_balance_plot(li))
    # composed: two panels sharing a projection, points marked, themed
    FP.activate!()
    fig = Figure(size = (900, 420))
    p1 = FP.fieldmap!(fig[1, 1], li; t = 1200.0, coords = :local, origin = :source, colorbar = false)
    p2 = FP.fieldmap!(fig[1, 2], swe; t = 1200.0, coords = p1.coords, colorrange = p1.range)
    FP.mark_points!(p1, [(SRC[1] + 0.0005, SRC[2])]; labels = ["gauge"])
    @test p1.coords === p2.coords && p2.range == p1.range
    check_png("composed", fig)
    set_theme!()
    # the extent map's counts are FA.compare's
    c = FA.compare(swe, li; threshold = 0.01)
    s = FP.extent_scores(FP.extent_classes(FA.depth(swe, c.k_ref), FA.depth(li, c.k_run); threshold = 0.01))
    @test (s.hits, s.misses, s.false_alarms, s.csi) == (c.hits, c.misses, c.false_alarms, c.csi)
    # storage of a closed domain: inflow volume 3 m³/s × t
    @test FP.storage_series(li)[end] ≈ 3.0 * li.t[end] rtol = 0.02
    mb = FP.mass_balance_series(li)
    @test mb.stored[end] ≈ mb.net rtol = 1e-6
    @test all(≈(3.0; rtol = 0.05), mb.rate[2:end])
end
println("figures in ", OUT)
