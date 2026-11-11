# test/test_plots.jl — the Makie-free logic of visualisation/FloodA5Plots.jl (FP)
#
#   julia --project=ci test/test_plots.jl
#
# CI has no CairoMakie: this loads FP without it (FLOODA5_PLOTS_NO_MAKIE=1)
# and tests colour ranges, frame selection, extent classification, map
# coordinates, sources and run styles on synthetic runs (no output files).
# Rendering every plot type is test/test_plots_render.jl (local only).

ENV["FLOODA5_PLOTS_NO_MAKIE"] = "1"
include(joinpath(@__DIR__, "..", "visualisation", "FloodA5Plots.jl"))
const FP = FloodA5Plots
const FA = FP.FA
using Test, Statistics

"A synthetic run (no file): an n × n grid of square cells of `dx` degrees at (lon0, lat0)."
function fake_run(; label = "fake.h5", t = collect(0.0:600.0:3600.0), n = 4, dx = 0.001,
                  lon0 = -2.9, lat0 = 54.9, scheme = "a5-v1", flow_model = "standard",
                  info = Dict{String,Any}(), z = nothing)
    polys = Matrix{Float64}[]; lon = Float64[]; lat = Float64[]
    for j in 1:n, i in 1:n
        x0 = lon0 + (i - 1) * dx; y0 = lat0 + (j - 1) * dx
        push!(polys, [x0 x0+dx x0+dx x0 x0; y0 y0 y0+dx y0+dx y0])
        push!(lon, x0 + dx / 2); push!(lat, y0 + dx / 2)
    end
    m = length(lon)
    mesh = FA.MeshGeom(["c$i" for i in 1:m], lon, lat, z === nothing ? zeros(m) : z, fill(100.0, m),
                       polys, zeros(Int32, 0, m), Dict{String,Vector{Float64}}(), scheme)
    return FA.Run("", label, 2, t, mesh, info, flow_model)
end

@testset "P1 — module loads without CairoMakie" begin
    @test !FP.HAS_MAKIE
    @test !isdefined(FP, :fieldmap)
    @test length(FP.SERIES) == 8 && all(startswith("#"), FP.SERIES)
    @test FP.MESH_STYLE["a5"].color == FP.SERIES[1] && FP.MESH_STYLE["cart"].color == FP.SERIES[2]
    @test FP.series_color(3) == FP.SERIES[3]
    @test_throws ErrorException FP.series_color(9)
end

@testset "P2 — colour ranges" begin
    @test FP.nice_ceil(3.7) == 4.0
    @test FP.nice_ceil(0.37) ≈ 0.4
    @test FP.nice_ceil(1.0) == 1.0
    @test FP.nice_ceil(1.01) == 1.2
    @test FP.nice_ceil(9.0) == 10.0
    @test FP.nice_ceil(4.2e6) == 5.0e6
    @test FP.nice_ceil(0.0) == 0.0 && FP.nice_ceil(NaN) == 0.0
    # depth with one huge outlier: the top is the 98th percentile, rounded up
    h = vcat(collect(range(0.0, 1.9; length = 999)), 1000.0)
    lo, hi = FP.robust_range(h; lo = 0.0)
    @test lo == 0.0 && hi == 2.0
    # mask: dry cells excluded; NaN ignored
    h2 = [0.0, 0.0, 0.5, NaN, 0.8]
    @test FP.robust_range(h2; mask = h2 .> 0.01, lo = 0.0, q = 1.0) == (0.0, 0.8)
    # symmetric difference
    d = vcat(fill(-0.3, 50), fill(0.25, 50), [5.0])
    @test FP.robust_range(d; symmetric = true) == (-0.3, 0.3)
    # quantile range (elevation), not rounded
    z = collect(1.0:100.0)
    a, b = FP.robust_range(z; nice = false)
    @test a ≈ quantile(z, 0.02) && b ≈ quantile(z, 0.98)
    # degenerate inputs still give valid ranges
    @test FP.robust_range(Float64[]) == (0.0, 1.0)
    @test FP.robust_range([NaN, NaN]; lo = 0.0) == (0.0, 1.0)
    a, b = FP.robust_range(fill(3.0, 10))
    @test a < 3.0 < b
    a, b = FP.robust_range(zeros(10); lo = 0.0)
    @test a == 0.0 && b > 0
    a, b = FP.robust_range(zeros(10); symmetric = true)
    @test a == -b && b > 0
    @test FP.default_range(:depth, [0.0, 1.3]; mask = nothing)[1] == 0.0
    @test FP.default_range(:elevation, collect(1.0:100.0))[1] > 1.0
    # colour index: clamped at both ends, NaN not drawn
    @test FP.color_index(0.0, (0.0, 1.0), 10) == 1
    @test FP.color_index(1.0, (0.0, 1.0), 10) == 10
    @test FP.color_index(0.55, (0.0, 1.0), 10) == 6
    @test FP.color_index(-5.0, (0.0, 1.0), 10) == 1
    @test FP.color_index(5.0, (0.0, 1.0), 10) == 10
    @test FP.color_index(NaN, (0.0, 1.0), 10) == 0
end

@testset "P3 — variables, labels, times" begin
    @test FP.var_label(:depth) == "Depth (m)"
    @test FP.var_label(:speed) == "Speed (m/s)"
    @test FP.var_label("saturation") == "saturation"
    @test FP.var_spec(:depth).masked && !FP.var_spec(:elevation).masked
    @test FP.hhmm(7210.0) == "02:00"
    @test FP.hhmm(1082.0) == "00:18"
    @test FP.hhmm(0.0) == "00:00"
    @test FP.hhmm(36 * 3600 + 29) == "36:00"
    @test FP.time_label(5400) == "t = 01:30"
    @test FP.time_axis(3600.0) == (60.0, "Time (min)")
    @test FP.time_axis(21600.0) == (3600.0, "Time (h)")
    @test FP.run_label(fake_run(label = "a5_li.h5")) == "a5_li"
    @test FP.run_label(fake_run(label = "my run")) == "my run"
    @test FP.volume_scale(5.9e6) == (1e6, "10⁶ m³")
    @test FP.volume_scale(2.1e4) == (1e3, "10³ m³")
    @test FP.volume_scale(900.0) == (1.0, "m³")
    @test FP.format_volume(5.92e6) == "5.92 × 10⁶ m³"
    @test FP.format_volume(212.0) == "212 m³"
    @test FP.format_volume(NaN) == "–"
end

@testset "P4 — frame selection" begin
    # frame times as Carlisle's: nominally every 900 s, actually a little later
    a = fake_run(t = [10.0; [900.0k + 10.0 + 2k for k in 1:24]])
    b = fake_run(t = collect(0.0:900.0:21600.0))
    @test FP.frame_tol(b) == 540.0
    K = FP.frame_grid([a, b], [2, 4, 6] .* 3600)
    @test size(K) == (2, 3)
    @test K[2, :] == [9, 17, 25]
    @test all(abs(a.t[K[1, j]] - [2, 4, 6][j] * 3600) < 100 for j in 1:3)
    @test_throws ErrorException FP.frame_grid([b], [30_000.0])         # past the end
    @test FP.frame_grid([b], [30_000.0]; tol = Inf) == fill(25, 1, 1)
    # matched pairs: by time, by the run's frame, or the last common time
    c = fake_run(t = collect(0.0:600.0:7200.0))
    @test FP.matched_pair(b, c; t = 3600.0) == (5, 7)
    @test FP.matched_pair(b, c; k = 4) == (3, 4)                        # c.t[4] = 1800
    @test FP.matched_pair(b, c) == (9, 13)                              # last common t = 7200
    @test FP.same_mesh(b, c)
    @test_throws ErrorException FP.same_mesh(b, fake_run(n = 3))
end

@testset "P5 — extent classes and differences" begin
    h_ref = [0.0, 0.5, 0.5, 0.0, NaN, 0.2, 0.05]
    h     = [0.0, 0.5, 0.0, 0.5, 0.5, NaN, 0.15]
    c = FP.extent_classes(h_ref, h; threshold = 0.1)
    @test c == Int8[0, 1, 2, 3, 3, 2, 3]
    s = FP.extent_scores(c)
    cf = FA.compare_fields(h_ref, h, ones(length(h)); threshold = 0.1)
    @test (s.hits, s.misses, s.false_alarms) == (cf.hits, cf.misses, cf.false_alarms)
    @test s.csi == cf.csi
    @test FP.extent_scores(zeros(Int8, 4)).csi == 0.0
    @test_throws ErrorException FP.extent_classes([0.0], [0.0, 1.0])
    d = FP.difference_values(h_ref, h, h_ref, h; wet = 0.01)
    @test isnan(d[1])                  # dry in both
    @test d[2] == 0.0 && d[3] == -0.5 && d[4] == 0.5
    @test isnan(d[5]) && isnan(d[6])   # a NaN in either
    @test d[7] ≈ 0.1
end

@testset "P6 — map coordinates" begin
    r = fake_run(n = 4, dx = 0.001, lon0 = 10.0, lat0 = 60.0)
    mc = FP.map_coords(r)                                   # :lonlat, bbox centre
    @test mc.coords === :lonlat && mc.lon0 ≈ 10.002 && mc.lat0 ≈ 60.002
    @test FP.project(mc, 1.0, 2.0) == (1.0, 2.0)
    @test FP.axis_labels(mc) == ("Longitude (°)", "Latitude (°)")
    P = FP.project_polygons(mc, r.mesh)
    @test P == r.mesh.polygons
    # true aspect: a 0.004° square at 60° N is twice as tall as wide
    xl, yl, asp = FP.map_extent(mc, P; pad = 0.0)
    @test xl == (10.0, 10.004) && yl == (60.0, 60.004)
    @test asp ≈ cosd(60.002) rtol = 1e-9
    xl, yl, asp = FP.map_extent(mc, P; pad = 0.05)
    @test asp ≈ (0.004 * cosd(60.002) + 0.1 * 0.004) / (0.004 + 0.1 * 0.004) rtol = 1e-9
    # local metres equal FA.local_xy about the origin
    ml = FP.map_coords(r; coords = :local, origin = (10.0, 60.0))
    x, y = FP.project(ml, r.mesh.lon, r.mesh.lat)
    xf, yf = FA.local_xy(r.mesh, 10.0, 60.0)
    @test x ≈ xf && y ≈ yf
    @test FP.axis_labels(ml) == ("x (m)", "y (m)")
    Pl = FP.project_polygons(ml, r.mesh; cells = [1])
    @test length(Pl) == 1
    @test Pl[1][:, 1] ≈ [0.0, 0.0] atol = 1e-9
    _, _, asp = FP.map_extent(ml, FP.project_polygons(ml, r.mesh); pad = 0.0)
    @test asp ≈ cosd(60.0) rtol = 1e-9
    # origin = :source uses the first source, :centre otherwise
    rs = fake_run(info = Dict{String,Any}("options" => Dict{String,Any}(
        "inflow_points" => [["60.003", "10.001", "C:\\x\\q10.csv", "q10"]])))
    m2 = FP.map_coords(rs; coords = :local, origin = :source)
    @test (m2.lon0, m2.lat0) == (10.001, 60.003)
    @test FP.map_coords(r; origin = :source).lon0 == mc.lon0      # no sources: the centre
    @test_throws ErrorException FP.map_coords(r; coords = :utm)
    @test_throws ErrorException FP.map_coords(r; origin = :corner)
    @test FP.map_coords(mc) === mc
end

@testset "P7 — sources and points" begin
    o = Dict{String,Any}(
        "inflow_points" => [["51.0", "-0.0284", "D:\\c\\inflow_q10.csv", "inflow_q10"],
                            ["51.1", "-0.03", "/tmp/hyd.csv"]],
        "injection_points" => [[0, 0, 10]],
        "rainpoints" => [[1.5, 2.5, 20.0]],
        "inflow_bci" => nothing)
    src = FP.run_sources(fake_run(info = Dict{String,Any}("options" => o)))
    @test length(src) == 4
    @test src[1] == (lon = -0.0284, lat = 51.0, label = "inflow_q10")
    @test src[2].label == "hyd" && src[2].lat == 51.1
    @test (src[3].lon, src[3].lat) == (0.0, 0.0)
    @test (src[4].lon, src[4].lat, src[4].label) == (1.5, 2.5, "rain point")
    @test isempty(FP.run_sources(fake_run()))                      # format 1 / no record
    bci = "P\t-2.895696316\t54.90901578\tQVAR\tmain\nP -2.9371 54.8838 QVAR cummersdale\nW\t54.89\t54.91\tFREE\t\n"
    pts = FP.parse_bci_points(bci)
    @test length(pts) == 2 && pts[1].label == "main" && pts[2].lon == -2.9371
    # a .bci on disk is read through the recorded path
    p = tempname() * ".bci"; write(p, bci)
    src2 = FP.run_sources(fake_run(info = Dict{String,Any}("options" => Dict{String,Any}("inflow_bci" => p))))
    @test length(src2) == 2 && src2[1].lat == 54.90901578
    rm(p)
    # points → cells: indices as given, coordinates → nearest cell centre
    r = fake_run(n = 3, dx = 0.001, lon0 = 0.0, lat0 = 0.0)
    @test FP.point_cells(r.mesh, [5, (0.0026, 0.0004), (lon = 0.0004, lat = 0.0027, label = "x")]) == [5, 3, 7]
    @test_throws ErrorException FP.point_cells(r.mesh, [10])
end

@testset "P8 — run styles" begin
    a5li  = fake_run(scheme = "a5-v1", flow_model = "standard")
    a5swe = fake_run(scheme = "a5-v0", flow_model = "swe")
    a5sgs = fake_run(scheme = "a5-v1", flow_model = "sgs")
    cli   = fake_run(scheme = "cartesian", flow_model = "standard")
    @test FP.mesh_family(a5li) == "a5" && FP.mesh_family(cli) == "cart"
    @test FP.mesh_family(fake_run(scheme = "")) == ""
    st = FP.run_styles([a5li, a5swe, a5sgs, cli])
    @test [s.color for s in st] == [FP.SERIES[1], FP.SERIES[1], FP.SERIES[1], FP.SERIES[2]]
    @test [s.linestyle for s in st] == [:solid, :dash, :dot, :solid]
    # duplicates (e.g. two resolutions of A5 local-inertial): series colours
    st2 = FP.run_styles([a5li, a5li, cli])
    @test [s.color for s in st2] == FP.SERIES[1:3]
    @test_throws ErrorException FP.run_styles([a5li, a5li]; style = :model)
    @test FP.run_styles([cli]; style = :series)[1].color == FP.SERIES[1]
    @test_throws ErrorException FP.run_styles(fill(a5li, 9); style = :series)
end
