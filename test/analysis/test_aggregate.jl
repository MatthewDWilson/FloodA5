# A fine run on a coarse mesh's cells (analysis/aggregate.jl): assignment,
# aggregation, SGS wetted area, and the like-for-like comparison.

"A (lon, lat) square polygon of side s (m) with its south-west corner at (x, y) m from (0°, 0°)."
function _agg_square(x, y, s)
    d = 180 / (π * FA.R_EARTH)                      # degrees per metre at the equator
    xs = [x, x + s, x + s, x, x] .* d; ys = [y, y, y + s, y + s, y] .* d
    return permutedims(hcat(xs, ys))
end

@testset "A-aggregate — 10 × 10 fine cells on a 2 × 2 coarse grid" begin
    d = 180 / (π * FA.R_EARTH)
    # coarse: 50 m squares, SW, SE, NW, NE
    coarse = [_agg_square(x, y, 50.0) for (x, y) in ((0, 0), (50, 0), (0, 50), (50, 50))]
    # fine: 10 m cells, centres at 5, 15, …, 95 m; depth 0.05·i along x (i = 1…10), bed 0
    I = [i for i in 1:10, j in 1:10]; J = [j for i in 1:10, j in 1:10]
    lon = vec((10.0 .* I .- 5.0) .* d); lat = vec((10.0 .* J .- 5.0) .* d)
    h = vec(0.05 .* I); η = copy(h); A = fill(100.0, 100); V = h .* A
    # three cells of the NE coarse cell take no part (no terrain): it is not covered
    gone = [k for k in 1:100 if I[k] == 10 && J[k] >= 8]
    V[gone] .= NaN; h[gone] .= NaN; η[gone] .= NaN
    # one more point outside every polygon
    own = FA.assign_to_cells(vcat(lon, 150d), vcat(lat, 5d), coarse; bucket = 20.0)
    @test own[end] == 0
    own = own[1:end-1]
    @test [count(==(k), own) for k in 1:4] == [25, 25, 25, 25]
    @test all(own[k] == (I[k] <= 5 ? 1 : 2) + (J[k] <= 5 ? 0 : 2) for k in 1:100)
    # the same with the default projection centre and bucket size
    @test FA.assign_to_cells(lon, lat, coarse) == own

    a = FA.aggregate_fine(own, fill(2500.0, 4), V, h, η, A; wet = 0.1, cover = 0.9)
    @test a.volume ≈ [375.0, 1000.0, 375.0, 850.0]
    @test a.fine_area == [2500.0, 2500.0, 2500.0, 2200.0]
    @test a.wet_area == [1500.0, 2500.0, 1500.0, 2200.0]     # depth > 0.1: i ≥ 3
    @test a.wet_frac ≈ [0.6, 1.0, 0.6, 1.0]
    @test a.level[1:3] ≈ [0.2, 0.4, 0.2]
    @test a.level[4] ≈ (5 * (0.3 + 0.35 + 0.4 + 0.45 + 0.5) - 3 * 0.5) / 22
    @test a.covered == [true, true, true, false]
    # a dry coarse cell: wet fraction 0, no level
    a0 = FA.aggregate_fine(own, fill(2500.0, 4), V, zeros(100), η, A)
    @test all(iszero, a0.wet_frac[1:3]) && all(isnan, a0.level)

    # a cell-mean model on the coarse cells
    Vm = [500.0, 1000.0, 0.0, 900.0]; Am = fill(2500.0, 4)
    ηm = [0.2, 0.45, 0.0, 0.36]
    fm = FA.wet_fraction_cells(Vm ./ Am; wet = 0.1)
    @test fm == [1.0, 1.0, 0.0, 1.0]
    c = FA.compare_aggregated(a, Vm, ηm, fm, Am)
    @test c.storage == 1500.0 && c.storage_fine ≈ 1750.0 && c.storage_ratio ≈ 1500 / 1750   # NE left out
    @test c.wet_area == 5000.0 && c.wet_area_fine ≈ 5500.0
    @test c.csi ≈ 1.6 / 2.6
    @test c.n_level == 2 && c.level_median ≈ 0.025 && c.level_mean ≈ 0.025
    @test c.level_rms ≈ sqrt(0.05^2 / 2)
    @test c.n_depth == 3 && c.depth_rms ≈ sqrt((0.05^2 + 0.15^2) / 3)
    @test_throws ErrorException FA.compare_aggregated(a, Vm[1:3], ηm, fm, Am)
end

@testset "A-aggregate — SGS wetted area" begin
    t = FA.SGSAreaTable([0.0, 1.0, 2.0], [0.0, 50.0, 100.0], 100.0, 0.0, 2.0)
    @test FA.sgs_wet_area(t, -1.0) == 0.0 && FA.sgs_wet_area(t, 0.0) == 0.0
    @test FA.sgs_wet_area(t, 0.5) ≈ 25.0 && FA.sgs_wet_area(t, 1.5) ≈ 75.0
    @test FA.sgs_wet_area(t, 3.0) == 100.0
    tn = FA.SGSAreaTable(Float64[], Float64[], NaN, NaN, NaN)          # a cell without terrain
    @test FA.wet_fraction_sgs([t, t, tn], [1.6, NaN, 1.0]; wet = 0.1) ≈ [0.75, 0.0, 0.0]
    @test FA.wet_fraction_cells([0.05, 0.2, NaN]) == [0.0, 1.0, 0.0]
end

@testset "A-aggregate — SGS tables from a mesh file, and the Run wrappers" begin
    # an SGS mesh of the fixture's AOI over a plane with a 4 m channel, and a short SGS run on it
    dem = joinpath(FIX_DIR, "agg_dem.tif"); sgsmesh = joinpath(FIX_DIR, "agg_sgs.parquet")
    sgsh5 = joinpath(FIX_DIR, "agg_sgs.h5"); aoi = joinpath(FIX_DIR, "aoi.geojson")
    hyd = joinpath(FIX_DIR, "hydro.csv")
    px = rad2deg(2.0 / FA.R_EARTH); nx = ny = 160
    z = Float32[0.001 * 2i + (abs(2j - 150) <= 2 ? -1.0 : 0.0) for i in 0:nx-1, j in 0:ny-1]
    ArchGDAL.create(dem; driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                    nbands = 1, dtype = Float32) do ds
        ArchGDAL.setgeotransform!(ds, [-0.0001, px, 0.0, 0.0028, 0.0, -px])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
        ArchGDAL.write!(ds, z, 1)
    end
    with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            main(["--meshgen-cartesian", aoi, "--dx", "25", "--dem", dem, "--meshout", sgsmesh,
                  "--flow-model", "sgs", "--mesh-only"])
            main(["--meshload", sgsmesh, "--flow-model", "sgs", "--closed-boundaries",
                  "--inflow-point", "$(FIX_SRC[2]),$(FIX_SRC[1]),$hyd",
                  "--sim-duration", "300", "--output", sgsh5, "--output-interval", "60"])
        end
    end
    # the tables as FA reads them (file row order) against the model's own lookup
    tabs = FA.load_sgs_tables(sgsmesh)
    m = with_logger(() -> load_mesh_geoparquet(sgsmesh), NullLogger())
    @test length(tabs) == length(m.cells)
    n_checked = 0
    for i in eachindex(tabs)
        isnan(tabs[i].z_min) && continue
        st = A5Grid.sgs_table(m, i)
        @test tabs[i].z_min == st.z_min && tabs[i].z_max == st.z_max && tabs[i].cell_area == st.cell_area
        for η in range(st.z_min - 0.2, st.z_max + 0.2; length = 23)
            FA.sgs_wet_area(tabs[i], η) == A5Grid.wetted_area_from_wse(st, η) || (@test false; break)
        end
        n_checked += 1
    end
    @test n_checked > 50
    # any list or scalar column, e.g. the per-face sill raise of connectivity, in row or id order
    cols = FA.load_mesh_columns(sgsmesh, ["sgs_face_raise", "sgs_z_min", "no_such_column"])
    @test !haskey(cols, "no_such_column")
    @test all(cols["sgs_face_raise"][i] == m.array_vars["sgs_face_raise"][:, i] for i in eachindex(m.cells))
    @test isequal(cols["sgs_z_min"], m.static_vars["sgs_z_min"])
    pq = FA.load_mesh_parquet(sgsmesh)
    rev = reverse(pq.ids)
    @test isequal(FA.load_mesh_columns(sgsmesh, ["sgs_z_min"]; ids = rev)["sgs_z_min"], reverse(cols["sgs_z_min"]))
    @test_throws ErrorException FA.load_mesh_columns(sgsmesh, ["sgs_z_min"]; ids = ["ffffffffffffffff"])
    # Run level: the tables in the output's cell order, checked against its sgs_z_min
    r = FA.load_run(sgsh5)
    @test r.flow_model == "sgs"
    t_run = FA.sgs_tables(r)
    @test [t.z_min for t in t_run] == r.mesh.static["sgs_z_min"]
    k = FA.nframes(r)
    f = FA.wetted_fraction(r, k)
    @test all(0 .<= f .<= 1) && any(>(0), f)
    @test f == FA.wet_fraction_sgs(t_run, FA.wse(r, k))
    @test_throws ErrorException FA.sgs_tables(r; parquet = FIX_MESH)      # not this mesh (no tables)
    # a run against itself, aggregated onto its own cells: identical
    s = FA.load_run(FIX_H5)
    c = FA.compare_to_fine(s, s, s.t[end]; wet = 0.01)       # the fixture is shallow: wet from 1 cm
    @test c.storage_ratio ≈ 1.0 && c.csi ≈ 1.0 && c.wet_area > 0
    @test c.storage ≈ FA.stored_volume(s, FA.nframes(s))
    @test c.n_level > 0 && abs(c.level_median) < 1e-12 && c.level_rms < 1e-12
    @test c.owner == collect(1:length(s.mesh))
    # the SGS run against the standard fixture (two 25 m grids of the same AOI): it runs end to end
    c2 = FA.compare_to_fine(s, r, 300.0; tables = t_run, tol = 60.0)
    @test c2.storage_fine <= FA.stored_volume(s, FA.frame_at(s, c2.t_fine)) * (1 + 1e-9)
    @test isfinite(c2.storage_ratio) && 0 <= c2.csi <= 1 && abs(c2.t_fine - 300.0) <= 60
end
