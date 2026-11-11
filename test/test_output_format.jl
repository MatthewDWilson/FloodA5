# test/test_output_format.jl — HDF5 output format 2 (0.4.0): run record,
# embedded mesh, results layout, --rerun, and the format-1 reader.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_output_format.jl
#
# O1  a CLI run writes format 2: results layout, run record (argv, options,
#     git, environment), embedded mesh and inputs, geometry, summary
# O2  --rerun reproduces every frame bit for bit from the file alone (the
#     original mesh and input files moved away)
# O3  a mesh generated in the run and not saved is embedded (re-saved copy),
#     and --rerun reproduces that run too
# O4  --embed-mesh geometry: no mesh bytes; --rerun uses the original file and
#     refuses it once it has changed
# O5  the reader still reads format-1 files (FloodA5 ≤ 0.3.x)

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
using Test, HDF5, Logging

const DIR = mktempdir()
const AOI = joinpath(DIR, "aoi.geojson")
write(AOI, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
  "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.0027,0.0],[0.0027,0.0027],[0.0,0.0027],[0.0,0.0]]]}}]}""")
const HYD = joinpath(DIR, "hydro.csv")
write(HYD, "t_s,Q\n0,0\n300,2\n600,1\n")
const MESH = joinpath(DIR, "m.parquet")

_run(args) = with_logger(NullLogger()) do
    redirect_stdout(devnull) do
        main(copy(args))
    end
end
_same_results(a, b) = FO.frame_times(a) == FO.frame_times(b) &&
    FO.frame_vars(a) == FO.frame_vars(b) &&
    all(v -> isequal(FO.read_var(a, v), FO.read_var(b, v)), FO.frame_vars(a))

const O1 = joinpath(DIR, "o1.h5")
const ARGS1 = ["--meshgen-cartesian", AOI, "--dx", "25", "--meshout", MESH,
               "--flow-model", "standard", "--rainfall", "40",
               "--inflow-point", "0.00135,0.00135,$HYD",
               "--sim-duration", "600", "--output", O1, "--output-interval", "60"]

@testset "O1 — format 2 written by a CLI run" begin
    _run(ARGS1)
    @test FO.format_version(O1) == 2
    nf = FO.n_frames(O1)
    @test nf >= 5
    @test FO.frame_times(O1)[end] == 600.0
    @test issubset(["water_depth", "volume", "saturation", "velocity", "vel_u", "vel_v"],
                   FO.frame_vars(O1))
    ct = FO.cell_table(O1); n = length(ct.ids)
    @test ct.elevation == zeros(n)                     # no DEM: recorded as the z = 0 the model used
    @test size(FO.read_var(O1, "volume")) == (n, nf)
    @test FO.read_frame(O1, "volume") == FO.read_var(O1, "volume")[:, end]
    @test FO.read_series(O1, "volume", 3) == FO.read_var(O1, "volume")[3, :]
    @test sum(FO.read_frame(O1, "volume")) > 0
    ri = FO.run_info(O1)
    @test ri["flooda5_version"] == FLOODA5_VERSION
    @test ri["argv"] == ARGS1
    @test ri["options"]["scheme"] == "recommended" && ri["options"]["eta_filter"] == 0.05
    @test ri["flow_model"] == "standard" && ri["threads"] == Threads.nthreads()
    @test !isempty(ri["git_commit"]) && !isempty(ri["project"])
    @test ri["summary"]["completed"] == true && ri["summary"]["n_frames"] == nf
    s = ri["summary"]                                  # the run's mass balance is recorded
    @test s["vol_added_m3"] > 0 && s["vol_domain_m3"] > 0
    @test s["vol_domain_m3"] ≈ sum(filter(isfinite, FO.read_frame(O1, "volume"))) rtol = 1e-12
    # conservation: what the stepping added is stored or left. mass_balance_error
    # also contains the sources' first-order time integration (vol_applied −
    # vol_added: Q(t)·dt against the hydrograph's exact integral)
    @test abs(s["conservation_error_m3"]) < 1e-9 * s["vol_applied_m3"] + 1e-9
    @test s["mass_balance_error_m3"] ≈ s["conservation_error_m3"] + s["vol_added_m3"] - s["vol_applied_m3"] atol = 1e-6
    ms = ri["mesh_source"]
    @test ms["embedded"] == "full" && ms["kind"] == "file" && ms["sha256"] == FO.sha256_file(MESH)
    @test length(ri["inputs"]) == 1 && ri["inputs"][1]["embedded"] == true &&
          ri["inputs"][1]["sha256"] == FO.sha256_file(HYD)
    # geometry: one polygon per cell; neighbour lists symmetric; areas ≈ dx²
    polys = FO.cell_polygons(O1)
    @test length(polys) == n && all(p -> size(p, 1) == 2 && size(p, 2) >= 4, polys)
    nb = FO.cell_neighbours(O1)
    @test all(1:n) do i
        all(j -> j == 0 || i in nb[:, j], nb[:, i])
    end
    m = isqrt(n)                                        # an m × m grid: 2m(m−1) faces
    @test m^2 == n && count(>(0), nb) == 2 * (2 * m * (m - 1))
    @test all(a -> isapprox(a, 625.0; rtol = 0.02), FO.cell_area(O1))
end

@testset "O2 — --rerun reproduces the run from the file alone" begin
    mv(MESH, MESH * ".away"); mv(HYD, HYD * ".away")
    O2 = joinpath(DIR, "o2.h5")
    _run(["--rerun", O1, "--output", O2])
    @test _same_results(O1, O2)
    ri = FO.run_info(O2)
    @test "--meshload" in ri["argv"] && !("--meshgen-cartesian" in ri["argv"])
    @test ri["mesh_source"]["sha256"] == FO.run_info(O1)["mesh_source"]["sha256"]
    # a rerun's own file is rerunnable, and an override replaces a recorded flag
    O3 = joinpath(DIR, "o3.h5")
    _run(["--rerun", O2, "--output", O3, "--sim-duration", "300"])
    @test FO.frame_times(O3)[end] == 300.0
    k = FO.n_frames(O3) - 1                 # frames before the final one (t = 300) agree
    @test FO.frame_times(O3)[1:k] == FO.frame_times(O1)[1:k]
    @test FO.read_var(O3, "volume")[:, 1:k] == FO.read_var(O1, "volume")[:, 1:k]
    mv(MESH * ".away", MESH); mv(HYD * ".away", HYD)
end

@testset "O3 — mesh generated in the run, not saved" begin
    O4 = joinpath(DIR, "o4.h5")
    _run(["--meshgen-cartesian", AOI, "--dx", "25", "--flow-model", "standard",
          "--rainfall", "40", "--sim-duration", "300", "--output", O4])
    ms = FO.run_info(O4)["mesh_source"]
    @test ms["kind"] == "resaved" && ms["embedded"] == "full"
    O5 = joinpath(DIR, "o5.h5")
    _run(["--rerun", O4, "--output", O5])
    @test _same_results(O4, O5)
end

@testset "O4 — --embed-mesh geometry" begin
    O6 = joinpath(DIR, "o6.h5")
    _run(["--meshload", MESH, "--flow-model", "standard", "--rainfall", "40",
          "--sim-duration", "300", "--output", O6, "--embed-mesh", "geometry"])
    ms = FO.run_info(O6)["mesh_source"]
    @test ms["embedded"] == "geometry" && ms["sha256"] == FO.sha256_file(MESH)
    @test h5open(f -> !haskey(f["mesh/source"], "parquet"), O6)
    @test length(FO.cell_polygons(O6)) == length(FO.cell_table(O6).ids)
    @test filesize(O6) < filesize(O1)
    O7 = joinpath(DIR, "o7.h5")
    _run(["--rerun", O6, "--output", O7])
    @test _same_results(O6, O7)
    open(MESH, "a") do io; write(io, UInt8[0]); end           # the mesh file changes
    @test_throws ErrorException _run(["--rerun", O6, "--output", joinpath(DIR, "o8.h5")])
end

@testset "O5 — format-1 files are still readable" begin
    p = joinpath(DIR, "v1.h5")
    h5open(p, "w") do f
        g = create_group(f, "mesh")
        g["cell_ids"] = ["a", "b", "c"]; g["center_lons"] = [0.0, 1.0, 2.0]
        g["center_lats"] = zeros(3); g["elevations"] = [1.0, 2.0, 3.0]
        fr = create_group(f, "frames")
        for k in 1:12                       # 12 frames: names sort numerically, not as text
            gk = create_group(fr, lpad(string(k), 6, '0'))
            gk["t"] = 10.0k; gk["water_depth"] = [k, 2k, 3k] .* 1.0
        end
    end
    @test FO.format_version(p) == 1 && FO.n_frames(p) == 12
    @test FO.frame_times(p) == 10.0 .* (1:12)
    @test FO.read_frame(p, "water_depth", 11) == [11.0, 22.0, 33.0]
    @test FO.read_series(p, "water_depth", 2) == 2.0 .* (1:12)
    @test size(FO.read_var(p, "water_depth")) == (3, 12)
    @test FO.frame_at(p, 52.0) == 5
    @test FO.cell_table(p).elevation == [1.0, 2.0, 3.0]
    @test_throws ErrorException FO.cell_polygons(p)
    @test FO.run_info(p)["format_version"] == 1
end

@testset "O6 — run_simulation! called directly creates its own output file" begin
    d = mktempdir()
    FO.extract_mesh(O1, d)
    m = with_logger(() -> A5Grid.load_mesh_geoparquet(joinpath(d, "mesh.parquet")), NullLogger())
    p = joinpath(DIR, "o9.h5")
    cp(joinpath(DIR, "v1.h5"), p)        # a stale format-1 file at the path is replaced …
    nf = Int[]
    for pass in 1:2                      # … and so is an earlier run's format-2 file
        s = with_logger(() -> initialise_flow_model(m, StandardFlow()), NullLogger())
        out = SimOutput(p, 60.0, -Inf, 0, true)
        _run_direct() = run_simulation!(s, m, 180.0, 10.0, nothing, :none, out;
                                        method = StandardFlow(), rainfall_rate = 40e-3 / 3600)
        with_logger(NullLogger()) do
            redirect_stdout(_run_direct, devnull)
        end
        @test FO.format_version(p) == 2
        @test FO.frame_times(p)[end] == 180.0 && issorted(FO.frame_times(p))
        @test length(FO.cell_table(p).ids) == length(m.cells)
        @test length(FO.cell_area(p)) == length(m.cells)
        push!(nf, FO.n_frames(p))
    end
    @test nf[1] == nf[2] >= 3
end

@testset "O7 — default flow model (standard since 0.5.0); --rerun pins the recorded one" begin
    d = mktempdir()
    mesh = FO.extract_mesh(O1, d)
    p = joinpath(DIR, "o10.h5")
    _run(["--meshload", mesh, "--rainfall", "40", "--sim-duration", "120",
          "--output", p, "--output-interval", "60"])                # no --flow-model
    ri = FO.run_info(p)
    @test ri["flow_model"] == "standard" && !("--flow-model" in ri["argv"])
    rv = FO.rerun_argv(p; output = joinpath(d, "r.h5"), workdir = d)
    i = findfirst(==("--flow-model"), rv)
    @test i !== nothing && rv[i+1] == "standard"                  # pinned
    j = findfirst(==("--scheme"), rv)
    @test j !== nothing && rv[j+1] == string(ri["options"]["scheme"])
    # a run that took the pre-0.5.0 default (sgs) reruns as sgs, not as today's default
    q = joinpath(DIR, "o10_sgs.h5"); cp(p, q)
    h5open(q, "r+") do f
        HDF5.delete_attribute(f["run"], "flow_model")
        attributes(f["run"])["flow_model"] = "sgs"
    end
    rv = FO.rerun_argv(q; output = joinpath(d, "r2.h5"), workdir = d)
    @test rv[findfirst(==("--flow-model"), rv) + 1] == "sgs"
    # an explicit flag is kept, not duplicated; an override wins
    rv = FO.rerun_argv(O1; output = joinpath(d, "r3.h5"), workdir = d)
    @test count(==("--flow-model"), rv) == 1
    rv = FO.rerun_argv(q; output = joinpath(d, "r4.h5"), workdir = d, overrides = ["--flow-model", "swe"])
    @test count(==("--flow-model"), rv) == 1 && rv[findfirst(==("--flow-model"), rv) + 1] == "swe"
end

println("\nOutput-format tests complete.")
