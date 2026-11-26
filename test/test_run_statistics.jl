# test/test_run_statistics.jl — statistics over the run written by the model
# (0.8.13: /statistics, --stats-interval, --flood-depth)
#
#   julia --project=ci test/test_run_statistics.jl
#
# - At the default cadence (the output frames) the model's statistics equal,
#   bit for bit, those derived from the frames written (FA.run_statistics, the
#   viewers' implementation): both go through the same accumulator with the
#   same Float32 values.
# - --stats-interval 0 (every step) changes no frame, counts every step, and
#   its maxima are at least those of the frames; --flood-depth sets the
#   threshold.
# - The accumulator added frame by frame equals run_statistics, and can be
#   read part-way.
# - Files without statistics (earlier versions) read as `nothing`.

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
isdefined(Main, :FloodA5Analysis) || include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
using Test, HDF5

const SDIR = mktempdir()
const AOI = joinpath(SDIR, "aoi.geojson")
write(AOI, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
  "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.003,0.0],[0.003,0.002],[0.0,0.002],[0.0,0.0]]]}}]}""")

function run_case(h5, extra = String[])
    args = ["--meshgen-cartesian", AOI, "--dx", "25", "--flow-model", "standard",
            "--injection-point", "0.001,0.0015,0.5", "--rainfall", "20",
            "--sim-duration", "600", "--output", h5, "--output-interval", "120", extra...]
    redirect_stdout(() -> main(args), devnull)
    return h5
end

# equal as stored: NaN = NaN, otherwise ===
same(a, b) = length(a) == length(b) && all(i -> isequal(a[i], b[i]), eachindex(a))

@testset "statistics over the run (model)" begin
    h1 = run_case(joinpath(SDIR, "default.h5"))
    h2 = run_case(joinpath(SDIR, "every_step.h5"), ["--stats-interval", "0", "--flood-depth", "0.05"])

    t = FO.frame_times(h1)
    D = FO.read_var(h1, "water_depth"); V = FO.read_var(h1, "velocity")
    n = size(D, 1)

    @testset "default cadence: equal to the frames' statistics" begin
        r = FO.read_statistics(h1)
        @test r !== nothing
        st, at = r
        @test at["cadence"] == "output" && at["interval_s"] == 120.0 && at["flood_depth"] == 0.1
        @test at["n_states"] == length(t) && at["t_last"] == t[end]
        ref = FA.run_statistics(t, k -> (Float32.(D[:, k]), Float32.(V[:, k])), n)
        @test Set(keys(st)) == Set(["max_water_depth", "time_max_water_depth", "time_first_flooded",
                                    "flood_duration", "max_velocity", "max_depth_velocity", "max_hazard_class_arr"])
        for (k, v) in st
            ok = same(v, ref[k])
            ok || @warn "statistic $k differs from the frames'"
            @test ok
        end
        @test any(>(0), st["max_water_depth"]) && any(isfinite, st["time_max_water_depth"])
        @test FO.run_info(h1)["statistics"]["cadence"] == "output"
    end

    @testset "--stats-interval 0, --flood-depth" begin
        # the statistics change no result
        @test FO.frame_times(h2) == t
        for v in ("water_depth", "volume", "velocity")
            @test FO.read_var(h2, v) == FO.read_var(h1, v)
        end
        st, at = FO.read_statistics(h2)
        @test at["cadence"] == "step" && at["interval_s"] == 0.0 && at["flood_depth"] == 0.05
        @test at["n_states"] > 2 * length(t)                 # every step (plus the start)
        ref = FA.run_statistics(t, k -> (Float32.(D[:, k]), Float32.(V[:, k])), n)
        @test all(st["max_water_depth"] .>= ref["max_water_depth"])
        @test all(st["max_velocity"] .>= ref["max_velocity"])
        tm = filter(isfinite, st["time_max_water_depth"])
        @test !isempty(tm) && all(0 .<= tm .<= t[end])
        # a lower threshold: flooded earlier and longer
        ff = st["time_first_flooded"]; fr = ref["time_first_flooded"]
        @test count(isfinite, ff) >= count(isfinite, fr)
        @test sum(st["flood_duration"]) >= sum(ref["flood_duration"])
        @test "--stats-interval" in FO.run_info(h2)["argv"]    # (so --rerun repeats it)
        @test FO.run_info(h2)["options"]["stats_interval"] == 0.0
    end

    @testset "device backend (ka-cpu): the statistics' syncs change no frame" begin
        k1 = run_case(joinpath(SDIR, "ka_default.h5"), ["--backend", "ka-cpu"])
        k2 = run_case(joinpath(SDIR, "ka_step.h5"), ["--backend", "ka-cpu", "--stats-interval", "0"])
        for v in ("water_depth", "volume", "velocity")
            @test FO.read_var(k2, v) == FO.read_var(k1, v)
        end
        s1, _ = FO.read_statistics(k1); s2, a2 = FO.read_statistics(k2)
        @test a2["n_states"] > 2 * length(FO.frame_times(k2))
        @test all(s2["max_water_depth"] .>= s1["max_water_depth"])
    end

    @testset "accumulator: incremental = one pass" begin
        p = FA.HazardParams(flood_depth = 0.02)
        ref = FA.run_statistics(t, k -> (D[:, k], V[:, k]), n, p)
        a = FA.StatsAccumulator(n, p)
        part = nothing
        for k in eachindex(t)
            FA.accumulate_frame!(a, t[k], D[:, k], V[:, k])
            k == 2 && (part = FA.accumulated_statistics(a))
        end
        full = FA.accumulated_statistics(a)
        @test Set(keys(full)) == Set(keys(ref))
        @test all(k -> same(full[k], ref[k]), keys(ref))
        # part-way: the statistics of the first two frames
        ref2 = FA.run_statistics(t[1:2], k -> (D[:, k], V[:, k]), n, p)
        @test all(k -> same(part[k], ref2[k]), keys(ref2))
        # without the parameter-dependent maxima
        b = FA.StatsAccumulator(n, p; hazards = false)
        FA.accumulate_frame!(b, 0.0, D[:, 1], V[:, 1])
        @test !haskey(FA.accumulated_statistics(b), "max_hazard_rating")
        @test haskey(FA.accumulated_statistics(b), "max_hazard_class_arr")
    end

    @testset "FA.write_run_statistics!: an older file gets the model's statistics" begin
        p = joinpath(SDIR, "older.h5")
        cp(h1, p)
        h5open(f -> delete_object(f, "statistics"), p, "r+")
        @test FO.read_statistics(p) === nothing
        st = FA.write_run_statistics!(p)
        @test st !== nothing
        s0, _ = FO.read_statistics(h1)
        s1, a1 = FO.read_statistics(p)
        @test Set(keys(s1)) == Set(keys(s0)) && all(k -> same(s1[k], s0[k]), keys(s0))
        @test a1["added_by"] == "FA.write_run_statistics!" && a1["n_states"] == length(t)
        @test FA.write_run_statistics!(p) === nothing              # already there: left alone
    end

    @testset "files without statistics" begin
        p = joinpath(SDIR, "old.h5")
        FO.create_output(p; cells = (ids = ["a", "b"], center_lons = [0.0, 1.0], center_lats = [0.0, 1.0],
                                     elevations = [1.0, 2.0], static_vars = Dict{String,Vector{Float64}}(),
                                     id_scheme = "test"))
        @test FO.read_statistics(p) === nothing
        @test !haskey(FO.run_info(p), "statistics")
    end
end
