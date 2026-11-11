# test/test_replay.jl — visualisation/replay_data.jl (FloodA5Replay), the
# Makie-free half of the viewer's replay mode (`MakieVisualiser.replay`).
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_replay.jl
#
# Two small format-2 fixture runs (Cartesian, 25 m, flat bed):
#   A  closed boundaries; inflow-point hydrograph (CSV), injection point,
#      rain point and uniform rainfall
#   B  open boundaries; --inflow-bci (P QVAR from a .bdy found next to the
#      .bci, P QFIX) and uniform rainfall
# plus a format-1 copy of A read with its GeoParquet mesh.
# Checks: the mesh adaptor, per-frame series (domain volume, wet cells),
# the input volume rebuilt from the run record against the model's own
# `total_cumulative_input` at every frame and the recorded /summary total,
# the removed-volume handling for closed/open runs, frame reads, and that a
# mismatch with the recorded total drops the rebuilt inputs.

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
include(joinpath(@__DIR__, "..", "visualisation", "replay_data.jl"))
const FA = FloodA5Analysis
const RP = FloodA5Replay
using Test, HDF5, Logging

const RDIR = mktempdir()
const MESH_A = joinpath(RDIR, "a.parquet")
const H5_A = joinpath(RDIR, "a.h5")
const H5_B = joinpath(RDIR, "b.h5")
const SRC_IN  = (0.0012, 0.0014)     # (lon, lat) of the inflow point
const SRC_INJ = (0.0005, 0.0021)
const SRC_RP  = (0.0021, 0.0006)

_quiet(f) = with_logger(NullLogger()) do
    redirect_stdout(f, devnull)
end

let aoi = joinpath(RDIR, "aoi.geojson"), hyd = joinpath(RDIR, "hydro.csv")
    write(aoi, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
      "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.0027,0.0],[0.0027,0.0027],[0.0,0.0027],[0.0,0.0]]]}}]}""")
    write(hyd, "t_s,Q\n0,0\n300,2\n450,1\n")          # ramp: the trapezoid matters
    _quiet() do
        main(["--meshgen-cartesian", aoi, "--dx", "25", "--meshout", MESH_A,
              "--flow-model", "standard", "--closed-boundaries",
              "--inflow-point", "$(SRC_IN[2]),$(SRC_IN[1]),$hyd",
              "--injection-point", "$(SRC_INJ[2]),$(SRC_INJ[1]),0.5",
              "--rainpoint", "$(SRC_RP[2]),$(SRC_RP[1]),200",
              "--rainfall", "20",
              "--sim-duration", "600", "--output", H5_A, "--output-interval", "60"])
    end
    bdir = joinpath(RDIR, "bci"); mkpath(bdir)
    write(joinpath(bdir, "in.bci"), "P\t$(SRC_IN[1])\t$(SRC_IN[2])\tQVAR\tq1\n" *
                                    "P\t$(SRC_INJ[1])\t$(SRC_INJ[2])\tQFIX\t0.3\n")
    write(joinpath(bdir, "in.bdy"), "fixture\nq1\n3\tseconds\n0\t0\n1.5\t200\n0.5\t500\n")
    _quiet() do
        main(["--meshload", MESH_A, "--flow-model", "standard",
              "--inflow-bci", joinpath(bdir, "in.bci"), "--rainfall", "10",
              "--sim-duration", "600", "--output", H5_B, "--output-interval", "60"])
    end
end

"The model's own input integral (boundaryinputs/sources.jl) for run A at time t."
function _model_input_A(r::FA.Run, t)
    g = r.mesh
    near(p) = RP.nearest_cell(g, p[1], p[2])
    ts, q = read_timeseries(TwoColumnCSVReader(joinpath(RDIR, "hydro.csv")))["hydro"]
    srcs = AbstractSource[
        InjectionPoint(near(SRC_INJ), "", 0.5, SRC_INJ...),
        RainPoint(near(SRC_RP), "", 200 / 3_600_000.0 * g.area[near(SRC_RP)], SRC_RP..., 200.0),
        InflowPoint(near(SRC_IN), "", SRC_IN..., ts, q, "hydro")]
    return total_cumulative_input(srcs, 20 / 3_600_000.0, g.area, Float64(t))
end

@testset "R1 — BCIEntry copy matches boundary_conditions.jl" begin
    @test fieldnames(RP.BCIEntry) == fieldnames(BCIEntry)
    @test fieldtypes(RP.BCIEntry) == fieldtypes(BCIEntry)
end

@testset "R2 — mesh adaptor, resolution, adjacency" begin
    rd = RP.open_replay(H5_A)
    g = rd.run.mesh
    @test length(rd.mesh) == length(g) == length(rd.mesh.cells)
    @test [c.id for c in rd.mesh.cells] == g.ids
    @test [c.center_lon for c in rd.mesh.cells] == g.lon
    @test all(c -> c.boundary[1] == c.boundary[end], rd.mesh.cells)          # closed rings
    @test all(i -> rd.mesh.cells[i].boundary == [g.polygons[i][:, j] for j in axes(g.polygons[i], 2)],
              1:length(g))
    @test rd.mesh.resolution == 0                                           # Cartesian
    for i in 1:length(g)
        @test Set(rd.adjacency[g.ids[i]]) == Set(g.ids[j] for j in g.neighbours[:, i] if j > 0)
    end
    # A5 level from the median cell area (level r: 1/(60·4^(r-1)) of the sphere)
    for lev in (12, 18, 20)
        a = 4π * FA.R_EARTH^2 / (60 * 4.0^(lev - 1))
        m = FA.MeshGeom(["x", "y"], [0.0, 0], [0.0, 0], [0.0, 0], [a, 1.1a], Matrix{Float64}[],
                        zeros(Int32, 0, 2), Dict{String,Vector{Float64}}(), "a5-v1")
        @test RP.a5_resolution(m) == lev
    end
    a15 = 4π * FA.R_EARTH^2 / (60 * 4.0^14)
    m = FA.MeshGeom(["x"], [0.0], [0.0], [0.0], [2 * a15], Matrix{Float64}[],
                    zeros(Int32, 0, 1), Dict{String,Vector{Float64}}(), "a5-v1")
    @test RP.a5_resolution(m) == 0                                          # between levels
    # sources: the cells the model uses (nearest centre in degrees), wet at frame 1
    @test rd.source_indices == [RP.nearest_cell(g, SRC_INJ...), RP.nearest_cell(g, SRC_RP...),
                                RP.nearest_cell(g, SRC_IN...)]
    @test all(i -> FA.depth(rd.run, 1)[i] > 0, rd.source_indices)
    close(rd)
end

@testset "R3 — closed run: series and rebuilt inputs" begin
    rd = RP.open_replay(H5_A)
    r = rd.run; s = rd.series; nf = length(r.t)
    summ = r.info["summary"]
    @test s.t == r.t
    for k in 1:nf
        @test s.vol_domain[k] == FA.stored_volume(r, k)
        @test s.n_wet[k] == count(>(1e-4), FA.depth(r, k))
        @test s.vol_added[k] ≈ _model_input_A(r, r.t[k]) rtol = 1e-12
    end
    @test s.inputs_known
    @test s.vol_added[end] ≈ summ["vol_added_m3"] rtol = 1e-9
    @test s.removed === :closed && s.end_mark === nothing
    @test all(==(0.0), s.vol_removed)
    @test s.mb_err == s.vol_added .- s.vol_domain
    @test s.mb_err[end] ≈ summ["mass_balance_error_m3"] atol = 1e-6 * summ["vol_added_m3"]
    @test any(n -> occursin("match the recorded total", n), rd.notes)
    # frames: as stored, cached
    for k in (1, nf)
        fr = rd.frame(k)
        @test fr.depth == FO.read_frame(H5_A, "water_depth", k)
        @test fr.saturation == FO.read_frame(H5_A, "saturation", k)
        @test fr.volume == FO.read_frame(H5_A, "volume", k)
        @test fr.velocity == FO.read_frame(H5_A, "velocity", k)
        @test fr.vel_u == FO.read_frame(H5_A, "vel_u", k)
        @test fr.vel_v == FO.read_frame(H5_A, "vel_v", k)
        @test rd.frame(k) === fr
    end
    @test_throws ErrorException rd.frame(nf + 1)
    @test isempty(rd.missing_vars)
    close(rd)
end

@testset "R4 — open run: BCI/BDY inputs, removed volume at the end only" begin
    rd = RP.open_replay(H5_B)
    r = rd.run; s = rd.series
    summ = r.info["summary"]
    @test s.inputs_known
    @test s.vol_added[end] ≈ summ["vol_added_m3"] rtol = 1e-9
    # QFIX 0.3 m³/s + QVAR (0→1.5 m³/s over 200 s, then →0.5 at 500 s) + 10 mm/h rain
    area = sum(a for a in r.mesh.area if a >= 1.0)
    qvar(t) = t <= 200 ? 0.75 / 200 * t^2 :
              t <= 500 ? 150 + (1.5 + (1.5 - (t - 200) / 300)) / 2 * (t - 200) :
              150 + 300 + 0.5 * (t - 500)
    for k in eachindex(r.t)
        t = r.t[k]
        @test s.vol_added[k] ≈ 0.3t + qvar(t) + 10 / 3_600_000 * t * area rtol = 1e-9
    end
    @test s.removed === :end_only
    @test isequal(s.vol_removed, s.vol_added .- s.vol_domain)      # removed + MB error
    @test all(isnan, s.mb_err)
    @test s.end_mark == (r.t[end], Float64(summ["vol_removed_m3"]),
                         Float64(summ["mass_balance_error_m3"]))
    @test s.removed_label == "Removed + MB err"
    @test length(rd.source_indices) == 2
    close(rd)
end

@testset "R5 — a recorded total that does not match drops the rebuilt inputs" begin
    p = joinpath(RDIR, "tampered.h5")
    cp(H5_A, p; force = true)
    h5open(p, "r+") do f
        a = attributes(f["summary"])
        v = read(a["vol_added_m3"])
        delete_attribute(f["summary"], "vol_added_m3")
        a["vol_added_m3"] = 1.01 * v
    end
    rd = RP.open_replay(p)
    @test !rd.series.inputs_known
    @test all(isnan, rd.series.vol_added) && all(isnan, rd.series.mb_err)
    @test any(n -> occursin("do not match", n), rd.notes)
    close(rd)
end

@testset "R6 — format-1 copy with its GeoParquet mesh" begin
    p = joinpath(RDIR, "v1.h5")
    r2 = FA.load_run(H5_A)
    h5open(p, "w") do f
        g = create_group(f, "mesh")
        g["cell_ids"] = r2.mesh.ids; g["center_lons"] = r2.mesh.lon
        g["center_lats"] = r2.mesh.lat; g["elevations"] = r2.mesh.elevation
        fr = create_group(f, "frames")
        for k in eachindex(r2.t)
            gk = create_group(fr, lpad(string(k), 6, '0'))
            gk["t"] = r2.t[k]
            for v in ("water_depth", "volume", "saturation", "velocity")   # no vel_u/vel_v
                gk[v] = FO.read_frame(H5_A, v, k)
            end
        end
    end
    @test_throws ErrorException RP.open_replay(p)          # no polygons without the mesh
    rd = RP.open_replay(p; parquet = MESH_A)
    rd2 = RP.open_replay(H5_A)
    @test [c.boundary for c in rd.mesh.cells] == [c.boundary for c in rd2.mesh.cells]
    @test rd.series.vol_domain == rd2.series.vol_domain
    @test rd.series.n_wet == rd2.series.n_wet
    @test all(isnan, rd.series.vol_added) && !rd.series.inputs_known
    @test rd.series.removed === :unknown
    @test sort(rd.missing_vars) == ["vel_u", "vel_v"]
    k = length(rd.t)
    @test rd.frame(k).depth == rd2.frame(k).depth
    @test all(==(0.0), rd.frame(k).vel_u)
    @test any(n -> occursin("format-1", n), rd.notes)
    @test isempty(rd.source_indices)
    close(rd); close(rd2)
end

@testset "R7 — --replay needs the viewer (stubbed here)" begin
    @test_throws ErrorException main(["--replay", H5_A])
end
