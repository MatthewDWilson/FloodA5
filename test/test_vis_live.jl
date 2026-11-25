# test/test_vis_live.jl — the web globe viewer's server end to end
#
#   julia --project=ci -t 2 test/test_vis_live.jl
#
# A small Cartesian run with `--vis cesium` (and an HDF5 output) while a
# WebSocket client follows it: the client sees the mesh, live frames and the
# end of the run; the last frame the server holds is the run's end state, and
# cell values (looked up through /api/cell, so in the server's cell order)
# equal the output file's. Then the same file is served as a replay
# (visualisation/web_replay.jl): frames and cell order as in the file.
# Also: a bare `--vis` now means the Makie viewer.

ENV["FLOODA5_NO_VIS"] = "1"
ENV["FLOODA5_VIS_INTERVAL"] = "0"     # keep every live frame (within the budget)
ENV["FLOODA5_VIS_LINGER"] = "3"       # serve 3 s after the run, then stop
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
include(joinpath(@__DIR__, "..", "visualisation", "web_replay.jl"))
const FA = FloodA5Analysis
using Test, HTTP, JSON3, HDF5, Logging

const LDIR = mktempdir()
const AOI = joinpath(LDIR, "aoi.geojson")
const H5 = joinpath(LDIR, "live.h5")
write(AOI, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
  "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.003,0.0],[0.003,0.002],[0.0,0.002],[0.0,0.0]]]}}]}""")

free_port() = (s = HTTP.Sockets.listen(HTTP.Sockets.IPv4("127.0.0.1"), 0);
               p = Int(HTTP.Sockets.getsockname(s)[2]); close(s); p)
get_(port, path) = HTTP.get("http://127.0.0.1:$port$path"; status_exception = false, retry = false)
json_(port, path) = JSON3.read(String(get_(port, path).body))

@testset "web viewer, live and replay" begin
    port = free_port()
    args = ["--meshgen-cartesian", AOI, "--dx", "25", "--flow-model", "standard",
            "--injection-point", "0.001,0.0015,0.5", "--rainfall", "20",
            "--sim-duration", "300", "--output", H5, "--output-interval", "60",
            "--vis", "cesium", "--vis-port", string(port)]
    args0 = copy(args)                      # (main may consume its argument vector)
    run_task = Threads.@spawn with_logger(NullLogger()) do
        redirect_stdout(() -> main(args), devnull)
    end

    msgs = Any[]; mesh_bytes = Ref{Vector{UInt8}}(UInt8[]); last_cells = Ref{Any}(nothing)
    live_run = Ref{Any}(nothing)
    client = Threads.@spawn begin
        # wait for the server
        t0 = time()
        while time() - t0 < 600
            istaskdone(run_task) && break
            ok = try get_(port, "/api/info").status == 200 catch; false end
            ok && break
            sleep(0.05)
        end
        try
            HTTP.WebSockets.open("ws://127.0.0.1:$port/api/live") do ws
                for m in ws
                    msg = JSON3.read(m); push!(msgs, msg)
                    if (msg.type == "complete") || (msg.type == "hello" && msg.complete)
                        r = get_(port, "/api/mesh"); r.status == 200 && (mesh_bytes[] = r.body)
                        fl = json_(port, "/api/frames")
                        last = fl.frames[end].id
                        last_cells[] = (frames = fl, cells = [json_(port, "/api/cell/$i?frame=$last") for i in 0:5:59])
                        live_run[] = json_(port, "/api/run")
                    end
                end
            end
        catch e
            e isa HTTP.WebSockets.WebSocketError || e isa EOFError || e isa Base.IOError || rethrow()
        end
    end
    fetch(run_task)
    wait(client)

    types = [m.type for m in msgs]
    @test !isempty(types) && types[1] == "hello"
    @test "complete" in types || msgs[1].complete
    nframe = count(==("frame"), types)
    @test nframe >= 2 || msgs[1].frames >= 2
    @test length(mesh_bytes[]) > 32 && mesh_bytes[][1:4] == b"FA5M"
    lc = last_cells[]
    @test lc !== nothing
    if lc !== nothing
        @test lc.frames.complete == true
        @test lc.frames.frames[end].t ≈ 300.0          # the end state is always pushed
        # values per cell id equal the output file's last frame
        run = FA.load_run(H5)
        h = h5open(f -> FloodA5Output.read_frame(f, "water_depth", length(run.t)), H5, "r")
        pos = Dict(id => k for (k, id) in enumerate(run.mesh.ids))
        @test length(lc.cells) == 12
        for c in lc.cells
            @test Float32(c.values.water_depth) == Float32(h[pos[c.id]])
        end
        @test any(c -> c.values.water_depth > 0, lc.cells)
    end
    # the live run's record (the viewer's Run tab), as the output file will have it
    lr = live_run[]
    @test lr !== nothing
    if lr !== nothing
        @test lr.flooda5_version == FLOODA5_VERSION && lr.flow_model == "standard"
        @test collect(String, lr.argv) == args0 && lr.options.manning_n isa Real
        @test lr.mesh.n_cells == length(FA.load_run(H5).mesh) && lr.mesh.id_scheme == "cartesian" && lr.mesh.cell_area_mean_m2 > 0
        @test lr.output == abspath(H5) && !haskey(lr, :project) && !haskey(lr, :manifest)
    end
    @test VisualisationServer._port_free("127.0.0.1", port)   # stopped after the linger

    # the same run as a replay
    rport = free_port()
    vis, f = web_replay_server(H5; port = rport, viewer = "giro3d", autoplay = true)
    try
        run = FA.load_run(H5)
        cfg = json_(rport, "/api/config")
        @test cfg.autoplay == true && !haskey(cfg, :vertical_offset)    # no --offset: estimated
        rr = json_(rport, "/api/ranges?wet=0")         # (a short test run: little is 5 cm deep)
        dr = rr.vars.water_depth
        @test rr.complete && rr.frames == min(length(run.t), 24) && dr.n > 0
        @test 0 <= dr.min <= dr.p02 <= dr.p98 <= dr.max
        info = json_(rport, "/api/info")
        @test info.mode == "replay" && info.complete && info.n_cells == length(run.mesh)
        @test "water_depth" in info.vars && "velocity" in info.vars
        fl = json_(rport, "/api/frames")
        @test [x.t for x in fl.frames] == run.t
        k = length(run.t)
        d = reinterpret(Float32, get_(rport, "/api/frames/$k/water_depth").body)
        @test d == Float32.(FloodA5Output.read_frame(f, "water_depth", k))
        @test json_(rport, "/api/cell/0").id == run.mesh.ids[1]          # file order
        @test json_(rport, "/api/cell/$(length(run.mesh) - 1)").id == run.mesh.ids[end]
        # a cell's series, read from the file in one go (time-series plots)
        c = length(run.mesh) ÷ 2
        sr = json_(rport, "/api/series/$(c - 1)?var=velocity")
        @test sr.id == run.mesh.ids[c] && collect(sr.t) == run.t && collect(sr.frames) == 1:k
        @test Float32.(sr.values) == Float32.(FloodA5Output.read_series(f, "velocity", c))
        @test Float32.(json_(rport, "/api/series/$(c - 1)").values)[k] == d[c]
        mb = get_(rport, "/api/mesh").body
        @test reinterpret(UInt32, mb[9:12])[1] == length(run.mesh)
        # the file's run record
        rec = json_(rport, "/api/run")
        @test rec.file == abspath(H5) && rec.label == basename(H5) && rec.flooda5_version == FLOODA5_VERSION
        @test collect(String, rec.argv) == args0 && rec.format_version == 2 && !isempty(rec.git_commit)
        @test rec.mesh.n_cells == length(run.mesh) && rec.frames.n == length(run.t) && rec.frames.t_last ≈ 300.0
        @test !haskey(rec, :project) && !haskey(rec, :manifest) && !haskey(rec, :git_diff)
        @test rec.mesh.cell_area_min_m2 ≈ rec.mesh.cell_area_max_m2 rtol = 0.05   # a 25 m Cartesian grid
        @test json_(rport, "/api/about").variables.water_depth.unit == "m"
    finally
        VisualisationServer.stop(vis); close(f)
    end

    # bare --vis is the Makie viewer now (stubbed here: it refuses to start)
    err = try
        with_logger(NullLogger()) do
            redirect_stdout(devnull) do
                main(["--meshgen-cartesian", AOI, "--dx", "50", "--flow-model", "standard",
                      "--sim-duration", "10", "--vis"])
            end
        end
        nothing
    catch e
        e
    end
    @test err isa ErrorException && occursin("--vis makie", err.msg)
end
