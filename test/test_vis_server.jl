# test/test_vis_server.jl — the web viewer's server (visualisation/VisualisationServer.jl)
#
#   julia --project=ci test/test_vis_server.jl
#
# Protocol 2 (docs/DATA_FORMATS.md §6): binary mesh layout, frame list and
# Float32 frames, the WebSocket hello (a client that connects after the run
# ended learns that it is complete), live-frame thinning and the memory
# budget, replay sources, static-file path checks, the Host-header guard and
# the viewer configuration. Also checks that the mesh encoding matches the
# fixture the JavaScript tests decode (visualisation/viewer-common/test/fixtures).

using Test
using HTTP, JSON3

isdefined(Main, :VisualisationServer) ||
    include(joinpath(@__DIR__, "..", "visualisation", "VisualisationServer.jl"))
const VS = VisualisationServer

const FIXTURE = joinpath(@__DIR__, "..", "visualisation", "viewer-common", "test", "fixtures", "mesh_small.bin")

# Three cells: a closed square ring (vector of points), an open triangle, and a
# pentagon given as a 2×k matrix (the HDF5 polygon layout).
function small_mesh()
    rings = Any[
        [[172.0, -43.0], [172.001, -43.0], [172.001, -42.999], [172.0, -42.999], [172.0, -43.0]],
        [[172.001, -43.0], [172.002, -43.0], [172.0015, -42.999]],
        [172.002 172.003 172.0035 172.0025 172.0015;
         -43.0  -43.0   -42.9995 -42.999   -42.9995],
    ]
    return VS.encode_mesh(rings, ["c1", "c2", "c3"], [10.0, 11.5, NaN], [9.5, 11.5, 12.0])
end

# Little-endian readers over the payload
r32(b, o) = ltoh(reinterpret(UInt32, b[o+1:o+4])[1])
f32(b, o) = ltoh(reinterpret(Float32, b[o+1:o+4])[1])
f64(b, o) = ltoh(reinterpret(Float64, b[o+1:o+8])[1])

free_port() = (s = HTTP.Sockets.listen(HTTP.Sockets.IPv4("127.0.0.1"), 0);
               p = Int(HTTP.Sockets.getsockname(s)[2]); close(s); p)
get_(port, path; headers = []) = HTTP.get("http://127.0.0.1:$port$path", headers;
                                         status_exception = false, retry = false)
json_(port, path) = JSON3.read(String(get_(port, path).body))

@testset "VisualisationServer" begin

    @testset "mesh encoding" begin
        m = small_mesh()
        b = m.bytes
        @test b[1:4] == b"FA5M"
        @test r32(b, 4) == 1 && r32(b, 8) == 3
        nv = r32(b, 12)
        @test nv == 4 + 3 + 5                       # closing repeat dropped
        olon = f64(b, 16); olat = f64(b, 24)
        @test olon ≈ (172.0 + 172.0035) / 2 && olat ≈ (-43.0 - 42.999) / 2
        offs = [r32(b, 32 + 4i) for i in 0:3]
        @test offs == [0, 4, 7, 12]
        o_dlon = 32 + 4 * 4; o_dlat = o_dlon + 4nv; o_elev = o_dlat + 4nv; o_base = o_elev + 4 * 3
        @test length(b) == o_base + 4 * 3
        @test olon + f32(b, o_dlon) ≈ 172.0 atol = 1e-7
        @test olat + f32(b, o_dlat + 4 * 6) ≈ -42.999 atol = 1e-7     # triangle apex
        @test f32(b, o_elev) == 10.0f0 && isnan(f32(b, o_elev + 8))
        @test f32(b, o_base + 4) == 11.5f0
        @test m.ids == ["c1", "c2", "c3"] && m.lon[1] ≈ 172.0005
        # the JavaScript tests decode this fixture: keep the two sides in step
        if get(ENV, "FLOODA5_UPDATE_FIXTURES", "") == "1"
            mkpath(dirname(FIXTURE)); write(FIXTURE, b)
        end
        @test isfile(FIXTURE) && read(FIXTURE) == b
        @test_throws ErrorException VS.encode_mesh([[[0.0, 0.0], [1.0, 0.0]]], ["x"], [0.0], [0.0])
        @test_throws ErrorException VS.encode_mesh([[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]], ["x", "y"], [0.0], [0.0])
    end

    @testset "depth datum" begin
        e = [1.0, 2.0, 3.0]
        @test VS.depth_datum(:standard, e, Dict()) == e
        @test VS.depth_datum("sgs", e, Dict("sgs_z_min" => [0.5, NaN, 2.0])) == [0.5, 2.0, 2.0]
        @test VS.depth_datum(:swe, e, Dict("sgs_z_min" => [0.0, 0.0, 0.0])) == e
    end

    @testset "frame budget" begin
        f = VS.LiveFrames(10 * 3 * 4)               # room for 10 frames of 3 Float32
        for k in 1:25
            push!(f.ids, k); push!(f.t, k); push!(f.data, Dict("h" => Float32[k, k, k]))
            f.bytes += 12
            VS._decimate!(f)
            @test f.bytes <= f.budget
            @test f.ids[1] == 1 && f.ids[end] == k  # first and latest always kept
            @test issorted(f.ids)
        end
    end

    @testset "live server" begin
        port = free_port()
        vis = VS.start(; port, min_interval = 0.2, info = Dict("label" => "test", "flow_model" => "standard"),
                       config = Dict("vertical_offset" => 1.5))
        try
            @test json_(port, "/api/info")["mesh_ready"] == false
            @test get_(port, "/api/mesh").status == 204
            m = small_mesh()
            VS.set_mesh!(vis, m)
            r = get_(port, "/api/mesh")
            @test r.status == 200 && r.body == m.bytes

            @test VS.wants_frame(vis)
            id1 = VS.push_frame!(vis, 10.0, Dict("water_depth" => [0.0, 0.5, 1.0], "velocity" => [0.0, 0.1, 0.2]))
            @test !VS.wants_frame(vis)              # thinned until min_interval has passed
            sleep(0.25)
            @test VS.wants_frame(vis)
            id2 = VS.push_frame!(vis, 20.0, Dict("water_depth" => [0.0, 0.6, 1.1], "velocity" => [0.0, 0.1, 0.3]))
            fl = json_(port, "/api/frames")
            @test fl["complete"] == false
            @test [x["id"] for x in fl["frames"]] == [id1, id2] && fl["frames"][2]["t"] == 20.0
            d = get_(port, "/api/frames/$id2/water_depth")
            @test d.status == 200 && reinterpret(Float32, d.body) == Float32[0.0, 0.6, 1.1]
            @test get_(port, "/api/frames/$id2/nonesuch").status == 404
            @test get_(port, "/api/frames/999/water_depth").status == 404
            @test get_(port, "/api/frames/x/water_depth").status == 400
            info = json_(port, "/api/info")
            @test info["protocol"] == 2 && info["mode"] == "live" && info["frames"] == 2
            @test sort(String.(info["vars"])) == ["velocity", "water_depth"] && info["n_cells"] == 3
            c = json_(port, "/api/cell/1?frame=$id2")
            @test c["id"] == "c2" && c["values"]["water_depth"] ≈ 0.6 && c["base"] == 11.5
            @test json_(port, "/api/cell/2")["elevation"] === nothing      # NaN → null
            @test get_(port, "/api/cell/3").status == 404
            cfg = json_(port, "/api/config")
            @test cfg["protocol"] == 2 && cfg["vertical_offset"] == 1.5
            @test endswith(cfg["cesium_url"], "/") && haskey(cfg, "ion_token")
            @test occursin("chart.js@", cfg["chartjs_url"])

            # a cell's series over the live frames (time-series plots)
            sr = json_(port, "/api/series/1?var=water_depth")
            @test sr["index"] == 1 && sr["id"] == "c2" && sr["var"] == "water_depth"
            @test sr["frames"] == [id1, id2] && sr["t"] == [10.0, 20.0]
            @test Float32.(sr["values"]) == Float32[0.5, 0.6]
            @test Float32.(json_(port, "/api/series/2?var=velocity")["values"]) == Float32[0.2, 0.3]
            @test json_(port, "/api/series/0")["var"] == "water_depth"           # the default variable
            @test get_(port, "/api/series/1?var=nonesuch").status == 404
            @test get_(port, "/api/series/3").status == 404
            @test get_(port, "/api/series/x").status == 404

            # the run record (Run tab) and the variable descriptions (About tab)
            @test isempty(json_(port, "/api/run"))                       # none set yet
            VS.set_run!(vis, Dict("flow_model" => "standard", "argv" => ["--dx", "25"],
                                  "mesh" => Dict("n_cells" => 3)))
            rec = json_(port, "/api/run")
            @test rec["flow_model"] == "standard" && rec["argv"] == ["--dx", "25"] && rec["mesh"]["n_cells"] == 3
            ab = json_(port, "/api/about")
            for v in ("water_depth", "level", "velocity", "saturation", "elevation")
                @test !isempty(ab["variables"][v]["label"]) && !isempty(ab["variables"][v]["definition"])
            end
            @test haskey(ab["topics"], "colour_range") && haskey(ab["topics"], "vertical_offset")

            # colour ranges over the run: wet cells (depth ≥ wet) of every frame so far
            rr = json_(port, "/api/ranges?wet=0.05")
            @test rr["complete"] == false && rr["frames"] == 2 && rr["of"] == 2 && rr["wet"] == 0.05
            dr = rr["vars"]["water_depth"]                      # wet depths 0.5, 1.0, 0.6, 1.1
            @test dr["n"] == 4 && isapprox(dr["min"], 0.5; atol = 1e-6) && isapprox(dr["max"], 1.1; atol = 1e-6)
            @test isapprox(dr["p02"], 0.506; atol = 1e-5) && isapprox(dr["p98"], 1.094; atol = 1e-5)
            @test rr["vars"]["velocity"]["n"] == 4 && isapprox(rr["vars"]["velocity"]["max"], 0.3; atol = 1e-6)
            lv = rr["vars"]["level"]                            # base + depth (cell 2: 11.5 + 0.5, …)
            @test lv["n"] == 4 && isapprox(lv["min"], 12.0; atol = 1e-5)
            @test json_(port, "/api/ranges?wet=0.55")["vars"]["water_depth"]["n"] == 3
            @test json_(port, "/api/ranges?wet=nonsense")["wet"] == 0.05
            @test json_(port, "/api/ranges?wet=0")["vars"]["water_depth"]["n"] == 6

            # WebSocket: frame and complete messages to a connected client …
            msgs = Channel{Any}(32)
            ws_task = @async HTTP.WebSockets.open("ws://127.0.0.1:$port/api/live") do ws
                for m in ws
                    put!(msgs, JSON3.read(m))
                    JSON3.read(m)["type"] == "complete" && break
                end
            end
            hello = take!(msgs)
            @test hello["type"] == "hello" && hello["complete"] == false && hello["frames"] == 2 &&
                  hello["mesh_ready"] == true
            sleep(0.25)
            id3 = VS.push_frame!(vis, 30.0, Dict("water_depth" => [0.1, 0.7, 1.2], "velocity" => [0.0, 0.2, 0.3]))
            fm = take!(msgs)
            @test fm["type"] == "frame" && fm["id"] == id3 && fm["t"] == 30.0
            VS.notify_complete!(vis)
            @test take!(msgs)["type"] == "complete"
            wait(ws_task)
            # … and a client that connects after the end is told so in its hello
            lateref = Ref{Any}(nothing)
            HTTP.WebSockets.open("ws://127.0.0.1:$port/api/live") do ws
                lateref[] = JSON3.read(HTTP.WebSockets.receive(ws))
            end
            late = lateref[]
            @test late["type"] == "hello"
            @test late["complete"] == true
            @test late["frames"] == 3
            @test json_(port, "/api/frames")["complete"] == true
            # the ranges follow the new frame (not the cached two-frame answer)
            rr = json_(port, "/api/ranges?wet=0.05")
            @test rr["complete"] == true && rr["frames"] == 3 && rr["vars"]["water_depth"]["n"] == 7

            # static files, path checks, Host guard
            for v in ("cesium", "giro3d")                 # both viewers are served
                r = get_(port, "/viz/$v/index.html")
                @test r.status == 200 && occursin("text/html", HTTP.header(r, "Content-Type"))
                @test get_(port, "/viz/$v/js/main.js").status == 200
            end
            r = get_(port, "/viz/viewer-common/js/app.js")
            @test r.status == 200 && occursin("javascript", HTTP.header(r, "Content-Type"))
            @test get_(port, "/").status == 200           # redirected to the viewer
            r = HTTP.get("http://127.0.0.1:$port/?basemap=osm"; redirect = false, status_exception = false)
            @test r.status == 302 && HTTP.header(r, "Location") == "/viz/cesium/index.html?basemap=osm"
            @test get_(port, "/viz/VisualisationServer.jl").status == 404     # only the viewer directories
            @test get_(port, "/viz/plots/makie.jl").status == 404
            @test get_(port, "/viz/../FloodModel.jl").status in (403, 404)
            @test get_(port, "/viz/cesium/../../FloodModel.jl").status in (403, 404)
            @test get_(port, "/viz/cesium/%2e%2e/%2e%2e/FloodModel.jl").status in (403, 404)
            @test get_(port, "/viz/cesium/..%5C..%5CFloodModel.jl").status in (403, 404)
            @test get_(port, "/viz/viewer-common/config.json").status == 404
            @test get_(port, "/viz/viewer-common/config.example.json").status == 200
            @test get_(port, "/viz/cesium/nonesuch.js").status == 404
            @test get_(port, "/api/info"; headers = ["Host" => "evil.example:$port"]).status == 403
            @test get_(port, "/api/info"; headers = ["Host" => "localhost:$port"]).status == 200
            @test HTTP.post("http://127.0.0.1:$port/api/info"; status_exception = false).status == 405

            # The server answers while this (main) task computes, as the time
            # loop does: in Julia 1.12 the main task runs on interactive thread 1,
            # the pool the server uses, so the loop yields once per "step"
            # (FloodModel.jl run_simulation! does the same with --vis cesium).
            t0 = time(); lat = Ref(NaN)
            probe = Threads.@spawn begin
                t1 = time(); get_(port, "/api/info"); lat[] = time() - t1
            end
            x = 0.0
            while time() - t0 < 2.0
                t_step = time()
                while time() - t_step < 0.01     # a 10 ms step without yielding
                    x += sin(x)
                end
                yield()
            end
            wait(probe)
            @test lat[] < 0.5
        finally
            VS.stop(vis)
        end
        sleep(0.2)
        @test VS._port_free("127.0.0.1", port)
    end

    @testset "replay server" begin
        port = free_port()
        vis = VS.start(; port)
        try
            VS.set_mesh!(vis, small_mesh())
            reads = Ref(0)
            VS.set_replay!(vis, [0.0, 60.0, 120.0], ["water_depth"],
                           (k, v) -> (reads[] += 1; Float32[k, k + 1, k + 2]))
            @test !VS.wants_frame(vis)
            @test_throws ErrorException VS.push_frame!(vis, 1.0, Dict("water_depth" => [0.0, 0, 0]))
            fl = json_(port, "/api/frames")
            @test fl["complete"] == true && [x["t"] for x in fl["frames"]] == [0.0, 60.0, 120.0]
            @test reinterpret(Float32, get_(port, "/api/frames/2/water_depth").body) == Float32[2, 3, 4]
            get_(port, "/api/frames/2/water_depth")
            @test reads[] == 1                       # cached
            @test get_(port, "/api/frames/4/water_depth").status == 404
            @test json_(port, "/api/info")["mode"] == "replay"
            # series without a series reader: frame by frame (k, k+1, k+2 at cell 2 → k+1)
            sr = json_(port, "/api/series/1?var=water_depth")
            @test sr["frames"] == [1, 2, 3] && sr["t"] == [0.0, 60.0, 120.0] && sr["values"] == [2, 3, 4]
            @test get_(port, "/api/series/1?var=velocity").status == 404
            # ranges over a replay: frames k = 1, 2, 3 hold k, k+1, k+2
            rr = json_(port, "/api/ranges")
            @test rr["complete"] == true && rr["frames"] == 3 && rr["wet"] == 0.05
            @test rr["vars"]["water_depth"]["n"] == 9 && rr["vars"]["water_depth"]["min"] == 1 &&
                  rr["vars"]["water_depth"]["max"] == 5
            n_reads = reads[]
            json_(port, "/api/ranges")
            @test reads[] == n_reads                  # cached
            @test VS.run_ranges(vis; max_frames = 2).frames == 2      # frames spread over the run
            @test VS.run_ranges(vis; max_values = 3).vars["water_depth"].n == 3   # cells subsampled
            # stop() releases connected clients (their read loops end)
            n_seen = Ref(0)
            client = @async try
                HTTP.WebSockets.open("ws://127.0.0.1:$port/api/live") do ws
                    for _ in ws
                        n_seen[] += 1
                    end
                end
            catch
            end
            t0 = time(); while n_seen[] == 0 && time() - t0 < 10; sleep(0.05); end
            @test n_seen[] >= 1                       # the hello
            VS.stop(vis)
            t0 = time(); while !istaskdone(client) && time() - t0 < 10; sleep(0.05); end
            @test istaskdone(client)
        finally
            VS.stop(vis)
        end
    end

    @testset "hazard variables and statistics" begin
        port = free_port()
        vis = VS.start(; port)
        try
            VS.set_mesh!(vis, small_mesh())
            # depth and speed over three frames (cells 1–3)
            D = Float32[0.0 0.3 1.0; 0.2 0.5 0.0; 0.4 2.0 NaN]       # frame k = column k
            V = Float32[0.0 1.0 2.5; 0.5 0.5 0.0; 0.1 0.5 0.0]
            VS.set_replay!(vis, [0.0, 60.0, 120.0], ["water_depth", "velocity"],
                           (k, v) -> v == "water_depth" ? D[:, k] : V[:, k])
            info = json_(port, "/api/info")
            @test sort(String.(info["derived"])) == ["depth_velocity", "hazard_class_arr", "hazard_degree", "hazard_rating"]
            @test String.(info["stats"]) == collect(VS.STAT_VARS)
            f32b(r) = reinterpret(Float32, r.body)
            dv = f32b(get_(port, "/api/frames/2/depth_velocity"))
            @test dv ≈ Float32[0.3, 0.25, 1.0]
            # Defra HR with the land use from the query (urban: DF 0 / 1 / 1)
            hr = f32b(get_(port, "/api/frames/2/hazard_rating?landuse=urban"))
            @test hr ≈ Float32[0.3 * 1.5 + 1, 0.5 * 1.0 + 1, 2.0 * 1.0 + 1]
            @test f32b(get_(port, "/api/frames/2/hazard_rating"))[1] ≈ 0.3 * 1.5 + 1      # fd2320
            @test f32b(get_(port, "/api/frames/2/hazard_rating?landuse=pasture"))[1] ≈ 0.45
            # the same code as FloodA5Analysis would run: analysis/hazard.jl
            p = VS.HazardParams(mass = 50.0, height = 1.6)
            hd = f32b(get_(port, "/api/frames/2/hazard_degree?mass=50&height=1.6"))
            @test hd ≈ VS.derive_hazard("hazard_degree", D[:, 2], V[:, 2], p)
            @test f32b(get_(port, "/api/frames/1/hazard_class_arr")) == Float32[0, 1, 2]
            @test isnan(f32b(get_(port, "/api/frames/3/depth_velocity"))[3])
            # a cell's series derives the same way as a frame (any length)
            @test VS.derived_values("depth_velocity", n -> n == "water_depth" ? D[1, :] : V[1, :], p) ≈ Float32[0, 0.3, 2.5]
            @test VS.derived_values("depth_velocity", n -> nothing, p) === nothing
            # … and over HTTP: /api/series of a hazard variable, with the parameters
            sr = json_(port, "/api/series/0?var=hazard_rating&landuse=urban")
            @test sr["t"] == [0.0, 60.0, 120.0]
            @test Float32.(something.(sr["values"], NaN)) ≈
                  VS.derive_hazard("hazard_rating", D[1, :], V[1, :], VS.HazardParams(landuse = "urban")) nans = true
            @test Float32.(sr["values"][2:3]) ≈ Float32[0.3 * 1.5 + 1, 1.0 * 3.0 + 1]
            @test Float32.(json_(port, "/api/series/0?var=depth_velocity")["values"]) ≈ Float32[0, 0.3, 2.5]
            @test Float32.(json_(port, "/api/series/1?var=hazard_degree&mass=50&height=1.6")["values"]) ≈
                  VS.derive_hazard("hazard_degree", D[2, :], V[2, :], p)
            @test get_(port, "/api/series/0?var=max_water_depth").status == 404          # statistics: no series
            c = json_(port, "/api/cell/0?frame=2&landuse=urban")
            @test c["values"]["hazard_rating"] ≈ 1.45 && c["values"]["depth_velocity"] ≈ 0.3
            # ranges include the hazard variables (with the parameters)
            rr = json_(port, "/api/ranges?wet=0.05&landuse=none")
            @test rr["vars"]["depth_velocity"]["max"] ≈ 2.5 && rr["vars"]["hazard_rating"]["max"] ≈ 3.0
            @test json_(port, "/api/ranges?wet=0.05&landuse=urban")["vars"]["hazard_rating"]["max"] ≈ 4.0
            # statistics over the run
            st(name, q = "") = f32b(get_(port, "/api/stats/$name$q"))
            @test st("max_water_depth") == Float32[1.0, 0.5, 2.0]
            @test st("max_velocity") == Float32[2.5, 0.5, 0.5]
            @test st("max_depth_velocity") ≈ Float32[2.5, 0.25, 1.0]
            t1 = st("time_first_flooded", "?flood_depth=0.25")
            @test t1[1:2] == Float32[60, 60] && t1[3] == 0
            @test st("flood_duration", "?flood_depth=0.25") == Float32[90, 60, 90]  # weights 30, 60, 30
            @test isnan(st("time_first_flooded", "?flood_depth=5")[1])
            @test get_(port, "/api/stats/nonesuch").status == 404
            @test length(vis.stats) == 3                                   # one entry per parameter set
            st("max_water_depth", "?flood_depth=0.25")
            @test length(vis.stats) == 3                                   # cached
        finally
            VS.stop(vis)
        end
        # live: statistics over the frames kept so far, recomputed as frames arrive
        port = free_port()
        vis = VS.start(; port, min_interval = 0.0)
        try
            VS.set_mesh!(vis, small_mesh())
            VS.push_frame!(vis, 0.0, Dict("water_depth" => [0.0, 0.2, 0.4]))
            @test isempty(json_(port, "/api/info")["derived"])            # no velocity: no hazard
            @test String.(json_(port, "/api/info")["stats"]) == ["max_water_depth", "time_max_water_depth", "time_first_flooded", "flood_duration"]
            @test reinterpret(Float32, get_(port, "/api/stats/max_water_depth").body) == Float32[0, 0.2, 0.4]
            @test get_(port, "/api/stats/max_velocity").status == 404
            VS.push_frame!(vis, 10.0, Dict("water_depth" => [0.5, 0.1, 0.4]))
            @test reinterpret(Float32, get_(port, "/api/stats/max_water_depth").body) == Float32[0.5, 0.2, 0.4]
        finally
            VS.stop(vis)
        end
    end

    @testset "port in use" begin
        port = free_port()
        vis = VS.start(; port)
        try
            @test_throws ErrorException VS.start(; port)
        finally
            VS.stop(vis)
        end
    end

    @testset "shared palettes (palettes.json)" begin
        good = Dict("flooda5_palettes" => 1,
                    "palettes" => [Dict("name" => " Team depth ",
                                        "stops" => [Dict("pos" => 0, "color" => "#C6EBFF"),
                                                    Dict("pos" => 1, "color" => "#0a1964")])])
        ps, err = VS.validate_palettes(JSON3.read(JSON3.write(good)))
        @test err === nothing && length(ps) == 1
        @test ps[1]["name"] == "Team depth" && ps[1]["stops"][1]["color"] == "#c6ebff"
        stops = good["palettes"][1]["stops"]
        bad_palettes = Any[
            Dict("name" => "x"),                                             # no stops
            Dict("name" => "x", "stops" => stops, "extra" => 1),             # unknown key
            Dict("name" => "", "stops" => stops),
            Dict("name" => "x"^41, "stops" => stops),
            Dict("name" => "a" * Char(7) * "b", "stops" => stops),           # control character
            Dict("name" => 42, "stops" => stops),
            Dict("name" => "x", "stops" => stops[1:1]),
            Dict("name" => "x", "stops" => [Dict("pos" => i / 32, "color" => "#000000") for i in 0:32]),
            Dict("name" => "x", "stops" => reverse(stops)),                  # unsorted
            Dict("name" => "x", "stops" => [stops[1], Dict("pos" => 1.5, "color" => "#000000")]),
            Dict("name" => "x", "stops" => [stops[1], Dict("pos" => "1", "color" => "#000000")]),
            Dict("name" => "x", "stops" => [stops[1], Dict("pos" => true, "color" => "#000000")]),
            Dict("name" => "x", "stops" => [stops[1], Dict("pos" => 1, "color" => "blue")]),
            Dict("name" => "x", "stops" => [stops[1], Dict("pos" => 1, "color" => "#fff")]),
            Dict("name" => "x", "stops" => [stops[1], Dict("pos" => 1, "color" => "#000000", "a" => 1)]),
        ]
        for p in bad_palettes
            d = Dict("flooda5_palettes" => 1, "palettes" => [p])
            @test VS.validate_palettes(JSON3.read(JSON3.write(d)))[1] === nothing
        end
        for d in Any[[], Dict("palettes" => [good["palettes"][1]]),
                     Dict("flooda5_palettes" => 2, "palettes" => good["palettes"]),
                     Dict("flooda5_palettes" => 1, "palettes" => good["palettes"], "x" => 1),
                     Dict("flooda5_palettes" => 1, "palettes" => []),
                     Dict("flooda5_palettes" => 1, "palettes" => [Dict("name" => "p$i", "stops" => stops) for i in 1:65]),
                     Dict("flooda5_palettes" => 1, "palettes" => [Dict("name" => "A", "stops" => stops),
                                                                  Dict("name" => "a", "stops" => stops)])]
            @test VS.validate_palettes(JSON3.read(JSON3.write(d)))[1] === nothing
        end

        # Over HTTP, from a viewer directory with palettes.json beside config.json
        viz = mktempdir()
        common = joinpath(viz, "viewer-common")
        mkpath(common)
        write(joinpath(common, "app.js"), "// app")
        write(joinpath(common, "config.json"), JSON3.write(Dict("cesium_ion_token" => "SECRET-CONFIG")))
        write(joinpath(common, "keys.json"), JSON3.write(Dict("FLOODA5_CESIUM_ION_TOKEN" => "SECRET-KEYS")))
        pal = joinpath(common, "palettes.json")
        port = free_port()
        vis = VS.start(; port, viz_dir = viz)
        try
            r = json_(port, "/api/palettes")                      # absent: nothing offered, no error
            @test r["flooda5_palettes"] == 1 && isempty(r["palettes"]) && !haskey(r, "error")
            write(pal, JSON3.write(good))
            r = json_(port, "/api/palettes")                      # read on each request
            @test length(r["palettes"]) == 1 && r["palettes"][1]["name"] == "Team depth" && !haskey(r, "error")
            @test r["palettes"][1]["stops"][2]["pos"] == 1
            write(pal, "{ not json")
            r = json_(port, "/api/palettes")
            @test isempty(r["palettes"]) && occursin("palettes.json", r["error"])
            write(pal, JSON3.write(Dict("flooda5_palettes" => 1, "palettes" => [bad_palettes[9]])))
            r = json_(port, "/api/palettes")
            @test isempty(r["palettes"]) && occursin("sorted", r["error"])
            big = Dict("flooda5_palettes" => 1, "palettes" => [Dict("name" => "x", "stops" => stops)])
            write(pal, JSON3.write(big) * " "^(300 * 1024))      # valid JSON, but over the size limit
            r = json_(port, "/api/palettes")
            @test isempty(r["palettes"]) && occursin("256 KiB", r["error"])
            write(pal, JSON3.write(good))
            # palettes.json only through /api/palettes; config.json and keys.json never
            @test get_(port, "/viz/viewer-common/app.js").status == 200
            for f in ("palettes.json", "config.json", "keys.json")
                r = get_(port, "/viz/viewer-common/$f")
                @test r.status == 404
                @test !occursin("SECRET", String(r.body))
            end
            @test !occursin("SECRET", String(get_(port, "/api/palettes").body))
            @test get_(port, "/viz/viewer-common/../viewer-common/keys.json").status in (403, 404)
        finally
            VS.stop(vis)
        end
    end

    @testset "viewer choice" begin
        port = free_port()
        vis = VS.start(; port, viewer = "giro3d")
        try
            r = HTTP.get("http://127.0.0.1:$port/"; redirect = false, status_exception = false)
            @test r.status == 302 && HTTP.header(r, "Location") == "/viz/giro3d/index.html"
            @test json_(port, "/api/info")["viewer"] == "giro3d"
            # the Giro3D page pins its packages in an import map, and the
            # engine's version constant matches the pinned Giro3D
            html = String(get_(port, "/viz/giro3d/index.html").body)
            @test occursin("\"importmap\"", html)
            pinned = match(r"@giro3d/giro3d@([0-9.]+)/", html)
            eng = String(get_(port, "/viz/giro3d/js/engine.js").body)
            @test pinned !== nothing && occursin("GIRO3D_VERSION = '$(pinned[1])'", eng)
            @test occursin("three@0.180.0", html)
        finally
            VS.stop(vis)
        end
        @test_throws ErrorException VS.start(; port = free_port(), viewer = "nonesuch")
    end
end
