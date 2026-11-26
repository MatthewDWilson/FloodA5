# test/test_vis_local.jl — the globe viewer's local files (0.8.10): an output
# file read in the browser must give what the server gives.
#
#   julia --project=ci test/test_vis_local.jl
#
# - The committed fixture (visualisation/viewer-common/test/fixtures/
#   local_run.json, from test/vis_local_fixture.jl) is still what the server
#   answers: regenerated and compared (the run's creation time and file paths
#   aside). If the server's answers changed on purpose, rerun
#   `julia --project=ci test/vis_local_fixture.jl --write` and update the
#   JavaScript port (viewer-common/js/local/).
# - The JavaScript port answers a freshly generated fixture as the server
#   (node, when installed: viewer-common/test/local.test.js).
# - The viewer without a server: about.json is io/variables.toml, defaults.json
#   is what the server's /api/config starts from (h5wasm_url, FLOODA5_H5WASM_URL).
# - web_replay.jl without a file: a server with no run (config no_data).

ENV["FLOODA5_NO_VIS"] = "1"
using Test, HTTP, JSON3, Logging
import TOML
include(joinpath(@__DIR__, "vis_local_fixture.jl"))
const VS = VisualisationServer
const VC = joinpath(@__DIR__, "..", "visualisation", "viewer-common")

# plain, order-independent data: dicts as sorted pairs, numbers as Float64
_canon(x::AbstractDict) = [String(k) => _canon(v) for (k, v) in sort!(collect(pairs(x)); by = p -> String(p[1]))]
_canon(x::AbstractVector) = Any[_canon(v) for v in x]
_canon(x::Real) = x isa Bool ? x : Float64(x)
_canon(x::AbstractString) = String(x)
_canon(x) = x

# what changes from one generation to the next: the times written, where the files were
const _VOLATILE = ("created_utc", "finished_utc", "file", "source")
_strip(x::AbstractDict) = Dict{String,Any}(String(k) => _strip(v) for (k, v) in x if !(String(k) in _VOLATILE))
_strip(x::AbstractVector) = Any[_strip(v) for v in x]
_strip(x) = x

free_port_() = (s = HTTP.Sockets.listen(HTTP.Sockets.IPv4("127.0.0.1"), 0);
                p = Int(HTTP.Sockets.getsockname(s)[2]); close(s); p)

@testset "viewer: local files" begin
    dir = mktempdir()
    fresh = make_local_fixture(dir)

    @testset "committed fixture is current" begin
        @test isfile(LOCAL_FIXTURE)
        old = JSON3.read(read(LOCAL_FIXTURE, String))
        new = JSON3.read(JSON3.write(fresh))
        @test _canon(_strip(old["picks"])) == _canon(_strip(new["picks"]))
        @test length(old["runs"]) == length(new["runs"])
        for (a, b) in zip(old["runs"], new["runs"])
            @test a["name"] == b["name"]
            same = _canon(_strip(a["h5"])) == _canon(_strip(b["h5"]))
            same || @warn "fixture contents changed: $(a["name"]) — regenerate with test/vis_local_fixture.jl --write"
            @test same
            @test length(a["requests"]) == length(b["requests"])
            for (qa, qb) in zip(a["requests"], b["requests"])
                ok = _canon(_strip(qa)) == _canon(_strip(qb))
                ok || @warn "fixture answer changed: $(qa["path"]) — regenerate with test/vis_local_fixture.jl --write"
                @test ok
            end
        end
    end

    @testset "the JavaScript port answers as the server (node)" begin
        node = Sys.which("node")
        if node === nothing
            @info "node not found: JavaScript check skipped"
        else
            fx = joinpath(dir, "local_run.json")
            write(fx, JSON3.write(fresh))
            cmd = Cmd(`$node --test test/local.test.js`; dir = VC)
            ok = success(pipeline(addenv(cmd, "FLOODA5_LOCAL_FIXTURE" => fx); stdout = devnull, stderr = devnull))
            ok || run(ignorestatus(addenv(cmd, "FLOODA5_LOCAL_FIXTURE" => fx)))    # show the failures
            @test ok
        end
    end

    @testset "no server: about.json and defaults.json" begin
        toml = TOML.parsefile(joinpath(@__DIR__, "..", "io", "variables.toml"))
        about = JSON3.read(read(joinpath(VC, "about.json"), String))
        same = _canon(about) == _canon(toml)
        same || @warn "about.json is out of date: julia --project=. visualisation/viewer-common/tools/make_about.jl"
        @test same
        d = JSON3.read(read(joinpath(VC, "defaults.json"), String))
        @test VS.DEFAULT_CESIUM_URL == d["cesium_url"]
        @test VS.DEFAULT_CHARTJS_URL == d["chartjs_url"]
        @test VS.DEFAULT_H5WASM_URL == d["h5wasm_url"]
        @test occursin("h5wasm@", d["h5wasm_url"]) && endswith(d["h5wasm_url"], ".js")
    end

    @testset "web_replay without a file; h5wasm_url" begin
        port = free_port_()
        withenv("FLOODA5_H5WASM_URL" => nothing) do
            vis = with_logger(NullLogger()) do
                web_replay_empty(; port, viewer = "giro3d")
            end
            try
                g(p) = HTTP.get("http://127.0.0.1:$port$p"; status_exception = false, retry = false)
                cfg = JSON3.read(String(g("/api/config").body))
                @test cfg["no_data"] == true
                @test cfg["h5wasm_url"] == VS.DEFAULT_H5WASM_URL
                @test g("/api/mesh").status == 204
                info = JSON3.read(String(g("/api/info").body))
                @test info["mesh_ready"] == false && info["frames"] == 0
                for f in ("defaults.json", "about.json", "js/local/worker.js", "js/local/localrun.js",
                          "img/GoogleMaps_Logo_WithDarkOutline.svg")       # (0.8.16: Giro3D's Google credit)
                    r = g("/viz/viewer-common/$f")
                    @test r.status == 200
                end
                @test startswith(HTTP.header(g("/viz/viewer-common/js/local/worker.js"), "Content-Type"), "text/javascript")
                @test HTTP.header(g("/viz/viewer-common/img/GoogleMaps_Logo_WithDarkOutline.svg"), "Content-Type") == "image/svg+xml"
                withenv("FLOODA5_H5WASM_URL" => "/viz/viewer-common/vendor/hdf5_hl.js") do
                    @test JSON3.read(String(g("/api/config").body))["h5wasm_url"] == "/viz/viewer-common/vendor/hdf5_hl.js"
                end
                @test HTTP.get("http://127.0.0.1:$port/"; redirect = false, status_exception = false, retry = false).status == 302
            finally
                VS.stop(vis)
            end
        end
    end
end
