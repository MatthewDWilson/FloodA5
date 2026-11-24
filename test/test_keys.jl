# test/test_keys.jl — access keys (io/FloodA5Keys.jl) and their use by the
# web viewers' server.
#
#   julia --project=ci test/test_keys.jl
#
# Lookup order (environment, FLOODA5_KEYS, the repository's keys.json, the
# user's ~/.flooda5/keys.json), placeholders and empty values, a malformed
# file, key_source without the value, the template, the server's
# /api/config, keys.json never served, and keys.json git-ignored. The tests
# pass their own files, so a developer's real keys never affect them.

using Test, Logging
using HTTP, JSON3

isdefined(Main, :VisualisationServer) ||
    include(joinpath(@__DIR__, "..", "visualisation", "VisualisationServer.jl"))
const VS = VisualisationServer
const K = VS.FloodA5Keys
const ROOT = abspath(joinpath(@__DIR__, ".."))
const NAME = "FLOODA5_CESIUM_ION_TOKEN"

writejson(path, d) = (mkpath(dirname(path)); write(path, JSON3.write(d)); path)
free_port() = (l = HTTP.Sockets.listen(HTTP.Sockets.IPv4("127.0.0.1"), 0);
               p = Int(HTTP.Sockets.getsockname(l)[2]); close(l); p)

@testset "access keys" begin
    dir = mktempdir()
    a = writejson(joinpath(dir, "a", "keys.json"), Dict(NAME => " token-a ", "FLOODA5_GOOGLE_TILES_KEY" => ""))
    b = writejson(joinpath(dir, "b", "keys.json"), Dict(NAME => "token-b", "FLOODA5_GOOGLE_TILES_KEY" => "google-b"))
    ph = writejson(joinpath(dir, "ph", "keys.json"), Dict(NAME => "YOUR_CESIUM_ION_TOKEN_HERE"))
    bad = joinpath(dir, "bad", "keys.json"); mkpath(dirname(bad)); write(bad, "{ not json")
    none = joinpath(dir, "none", "keys.json")

    withenv(NAME => nothing, "FLOODA5_GOOGLE_TILES_KEY" => nothing) do
        @testset "files in order" begin
            @test K.flooda5_key(NAME; files = [a, b]) == "token-a"          # trimmed, first file wins
            @test K.flooda5_key(NAME; files = [b, a]) == "token-b"
            @test K.key_source(NAME; files = [a, b]) == a
            # an empty value or a placeholder counts as unset: the next file is used
            @test K.flooda5_key("FLOODA5_GOOGLE_TILES_KEY"; files = [a, b]) == "google-b"
            @test K.flooda5_key(NAME; files = [ph, b]) == "token-b"
            @test K.flooda5_key(NAME; files = [none, b]) == "token-b"        # missing file skipped
            @test K.flooda5_key(NAME; files = [none]) == ""
            @test K.flooda5_key(NAME; files = [none], default = "d") == "d"
            @test K.key_source(NAME; files = [none]) === nothing
        end
        @testset "malformed file" begin
            v = with_logger(NullLogger()) do
                K.flooda5_key(NAME; files = [bad, b])
            end
            @test v == "token-b"
            @test_logs (:warn, r"not a JSON object of keys") match_mode = :any begin
                empty!(K._WARNED); K.flooda5_key(NAME; files = [bad])
            end
        end
        @testset "default file list" begin
            fs = withenv("FLOODA5_KEYS" => a) do; K.keys_files(); end
            @test fs[1] == a
            @test fs[2] == joinpath(ROOT, "keys.json")
            @test fs[3] == joinpath(homedir(), ".flooda5", "keys.json")
            @test withenv(K.keys_files, "FLOODA5_KEYS" => nothing) == fs[2:3]
            # FLOODA5_KEYS outranks the repository's and the user's files
            @test withenv(() -> K.flooda5_key(NAME), "FLOODA5_KEYS" => b) == "token-b"
        end
    end

    @testset "environment first" begin
        withenv(NAME => "token-env") do
            @test K.flooda5_key(NAME; files = [a, b]) == "token-env"
            @test K.key_source(NAME; files = [a, b]) == "environment"
        end
        withenv(NAME => "  ") do                                           # blank: falls through
            @test K.flooda5_key(NAME; files = [b]) == "token-b"
        end
    end

    @testset "template" begin
        t = JSON3.read(read(joinpath(ROOT, "keys.example.json"), String), Dict{String,Any})
        @test sort([k for k in keys(t) if !startswith(k, "_")]) == sort(collect(K.KEY_NAMES))
        @test all(v -> v == "", (t[k] for k in K.KEY_NAMES))                 # no key in the template
    end

    @testset "git-ignored" begin
        # (.git is a directory in a clone, a file in a worktree)
        if Sys.which("git") !== nothing && ispath(joinpath(ROOT, ".git"))
            ig(p) = success(Cmd(`git check-ignore -q $p`; dir = ROOT))
            @test ig("keys.json")
            @test ig("visualisation/viewer-common/keys.json")
            @test !ig("keys.example.json")
        end
    end

    @testset "server" begin
        s = VS.start(; port = free_port())
        try
            withenv(NAME => nothing, "FLOODA5_GOOGLE_TILES_KEY" => nothing, "FLOODA5_KEYS" => b) do
                c = VS.viewer_config(s)
                @test c["ion_token"] == "token-b" && c["google_key"] == "google-b"
            end
            withenv(NAME => "token-env", "FLOODA5_KEYS" => b) do
                @test VS.viewer_config(s)["ion_token"] == "token-env"
            end
            # a keys.json placed among the viewer files is never served
            f = joinpath(s.viz_dir, "viewer-common", "keys.json")
            existed = isfile(f)
            existed || write(f, "{}")
            try
                r = HTTP.get("http://127.0.0.1:$(s.port)/viz/viewer-common/keys.json"; status_exception = false)
                @test r.status == 404
            finally
                existed || rm(f)
            end
        finally
            VS.stop(s)
        end
    end
end
