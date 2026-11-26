# test/vis_local_fixture.jl — the fixture that ties the viewer's local-file
# reader (visualisation/viewer-common/js/local/: the /api answered in the
# browser) to the server (VisualisationServer.jl + web_replay.jl).
#
#   julia --project=ci test/vis_local_fixture.jl --write
#
# Builds two small output files (an SGS run with NaNs, non-finite options and
# closed and open polygon rings; a depth-only standard run), serves each with
# web_replay_server and records, as JSON: the files' contents (every group,
# dataset and attribute, as the JavaScript test reader serves them) and the
# server's answer to a list of requests. The JavaScript tests
# (viewer-common/test/local.test.js) answer the same requests from the same
# contents and must agree; test/test_vis_local.jl checks that this committed
# fixture is still what the server answers. After a deliberate change to the
# server's answers: rerun with --write and update the JavaScript port.

isdefined(Main, :FloodA5Output) || include(joinpath(@__DIR__, "..", "io", "FloodA5Output.jl"))
isdefined(Main, :web_replay_server) || include(joinpath(@__DIR__, "..", "visualisation", "web_replay.jl"))
using HDF5, HTTP, JSON3, Base64, Logging
isdefined(Main, :FO) || (const FO = FloodA5Output)

const LOCAL_FIXTURE = joinpath(@__DIR__, "..", "visualisation", "viewer-common", "test", "fixtures", "local_run.json")

# ---------------------------------------------------------------------------
# The runs
# ---------------------------------------------------------------------------

# a square cell of 0.001° at column c, row r; `closed`: repeat the first vertex
_square(c, r; closed = false) = (p = [[172.6 + 0.001c, -43.5 + 0.001r], [172.601 + 0.001c, -43.5 + 0.001r],
                                       [172.601 + 0.001c, -43.499 + 0.001r], [172.6 + 0.001c, -43.499 + 0.001r]];
                                  closed ? [p; [p[1]]] : p)

function _sgs_run(path)
    polys = Any[_square(0, 0; closed = true), _square(1, 0),
                [[172.602, -43.5], [172.603, -43.5], [172.6025, -43.499]],                       # triangle, open
                [[172.6, -43.499], [172.601, -43.499], [172.6005, -43.498], [172.6, -43.499]],  # triangle, closed
                [[172.601, -43.499], [172.602, -43.499], [172.6022, -43.4985], [172.6015, -43.498], [172.6008, -43.4985]],
                _square(2, 1)]
    n = length(polys)
    clon = [sum(p[k][1] for k in 1:3) / 3 for p in polys]; clat = [sum(p[k][2] for k in 1:3) / 3 for p in polys]
    elev = [10.0, 10.5, 11.25, 9.75, NaN, 12.0]
    cells = (ids = ["c$i" for i in 1:n], center_lons = clon, center_lats = clat, elevations = elev,
             static_vars = Dict("sgs_z_min" => [9.5, 10.0, NaN, 9.0, 11.0, 11.5], "friction" => fill(0.035, n)),
             id_scheme = "polygon", cell_area = [8000.0, 8000.0, 4000.0, 4000.0, NaN, 8000.0], polygons = polys,
             mesh_meta = Dict{String,Any}("mesh_family" => "fixture", "dggs_equal_area" => true, "target_dx_m" => 90.0))
    run = Dict{String,Any}("argv" => ["--flow-model", "sgs", "--output", "x.h5"],
                           "options" => Dict{String,Any}("flow_model" => "sgs", "courant" => 0.5, "t_end" => Inf,
                                                         "eta_filter" => NaN, "floor" => -Inf, "name" => "Infinity NaN"),
                           "flow_model" => "sgs", "git_commit" => "0123abc", "git_dirty" => true,
                           "git_diff" => "line 1\nline 2\n", "threads" => 4, "backend" => "cpu",
                           "precision" => "float64", "project" => "[deps]\n")
    FO.create_output(path; flooda5_version = "0.8.10-fixture", cells, run)
    ts = [0.0, 60.0, 120.0, 300.0, 600.0, 900.0, 1800.0]
    for (k, t) in enumerate(ts)
        d = [round(max(0.0, 1.2 * sin(0.7i + 0.9k) + 0.4 * (k - 2) / 3), digits = 3) for i in 1:n]
        v = [round(abs(2.6 * sin(1.3i - 0.5k)) + (i == 5 ? 1.9 : 0.0), digits = 3) for i in 1:n]
        k <= 2 && (d[6] = NaN)
        k == 4 && (v[3] = NaN)
        k == 5 && (d[2] = 2.1)                     # deeper than the person: floating
        FO.append_frame!(path, t, ["water_depth" => d, "velocity" => v,
                                   "vel_u" => v .* cos.(1:n), "vel_v" => v .* sin.(1:n),
                                   "saturation" => min.(1.0, d ./ 0.3), "volume" => d .* 100])
    end
    FO.write_summary!(path, Dict{String,Any}("completed" => true, "n_frames" => length(ts),
                                             "mass_balance_error_m3" => 1.5e-9, "vol_added_m3" => 12.0))
    # the model's statistics (0.8.13), with values of their own: the answers show which were used
    FO.write_statistics!(path, Dict(
        "max_water_depth" => [1.5, 2.5, 0.9, 1.8, 1.25, NaN], "time_max_water_depth" => [900, 600, 1800, 300, 900, NaN],
        "time_first_flooded" => [60, 0, 120, 60, 300, NaN], "flood_duration" => [1500, 1800, 1200, 1700, 900, 0],
        "max_velocity" => [2.5, 3.1, 0.7, 2.9, 4.6, NaN], "max_depth_velocity" => [2.1, 4.0, 0.5, 3.0, 4.4, NaN],
        "max_hazard_class_arr" => [5, 5, 2, 5, 6, 0]),
        Dict("cadence" => "step", "interval_s" => 0.0, "flood_depth" => 0.1, "n_states" => 250,
             "t_last" => 1800.0, "flooda5_version" => "0.8.13-fixture"))
    return path
end

function _depth_only_run(path)
    polys = Any[_square(0, 0), _square(1, 0), _square(0, 1), _square(1, 1; closed = true)]
    n = length(polys)
    cells = (ids = ["d$i" for i in 1:n], center_lons = [172.6005, 172.6015, 172.6005, 172.6015],
             center_lats = [-43.4995, -43.4995, -43.4985, -43.4985], elevations = [5.0, 5.5, 6.0, 6.5],
             static_vars = Dict("sgs_z_min" => [4.0, 4.0, 4.0, 4.0]),     # (not the datum: not an SGS run)
             id_scheme = "cartesian", polygons = polys)
    run = Dict{String,Any}("argv" => String[], "options" => Dict{String,Any}("flow_model" => "standard"),
                           "flow_model" => "standard")
    FO.create_output(path; flooda5_version = "0.8.10-fixture", cells, run)
    for (k, t) in enumerate([0.0, 30.0, 90.0])
        FO.append_frame!(path, t, ["water_depth" => [0.0, 0.05 * k, 0.11 * k, 0.3 * (k - 1)]])
    end
    return path
end

const _HZ = "landuse=woodland&mass=40&height=1.4&theta=70&alpha=3&beta=0.3&flood_depth=0.3"
const _HAZ = ("depth_velocity", "hazard_rating", "hazard_degree", "hazard_class_arr")
const _STATS = ("max_water_depth", "max_velocity", "max_depth_velocity", "max_hazard_rating",
                "max_hazard_degree", "max_hazard_class_arr", "time_max_water_depth", "time_first_flooded", "flood_duration")

function _sgs_requests()
    r = ["/api/info", "/api/frames", "/api/mesh", "/api/run"]
    for k in 1:7, v in ("water_depth", "velocity", "vel_u", "vel_v", "saturation")
        push!(r, "/api/frames/$k/$v")
    end
    for k in (2, 4, 5, 7), v in _HAZ
        push!(r, "/api/frames/$k/$v"); push!(r, "/api/frames/$k/$v?$_HZ"); push!(r, "/api/frames/$k/$v?landuse=pasture")
    end
    append!(r, ["/api/frames/0/water_depth", "/api/frames/8/water_depth", "/api/frames/x/water_depth",
                "/api/frames/3/volume", "/api/frames/3/nonesuch", "/api/frames/3/hazard_rating?landuse=bogus&mass=abc"])
    append!(r, ["/api/ranges", "/api/ranges?wet=0.05", "/api/ranges?wet=0&$_HZ", "/api/ranges?wet=abc",
                "/api/ranges?wet=-1", "/api/ranges?wet=0.5&landuse=none"])
    for s in _STATS
        push!(r, "/api/stats/$s"); push!(r, "/api/stats/$s?$_HZ")
    end
    push!(r, "/api/stats/nonesuch")
    for i in 0:5
        push!(r, "/api/cell/$i?frame=4&$_HZ")
    end
    append!(r, ["/api/cell/2", "/api/cell/1?frame=99", "/api/cell/6", "/api/cell/-1", "/api/cell/x"])
    for i in (0, 2, 4), v in ("water_depth", "velocity", "saturation", _HAZ...)
        push!(r, "/api/series/$i?var=$v&$_HZ")
    end
    append!(r, ["/api/series/3", "/api/series/1?var=volume", "/api/series/6", "/api/nonesuch"])
    return r
end

_depth_requests() = ["/api/info", "/api/frames", "/api/mesh", "/api/run", "/api/frames/1/water_depth",
                     "/api/frames/3/water_depth", "/api/frames/1/hazard_rating", "/api/ranges", "/api/ranges?wet=0",
                     "/api/stats/max_water_depth", "/api/stats/max_velocity", "/api/stats/time_first_flooded",
                     "/api/stats/time_max_water_depth",
                     "/api/stats/flood_duration?flood_depth=0.05", "/api/cell/3?frame=2", "/api/series/3",
                     "/api/series/1?var=depth_velocity"]

# ---------------------------------------------------------------------------
# Recording
# ---------------------------------------------------------------------------

# numbers JSON cannot hold, as the JavaScript reader decodes them
_num(x::AbstractFloat) = isfinite(x) ? x : string(x)            # "NaN", "Inf", "-Inf"
_num(x) = x
_enc(v::AbstractArray{<:AbstractFloat}) = Any[_num(x) for x in v]
_enc(v::AbstractArray) = collect(v)
_enc(v::AbstractFloat) = _num(v)
_enc(v) = v
_attr(v::AbstractFloat) = isfinite(v) ? v : Dict("float" => string(v))
_attr(v::AbstractArray{<:AbstractFloat}) = Any[_attr(x) for x in v]
_attr(v) = v

"Every group, dataset and attribute of `f`: path → {kind, keys | shape + dtype + data, attrs} (C order)."
function dump_h5(f)
    objs = Dict{String,Any}()
    attrs(o) = Dict{String,Any}(k => _attr(read(attributes(o)[k])) for k in keys(attributes(o)))
    function walk(g, path)
        objs[path] = Dict{String,Any}("kind" => "group", "keys" => sort!(collect(keys(g))), "attrs" => attrs(g))
        for k in keys(g)
            o = g[k]; p = path == "/" ? "/$k" : "$path/$k"
            if o isa HDF5.Group
                walk(o, p)
            else
                v = read(o)
                arr = v isa AbstractArray
                objs[p] = Dict{String,Any}("kind" => "dataset", "shape" => arr ? reverse(collect(size(v))) : Int[],
                                           "dtype" => string(eltype(v)), "data" => _enc(arr ? vec(v) : v),
                                           "attrs" => attrs(o))
            end
        end
    end
    walk(f, "/")
    return objs
end

function _record(port, path)
    r = HTTP.get("http://127.0.0.1:$port$path"; status_exception = false, retry = false)
    ct = HTTP.header(r, "Content-Type", "")
    d = Dict{String,Any}("path" => path, "status" => r.status)
    if r.status != 200
        d["text"] = String(r.body)
    elseif startswith(ct, "application/json")
        d["json"] = JSON3.read(String(r.body))
    else
        d["base64"] = base64encode(r.body)
    end
    return d
end

_free_port() = (s = HTTP.Sockets.listen(HTTP.Sockets.IPv4("127.0.0.1"), 0);
                p = Int(HTTP.Sockets.getsockname(s)[2]); close(s); p)

"Julia's frame picks for run_ranges (`unique(round.(Int, range(1, n; length = m)))`), for the JavaScript port."
_picks() = [Dict("n" => n, "m" => m, "picks" => (n <= m ? collect(1:n) : unique(round.(Int, range(1, n; length = m)))))
            for (n, m) in ((5, 24), (24, 24), (25, 24), (112, 24), (1000, 24), (31, 25), (100, 25), (51, 5), (7, 3))]

"""
    make_local_fixture(dir = mktempdir()) → Dict

Build the runs in `dir`, serve each, record the contents and the answers.
"""
function make_local_fixture(dir::AbstractString = mktempdir())
    runs = Any[]
    for (name, build, reqs) in (("sgs_fixture.h5", _sgs_run, _sgs_requests()),
                                ("depth_fixture.h5", _depth_only_run, _depth_requests()))
        p = build(joinpath(dir, name))
        port = _free_port()
        vis, f = with_logger(NullLogger()) do
            web_replay_server(p; port)
        end
        try
            answers = [_record(port, q) for q in reqs]
            contents = h5open(dump_h5, p)
            push!(runs, Dict{String,Any}("name" => name, "h5" => contents, "requests" => answers))
        finally
            VisualisationServer.stop(vis); close(f)
        end
    end
    return Dict{String,Any}("_comment" => "Generated by test/vis_local_fixture.jl (--write): do not edit.",
                            "runs" => runs, "picks" => _picks())
end

write_local_fixture(path = LOCAL_FIXTURE; dir = mktempdir()) =
    (write(path, JSON3.write(make_local_fixture(dir)) * "\n"); path)

if abspath(PROGRAM_FILE) == @__FILE__
    "--write" in ARGS || (println("usage: julia --project=ci test/vis_local_fixture.jl --write"); exit(1))
    println("wrote ", write_local_fixture())
end
