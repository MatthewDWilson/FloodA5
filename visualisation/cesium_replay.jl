# visualisation/cesium_replay.jl — replay an output file in the web globe viewer
#
#   julia --project=. visualisation/cesium_replay.jl RUN.h5 [options]
#   julia --project=. FloodModel.jl --replay RUN.h5 --vis cesium [--vis-port PORT]
#
# Serves a recorded run (FloodA5 output, format 2: the mesh is in the file) to
# the CesiumJS viewer (visualisation/cesium/) and waits; open the printed URL.
# Frames are read from the file on demand, so files of any length work, and the
# file can still be growing (frames written after start-up are not picked up:
# restart to see them).
#
# Options:
#   --port N             HTTP port (default 8080; --vis-port is accepted too)
#   --mesh FILE.parquet  the run's GeoParquet mesh (format-1 files only)
#   --offset M           default vertical offset in metres (model heights →
#                        ellipsoidal heights); the viewer remembers its own
#   --open               open the viewer in the default browser
#
# Cesium ion: set FLOODA5_CESIUM_ION_TOKEN for Google Photorealistic 3D Tiles
# and Cesium World Terrain (docs/USER_GUIDE.md §6.2).

isdefined(Main, :FloodA5Output) || include(joinpath(@__DIR__, "..", "io", "FloodA5Output.jl"))
isdefined(Main, :FloodA5Analysis) || include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
isdefined(Main, :VisualisationServer) || include(joinpath(@__DIR__, "VisualisationServer.jl"))

using HDF5

const _CESIUM_REPLAY_VARS = ("water_depth", "velocity", "vel_u", "vel_v", "saturation")

"""
    cesium_replay_server(h5; port = 8080, parquet = nothing, offset = nothing)
        → (server, file)

Start a `VisualisationServer` serving the run recorded in `h5`. The HDF5 file
stays open for reading (close it after `VisualisationServer.stop`).
"""
function cesium_replay_server(h5::AbstractString; port::Integer = 8080, parquet = nothing,
                              offset = nothing)
    VS = VisualisationServer; FO = FloodA5Output; FA = FloodA5Analysis
    run = FA.load_run(h5; parquet)
    m = run.mesh
    isempty(m.polygons) &&
        error("$h5 has no cell polygons: a format-1 file needs --mesh FILE.parquet")
    isempty(run.t) && error("$h5 has no frames yet")
    base = VS.depth_datum(run.flow_model, m.elevation, m.static)
    payload = VS.encode_mesh(m.polygons, m.ids, m.elevation, base; lon = m.lon, lat = m.lat)
    f = h5open(h5, "r")
    vars = [v for v in _CESIUM_REPLAY_VARS if FO.has_var(f, v)]
    "water_depth" in vars || (close(f); error("$h5 has no water_depth frames"))
    iolock = ReentrantLock()                       # HDF5.jl is not thread-safe
    readvar = (k, v) -> lock(() -> Float32.(FO.read_frame(f, v, k)), iolock)
    ri = FO.run_info(h5)
    info = Dict{String,Any}("label" => run.label, "flow_model" => run.flow_model,
                            "source" => abspath(h5), "id_scheme" => m.id_scheme,
                            "flooda5_version" => string(get(ri, "flooda5_version", "")))
    config = offset === nothing ? Dict{String,Any}() : Dict{String,Any}("vertical_offset" => Float64(offset))
    vis = try
        VS.start(; port = Int(port), info, config)
    catch
        close(f); rethrow()
    end
    VS.set_mesh!(vis, payload)
    VS.set_replay!(vis, run.t, vars, readvar)
    return vis, f
end

function _cesium_replay_usage(code::Int)
    println("""
Usage: julia --project=. visualisation/cesium_replay.jl RUN.h5 [options]

  --port N             HTTP port (default 8080)
  --mesh FILE.parquet  the run's GeoParquet mesh (format-1 files only)
  --offset M           default vertical offset (m) to ellipsoidal heights
  --open               open the viewer in the default browser
""")
    code == 0 || exit(code)
    return nothing
end

function _open_browser(url)
    try
        Sys.iswindows() ? run(`cmd /c start "" $url`; wait = false) :
        Sys.isapple()   ? run(`open $url`; wait = false) : run(`xdg-open $url`; wait = false)
    catch e
        @warn "Could not open a browser" exception = e
    end
end

"""
    cesium_replay_main(args) → nothing

The command line of visualisation/cesium_replay.jl (also `FloodModel.jl
--replay FILE --vis cesium`). Blocks until interrupted.
"""
function cesium_replay_main(args::Vector{String})
    ("--help" in args || "-h" in args) && return _cesium_replay_usage(0)
    args = copy(args)
    opt(flag) = (i = findfirst(==(flag), args);
                 i === nothing ? nothing :
                 (i == length(args) ? error("$flag needs a value") :
                  (v = args[i+1]; deleteat!(args, [i, i+1]); v)))
    flag(fl) = (i = findfirst(==(fl), args); i === nothing ? false : (deleteat!(args, i); true))
    port = parse(Int, something(opt("--port"), opt("--vis-port"), "8080"))
    mesh = opt("--mesh")
    off = opt("--offset")
    open_it = flag("--open")
    files = filter(a -> !startswith(a, "--"), args)
    length(files) == 1 || (println("ERROR: give one output file"); return _cesium_replay_usage(1))
    length(args) == 1 || error("unknown option(s): $(join(filter(startswith("--"), args), " "))")
    vis, f = cesium_replay_server(files[1]; port, parquet = mesh,
                                  offset = off === nothing ? nothing : parse(Float64, off))
    url = "http://localhost:$port/"
    @info "Replay: open $url in a browser (Ctrl-C to stop)"
    open_it && _open_browser(url)
    try
        wait(vis.shutdown)
    catch e
        e isa InterruptException || rethrow()
    finally
        VisualisationServer.stop(vis)
        close(f)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    cesium_replay_main(copy(ARGS))
end
