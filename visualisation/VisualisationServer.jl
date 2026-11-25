# VisualisationServer.jl
# ----------------------
# HTTP + WebSocket server for the web globe viewers — visualisation/cesium/
# (CesiumJS) and visualisation/giro3d/ (Giro3D), which share
# visualisation/viewer-common/ — for a running simulation (`--vis cesium`,
# `--vis giro3d`) or a recorded one (replay of an HDF5 output,
# visualisation/web_replay.jl). Plain HTTP.jl, no framework.
#
# The protocol (version 2, docs/DATA_FORMATS.md §6) is independent of the
# rendering engine: a binary mesh once, then per-frame Float32 arrays fetched
# on demand, one variable at a time, plus a WebSocket for notifications.
#
#   GET /                          → redirect to /viz/<viewer>/index.html
#   GET /viz/<dir>/<path>          → static files of the viewers (<dir>: cesium,
#                                    giro3d, viewer-common; no "..", no absolute paths)
#   GET /api/config                → viewer configuration (ion token, Cesium URL, …)
#   GET /api/info                  → run description, protocol version, state
#   GET /api/mesh                  → binary mesh (see `encode_mesh`); 204 until set
#   GET /api/frames                → {complete, frames: [{id, t}]}
#   GET /api/frames/{id}/{var}     → raw little-endian Float32[n_cells]; {var} may be a
#                                    hazard variable (HAZARD_VARS), with its parameters
#                                    in the query (`hazard_params`: landuse, mass, …)
#   GET /api/ranges?wet=0.05       → each variable's 2nd/98th percentiles over the run (`run_ranges`;
#                                    hazard variables too, with the hazard parameters)
#   GET /api/stats/{name}?…        → a statistic over the run (STAT_VARS: max_water_depth,
#                                    time_first_flooded, …), Float32[n_cells] (`run_stats`)
#   GET /api/cell/{i}?frame={id}   → one cell (0-based index): id, centre, values
#   GET /api/series/{i}?var=V      → one cell's values of V over the frames kept (`series_data`)
#   GET /api/run                   → the run record (`set_run!`; {} until set)
#   GET /api/about                 → variable descriptions (io/variables.toml)
#
# The viewer can also read an output file from the user's computer itself
# (viewer-common/js/local/: h5wasm in a worker answers the same /api), and
# runs without this server from a static web host (defaults.json, about.json).
#   WS  /api/live                  → hello, mesh, frame, reindex, complete (JSON)
#
# Every request must carry a Host header naming localhost / 127.0.0.1 (a
# guard against DNS rebinding: /api/config can hold a Cesium ion token), and
# the server listens on 127.0.0.1 only.
#
# Live frames are thinned by wall-clock time (`wants_frame`, `min_interval`)
# and kept within a memory budget: when it is exceeded every second frame is
# dropped (the first and last are kept), so a long run stays viewable end to
# end at a coarser time resolution. Frames have stable ids; clients refresh
# the list on "reindex".
#
# Usage (live):
#   vis = VisualisationServer.start(port = 8080, info = Dict("label" => "run"))
#   VisualisationServer.set_mesh!(vis, VisualisationServer.encode_mesh(rings, ids, elev, base))
#   VisualisationServer.wants_frame(vis) &&
#       VisualisationServer.push_frame!(vis, t, Dict("water_depth" => h, ...))
#   VisualisationServer.notify_complete!(vis)
#   VisualisationServer.stop(vis)
#
# Usage (replay): set_replay!(vis, times, vars, (k, var) -> Vector{Float32};
#                             series = (var, cell) -> Vector{Float32}).

module VisualisationServer

using HTTP
using JSON3
using Sockets
import TOML

# Access keys (Cesium ion token, Google key): environment, then key files
include(joinpath(@__DIR__, "..", "io", "FloodA5Keys.jl"))
using .FloodA5Keys: flooda5_key, key_source

# Hazard variables and statistics over a run: the same pure functions as
# FloodA5Analysis (FA.hazard_rating, FA.run_statistics, …)
include(joinpath(@__DIR__, "..", "analysis", "hazard.jl"))

# Call everything qualified (VisualisationServer.start, …): FloodModel.jl
# also loads MakieVisualiser, which has its own start/stop/push_frame!.
export VisServer, MeshPayload

const PROTOCOL = 2
const MESH_MAGIC = b"FA5M"
const MESH_VERSION = UInt32(1)
# Defaults of /api/config, shared with the viewer without a server (viewer-common/defaults.json)
const _DEFAULTS = Dict{String,Any}(String(k) => v for (k, v) in
    JSON3.read(read(joinpath(@__DIR__, "viewer-common", "defaults.json"), String)) if k != :_comment)
const DEFAULT_CESIUM_URL = String(_DEFAULTS["cesium_url"])
# Chart.js (the time-series plots), an ES module loaded when the first plot opens
const DEFAULT_CHARTJS_URL = String(_DEFAULTS["chartjs_url"])
"The h5wasm ES module (HDF5 in WebAssembly) the viewer reads local files with (0.8.10)."
const DEFAULT_H5WASM_URL = String(_DEFAULTS["h5wasm_url"])

# ---------------------------------------------------------------------------
# Mesh payload
# ---------------------------------------------------------------------------

"""
A mesh encoded for the viewer, plus the per-cell values `/api/cell` reports.

Binary layout (little-endian; every section 4-byte aligned):

    0   "FA5M"                 magic
    4   UInt32  version (1)
    8   UInt32  n_cells
    12  UInt32  n_vertices      (ring vertices; the closing repeat is dropped)
    16  Float64 origin_lon      (degrees; the mesh's bounding-box centre)
    24  Float64 origin_lat
    32  UInt32  offsets[n_cells + 1]     cell i (0-based) has vertices offsets[i] … offsets[i+1]-1
        Float32 dlon[n_vertices]         vertex lon − origin_lon (degrees)
        Float32 dlat[n_vertices]         vertex lat − origin_lat (degrees)
        Float32 elevation[n_cells]       cell bed elevation (m; NaN where unknown)
        Float32 base[n_cells]            depth datum: water level = base + depth (m)

Offsets from the origin in Float32 keep vertices to about 1 mm over a few
tens of kilometres.
"""
struct MeshPayload
    bytes     :: Vector{UInt8}
    n_cells   :: Int
    ids       :: Vector{String}
    lon       :: Vector{Float64}     # cell centres
    lat       :: Vector{Float64}
    elevation :: Vector{Float32}
    base      :: Vector{Float32}
end

# A cell ring as (lon, lat) vectors, without the closing repeat.
function _ring(b::AbstractVector{<:AbstractVector})
    return _open_ring(Float64[p[1] for p in b], Float64[p[2] for p in b])
end
_ring(m::AbstractMatrix) = size(m, 1) == 2 ? _open_ring(Float64.(m[1, :]), Float64.(m[2, :])) :
                                             _open_ring(Float64.(m[:, 1]), Float64.(m[:, 2]))
function _open_ring(lon::Vector{Float64}, lat::Vector{Float64})
    n = length(lon)
    (n > 1 && lon[1] == lon[n] && lat[1] == lat[n]) && (n -= 1)
    return lon[1:n], lat[1:n]
end

"""
    encode_mesh(rings, ids, elevation, base; lon = nothing, lat = nothing) → MeshPayload

`rings[i]` is cell i's boundary: a vector of `[lon, lat]` points (an A5Mesh
cell's `boundary`) or a 2×k matrix (`FA.MeshGeom.polygons`). `base` is the
depth datum (see `depth_datum`). Centres default to the ring vertex means.
"""
function encode_mesh(rings::AbstractVector, ids::AbstractVector, elevation::AbstractVector,
                     base::AbstractVector; lon = nothing, lat = nothing)
    n = length(rings)
    length(ids) == n == length(elevation) == length(base) ||
        error("encode_mesh: rings, ids, elevation and base must have the same length")
    n > 0 || error("encode_mesh: empty mesh")
    rs = [_ring(r) for r in rings]
    any(r -> length(r[1]) < 3, rs) && error("encode_mesh: a cell has fewer than 3 vertices")
    nv = sum(r -> length(r[1]), rs)
    lo_min = minimum(r -> minimum(r[1]), rs); lo_max = maximum(r -> maximum(r[1]), rs)
    la_min = minimum(r -> minimum(r[2]), rs); la_max = maximum(r -> maximum(r[2]), rs)
    olon = (lo_min + lo_max) / 2; olat = (la_min + la_max) / 2
    io = IOBuffer()
    write(io, MESH_MAGIC, htol(MESH_VERSION), htol(UInt32(n)), htol(UInt32(nv)),
          htol(Float64(olon)), htol(Float64(olat)))
    off = UInt32(0)
    write(io, htol(off))
    for r in rs
        off += UInt32(length(r[1])); write(io, htol(off))
    end
    for r in rs, x in r[1]
        write(io, htol(Float32(x - olon)))
    end
    for r in rs, y in r[2]
        write(io, htol(Float32(y - olat)))
    end
    elev32 = Float32.(elevation); base32 = Float32.(base)
    for v in elev32; write(io, htol(v)); end
    for v in base32; write(io, htol(v)); end
    clon = lon === nothing ? [sum(r[1]) / length(r[1]) for r in rs] : Float64.(lon)
    clat = lat === nothing ? [sum(r[2]) / length(r[2]) for r in rs] : Float64.(lat)
    return MeshPayload(take!(io), n, String.(ids), clon, clat, elev32, base32)
end

"""
    depth_datum(flow_model, elevation, static) → Vector{Float64}

The level the model's water depth is measured from: SGS depth is above the
cell's lowest sub-grid elevation (`sgs_z_min`), the other models' above the
cell elevation. Water level = datum + depth.
"""
function depth_datum(flow_model, elevation::AbstractVector, static::AbstractDict)
    if string(flow_model) == "sgs" && haskey(static, "sgs_z_min")
        z = Float64.(static["sgs_z_min"])
        return [isfinite(z[i]) ? z[i] : Float64(elevation[i]) for i in eachindex(z)]
    end
    return Float64.(elevation)
end

# ---------------------------------------------------------------------------
# Frame stores
# ---------------------------------------------------------------------------

mutable struct LiveFrames
    ids    :: Vector{Int}
    t      :: Vector{Float64}
    data   :: Vector{Dict{String,Vector{Float32}}}
    bytes  :: Int
    budget :: Int          # bytes
    next   :: Int
end
LiveFrames(budget) = LiveFrames(Int[], Float64[], Dict{String,Vector{Float32}}[], 0, budget, 1)

struct ReplayFrames
    t      :: Vector{Float64}
    vars   :: Vector{String}
    read   :: Function                  # (k, var) → Vector{Float32}
    series :: Union{Function,Nothing}   # (var, cell) → Vector{Float32} over all frames, or nothing
end
ReplayFrames(t, vars, read) = ReplayFrames(t, vars, read, nothing)

_frame_bytes(d) = sum(v -> sizeof(v), values(d); init = 0)

"Drop every second frame (keeping the first and the last) until within budget."
function _decimate!(f::LiveFrames)
    changed = false
    while f.bytes > f.budget && length(f.ids) > 2
        keep = [i for i in eachindex(f.ids) if isodd(i) || i == length(f.ids)]
        f.ids = f.ids[keep]; f.t = f.t[keep]; f.data = f.data[keep]
        f.bytes = sum(_frame_bytes, f.data; init = 0)
        changed = true
    end
    return changed
end

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

mutable struct VisServer
    host          :: String
    port          :: Int
    viz_dir       :: String
    lock          :: ReentrantLock
    mesh          :: Union{MeshPayload,Nothing}
    frames        :: Union{LiveFrames,ReplayFrames}
    complete      :: Bool
    info          :: Dict{String,Any}
    config        :: Dict{String,Any}
    clients       :: Vector{Channel{String}}
    min_interval  :: Float64
    last_push     :: Float64
    cache         :: Dict{Tuple{Int,String},Vector{Float32}}   # replay reads
    cache_order   :: Vector{Tuple{Int,String}}
    cache_size    :: Int
    server        :: Any
    task          :: Union{Task,Nothing}
    shutdown      :: Base.Event
    ranges        :: Dict{Any,Any}          # run_ranges, by (wet, frames pushed, limits, hazard parameters)
    run           :: Dict{String,Any}       # /api/run (`set_run!`)
    stats         :: Dict{Any,Any}          # run_stats: (hazard parameters, frames pushed) → Task
end

_envnum(name, default) = something(tryparse(Float64, get(ENV, name, "")), default)

"The web viewers (directories of `visualisation/`), and the directory they share."
const VIEWERS = ("cesium", "giro3d")
const _STATIC_DIRS = (VIEWERS..., "viewer-common")

"""
    start(; port = 8080, host = "127.0.0.1", viewer = "cesium", viz_dir, info = Dict(),
          config = Dict(), min_interval = 1.0, frame_budget_mb = 1024) → VisServer

Start the server without blocking. `viewer` (`cesium` or `giro3d`) is the one
`/` opens; both are served from `viz_dir` (`visualisation/`). `min_interval` (s, wall clock) thins live
frames (`wants_frame`); `frame_budget_mb` bounds the live frames kept in
memory. Their defaults can be set with the environment variables
`FLOODA5_VIS_INTERVAL` and `FLOODA5_VIS_BUDGET_MB`. `info` describes the run
(`label`, `flow_model`, …) for `/api/info`;
`config` adds viewer settings to `/api/config` (see `viewer_config`). The
handlers run on the interactive thread pool when there is one (Julia 1.12
starts one by default), so the viewer stays responsive while the time loop
keeps the default threads busy.
"""
function start(; port::Int = 8080, host::AbstractString = "127.0.0.1",
               viewer::AbstractString = "cesium",
               viz_dir::AbstractString = @__DIR__,
               info::AbstractDict = Dict{String,Any}(), config::AbstractDict = Dict{String,Any}(),
               min_interval::Real = _envnum("FLOODA5_VIS_INTERVAL", 1.0),
               frame_budget_mb::Real = _envnum("FLOODA5_VIS_BUDGET_MB", 1024))
    viewer in VIEWERS || error("VisServer: viewer must be one of $(join(VIEWERS, ", ")) (got $viewer)")
    s = VisServer(String(host), port, String(viz_dir), ReentrantLock(), nothing,
                  LiveFrames(round(Int, frame_budget_mb * 2^20)), false,
                  merge(Dict{String,Any}(String(k) => v for (k, v) in info), Dict{String,Any}("viewer" => viewer)),
                  Dict{String,Any}(String(k) => v for (k, v) in config),
                  Channel{String}[], Float64(min_interval), -Inf,
                  Dict{Tuple{Int,String},Vector{Float32}}(), Tuple{Int,String}[], 48,
                  nothing, nothing, Base.Event(), Dict{Any,Any}(),
                  Dict{String,Any}(), Dict{Any,Any}())
    _port_free(s.host, port) ||
        error("VisServer: port $port is already in use. Stop the other process or use --vis-port.")
    ready = Channel{Any}(1)
    body = () -> begin
        try
            srv = HTTP.serve!(http -> _handle_stream(s, http), s.host, port;
                              stream = true, verbose = false, access_log = nothing)
            put!(ready, srv)
            wait(srv)
        catch e
            isready(ready) || put!(ready, e)
            e isa InterruptException || e isa Base.IOError || @warn "VisServer stopped" exception = e
        finally
            notify(s.shutdown)
        end
    end
    s.task = Threads.nthreads(:interactive) > 0 ? Threads.@spawn(:interactive, body()) :
                                                  Threads.@spawn(body())
    r = take!(ready)
    r isa Exception && throw(r)
    s.server = r
    @info "VisServer: listening on http://localhost:$port"
    src = key_source("FLOODA5_CESIUM_ION_TOKEN")
    src === nothing ?
        @info("VisServer: no Cesium ion token (FLOODA5_CESIUM_ION_TOKEN in the environment, keys.json " *
              "or ~/.flooda5/keys.json; template keys.example.json): Google 3D tiles unavailable") :
        @info("VisServer: Cesium ion token from $src")
    return s
end

function _port_free(host, port)
    try
        close(Sockets.listen(Sockets.IPv4(host), port)); return true
    catch
        return false
    end
end

"Stop the server and close client connections."
function stop(s::VisServer)
    lock(s.lock) do
        foreach(close, s.clients); empty!(s.clients)
    end
    s.server === nothing || try close(s.server) catch end
    notify(s.shutdown)
    return nothing
end

"""
    set_run!(s, record) → s

The run record `/api/run` serves (the viewer's *Run* tab):
`FloodA5Output.run_record` of the output file (replay) or of the run in
progress. Plain JSON-able data (`run_record` makes it so).
"""
function set_run!(s::VisServer, record::AbstractDict)
    lock(s.lock) do
        s.run = Dict{String,Any}(String(k) => v for (k, v) in record)
    end
    return s
end

const _VARIABLES_TOML = joinpath(@__DIR__, "..", "io", "variables.toml")
const _ABOUT = Ref{Any}(nothing)

"The variable descriptions (io/variables.toml; FloodA5Output.variable_docs reads the same file)."
function about()
    a = _ABOUT[]
    a === nothing || return a
    return (_ABOUT[] = TOML.parsefile(_VARIABLES_TOML))
end

"Set (or replace) the mesh and tell connected clients."
function set_mesh!(s::VisServer, m::MeshPayload)
    lock(s.lock) do
        s.mesh = m
        empty!(s.ranges); empty!(s.stats)
        s.info["n_cells"] = m.n_cells
        _broadcast!(s, (type = "mesh", n_cells = m.n_cells))
    end
    return s
end

"""
    set_replay!(s, times, vars, read; series = nothing) → s

Serve a recorded run: `read(k, var)` returns frame k's `var` (in mesh order).
`series(var, cell)` (optional; cell 1-based, mesh order) returns one cell's
values over all frames, for `/api/series` — without it a series is read
frame by frame. The run is complete from the start.
"""
function set_replay!(s::VisServer, times::AbstractVector, vars::AbstractVector, read::Function;
                     series::Union{Function,Nothing} = nothing)
    lock(s.lock) do
        s.frames = ReplayFrames(Float64.(times), String.(vars), read, series)
        s.complete = true
        empty!(s.cache); empty!(s.cache_order); empty!(s.ranges); empty!(s.stats)
        _broadcast!(s, (type = "reindex",))
        _broadcast!(s, (type = "complete", frames = length(times)))
    end
    return s
end

"""
    wants_frame(s) → Bool

Whether a live frame pushed now would be kept: the first frame, then one per
`min_interval` seconds of wall-clock time. Check it before building the
arrays for `push_frame!`.
"""
wants_frame(s::VisServer) = s.frames isa LiveFrames &&
    (isempty(s.frames.ids) || time() - s.last_push >= s.min_interval)

"""
    push_frame!(s, t, vars) → frame id

Store a live frame (`vars`: name → values in mesh order; copied as Float32)
and notify clients.
"""
function push_frame!(s::VisServer, t::Real, vars::AbstractDict)
    d = Dict{String,Vector{Float32}}(String(k) => Vector{Float32}(v) for (k, v) in vars)
    return lock(s.lock) do
        f = s.frames
        f isa LiveFrames || error("push_frame!: the server is serving a replay")
        id = f.next; f.next += 1
        push!(f.ids, id); push!(f.t, Float64(t)); push!(f.data, d)
        f.bytes += _frame_bytes(d)
        s.last_push = time()
        s.info["vars"] = sort!(collect(keys(d)))
        _decimate!(f) && _broadcast!(s, (type = "reindex",))
        _broadcast!(s, (type = "frame", id = id, t = Float64(t), frames = length(f.ids)))
        id
    end
end

"Mark the run complete (no more frames) and tell clients."
function notify_complete!(s::VisServer)
    n = lock(s.lock) do
        s.complete = true
        _broadcast!(s, (type = "complete", frames = _nframes(s.frames)))
        _nframes(s.frames)
    end
    @info "VisServer: run complete ($n frames kept)"
    return nothing
end

_nframes(f::LiveFrames) = length(f.ids)
_nframes(f::ReplayFrames) = length(f.t)
_frame_list(f::LiveFrames) = [(id = f.ids[i], t = f.t[i]) for i in eachindex(f.ids)]
_frame_list(f::ReplayFrames) = [(id = k, t = f.t[k]) for k in eachindex(f.t)]
_vars(s) = s.frames isa ReplayFrames ? s.frames.vars : get(s.info, "vars", String[])

# ---------------------------------------------------------------------------
# Derived (hazard) variables
# ---------------------------------------------------------------------------

"The frame variables the hazard variables are derived from."
const DERIVED_INPUTS = ("water_depth", "velocity")

"Whether `var` is a hazard variable derived from depth and velocity (`HAZARD_VARS`)."
is_derived(var::AbstractString) = haskey(HAZARD_VARS, var)

"""
    derived_values(var, get, p::HazardParams) → Vector{Float32} or nothing

Hazard variable `var` from its inputs: `get(name)` returns the input variable
"water_depth" or "velocity" (or nothing). The inputs may be of any length —
one frame's cells, or one cell's values over the frames — so a time series
derives the same way as a frame.
"""
function derived_values(var::AbstractString, get::Function, p::HazardParams)
    d = get("water_depth"); d === nothing && return nothing
    v = get("velocity");    v === nothing && return nothing
    return derive_hazard(var, d, v, p)
end

"Whether the run has the inputs of the hazard variables."
_has_derived(s::VisServer) = all(in(_vars(s)), DERIVED_INPUTS)

"Names of the hazard variables and statistics the run supports (for /api/info)."
_derived_names(s::VisServer) = _has_derived(s) ? sort!(collect(keys(HAZARD_VARS))) : String[]
_stat_names(s::VisServer) = _has_derived(s) ? collect(STAT_VARS) :
    ("water_depth" in _vars(s) ? ["max_water_depth", "time_first_flooded", "flood_duration"] : String[])

# Frame data (Float32, mesh order), or nothing. Hazard variables are derived
# (not cached: elementwise, cheap) from the cached inputs.
function _frame_data(s::VisServer, id::Int, var::String; params::HazardParams = HazardParams())
    is_derived(var) && return derived_values(var, v -> _frame_data(s, id, v), params)
    f = s.frames
    if f isa LiveFrames
        return lock(s.lock) do
            i = searchsortedfirst(f.ids, id)
            (i <= length(f.ids) && f.ids[i] == id) ? get(f.data[i], var, nothing) : nothing
        end
    end
    (1 <= id <= length(f.t) && var in f.vars) || return nothing
    key = (id, var)
    v = lock(s.lock) do; get(s.cache, key, nothing); end
    v === nothing || return v
    v = Vector{Float32}(f.read(id, var))
    lock(s.lock) do
        if !haskey(s.cache, key)
            s.cache[key] = v; push!(s.cache_order, key)
            while length(s.cache_order) > s.cache_size
                delete!(s.cache, popfirst!(s.cache_order))
            end
        end
    end
    return v
end

"""
    series_data(s, cell, var; params = HazardParams()) → (ids, t, values) or nothing

One cell's values of frame variable `var` over the frames kept (`cell`:
1-based, mesh order): frame ids, times (s) and Float32 values. `nothing` if
the variable is not served. A replay with a `series` reader reads the cell's
row once; otherwise (live, or no reader) the values come from each frame.
Hazard variables are derived from the depth and velocity series with
`params` (`derived_values`, as for a frame).
"""
function series_data(s::VisServer, cell::Int, var::String; params::HazardParams = HazardParams())
    if is_derived(var)
        base = Dict{String,Any}()
        for v in DERIVED_INPUTS
            r = series_data(s, cell, String(v))
            r === nothing && return nothing
            base[String(v)] = r
        end
        ids, t, _ = base["water_depth"]
        vals = derived_values(var, v -> base[v][3], params)
        return vals === nothing ? nothing : (ids, t, Vector{Float32}(vals))
    end
    f = s.frames
    if f isa LiveFrames
        return lock(s.lock) do
            any(d -> haskey(d, var), f.data) || return nothing
            vals = Float32[(v = get(d, var, nothing); v === nothing ? NaN32 : v[cell]) for d in f.data]
            (copy(f.ids), copy(f.t), vals)
        end
    end
    var in f.vars || return nothing
    ids = collect(1:length(f.t))
    if f.series !== nothing
        return (ids, copy(f.t), Vector{Float32}(f.series(var, cell)))
    end
    vals = Float32[(d = _frame_data(s, k, var); d === nothing ? NaN32 : d[cell]) for k in ids]
    return (ids, copy(f.t), vals)
end

# ---------------------------------------------------------------------------
# Colour ranges over the run
# ---------------------------------------------------------------------------

# q-quantile (0–1) of sorted values, linear between order statistics (NaN if empty)
function _sorted_quantile(a::AbstractVector, q::Real)
    isempty(a) && return NaN
    h = (length(a) - 1) * q + 1
    lo = floor(Int, h); hi = min(lo + 1, length(a))
    return Float64(a[lo]) + (h - lo) * (Float64(a[hi]) - Float64(a[lo]))
end

"""
    run_ranges(s; wet = 0.05, max_frames = 24, max_values = 4_000_000) → NamedTuple

Colour ranges fixed for the whole run: for each frame variable, and `level`
(depth datum + depth), the 2nd and 98th percentiles (plus min, max and the
count) over the wet cells — water depth ≥ `wet` — of up to `max_frames`
frames spread evenly over the run (a live run: over its frames so far). Cells
are subsampled evenly when frames × cells exceeds `max_values`. Cached until
the mesh, the replay or the live frames change.
"""
function run_ranges(s::VisServer; wet::Real = 0.05, max_frames::Int = 24, max_values::Int = 4_000_000,
                    params::HazardParams = HazardParams())
    m = s.mesh
    m === nothing && return (wet = Float64(wet), frames = 0, of = 0, complete = false, vars = Dict{String,Any}())
    # what to scan (under the lock), read outside it
    picks, getter, stamp, n, complete = lock(s.lock) do
        f = s.frames
        n = _nframes(f)
        ks = n <= max_frames ? collect(1:n) : unique(round.(Int, range(1, n; length = max_frames)))
        if f isa LiveFrames
            ([f.data[k] for k in ks], (d, v) -> get(d, v, nothing), f.next, n, s.complete)
        else
            (ks, (k, v) -> v in f.vars ? f.read(k, v) : nothing, 0, n, true)
        end
    end
    key = (Float64(wet), stamp, max_frames, max_values, hazard_key(params))
    hit = lock(s.lock) do; get(s.ranges, key, nothing); end
    hit === nothing || return hit
    names = String[v for v in _vars(s)]
    derived = _derived_names(s)
    vals = Dict{String,Vector{Float32}}(v => Float32[] for v in [names; derived])
    vals["level"] = Float32[]
    stride = max(1, cld(m.n_cells * max(length(picks), 1), max_values))
    for p in picks
        d = getter(p, "water_depth")
        d === nothing && continue
        idx = [i for i in 1:stride:min(length(d), m.n_cells) if d[i] >= wet]
        vel = nothing
        for v in names
            x = v == "water_depth" ? d : getter(p, v)
            v == "velocity" && (vel = x)
            x === nothing || append!(vals[v], (x[i] for i in idx if isfinite(x[i])))
        end
        if !isempty(derived) && vel !== nothing
            dw = Float32[d[i] for i in idx]; vw = Float32[vel[i] for i in idx]
            for v in derived
                append!(vals[v], filter(isfinite, derive_hazard(v, dw, vw, params)))
            end
        end
        append!(vals["level"], (m.base[i] + d[i] for i in idx if isfinite(m.base[i] + d[i])))
    end
    fl(x) = isfinite(x) ? x : nothing           # (JSON has no NaN)
    out = Dict{String,Any}()
    for (v, a) in vals
        sort!(a)
        out[v] = (p02 = fl(_sorted_quantile(a, 0.02)), p98 = fl(_sorted_quantile(a, 0.98)),
                  min = isempty(a) ? nothing : fl(Float64(a[1])), max = isempty(a) ? nothing : fl(Float64(a[end])),
                  n = length(a))
    end
    r = (wet = Float64(wet), frames = length(picks), of = n, complete = complete, vars = out)
    lock(s.lock) do
        filter!(kv -> kv[1][2] == stamp, s.ranges)     # older live states are stale
        s.ranges[key] = r
    end
    return r
end

# ---------------------------------------------------------------------------
# Statistics over the run (static maps)
# ---------------------------------------------------------------------------

"""
    run_stats(s, p::HazardParams) → Dict{String,Vector{Float32}}

`run_statistics` (maxima of depth, velocity and the hazard variables, time of
first flooding, duration of flooding) over every frame of the replay, or over
the live frames kept so far — at the resolution of those frames. Computed
once per parameter set and set of frames (concurrent requests share the
work), off the server's thread.
"""
function run_stats(s::VisServer, p::HazardParams)
    m = s.mesh
    m === nothing && return Dict{String,Vector{Float32}}()
    t = lock(s.lock) do
        f = s.frames
        stamp = f isa LiveFrames ? f.next : 0
        key = (hazard_key(p), stamp)
        task = get(s.stats, key, nothing)
        if task === nothing
            if f isa LiveFrames
                times = copy(f.t)
                data = [(get(d, "water_depth", nothing), get(d, "velocity", nothing)) for d in f.data]
                frame = k -> data[k]
            else
                times = f.t; rd = f.read; hv = "velocity" in f.vars
                frame = k -> (rd(k, "water_depth"), hv ? rd(k, "velocity") : nothing)
            end
            n = m.n_cells
            task = Threads.@spawn run_statistics(times, frame, n, p)
            filter!(kv -> kv[1][2] == stamp, s.stats)          # older live states are stale
            length(s.stats) >= 4 && empty!(s.stats)              # a few parameter sets at most
            s.stats[key] = task
        end
        task
    end
    try
        return fetch(t)
    catch
        lock(s.lock) do; filter!(kv -> kv[2] !== t, s.stats); end
        rethrow()
    end
end

function _broadcast!(s::VisServer, msg)
    str = JSON3.write(msg)
    for ch in s.clients
        isopen(ch) && put!(ch, str)
    end
end

# ---------------------------------------------------------------------------
# Configuration for the viewer
# ---------------------------------------------------------------------------

"""
    viewer_config(s) → Dict

`/api/config`: `visualisation/viewer-common/config.json` (git-ignored) if present,
then the server's `config`, then the keys — `FLOODA5_CESIUM_ION_TOKEN`
(Cesium ion access token: Google Photorealistic 3D Tiles, Cesium World
Terrain) and `FLOODA5_GOOGLE_TILES_KEY` (Google Map Tiles API key, the other
route to the Google tiles) — from the environment or a key file
(`io/FloodA5Keys.jl`: keys.json in the repository root or
~/.flooda5/keys.json), then `FLOODA5_CESIUM_URL` (base URL of a CesiumJS
build, e.g. a local copy for offline use), `FLOODA5_CHARTJS_URL` (the
Chart.js ES module) and `FLOODA5_H5WASM_URL` (the h5wasm ES module, for
local files) from the environment. The defaults are viewer-common/defaults.json.
"""
function viewer_config(s::VisServer)
    c = merge(Dict{String,Any}(_DEFAULTS), Dict{String,Any}("protocol" => PROTOCOL, "ion_token" => "", "google_key" => ""))
    f = joinpath(s.viz_dir, "viewer-common", "config.json")
    if isfile(f)
        try
            for (k, v) in JSON3.read(read(f, String), Dict{String,Any})
                k == "cesium_ion_token" ? (c["ion_token"] = v) : (c[String(k)] = v)
            end
        catch e
            @warn "VisServer: cannot read $f" exception = e
        end
    end
    merge!(c, s.config)
    for (k, e) in (("ion_token", "FLOODA5_CESIUM_ION_TOKEN"), ("google_key", "FLOODA5_GOOGLE_TILES_KEY"))
        v = flooda5_key(e)
        isempty(v) || (c[k] = v)
    end
    v = get(ENV, "FLOODA5_CESIUM_URL", "")
    isempty(v) || (c["cesium_url"] = v)
    v = get(ENV, "FLOODA5_CHARTJS_URL", "")
    isempty(v) || (c["chartjs_url"] = v)
    v = get(ENV, "FLOODA5_H5WASM_URL", "")
    isempty(v) || (c["h5wasm_url"] = v)
    occursin("YOUR_", string(c["ion_token"])) && (c["ion_token"] = "")
    endswith(string(c["cesium_url"]), "/") || (c["cesium_url"] = string(c["cesium_url"], "/"))
    return c
end

# ---------------------------------------------------------------------------
# HTTP handling
# ---------------------------------------------------------------------------

const _MIME = Dict(".html" => "text/html; charset=utf-8", ".js" => "text/javascript; charset=utf-8",
                   ".mjs" => "text/javascript; charset=utf-8", ".css" => "text/css; charset=utf-8",
                   ".json" => "application/json", ".png" => "image/png", ".svg" => "image/svg+xml",
                   ".ico" => "image/x-icon", ".wasm" => "application/wasm",
                   ".glb" => "model/gltf-binary", ".txt" => "text/plain; charset=utf-8")

_json(x; status = 200) = HTTP.Response(status, ["Content-Type" => "application/json",
                                                "Cache-Control" => "no-store"], JSON3.write(x))
_text(status, msg) = HTTP.Response(status, ["Content-Type" => "text/plain; charset=utf-8"], msg)
_binary(bytes) = HTTP.Response(200, ["Content-Type" => "application/octet-stream",
                                     "Cache-Control" => "no-store"], bytes)

function _host_ok(req::HTTP.Request)
    h = lowercase(HTTP.header(req, "Host", ""))
    name = startswith(h, "[") ? h[1:something(findfirst(']', h), length(h))] :
           first(split(h, ':'))
    return name in ("localhost", "127.0.0.1", "[::1]")
end

function _handle_stream(s::VisServer, http::HTTP.Stream)
    req = http.message
    if !_host_ok(req)
        HTTP.setstatus(http, 403); HTTP.startwrite(http); write(http, "Forbidden host")
        return
    end
    if HTTP.WebSockets.isupgrade(req) && HTTP.URI(req.target).path == "/api/live"
        HTTP.WebSockets.upgrade(ws -> _serve_ws(s, ws), http)
        return
    end
    resp = try
        _route(s, req)
    catch e
        @warn "VisServer: error handling $(req.method) $(req.target)" exception = (e, catch_backtrace())
        _text(500, "Internal error")
    end
    HTTP.setstatus(http, resp.status)
    for (k, v) in resp.headers
        HTTP.setheader(http, k => v)
    end
    HTTP.setheader(http, "Content-Length" => string(length(resp.body)))
    HTTP.startwrite(http)
    req.method == "HEAD" || write(http, resp.body)
    return
end

function _serve_ws(s::VisServer, ws)
    ch = Channel{String}(Inf)
    lock(s.lock) do
        put!(ch, JSON3.write((type = "hello", protocol = PROTOCOL,
                              mode = s.frames isa ReplayFrames ? "replay" : "live",
                              complete = s.complete, frames = _nframes(s.frames),
                              mesh_ready = s.mesh !== nothing)))
        push!(s.clients, ch)
    end
    # The writer sends queued messages; when the queue is closed (`stop`, or
    # the client left) it closes the WebSocket, which ends the reader below
    # and the client's side.
    writer = @async begin
        try
            for msg in ch
                HTTP.WebSockets.send(ws, msg)
            end
        catch
        end
        try close(ws) catch end
    end
    try
        for _ in ws      # the client sends nothing we need; read until it closes
        end
    catch
    finally
        lock(s.lock) do; filter!(c -> c !== ch, s.clients); end
        close(ch)
        try wait(writer) catch end
    end
    return nothing
end

function _route(s::VisServer, req::HTTP.Request)
    req.method in ("GET", "HEAD") || return _text(405, "Method not allowed")
    uri = HTTP.URI(req.target)
    path = HTTP.unescapeuri(uri.path)
    parts = filter(!isempty, split(path, '/'))
    isempty(parts) && return HTTP.Response(302, ["Location" => "/viz/$(s.info["viewer"])/index.html" *
                                                 (isempty(uri.query) ? "" : "?" * uri.query)], "")
    if parts[1] == "viz"
        return _static(s, parts[2:end])
    elseif parts[1] == "api" && length(parts) >= 2
        a = parts[2]
        a == "config" && length(parts) == 2 && return _json(viewer_config(s))
        a == "info"   && length(parts) == 2 && return _json(_info(s))
        a == "run"    && length(parts) == 2 && return _json(lock(() -> copy(s.run), s.lock))
        a == "about"  && length(parts) == 2 && return _json(about())
        if a == "mesh" && length(parts) == 2
            m = s.mesh
            return m === nothing ? HTTP.Response(204, "") : _binary(m.bytes)
        end
        if a == "frames"
            length(parts) == 2 && return _json(lock(s.lock) do
                (complete = s.complete, frames = _frame_list(s.frames))
            end)
            if length(parts) == 4
                id = tryparse(Int, parts[3])
                id === nothing && return _text(400, "Bad frame id")
                v = _frame_data(s, id, String(parts[4]); params = hazard_params(HTTP.queryparams(uri)))
                v === nothing && return _text(404, "No frame $id variable $(parts[4])")
                return _binary(collect(reinterpret(UInt8, htol.(v))))
            end
        end
        if a == "ranges" && length(parts) == 2
            q = HTTP.queryparams(uri)
            w = something(tryparse(Float64, get(q, "wet", "")), 0.05)
            w = isfinite(w) ? max(w, 0.0) : 0.05
            # (scanning frames takes a moment: off the server's thread)
            return _json(fetch(Threads.@spawn run_ranges(s; wet = w, params = hazard_params(q))))
        end
        if a == "stats" && length(parts) == 3
            name = String(parts[3])
            name in STAT_VARS || return _text(404, "No statistic $name")
            st = run_stats(s, hazard_params(HTTP.queryparams(uri)))
            v = get(st, name, nothing)
            v === nothing && return _text(404, "No statistic $name for this run")
            return _binary(collect(reinterpret(UInt8, htol.(v))))
        end
        if a == "cell" && length(parts) == 3
            return _cell(s, tryparse(Int, parts[3]), HTTP.queryparams(uri))
        end
        if a == "series" && length(parts) == 3
            q = HTTP.queryparams(uri)
            return _series(s, tryparse(Int, parts[3]), get(q, "var", "water_depth"), hazard_params(q))
        end
    end
    return _text(404, "Not found")
end

function _info(s::VisServer)
    lock(s.lock) do
        d = Dict{String,Any}(s.info)
        d["protocol"] = PROTOCOL
        d["mode"] = s.frames isa ReplayFrames ? "replay" : "live"
        d["complete"] = s.complete
        d["frames"] = _nframes(s.frames)
        d["vars"] = _vars(s)
        d["derived"] = _derived_names(s)
        d["stats"] = _stat_names(s)
        d["mesh_ready"] = s.mesh !== nothing
        s.mesh === nothing || (d["n_cells"] = s.mesh.n_cells)
        return d
    end
end

function _cell(s::VisServer, i, q)
    m = s.mesh
    (m === nothing || i === nothing || !(0 <= i < m.n_cells)) && return _text(404, "No such cell")
    k = i + 1          # the protocol's cell indices are 0-based (mesh order)
    fl = (x -> isfinite(x) ? x : nothing)
    vals = Dict{String,Any}()
    fid = tryparse(Int, get(q, "frame", ""))
    if fid !== nothing
        for v in _vars(s)
            d = _frame_data(s, fid, v)
            d === nothing || (vals[v] = fl(d[k]))
        end
        p = hazard_params(q)
        for v in _derived_names(s)
            d = derived_values(v, n -> (x = _frame_data(s, fid, n); x === nothing ? nothing : x[k:k]), p)
            d === nothing || (vals[v] = fl(d[1]))
        end
    end
    return _json((index = i, id = m.ids[k], lon = m.lon[k], lat = m.lat[k],
                  elevation = fl(m.elevation[k]), base = fl(m.base[k]), frame = fid, values = vals))
end

function _series(s::VisServer, i, var, params::HazardParams = HazardParams())
    m = s.mesh
    (m === nothing || i === nothing || !(0 <= i < m.n_cells)) && return _text(404, "No such cell")
    r = series_data(s, i + 1, String(var); params)
    r === nothing && return _text(404, "No variable $var")
    ids, t, vals = r
    return _json((index = i, id = m.ids[i + 1], var = var, frames = ids, t = t,
                  values = [isfinite(x) ? x : nothing for x in vals]))
end

function _static(s::VisServer, parts)
    (length(parts) >= 2 && parts[1] in _STATIC_DIRS) || return _text(404, "Not found")
    for p in parts
        (p in (".", "..") || occursin('\\', p) || occursin(':', p) || startswith(p, ".")) &&
            return _text(403, "Forbidden")
    end
    root = abspath(s.viz_dir)
    f = abspath(joinpath(root, parts...))
    (startswith(f, root) && isfile(f)) || return _text(404, "Not found")
    basename(f) in ("config.json", "keys.json") && return _text(404, "Not found")   # never served as files
    ext = lowercase(splitext(f)[2])
    return HTTP.Response(200, ["Content-Type" => get(_MIME, ext, "application/octet-stream"),
                               "Cache-Control" => "no-cache"], read(f))
end

end # module VisualisationServer
