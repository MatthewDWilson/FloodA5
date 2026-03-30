# VisualisationServer.jl
# ----------------------
# Flood visualisation server: HTTP + WebSocket via Oxygen.jl
#
# Frame protocol (Option B — binary per-variable, on-demand fetch)
# ----------------------------------------------------------------
# Each frame stores a dense Float32 array per variable, values in mesh-cell
# index order (same order as the cell_ids array sent with the mesh).
# The client fetches one variable at a time via HTTP:
#
#   GET /frames/{idx}/{varname}   → application/octet-stream
#                                   raw little-endian Float32 array, n_cells × 4 bytes
#
# This scales to 1M+ cells: one variable = 4 MB vs ~25 MB JSON.
# Adding a new variable (e.g. contaminant) requires no structural changes —
# just include it in the vars dict passed to push_frame!.
#
# WebSocket messages
# ------------------
#   Server → Client  type="mesh"       data: GeoJSON FeatureCollection
#                                       cell_order: [id, id, ...] (mesh index order)
#   Server → Client  type="framecount"   count: N
#   Server → Client  type="newframe"     idx: N, t: Float64, vars: [name, ...]
#   Server → Client  type="simcomplete"  frames: N  (sent once when simulation ends)
#   Client → Server  (no messages expected; connection kept alive)
#
# HTTP endpoints
# --------------
#   GET /                         → redirect to /viz/index.html
#   GET /viz/{file}               → static files
#   GET /mesh                     → mesh GeoJSON + cell_order array (JSON)
#   GET /frames/count             → {count: N}
#   GET /frames/{idx}             → {t, vars: [name, ...]}  (metadata only)
#   GET /frames/{idx}/{varname}   → raw Float32 binary (application/octet-stream)
#   GET /status                   → server diagnostics
#
# Usage
# -----
#   server = VisualisationServer.start(port=8080)
#   VisualisationServer.set_mesh!(server, geojson_string, cell_ids)
#   VisualisationServer.push_frame!(server, t, Dict(
#       "depth"      => Float32.(state.water_depth),
#       "saturation" => Float32.(sat),
#       "volume"     => Float32.(state.volume),
#       "velocity"   => Float32.(state.velocity),
#   ))
#   VisualisationServer.stop(server)

module VisualisationServer

using Oxygen
using HTTP
using JSON3
using Sockets

export VisServer, start, stop, push_frame!, set_mesh!, notify_complete!

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
A single simulation snapshot.

`vars` maps variable name → dense Float32 vector in mesh-cell index order.
Using Float32 halves wire size vs Float64 with negligible precision loss for
visualisation purposes.
"""
struct Frame
    t    :: Float64
    vars :: Dict{String, Vector{Float32}}   # varname → values[n_cells]
end

mutable struct VisServer
    port         :: Int
    mesh_geojson :: Union{String, Nothing}
    cell_order   :: Vector{String}          # ordered cell IDs, set by set_mesh!
    frames       :: Vector{Frame}
    clients      :: Vector{HTTP.WebSockets.WebSocket}
    client_lock  :: ReentrantLock
    frame_lock   :: ReentrantLock
    task         :: Union{Task, Nothing}
    shutdown     :: Base.Event
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    start(; port=8080, viz_dir=...) → VisServer

Start the HTTP + WebSocket server non-blocking. Returns a VisServer handle.
"""
function start(; port::Int=8080, viz_dir::String=joinpath(@__DIR__, "viz"))
    server = VisServer(port, nothing, String[], Frame[],
                       HTTP.WebSockets.WebSocket[],
                       ReentrantLock(), ReentrantLock(), nothing, Base.Event())

    # --- Static files -------------------------------------------------------
    @get "/viz/{file}" function(req::HTTP.Request, file::String)
        path = joinpath(viz_dir, file)
        isfile(path) || return HTTP.Response(404, "Not found")
        return HTTP.Response(200, ["Content-Type" => _mime_type(file)], read(path))
    end

    @get "/" function(req::HTTP.Request)
        HTTP.Response(302, ["Location" => "/viz/index.html"], "")
    end

    # --- Mesh ---------------------------------------------------------------
    @get "/mesh" function(req::HTTP.Request)
        isnothing(server.mesh_geojson) &&
            return HTTP.Response(204, "No mesh loaded yet")
        # Return GeoJSON plus the ordered cell ID list so the client can
        # build its index mapping without parsing the full GeoJSON again.
        payload = JSON3.write((
            geojson     = JSON3.read(server.mesh_geojson),
            cell_order  = server.cell_order,
        ))
        HTTP.Response(200, ["Content-Type" => "application/json"], payload)
    end

    # --- Frame metadata -----------------------------------------------------
    @get "/frames/count" function(req::HTTP.Request)
        n = lock(server.frame_lock) do; length(server.frames); end
        HTTP.Response(200, ["Content-Type" => "application/json"],
                      JSON3.write((count=n,)))
    end

    # Frame metadata: time + list of available variable names
    @get "/frames/{idx}" function(req::HTTP.Request, idx::Int)
        lock(server.frame_lock) do
            (idx < 1 || idx > length(server.frames)) &&
                return HTTP.Response(404, "Frame out of range")
            f = server.frames[idx]
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write((t=f.t, vars=collect(keys(f.vars)))))
        end
    end

    # --- Binary variable data -----------------------------------------------
    # Returns a raw little-endian Float32 array in mesh-cell index order.
    # n_cells values × 4 bytes each.  The client decodes with:
    #   new Float32Array(await resp.arrayBuffer())
    @get "/frames/{idx}/{varname}" function(req::HTTP.Request,
                                             idx::Int, varname::String)
        lock(server.frame_lock) do
            (idx < 1 || idx > length(server.frames)) &&
                return HTTP.Response(404, "Frame out of range")
            f = server.frames[idx]
            haskey(f.vars, varname) ||
                return HTTP.Response(404, "Variable '$varname' not in frame")
            data = f.vars[varname]
            # Serialise as raw little-endian bytes.
            # Vector{UInt8}(reinterpret(...)) forces a concrete copy — needed
            # because HTTP.jl serialises the raw buffer and a lazy reinterpret
            # view can produce incorrect output depending on the runtime layout.
            bytes = Vector{UInt8}(reinterpret(UInt8, data))
            HTTP.Response(200,
                ["Content-Type"   => "application/octet-stream",
                 "Content-Length" => string(length(bytes)),
                 "X-Cell-Count"   => string(length(data)),
                 "X-Var-Name"     => varname],
                bytes)
        end
    end

    # --- Status -------------------------------------------------------------
    @get "/status" function(req::HTTP.Request)
        n_frames  = lock(server.frame_lock)  do; length(server.frames);  end
        n_clients = lock(server.client_lock) do; length(server.clients); end
        vars = n_frames > 0 ?
            lock(server.frame_lock) do; collect(keys(server.frames[end].vars)); end :
            String[]
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write((status="running", port=server.port,
                         mesh_loaded=!isnothing(server.mesh_geojson),
                         n_cells=length(server.cell_order),
                         frames=n_frames, live_clients=n_clients,
                         variables=vars)))
    end

    # --- WebSocket ----------------------------------------------------------
    @websocket "/live" function(ws::HTTP.WebSockets.WebSocket)
        lock(server.client_lock) do; push!(server.clients, ws); end
        n = lock(server.client_lock) do; length(server.clients); end
        @info "VisServer: client connected ($n total)"

        # Send mesh on connect (GeoJSON + cell_order)
        if !isnothing(server.mesh_geojson)
            payload = JSON3.write((
                type       = "mesh",
                data       = JSON3.read(server.mesh_geojson),
                cell_order = server.cell_order,
            ))
            HTTP.WebSockets.send(ws, payload)
        end

        # Tell client how many historical frames exist
        n_frames = lock(server.frame_lock) do; length(server.frames); end
        vars = n_frames > 0 ?
            lock(server.frame_lock) do; collect(keys(server.frames[end].vars)); end :
            String[]
        HTTP.WebSockets.send(ws, JSON3.write((
            type  = "framecount",
            count = n_frames,
            vars  = vars,
        )))

        try
            while !eof(ws)
                HTTP.WebSockets.receive(ws)
            end
        catch
        finally
            lock(server.client_lock) do
                filter!(c -> c !== ws, server.clients)
            end
            @info "VisServer: client disconnected"
        end
    end

    # --- Start server -------------------------------------------------------
    try
        test_sock = Sockets.listen(Sockets.IPv4(0), port)
        close(test_sock)
    catch
        error("VisServer: port $port is already in use. " *
              "Kill the previous Julia process or use --vis-port to choose another.")
    end

    server.task = @async begin
        @info "VisServer: listening on http://localhost:$port"
        try
            serve(port=port, async=false)
        catch e
            e isa InterruptException || @warn "VisServer task ended: $e"
        finally
            notify(server.shutdown)
        end
    end

    for _ in 1:20
        sleep(0.1)
        (istaskfailed(server.task) || istaskstarted(server.task)) && break
    end
    istaskfailed(server.task) &&
        error("VisServer failed to start: $(sprint(showerror, task_exception(server.task)))")

    return server
end

"""
    set_mesh!(server, geojson_string, cell_ids)

Store the mesh GeoJSON and the ordered cell ID list, then broadcast to any
already-connected clients.

`cell_ids` must be in the same order as the mesh cells — this order is used
as the index into all subsequent binary frame arrays.
"""
function set_mesh!(server::VisServer, geojson::String, cell_ids::Vector{String})
    server.mesh_geojson = geojson
    server.cell_order   = cell_ids
    payload = JSON3.write((
        type       = "mesh",
        data       = JSON3.read(geojson),
        cell_order = cell_ids,
    ))
    _broadcast!(server, payload)
end

"""
    push_frame!(server, t, vars)

Store a simulation timestep and notify live clients.

Arguments
---------
  t     Simulation time in seconds.
  vars  Dict mapping variable name → Float32 vector in mesh-cell index order.
        Example:
          Dict(
            "depth"      => Float32.(state.water_depth),
            "saturation" => Float32.(sat),
            "volume"     => Float32.(state.volume),
            "velocity"   => Float32.(state.velocity),
          )
        Any number of variables can be included; the client fetches whichever
        one it is currently displaying via GET /frames/{idx}/{varname}.

The WebSocket notification is a small JSON message listing the available
variable names — the client fetches binary data only on demand.
"""
function push_frame!(server::VisServer, t::Float64,
                     vars::Dict{String, Vector{Float32}})
    frame = Frame(t, vars)
    idx = lock(server.frame_lock) do
        push!(server.frames, frame)
        length(server.frames)
    end
    # Lightweight notification — no payload data over WebSocket
    _broadcast!(server, JSON3.write((
        type = "newframe",
        idx  = idx,
        t    = t,
        vars = collect(keys(vars)),
    )))
end

"""
    notify_complete!(server)

Broadcast a `simcomplete` message to all connected clients, signalling that
the simulation has finished and no further `newframe` messages will arrive.
The server continues running so the client can replay frames via HTTP.
Call this immediately after `run_simulation!` returns, before the keep-alive
loop.
"""
function notify_complete!(server::VisServer)
    n_frames = lock(server.frame_lock) do; length(server.frames); end
    _broadcast!(server, JSON3.write((
        type   = "simcomplete",
        frames = n_frames,
    )))
    @info "VisServer: simulation complete notification sent ($(n_frames) frames)"
end

"""
    stop(server)

Cleanly shut down the Oxygen server and signal the shutdown event.
"""
function stop(server::VisServer)
    try
        Oxygen.terminate()
    catch
    end
    notify(server.shutdown)
    @info "VisServer: stopped"
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _broadcast!(server::VisServer, msg::String)
    dead = HTTP.WebSockets.WebSocket[]
    lock(server.client_lock) do
        for ws in server.clients
            try
                HTTP.WebSockets.send(ws, msg)
            catch
                push!(dead, ws)
            end
        end
        filter!(c -> c ∉ dead, server.clients)
    end
end

function _mime_type(filename::String)::String
    ext = lowercase(splitext(filename)[2])
    get(Dict(
        ".html" => "text/html",
        ".js"   => "application/javascript",
        ".css"  => "text/css",
        ".json" => "application/json",
    ), ext, "text/plain")
end

end # module VisualisationServer
