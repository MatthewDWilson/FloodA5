"""
FloodModel.jl
-------------
Flood modelling application based on the A5 pentagonal DGGS.

Usage:
    julia --threads auto FloodModel.jl  --meshgen <aoi.geojson>  --meshres <N>
                                        [--meshout <file>]
                                        [--dem <file.tif>]  [--dem-strict]
                                        [--vis [mode]]  [--vis-port PORT]
                                        [--output <file.h5>]  [--output-interval SECS]

    julia --threads auto FloodModel.jl  --meshload <mesh.parquet>
                                        [--dem <file.tif>]  [--dem-strict]
                                        [--vis [mode]]  [--vis-port PORT]
                                        [--output <file.h5>]  [--output-interval SECS]

    julia FloodModel.jl --help | -h

DEM flags:
    --dem FILE          GeoTIFF elevation file. Samples the DEM onto the mesh
                        and stores elevation as a static variable. When used
                        with --meshout / --meshload the sampled elevation is
                        persisted in the parquet and does not need to be
                        re-sampled on future runs.
    --dem-strict        Error if any cell centre falls outside the DEM extent.
                        Default (no flag): assign NaN for out-of-bounds cells
                        and continue with a warning.

Output flags:
    --output FILE       Write simulation output to an HDF5 file (default: off).
                        Extension should be .h5 or .hdf5.
    --output-interval N Write a snapshot every N seconds of simulation time
                        (default: 60.0).
"""

push!(LOAD_PATH, @__DIR__)
# Guard against double-include when FloodModel.jl is included from a test
# harness that has already loaded A5Grid or stubbed the vis modules.
if !isdefined(Main, :A5Grid)
    include(joinpath(@__DIR__, "mesh", "A5Grid.jl"))
end
if !isdefined(Main, :VisualisationServer)
    include(joinpath(@__DIR__, "visualisation", "VisualisationServer.jl"))
end
if !isdefined(Main, :MakieVisualiser)
    include(joinpath(@__DIR__, "visualisation", "MakieVisualiser.jl"))
end

using .A5Grid
using .A5Grid: SGSTable, wse_from_volume, wetted_area_from_wse,
               flow_area_from_wse, wetted_perim_from_wse, hydraulic_radius_from_wse,
               sgs_table, build_sgs_tables!,
               grid_disk_neighbours, grid_disk_neighbours_batch,
               _polygon_area_m2, _shared_edge, _haversine_m
using .VisualisationServer
using .MakieVisualiser
using JSON3
using Dates
using Statistics: mean
using HDF5

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Push a frame to the visualiser every N simulation steps.
const VIS_INTERVAL = 1

# Supported visualisation modes.
const VIS_MODES = (:cesium, :makie)

# Number of sides (neighbours) of an A5 pentagon interior cell.
# The per-edge donor limiter caps each edge transfer at volume[donor] / DONOR_EDGE_DIVISOR.
# Divisor = 2*N_SIDES = 10: even if all 5 edges of a cell fire simultaneously,
# the maximum total outflow is 5 * V/10 = V/2 = 50%, preserving the intended
# half-step stability criterion.  (Divisor = N_SIDES = 5 allowed 100% drainage
# when all edges fired together, driving oscillations — Bug 49 Fix A.)
const N_SIDES           = 5
const DONOR_EDGE_DIVISOR = 2 * N_SIDES   # = 10

# ---------------------------------------------------------------------------
# Flow model types
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Flow model types
# ---------------------------------------------------------------------------

"""
    FlowMethod — abstract supertype for solver selection.

  StandardFlow(; manning_n=0.03)
      Diffusive-wave solver on mean cell elevation.  WSE = bed + depth.
      Fast; suitable for initial runs and validation.

  SGSFlow(; manning_n=0.03)
      Diffusive-wave solver with sub-grid hypsometric storage tables.
      Allows partial cell wetting.  Requires `build_sgs_tables!` to have
      been run; tables are loaded from the `.sgs.h5` companion file and
      stored in `FlowState.sgs_tables` at initialisation time.
"""
abstract type FlowMethod end

struct StandardFlow <: FlowMethod
    manning_n :: Float64
end
StandardFlow(; manning_n::Float64 = 0.03) = StandardFlow(manning_n)

struct SGSFlow <: FlowMethod
    manning_n :: Float64
end
SGSFlow(; manning_n::Float64 = 0.03) = SGSFlow(manning_n)

"""
    EdgeList

All undirected edges in the mesh, indexed 1:n_edges.  Each edge is stored
exactly once — from the perspective of the cell with the lower array index
(`cell_i < cell_j` always).  This canonical ordering is enforced by
`_build_edge_list` and is the basis for correct single-pass flux computation.

Sign convention for `flux`:
  positive  → flow from cell_i to cell_j
  negative  → flow from cell_j to cell_i

Forward-compatibility for multi-resolution (Phase 3):
  - `res_i`, `res_j` — A5 resolution levels (added in Phase 3)
  - `active`         — bitmask for AMR active/inactive edges (added in Phase 3)
  - The struct is sized to `max_edges` at initialisation; `n_edges` is the
    current active count (for single-resolution, n_edges == max_edges).
"""
struct EdgeList
    n_edges   :: Int
    cell_i    :: Vector{Int}        # lower cell index  (cell_i < cell_j always)
    cell_j    :: Vector{Int}        # higher cell index
    width     :: Vector{Float64}    # shared edge length (m)
    L         :: Vector{Float64}    # centre-to-centre haversine distance (m)
    cos_theta :: Vector{Float64}    # non-orthogonality correction (1.0 = orthogonal)
    sill      :: Vector{Float64}    # sill elevation (m) — bed or SGS minimum
    flux      :: Vector{Float64}    # q (m²/s) at t-dt, signed cell_i → cell_j
                                    #   used by standard flow and SGS Bates kernel
    flux_Q    :: Vector{Float64}    # Q (m³/s) at t-dt, signed cell_i → cell_j
                                    #   used by SGS R-A kernel only; zero otherwise
end

"""
Hydrodynamic state of the flood model at a single timestep.

Primary state variable is `volume` (m³ stored water per cell).  All other
fields are either static mesh geometry or diagnostics updated each step.
Fields are in mesh.cells order (index i ↔ cell i).

`edges` holds all undirected mesh edges; `adj_matrix` is retained for
neighbour queries (wet/dry triggers, BC injection) independently of the
flux computation.
"""
mutable struct FlowState
    cell_ids    :: Vector{String}
    water_depth :: Vector{Float64}   # m above local bed (diagnostic)
    volume      :: Vector{Float64}   # m³ stored — primary state variable
    velocity    :: Vector{Float64}   # m/s scalar magnitude (diagnostic)
    vel_u       :: Vector{Float64}   # m/s eastward  velocity component (lon direction)
    vel_v       :: Vector{Float64}   # m/s northward velocity component (lat direction)
    elevation   :: Vector{Float64}   # bed elevation above datum (m)
    manning_n   :: Vector{Float64}   # Manning's roughness per cell
    cell_area   :: Vector{Float64}   # plan area (m²)
    cell_lons   :: Vector{Float64}   # cell centre longitude (degrees, EPSG:4326)
    cell_lats   :: Vector{Float64}   # cell centre latitude  (degrees, EPSG:4326)
    adjacency   :: Dict{String, Vector{String}}  # cell_id → neighbour ids
    adj_matrix  :: Matrix{Int}       # (max_nb × n_cells) index matrix, 0=none
    edges       :: EdgeList          # all undirected edges, computed once
    sgs_tables  :: Vector{Any}       # Vector{SGSTable} (SGSFlow) or empty
end

# ---------------------------------------------------------------------------
# HDF5 output
# ---------------------------------------------------------------------------

"""
    SimOutput

Manages HDF5 file output for dynamic simulation state.

Structure of the HDF5 file
--------------------------
/mesh/
    cell_ids          String dataset  (n_cells,)
    elevations        Float64 dataset (n_cells,)
    center_lons       Float64 dataset (n_cells,)
    center_lats       Float64 dataset (n_cells,)

/frames/
    /0001/
        t             scalar Float64 — simulation time in seconds
        water_depth   Float64 dataset (n_cells,) — depth above cell z_min (m)
        volume        Float64 dataset (n_cells,) — stored volume (m³)
        saturation    Float64 dataset (n_cells,) — fractional wetted area [0-1]
        velocity      Float64 dataset (n_cells,) — scalar velocity (m/s)
    /0002/ ...

Rationale for HDF5
------------------
Dynamic outputs are stored separately from the mesh parquet because:
  - Frame counts can reach thousands; parquet is not optimised for this pattern
  - HDF5 supports chunked / compressed storage for large arrays
  - Post-processing tools (Python xarray, Julia HDF5.jl, QGIS) read HDF5 natively
  - The mesh can be reloaded from parquet independently of any output file

Chunking and deflate (gzip level 4) compression are applied to all
per-cell datasets. chunk_size = min(n_cells, 4096).
"""
mutable struct SimOutput
    path             :: String
    output_interval  :: Float64    # seconds of sim time between writes
    last_write_t     :: Float64    # sim time of last write
    frame_count      :: Int
    enabled          :: Bool
end

SimOutput(; path="", interval=60.0, enabled=false) =
    SimOutput(path, interval, -Inf, 0, enabled)

"""
    _write_mesh_metadata!(output, mesh)

Write the static mesh metadata to the HDF5 file's /mesh/ group.
Called once when the output file is first opened.
"""
function _write_mesh_metadata!(output::SimOutput, mesh::A5Mesh)
    output.enabled || return
    HDF5.h5open(output.path, "w") do fid
        g = HDF5.create_group(fid, "mesh")
        g["cell_ids"]    = [c.id         for c in mesh.cells]
        g["center_lons"] = [c.center_lon for c in mesh.cells]
        g["center_lats"] = [c.center_lat for c in mesh.cells]
        g["elevations"]  = get(mesh.static_vars, "elevation",
                               fill(NaN, length(mesh.cells)))
        for (varname, vals) in mesh.static_vars
            varname == "elevation" && continue
            g[varname] = vals
        end
        HDF5.create_group(fid, "frames")
    end
    @info "HDF5 output file created: $(output.path)"
end

"""
    _write_frame!(output, state, t)

Append a simulation snapshot to /frames/ in the HDF5 file.
Frame group names are zero-padded to 6 digits for lexicographic ordering.

Datasets per frame:
  t            — simulation time (s)
  water_depth  — depth above cell minimum (m), diagnostic
  volume       — stored volume per cell (m³), primary state variable
  saturation   — fractional wetted area [0–1] (SGS only; 1.0 where wet otherwise)
  velocity     — scalar velocity magnitude (m/s)
"""
function _write_frame!(output::SimOutput, state::FlowState, t::Float64)
    output.enabled || return
    output.frame_count += 1
    frame_name = lpad(string(output.frame_count), 6, '0')
    sat   = saturation_fraction(state)
    n     = length(state.cell_ids)
    chunk = (min(n, 4096),)   # chunk size: one contiguous read ≤ 4096 cells
    HDF5.h5open(output.path, "r+") do fid
        g = HDF5.create_group(fid["frames"], frame_name)
        # Scalar — no chunking needed
        g["t"] = t
        # Per-cell arrays — chunked + deflate (gzip level 4)
        for (name, data) in (("water_depth", state.water_depth),
                              ("volume",      state.volume),
                              ("saturation",  sat),
                              ("velocity",    state.velocity))
            ds = HDF5.create_dataset(g, name, eltype(data), (n,);
                                     chunk=chunk, deflate=4)
            write(ds, data)
        end
        # Write flux_Q (m³/s, SGS R-A only) when non-zero.
        # Written per-edge rather than per-cell; only present in SGS R-A runs.
        ne = state.edges.n_edges
        if ne > 0 && any(!=(0.0), state.edges.flux_Q)
            eq_chunk = (min(ne, 4096),)
            ds = HDF5.create_dataset(g, "flux_Q", Float64, (ne,);
                                     chunk=eq_chunk, deflate=4)
            write(ds, state.edges.flux_Q)
        end
    end
    output.last_write_t = t
end

"""
    _should_write_frame(output, t) → Bool

Return true if enough sim time has elapsed since the last write.
"""
function _should_write_frame(output::SimOutput, t::Float64)::Bool
    output.enabled || return false
    return (t - output.last_write_t) >= output.output_interval
end

# ---------------------------------------------------------------------------
# Flow model helpers
# ---------------------------------------------------------------------------

"""Haversine-based angular distance in degrees — fast scalar approximation."""
function _haversine_deg(lon1, lat1, lon2, lat2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    return sqrt(dlat^2 + (dlon * cosd((lat1 + lat2) / 2))^2)
end

"""
    _build_adjacency_shared_vertices(mesh) → Dict{String,Vector{String}}

Build exact edge-sharing adjacency from cell boundary vertices.
Two cells are edge-sharing neighbours iff their boundaries share exactly
2 vertices (the endpoints of the shared edge).

This bypasses pya5's grid_disk, which returns compact/mixed-resolution
results that cause some cells to report fewer than 5 neighbours.
Vertex coordinates are rounded to 7 decimal places (~1 cm) to handle
floating-point differences between pya5's boundary representations for
adjacent cells.
"""
function _build_adjacency_shared_vertices(mesh::A5Mesh)::Dict{String,Vector{String}}
    PREC = 1e7   # multiply then round to get 7 decimal places

    norm = id -> A5Grid._to_hex(parse(UInt64, id, base=16))
    ids  = [norm(c.id) for c in mesh.cells]
    n    = length(ids)

    # Build vertex → cell-index index
    vertex_to_cells = Dict{Tuple{Int64,Int64}, Vector{Int}}()
    for i in 1:n
        bnd = mesh.cells[i].boundary
        # Exclude closing vertex (duplicate of first) if present
        verts = (length(bnd) > 1 && bnd[1] ≈ bnd[end]) ? bnd[1:end-1] : bnd
        for v in verts
            key = (round(Int64, v[1] * PREC), round(Int64, v[2] * PREC))
            push!(get!(vertex_to_cells, key, Int[]), i)
        end
    end

    # For each cell count shared vertices with every other cell
    adj = Dict{String,Vector{String}}(id => String[] for id in ids)
    for i in 1:n
        bnd = mesh.cells[i].boundary
        verts = (length(bnd) > 1 && bnd[1] ≈ bnd[end]) ? bnd[1:end-1] : bnd
        shared = Dict{Int,Int}()
        for v in verts
            key = (round(Int64, v[1] * PREC), round(Int64, v[2] * PREC))
            for j in get(vertex_to_cells, key, Int[])
                j == i && continue
                shared[j] = get(shared, j, 0) + 1
            end
        end
        adj[ids[i]] = [ids[j] for (j, cnt) in shared if cnt == 2]
    end

    n_edges = sum(length(v) for v in values(adj)) ÷ 2
    @info "Adjacency built: $n cells, $n_edges undirected edges (shared-vertex method)"
    return adj
end

"""
    _edge_length_m(bnd_i, bnd_j) → Float64

Haversine length in metres of the shared edge between two A5 pentagon cells.
Returns NaN if no shared edge is found.
"""
function _edge_length_m(bnd_i::Vector{Vector{Float64}},
                         bnd_j::Vector{Vector{Float64}})::Float64
    edge = A5Grid._shared_edge(bnd_i, bnd_j)
    edge === nothing && return NaN
    lon1, lat1, lon2, lat2 = edge
    return A5Grid._haversine_m(lon1, lat1, lon2, lat2)
end

"""
    _cell_area_m2(cell) → Float64

Plan area of an A5 cell in m² (equirectangular approximation, ~0.1% error
at resolution 14).
"""
_cell_area_m2(cell::A5Cell) = A5Grid._polygon_area_m2(cell.boundary)

"""
    _build_adjacency_matrix!(adj_matrix, cells, id_idx, adj)

Populate the adjacency index matrix in-place.  `adj_matrix[slot, i]` holds
the array index of the slot-th neighbour of cell i (0 = no neighbour).
Kept separate from edge geometry so it can be updated cheaply during AMR
without rebuilding the full EdgeList.
"""
function _build_adjacency_matrix!(adj_matrix :: Matrix{Int},
                                   cells      :: Vector{A5Cell},
                                   id_idx     :: Dict{String,Int},
                                   adj        :: Dict{String,Vector{String}})
    n      = length(cells)
    _norm(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
    ids    = [_norm(c.id) for c in cells]
    max_nb = size(adj_matrix, 1)
    for i in 1:n
        nbrs = get(adj, ids[i], String[])
        for (slot, nb_id) in enumerate(nbrs)
            slot > max_nb && break
            j = get(id_idx, nb_id, 0)
            j == 0 && continue
            adj_matrix[slot, i] = j
        end
    end
end

"""
    _build_edge_list(cells, id_idx, adj, areas, sill_matrix) → EdgeList

Build a flat list of all undirected edges in the mesh.  Each edge is stored
exactly once with `cell_i < cell_j` (canonical lower-index convention).

Accepts an optional pre-computed `sill_matrix` (max_nb × n) from SGS
pre-processing.  If `nothing`, sills default to `max(elev_i, elev_j)` and
must be overridden by the caller for SGS runs.

Works for mixed A5 resolution levels — cell geometry is read directly from
`A5Cell.boundary` and centre coordinates, making no assumption of uniform
cell size.  This makes the function valid for Phase 3 multi-resolution meshes
without modification.

Logs the edge count, and the min/mean/max of cos θ for the non-orthogonality
distribution of the loaded mesh.
"""

"""
    _check_mesh_connectivity(edges, n_cells) → Vector{Int}

Finds connected components in the cell graph using BFS.
Returns a component label vector (length n_cells) where cells with the
same label are mutually reachable via edges.  Logs a warning if the mesh
has more than one component, as isolated cells will never receive flux.
"""
function _check_mesh_connectivity(edges::EdgeList, n_cells::Int)::Vector{Int}
    # Build adjacency list from EdgeList
    adj = [Int[] for _ in 1:n_cells]
    for e in 1:edges.n_edges
        push!(adj[edges.cell_i[e]], edges.cell_j[e])
        push!(adj[edges.cell_j[e]], edges.cell_i[e])
    end

    component = zeros(Int, n_cells)
    comp_id   = 0

    for start in 1:n_cells
        component[start] != 0 && continue
        comp_id += 1
        queue = [start]
        component[start] = comp_id
        while !isempty(queue)
            node = popfirst!(queue)
            for nb in adj[node]
                if component[nb] == 0
                    component[nb] = comp_id
                    push!(queue, nb)
                end
            end
        end
    end

    if comp_id > 1
        comp_sizes = [count(==(k), component) for k in 1:comp_id]
        largest    = argmax(comp_sizes)
        isolated   = n_cells - comp_sizes[largest]
        @warn "Mesh has $comp_id disconnected components. " *
              "Largest: $(comp_sizes[largest]) cells. " *
              "Isolated (unreachable from largest): $isolated cells. " *
              "These cells will never receive flux — consider a larger AOI."
        for k in 1:comp_id
            sz = comp_sizes[k]
            mark = k == largest ? " (largest — source should be here)" : " (isolated)"
            @info "  Component $k: $sz cells$mark"
        end
    else
        @info "Mesh connectivity: fully connected ($n_cells cells, 1 component)"
    end

    return component
end

function _build_edge_list(cells       :: Vector{A5Cell},
                           id_idx      :: Dict{String,Int},
                           adj         :: Dict{String,Vector{String}},
                           areas       :: Vector{Float64},
                           sill_matrix :: Union{Matrix{Float64}, Nothing} = nothing,
                           elevations  :: Union{Vector{Float64}, Nothing} = nothing)::EdgeList
    n    = length(cells)
    _norm(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
    ids  = [_norm(c.id) for c in cells]

    # Pre-size to worst case (5 edges per pentagon, each shared once = 5n/2)
    # +10% margin for boundary effects.  We compact at the end.
    max_e = n * 3   # conservative upper bound: 5n/2 rounded up with margin
    ci    = Vector{Int}(undef, max_e)
    cj    = Vector{Int}(undef, max_e)
    ws    = Vector{Float64}(undef, max_e)
    Ls    = Vector{Float64}(undef, max_e)
    cts   = Vector{Float64}(undef, max_e)
    sls   = Vector{Float64}(undef, max_e)

    e = 0   # edge counter
    seen = Set{Tuple{Int,Int}}()   # (min,max) pairs already added

    for i in 1:n
        cell_i = cells[i]
        nbrs   = get(adj, ids[i], String[])
        for nb_id in nbrs
            j = get(id_idx, nb_id, 0)
            j == 0 && continue           # external / out-of-mesh neighbour

            lo, hi = minmax(i, j)
            (lo, hi) in seen && continue # edge already recorded
            push!(seen, (lo, hi))

            e += 1
            if e > max_e                 # grow if needed (shouldn't happen)
                append!(ci,  zeros(Int, n))
                append!(cj,  zeros(Int, n))
                append!(ws,  zeros(Float64, n))
                append!(Ls,  zeros(Float64, n))
                append!(cts, zeros(Float64, n))
                append!(sls, zeros(Float64, n))
                max_e += n
            end

            cell_j = cells[j]
            ci[e] = lo
            cj[e] = hi

            # ── Edge width ────────────────────────────────────────────────
            w = _edge_length_m(cell_i.boundary, cell_j.boundary)
            ws[e] = isnan(w) ? sqrt(areas[i]) : w

            # ── Centre-to-centre distance ─────────────────────────────────
            Ls[e] = A5Grid._haversine_m(
                cell_i.center_lon, cell_i.center_lat,
                cell_j.center_lon, cell_j.center_lat)

            # ── Non-orthogonality correction ──────────────────────────────
            ct = A5Grid._edge_cos_theta(
                cell_i.boundary, cell_j.boundary,
                cell_i.center_lon, cell_i.center_lat,
                cell_j.center_lon, cell_j.center_lat)
            cts[e] = isnan(ct) ? 1.0 : ct  # fallback: assume orthogonal

            # ── Sill elevation ────────────────────────────────────────────
            # Prefer SGS pre-computed sill; fall back to max(elev_i, elev_j).
            # The caller is responsible for passing the correct sill_matrix
            # for SGS runs; for standard runs this defaults to bed elevation.
            #
            # SGS sill lookup: sill_matrix[slot, cell] stores the minimum
            # DEM elevation along the shared edge, indexed by adjacency slot.
            # The slot ordering in sill_matrix was set during build_sgs_tables!
            # using grid_disk_neighbours_batch(), which may differ from the
            # adjacency ordering stored in the parquet adj columns and loaded
            # into the adj dict here.  To guard against ordering mismatches
            # we search BOTH the lo-cell and hi-cell perspectives and take
            # whichever returns a valid (non-NaN) result first.
            sls[e] = if sill_matrix !== nothing
                s = NaN
                # Try lo → hi
                lo_nbrs = get(adj, ids[lo], String[])
                for (slot, nb_id2) in enumerate(lo_nbrs)
                    slot > size(sill_matrix, 1) && break
                    if get(id_idx, nb_id2, 0) == hi
                        s = sill_matrix[slot, lo]
                        break
                    end
                end
                # If lo lookup failed, try hi → lo
                if isnan(s)
                    hi_nbrs = get(adj, ids[hi], String[])
                    for (slot, nb_id2) in enumerate(hi_nbrs)
                        slot > size(sill_matrix, 1) && break
                        if get(id_idx, nb_id2, 0) == lo
                            s = sill_matrix[slot, hi]
                            break
                        end
                    end
                end
                isnan(s) ? (elevations !== nothing ?
                    (isnan(elevations[lo]) || isnan(elevations[hi]) ? NaN :
                    max(elevations[lo], elevations[hi])) : 0.0) : s
            else
                elevations !== nothing ?
                    (isnan(elevations[lo]) || isnan(elevations[hi]) ? NaN :
                    max(elevations[lo], elevations[hi])) : 0.0
            end
        end
    end

    n_edges = e
    @info "Edge list built: $n_edges edges for $n cells ($(round(n_edges/max(n,1), digits=2)) edges/cell)"
                   
    valid_ct = filter(isfinite, cts[1:n_edges])
    if !isempty(valid_ct)
        @info "Edge non-orthogonality (cos θ):  min=$(round(minimum(valid_ct),digits=3))  mean=$(round(mean(valid_ct),digits=3))  max=$(round(maximum(valid_ct),digits=3))"
    end

    if sill_matrix !== nothing
        # Diagnostic: count how many edges used the SGS sill vs the fallback.
        # Any edge where sls[e] == max(elev_lo, elev_hi) AND elevations are valid
        # may have fallen back.  We count NaN sills (geometry fallback) separately.
        n_nan_sill  = count(isnan,  sls[1:n_edges])
        n_sgs_sill  = n_edges - n_nan_sill   # all finite sills came from somewhere
        # Warn if a significant fraction are NaN — indicates sill lookup failures
        if n_nan_sill > 0
            @warn "SGS edge sills: $n_nan_sill / $n_edges edges have NaN sill " *
                  "(no DEM data along edge). These edges will have no flux. " *
                  "Check DEM coverage."
        else
            @info "SGS edge sills: all $n_edges edges have valid sill elevations ✓"
        end
    end

    return EdgeList(
        n_edges,
        ci[1:n_edges],
        cj[1:n_edges],
        ws[1:n_edges],
        Ls[1:n_edges],
        cts[1:n_edges],
        sls[1:n_edges],
        zeros(Float64, n_edges),   # flux  (m²/s) — standard flow and SGS Bates
        zeros(Float64, n_edges),   # flux_Q (m³/s) — SGS R-A kernel only
    )
end

# ---------------------------------------------------------------------------
# Flow model initialisation
# ---------------------------------------------------------------------------

"""
    initialise_flow_model(mesh, method;
                          manning_n=0.03, friction_raster=nothing) → FlowState

Initialise the flow model on the A5 mesh.

- `method`          — `StandardFlow()` or `SGSFlow()`
- `manning_n`       — global Manning's roughness (default 0.03)
- `friction_raster` — path to a friction GeoTIFF; per-cell n overrides global value

For SGSFlow, `build_sgs_tables!` must have been called on the mesh first
(tables are stored in `mesh.array_vars`).
"""
function initialise_flow_model(mesh::A5Mesh,
                                method::FlowMethod = StandardFlow();
                                manning_n::Float64 = 0.03,
                                friction_raster    = nothing)::FlowState
    n       = length(mesh)
    # Normalise cell IDs to 16-char zero-padded hex throughout — ensures
    # consistency between parquet-stored IDs (via pya5 u64_to_hex, may omit
    # leading zeros) and IDs returned by grid_disk_neighbours (always 16 chars).
    _norm_id(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
    ids     = [_norm_id(c.id) for c in mesh.cells]
    id_idx  = Dict{String,Int}(ids[i] => i for i in 1:n)

    # ── Elevation ──────────────────────────────────────────────────────────
    elevations = if haskey(mesh.static_vars, "elevation")
        copy(mesh.static_vars["elevation"])
    else
        @warn "No elevation data — all cells at z=0. Use --dem before running."
        zeros(Float64, n)
    end

    nan_elev_idx = findall(isnan, elevations)
    if !isempty(nan_elev_idx)
        @warn "$(length(nan_elev_idx)) cells have NaN elevation and will be hydraulically inert (no flux on any adjacent edge). Indices: $(nan_elev_idx[1:min(10,end)])"
    end

    # ── Manning's n ────────────────────────────────────────────────────────
    n_vec = fill(manning_n, n)
    if friction_raster !== nothing
        @info "Sampling friction raster onto mesh..."
        sample_dem_centroid!(mesh, FileDEM(friction_raster); var_name="friction")
        fr = mesh.static_vars["friction"]
        for i in 1:n
            isfinite(fr[i]) && fr[i] > 0.0 && (n_vec[i] = fr[i])
        end
    elseif haskey(mesh.static_vars, "friction")
        fr = mesh.static_vars["friction"]
        for i in 1:n
            isfinite(fr[i]) && fr[i] > 0.0 && (n_vec[i] = fr[i])
        end
    end

    # ── Cell areas ─────────────────────────────────────────────────────────
    # Always recompute cell areas from polygon boundaries rather than trusting
    # the parquet value. The Python bridge stores a "cell_area" column but it
    # may reflect a total-AOI or approximate value rather than the exact geodetic
    # area of each individual pentagon. _polygon_area_m2 is fast (< 1ms for
    # typical mesh sizes) and guarantees correctness.
    # SGS builds store "sgs_cell_area" computed the same way during build_sgs_tables!
    # so we keep that path as it is already correct and consistent with the tables.
    areas = if haskey(mesh.static_vars, "sgs_cell_area")
        copy(mesh.static_vars["sgs_cell_area"])
    else
        [_polygon_area_m2(c.boundary) for c in mesh.cells]
    end
    # Persist cell areas into static_vars so they are written into the parquet
    # on --meshout and available without recomputation on --meshload.
    # "sgs_cell_area" is the canonical key; for non-SGS meshes we write it
    # here so downstream tools (and the parquet) always have a cell_area column.
    if !haskey(mesh.static_vars, "sgs_cell_area")
        mesh.static_vars["sgs_cell_area"] = copy(areas)
    end

    # ── Debug: cell area sanity check ──────────────────────────────────────
    n_nan_area = count(isnan, areas)
    n_zero_area = count(a -> a < 1.0, areas)
    @info "Cell areas: min=$(round(minimum(areas), sigdigits=4)) max=$(round(maximum(areas), sigdigits=4)) " *
          "mean=$(round(sum(areas)/length(areas), sigdigits=4)) NaN=$n_nan_area zero/small=$n_zero_area"
    if n_nan_area > 0 || n_zero_area > 0
        @warn "Cell area problems detected — NaN or near-zero areas will cause NaN/Inf volumes"
    end

    # ── Adjacency ───────────────────────────────────────────────────────────
    # Prefer adjacency pre-computed during mesh generation and stored in the
    # parquet (mesh.adjacency).  If absent (old parquets), fall back to the
    # shared-vertex method computed from boundary geometry — which is
    # geometrically exact and does not call pya5 at all.
    t0  = time()
    adj = if !isempty(mesh.adjacency)
        @info "Using pre-computed adjacency from mesh parquet..."
        # Normalise keys to match ids (16-char zero-padded hex)
        norm = id -> A5Grid._to_hex(parse(UInt64, id, base=16))
        Dict{String,Vector{String}}(
            norm(k) => [norm(nb) for nb in v]
            for (k, v) in mesh.adjacency
        )
    else
        @info "No adjacency in parquet — computing via shared-vertex detection..."
        _build_adjacency_shared_vertices(mesh)
    end
    @info "  Adjacency ready in $(round(time()-t0, digits=1))s ($(length(adj)) cells)"

    max_nb     = 5
    adj_matrix = zeros(Int, max_nb, n)

    _build_adjacency_matrix!(adj_matrix, mesh.cells, id_idx, adj)

    # ── SGS tables ─────────────────────────────────────────────────────────
    sgs_tables  = SGSTable[]
    sill_matrix = nothing   # will be set below for SGS runs

    if method isa SGSFlow
        haskey(mesh.array_vars, "sgs_elev_bins") ||
            error("SGS tables not found in mesh — call build_sgs_tables! first.")
        @info "Loading SGS tables from mesh array_vars..."
        sgs_tables = [sgs_table(mesh, i) for i in 1:n]

        if haskey(mesh.array_vars, "sgs_edge_sills")
            sill_matrix = copy(mesh.array_vars["sgs_edge_sills"])
        else
            @warn "sgs_edge_sills not found — using cell z_min fallback."
            z_mins = get(mesh.static_vars, "sgs_z_min", elevations)
            sill_matrix = fill(NaN, max_nb, n)
            for i in 1:n, slot in 1:max_nb
                j = adj_matrix[slot, i]
                j == 0 && continue
                sill_matrix[slot, i] = min(z_mins[i], z_mins[j])
            end
        end
    end

    # ── Edge list ───────────────────────────────────────────────────────────
    @info "Building edge list..."
    t1 = time()
    edges = _build_edge_list(mesh.cells, id_idx, adj, areas,
                              sill_matrix, elevations)
    @info "  Edge list built in $(round(time()-t1, digits=1))s"

    # Check graph connectivity -- warn if mesh has isolated components
    _check_mesh_connectivity(edges, n)

    volumes = zeros(Float64, n)

    lons = [c.center_lon for c in mesh.cells]
    lats = [c.center_lat for c in mesh.cells]

    return FlowState(
        ids,
        zeros(Float64, n),
        volumes,
        zeros(Float64, n),
        zeros(Float64, n),
        zeros(Float64, n),
        elevations,
        n_vec,
        areas,
        lons,
        lats,
        adj,
        adj_matrix,
        edges,
        sgs_tables,
    )
end

# Backwards-compatible no-method overload
initialise_flow_model(mesh::A5Mesh) = initialise_flow_model(mesh, StandardFlow())


# ---------------------------------------------------------------------------
# Diffusive wave physics — Bates et al. (2010) inertial formulation
# ---------------------------------------------------------------------------
# Pure physics kernels (_bates_flux, _bates_flux_limited, _manning_flux_ra,
# _adjacency_slot, _cfl_dt) and their constants (_G, HFLOW_THRESHOLD,
# FROUDE_LIMIT) are defined in surfacewater/flow2d.jl.
#
# That file has no FlowState dependency and can be included independently by
# external code (GPU kernel wrappers, 3DGeo plugins, unit tests) without
# loading the full application stack.
#
# ENV note: this include contains only definitions — no top-level execution.

include(joinpath(@__DIR__, "surfacewater", "flow2d.jl"))

"""
    step_standard!(state, dt)

Advance the standard (non-SGS) model by one timestep `dt` using the
Bates et al. (2010) inertial formulation.

WSE is derived from `volume / cell_area + elevation`.  Volume is the
primary state variable; water_depth is a derived diagnostic.  The sill
between adjacent cells is the maximum of their two bed elevations.  The
slope distance L is stored in `state.edges.L`.
"""

# ---------------------------------------------------------------------------
# Volume/depth update helpers — called by both step functions
# ---------------------------------------------------------------------------

"""
    _apply_dV_standard!(state, dV)

Apply a vector of net volume increments `dV` (m³) to a `StandardFlow` state.

`volume` is the primary state variable (mirrors SGS).  `water_depth` is
derived as `volume / cell_area` after each update.

The per-edge donor limiter in `step_standard!` already ensures no cell loses
more than `volume / N_SIDES` per edge, so net outflow across all edges is
bounded at 50%.  The `max(0.0, …)` floor here is a last-resort guard only.
"""
function _apply_dV_standard!(state::FlowState, dV::Vector{Float64})
    n = length(state.cell_ids)
    Threads.@threads for i in 1:n
        A_i = state.cell_area[i]
        A_i < 1.0 && continue
        state.volume[i]      = max(0.0, state.volume[i] + dV[i])
        state.water_depth[i] = state.volume[i] / A_i
    end
end

"""
    _apply_dV_sgs!(state, dV)

Apply a vector of net volume increments `dV` (m³) to an `SGSFlow` state.

Mirrors `_apply_dV_standard!` but operates on `state.volume` directly (the
SGS primary state variable) and derives `water_depth` from the hypsometric
lookup.

The per-edge donor limiter in `step_sgs!` already ensures no cell loses more
than `volume / N_SIDES` per edge.  The `max(0.0, …)` floor here is a
last-resort guard only.
"""
function _apply_dV_sgs!(state::FlowState, dV::Vector{Float64})
    n = length(state.cell_ids)
    Threads.@threads for i in 1:n
        state.volume[i] = max(0.0, state.volume[i] + dV[i])
        tbl = state.sgs_tables[i]
        if isnan(tbl.z_min)
            state.water_depth[i] = 0.0
            continue
        end
        new_wse = wse_from_volume(tbl, state.volume[i])
        # water_depth is the physically meaningful depth above the cell thalweg,
        # used for display, HDF5 output, and the velocity computation.
        # For overfull cells (volume > vol_curve[end]), wse_from_volume extrapolates
        # linearly above z_max (Bug 52 fix).  We report the true depth here so that
        # channel cells show their actual inundation depth in the visualiser.
        # The CFL timestep uses a separate percentile-capped depth (see _cfl_dt)
        # to prevent a single deep cell from forcing the whole domain to tiny dt.
        state.water_depth[i] = max(0.0, new_wse - tbl.z_min)
    end
end

# ---------------------------------------------------------------------------
# Velocity computation — called at end of every step_standard! / step_sgs!
# ---------------------------------------------------------------------------

"""
    _compute_velocity!(state)

Compute per-cell scalar speed and (u, v) velocity components from the edge
flux array stored in `state.edges.flux`.

Each undirected edge `e` carries unit discharge `q` (m²/s) with sign
convention: `q > 0` means flow from `cell_j` to `cell_i`.  The volumetric
flux is `Q = q * width` (m³/s).

The (u, v) components are accumulated in geographic coordinates:
  - u  eastward  (positive = flow toward increasing longitude)
  - v  northward (positive = flow toward increasing latitude)

The unit direction vector from cell_i to cell_j is:
  dx = lon_j - lon_i,  dy = lat_j - lat_i   (un-normalised geographic delta)
normalised to a unit vector.  When `Q > 0` (flow j→i), the momentum
contribution to cell_i is in the −(dx,dy) direction and to cell_j in the
+(dx,dy) direction, and vice versa when `Q < 0`.

After accumulation each cell's speed is:
  speed = sqrt(u² + v²)
  u, v are normalised to m/s by dividing by `cell_area * max(depth, h_min)`.

Dry cells (depth < h_min) have velocity zeroed.
"""
function _compute_velocity!(state::FlowState)
    n      = length(state.cell_ids)
    edges  = state.edges
    h_min  = 1e-4   # m — depth threshold below which velocity = 0

    # Flux-weighted momentum scatter (serial — ci/cj not unique across edges)
    sum_u  = zeros(Float64, n)
    sum_v  = zeros(Float64, n)

    @inbounds for e in 1:edges.n_edges
        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        dlon = state.cell_lons[cj] - state.cell_lons[ci]
        dlat = state.cell_lats[cj] - state.cell_lats[ci]
        dist = sqrt(dlon*dlon + dlat*dlat)
        dist < 1e-12 && continue
        ux = dlon / dist
        uy = dlat / dist

        # Use flux_Q (m³/s) for the R-A SGS kernel; fall back to flux*width for
        # standard flow (flux is unit discharge m²/s, flux_Q is zero).
        Q = edges.flux_Q[e] != 0.0 ? edges.flux_Q[e] :
                                      edges.flux[e] * edges.width[e]
        absQ = abs(Q)
        if Q < 0.0
            sum_u[ci] += absQ *  ux;  sum_v[ci] += absQ *  uy
            sum_u[cj] += absQ * -ux;  sum_v[cj] += absQ * -uy
        else
            sum_u[ci] += absQ * -ux;  sum_v[ci] += absQ * -uy
            sum_u[cj] += absQ *  ux;  sum_v[cj] += absQ *  uy
        end
    end

    # Per-cell normalisation (parallel — each i is independent)
    Threads.@threads for i in 1:n
        h = state.water_depth[i]
        if h < h_min
            state.velocity[i] = 0.0
            state.vel_u[i]    = 0.0
            state.vel_v[i]    = 0.0
            continue
        end
        denom = state.cell_area[i] * h
        u = sum_u[i] / denom
        v = sum_v[i] / denom
        state.vel_u[i]    = u
        state.vel_v[i]    = v
        state.velocity[i] = sqrt(u*u + v*v)
    end
end

# Debug counter for step_standard! — counts calls, logs first 2
const _step_debug_count = Ref{Int}(0)

"""
    step_standard!(state, dt)

Advance the standard diffusive-wave model by one timestep `dt`.

## Parallel two-phase design

The combined read-compute-scatter loop is split into three phases so that
the expensive flux computation is embarrassingly parallel while the
write-hazard (multiple edges sharing a cell) is confined to a cheap serial
scatter pass.

**Phase A — edge flux computation** (`Threads.@threads`, one task per edge):
  - Reads `state.volume`, `elevation`, `cell_area`, `manning_n` — all
    read-only during this phase; no writes to cell state.
  - Writes only to `edge_vol[e]` and `edges.flux[e]`, each indexed by `e`
    (unique per thread — no data race).
  - `_bates_flux` is a pure function with no side effects.

**Phase B — dV scatter with donor limiter** (serial):
  - Applies the Bug-46 per-edge donor limit using the start-of-step volume
    (still unmodified) then accumulates into `dV[ci]` / `dV[cj]`.
  - Must be serial: `ci` and `cj` are not unique across edges, so
    concurrent `+=` would race.  Cost is O(n_edges) additions — negligible.

**Phase C — cell volume update** (`Threads.@threads`, one task per cell):
  - Each thread works on a unique index `i` — no data race.

## GPU roadmap
When `SOLVER_BACKEND[] == :gpu` (Phase 5):
  - Phase A → CUDA kernel, one thread per edge.
  - Phase B → `CUDA.@atomic` scatter-add kernel, or prefix-sum reduce.
  - Phase C → CUDA kernel, one thread per cell.
  The interface between phases (edge_vol, dV vectors) is the same for
  both CPU and GPU paths; only the dispatch changes.
"""
function step_standard!(state::FlowState, dt::Float64)
    n     = length(state.cell_ids)
    edges = state.edges
    ne    = edges.n_edges

    # ── Debug: log NaN/depth state before flux loop (first 2 calls) ───────
    if _step_debug_count[] < 2
        _step_debug_count[] += 1
        k = _step_debug_count[]
        n_nan_vol   = count(isnan, state.volume)
        n_nan_elev  = count(isnan, state.elevation)
        n_nan_area  = count(isnan, state.cell_area)
        n_zero_area = count(a -> a < 1.0, state.cell_area)
        vol_sum     = sum(v for v in state.volume if isfinite(v); init=0.0)
        max_depth   = maximum(state.water_depth; init=0.0)
        @info "step_standard! call $k: NaN_vol=$n_nan_vol  NaN_elev=$n_nan_elev" *
              " NaN_area=$n_nan_area  zero_area=$n_zero_area" *
              " vol_sum=$(round(vol_sum,sigdigits=4))" *
              "  max_water_depth=$(round(max_depth,sigdigits=4))"
        if n_nan_vol > 0
            bad = findall(isnan, state.volume)
            @warn "  NaN volumes at indices: $(bad[1:min(5,end)])"
        end
        for i in [1, 2, 3, argmax(state.volume)]
            wse = state.elevation[i] + state.volume[i] / max(state.cell_area[i], 1.0)
            @info "  cell[$i]: vol=$(round(state.volume[i],sigdigits=4))" *
                  "  depth=$(round(state.water_depth[i],sigdigits=4))" *
                  "  wse=$(round(wse,sigdigits=6))"
        end
    end

    # ── Phase A: parallel edge flux computation ─────────────────────────────────────
    # Each iteration reads cell state (read-only) and writes to unique edge
    # slots edge_vol[e] and edges.flux[e].  No two threads share a write
    # target — this is provably race-free.
    #
    # _bates_flux_limited applies three stability fixes vs the earlier call:
    #   • Froude limiter      (CAESAR froude_limit=0.8: caps supercritical q)
    #   • Volume limiter      (CAESAR depth/5 cap: ≤ 1/5 donor depth per edge)
    #   • Consistent q_prev  (stores post-limiting q, not raw Bates q)
    # See _bates_flux_limited docstring for full rationale.
    edge_vol = Vector{Float64}(undef, ne)

    Threads.@threads for e in 1:ne
        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        # Degenerate edge guard (NaN in geometry or elevation).
        if (isnan(edges.width[e])      || isnan(edges.L[e])          ||
            isnan(edges.cos_theta[e])  || isnan(edges.sill[e])       ||
            isnan(state.elevation[ci]) || isnan(state.elevation[cj]))
            edge_vol[e] = 0.0
            edges.flux[e] = 0.0   # clear stale momentum on degenerate edges
            continue
        end

        wse_ci = state.elevation[ci] + state.volume[ci] / max(state.cell_area[ci], 1.0)
        wse_cj = state.elevation[cj] + state.volume[cj] / max(state.cell_area[cj], 1.0)

        # Identify donor depth for the volume limiter (higher-WSE side donates).
        depth_donor = wse_ci >= wse_cj ? state.water_depth[ci] : state.water_depth[cj]

        Q, q_stored = _bates_flux_limited(
            edges.flux[e], wse_ci, wse_cj, edges.sill[e],
            edges.width[e], edges.L[e], edges.cos_theta[e],
            min(state.manning_n[ci], state.manning_n[cj]), dt,
            depth_donor)

        edges.flux[e] = q_stored   # Fix C: post-limiting q → consistent q_prev
        edge_vol[e]   = Q * dt     # signed volume (m³)
    end

    # ── Phase B: serial dV scatter ─────────────────────────────────────────────
    # _bates_flux_limited already bounds the transfer via the Froude and volume
    # limiters.  The DONOR_EDGE_DIVISOR cap is retained as a last-resort
    # mass-conservation guard (should rarely trigger).
    # state.volume[] is still the start-of-step snapshot (Phase A did not
    # modify it), so the donor-limit comparison is consistent.
    dV = zeros(Float64, n)
    @inbounds for e in 1:ne
        ev = edge_vol[e]
        iszero(ev) && continue

        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        if ev > 0.0
            ev = min(ev,  state.volume[cj] / DONOR_EDGE_DIVISOR)   # cj is donor
        else
            ev = max(ev, -state.volume[ci] / DONOR_EDGE_DIVISOR)   # ci is donor
        end

        dV[ci] += ev   # ci gains when ev > 0
        dV[cj] -= ev   # cj loses when ev > 0
    end

    # ── Phase C: parallel cell volume update ──────────────────────────────
    _apply_dV_standard!(state, dV)

    # ── Phase D: velocity ─────────────────────────────────────────────────
    _compute_velocity!(state)
end

"""
    step_sgs!(state, dt)

Advance the SGS diffusive-wave model by one timestep `dt`.

Uses the same parallel A/B/C/D phase structure as `step_standard!` — see
that function's docstring for the thread-safety rationale.  The SGS-specific
differences are:
  - WSE is derived from the hypsometric lookup (not elevation + depth).
  - The edge sill is the pre-computed SGS minimum along the shared boundary.
  - Phase A reads the pre-computed `wse[]` array (computed serially in
    the SGS-specific Step 0) rather than computing WSE inline.

Stability fixes applied in this function (branch: sgs_flow_fixes):

**Fix A — Froude limiter**
  After the Bates eq. 9 call, |q| is capped at h_flow_eff × √(g × h_flow_eff) × FROUDE_LIMIT.
  h_flow_eff is the flow depth that _bates_flux saw (= max(wse_i_eff, wse_j_eff) - z_sill_eff),
  which is already bounded by max(depth_ci, depth_cj) via the Bug 49 z_sill_eff correction.
  LISFLOOD-FP SGC achieves equivalent stability through its R-A formulation; the explicit
  Froude cap is the appropriate mechanism given FloodA5's Bates-based SGS kernel.

**Fix B — Volume limiter (added for high-resolution stability)**
  Caps |Q×dt| at `depth_donor × width / 5`, preventing more than ~20% of the donor
  cell's water column from leaving via one edge per step.  Mirrors `_bates_flux_limited`
  Fix B used in `step_standard!`.  At high resolution (small dx, small dt) the Froude cap
  alone is insufficient — the momentum term drives checkerboard oscillation between adjacent
  cells.  `state.water_depth` is capped at `(z_max - z_min)` in `_apply_dV_sgs!` (Bug 58),
  making it a bounded physical quantity safe to use as the depth-based limiter scale.

**Fix C — Consistent q_stored (primary + full)**
  Primary: edges.flux[e] is set to the post-Froude q_stored (not the raw Bates q).
  Full: if the Phase B DONOR_EDGE_DIVISOR cap further clips the volume, edges.flux[e]
  is updated again to match, so q_prev next step always reflects what was transferred.
  This eliminates the divergence between stored momentum and actual hydraulic state
  that drove SGS oscillations analogously to the standard-flow checkerboard instability.

**Manning's n — arithmetic mean**
  Changed from min(n_ci, n_cj) to 0.5*(n_ci + n_cj) per edge, matching LISFLOOD-FP
  and the standard shallow-water literature.  The minimum is marginally non-standard
  and slightly over-conductive.
"""
function step_sgs!(state::FlowState, dt::Float64)
    n     = length(state.cell_ids)
    edges = state.edges
    ne    = edges.n_edges

    # ── Step 0: WSE and wetted area from hypsometric lookup (serial) ──────
    # Must be serial: wse_from_volume / wetted_area_from_wse each read a
    # cell's SGS table — no shared writes, but the lookup is not @inline
    # and benefits from sequential cache access patterns.
    # Can be made @threads if profiling shows it is a bottleneck.
    wse   = Vector{Float64}(undef, n)
    A_wet = Vector{Float64}(undef, n)
    for i in 1:n
        tbl      = state.sgs_tables[i]
        wse[i]   = wse_from_volume(tbl, state.volume[i])
        A_wet[i] = wetted_area_from_wse(tbl, wse[i])
    end

    # ── Phase A: parallel edge flux computation ────────────────────────────
    # Reads: wse[] (just computed, read-only), state.elevation, manning_n,
    #        edges geometry.  Writes only to edge_vol[e] and edges.flux[e]
    #        (unique per thread — no data race).
    #
    # Dry-cell WSE correction (Bug 48):
    # wse_from_volume(V=0) returns tbl.z_min (the cell thalweg elevation).
    # For a high-elevation dry cell, z_min can be >> z_sill of an adjacent edge,
    # creating a spurious driving head that drives large flux from a dry cell.
    # The fix: for flux computation only, cap a dry cell's effective WSE at z_sill,
    # so it contributes zero head above the sill — matching the standard solver's
    # behaviour where a dry cell never drives flux.
    edge_vol = Vector{Float64}(undef, ne)

    Threads.@threads for e in 1:ne
        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        (isnan(state.elevation[ci]) || isnan(state.elevation[cj])) &&
            (edge_vol[e] = 0.0; continue)

        z_sill = edges.sill[e]
        if isnan(z_sill)
            z_sill = min(state.sgs_tables[ci].z_min, state.sgs_tables[cj].z_min)
        end

        if (isnan(edges.width[e]) || isnan(edges.L[e]) || isnan(edges.cos_theta[e]))
            edge_vol[e] = 0.0
            continue
        end

        # Bug 48 / Bug 49 refinement: effective WSE for dry cells.
        # A dry cell (volume=0) contributes zero pooling head.  Its effective WSE
        # must clear TWO physical barriers before it can receive water:
        #   1. z_sill   — the channel bed (lowest DEM along the shared edge arc)
        #   2. tbl.z_min — the basin floor (lowest terrain in the receiving cell)
        wse_ci_eff = state.volume[ci] > 0.0 ? wse[ci] :
                     max(z_sill, state.sgs_tables[ci].z_min)
        wse_cj_eff = state.volume[cj] > 0.0 ? wse[cj] :
                     max(z_sill, state.sgs_tables[cj].z_min)

        # If both cells are dry, skip
        if state.volume[ci] <= 0.0 && state.volume[cj] <= 0.0
            edge_vol[e] = 0.0
            edges.flux_Q[e] = 0.0
            continue
        end

        # ── R-A flux kernel ──────────────────────────────────────────────────
        # Look up cross-sectional flow area A and hydraulic radius R from the
        # pre-computed edge hydraulic curves in each cell's SGSTable.
        # Both perspectives (ci and cj) are averaged for a symmetric treatment.
        # The adjacency slot is the position of cj in ci's neighbour list,
        # looked up via adj_matrix (O(5) per edge, negligible).
        wse_flow = max(wse_ci_eff, wse_cj_eff)

        slot_i = _adjacency_slot(state.adj_matrix, ci, cj)
        slot_j = _adjacency_slot(state.adj_matrix, cj, ci)
        A_i = A5Grid.flow_area_from_wse(state.sgs_tables[ci], slot_i, wse_flow)
        R_i = A5Grid.hydraulic_radius_from_wse(state.sgs_tables[ci], slot_i, wse_flow)
        A_j = A5Grid.flow_area_from_wse(state.sgs_tables[cj], slot_j, wse_flow)
        R_j = A5Grid.hydraulic_radius_from_wse(state.sgs_tables[cj], slot_j, wse_flow)
        A_edge = 0.5 * (A_i + A_j)
        R_edge = 0.5 * (R_i + R_j)

        Q_new = _manning_flux_ra(edges.flux_Q[e], wse_ci_eff, wse_cj_eff, z_sill,
                                  A_edge, R_edge,
                                  edges.L[e], edges.cos_theta[e],
                                  0.5 * (state.manning_n[ci] + state.manning_n[cj]),
                                  dt)

        # Fix C: store Q (m³/s) as momentum state for next step.
        # The R-A form is self-stabilising — no Froude or volume limiter needed here.
        # Phase B DONOR_EDGE_DIVISOR cap is retained as a last-resort mass guard.
        edges.flux_Q[e] = Q_new
        edge_vol[e]     = Q_new * dt
    end

    # ── Phase B: serial dV scatter with per-edge donor limiter ────────────
    dV = zeros(Float64, n)
    @inbounds for e in 1:ne
        ev = edge_vol[e]
        iszero(ev) && continue

        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        # Last-resort per-edge donor cap: V/DONOR_EDGE_DIVISOR = V/10.
        # The R-A form is self-stabilising so this should rarely bind.
        # Fix C (full): back-propagate clipped Q to flux_Q so q_prev next
        # step reflects what was actually transferred.
        if ev > 0.0
            ev_capped = min(ev, state.volume[cj] / DONOR_EDGE_DIVISOR)
            if ev_capped < ev
                edges.flux_Q[e] = ev_capped / dt
            end
            ev = ev_capped
        else
            ev_capped = max(ev, -state.volume[ci] / DONOR_EDGE_DIVISOR)
            if ev_capped > ev
                edges.flux_Q[e] = ev_capped / dt
            end
            ev = ev_capped
        end

        dV[ci] += ev
        dV[cj] -= ev
    end

    # ── Phase C: parallel cell volume update ──────────────────────────────
    _apply_dV_sgs!(state, dV)

    # ── Phase D: velocity ─────────────────────────────────────────────────
    _compute_velocity!(state)
end

# ---------------------------------------------------------------------------
# Wetted-area saturation (for visualisation)
# ---------------------------------------------------------------------------

"""
    saturation_fraction(state) → Vector{Float64}

Return the fractional wetted area (0–1) for each cell.
For SGS this is meaningful sub-grid information.
For standard flow it is 1.0 wherever depth > 0.
"""
function saturation_fraction(state::FlowState)::Vector{Float64}
    n = length(state.cell_ids)
    sat = Vector{Float64}(undef, n)
    if isempty(state.sgs_tables)
        # Standard: binary wet/dry
        for i in 1:n
            sat[i] = state.water_depth[i] > 1e-4 ? 1.0 : 0.0
        end
    else
        for i in 1:n
            tbl    = state.sgs_tables[i]
            w_area = wetted_area_from_wse(tbl, wse_from_volume(tbl, state.volume[i]))
            sat[i] = tbl.cell_area > 0.0 ? clamp(w_area / tbl.cell_area, 0.0, 1.0) : 0.0
        end
    end
    return sat
end

# ---------------------------------------------------------------------------
# Simulation loop
# ---------------------------------------------------------------------------

"""
    InjectionPoint

A fixed-rate point source: water added to the nearest mesh cell at a
constant volumetric flow rate (m³/s) for the duration of the simulation.
"""
struct InjectionPoint
    cell_index :: Int       # index into state.cell_ids
    cell_id    :: String    # hex cell ID (for logging)
    rate_m3s   :: Float64   # volumetric flow rate (m³/s)
    lon        :: Float64   # source longitude (degrees)
    lat        :: Float64   # source latitude (degrees)
end

"""
    RainPoint

A localised rainfall source: water added to the nearest mesh cell at a
rate equivalent to a given rainfall intensity (mm/hr) applied over that
cell's plan area.  Unlike `--rainfall` (which applies to every cell),
`--rainpoint` applies only to the single nearest cell.

`rate_m3s` is pre-computed as `rainfall_mm_hr / 3_600_000 × cell_area_m2`
and stored so the simulation loop is identical to the InjectionPoint path.
"""
struct RainPoint
    cell_index     :: Int       # index into state.cell_ids
    cell_id        :: String    # hex cell ID (for logging)
    rate_m3s       :: Float64   # volumetric flow rate (m³/s) = mm_hr/3.6e6 × area
    lon            :: Float64   # requested longitude (degrees)
    lat            :: Float64   # requested latitude (degrees)
    rainfall_mm_hr :: Float64   # original user input (mm/hr) — for logging
end

"""
    _find_nearest_cell(mesh, lon, lat) → (index, cell_id, dist_m)

Return the index and ID of the mesh cell whose centre is closest to (lon, lat).
Uses Euclidean distance in degree-space (sufficient for small AOIs).
"""
function _find_nearest_cell(mesh::A5Mesh, lon::Float64, lat::Float64)
    best_i    = 1
    best_dist = Inf
    for (i, c) in enumerate(mesh.cells)
        dlon = c.center_lon - lon
        dlat = c.center_lat - lat
        d    = sqrt(dlon^2 + dlat^2)
        if d < best_dist
            best_dist = d
            best_i    = i
        end
    end
    dist_m = best_dist * 111_000.0   # rough degrees → metres
    # Warn if the nearest cell is far from the requested point — likely outside the mesh
    if dist_m > 2000.0
        @warn "_find_nearest_cell: requested point (lon=$lon, lat=$lat) is $(round(dist_m/1000,digits=1)) km from nearest mesh cell. " *
              "Point may be outside the mesh AOI. Nearest cell: $(mesh.cells[best_i].id)"
    end
    return best_i, mesh.cells[best_i].id, dist_m
end

"""
    run_simulation!(state, mesh, sim_duration, dt_max, vis, vis_mode, output,
                    method, rainfall_rate)

Time-stepping loop.  Advances the model for `sim_duration` seconds using
adaptive dt (CFL-limited, capped at `dt_max`).

Dispatches `push_frame!` to the active visualiser every `VIS_INTERVAL` steps.
Writes HDF5 snapshots according to `output.output_interval`.

Arguments
---------
  method           — `StandardFlow()` or `SGSFlow()`
  rainfall_rate    — uniform rainfall in m/s applied to every cell (default 0)
  injection_points — fixed volumetric point sources (m³/s)
  rain_points      — localised rainfall sources (single nearest cell, mm/hr × area)
"""
function run_simulation!(state            :: FlowState,
                         mesh             :: A5Mesh,
                         sim_duration     :: Float64,
                         dt_max           :: Float64,
                         vis,
                         vis_mode         :: Symbol = :none,
                         output           :: SimOutput = SimOutput();
                         method           :: FlowMethod = StandardFlow(),
                         rainfall_rate    :: Float64   = 0.0,
                         injection_points :: Vector{InjectionPoint} = InjectionPoint[],
                         rain_points      :: Vector{RainPoint}      = RainPoint[],
                         sgs_diag                                   = nothing)
    t    = 0.0
    step = 0
    use_sgs = method isa SGSFlow

    # ── Debug: pre-simulation state check ──────────────────────────────────
    @info "Pre-sim check: n_cells=$(length(state.volume))  " *
          "initial_vol_sum=$(sum(state.volume))  " *
          "rainfall_rate=$rainfall_rate  " *
          "n_injection=$(length(injection_points))  n_rainpoints=$(length(rain_points))"
    if !isempty(rain_points)
        for (k, rp) in enumerate(rain_points)
            @info "  rain_point[$k]: idx=$(rp.cell_index)  id=$(rp.cell_id)  " *
                  "rate=$(rp.rate_m3s) m3/s  " *
                  "valid_idx=$(1 <= rp.cell_index <= length(state.volume))"
        end
    end

    while t < sim_duration
        # ── Pause polling: if Makie pause button pressed, sleep until resumed ──
        # Checks vis.paused (Threads.Atomic{Bool}) between steps — safe with GPU
        # because CUDA kernels complete synchronously before this point.
        if vis_mode === :makie && vis !== nothing
            while vis.paused[]
                sleep(0.05)
            end
        end

        # Adaptive dt
        dt = min(_cfl_dt(state, method), dt_max, sim_duration - t)
        dt = max(dt, 0.1)   # floor: 0.1s to prevent infinite loops on dry mesh

        # Apply rainfall source (uniform, before routing).
        # Volume is primary state for both standard and SGS.
        # dV = rate (m/s) × dt (s) × cell_area (m²)
        if rainfall_rate > 0.0
            for i in eachindex(state.cell_ids)
                state.volume[i] += rainfall_rate * dt * state.cell_area[i]
            end
        end

        # Apply injection point sources (fixed volumetric rate m³/s).
        # Volume is primary state for both methods.
        for inj in injection_points
            state.volume[inj.cell_index] += inj.rate_m3s * dt
        end

        # Apply localised rainfall point sources.
        # Each RainPoint pre-stores the effective m³/s for its cell
        # (mm/hr converted to m/s × cell area), so the loop is identical
        # to the injection-point path.
        for rp in rain_points
            state.volume[rp.cell_index] += rp.rate_m3s * dt
        end

        # ── Debug: post-source volume check (first 5 steps + every 10th) ──
        if step <= 5 || step % 10 == 0
            src_vols = isempty(rain_points) ? Float64[] :
                       [state.volume[rp.cell_index] for rp in rain_points]
            wet_now  = count(>(1e-4), state.water_depth)
            non_src_wet = count(i -> state.water_depth[i] > 1e-4 &&
                all(rp.cell_index != i for rp in rain_points), eachindex(state.volume))
            @info "  Step $step: vol_sum=$(round(sum(state.volume),sigdigits=5))  " *
                  "src_vol=$(round.(src_vols,sigdigits=5))  " *
                  "wet=$wet_now (non-src=$non_src_wet)  dt=$dt"
        end

        # Sync water_depth from volume after all sources have been applied.
        # Sources (rainfall, injection, rainpoint) write directly to state.volume
        # but state.water_depth is only updated inside _apply_dV_standard!/sgs!.
        # Without this sync, the first-step flux loop and the progress log both
        # read stale water_depth = 0, giving wet=0 and wrong CFL even though
        # volume is non-zero.
        if !use_sgs
            for i in eachindex(state.cell_ids)
                state.cell_area[i] >= 1.0 &&
                    (state.water_depth[i] = state.volume[i] / state.cell_area[i])
            end
        end
        # SGS: water_depth is derived from the hypsometric curve, which is
        # evaluated inside step_sgs! — no pre-sync needed for SGS.

        # Physics step
        if use_sgs
            step_sgs!(state, dt)
        else
            step_standard!(state, dt)
        end

        t    += dt
        step += 1

        # Visualiser push
        if vis !== nothing && step % VIS_INTERVAL == 0
            sat = saturation_fraction(state)
            if vis_mode === :cesium
                VisualisationServer.push_frame!(vis, t, Dict(
                    "depth"      => Float32.(state.water_depth),
                    "saturation" => Float32.(sat),
                    "volume"     => Float32.(state.volume),
                    "velocity"   => Float32.(state.velocity),
                ))
            elseif vis_mode === :makie
                # Compute cumulative volume budget for the mass-balance plots.
                # vol_added  = all water injected since t=0 (rainfall + point sources).
                # vol_domain = water currently in the domain (primary state sum).
                # vol_removed is 0 until Phase 2 adds open outflow BCs.
                _vis_vol_added = rainfall_rate * t *
                    sum(a for a in state.cell_area if a >= 1.0; init=0.0) +
                    sum(inj.rate_m3s * t for inj in injection_points; init=0.0) +
                    sum(rp.rate_m3s  * t for rp  in rain_points;      init=0.0)
                _vis_vol_domain = sum(state.volume)
                _vis_n_wet      = count(>(1e-4), state.water_depth)
                MakieVisualiser.push_frame!(
                    vis, state.cell_ids,
                    state.water_depth,
                    sat,
                    state.volume,
                    state.velocity,
                    state.vel_u,
                    state.vel_v,
                    t;
                    vol_added   = _vis_vol_added,
                    vol_domain  = _vis_vol_domain,
                    vol_removed = 0.0,
                    n_wet       = _vis_n_wet,
                    sim_step    = step,
                    sim_dt      = dt)
                # Update SGS diagnostic window if open
                if sgs_diag !== nothing
                    _update_sgs_diagnostic!(sgs_diag, state, t)
                    yield()   # let GLMakie render thread redraw the diagnostic window
                end
            end
        end

        # HDF5 output
        if _should_write_frame(output, t)
            _write_frame!(output, state, t)
        end

        if step % 50 == 0
            n_wet     = count(>(1e-4), state.water_depth)
            max_depth = isempty(state.water_depth) ? 0.0 :
                        something(maximum(v for v in state.water_depth if isfinite(v); init=0.0), 0.0)
            # Mass balance: cumulative rainfall input vs total domain volume.
            # With the per-edge donor limiter (Bug 46 fix), flux is mass-conserving,
            # so mb_err should be ~0 for a closed domain with no outflow BCs.
            # Small residuals (<< 1 m³) from the max(0,…) floor are acceptable.
            # mb_err < 0 means domain_vol > input_vol (mass creation — should not occur).
            # mb_err > 0 means volume has left the domain (open boundaries or floor clips).
            input_vol  = rainfall_rate * t *
                sum(a for a in state.cell_area if a >= 1.0; init=0.0) +
                sum(inj.rate_m3s * t for inj in injection_points; init=0.0) +
                sum(rp.rate_m3s  * t for rp  in rain_points;      init=0.0)
            domain_vol = sum(state.volume)
            mb_err     = input_vol - domain_vol
            @info "  step=$(lpad(step,5))  t=$(round(t,digits=1))s  dt=$(round(dt,digits=2))s  wet=$n_wet  max_depth=$(round(max_depth,digits=3))m  domain_vol=$(round(domain_vol,digits=1))m³  mb_err=$(round(mb_err,digits=1))m³"
        end
    end
    @info "Simulation finished at t=$(round(t,digits=1))s  ($(step) steps)"

    if output.enabled && t > output.last_write_t
        _write_frame!(output, state, t)
    end
    output.enabled &&
        @info "HDF5 output: $(output.frame_count) frames → $(output.path)"
end

# ---------------------------------------------------------------------------
# Application entry point
# ---------------------------------------------------------------------------

"""
    run_flood_model(; mesh_source, vis_mode, vis_port, dem_source, dem_strict,
                      dem_method, dem_samples, halton_seed,
                      flow_method, sgs_bins, manning_n, friction_source,
                      sim_duration, dt_max, rainfall_rate,
                      output_path, output_interval)

Main application entry point.

`mesh_source` is one of:
  - `(:generate, aoi_path, resolution, output_path_or_nothing)`
  - `(:load, parquet_path)`

`dem_source` is one of:
  - `nothing`              — no DEM; elevation stays zero
  - `FileDEM(path)`        — sample from a local GeoTIFF

`flow_method` is one of:
  - `:standard`  — diffusive wave on mean cell elevation (fast)
  - `:sgs`       — diffusive wave with sub-grid hypsometric lookup (accurate)
"""
function run_flood_model(;
    mesh_source      :: Tuple,
    mesh_only        :: Bool    = false,
    vis_mode         :: Symbol  = :none,
    vis_port         :: Int     = 8080,
    dem_source                  = nothing,
    dem_strict       :: Bool    = false,
    dem_method       :: Symbol  = :mean,
    dem_samples      :: Int     = 256,
    halton_seed      :: Int     = 0,
    flow_method      :: Symbol  = :sgs,
    sgs_bins         :: Int     = 100,
    sgs_samples      :: Int     = 512,
    manning_n        :: Float64 = 0.03,
    friction_source             = nothing,
    sim_duration      :: Float64 = 3600.0,
    dt_max            :: Float64 = 60.0,
    rainfall_rate     :: Float64 = 0.0,
    injection_specs   :: Vector{Tuple{Float64,Float64,Float64}} = Tuple{Float64,Float64,Float64}[],
    rainpoint_specs   :: Vector{Tuple{Float64,Float64,Float64}} = Tuple{Float64,Float64,Float64}[],
    output_path       :: Union{String,Nothing} = nothing,
    output_interval   :: Float64 = 60.0)

    @info "=== A5 Flood Model ===" Dates.now()
    @info "Vis mode    : $vis_mode"
    if !mesh_only
        @info "Flow method : $flow_method"
    end

    # 1. Start Cesium server early
    vis      = if vis_mode === :cesium
        VisualisationServer.start(port = vis_port,
                                  viz_dir = joinpath(@__DIR__, "visualisation", "cesium"))
    else
        nothing
    end
    sgs_diag = nothing   # SGS hypsometric diagnostic window; assigned later if applicable
    vis_mode === :cesium &&
        @info "Cesium viewer: open http://localhost:$vis_port in your browser"

    # 2. Acquire mesh
    mesh = if mesh_source[1] === :generate
        _, aoi_path, resolution, mesh_out = mesh_source
        @info "AOI        : $aoi_path"
        @info "Resolution : $resolution"
        @info "Mesh out   : $(mesh_out === nothing ? "(not saved)" : mesh_out)"
        @info "Generating A5 pentagon mesh..."
        t0 = time()
        m = if mesh_out !== nothing
            fmt = endswith(mesh_out, ".parquet") ? :geoparquet : :geojson
            mesh_for_aoi(aoi_path, resolution; output_path=mesh_out, format=fmt)
        else
            mesh_for_aoi(aoi_path, resolution)
        end
        @info "Mesh generated in $(round(time()-t0, digits=1))s — $(length(m)) cells"
        m
    else
        _, parquet_path, mesh_out = mesh_source
        @info "Loading mesh from $parquet_path..."
        t0 = time()
        m = load_mesh_geoparquet(parquet_path)
        @info "Mesh loaded in $(round(time()-t0, digits=1))s — $(length(m)) cells"

        # ── Compatibility check ──────────────────────────────────────────
        if flow_method === :sgs && !haskey(m.array_vars, "sgs_elev_bins")
            error(
                "Mesh at '$parquet_path' has no SGS hypsometric tables, " *
                "but --flow-model sgs was requested.\n" *
                "  Re-generate the mesh with:\n" *
                "    --meshgen <aoi.geojson> --meshres <N> --dem <dem.tif> " *
                "--flow-model sgs --mesh-only --meshout <mesh.parquet>"
            )
        end
        if flow_method !== :sgs && haskey(m.array_vars, "sgs_elev_bins")
            @info "Note: mesh contains SGS tables but --flow-model standard was " *
                  "requested — SGS data will be ignored."
        end
        if isempty(m.adjacency)
            @warn "Mesh has no pre-computed adjacency (old format). " *
                  "Adjacency will be computed from boundary geometry at initialisation."
        end
        m
    end

    # 3. DEM sampling (if requested)
    if dem_source !== nothing
        if mesh_source[1] === :load && haskey(mesh.static_vars, "elevation")
            n_valid = count(isfinite, mesh.static_vars["elevation"])
            n_total = length(mesh.cells)
            error(
                "The loaded mesh already contains elevation data " *
                "($n_valid / $n_total cells with valid values).\n" *
                "  To re-sample, first remove the elevation column and re-save."
            )
        end

        @info "Checking DEM coverage..."
        cov = check_dem_coverage(mesh, dem_source)
        @info "DEM coverage: $(round(cov.coverage_pct, digits=1))% " *
              "($(cov.cells_covered)/$(cov.total_cells) cells), CRS: $(cov.dem_crs)"
        if cov.cells_outside > 0
            msg = "$(cov.cells_outside) cell(s) outside DEM extent"
            dem_strict ? error("$msg.") : @warn "$msg — NaN assigned."
        end

        @info "Sampling DEM onto mesh (method=:$(dem_method))..."
        t0 = time()
        if dem_method === :mean
            sample_dem_mean!(mesh, dem_source; strict=dem_strict,
                             n_samples=dem_samples, halton_seed=halton_seed)
        elseif dem_method === :centroid
            sample_dem_centroid!(mesh, dem_source; strict=dem_strict)
        else
            error("Unknown dem_method :$dem_method")
        end
        @info "DEM sampled in $(round(time()-t0, digits=1))s"

        mesh_save_path = if mesh_source[1] === :generate
            mesh_source[4]
        else
            something(mesh_source[3], mesh_source[2])
        end
        if mesh_save_path !== nothing && endswith(mesh_save_path, ".parquet")
            @info "Re-saving mesh with elevation to $mesh_save_path..."
            save_mesh_geoparquet(mesh, mesh_save_path)
        end
    else
        if !haskey(mesh.static_vars, "elevation")
            @warn "No DEM provided and mesh has no saved elevation data. " *
                  "The model will run with flat terrain (z=0 everywhere). " *
                  "Use --dem <file.tif> to sample elevation data."
        end
    end

    # 4. Log mesh summary (before run guard — always shown)
    @info mesh_summary(mesh)

    # 4b. Run guard (first pass): if mesh_only but NOT sgs, exit now.
    # If mesh_only AND sgs, we continue to step 5 to build SGS tables first,
    # then exit after saving. If no water source and not mesh_only, also exit.
    has_water = rainfall_rate > 0.0 ||
                !isempty(injection_specs) ||
                !isempty(rainpoint_specs)   # extend here for inflow/BC in Phase 2
    if !mesh_only && !has_water
        @info "No water source provided (rainfall=0, no inflow). " *
              "Mesh saved and ready. Re-run with --rainfall <mm/hr> to simulate."
        return nothing
    end
    if mesh_only && flow_method !== :sgs
        @info "Mesh-only run complete (--mesh-only flag set). Exiting without simulation."
        return nothing
    end
    # mesh_only+sgs: fall through to SGS pre-processing, then exit after.

    # 5. SGS pre-processing — runs when flow_method === :sgs, regardless of
    #    mesh_only flag.  Rationale: a mesh built with --flow-model sgs should
    #    be fully ready for SGS simulation (including hypsometric tables) so the
    #    user can inspect and validate the mesh before running the model.
    #    A mesh built without --flow-model sgs will fail a compatibility check
    #    at --meshload time if sgs is requested (see step 2 load guard).
    if flow_method === :sgs
        if !haskey(mesh.array_vars, "sgs_elev_bins")
            dem_source !== nothing ||
                error("SGS flow requires a DEM — provide --dem <file.tif>.")
            @info "Building SGS tables ($sgs_bins bins, $sgs_samples pts/cell)..."
            t0 = time()
            build_sgs_tables!(mesh, dem_source;
                              n_bins=sgs_bins, n_samples=sgs_samples,
                              halton_seed=halton_seed)
            @info "SGS tables built in $(round(time()-t0, digits=1))s"

            mesh_save_path = mesh_source[1] === :generate ?
                mesh_source[4] : something(mesh_source[3], mesh_source[2])
            if mesh_save_path !== nothing && endswith(mesh_save_path, ".parquet")
                @info "Re-saving mesh with SGS tables to $mesh_save_path..."
                save_mesh_geoparquet(mesh, mesh_save_path)
            end
        else
            n_bins_loaded = Int(mesh.static_vars["sgs_n_bins"][1])
            @info "SGS tables loaded from mesh ($n_bins_loaded bins/cell)"
        end
    end

    # 5c. Deferred mesh-only exit (after SGS tables built and saved)
    if mesh_only
        @info "Mesh-only run complete (--mesh-only --flow-model sgs). " *
              "SGS tables built and saved. Exiting without simulation."
        return nothing
    end

    # 6. Hand mesh to Cesium visualiser.
    # Makie is deferred to step 7b so source cell indices can be highlighted.
    if vis_mode === :cesium
        @info "Pushing mesh to Cesium server ($(length(mesh)) cells)..."
        VisualisationServer.set_mesh!(vis, mesh_to_geojson_string(mesh),
                                    [c.id for c in mesh.cells])
        @info "Mesh pushed — viewer ready"
    end

    # 7. Initialise flow model
    method_obj = flow_method === :sgs ? SGSFlow() : StandardFlow()
    @info "Initialising flow model ($(flow_method)) on $(length(mesh)) cells..."
    t0 = time()
    flow_state = initialise_flow_model(mesh, method_obj;
                                        manning_n       = manning_n,
                                        friction_raster = friction_source)
    @info "Flow model ready in $(round(time()-t0, digits=1))s  " *
          "(adjacency: $(length(flow_state.adjacency)) cells, " *
          "Manning n = $(manning_n))"

    # Resolve injection point specs (lon, lat, rate_m3s) → InjectionPoint structs
    injection_points = InjectionPoint[]
    for (lon, lat, rate) in injection_specs
        idx, cid, dist_m = _find_nearest_cell(mesh, lon, lat)
        push!(injection_points, InjectionPoint(idx, cid, rate, lon, lat))
        @info "Injection point: ($(round(lon,digits=5)), $(round(lat,digits=5))) → cell $cid  (dist=$(round(dist_m,digits=0))m)  rate=$(round(rate,digits=4)) m³/s"
    end

    # Resolve rainpoint specs (lon, lat, mm_hr) → RainPoint structs.
    # Rate is converted from mm/hr to m/s, then multiplied by the cell area
    # so the simulation loop treats it identically to an InjectionPoint.
    rain_points = RainPoint[]
    for (lon, lat, mm_hr) in rainpoint_specs
        idx, cid, dist_m = _find_nearest_cell(mesh, lon, lat)
        area_m2  = flow_state.cell_area[idx]
        rate_m3s = (mm_hr / 3_600_000.0) * area_m2
        push!(rain_points, RainPoint(idx, cid, rate_m3s, lon, lat, mm_hr))
        @info "Rain point: ($(round(lon,digits=5)), $(round(lat,digits=5))) → cell $cid  " *
              "(dist=$(round(dist_m,digits=0))m)  $(round(mm_hr,digits=2)) mm/hr  " *
              "= $(round(rate_m3s, sigdigits=4)) m³/s  (cell area $(round(area_m2,digits=0)) m²)"
        @info "  RainPoint debug: idx=$idx  n_cells=$(length(flow_state.cell_area))  " *
              "cell_area[idx]=$(flow_state.cell_area[idx])  " *
              "isnan(area)=$(isnan(area_m2))  isnan(rate)=$(isnan(rate_m3s))"
    end

    # 7b. Open Makie viewer now that source indices and the flow adjacency dict
    # are both known.  Passing adjacency activates ring mode: a BFS from the
    # source cells assigns each cell a ring index, enabling the "Ring index"
    # map overlay and the per-ring volume bar chart in the bottom strip.
    if vis_mode === :makie
        src_label      = basename(mesh_source[2])
        source_indices = vcat(
            [inj.cell_index for inj in injection_points],
            [rp.cell_index  for rp  in rain_points],
        )
        @info "Opening Makie viewer ($(length(mesh)) cells, " *
              "$(length(source_indices)) source cell(s), ring mode active)..."
        vis = MakieVisualiser.start(mesh;
                  source_indices = source_indices,
                  adjacency      = flow_state.adjacency,
                  title          = "FloodA5  res=$(mesh.resolution)  $src_label")

        # SGS diagnostic window — opened below if applicable
        if flow_method === :sgs && !isempty(source_indices)
            sgs_diag = _open_sgs_diagnostic(flow_state, source_indices[1])
        end
    end

    # 8. Prepare HDF5 output
    sim_output = SimOutput(
        path     = something(output_path, ""),
        interval = output_interval,
        enabled  = output_path !== nothing,
    )
    sim_output.enabled && _write_mesh_metadata!(sim_output, mesh)

    # 9. Simulation loop
    @info "Starting simulation: duration=$(sim_duration)s  dt_max=$(dt_max)s  " *
          "rainfall=$(round(rainfall_rate*3600*1000, digits=2)) mm/hr"
    run_simulation!(flow_state, mesh, sim_duration, dt_max,
                    vis, vis_mode, sim_output;
                    method           = method_obj,
                    rainfall_rate    = rainfall_rate,
                    injection_points = injection_points,
                    rain_points      = rain_points,
                    sgs_diag         = sgs_diag)

    # 10. Notify visualiser that the simulation is done, then keep alive for replay
    if vis !== nothing && vis_mode === :cesium
        VisualisationServer.notify_complete!(vis)
    end

    if vis !== nothing
        try
            if vis_mode === :cesium
                @info "Simulation complete. Cesium viewer at http://localhost:$vis_port"
                @info "Press Ctrl-C to stop."
                while !istaskdone(vis.task)
                    sleep(1)
                end
            elseif vis_mode === :makie
                @info "Simulation complete. Close the Makie window or press Ctrl-C."
                while vis.running[]
                    sleep(1)
                end
            end
        catch e
            e isa InterruptException || rethrow(e)
        finally
            @info "Shutting down visualiser..."
            if vis_mode === :cesium
                VisualisationServer.stop(vis)
                sleep(0.5)
            elseif vis_mode === :makie
                MakieVisualiser.stop(vis)
            end
        end
    end

    return mesh, flow_state
end

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

function print_help(exit_code::Int = 0)
    modes_str = join(string.(VIS_MODES), " | ")
    println("""
FloodA5 — A5 Pentagon Flood Model
==================================
Usage:
  julia [--threads auto] FloodModel.jl  --meshgen <aoi.geojson>  --meshres <N>
                                        [--meshout <file>]
                                        [--dem <file.tif>]  [--dem-strict]
                                        [--vis [mode]]  [--vis-port PORT]
                                        [--output <file.h5>]
                                        [--output-interval SECS]

  julia [--threads auto] FloodModel.jl  --meshload <mesh.parquet>
                                        [--dem <file.tif>]  [--dem-strict]
                                        [--vis [mode]]  [--vis-port PORT]
                                        [--output <file.h5>]
                                        [--output-interval SECS]

  julia FloodModel.jl --help

Mesh options (choose one):
  --meshgen  FILE    GeoJSON area of interest. Triggers mesh generation via pya5.
  --meshres  N       A5 resolution level (required with --meshgen).
  --mesh-only        Generate (or load) and save the mesh, then exit without
                     simulating. Implied when no water source is provided.
  --meshout  FILE    Save generated mesh (optional). Extension controls format:
                       .parquet  → GeoParquet  (recommended)
                       .geojson  → GeoJSON
  --meshload FILE    Load a previously saved GeoParquet mesh.

DEM options:
  --dem      FILE    GeoTIFF elevation file. Samples DEM and bakes elevation
                     into the mesh parquet for future reuse. When used with
                     --meshload: errors if elevation is already present
                     (prevents silent overwrites).
  --dem-strict       Error if any sample point falls outside the DEM extent.
                     Default: assign NaN and continue with a warning.
  --dem-method M     Sampling method: mean (default) or centroid.
                       mean      arithmetic mean of Halton-distributed points
                                 within the cell polygon (better physics).
                       centroid  single bilinear sample at cell centre (fast).
  --dem-samples N    Halton candidate points per cell for --dem-method mean
                     (default: 256).
  --dem-seed N       Halton sequence offset (default: 0). Use different values
                     for independent Monte Carlo uncertainty runs.

Flow model options:
  --flow-model M     Flow routing method: sgs (default) or standard.
                       sgs       Sub-Grid Sampling — hypsometric volume lookup.
                                 Captures partial wetting, ditches, channels.
                                 Requires DEM for pre-processing (auto-run).
                       standard  Diffusive wave on mean cell elevation (fast).
  --sgs-bins N       Elevation bins for SGS hypsometric curves (default: 100).
  --sgs-samples N    Halton points per cell for SGS pre-processing (default: 512).
  --manning-n N      Global Manning's roughness coefficient (default: 0.03).
  --friction FILE    GeoTIFF friction raster. Per-cell Manning's n sampled at
                     cell centres; overrides --manning-n where finite.
  --sim-duration S   Simulation duration in seconds (default: 3600).
  --dt-max S         Maximum adaptive timestep in seconds (default: 60).
  --rainfall R       Uniform rainfall rate in mm/hr (default: 0).
  --injection-point LAT,LON,RATE
                     Add a point source at (lat, lon) with volumetric flow rate
                     RATE (m³/s). Repeatable for multiple sources.
                     Example: --injection-point -43.386,172.648,0.5
  --rainpoint LAT,LON,RATE_MM_HR
                     Add a localised rainfall source at (lat, lon) with intensity
                     RATE_MM_HR (mm/hr) applied to the area of the single nearest
                     cell only. Unlike --rainfall (which is applied to every cell),
                     this injects water into one cell — useful for point-source
                     testing and validation on flat meshes. Repeatable.
                     Example: --rainpoint -43.531,172.636,50.0

Visualisation options:
  --vis [MODE]       Enable visualisation (off by default).
                     Available modes: $modes_str
                     Omitting MODE defaults to 'cesium'.
  --vis-port PORT    Port for the Cesium HTTP server (default: 8080).

Output options:
  --output FILE      Write HDF5 simulation output (default: off).
  --output-interval  Seconds of simulation time between output snapshots
                     (default: 60.0).

Examples:
  # Generate mesh, sample DEM, build SGS tables, run 1hr SGS simulation
  julia --threads auto FloodModel.jl \\
      --meshgen christchurch_aoi.geojson --meshres 14 --meshout mesh_sgs.parquet \\
      --dem linz_dem.tif --mesh-only

  # Reload mesh with pre-built SGS tables, run with friction raster
  julia --threads auto FloodModel.jl \\
      --meshload mesh_sgs.parquet --friction land_use_n.tif \\
      --rainfall 10 --sim-duration 7200 --output sim_out.h5

  # Standard (non-SGS) run for quick testing
  julia --threads auto FloodModel.jl \\
      --meshload mesh_sgs.parquet --flow-model standard --sim-duration 1800

Resolution guide (approximate cell area — equal-area globally):
  Level  5  ~33,100 km²  Continental / regional
  Level  8  ~518 km²     Large catchment
  Level 10  ~32 km²      Medium catchment
  Level 12  ~2.02 km²    Small catchment
  Level 14  ~12.6 ha     Urban / detailed
  Level 17  ~1,976 m²    High-resolution modelling
  Level 18  ~494 m²      High-resolution channel
""")
    exit(exit_code)
end

function _pop_flag(args::Vector{String}, flag::String)
    idx = findfirst(==(flag), args)
    idx === nothing && return (nothing, args)
    val = get(args, idx + 1, nothing)
    new_args = deleteat!(copy(args), val !== nothing ? [idx, idx+1] : [idx])
    return (val, new_args)
end

function _pop_bool(args::Vector{String}, flag::String)
    idx = findfirst(==(flag), args)
    idx === nothing && return (false, args)
    return (true, deleteat!(copy(args), idx))
end


# ---------------------------------------------------------------------------
# SGS hypsometric diagnostic window
# ---------------------------------------------------------------------------
# Opens a second GLMakie figure showing depth/area curves for the source cell
# and its immediate neighbours.  Updates every visualisation interval so you
# can watch the water surface rise and fall within each cell's storage curve.
# ---------------------------------------------------------------------------

"""
    SGSDiagnostic

Holds the Observables and axis references for the SGS diagnostic figure so
that `_update_sgs_diagnostic!` can update them without recreating the figure.
"""
struct SGSDiagnostic
    fig          :: Any
    cell_indices :: Vector{Int}
    cell_labels  :: Vector{String}
    wse_obs      :: Vector{Any}    # Observable{Vector{Float64}} — current WSE per cell
    axes         :: Vector{Any}
    time_obs     :: Any            # Observable{String} — sim time label
    label_obs    :: Vector{Any}    # Observable{String} — stats text per cell
end

"""
    _open_sgs_diagnostic(state, source_idx) → SGSDiagnostic | nothing

Open a GLMakie figure showing the hypsometric storage curve (bed area vs
elevation) for the source cell and each of its neighbours.  A horizontal
line on each subplot shows the current water surface elevation.

Returns `nothing` if GLMakie is not loaded or SGS tables are absent.
"""
function _open_sgs_diagnostic(state::FlowState, source_idx::Int)
    isempty(state.sgs_tables) && return nothing
    isdefined(MakieVisualiser, :GLMakie) || return nothing
    GM = MakieVisualiser.GLMakie

    # Collect source + up to 5 neighbours from adj_matrix
    nbr_indices = Int[]
    for slot in 1:5
        j = state.adj_matrix[slot, source_idx]
        j > 0 && push!(nbr_indices, j)
    end
    cell_indices = vcat([source_idx], nbr_indices)
    n_panels     = length(cell_indices)

    # Find the edge sill between source and each neighbour from the EdgeList.
    # edge_sill[k] = sill elevation for the edge between source and cell_indices[k]
    #              = NaN for the source cell itself (k=1) or if edge not found.
    edge_sills = fill(NaN, n_panels)
    for e in 1:state.edges.n_edges
        ci = state.edges.cell_i[e]
        cj = state.edges.cell_j[e]
        for (k, idx) in enumerate(cell_indices)
            k == 1 && continue   # skip source-vs-source
            if (ci == source_idx && cj == idx) || (cj == source_idx && ci == idx)
                edge_sills[k] = state.edges.sill[e]
            end
        end
    end

    cell_labels = String[]
    for (k, ci) in enumerate(cell_indices)
        push!(cell_labels, k == 1 ? "source [$(ci)]" : "nbr$(k-1) [$(ci)]")
    end

    ncols = max(2, ceil(Int, n_panels / 2))
    nrows = ceil(Int, n_panels / ncols)
    fig   = GM.Figure(size=(320*ncols, 300*nrows + 30),
                      title="SGS Hypsometric Diagnostic")

    time_obs = GM.Observable("t = 0.00 h")
    GM.Label(fig[0, 1:ncols], time_obs;
             fontsize=13, halign=:center, tellwidth=false)

    axes      = Any[]
    wse_obs   = Any[]    # Observable{Vector{Float64}} — current WSE per cell
    label_obs = Any[]    # Observable{String} — stats text per cell

    for (k, ci) in enumerate(cell_indices)
        tbl  = state.sgs_tables[ci]
        row  = ceil(Int, k / ncols)
        col  = ((k-1) % ncols) + 1
        ax   = GM.Axis(fig[row, col];
                   title     = cell_labels[k],
                   xlabel    = "Wetted area (m²)",
                   ylabel    = "Elevation (m)",
                   titlesize = 11, xlabelsize = 9, ylabelsize = 9)

        # Static hypsometric curve
        GM.lines!(ax, tbl.area_curve, tbl.elev_bins; color=:steelblue, linewidth=2)
        GM.band!(ax, tbl.area_curve,
                 fill(tbl.z_min, length(tbl.elev_bins)),
                 tbl.elev_bins; color=(:steelblue, 0.15))

        # z_min / z_max reference lines (dashed grey)
        GM.hlines!(ax, [tbl.z_min, tbl.z_max]; color=:gray60, linewidth=1, linestyle=:dash)

        # Edge sill line (dotted green) — shows whether terrain blocks flow from source
        if k > 1 && isfinite(edge_sills[k])
            GM.hlines!(ax, [edge_sills[k]]; color=:forestgreen, linewidth=1.5, linestyle=:dot)
        end

        # Live WSE line (red)
        wse_init = A5Grid.wse_from_volume(tbl, state.volume[ci])
        obs_wse  = GM.Observable([wse_init])
        push!(wse_obs, obs_wse)
        GM.hlines!(ax, obs_wse; color=:tomato, linewidth=2)

        # Text overlay: WSE, volume, depth above z_min, and (for nbrs) edge sill
        depth_init = wse_init - tbl.z_min
        sill_str   = (k > 1 && isfinite(edge_sills[k])) ?
                     "\nsill=$(round(edge_sills[k], digits=1))m" : ""
        obs_txt = GM.Observable(
            "WSE=$(round(wse_init,digits=2))m  d=$(round(depth_init,digits=3))m$(sill_str)")
        push!(label_obs, obs_txt)
        GM.text!(ax, 0.02, 0.97; text=obs_txt, space=:relative,
                 align=(:left, :top), fontsize=8, color=:black)

        push!(axes, ax)
    end

    GM.display(GM.Screen(), fig)
    @info "SGS diagnostic window opened: $(n_panels) cells (source + neighbours)"
    for (k, ci) in enumerate(cell_indices)
        k == 1 && continue
        sill_info = isfinite(edge_sills[k]) ?
            "sill=$(round(edge_sills[k],digits=2))m" : "sill=NaN (edge not in EdgeList!)"
        tbl = state.sgs_tables[ci]
        @info "  $(cell_labels[k]): z_min=$(round(tbl.z_min,digits=2))m  $sill_info"
    end

    return SGSDiagnostic(fig, cell_indices, cell_labels, wse_obs, axes, time_obs, label_obs)
end

"""
    _update_sgs_diagnostic!(diag, state, t)

Update the WSE horizontal lines and time label on the SGS diagnostic figure.
Called every visualisation frame alongside `push_frame!`.
`t` is the current simulation time in seconds.
"""
function _update_sgs_diagnostic!(diag::SGSDiagnostic, state::FlowState, t::Float64)
    diag === nothing && return
    isempty(state.sgs_tables) && return
    for (k, ci) in enumerate(diag.cell_indices)
        ci > length(state.sgs_tables) && continue
        tbl   = state.sgs_tables[ci]
        wse   = A5Grid.wse_from_volume(tbl, state.volume[ci])
        depth = wse - tbl.z_min
        diag.wse_obs[k][] = [wse]
        if k <= length(diag.label_obs)
            diag.label_obs[k][] = "WSE=$(round(wse,digits=2))m  d=$(round(depth,digits=3))m"
        end
    end
    diag.time_obs[] = "t = $(round(t / 3600.0, digits=2)) h"
end

function main(args=String[])
    profile = get(ENV, "DEBUG_PROFILE_NAME", "Production/Shell")
    @debug "Model profile:" profile=profile args=args

    #args = copy(ARGS)
    if isinteractive()
        # Use abspath to avoid confusing the VS Code terminal server
        target_dir = "F:/OneDrive - University of Canterbury/Julia/FloodA5"
        if pwd() != target_dir
            cd(target_dir)
        end
    end
    @info "Arguments passed: $(args)"

    ("--help" in args || "-h" in args) && print_help(0)

    # --vis [mode]
    vis_mode = :none
    vis_idx  = findfirst(==("--vis"), args)
    if vis_idx !== nothing
        next = get(args, vis_idx + 1, "")
        if !isempty(next) && !startswith(next, "-") && Symbol(next) in VIS_MODES
            vis_mode = Symbol(next)
            deleteat!(args, [vis_idx, vis_idx + 1])
        else
            vis_mode = :cesium
            deleteat!(args, vis_idx)
        end
    end

    # --vis-port
    vis_port_val, args = _pop_flag(args, "--vis-port")
    vis_port = vis_port_val !== nothing ? parse(Int, vis_port_val) : 8080

    # --meshgen / --meshres / --meshout / --meshload / --mesh-only
    meshgen_val,  args = _pop_flag(args, "--meshgen")
    meshres_val,  args = _pop_flag(args, "--meshres")
    meshout_val,  args = _pop_flag(args, "--meshout")
    meshload_val, args = _pop_flag(args, "--meshload")
    mesh_only,    args = _pop_bool(args, "--mesh-only")

    # --dem / --dem-strict / --dem-method / --dem-samples / --dem-seed
    dem_val,         args = _pop_flag(args, "--dem")
    dem_strict,      args = _pop_bool(args, "--dem-strict")
    dem_method_val,  args = _pop_flag(args, "--dem-method")
    dem_samples_val, args = _pop_flag(args, "--dem-samples")
    dem_seed_val,    args = _pop_flag(args, "--dem-seed")

    dem_method  = dem_method_val  !== nothing ? Symbol(dem_method_val)      : :mean
    dem_samples = dem_samples_val !== nothing ? parse(Int, dem_samples_val) : 256
    halton_seed = dem_seed_val    !== nothing ? parse(Int, dem_seed_val)    : 0

    dem_method ∉ (:centroid, :mean) &&
        (println("ERROR: --dem-method must be 'centroid' or 'mean'\n"); print_help(1))

    # --flow-model / --sgs-bins / --sgs-samples / --manning-n / --friction
    flow_method_val,  args = _pop_flag(args, "--flow-model")
    sgs_bins_val,     args = _pop_flag(args, "--sgs-bins")
    sgs_samples_val,  args = _pop_flag(args, "--sgs-samples")
    manning_n_val,    args = _pop_flag(args, "--manning-n")
    friction_val,     args = _pop_flag(args, "--friction")

    # Default flow method: standard for mesh-only runs (no simulation needed),
    # sgs for simulation runs (full accuracy by default).
    flow_method  = flow_method_val !== nothing ? Symbol(flow_method_val) :
                   mesh_only ? :standard : :sgs
    sgs_bins     = sgs_bins_val     !== nothing ? parse(Int, sgs_bins_val)      : 100
    sgs_samples  = sgs_samples_val  !== nothing ? parse(Int, sgs_samples_val)   : 512
    manning_n    = manning_n_val    !== nothing ? parse(Float64, manning_n_val) : 0.03

    flow_method ∉ (:sgs, :standard) &&
        (println("ERROR: --flow-model must be 'sgs' or 'standard'\n"); print_help(1))

    # --sim-duration / --dt-max / --rainfall
    sim_dur_val,  args = _pop_flag(args, "--sim-duration")
    dt_max_val,   args = _pop_flag(args, "--dt-max")
    rainfall_val, args = _pop_flag(args, "--rainfall")

    sim_duration  = sim_dur_val  !== nothing ? parse(Float64, sim_dur_val)  : 3600.0
    dt_max        = dt_max_val   !== nothing ? parse(Float64, dt_max_val)   : 60.0
    # Rainfall: user gives mm/hr, convert to m/s
    rainfall_rate = rainfall_val !== nothing ?
                    parse(Float64, rainfall_val) / 3_600_000.0 : 0.0

    # --output / --output-interval
    output_val,     args = _pop_flag(args, "--output")
    output_int_val, args = _pop_flag(args, "--output-interval")
    output_interval = output_int_val !== nothing ? parse(Float64, output_int_val) : 60.0

    # --injection-point lat,lon,rate  (repeatable)
    # Format: --injection-point -43.386,172.648,0.5  (lat, lon, m³/s)
    injection_specs = Tuple{Float64,Float64,Float64}[]
    while true
        inj_val, args = _pop_flag(args, "--injection-point")
        inj_val === nothing && break
        parts = split(inj_val, ",")
        length(parts) == 3 ||
            (println("ERROR: --injection-point must be lat,lon,rate_m3s\n"); print_help(1))
        lat_inj = parse(Float64, strip(parts[1]))
        lon_inj = parse(Float64, strip(parts[2]))
        rate    = parse(Float64, strip(parts[3]))
        push!(injection_specs, (lon_inj, lat_inj, rate))   # stored as (lon,lat,rate)
    end

    # --rainpoint lat,lon,mm_hr  (repeatable)
    # Format: --rainpoint -43.531,172.636,50.0  (lat, lon, mm/hr)
    rainpoint_specs = Tuple{Float64,Float64,Float64}[]
    while true
        rp_val, args = _pop_flag(args, "--rainpoint")
        rp_val === nothing && break
        parts = split(rp_val, ",")
        length(parts) == 3 ||
            (println("ERROR: --rainpoint must be lat,lon,mm_hr\n"); print_help(1))
        lat_rp  = parse(Float64, strip(parts[1]))
        lon_rp  = parse(Float64, strip(parts[2]))
        mm_hr   = parse(Float64, strip(parts[3]))
        push!(rainpoint_specs, (lon_rp, lat_rp, mm_hr))   # stored as (lon,lat,mm_hr)
    end

    # Validate mesh source
    has_gen  = meshgen_val !== nothing
    has_load = meshload_val !== nothing

    !has_gen && !has_load &&
        (println("ERROR: provide --meshgen or --meshload\n"); print_help(1))
    has_gen && has_load &&
        (println("ERROR: --meshgen and --meshload are mutually exclusive\n"); print_help(1))
    has_gen && meshres_val === nothing &&
        (println("ERROR: --meshgen requires --meshres\n"); print_help(1))
    !isempty(args) &&
        (println("ERROR: unexpected arguments: $(join(args, " "))\n"); print_help(1))
    
    

    mesh_source = has_gen ?
        (:generate, meshgen_val, parse(Int, meshres_val), meshout_val) :
        (:load, meshload_val, meshout_val)

    dem_source      = dem_val      !== nothing ? FileDEM(dem_val)      : nothing
    friction_source = friction_val !== nothing ? friction_val          : nothing
    
    run_flood_model(;
        mesh_source     = mesh_source,
        mesh_only       = mesh_only,
        vis_mode        = vis_mode,
        vis_port        = vis_port,
        dem_source      = dem_source,
        dem_strict      = dem_strict,
        dem_method      = dem_method,
        dem_samples     = dem_samples,
        halton_seed     = halton_seed,
        flow_method     = flow_method,
        sgs_bins        = sgs_bins,
        sgs_samples     = sgs_samples,
        manning_n       = manning_n,
        friction_source = friction_source,
        sim_duration    = sim_duration,
        dt_max          = dt_max,
        rainfall_rate    = rainfall_rate,
        injection_specs  = injection_specs,
        rainpoint_specs  = rainpoint_specs,
        output_path      = output_val,
        output_interval  = output_interval,
    )
end


@info "Starting FloodA5 model..." Dates.now()
let
    # invokelatest is required under Julia 1.12+ strict world-age semantics.
    # All top-level definitions (main, run_flood_model, etc.) are defined in
    # a prior world relative to this let-block; invokelatest ensures we always
    # call the most recent definition, which also suppresses the world-age warnings.

    # 0. Skip execution entirely when included by an external script (e.g.
    #    benchmark_sim.jl).  The caller sets this env var before include() and
    #    clears it immediately after; we just need to not call main() in that case.
    if get(ENV, "FLOODMODEL_INCLUDE_ONLY", "") == "1"
        # nothing — symbols are now defined; caller drives execution

    # 1. REPL / interactive session with no args: apply default test set
    elseif isinteractive() && isempty(ARGS)
        @info "REPL detected. Applying default test set."
        REPL_ARGS = [
            "--meshload", "test/kaiapoi_mesh16_sgs.parquet", 
            "--rainfall", "50", 
            "--output", "test/kaiapoi_test4_repl.h5",
            "--output-interval", "3600",
            "--sim-duration", "72000",
            "--flow-model", "sgs"
        ]
        Base.invokelatest(main, REPL_ARGS)

    # 2. Shell or VS Code debugger: run with the supplied ARGS
    else
        is_shell    = abspath(PROGRAM_FILE) == @__FILE__
        is_debugger = occursin("run_debugger.jl", PROGRAM_FILE)
        if is_shell || is_debugger
            Base.invokelatest(main, ARGS)
        end
    end
end
@info "FloodA5 model finished." Dates.now()

