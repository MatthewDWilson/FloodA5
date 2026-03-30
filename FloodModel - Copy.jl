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
include(joinpath(@__DIR__, "A5Grid.jl"))
include(joinpath(@__DIR__, "VisualisationServer.jl"))
include(joinpath(@__DIR__, "MakieVisualiser.jl"))

using .A5Grid
using .A5Grid: SGSTable, wse_from_volume, wetted_area_from_wse,
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
Hydrodynamic state of the flood model at a single timestep.

Primary state variable is `volume` (m³ stored water per cell).  All other
fields are either static mesh geometry or diagnostics updated each step.
Fields are in mesh.cells order (index i ↔ cell i).
"""
mutable struct FlowState
    cell_ids    :: Vector{String}
    water_depth :: Vector{Float64}   # m above local bed (diagnostic)
    volume      :: Vector{Float64}   # m³ stored — primary state variable
    velocity    :: Vector{Float64}   # m/s scalar magnitude (diagnostic)
    elevation   :: Vector{Float64}   # bed elevation above datum (m)
    manning_n   :: Vector{Float64}   # Manning's roughness per cell
    cell_area   :: Vector{Float64}   # plan area (m²)
    adjacency   :: Dict{String, Vector{String}}  # cell_id → neighbour ids
    adj_matrix  :: Matrix{Int}       # (max_nb × n_cells) index matrix, 0=none
    edge_width  :: Matrix{Float64}   # (max_nb × n_cells) shared edge length (m)
    edge_sill   :: Matrix{Float64}   # (max_nb × n_cells) sill elevation (m)
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

TODO: add chunking and compression (HDF5.jl chunk/deflate options) once the
solver produces real data and representative frame sizes are known.
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
    sat = saturation_fraction(state)
    HDF5.h5open(output.path, "r+") do fid
        g = HDF5.create_group(fid["frames"], frame_name)
        g["t"]           = t
        g["water_depth"] = state.water_depth
        g["volume"]      = state.volume
        g["saturation"]  = sat
        g["velocity"]    = state.velocity
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
    _build_adjacency_grid_disk(mesh) → Dict{String, Vector{String}}

Build exact topological adjacency using pya5.grid_disk(cell, 1).
Returns only cells that are present in the mesh (inter-mesh boundary cells
returned by grid_disk but not in the AOI are silently dropped).
"""
function _build_adjacency_grid_disk(mesh::A5Mesh)::Dict{String,Vector{String}}
    all_ids = Set(c.id for c in mesh.cells)
    adj     = Dict{String,Vector{String}}()
    for cell in mesh.cells
        nbrs = grid_disk_neighbours(cell.id)
        adj[cell.id] = filter(id -> id in all_ids, nbrs)
    end
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
    ids     = [c.id for c in mesh.cells]
    id_idx  = Dict{String,Int}(ids[i] => i for i in 1:n)

    # ── Elevation ──────────────────────────────────────────────────────────
    elevations = if haskey(mesh.static_vars, "elevation")
        copy(mesh.static_vars["elevation"])
    else
        @warn "No elevation data — all cells at z=0. Use --dem before running."
        zeros(Float64, n)
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
    areas = if haskey(mesh.static_vars, "sgs_cell_area")
        copy(mesh.static_vars["sgs_cell_area"])
    else
        [_polygon_area_m2(c.boundary) for c in mesh.cells]
    end

    # ── Exact adjacency via grid_disk ──────────────────────────────────────
    @info "Building exact topological adjacency via grid_disk..."
    t0  = time()
    adj = _build_adjacency_grid_disk(mesh)
    @info "  Adjacency built in $(round(time()-t0, digits=1))s"

    max_nb     = 5
    adj_matrix = zeros(Int, max_nb, n)
    edge_width = fill(NaN, max_nb, n)

    for i in 1:n
        cell_i = mesh.cells[i]
        nbrs   = get(adj, ids[i], String[])
        for (slot, nb_id) in enumerate(nbrs)
            slot > max_nb && break
            j = get(id_idx, nb_id, 0)
            j == 0 && continue
            adj_matrix[slot, i] = j
            w = _edge_length_m(cell_i.boundary, mesh.cells[j].boundary)
            edge_width[slot, i] = isnan(w) ? sqrt(areas[i]) : w
        end
    end

    # ── SGS tables ─────────────────────────────────────────────────────────
    sgs_tables = SGSTable[]
    edge_sill  = fill(NaN, max_nb, n)

    if method isa SGSFlow
        haskey(mesh.array_vars, "sgs_elev_bins") ||
            error("SGS tables not found in mesh — call build_sgs_tables! first.")
        @info "Loading SGS tables from mesh array_vars..."
        sgs_tables = [sgs_table(mesh, i) for i in 1:n]

        if haskey(mesh.array_vars, "sgs_edge_sills")
            edge_sill = copy(mesh.array_vars["sgs_edge_sills"])
        else
            @warn "sgs_edge_sills not found — using cell z_min fallback."
            z_mins = get(mesh.static_vars, "sgs_z_min", elevations)
            for i in 1:n, slot in 1:max_nb
                j = adj_matrix[slot, i]
                j == 0 && continue
                edge_sill[slot, i] = min(z_mins[i], z_mins[j])
            end
        end
    end

    volumes = zeros(Float64, n)

    return FlowState(
        ids,
        zeros(Float64, n),
        volumes,
        zeros(Float64, n),
        elevations,
        n_vec,
        areas,
        adj,
        adj_matrix,
        edge_width,
        edge_sill,
        sgs_tables,
    )
end

# Backwards-compatible no-method overload
initialise_flow_model(mesh::A5Mesh) = initialise_flow_model(mesh, StandardFlow())

# ---------------------------------------------------------------------------
# Diffusive wave physics
# ---------------------------------------------------------------------------

const _G = 9.81   # m s⁻²

"""
    _diffusive_flux(wse_i, wse_j, z_sill, A_wet_i, width, L, n_i) → Float64

Compute the diffusive-wave volumetric flux (m³/s) from cell i to cell j.

Uses Manning's equation for the cross-section:
    Q = (A_wet / n) × R^(2/3) × S^(1/2)

where:
  S   = |WSE_i - WSE_j| / L          (water surface slope)
  L   = distance between cell centres (m), passed in as `width` here
       [Note: caller passes the shared edge width as L per user decision] 
       #### TODO: This should be computed between cell centres on the sphere using haversine (great circle) distances
       #### TODO: However, the grid is non-orthogonal and a slope correction will be needed because of this. 

  R   = hydraulic radius ≈ A_wet / width  (wide-channel approximation)
  A_wet = wetted cross-sectional area ≈ wet_depth × edge_width

For SGS, wet_depth at the edge = WSE_i - z_sill (overtopping depth).
For standard flow, wet_depth = water_depth at cell i.

Returns a positive value if flow is from i → j (wse_i > wse_j).
Returns 0.0 if no overtopping (wse_i ≤ z_sill).
"""
@inline function _diffusive_flux(wse_i::Float64, wse_j::Float64,
                                  z_sill::Float64,
                                  A_wet_i::Float64,
                                  width::Float64,
                                  n_mann::Float64)::Float64
    wse_i <= z_sill && return 0.0
    dh = wse_i - max(wse_j, z_sill)
    dh <= 0.0 && return 0.0

    # Overtopping depth at the sill
    h_sill = wse_i - z_sill
    # Cross-sectional area at edge: rectangular section, width = edge width
    A_cross = h_sill * width
    # Hydraulic radius (wide channel): R ≈ h_sill
    R = h_sill
    # Slope
    S = dh / max(width, 1.0)
    # Manning's Q (m³/s)
    return (A_cross / n_mann) * R^(2.0/3.0) * sqrt(S)
    #### TODO: change this formula to use a modified version of Bates et al. 2010 inertial formulation
    #### Here is a tex of this equation: Q_x^t = \frac{q_x^{t-\Delta t} - g h^t_{[\text{flow}]} \Delta t \frac{\Delta (h^t + z)}{\Delta x} } {(1 + g h^t_{[\text{flow}]} \Delta t n^2 |q_x^{t-\Delta t}| / (h^t_{[\text{flow}]} )^{10/3})} \Delta y 
 end

"""
    _cfl_dt(state, method) → Float64

Compute the maximum stable timestep from the CFL condition for the
diffusive wave equation.  The diffusive stability criterion is:

    dt ≤ CFL × dx² / (2D)

where D = Q/A is the diffusivity, dx is the cell length scale,
and CFL = 0.5 is a safety factor.

We approximate D conservatively as g × h_max × dx / n_min for each
active cell, then take the minimum over all cells.
"""
function _cfl_dt(state::FlowState, method::FlowMethod; cfl::Float64=0.5)::Float64
    n        = length(state.cell_ids)
    dt_min   = Inf
    h_thresh = 1e-4   # ignore nearly-dry cells

    for i in 1:n
        h = state.water_depth[i]
        h < h_thresh && continue
        # Length scale: sqrt of cell area
        dx  = sqrt(state.cell_area[i])
        n_i = state.manning_n[i]
        # Diffusivity: D ≈ (1/n) × h^(5/3) / S^(1/2)
        # Bound S away from zero for stability; use a representative 0.001 if flat
        D   = (1.0 / n_i) * h^(5.0/3.0) * 0.032   # √0.001 ≈ 0.032 typical slope
        D   = max(D, 0.001)
        dt  = cfl * dx^2 / (2.0 * D)
        dt < dt_min && (dt_min = dt)
    end

    return isfinite(dt_min) ? dt_min : 60.0   # default 60s if all cells dry
end

# ---------------------------------------------------------------------------
# Standard flow step (mean elevation, no SGS)
# ---------------------------------------------------------------------------

"""
    step_standard!(state, dt)

Advance the standard (non-SGS) diffusive-wave model by one timestep `dt`.

WSE = elevation + water_depth.  Flow uses the mean cell bed elevation as
the sill between adjacent cells.
"""
function step_standard!(state::FlowState, dt::Float64)
    n       = length(state.cell_ids)
    dV      = zeros(Float64, n)   # volume increments this step

    for i in 1:n
        state.water_depth[i] < 1e-6 && continue
        wse_i = state.elevation[i] + state.water_depth[i]
        A_i   = state.cell_area[i]

        for slot in 1:5
            j = state.adj_matrix[slot, i]
            j == 0 && break
            wse_j  = state.elevation[j] + state.water_depth[j]
            wse_i <= wse_j && continue   # only route downhill

            # Sill = max of the two bed elevations (weir-like)
            z_sill = max(state.elevation[i], state.elevation[j])
            width  = state.edge_width[slot, i]
            isnan(width) && continue

            Q = _diffusive_flux(wse_i, wse_j, z_sill, A_i,
                                 width, state.manning_n[i])
            vol = Q * dt
            # Don't remove more than is available
            vol = min(vol, state.water_depth[i] * A_i * 0.5)
            dV[i] -= vol
            dV[j] += vol
        end
    end

    # Update depths from volumes
    for i in 1:n
        A_i = state.cell_area[i]
        A_i < 1.0 && continue
        new_depth = state.water_depth[i] + dV[i] / A_i
        state.water_depth[i] = max(0.0, new_depth)
    end

    # Scalar velocity magnitude (diagnostic only)
    for i in 1:n
        state.volume[i] = state.water_depth[i] * state.cell_area[i]
    end
end

# ---------------------------------------------------------------------------
# SGS flow step (hypsometric volume lookup)
# ---------------------------------------------------------------------------

"""
    step_sgs!(state, dt)

Advance the SGS diffusive-wave model by one timestep `dt`.

WSE for each cell is derived from its stored volume via the hypsometric
curve (inverse lookup of vol_curve).  Wetted area is obtained from the
area_curve.  The edge sill is the minimum DEM elevation along the shared
boundary, sampled during pre-processing.

This allows water to be "partially wet" — a cell with volume below the
sill still exerts a real WSE, and flow only begins once WSE exceeds the
sill.  Channels and ditches sub-grid are captured correctly.
"""
function step_sgs!(state::FlowState, dt::Float64)
    n     = length(state.cell_ids)
    dV    = zeros(Float64, n)

    # Step 1: compute WSE and wetted area for every cell
    wse      = Vector{Float64}(undef, n)
    A_wet    = Vector{Float64}(undef, n)
    for i in 1:n
        tbl      = state.sgs_tables[i]
        wse[i]   = wse_from_volume(tbl, state.volume[i])
        A_wet[i] = wetted_area_from_wse(tbl, wse[i])
    end

    # Step 2: route flow across edges
    for i in 1:n
        state.volume[i] < 1e-6 && continue

        for slot in 1:5
            j = state.adj_matrix[slot, i]
            j == 0 && break
            wse[i] <= wse[j] && continue   # only route downhill

            z_sill = state.edge_sill[slot, i]
            isnan(z_sill) && (z_sill = min(state.sgs_tables[i].z_min,
                                            state.sgs_tables[j].z_min))
            width  = state.edge_width[slot, i]
            isnan(width) && continue

            Q   = _diffusive_flux(wse[i], wse[j], z_sill, A_wet[i],
                                   width, state.manning_n[i])
            vol = Q * dt
            # Volume limiter: don't route more than half the available in one step
            vol = min(vol, state.volume[i] * 0.5)
            dV[i] -= vol
            dV[j] += vol
        end
    end

    # Step 3: update volumes and diagnostic fields
    for i in 1:n
        state.volume[i] = max(0.0, state.volume[i] + dV[i])
        tbl = state.sgs_tables[i]
        # Skip degenerate cells (no valid DEM samples → NaN z_min)
        if isnan(tbl.z_min)
            state.water_depth[i] = 0.0
            continue
        end
        new_wse = wse_from_volume(tbl, state.volume[i])
        state.water_depth[i] = max(0.0, new_wse - tbl.z_min)
    end
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
    run_simulation!(state, mesh, sim_duration, dt_max, vis, vis_mode, output,
                    method, rainfall_rate)

Time-stepping loop.  Advances the model for `sim_duration` seconds using
adaptive dt (CFL-limited, capped at `dt_max`).

Dispatches `push_frame!` to the active visualiser every `VIS_INTERVAL` steps.
Writes HDF5 snapshots according to `output.output_interval`.

Arguments
---------
  method         — `StandardFlow()` or `SGSFlow()`
  rainfall_rate  — uniform rainfall in m/s (default 0; e.g. 1e-5 for ~36 mm/hr)
"""
function run_simulation!(state         :: FlowState,
                         mesh          :: A5Mesh,
                         sim_duration  :: Float64,
                         dt_max        :: Float64,
                         vis,
                         vis_mode      :: Symbol = :none,
                         output        :: SimOutput = SimOutput();
                         method        :: FlowMethod = StandardFlow(),
                         rainfall_rate :: Float64   = 0.0)
    t    = 0.0
    step = 0
    use_sgs = method isa SGSFlow

    while t < sim_duration
        # Adaptive dt
        dt = min(_cfl_dt(state, method), dt_max, sim_duration - t)
        dt = max(dt, 0.1)   # floor: 0.1s to prevent infinite loops on dry mesh

        # Apply rainfall source (uniform, before routing)
        if rainfall_rate > 0.0
            for i in eachindex(state.cell_ids)
                if use_sgs
                    state.volume[i] += rainfall_rate * dt * state.cell_area[i]
                else
                    state.water_depth[i] += rainfall_rate * dt
                end
            end
        end

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
                VisualisationServer.push_frame!(
                    vis, state.cell_ids, state.water_depth, t)
            elseif vis_mode === :makie
                MakieVisualiser.push_frame!(
                    vis, state.cell_ids, state.water_depth, t)
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
            @info "  step=$(lpad(step,5))  t=$(round(t,digits=1))s  " *
                  "dt=$(round(dt,digits=2))s  " *
                  "wet=$n_wet  " *
                  "max_depth=$(round(max_depth,digits=3))m"
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
    sim_duration     :: Float64 = 3600.0,
    dt_max           :: Float64 = 60.0,
    rainfall_rate    :: Float64 = 0.0,
    output_path      :: Union{String,Nothing} = nothing,
    output_interval  :: Float64 = 60.0)

    @info "=== A5 Flood Model ===" Dates.now()
    @info "Vis mode    : $vis_mode"
    @info "Flow method : $flow_method"

    # 1. Start Cesium server early
    vis = if vis_mode === :cesium
        VisualisationServer.start(port = vis_port,
                                  viz_dir = joinpath(@__DIR__, "viz"))
    else
        nothing
    end
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

    # 4. SGS pre-processing (if using SGS flow)
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

    # 5. Log mesh summary
    @info mesh_summary(mesh)

    # 6. Hand mesh to visualiser
    if vis_mode === :makie
        src_label = basename(mesh_source[2])
        @info "Opening Makie viewer ($(length(mesh)) cells)..."
        vis = MakieVisualiser.start(mesh;
                  title = "FloodA5  res=$(mesh.resolution)  $src_label")
    elseif vis_mode === :cesium
        @info "Pushing mesh to Cesium server ($(length(mesh)) cells)..."
        VisualisationServer.set_mesh!(vis, mesh_to_geojson_string(mesh))
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
                    method        = method_obj,
                    rainfall_rate = rainfall_rate)

    # 10. Keep visualiser alive
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
      --dem linz_dem.tif --flow-model sgs --rainfall 30 --sim-duration 3600 --vis

  # Reload mesh with pre-built SGS tables, run with friction raster
  julia --threads auto FloodModel.jl \\
      --meshload mesh_sgs.parquet --friction land_use_n.tif \\
      --rainfall 10 --sim-duration 7200 --output sim_out.h5

  # Standard (non-SGS) run for quick testing
  julia --threads auto FloodModel.jl \\
      --meshload mesh_sgs.parquet --flow-model standard --sim-duration 1800

Resolution guide (approximate cell area):
  Level  5  ~5 000 km²  Continental / regional
  Level  8  ~250 km²    Large catchment
  Level 10  ~50 km²     Medium catchment
  Level 12  ~10 km²     Small catchment
  Level 14  ~2 km²      Urban / detailed
  Level 17  ~0.1 km²    High-resolution modelling
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

function main()
    args = copy(ARGS)

    "--help" in args || "-h" in args && print_help(0)

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

    # --meshgen / --meshres / --meshout / --meshload
    meshgen_val,  args = _pop_flag(args, "--meshgen")
    meshres_val,  args = _pop_flag(args, "--meshres")
    meshout_val,  args = _pop_flag(args, "--meshout")
    meshload_val, args = _pop_flag(args, "--meshload")

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

    flow_method  = flow_method_val  !== nothing ? Symbol(flow_method_val)       : :sgs
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
        rainfall_rate   = rainfall_rate,
        output_path     = output_val,
        output_interval = output_interval,
    )
end

main()
