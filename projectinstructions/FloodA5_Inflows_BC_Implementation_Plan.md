# FloodA5 — Dynamic Inflows & Open Boundary Conditions: Implementation Plan

**Branch:** `inflows_and_bcs`  
**Date:** 2026-06-10  
**Status:** Design confirmed — ready for implementation

---

## 1. Scope

This document covers two features to be implemented together in branch `inflows_and_bcs`:

**Feature 1 — Dynamic fluvial inflows:** Time-varying volumetric inflow (m³/s) at one or more upstream boundary cells, driven by a hydrograph in a text file. Initial format: LISFLOOD-FP `.bci`/`.bdy` style for cross-portability, with a clean extensibility path to HydroML and CF-NetCDF in future sessions.

**Feature 2 — Open outflow boundaries:** Allow water to exit the domain at domain-edge cells by default, using a zero-gradient (kinematic slope) approximation. User can override per-boundary segment via a GeoJSON vector file to specify closed or custom BC types per segment.

Both features are designed with the full future development plan in mind (tide BCs, spatial rainfall, infiltration, soil/friction inputs, Phase 3 AMR). No shortcuts that would require a structural rewrite.

---

## 2. Design Principles

- **No FlowState surgery:** All new logic slots into existing patterns — source application before the flux step, BC application after the flux step.
- **Extensibility-first:** `AbstractSource` and `AbstractBC` hierarchies allow new types (tide, infiltration, spatial rainfall) to be added as new concrete subtypes without touching the simulation loop.
- **Boundary cell detection is mesh-level, not simulation-level:** Which cells are domain edges is a property of the mesh, computed once at `initialise_flow_model` time.
- **Input format abstraction:** A `TimeSeriesReader` abstraction isolates file I/O from the rest of the model. LISFLOOD-FP format is the first implementation; adding HydroML or CSV later requires only a new reader subtype.
- **Mass budget accounting:** Inflows add to `input_vol`; outflow losses accumulate in `vol_removed`. The existing `vol_removed` placeholder in the Makie visualiser is already wired for this.
- **No global state:** All new structs are passed explicitly into `run_simulation!` and `run_flood_model`. The CLI parses new flags and builds the appropriate structs.

---

## 3. Feature 1 — Dynamic Fluvial Inflows

### 3.1 Concept

A **fluvial inflow** is a point source whose flow rate varies in time according to a user-supplied hydrograph. It is injected into the nearest mesh cell, exactly like the existing `InjectionPoint`, but with rate `f(t)` rather than a constant. The cell receiving the inflow is typically at or near the upstream boundary of the domain — a river channel cell where the river enters the flood model extent.

Physically, this represents river discharge passing an upstream gauge point. The model does not route water through a channel network internally; it receives a pre-computed hydrograph (from a rainfall-runoff model or gauge observation) and injects the corresponding volume flux.

Multiple inflow points are fully composable with each other and with existing `--injection-point`, `--rainpoint`, and `--rainfall` inputs. If a `.bci`/`--inflow-point` entry and a `--injection-point` entry both target the same cell, their volumes are summed before flux routing — this is by design and should be documented.

### 3.2 Data structures

```julia
# boundaryinputs/sources.jl

"""
    AbstractSource

Supertype for all water-adding boundary inputs.  Concrete subtypes must implement:
  apply_source!(state, source, t, dt)  — add dV to state.volume
  cumulative_volume(source, t)         — total volume injected up to time t (m³)

Existing types InjectionPoint and RainPoint are promoted to AbstractSource subtypes.
New types (SpatialRainfall, etc.) are added here in future sessions.
"""
abstract type AbstractSource end

"""
    InflowPoint <: AbstractSource

Time-varying volumetric inflow at a single mesh cell.
Rate is linearly interpolated from a hydrograph (t_s, Q_m3s) at each timestep.
Flat extrapolation beyond the first and last knot points (not clipped to zero —
a hydrograph starting mid-storm still carries its first-knot value before t=0).
"""
struct InflowPoint <: AbstractSource
    cell_index   :: Int
    cell_id      :: String
    lon          :: Float64            # source longitude (degrees)
    lat          :: Float64            # source latitude (degrees)
    t_s          :: Vector{Float64}    # times (seconds from sim start), monotonically increasing
    Q_m3s        :: Vector{Float64}    # discharge (m³/s), same length as t_s
    label        :: String             # gauge ID / name for logging and mass balance output
end
```

`InjectionPoint` and `RainPoint` gain `<: AbstractSource` as their only structural change.

### 3.3 Hydrograph interpolation

```julia
"""
    _interp_hydrograph(t_s, Q_m3s, t) → Float64

Linear interpolation of a (time, discharge) hydrograph at simulation time t.
Flat extrapolation beyond first/last knot.  Binary search for O(log n) bracket.
"""
@inline function _interp_hydrograph(t_s   :: Vector{Float64},
                                     Q_m3s :: Vector{Float64},
                                     t     :: Float64)::Float64
    isempty(t_s)      && return 0.0
    t <= t_s[1]       && return Q_m3s[1]
    t >= t_s[end]     && return Q_m3s[end]
    lo = searchsortedlast(t_s, t)
    frac = (t - t_s[lo]) / (t_s[lo+1] - t_s[lo])
    return Q_m3s[lo] + frac * (Q_m3s[lo+1] - Q_m3s[lo])
end
```

### 3.4 Cumulative volume (for mass balance)

```julia
"""
    cumulative_volume(src::InflowPoint, t) → Float64

Cumulative volume injected from t=0 to t, computed by trapezoidal integration
over the stored hydrograph knots, then interpolating the final partial interval.
Used by the mass balance logger.
"""
function cumulative_volume(src::InflowPoint, t::Float64)::Float64
    isempty(src.t_s) && return 0.0
    t <= src.t_s[1]  && return src.Q_m3s[1] * t
    vol = 0.0
    for k in 1:length(src.t_s)-1
        t1, t2 = src.t_s[k], src.t_s[k+1]
        Q1, Q2 = src.Q_m3s[k], src.Q_m3s[k+1]
        if t2 <= t
            vol += 0.5 * (Q1 + Q2) * (t2 - t1)
        else
            frac = (t - t1) / (t2 - t1)
            Q_t  = Q1 + frac * (Q2 - Q1)
            vol += 0.5 * (Q1 + Q_t) * (t - t1)
            break
        end
    end
    t >= src.t_s[end] && (vol += src.Q_m3s[end] * (t - src.t_s[end]))
    return vol
end

# Constant-rate sources — trivial implementations
cumulative_volume(src::InjectionPoint, t::Float64) = src.rate_m3s * t
cumulative_volume(src::RainPoint,      t::Float64) = src.rate_m3s * t
```

### 3.5 LISFLOOD-FP file format support

#### 3.5.1 `.bdy` file format

The `.bdy` file contains one or more named time series. Multiple named series can appear in one file. Time units are `seconds`, `hours`, or `days` (converted to seconds on read). FloodA5 reads this format without modification — a `.bdy` produced for LISFLOOD-FP runs directly in FloodA5.

```
GAUGE_UPSTREAM
12                    # number of time steps
hours                 # time unit: seconds | hours | days
0.0    0.0
1.0    5.2
2.0    18.1
3.0    42.7
...
GAUGE_DOWNSTREAM
8
hours
0.0    0.0
...
```

Comments (lines beginning `#`) and blank lines are ignored. The reader returns `Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}` mapping series name to `(t_seconds, Q_m3s)`.

#### 3.5.2 `.bci` file format

The `.bci` file maps boundary locations to BC types. FloodA5 implements the LISFLOOD-FP format with the following columns:

| Col | Description |
|-----|-------------|
| 1 | Boundary type: `N`, `E`, `S`, `W` (edge), `P` (point source), `F` (internal free boundary) |
| 2 | Start of segment (easting/lon) — for `P`: longitude or easting of the point |
| 3 | End of segment (northing/lat) — for `P`: latitude or northing of the point |
| 4 | BC type: `QVAR`, `QFIX`, `HFIX`, `FREE` (see Table 1 below) |
| 5 | BC value: series name (for `QVAR`), constant discharge m³/s (for `QFIX`), WSE (for `HFIX`) |

**Table 1 — BC type codes**

| Code | Meaning | FloodA5 status |
|------|---------|----------------|
| `QVAR` | Time-varying discharge from `.bdy` series | ✅ Implemented (InflowPoint) |
| `QFIX` | Fixed discharge (m³/s) | ✅ Implemented (maps to InjectionPoint) |
| `FREE` | Zero-gradient outflow | ✅ Implemented (ZeroGradient BC) |
| `HFIX` | Fixed water surface elevation | ⏳ Future (tide session) |

**Coordinate system:** FloodA5 defaults to **decimal degrees (WGS 84)** matching LISFLOOD-FP's `latlong` option. If the user supplies a projected CRS EPSG code (via `--bc-epsg` CLI flag, see §3.7), coordinates are treated as metres (easting, northing) and converted to lon/lat internally using `Proj.jl` or the `ArchGDAL` coordinate transform already available in the stack. The default (no `--bc-epsg`) is always geographic, making it consistent with FloodA5's coordinate system throughout.

**`P` (point source) entries:** columns 2 and 3 are lon and lat (or easting and northing). The nearest mesh cell is found with `_find_nearest_cell`, with a warning logged if the distance exceeds 2 km. There must not be more than one `.bci` `P` entry per cell — this is checked at parse time and raises an error if violated.

**`N/E/S/W` (edge) entries:** columns 2 and 3 are the start and end coordinates of the boundary segment. All boundary cells whose centres project onto this segment (within spatial tolerance) are assigned the specified BC type. This is used primarily for `FREE` (open outflow) and future `HFIX` (tide) segments. For `QVAR`/`QFIX` on edge entries, the discharge is distributed uniformly across the matched cells — this is the LISFLOOD-FP convention for distributed inflow along a river boundary.

**`F` (internal free boundary):** not implemented in this session; parsed without error and logged as unsupported.

#### 3.5.3 Two-column CSV fallback

When `--inflow-point` is used without a `.bci` file, the hydrograph file can be:
- A plain two-column CSV (no header) with columns `t_s, Q_m3s`
- A two-column CSV with a header row (auto-detected: if first row non-numeric, skip it)
- A LISFLOOD-FP `.bdy` file with a single series (series name is ignored; the one series is used)

The reader auto-detects format from the file content, not the extension.

### 3.6 Reader abstraction

```julia
# boundaryinputs/timeseries_io.jl

abstract type AbstractTimeSeriesReader end

struct LisfloodBDYReader <: AbstractTimeSeriesReader
    path :: String
end

struct TwoColumnCSVReader <: AbstractTimeSeriesReader
    path  :: String
    t_col :: Int   # default 1
    Q_col :: Int   # default 2
end

# Future readers (stubs):
# struct WaterML2Reader <: AbstractTimeSeriesReader; path::String; end
# struct CFNetCDFGaugeReader <: AbstractTimeSeriesReader; path::String; station_id::String; end

"""
    read_timeseries(reader) → Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}

Returns a Dict mapping series name → (t_seconds, Q_m3s).
For TwoColumnCSVReader and single-series .bdy files, the dict has one entry
with a synthetic key equal to the filename (without extension).
"""
function read_timeseries(r::LisfloodBDYReader) :: Dict{String, Tuple{Vector{Float64},Vector{Float64}}}
    # parse .bdy; handle time units seconds/hours/days; skip comments and blanks
    # ...
end

function read_timeseries(r::TwoColumnCSVReader) :: Dict{String, Tuple{Vector{Float64},Vector{Float64}}}
    # auto-detect header; parse two columns
    # ...
end
```

### 3.7 CLI flags

```
--inflow-point LAT,LON,FILE[,LABEL]
    Add a time-varying inflow at (lat, lon) driven by FILE.
    FILE is a two-column CSV (t_s, Q_m3s), a single-series .bdy file,
    or a .bdy file — in which case LABEL must match a series name in the file.
    LABEL is optional; defaults to the filename stem. Used in log output.
    Repeatable for multiple inflow points.
    Example: --inflow-point -43.386,172.648,waimakariri.csv,Waimakariri

--inflow-bci FILE
    Load inflow and internal boundary configuration from a LISFLOOD-FP .bci file.
    QVAR entries reference series in a .bdy file in the same directory as FILE.
    QFIX entries create constant-rate InjectionPoints.
    FREE entries override the default open-boundary type to ZeroGradient.
    Repeatable; entries accumulate.
    Example: --inflow-bci carlisle_inflows.bci

--bc-epsg CODE
    Treat coordinates in --inflow-bci files as projected (easting/northing)
    in the given EPSG coordinate reference system, converting to lon/lat internally.
    Default (omitted): coordinates are decimal degrees WGS 84.
    Example: --bc-epsg 27700   (British National Grid)
```

### 3.8 Time reference

All hydrograph time values are **seconds from simulation start** (relative), consistent with `--sim-duration`. This matches LISFLOOD-FP's convention.

For future absolute-time support (e.g., HydroML ISO 8601 timestamps), the `_interp_hydrograph` function accepts a `t_offset::Float64 = 0.0` keyword that is subtracted from the hydrograph times before interpolation. The CLI and JSON config paths will expose this as `--inflow-t-start "2025-01-15T06:00:00"` in a later session.

---

## 4. Feature 2 — Open Outflow Boundaries

### 4.1 Concept

**Boundary cells** are mesh cells with fewer than 5 edge-sharing neighbours — their missing edges face the exterior of the domain. Currently these are implicitly treated as closed walls (no flux). The default is changed to **open (ZeroGradient)**: water that reaches a domain edge flows out freely.

The zero-gradient (kinematic) approximation sets the ghost-cell WSE equal to the boundary cell's WSE, so `dWSE = 0` across the ghost edge. This means the momentum term carries existing flow off the edge without reflection, but there is no additional pressure gradient driving flow out. This is the standard "non-reflective" or "transmissive" outflow condition used in most flood models.

### 4.2 Boundary cell detection

At `initialise_flow_model` time, after the EdgeList is built, identify domain-edge cells:

```julia
n_neighbours = zeros(Int, n)
for e in 1:edges.n_edges
    n_neighbours[edges.cell_i[e]] += 1
    n_neighbours[edges.cell_j[e]] += 1
end
boundary_mask = BitVector(n_neighbours[i] < N_SIDES for i in 1:n)
```

This is O(n_edges) and negligibly fast. The result is stored in `FlowState.boundary_mask`.

### 4.3 BC type enum and BoundarySegment

```julia
# boundaryinputs/boundary_conditions.jl

"""
    BCType

Supported boundary condition types.  Closed is the legacy implicit behaviour.
ZeroGradient is the new default.  Critical and FixedWSE are planned for
subsequent sessions (free outfall and tide respectively).
"""
@enum BCType begin
    Closed         # no flux — legacy implicit behaviour
    ZeroGradient   # ghost WSE = boundary cell WSE; transmissive outflow
    Critical       # Q = width × √(g × h³); free outfall condition
    FixedWSE       # fixed water surface — tide/downstream stage (future)
    FixedQ         # fixed outflow discharge (future)
end

"""
    BoundarySegment

Associates a BC type with a set of boundary cell indices.
Constructed at init time from the --bc-file GeoJSON or from .bci FREE entries.
"""
struct BoundarySegment
    cell_indices :: Vector{Int}
    bc_type      :: BCType
    label        :: String
    # For FixedWSE (tide): time series loaded here, queried at runtime
    wse_series   :: Union{Tuple{Vector{Float64},Vector{Float64}}, Nothing}
end
```

### 4.4 Default behaviour

**Default (no `--bc-file` provided):** all boundary cells get `ZeroGradient`. This is a breaking change from the current implicit closed-boundary behaviour. The `--closed-boundaries` flag restores the old behaviour for mass balance benchmarking and catchment studies where no water should leave the domain.

The synthetic DEM test T4 mass balance criterion is updated from:
```
|input_vol - domain_vol| / input_vol < 0.01%
```
to:
```
|input_vol - (domain_vol + vol_removed)| / input_vol < 0.01%
```

### 4.5 GeoJSON boundary control file

Users can override the default per segment:

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": { "type": "LineString",
                    "coordinates": [[172.55, -43.60], [172.60, -43.60]] },
      "properties": { "bc_type": "Closed", "label": "Northern levee" }
    },
    {
      "type": "Feature",
      "geometry": { "type": "Polygon", ... },
      "properties": { "bc_type": "ZeroGradient", "label": "Southern outlet" }
    }
  ]
}
```

**Cell matching:** a boundary cell is assigned to a feature if its centre lies within `1.5 × cell_diameter` of the feature's geometry (haversine distance to nearest point on the linestring or polygon perimeter). Cell diameter is `2 × sqrt(cell_area / π)` — the equivalent circle radius. Cells not matched by any feature use the default BC type.

**CLI flag:**
```
--bc-file FILE    GeoJSON file specifying boundary condition types per segment.
                  Default if omitted: ZeroGradient for all boundary cells.
--closed-boundaries
                  Treat all boundary cells as closed walls.
                  Equivalent to --bc-file with all segments set to Closed.
                  Useful for mass balance benchmarking.
```

### 4.6 Ghost edge geometry

Each boundary cell has `N_SIDES - n_neighbours[ci]` missing neighbour slots. For each missing slot, a **ghost edge** is pre-computed at `initialise_flow_model` time.

**Ghost edge geometry fields:**

| Field | Value |
|-------|-------|
| `width` | **Actual estimated length** of the missing edge, computed from the boundary cell's polygon. The cell polygon has `N_SIDES` vertices; the edge length for each side is computed from the two adjacent vertices using `_haversine_m`. Missing slots (edges facing the exterior) are identified by their vertex pair and their length measured directly. This avoids the inaccuracy of the "mean edge width" approximation noted in review. |
| `L` | Distance from the boundary cell centre to a virtual point one cell-diameter outward along the outward-facing normal of the ghost edge. Computed from the midpoint of the ghost edge arc and the cell centre. |
| `cos_theta` | 1.0 — ghost edges are defined as orthogonal to the outward normal by construction. |
| `sill` (standard) | `elevation[ci]` — the boundary cell's own bed elevation. |
| `sill` (SGS) | `tbl.z_min` from the SGS table, consistent with how interior SGS sills are handled. |
| `flux_prev` | `0.0` at init time; updated each step via Fix C write-back. |

**Ghost edge identification from cell polygon:** the boundary cell polygon has `N_SIDES` sides. Sides that are *not* shared with any mesh neighbour are the ghost edges. These are identified by checking whether each edge's midpoint vertex pair appears in the `_shared_edge` lookup for any of the cell's neighbours. This is done once at init time.

```julia
struct GhostEdge
    cell_index :: Int        # boundary cell index
    width      :: Float64    # actual edge length (m), from polygon vertices
    L          :: Float64    # cell-centre to ghost-cell-centre distance (m)
    cos_theta  :: Float64    # always 1.0
    sill       :: Float64    # bed elevation at this edge (m)
    flux_prev  :: Float64    # q or Q from previous step (momentum state)
    is_Q_flux  :: Bool       # true = SGS R-A kernel (Q in m³/s); false = q (m²/s)
end
```

`GhostEdge` instances are stored as `state.ghost_edges :: Vector{GhostEdge}` (flat list, all boundary cells concatenated). A companion `state.ghost_cell_bc :: Vector{BCType}` stores the BC type per ghost edge (set at init from `BoundarySegment` matching).

### 4.7 Ghost-cell WSE by BC type

At each step, the ghost-cell WSE is set before calling the flux kernel:

```julia
function _ghost_wse(wse_ci::Float64, sill::Float64, bc::BCType,
                    wse_series_val::Float64=0.0)::Float64
    bc === ZeroGradient && return wse_ci
    bc === Critical     && return sill + (2.0/3.0) * max(0.0, wse_ci - sill)
    bc === FixedWSE     && return wse_series_val   # future
    bc === Closed       && return -Inf             # sentinel: skip this edge
    return wse_ci   # safe fallback
end
```

For `ZeroGradient`, `wse_ghost = wse_ci` → `dWSE = 0` → the numerator of Bates eq. 9 reduces to `q_prev` (pure momentum carry-out) → water exits at the rate it arrived, with no reflection.

### 4.8 Ghost-edge flux function (`flow2d.jl`)

```julia
"""
    _bates_ghost_flux(q_prev, wse_ci, wse_ghost, z_sill, width, L,
                      n_mann, dt, depth_ci) → (Q_out, q_stored)

Outflow flux across a ghost-cell boundary edge.

wse_ghost is computed by _ghost_wse() based on the BC type.
Q_out is always ≥ 0 (ghost cells never inject water — inflow from ghost
cells is physically meaningless and numerically destabilising).

Sign convention: Q_out > 0 means water leaves the domain through cell ci.
Caller applies: dV[ci] -= Q_out * dt

Uses _bates_flux_limited internally with Froude and volume limiters.
"""
@inline function _bates_ghost_flux(q_prev    :: Float64,
                                    wse_ci    :: Float64,
                                    wse_ghost :: Float64,
                                    z_sill    :: Float64,
                                    width     :: Float64,
                                    L         :: Float64,
                                    n_mann    :: Float64,
                                    dt        :: Float64,
                                    depth_ci  :: Float64)::Tuple{Float64,Float64}
    wse_ghost == -Inf && return (0.0, 0.0)   # Closed BC sentinel
    Q, q_s = _bates_flux_limited(q_prev, wse_ghost, wse_ci, z_sill,
                                   width, L, 1.0, n_mann, dt, depth_ci)
    # One-way enforcement: ghost cells never push water in
    Q < 0.0 && return (0.0, 0.0)
    return (Q, q_s)
end
```

Note the argument order `(wse_ghost, wse_ci)` — ghost is "cell_i" in the EdgeList sense. Since `wse_ghost ≤ wse_ci` by construction for ZeroGradient, the flux will be ≥ 0 (flow from `wse_ci` outward). The one-way guard catches edge cases where floating-point rounding produces a tiny negative flux.

For the SGS R-A kernel, a parallel `_manning_ghost_flux` is added using the same ghost-WSE logic but calling `_manning_flux_ra` with `A` and `R` computed from the ghost edge's SGS area/perimeter curves (using `wse_ci` as the reference WSE, since the ghost cell has no terrain data).

### 4.9 Phase D — boundary outflow pass in step functions

Both `step_standard!` and `step_sgs!` gain a **Phase D** after the existing Phase C (dV scatter):

```julia
# Phase D — boundary outflow (serial; ghost edges are per-cell, not shared)
vol_out_this_step = 0.0
for ge in state.ghost_edges
    state.ghost_cell_bc[ge] === Closed && continue
    ci = ge.cell_index
    wse_ci = use_sgs ? wse_from_volume(state.sgs_tables[ci], state.volume[ci]) :
                       state.elevation[ci] + state.volume[ci] / state.cell_area[ci]
    wse_ghost = _ghost_wse(wse_ci, ge.sill, state.ghost_cell_bc[ge])
    depth_ci  = state.water_depth[ci]
    Q_out, q_s = _bates_ghost_flux(ge.flux_prev, wse_ci, wse_ghost,
                                    ge.sill, ge.width, ge.L,
                                    0.5 * (state.manning_n[ci] + state.manning_n[ci]),  # ghost has same n
                                    dt, depth_ci)
    dV_ghost = Q_out * dt
    # Donor cap: ghost edge counts toward the boundary cell's DONOR_EDGE_DIVISOR budget
    max_out = state.volume[ci] / DONOR_EDGE_DIVISOR
    dV_ghost = min(dV_ghost, max_out)
    state.volume[ci]  -= dV_ghost
    state.water_depth[ci] = max(0.0, state.water_depth[ci] - dV_ghost / state.cell_area[ci])
    # Fix C: update ghost edge momentum state
    ge_new = GhostEdge(ge.cell_index, ge.width, ge.L, ge.cos_theta,
                        ge.sill, q_s, ge.is_Q_flux)
    state.ghost_edges[g_idx] = ge_new   # update in-place
    vol_out_this_step += dV_ghost
end
# Accumulate in state for mass balance reporting
state.vol_removed += vol_out_this_step
```

`vol_removed` is added to `FlowState` as a new running accumulator (Float64, initialised to 0.0).

Phase D is serial because ghost edges are owned by single cells (no write hazard), and the number of boundary cells is small relative to interior cells.

### 4.10 Mass balance accounting

`FlowState` gains:
```julia
mutable struct FlowState
    # ... existing fields ...
    boundary_mask  :: BitVector
    ghost_edges    :: Vector{GhostEdge}
    ghost_cell_bc  :: Vector{BCType}
    vol_removed    :: Float64           # cumulative outflow through ghost edges (m³)
end
```

The mass balance log and Makie visualiser both use `vol_removed`:
```julia
input_vol  = cumulative_input(all_sources, t) + rainfall_rate * t * sum_valid_areas
domain_vol = sum(state.volume)
mb_err     = input_vol - domain_vol - state.vol_removed   # should be ~0
```

---

## 5. New File Structure

```
FloodA5/
├── FloodModel.jl                    # modified (CLI flags, unified source dispatch,
│                                    #   boundary init, Phase D wiring)
├── surfacewater/
│   └── flow2d.jl                    # modified (_bates_ghost_flux, _manning_ghost_flux,
│                                    #   _ghost_wse)
├── boundaryinputs/                  # NEW directory
│   ├── sources.jl                   # AbstractSource, InflowPoint, UniformRainfall,
│   │                                #   _interp_hydrograph, cumulative_volume dispatch
│   ├── boundary_conditions.jl       # BCType, BoundarySegment, GhostEdge, GhostEdgeList,
│   │                                #   _build_ghost_edges, _match_geojson_segments
│   └── timeseries_io.jl             # LisfloodBDYReader, TwoColumnCSVReader,
│                                    #   read_timeseries dispatch, parse_bci_file
└── test/
    ├── test_inflow_point.jl         # NEW — unit tests T-IP1..T-IP7
    ├── test_open_boundary.jl        # NEW — unit tests T-BC1..T-BC6
    └── (existing tests unchanged)
```

Included from `FloodModel.jl`:
```julia
include(joinpath(@__DIR__, "boundaryinputs", "sources.jl"))
include(joinpath(@__DIR__, "boundaryinputs", "boundary_conditions.jl"))
include(joinpath(@__DIR__, "boundaryinputs", "timeseries_io.jl"))
```

---

## 6. Changes to Existing Files — Summary

### `FlowState` struct (FloodModel.jl)

New fields: `boundary_mask`, `ghost_edges`, `ghost_cell_bc`, `vol_removed`.

### `initialise_flow_model` (FloodModel.jl)

After `_build_edge_list`:
```julia
boundary_mask, ghost_edges, ghost_cell_bc =
    _build_ghost_edges(mesh, edges, elevations, sgs_tables, n, default_bc_type)
```

### `run_flood_model` (FloodModel.jl)

New keyword arguments: `inflow_specs`, `inflow_bci_path`, `bc_file_path`, `closed_boundaries`, `bc_epsg`.

The existing `injection_specs` and `rainpoint_specs` are retained (same CLI flags) but their types are unified into `all_sources :: Vector{AbstractSource}` before being passed to `run_simulation!`.

### `run_simulation!` (FloodModel.jl)

- `injection_points` and `rain_points` parameters replaced by `all_sources :: Vector{AbstractSource}`
- New: `default_bc_type :: BCType = ZeroGradient`
- Phase D boundary outflow pass added after Phase C
- Mass balance logging updated to include `state.vol_removed`
- `vol_removed` tracking in Makie visualiser wired to `state.vol_removed`

### `has_water` guard (FloodModel.jl ~line 1757)

```julia
has_water = rainfall_rate > 0.0 ||
            !isempty(all_sources) ||
            default_bc_type != Closed
```

### `flow2d.jl`

New functions: `_ghost_wse`, `_bates_ghost_flux`, `_manning_ghost_flux`.

---

## 7. Testing Plan

### Unit tests: `test/test_inflow_point.jl`

| ID | Test |
|----|------|
| T-IP1 | `_interp_hydrograph` at knot points returns exact values |
| T-IP2 | Linear interpolation between knots: midpoint gives arithmetic mean |
| T-IP3 | Flat extrapolation before first knot: returns Q[1] |
| T-IP4 | Flat extrapolation after last knot: returns Q[end] |
| T-IP5 | `apply_source!` on `InflowPoint` adds `Q(t) × dt` to correct cell |
| T-IP6 | `cumulative_volume` matches trapezoidal integral of hydrograph |
| T-IP7 | `.bdy` file with `hours` time unit: times correctly converted to seconds |
| T-IP8 | `.bdy` file with two series: both parsed, names correct |
| T-IP9 | `.bci` file with `QVAR` entry: resolves cell, links to correct series |
| T-IP10 | `.bci` file with `QFIX` entry: creates InjectionPoint with correct rate |
| T-IP11 | `.bci` with `--bc-epsg 27700`: projected coordinates converted to lon/lat |
| T-IP12 | Two InflowPoints on same cell: volumes sum correctly |

### Unit tests: `test/test_open_boundary.jl`

| ID | Test |
|----|------|
| T-BC1 | `_build_ghost_edges`: correct count of ghost edges for known mesh |
| T-BC2 | Ghost edge widths match actual polygon edge lengths (not mean) |
| T-BC3 | `_ghost_wse` ZeroGradient: returns `wse_ci` |
| T-BC4 | `_ghost_wse` Critical: returns `sill + (2/3)(wse_ci - sill)` |
| T-BC5 | `_bates_ghost_flux` never returns negative Q (no inflow from ghost) |
| T-BC6 | Closed BC: `_bates_ghost_flux` returns `(0.0, 0.0)` |
| T-BC7 | Mass balance on open flat domain: `input = domain + vol_removed` to <0.01% |
| T-BC8 | `--closed-boundaries` flag: `vol_removed = 0` throughout |
| T-BC9 | GeoJSON bc-file: Closed segment overrides default ZeroGradient correctly |
| T-BC10 | Momentum persistence: `flux_prev` on ghost edge non-zero after first outflow step |

### Regression: synthetic DEM T4 (updated)

T4 mass balance criterion updated to:
```julia
@test abs(injected - (dn_vol + up_vol + vol_removed)) / injected < 1e-4
```
Both SGS and standard flow paths tested.

### Integration: Carlisle res 14

- 1h 50mm/hr rainfall, open boundaries: `vol_removed` grows monotonically, no spurious reflections at edges, `mb_err < 0.1%`
- Same run with `--closed-boundaries`: `vol_removed = 0`, `mb_err < 0.01%`
- Inflow test: `--inflow-point` with a synthetic ramp hydrograph (0→10→0 m³/s over 3h): total injected volume ≈ `cumulative_volume` to <0.1%

---

## 8. Implementation Sequence

1. Create `boundaryinputs/` directory; add stub includes to `FloodModel.jl` (no behaviour change yet)
2. Add `AbstractSource` supertype; promote `InjectionPoint` and `RainPoint`; unify source loop in `run_simulation!`
3. Implement `InflowPoint`, `_interp_hydrograph`, `cumulative_volume` — unit test T-IP1..T-IP6
4. Implement `LisfloodBDYReader` and `TwoColumnCSVReader` — unit test T-IP7..T-IP8
5. Implement `parse_bci_file`; wire `--inflow-point` and `--inflow-bci` CLI flags — T-IP9..T-IP12
6. Add `boundary_mask`, `ghost_edges`, `ghost_cell_bc`, `vol_removed` to `FlowState`
7. Implement `_build_ghost_edges` with actual polygon edge-length geometry — T-BC1..T-BC2
8. Add `_ghost_wse`, `_bates_ghost_flux`, `_manning_ghost_flux` to `flow2d.jl` — T-BC3..T-BC6
9. Wire Phase D into `step_standard!` and `step_sgs!`; update mass balance logging — T-BC7..T-BC8
10. Implement GeoJSON `--bc-file` parsing and `_match_geojson_segments` — T-BC9
11. Implement `--closed-boundaries` flag and `--bc-epsg` flag
12. Update T4 criterion; run full test suite including Carlisle integration

---

## 9. Future Development Hooks

| Future feature | Hook in this design |
|---|---|
| Tidal WSE BC | `BCType.FixedWSE`; `BoundarySegment.wse_series`; `_ghost_wse` FixedWSE branch |
| Tide `.bci` `HFIX` entries | `parse_bci_file` already stubs this; reader needs a `.bdy`-style WSE series |
| Spatial rainfall fields | New `SpatialRainfall <: AbstractSource` in `sources.jl` |
| Infiltration | Separate `surfacewater/infiltration.jl`; volume sink, not a source |
| Evaporation | Same pattern |
| Soil/friction input | Static raster, sampled at init like existing friction raster |
| HydroML / WaterML 2.0 | New `WaterML2Reader <: AbstractTimeSeriesReader` in `timeseries_io.jl` |
| CF-NetCDF gauge input | Same |
| Absolute-time hydrographs | `_interp_hydrograph` `t_offset` kwarg; ISO 8601 CLI flag |
| Phase 3 AMR | Boundary detection uses EdgeList `n_neighbours` — already multi-resolution safe |
| Critical-depth BC | `BCType.Critical` already in enum; `_ghost_wse` branch ready |
| QVAR on N/E/S/W segments | `parse_bci_file` distributes Q uniformly across matched cells |

---

*FloodA5 — University of Canterbury | Branch: inflows_and_bcs | 2026-06-10*
