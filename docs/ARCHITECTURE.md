# FloodA5 — Architecture

This document is for developers who want to understand how FloodA5 is structured,
where to find things, and how to extend the model. For the physics see
[METHODS.md](METHODS.md); for running the model see [USER_GUIDE.md](USER_GUIDE.md).

---

## 1. Overview

FloodA5 is a Julia application with a Python subprocess for grid operations. The
split is deliberate: Julia provides multithreaded performance for the simulation
loop; Python provides access to the `pya5` A5 DGGS library, which has no Julia
binding.

```
FloodModel.jl          Entry point, CLI, simulation loop, FlowState
mesh/A5Grid.jl         Julia module: mesh struct, DEM sampling, SGS tables
mesh/a5_bridge.py      Python subprocess: pya5 calls, GeoParquet I/O
surfacewater/flow2d.jl Pure physics kernels (Bates, Manning R-A, CFL)
boundaryinputs/        Source types, BC types, hydrograph I/O
visualisation/         CesiumJS and GLMakie backends, post-processing tools
```

All Julia modules are included (not imported) by `FloodModel.jl` using `include()`.
There are no registered Julia packages — the entire application lives in the
repository and is run directly with `julia FloodModel.jl`.

---

## 2. Module Responsibilities

### `FloodModel.jl`

The entry point and main application file. Responsibilities:

- **Type definitions** — `FlowState`, `EdgeList`, `InjectionPoint`, `RainPoint`,
  `AbstractSource`, `GhostEdge`, `BCType`, and the `FlowMethod` hierarchy
  (`StandardFlow`, `SGSFlow`). These are defined here rather than in their
  respective modules because `FlowState` references them and Julia requires
  types to be defined before use.
- **CLI parsing** — `_pop_flag` / `_pop_bool` helpers and the `main()` function.
  The `print_help()` function is the authoritative CLI reference.
- **Mesh initialisation** — `initialise_flow_model(mesh, method)` builds the
  `EdgeList`, assigns SGS tables, detects boundary cells, builds ghost edges, and
  resolves source locations.
- **Simulation loop** — `run_simulation!(state, method, sources, ...)` drives
  the timestep loop. Per-step phases:
  - Phase A: apply water sources
  - Phase B: sync `water_depth` from `volume` (standard solver)
  - Phase C: compute edge fluxes and accumulate `dV`
  - Phase D: apply volume limiter and scatter `dV` to `volume`
  - Phase E: apply ghost-edge (boundary outflow) fluxes
  - Phase F: log progress and write output frames

The step functions `step_standard!` and `step_sgs!` implement Phases C–D for
their respective solvers.

**Include order matters.** `flow2d.jl` is included after type definitions but
before the step functions, because the step functions call kernels defined there.
`boundaryinputs/` files are included after `InjectionPoint` and `RainPoint` are
defined (sources dispatch on these types). The include order is:

```julia
include("mesh/A5Grid.jl")
include("visualisation/VisualisationServer.jl")
include("visualisation/MakieVisualiser.jl")
# ... type definitions ...
include("surfacewater/flow2d.jl")
include("boundaryinputs/boundary_conditions.jl")
include("boundaryinputs/sources.jl")
include("boundaryinputs/timeseries_io.jl")
```

### `mesh/A5Grid.jl`

Julia module (`module A5Grid ... end`). Provides the mesh API used by
`FloodModel.jl`. Key responsibilities:

- **`A5Mesh` and `A5Cell` structs** — the in-memory mesh representation.
  `A5Mesh.cells` is a `Vector{A5Cell}`; static variables (elevation, friction)
  live in `A5Mesh.static_vars :: Dict{String, Vector{Float64}}`. Array variables
  (SGS curves) live in `A5Mesh.array_vars :: Dict{String, Matrix{Float64}}`.
- **Mesh generation** — `mesh_for_aoi(geojson, resolution)` delegates to
  `a5_bridge.py` via subprocess. Python handles all `pya5` calls and writes the
  GeoParquet file; Julia reads it back.
- **DEM sampling** — `sample_dem_mean!` and `sample_dem_centroid!` use ArchGDAL
  for GDAL-based raster access and CRS reprojection. Halton points are generated
  in `_halton(base, n, offset)`.
- **SGS pre-processing** — `build_sgs_tables!` builds hypsometric curves for
  cells and edges. The cell PIP step runs on GPU (CUDA) when available.
- **SGS runtime lookups** — `wse_from_volume`, `wetted_area_from_wse`,
  `flow_area_from_wse`, `wetted_perim_from_wse`, `hydraulic_radius_from_wse`.
  All are `@inline` for zero-overhead calls in the inner loop.
- **Cell query API** — `lonlat_to_cell`, `cell_to_lonlat`, `cell_to_boundary`,
  `cell_to_children`, `cell_to_parent`, `get_resolution`. These call `pya5` via
  PyCall with a `ReentrantLock` for thread safety.

### `mesh/a5_bridge.py`

Python subprocess called by `A5Grid.jl` for mesh generation. It:

1. Calls `pya5.fill_polygon` + `pya5.uncompact` to enumerate cells at the target
   resolution, covering both A5 sublattices by adding `grid_disk` neighbours.
2. Computes boundary polygons via `pya5.cell_to_boundary` for each cell.
3. Normalises all cell IDs to 16-character zero-padded hex and all longitudes to
   [-180, 180].
4. Writes the mesh as GeoParquet using geopandas and pyarrow.

The bridge is invoked as `python mesh/a5_bridge.py mesh_for_aoi <aoi_wkt> <res> <out>`.
A smoke-test mode is available: `python mesh/a5_bridge.py check`.

**Why subprocess and not PyCall?** PyCall can call `pya5` for single-cell queries
(with a GIL lock), but `pya5`'s batch operations interact poorly with PyCall's
garbage collector in multithreaded Julia contexts. The subprocess boundary provides
complete isolation.

### `surfacewater/flow2d.jl`

Pure physics kernels — no `FlowState`, no mesh, no I/O. All functions are
`@inline` and operate on scalar `Float64` arguments only.

Key functions:

| Function | Description |
|---|---|
| `_bates_flux` | Bates eq. 9 inertial formulation. Base flux kernel. |
| `_bates_flux_limited` | `_bates_flux` + Froude limiter + volume limiter + consistent `q_stored`. Used by `step_standard!`. |
| `_manning_flux_ra` | Manning R-A formulation (`A`, `R` based). Used by `step_sgs!`. |
| `_ghost_wse` | Ghost-cell WSE for each BC type. |
| `_bates_ghost_flux` | Outflow flux across a ghost edge (standard solver). |
| `_manning_ghost_flux` | Outflow flux across a ghost edge (SGS solver). |
| `_cfl_dt` | Wave-speed CFL timestep. |
| `_q_centred` | Q-centred (spatially smoothed) unit discharge. |

Constants: `_G = 9.81`, `HFLOW_THRESHOLD = 0.001`, `FROUDE_LIMIT = 0.8`,
`Q_CENTRE_THETA = 0.9`, `N_SIDES = 5`, `DONOR_EDGE_DIVISOR = 10`.

This file is designed to be includable independently (e.g. by a GPU kernel
wrapper or a 3D coupling plugin) without pulling in the full application.

### `boundaryinputs/`

Three files, included after type definitions in `FloodModel.jl`:

**`sources.jl`** — `AbstractSource` abstract type and `InflowPoint` concrete type.
Dispatch methods `apply_source!(state, src, t, dt)` and `cumulative_volume(src, t)`
for all source types. `InjectionPoint` and `RainPoint` are defined in `FloodModel.jl`
as `AbstractSource` subtypes; `InflowPoint` is defined here.

**`boundary_conditions.jl`** — `BCType` enum, `GhostEdge` struct, `BoundarySegment`
struct, `_build_ghost_edges(mesh, edges, ...)`, `_ghost_wse`, `_bates_ghost_flux`,
`_manning_ghost_flux`, `load_bc_file` (GeoJSON reader), `apply_bci_free_entries!`.

**`timeseries_io.jl`** — `AbstractTimeSeriesReader`, `LisfloodBDYReader`,
`TwoColumnCSVReader`, `read_timeseries` dispatch, `parse_bci_file`.

### `visualisation/`

Four scripts and one static directory:

| File | Description |
|---|---|
| `MakieVisualiser.jl` | GLMakie native desktop window. Module `MakieVisualiser`. |
| `VisualisationServer.jl` | Oxygen.jl HTTP + WebSocket server. Module `VisualisationServer`. |
| `view_h5_output.jl` | Standalone post-processing viewer for HDF5 output. |
| `visualise_mesh.jl` | Publication-quality mesh and SGS figures (wire-frame, cross-section, hypsometric curve). |
| `cesium/` | CesiumJS static files (`index.html`, `config.example.json`). |

---

## 3. Key Data Structures

### `FlowState`

The central mutable struct, passed through the entire simulation loop:

```julia
mutable struct FlowState
    cell_ids    :: Vector{String}
    water_depth :: Vector{Float64}   # m — diagnostic; derived from volume each step
    volume      :: Vector{Float64}   # m³ — PRIMARY STATE VARIABLE
    velocity    :: Vector{Float64}   # m/s scalar magnitude
    elevation   :: Vector{Float64}   # bed elevation (m)
    manning_n   :: Vector{Float64}   # Manning's roughness per cell
    cell_area   :: Vector{Float64}   # plan area (m²)
    cell_lons   :: Vector{Float64}   # cell centre longitude (degrees)
    cell_lats   :: Vector{Float64}   # cell centre latitude (degrees)
    adjacency   :: Dict{String, Vector{String}}
    adj_matrix  :: Matrix{Int}       # (5 × n_cells) — for neighbour queries
    edges       :: EdgeList
    sgs_tables  :: Vector{Any}       # Vector{SGSTable} for SGS; empty for standard
    boundary_mask  :: BitVector      # true = domain-edge cell
    ghost_edges    :: Vector{GhostEdge}
    ghost_cell_bc  :: Vector{BCType}
    vol_removed    :: Float64        # cumulative outflow through ghost edges (m³)
end
```

`volume` is the primary state variable. `water_depth` is a derived diagnostic,
recomputed each step.

### `EdgeList`

Undirected edge list — each edge stored once, with `cell_i < cell_j`:

```julia
struct EdgeList
    n_edges   :: Int
    cell_i    :: Vector{Int}       # lower-index cell
    cell_j    :: Vector{Int}       # higher-index cell
    width     :: Vector{Float64}   # shared edge length (m)
    L         :: Vector{Float64}   # centre-to-centre haversine distance (m)
    cos_theta :: Vector{Float64}   # non-orthogonality correction (1.0 = orthogonal)
    sill      :: Vector{Float64}   # sill elevation (m)
    flux      :: Vector{Float64}   # q (m²/s) at t-dt — standard solver
    flux_Q    :: Vector{Float64}   # Q (m³/s) at t-dt — SGS R-A solver
    collinear_i :: Vector{Int}     # index of most collinear edge for cell_i (Q-centred)
    collinear_j :: Vector{Int}     # index of most collinear edge for cell_j
end
```

**Sign convention:** `flux > 0` means flow from `cell_j` to `cell_i` (j is higher).
When `WSE_i > WSE_j`: `dWSE > 0` → `q_new < 0` → `Q < 0` → `dV[i] += Q·dt` (i
loses volume) → correct.

### `SGSTable`

Pre-computed hypsometric curves for one cell:

```julia
struct SGSTable
    elev_bins        :: Vector{Float64}   # WSE knots (n_bins)
    vol_curve        :: Vector{Float64}   # cumulative volume (m³)
    area_curve       :: Vector{Float64}   # wetted plan area (m²)
    cell_area        :: Float64
    z_min            :: Float64
    z_max            :: Float64
    edge_area_curves :: Matrix{Float64}   # (n_bins × 5) flow area per slot (m²)
    edge_perim_curves:: Matrix{Float64}   # (n_bins × 5) wetted perimeter per slot (m)
end
```

---

## 4. Extension Points

### Adding a new water source type

1. Define a new struct `<: AbstractSource` in `boundaryinputs/sources.jl`.
2. Implement `apply_source!(state, src::MySource, t, dt)` — add volume to
   `state.volume[src.cell_index]`.
3. Implement `cumulative_volume(src::MySource, t)` — return the cumulative volume
   added from t=0 to t (for mass balance logging).
4. Add a CLI flag in `FloodModel.jl`'s `main()` and construct the source there.

The simulation loop calls `apply_source!` for every element of `all_sources ::
Vector{AbstractSource}` each step, so no further changes are needed in the loop.

### Adding a new boundary condition type

1. Add a new variant to the `BCType` enum in `boundaryinputs/boundary_conditions.jl`.
2. Add a branch in `_ghost_wse(wse_ci, sill, bc, ...)` to return the appropriate
   ghost-cell WSE.
3. Add the string name to `load_bc_file` and `apply_bci_free_entries!` for CLI/file
   recognition.
4. The flux kernel (`_bates_ghost_flux` / `_manning_ghost_flux`) already passes
   `wse_ghost` from `_ghost_wse` — no change needed there unless the new BC type
   requires a fundamentally different flux computation.

`FixedWSE` (tidal BC) and `FixedQ` stubs are already present in `BCType`; their
`_ghost_wse` branches return `wse_series_val` and are awaiting implementation of the
time-series lookup.

### Adding a new physics kernel

`surfacewater/flow2d.jl` is designed to hold all scalar physics. New kernels
(infiltration, evaporation, a Riemann solver option) should be added here as
`@inline` functions taking only `Float64` scalar arguments. The step functions
in `FloodModel.jl` call these in their inner loops.

A companion file (e.g. `surfacewater/infiltration.jl`) is preferable for larger
additions, following the same pattern: pure scalar kernels, no `FlowState`
dependency, included by `FloodModel.jl`.

### Adding a new time-series format

1. Define a new struct `<: AbstractTimeSeriesReader` in `boundaryinputs/timeseries_io.jl`.
2. Implement `read_timeseries(reader::MyReader) →
   Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}`, where the dict maps
   series name → `(t_seconds, Q_m3s)`.
3. The `InflowPoint` constructor and `parse_bci_file` call `read_timeseries` through
   dispatch — no other changes are needed.

Planned future readers: `WaterML2Reader`, `CFNetCDFGaugeReader`, `HydroMLReader`.

### Adding a new visualisation variable

1. Add the variable to `FlowState` and compute it in the step function.
2. Add the dataset to `_write_frame!` in `FloodModel.jl`.
3. Add the variable to `VAR_DEFS` in `visualisation/cesium/index.html`:

```javascript
const VAR_DEFS = {
    depth:      { colormap: "turbo",  fixedMax: null },
    saturation: { colormap: "blues",  fixedMax: 1.0  },
    // add new variable here
    myvar:      { colormap: "plasma", fixedMax: null },
};
```

Unknown variables sent by the server are handled by a fallback in the CesiumJS
viewer, so new variables will appear in the dropdown automatically. Explicit
entries in `VAR_DEFS` allow you to set a specific colormap and fixed maximum.

---

## 5. Python Bridge Details

### A5 cell ID conventions

- pya5 returns cell IDs as Python `int`. The bridge converts to hex via
  `hex(cell_id).lstrip("0x").zfill(16)` — always 16 characters, zero-padded.
- Julia normalises all IDs to 16-char hex on load: `_norm_id(id) = A5Grid._to_hex(parse(UInt64, id, base=16))`.
- Never compare cell IDs without normalising both sides to 16-char hex.

### Dual sublattice

`pya5.fill_polygon + uncompact` returns cells from only one A5 sublattice. The
bridge adds neighbours of these primary cells (via `grid_disk(cell, 1)`) and
includes those whose centres fall inside the AOI — capturing both sublattices.
Any mesh generated without this step will have zero edges and produce no flux.

### `grid_disk` returns compact results

`pya5.grid_disk(cell, 1)` returns a compact result — a mix of cells at different
resolution levels. Always follow with `pya5.uncompact(disk, target_resolution)`.

---

## 6. Thread Safety

- **Mesh generation** — use `julia --threads 1`. The Python bridge makes PyCall
  calls inside loops; garbage collection on non-main threads can trigger PyCall
  finalisers and crash. This is a known limitation.
- **Simulation loop** — safe with any thread count. The edge flux loop uses
  `Threads.@threads` over the `EdgeList`; there are no write hazards because
  the loop accumulates into per-edge `edge_vol` arrays before scattering to cells.
- **Ghost-edge Phase E** — serial (ghost edges are owned by single cells).
- **Cell queries** (`lonlat_to_cell` etc.) — thread-safe via `ReentrantLock(_pycall_lock)`.

---

## 7. Testing

Tests live in `test/`. Run individual test files with:

```bash
julia --threads auto --project=. test/test_sgs_unit.jl
```

Key test files:

| File | What it tests |
|---|---|
| `test_a5grid.jl` | A5Grid cell API |
| `test_edge_geometry.jl` | `_shared_edge`, `_edge_cos_theta`, `_edge_length_m` |
| `test_sgs_unit.jl` | SGS hypsometric lookups, 5-cell chain routing |
| `test_sgs_edge_ra_tables.jl` | SGS edge flow area / wetted perimeter curves |
| `test_inflow_point.jl` | `InflowPoint`, `_interp_hydrograph`, `.bdy` reader, `.bci` parser |
| `test_open_boundary.jl` | Ghost edge geometry, BC types, outflow flux, mass balance |
| `test_flat_rainpoint.jl` | Flat-terrain point-source integration test |
| `synthetic_dem/test_sgs_synthetic.jl` | Full SGS validation: mesh, routing, sub-cell notch, mass balance |

The synthetic DEM test (`test_sgs_synthetic.jl`) is the most comprehensive
validation: it uses a purpose-built DEM with a Gaussian embankment containing
a sub-cell notch, and confirms that the SGS solver routes 5–14× more water through
the notch than the standard solver at the same injection history.
