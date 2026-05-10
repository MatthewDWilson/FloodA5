# FloodA5 — Project State Summary
_Last updated: 2026-03-31 (sublattice fix). Paste this into a new conversation to resume._

---

## Project Overview

A Julia flood model using the **A5 pentagonal DGGS**. Cells are pentagons
(5 neighbours). No native Julia A5 package — Python (`pya5`) used via subprocess.

Two visualisation backends: **Cesium** (CesiumJS web viewer, binary wire protocol)
and **Makie** (GLMakie native window). Both are opt-in via `--vis [mode]`.

**Location:** `F:\OneDrive - University of Canterbury\Julia\FloodA5\`

---

## File Structure

```
FloodA5/
├── FloodModel.jl           # Entry point, flow model, CLI
├── A5Grid.jl               # Mesh generation + cell query API
├── VisualisationServer.jl  # HTTP + WebSocket server (Oxygen.jl)
├── MakieVisualiser.jl      # GLMakie native desktop viewer
├── a5_bridge.py            # Python bridge: pya5, geoparquet
├── benchmark_pip.jl        # PIP benchmark
├── viz/
│   ├── index.html          # CesiumJS viewer
│   ├── config.json         # Runtime config — gitignored, create from example
│   └── config.example.json # Committed template
└── data/
    ├── christchurch_aoi.geojson
    └── example_aoi.geojson
```

---

## Architecture

```
Julia (FloodModel.jl)
  ├── A5Grid.jl
  │     ├── Phase 1: PIP sampling [GPU (CUDA) or CPU threads]
  │     └── Phase 2: subprocess → a5_bridge.py
  │
  ├── VisualisationServer.jl  (--vis cesium)
  │     Endpoints:
  │       GET /viz/{file}           static files
  │       GET /mesh                 GeoJSON + cell_order array
  │       GET /frames/count         {count}
  │       GET /frames/{idx}         {t, vars:[names]}  metadata only
  │       GET /frames/{idx}/{var}   raw Float32 binary, n_cells×4 bytes
  │       GET /status               diagnostics
  │       WS  /live                 JSON notifications
  │     WS message types:
  │       mesh        {type, data:GeoJSON, cell_order:[ids]}
  │       framecount  {type, count, vars:[names]}
  │       newframe    {type, idx, t, vars:[names]}
  │       simcomplete {type, frames}
  │
  └── MakieVisualiser.jl  (--vis makie)
        GLMakie Figure: poly!, Colorbar, Menu dropdown, sidebar text
```

---

## CLI — Full Reference

```
julia [--threads auto] FloodModel.jl
    --meshgen FILE  --meshres N  [--meshout FILE]   # generate mesh
    OR
    --meshload FILE                                  # load saved mesh

    [--dem FILE]  [--dem-strict]
    [--dem-method mean|centroid]  [--dem-samples N]  [--dem-seed N]

    [--flow-model sgs|standard]
    [--sgs-bins N]  [--sgs-samples N]
    [--manning-n N]  [--friction FILE]
    [--sim-duration S]  [--dt-max S]  [--rainfall R]

    [--output FILE]  [--output-interval S]

    [--vis [cesium|makie]]  [--vis-port PORT]

    [--help | -h]
```

**Defaults:** `--flow-model sgs`, `--sim-duration 3600`, `--dt-max 60`,
`--manning-n 0.03`, `--dem-method mean`, `--dem-samples 256`,
`--vis-port 8080`, `--output-interval 60`, rainfall `0`.

---

## VisualisationServer.jl

### Option B wire protocol (binary per-variable)

Frames store a `Dict{String, Vector{Float32}}` — one dense array per variable
in mesh-cell index order. The client fetches one variable at a time:

```
GET /frames/{idx}/{varname}  →  application/octet-stream
                                Float32 LE, n_cells × 4 bytes
```

Client decodes: `new Float32Array(await resp.arrayBuffer())`

Scaling: at 1M cells, one variable = ~4 MB. All four variables = ~16 MB.
Only the actively displayed variable is fetched.

### push_frame! signature

```julia
VisualisationServer.push_frame!(server, t, Dict(
    "depth"      => Float32.(state.water_depth),
    "saturation" => Float32.(sat),
    "volume"     => Float32.(state.volume),
    "velocity"   => Float32.(state.velocity),
))
```

Adding a new variable (e.g. contaminant) requires only adding it to this dict —
no other server or client changes needed. Unknown variable names auto-appear
in the browser dropdown.

### set_mesh! signature

```julia
VisualisationServer.set_mesh!(server, geojson_string, [c.id for c in mesh.cells])
```

The cell ID list is sent to the client as `cell_order` so frames can be
decoded as dense arrays without re-parsing GeoJSON.

### notify_complete!

```julia
VisualisationServer.notify_complete!(vis)
```

Broadcasts `{type:"simcomplete", frames:N}`. Called from `run_flood_model`
immediately after `run_simulation!` returns. Client disables live button,
sets status to "complete", clears `pendingFrame`.

### Exports

```julia
export VisServer, start, stop, push_frame!, set_mesh!, notify_complete!
```

---

## CesiumJS Viewer (index.html)

### Config file

`viz/config.json` (gitignored). Copy from `config.example.json`:

```json
{
  "cesium_ion_token": "YOUR_TOKEN",
  "home_lon": 172.636,
  "home_lat": -43.531,
  "home_alt": 25000,
  "default_basemap": "google3d"
}
```

Loaded at startup via `loadConfig()`. Falls back to defaults if missing.
Token validated against `api.cesium.com/v1/assets` before init.

### Basemap modes

| Value | Description |
|-------|-------------|
| `google3d` | Google Photorealistic 3D Tiles (Ion asset 2275207) |
| `osm` | OpenStreetMap imagery on ellipsoid globe |
| `none` | No imagery (dark background) |

Switched live via panel buttons — `setBasemap(name)` tears down the current
tileset/imagery layer and loads the new one without page reload.

### Variable selector

Dropdown populated from variable names in `framecount` / `newframe` WS messages.
New server-side variables auto-appear via `varDef(name)` fallback.

Known variables and their colormaps:

```javascript
const VAR_DEFS = {
  depth:      { colormap: 'turbo',  fixedMax: null },
  saturation: { colormap: 'blues',  fixedMax: 1.0  },
  volume:     { colormap: 'turbo',  fixedMax: null },
  velocity:   { colormap: 'plasma', fixedMax: null },
};
```

### Render pipeline (flash-free design)

All async work (meta fetch + binary fetch) completes before any scene update.
Single combined render pass — no intermediate grey flash.

Key implementation details:
- Each entity has `_colorProp` (a `ColorMaterialProperty`) created once at mesh load
- Per-frame updates call `ent._colorProp.color.setValue(col)` — mutates in place, zero allocation
- `DRY_COLOR` and `SCRATCH_COLOR` are module-level pre-allocated `Cesium.Color` objects
- Colourmap functions use `Cesium.Color.fromBytes(..., SCRATCH_COLOR)` — zero heap allocation per cell
- `cellPrevValues` Float32Array tracks last-rendered value per cell — unchanged cells are skipped
- `displaySettingsChanged = true` (starts true) forces full redraw after settings changes or mesh load

**Known outstanding issues (in progress):**
- Flash on frame updates still present in some conditions
- `simcomplete` signal not always stopping live frame updates

### State variables

```javascript
let renderInProgress    = false;  // async render in flight
let pendingFrame        = null;   // {idx, varName} latest-wins queue
let simComplete         = false;  // set permanently on simcomplete
let cellPrevValues      = [];     // Float32Array, one per cell
let displaySettingsChanged = true; // force full redraw on next applyFrame
```

---

## MakieVisualiser.jl

### Public API

```julia
vis = MakieVisualiser.start(mesh; title="FloodA5 …")
MakieVisualiser.push_frame!(vis, cell_ids, depths, saturations, volumes, velocities, t)
MakieVisualiser.stop(vis)
```

### Layout

1200 × 760 window:
- `fig[1,1]` — map axis with `poly!`, depth Observable, `:turbo` colormap
- `fig[2,1]` — horizontal colorbar, auto-scaling
- `fig[3,1]` — variable `Menu` dropdown + label
- `fig[1:3,2]` — monospace diagnostics sidebar

### Key implementation details

- `colsize!(fig.layout, ...)` — must target `.layout` not `fig`
- `poly!` takes `Vector{Vector{Point2f}}` — one closed ring per cell
- All mutable state in `Observable`s — `push_frame!` writes, GLMakie re-renders
- `display(fig)` is non-blocking
- `vis.running[]` polled by keep-alive loop; also set `false` on window close

---

## FlowState struct

```julia
mutable struct FlowState
    cell_ids    :: Vector{String}
    water_depth :: Vector{Float64}   # m above local bed (diagnostic)
    volume      :: Vector{Float64}   # m³ stored — primary state variable
    velocity    :: Vector{Float64}   # m/s scalar magnitude (diagnostic)
    elevation   :: Vector{Float64}   # bed elevation above datum (m)
    manning_n   :: Vector{Float64}   # Manning's roughness per cell
    cell_area   :: Vector{Float64}   # plan area (m²)
    adjacency   :: Dict{String, Vector{String}}
    adj_matrix  :: Matrix{Int}       # (max_nb × n_cells)
    edge_width  :: Matrix{Float64}   # shared edge length (m)
    edge_sill   :: Matrix{Float64}   # sill elevation (m)
    centre_dist :: Matrix{Float64}   # centre-to-centre haversine distance (m)
    edge_flux   :: Matrix{Float64}   # unit discharge q (m²/s) at t-dt
    sgs_tables  :: Vector{Any}       # SGSTable per cell (SGSFlow) or empty
end
```

**Note:** the `sgs_tables` line must end with `empty` then `end` on the next line.
The closing `end` has previously been accidentally concatenated (`or emptyend`) —
check line ~118 after any manual edit.

---

## Physics — Bates et al. (2010) Inertial Formulation

Both `step_standard!` and `step_sgs!` use `_bates_flux`:

```
q^t = [ q^{t-dt} - g·h_flow·dt·dWSE/L ]
      / [ 1 + g·h_flow·dt·n²·|q^{t-dt}| / h_flow^(10/3) ]
Q^t = q^t · width
```

- `h_flow = max(WSE_i, WSE_j) - z_sill`
- `L` = `centre_dist[slot,i]` (haversine, metres)
- `q^{t-dt}` = `edge_flux[slot,i]` (stored and updated each step)

**Known issue:** A5 pentagons are non-orthogonal — the shared edge is not
perpendicular to the centre-to-centre vector. A projection correction
(cos θ factor) would improve accuracy on skewed cell pairs. Deferred pending
analysis of typical skew at resolution 14.

---

## HDF5 Output Structure

```
/mesh/   cell_ids, elevations, center_lons, center_lats, [static vars]
/frames/ /000001/  t, water_depth, volume, saturation, velocity
         /000002/  ...
```

Datasets: `chunk=(min(n_cells, 4096),)`, `deflate=4`.

---

## Key Bugs Fixed (cumulative)

| # | Bug | Fix |
|---|-----|-----|
| 1–8 | Various coordinate, mesh, and Python bridge issues | See earlier sessions |
| 9 | `FlowState` not defined | Struct before functions |
| 10 | `EADDRINUSE` | Pre-check port |
| 11–17 | Cesium rendering issues | Globe/tiles/depth-test fixes |
| 18 | `colsize!(fig,…)` MethodError | Use `fig.layout` |
| 19 | `emptyend` parse error line ~118 | `sgs_tables …or empty\nend` |
| 20 | Binary frame zeros | `Vector{UInt8}(reinterpret(UInt8, data))` |
| 21 | Wrong cell colours | `_colorProp.color.setValue()` not raw assignment |
| 22 | `set_mesh!` MethodError | Pass `cell_ids` as third argument |
| 23 | `push_frame!` MethodError | New signature: `(server, t, Dict(...))` |
| 24 | **No inter-cell flow** — `j == 0 && break` in `step_standard!` and `step_sgs!` exited the slot loop at the first empty slot, silently skipping all subsequent valid neighbours (boundary cells produce sparse `adj_matrix` slots) | Changed to `j == 0 && continue` in both step functions |
| 25 | `cos_theta` NaN not guarded — a missing shared edge yields NaN from `_edge_cos_theta`, which was not caught before being passed to `_bates_flux` | Added `isnan(cos_theta)` to the NaN guard in both step functions |
| 26 | **All flux suppressed by cos_theta NaN** — `_edge_cos_theta` returned `NaN` (instead of a safe fallback) when `_shared_edge` found no shared vertices; this triggered the NaN guard in both step functions and silently zeroed all inter-cell flux. Combined with bug 24, this caused the no-flow symptom. | `_edge_cos_theta` now returns `1.0` (orthogonal fallback) instead of `NaN`; `_build_edge_geometry!` adds a second safety net converting any residual NaN to `1.0` |
| 27 | `main()` called unconditionally — `include("FloodModel.jl")` in tests triggered the CLI entry point, crashing with `--meshgen or --meshload` required | Wrapped `main()` in `if abspath(PROGRAM_FILE) == @__FILE__` guard |
| 28 | **Double-counted edge fluxes** — both step functions visited each edge twice (once from i, once from j), producing nearly-cancelling contributions with near-uniform WSE | Added `j <= i && continue` so each edge is processed exactly once from the lower-index cell |
| 31 | **`j<=i` guard order-dependent** — cells whose A5 neighbours all have smaller array indices had all their edges silently skipped | Removed `j<=i` guard; implemented **Option D EdgeList** — flat array of undirected edges, each stored once with `cell_i < cell_j`. Flux loop iterates `1:n_edges` exactly once per edge. No halving factor needed. `FlowState` gains `edges::EdgeList`; old per-slot matrices (`edge_width`, `edge_sill`, `centre_dist`, `edge_cos_theta`, `edge_flux`) removed. `adj_matrix` retained for neighbour queries. |
| 32 | **Zero edges built for real mesh** — cell IDs in parquet are written by pya5's `u64_to_hex` which may omit leading zeros (e.g. `"8a2a1072b59ffff"` 15 chars), while `grid_disk_neighbours` returns IDs via Julia's `_to_hex` always padded to 16 chars (`"08a2a1072b59ffff"`). The `filter(id -> id in all_ids, nbrs)` check in `_build_adjacency_grid_disk` found no matches, giving empty `adj` and zero edges. Model ran but produced no inter-cell flux. | Normalise all cell IDs to 16-char zero-padded hex via `A5Grid._to_hex` at the start of `initialise_flow_model` and inside `_build_adjacency_grid_disk`, `_build_adjacency_matrix!`, and `_build_edge_list`. |
| 33 | **`test_a5grid.jl` name mismatches** — test called `lon_lat_to_cell` and `cell_to_lon_lat`; exported names are `lonlat_to_cell` and `cell_to_lonlat` | Fixed in test file |
| 34 | **`grid_disk_neighbours` returns wrong resolution** — `pya5.grid_disk(cell, 1)` returns a *compact* mixed-resolution result; neighbours were at resolution 15 (`8e...`) while mesh cells are at resolution 16 (`9e...`). The `filter(id -> id in all_ids_norm, nbrs)` check found zero matches in every case, giving empty adjacency and zero edges even after the ID normalisation fix. | Fixed `grid_disk_neighbours` to call `pya5.uncompact(disk, resolution)` before filtering, expanding the compact disk to the target resolution. This is the same pattern used in `a5_bridge.py` for `fill_polygon`. |
| 35 | **A5 dual-sublattice: mesh missing half its cells** — A5 pentagons tile in TWO interleaved sublattices at each resolution. `fill_polygon+uncompact` returns only one sublattice (e.g. `9e...` cells). Their edge neighbours belong to the complementary sublattice (`8e...` cells) which `fill_polygon` never returns. No `9e...` cell has any `9e...` neighbour — `adj` is always empty and the edge list has zero entries. This has been the root cause of zero flux throughout all model runs. | Fixed `generate_mesh` in `a5_bridge.py`: after `fill_polygon+uncompact`, collect all `grid_disk(cell,1)` neighbours of the primary cells, filter to those whose centres lie inside the AOI, and add them. This includes both sublattices in every mesh. Existing meshes in parquet files must be regenerated with `--meshgen` to pick up the fix. |
| 29 | **Inverted continuity sign** — `dV[i] -= vol; dV[j] += vol` with `vol = Q*dt` gave the wrong direction: `_bates_flux` returns Q < 0 when flow goes i→j (wse_i > wse_j), so `dV[i] -= negative = dV[i] increases` (cell i gained volume when it should lose it). Result: water flowed uphill. | Changed to `dV[i] += vol; dV[j] -= vol` matching the sign convention |
| 30 | **Per-edge volume limiter** — limiter applied 5× per cell (once per edge) allowed up to 250% drainage per step | Moved limiter into `_apply_dV_standard!` / `_apply_dV_sgs!` as a single cell-level post-hoc cap at 50% of current volume, applied after all edge fluxes accumulate |

---

## Python Environment

```
Python:   C:\Python311\python.exe
Packages: pya5, geopandas, pyarrow, shapely≥2.0, numpy
```

---

## Julia Packages

```julia
Pkg.add(["PyCall", "JSON3", "DataFrames", "CUDA", "BenchmarkTools",
         "Oxygen", "HTTP", "GLMakie", "HDF5", "Statistics"])
```

---

## Quick Reference

```bash
# Full run — generate mesh, SGS, Cesium viewer
julia --threads auto FloodModel.jl \
    --meshgen christchurch_aoi.geojson --meshres 14 --meshout mesh_sgs.parquet \
    --dem linz_dem.tif --rainfall 30 --sim-duration 3600 --vis

# Repeat run — load mesh, Makie window
julia --threads auto FloodModel.jl \
    --meshload mesh_sgs.parquet --rainfall 10 --vis makie

# No vis
julia --threads auto FloodModel.jl --meshload mesh_sgs.parquet

# Help
julia FloodModel.jl --help

# Check Python
python a5_bridge.py check
```

---

## Pending Work

### Visualisation (deferred — Makie available as fallback)
- CesiumJS flash on frame updates — flagged for future session
- `simcomplete` signal not reliably stopping live updates — flagged for future session

### Phase 2 — Boundary conditions
- `TimeSeries` struct: tiered input (CSV native, WaterML 2.0, CF-NetCDF)
- Upstream hydrograph: time-varying inflow injected as volume source
- Open/tidal BC: fixed WSE at boundary cells, ghost-cell enforced

### Phase 3 — Multi-resolution AMR
- GPU data layout review at Phase 3 start — Strategy 1 preferred (solver on GPU,
  AMR bookkeeping on CPU), with conservative refinement to minimise sync events.
- MultiResMesh data model: multi-level cell store, parent/child A5 IDs
- Coarse/fine flux interface: ghost-cell projection + buffer zone
- Static MR mesh (user-specified regions) before dynamic adaptation
- Refinement triggers: wet/dry front + WSE gradient + momentum threshold
- Coarsening criterion: all-5-children stable for N timesteps (dry = sufficient alone)
- Dynamic AMR: runtime refine/coarsen with CFL-aware dt update

### Phase 4 — Output and validation
- NetCDF output: CF-1.8 + UGRID unstructured mesh convention
- MR-aware visualisation: show cell resolution level in Cesium and Makie
- Validation suite: Thacker planar surface, circular dam-break benchmarks

### Phase 5 — Future physics and performance
- Riemann solver: swap `_bates_flux` kernel; ghost-cell interface already compatible.
  Full over-relaxed non-orthogonality correction as upgrade path (see Weller 2014,
  Gemini conversation 2026-03-25).  LSQ gradient reconstruction also deferred here.
- GPU solver kernel for `step_standard!` / `step_sgs!`
- Ensemble / calibration: Manning sweep, skill scores
