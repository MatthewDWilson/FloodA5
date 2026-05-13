# FloodA5 — Project State Summary

_Last updated: 2026-05-12 (hydraulics validation session). Paste this into a new conversation to resume._

---

## Project Overview

A Julia flood model using the **A5 pentagonal DGGS**. Cells are pentagons (5 neighbours each for interior cells). No native Julia A5 package — Python (`pya5`) used via subprocess through `a5_bridge.py`.

Two visualisation backends: **Cesium** (CesiumJS web viewer, binary wire protocol) and **Makie** (GLMakie native window). Both opt-in via `--vis [mode]`.

**Location:** `F:\OneDrive - University of Canterbury\Julia\FloodA5\`
**Julia version:** 1.12  **Python:** 3.11 (`pya5`, `geopandas`, `pyarrow`, `shapely≥2.0`, `numpy`)
**GPU:** CUDA detected — PIP sampling runs on GPU above 500,000 sample-point threshold.

---

## File Structure

```
FloodA5/
├── FloodModel.jl              # Entry point, flow model, CLI, simulation loop
├── A5Grid.jl                  # Mesh generation, cell query API, SGS tables
├── VisualisationServer.jl     # HTTP + WebSocket server (Oxygen.jl)
├── MakieVisualiser.jl         # GLMakie native desktop viewer
├── a5_bridge.py               # Python bridge: pya5, GeoParquet I/O
├── a5_mesh_diagnostic.py      # Mesh diagnostic utilities
├── test_flat_rainpoint.jl     # Flat-terrain point-source validation test
├── viz/
│   ├── index.html             # CesiumJS viewer
│   ├── config.json            # Runtime config — gitignored
│   └── config.example.json   # Committed template
├── test/
│   ├── flat_test_aoi.geojson  # Small AOI for point-source tests (~5.5 km²)
│   ├── flat_mesh_res14.parquet # 61-cell flat mesh (no DEM)
│   ├── kaiapoi_aoi.geojson
│   └── kaiapoi_dem.tif
└── .vscode/
    └── launch.json            # VS Code debug configurations [00]–[10]
```

---

## Architecture

```
Julia (FloodModel.jl)
  ├── A5Grid.jl
  │     ├── Phase 1: PIP sampling [GPU (CUDA) or CPU threads]
  │     └── Phase 2: subprocess → a5_bridge.py
  │                   ├── pya5: fill_polygon + uncompact (both sublattices)
  │                   ├── grid_disk neighbours + uncompact
  │                   └── GeoParquet output
  │
  ├── VisualisationServer.jl  (--vis cesium)
  │     GET /viz/{file}, /mesh, /frames/count, /frames/{idx}, /frames/{idx}/{var}
  │     WS /live — mesh, framecount, newframe, simcomplete
  │
  └── MakieVisualiser.jl  (--vis makie)
        GLMakie Figure: poly!, Colorbar, Menu dropdown, diagnostics sidebar
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
    [--injection-point LAT,LON,RATE_M3S]   # repeatable
    [--rainpoint LAT,LON,RATE_MM_HR]       # repeatable — NEW

    [--output FILE]  [--output-interval S]

    [--vis [cesium|makie]]  [--vis-port PORT]

    [--help | -h]
```

**Defaults:** `--flow-model sgs`, `--sim-duration 3600`, `--dt-max 60`, `--manning-n 0.03`,
`--dem-method mean`, `--dem-samples 256`, `--vis-port 8080`, `--output-interval 60`, rainfall 0.

---

## VS Code Launch Configurations

| # | Name | Purpose |
|---|------|---------|
| [00] | Write flat test AOI GeoJSON | Run once; writes `test/flat_test_aoi.geojson` |
| [01] | Standard mesh + elevation | Kaiapoi, res 16, with DEM |
| [02] | SGS-ready mesh | Kaiapoi, res 16, with DEM + hypsometric tables |
| [03] | Standard simulation | Load kaiapoi standard mesh, Makie viewer |
| [04] | SGS simulation | Load kaiapoi SGS mesh, Makie viewer |
| [05] | Standard mesh, no DEM | Kaiapoi, res 16, no elevation |
| [06] | Mesh + simulation in one run | Kaiapoi standard |
| [07] | Flat mesh gen, res 14 | Generate `test/flat_mesh_res14.parquet` (61 cells) |
| [08] | Flat rainpoint 1hr 50mm/hr | Baseline validation run |
| [09] | Flat rainpoint + Makie | Same as [08] with live viewer |
| [08b] | Flat rainpoint 2hr 50mm/hr | Longer run for spread validation |
| [08c] | Flat rainpoint 1hr 200mm/hr | Accelerated spread test |
| [08d] | Flat rainpoint 1hr 1000mm/hr | Full-domain stress test |
| [10] | Validate HDF5 output | Runs `test_flat_rainpoint.jl --analyse [--file 2hr\|200\|1000]` |
| Python | Mesh diagnostic | `a5_mesh_diagnostic.py` on active file |

**Important:** all Julia configs use `"program": "${workspaceFolder}/FloodModel.jl"` (hardcoded) or the specific script, NOT `"${file}"`. Using `${file}` caused silent failures when the wrong file was active.

---

## FlowState Struct

```julia
mutable struct FlowState
    cell_ids    :: Vector{String}
    water_depth :: Vector{Float64}   # m above local bed (diagnostic only)
    volume      :: Vector{Float64}   # m³ stored — PRIMARY STATE VARIABLE
    velocity    :: Vector{Float64}   # m/s scalar magnitude (NOT YET COMPUTED — always 0)
    elevation   :: Vector{Float64}   # bed elevation above datum (m)
    manning_n   :: Vector{Float64}   # Manning's roughness per cell
    cell_area   :: Vector{Float64}   # plan area (m²) — recomputed from polygon at init
    adjacency   :: Dict{String, Vector{String}}
    adj_matrix  :: Matrix{Int}       # (max_nb × n_cells) — retained for neighbour queries
    edges       :: EdgeList          # undirected edge list (primary flux data structure)
    sgs_tables  :: Vector{Any}       # SGSTable per cell (SGSFlow only) or empty
end
```

**Critical:** `velocity` is always zero — never computed. See Bug 36 below.

---

## EdgeList Struct

```julia
struct EdgeList
    n_edges   :: Int
    cell_i    :: Vector{Int}       # lower-index cell (cell_i < cell_j always)
    cell_j    :: Vector{Int}       # higher-index cell
    width     :: Vector{Float64}   # shared edge length (m)
    L         :: Vector{Float64}   # centre-to-centre haversine distance (m)
    cos_theta :: Vector{Float64}   # non-orthogonality correction (1.0 = orthogonal)
    sill      :: Vector{Float64}   # sill elevation (m)
    flux      :: Vector{Float64}   # q (m²/s) at t-dt, signed: +ve = flow i→j
end
```

**Sign convention:** `flux > 0` → flow from `cell_j` to `cell_i` (j is higher).
When `WSE_i > WSE_j` (i is higher): `dWSE > 0` → `q_new < 0` → `Q < 0` →
`dV[ci] += Q*dt` (ci loses volume) ✓, `dV[cj] -= Q*dt` (cj gains volume) ✓.

---

## Physics — Bates et al. (2010) Inertial Formulation

```
q^t = [ q^{t-dt} - g·h_flow·dt·(dWSE/L_eff) ]
      / [ 1 + g·h_flow·dt·n²·|q^{t-dt}| / h_flow^(10/3) ]
Q^t = q^t · width
```

- `h_flow = max(WSE_i, WSE_j) - z_sill`,  floored at `max(h_flow, 1e-6)` (Bug 38)
- `L_eff  = L × cos_theta`,  `cos_theta` floored at 0.10
- `z_sill = max(elev_i, elev_j)` for standard flow; pre-computed edge minimum for SGS
- `dWSE   = WSE_i - WSE_j`  (positive when i is higher → flow from i to j → Q negative)

**Water depth pre-sync (Bug 37):** After source injection (rainfall/rainpoint/injection),
`water_depth` is synced from `volume / cell_area` before the flux loop, so the first
timestep sees the correct WSE rather than stale zeros.

---

## Source Types

### `--rainfall R` (mm/hr)
Applied uniformly to every cell every timestep. `dV = rate_m_s × dt × cell_area`.

### `--injection-point LAT,LON,RATE` (m³/s)
Fixed volumetric point source at nearest mesh cell. Rate in m³/s.

### `--rainpoint LAT,LON,RATE_MM_HR` (NEW)
Localised rainfall at nearest mesh cell only. Rate in mm/hr, converted to m³/s
using the individual cell area: `rate_m3s = (mm_hr / 3_600_000) × cell_area`.
Unlike `--rainfall`, does not apply to other cells. Repeatable.
`_find_nearest_cell` warns if the requested point is >2 km from the nearest cell.

---

## Mesh Connectivity Check (NEW)

`_check_mesh_connectivity(edges, n_cells)` runs after `_build_edge_list` using BFS.
Logs a warning if the mesh has >1 connected component (isolated cells that will never
receive flux). Logs component sizes and marks the largest as the expected source location.

**Key finding from validation:** Small AOIs at res 14 can produce disconnected meshes
where the two A5 sublattices form separate graph components. Both sublattices are
captured by the mesh generator, but at very small AOIs (~2 km²) the sublattices may
not share physical boundaries and form two isolated clusters. Use AOIs ≥5 km² at
res 14 to ensure full connectivity. The test AOI (~5.5 km²) gives 61 fully-connected cells.

---

## Key Bugs Fixed (cumulative)

| # | Bug | Fix |
|---|-----|-----|
| 1–8 | Various coordinate, mesh, and Python bridge issues | Earlier sessions |
| 9 | `FlowState` not defined | Struct before functions |
| 10 | `EADDRINUSE` | Pre-check port |
| 11–17 | Cesium rendering issues | Globe/tiles/depth-test fixes |
| 18 | `colsize!(fig,…)` MethodError | Use `fig.layout` |
| 19 | `emptyend` parse error | `sgs_tables …or empty\nend` |
| 20 | Binary frame zeros | `Vector{UInt8}(reinterpret(UInt8, data))` |
| 21 | Wrong cell colours | `_colorProp.color.setValue()` |
| 22 | `set_mesh!` MethodError | Pass `cell_ids` as third argument |
| 23 | `push_frame!` MethodError | New signature with Dict |
| 24 | No inter-cell flow — `j==0 && break` | Changed to `continue` |
| 25 | `cos_theta` NaN not guarded | Added `isnan` check |
| 26 | All flux suppressed by cos_theta NaN | `_edge_cos_theta` returns 1.0 fallback |
| 27 | `main()` called on `include` | `abspath(PROGRAM_FILE) == @__FILE__` guard |
| 28 | Double-counted edge fluxes | `j <= i && continue` guard |
| 29 | Inverted continuity sign — water flowed uphill | `dV[i] += vol; dV[j] -= vol` |
| 30 | Per-edge volume limiter (5× over-drainage) | Cell-level post-hoc cap in `_apply_dV_*` |
| 31 | `j<=i` guard order-dependent — edges skipped | EdgeList (Option D): undirected, `cell_i < cell_j` |
| 32 | Zero edges — ID padding mismatch (15 vs 16 chars) | Normalise all IDs to 16-char hex |
| 33 | Test name mismatches | Fixed in test file |
| 34 | `grid_disk` returns wrong resolution | `uncompact(disk, resolution)` after `grid_disk` |
| 35 | Dual-sublattice — mesh missing half its cells | Bridge adds neighbours of primary cells (both sublattices) |
| 36 | **`velocity` never computed — always zero** | Not yet fixed. See Pending Work. |
| 37 | **`water_depth` not synced after source injection** — `wet=0` on first step despite volume > 0; CFL reads stale depth=0 | Added pre-physics-step sync: `water_depth[i] = volume[i] / cell_area[i]` for all cells (standard flow) |
| 38 | **NaN in domain_vol from `h_flow^(10/3)` underflow** — at tiny depths (~4×10⁻⁴ m), `h_flow^(10/3) ≈ 1.4×10⁻¹³`; when `q_prev` was non-zero the denominator became `Inf`, giving `NaN` flux and NaN volumes by step ~50 | Added `h_flow = max(h_flow, 1e-6)` floor in `_bates_flux` |
| 39 | **Wrong cell areas from parquet** — the Python bridge stores `cell_area` column but at small AOI sizes this could be unreliable; `initialise_flow_model` now recognises both `cell_area` and `sgs_cell_area` keys, and for non-SGS meshes recomputes areas from `_polygon_area_m2(c.boundary)` | Recompute from polygon for non-SGS meshes |
| 40 | **`_find_nearest_cell` silent out-of-mesh** — when `--rainpoint` coordinates are outside the mesh AOI, the nearest mesh cell (possibly boundary) is used silently, causing source at wrong location | Added `@warn` when nearest cell is >2 km from requested point |
| 41 | **`launch.json` used `${file}`** — VS Code debug configs used the active file as program, so running [07] while `test_flat_rainpoint.jl` was open ran the wrong script silently | All Julia configs now use `${workspaceFolder}/FloodModel.jl` (hardcoded) |
| 42 | **Julia 1.12 rejects `using` inside function body** — `test_flat_rainpoint.jl` tried lazy `import HDF5` inside `validate_output()`, which the JuliaInterpreter debugger rejects | Moved to conditional top-level `if DO_ANAL; using HDF5; end` |
| 43 | **Julia 1.12 rejects `Main.varname =` in function** — debug code used `Main._step_std_debug_done = true`; Julia 1.12 forbids assigning new globals in `Main` from within a function | Replaced with module-level `const _step_debug_count = Ref{Int}(0)` |
| 44 | **Disconnected mesh components** — small AOIs at res 14 produce two A5 sublattices that form isolated graph components; 12/29 cells permanently unreachable | Added `_check_mesh_connectivity` (BFS); warns with component sizes. Use AOI ≥5 km² at res 14. |
| 45 | **Wrong sign convention comment** — `step_standard!` had inline comment "Q>0 means flow from cell_i to cell_j" which is backwards | Corrected: Q<0 when i is higher (flow i→j), Q>0 when j is higher (flow j→i) |

---

## Validation Results (2026-05-12)

Flat-terrain point-source test (standard flow, no DEM, single `--rainpoint`):

| Run | Mesh | Rate | Duration | Mass balance | Max wet cells | Notes |
|-----|------|------|----------|-------------|---------------|-------|
| [08] | 29 cells (disconnected) | 50 mm/hr | 1 hr | Perfect | 17/29 (59%) | 12 cells isolated |
| [08b] | 29 cells | 50 mm/hr | 2 hr | Perfect | 17/29 | Same isolation issue |
| [08c] | 29 cells→61 cells | 200 mm/hr | 1 hr | Perfect | 17/61 (28%) | Sim too short |
| [08c] | 61 cells (connected) | 200 mm/hr | 1 hr | Perfect | 17/61 | Needs ~1.7hr for ring 3 |
| [08d] | 61 cells (connected) | 1000 mm/hr | 1 hr | Perfect | 33/61 (54%) | Ring 4 reached |

**Mass balance is exact** at all checkpoints across all runs (`domain_vol = rate × t` to <0.002%).
**Flux direction is correct** — water flows from source outward in all runs.
**Spread is rate/time limited** — the wet-cell frontier advances predictably; stalls are due to
insufficient depth at outer ring cells rather than any hydraulic bug.

**Ring cascade timing at 1000 mm/hr:**
- Ring 1 (5 cells): t=120s (step 4)
- Ring 2 (11 cells): t=900s (step 30)
- Ring 3 (~6 cells): t=2100s (step 70)
- Ring 4 (~6 cells): t=2700s (step 90)
- Full 61 cells: estimated ~3 hours at 1000 mm/hr

---

## Pending Work

### Immediate (before further development)
- **Bug 36: implement velocity computation** — `state.velocity` always zero. Should be
  computed from edge fluxes: per-cell magnitude from `Σ|Q|_edges / (cell_area × depth)`.
- **Analytical validation** — compare against Thacker planar surface (exact solution for
  oscillating water surface in parabolic bowl). See Validation Strategy section below.
- **LINZ DEM ingestion** — Kaiapoi domain at res 14 with 1m LiDAR.

### Visualisation (deferred)
- CesiumJS flash on frame updates
- `simcomplete` signal not reliably stopping live updates

### Phase 2 — Boundary conditions
- `TimeSeries` struct: tiered input (CSV, WaterML 2.0, CF-NetCDF)
- Upstream hydrograph: time-varying inflow as volume source
- Open/tidal BC: fixed WSE at boundary cells, ghost-cell enforced
- `--source-duration` flag: stop injection after N seconds

### Phase 3 — Multi-resolution AMR
- GPU data layout review — Strategy 1 preferred (solver on GPU, AMR bookkeeping on CPU)
- `MultiResMesh` data model: multi-level cell store, parent/child A5 IDs
- Coarse/fine flux interface: ghost-cell projection + buffer zone
- Static MR mesh before dynamic adaptation
- Refinement triggers: wet/dry front, WSE gradient, momentum threshold

### Phase 4 — Output and validation
- NetCDF output: CF-1.8 + UGRID
- MR-aware visualisation

### Phase 5 — Physics and performance
- Full over-relaxed non-orthogonality correction (Weller 2014)
- Riemann solver option
- GPU solver kernel
- Ensemble/calibration: Manning sweep, skill scores

---

## Validation Strategy (next session)

### Option A — Thacker (1981) planar surface
**Best first benchmark.** An exact analytical solution exists for a planar water surface
oscillating in a parabolic bowl. Simpler variant: sloped planar surface (uniform slope,
initially stationary free surface, released). Requires a synthetic sloped DEM (generate
with `generate_slope_dem.py` in `test/`). Tests correct flow direction, approximate
spreading rate, and mass conservation simultaneously.
- Generate 0.1% east-west slope GeoTIFF covering the test AOI
- Set initial depth = 0 on high side, uniform depth on low side
- Run standard model; compare WSE profile across cells to analytical solution
- Success criterion: flow is in the correct direction (eastward), depth profiles within
  ~10% of analytical at t=300s, t=600s, t=900s

### Option B — Symmetric dam-break (circular)
Release a circular water column (high depth at centre, zero elsewhere) on flat terrain.
The Bates (2010) diffusive wave approximation damps the shock front relative to full
Saint-Venant, but the spreading rate and radial symmetry are still testable.
- Set initial volume in source cell (or cluster) rather than using `--rainpoint`
  (requires an `--initial-depth` flag to be added)
- Measure angular variance of the wet-cell front radius — A5 grid should show lower
  directional bias than a rectangular grid at the same resolution
- This is the key scientific validation for the A5 grid paper

### Option C — Mass-conservative equilibrium (simplest)
On a flat domain, inject constant rainfall everywhere (`--rainfall R`) for long enough
that the entire domain reaches a uniform depth. The equilibrium depth should equal
exactly `rate × t / (n_cells × cell_area)`. This is already essentially confirmed by
the point-source tests (mass balance is exact). Extend to uniform rainfall to verify
no preferred-direction bias.

### Option D — Comparison with LISFLOOD-FP or HEC-RAS
For the Kaiapoi domain with real topography: compare inundation extent and timing
against a published or reference model run. This is the ultimate validation but requires
a calibrated event (observed flood data) or published benchmark output.

### Immediate practical step
Add `--initial-depth DEPTH_M` flag: sets all cells to a specified depth at t=0.
This enables dam-break tests without needing a separate mesh preparation step.
Combined with `--rainpoint` turned off and a sloped DEM, this is the cleanest
Thacker-style test setup.

---

## Python Environment

```
Python: C:\Python311\python.exe
Packages: pya5, geopandas, pyarrow, shapely≥2.0, numpy, rasterio (for DEM work)
```

## Julia Packages

```julia
Pkg.add(["PyCall", "JSON3", "DataFrames", "CUDA", "BenchmarkTools",
         "Oxygen", "HTTP", "GLMakie", "HDF5", "Statistics", "ArchGDAL", "Arrow"])
```

## Quick Reference

```bash
# Point-source test (flat terrain, no DEM)
julia --threads auto FloodModel.jl \
    --meshload test/flat_mesh_res14.parquet \
    --flow-model standard \
    --rainpoint -43.4043,172.6644,1000.0 \
    --sim-duration 3600 --dt-max 30 \
    --output test/out.h5 --output-interval 300

# Kaiapoi full run
julia --threads auto FloodModel.jl \
    --meshgen test/kaiapoi_aoi.geojson --meshres 16 \
    --meshout test/kaiapoi_mesh.parquet \
    --dem test/kaiapoi_dem.tif \
    --rainfall 50 --sim-duration 7200 \
    --output test/kaiapoi_out.h5 --vis makie

# Validate test output
julia --project=. test_flat_rainpoint.jl --analyse --file 1000
```
