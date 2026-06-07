# FloodA5 — Project State Summary

_Last updated: 2026-06-07 (Bugs 51–52; SGS synthetic DEM validation suite T0–T4 all passing). Paste this into a new conversation to resume._

---

## Project Overview

A Julia flood model using the **A5 pentagonal DGGS**. Cells are pentagons (5 neighbours each for interior cells). No native Julia A5 package — Python (`pya5`) used via subprocess through `a5_bridge.py`.

Two visualisation backends: **Cesium** (CesiumJS web viewer, binary wire protocol) and **Makie** (GLMakie native window). Both opt-in via `--vis [mode]`.

**Location:** `D:\FloodA5\`  (also synced from `F:\OneDrive - University of Canterbury\Julia\FloodA5\`)
**Julia version:** 1.12  **Python:** conda env `flooda5`, Python 3.13 (`pya5`, `geopandas`, `pyarrow`, `shapely≥2.0`, `numpy`)
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
├── test_sgs_unit.jl           # SGS 5-cell unit test (Bug 48 regression)
├── benchmark_sim.jl           # Standalone simulation performance benchmark
├── run_benchmarks.ps1         # PowerShell thread-sweep benchmark runner
├── run_benchmarks.sh          # Bash equivalent
├── viz/
│   ├── index.html             # CesiumJS viewer
│   ├── config.json            # Runtime config — gitignored
│   └── config.example.json   # Committed template
├── test/
│   ├── flat_test_aoi.geojson  # Small AOI for point-source tests (~5.5 km²)
│   ├── flat_mesh_res14.parquet # 61-cell flat mesh (no DEM)
│   ├── kaiapoi_aoi.geojson
│   ├── kaiapoi_dem.tif
│   ├── carlisle/                  # Carlisle test meshes and results
│   └── synthetic_dem/             # Synthetic DEM validation suite
│       ├── generate_synthetic_dem.py
│       ├── test_sgs_synthetic.jl
│       ├── synthetic_dem.tif      # generated, gitignored
│       ├── synthetic_dem_params.json
│       ├── synthetic_aoi.geojson
│       ├── synthetic_mesh_sgs.parquet  # generated, gitignored
│       └── synthetic_mesh_std.parquet  # generated, gitignored
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
    [--rainpoint LAT,LON,RATE_MM_HR]       # repeatable

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
    vel_u       :: Vector{Float64}   # east velocity component (NOT YET COMPUTED)
    vel_v       :: Vector{Float64}   # north velocity component (NOT YET COMPUTED)
    elevation   :: Vector{Float64}   # bed elevation above datum (m)
    manning_n   :: Vector{Float64}   # Manning's roughness per cell
    cell_area   :: Vector{Float64}   # plan area (m²) — recomputed from polygon at init
    cell_lons   :: Vector{Float64}   # cell centre longitude (degrees)
    cell_lats   :: Vector{Float64}   # cell centre latitude (degrees)
    adjacency   :: Dict{String, Vector{String}}
    adj_matrix  :: Matrix{Int}       # (max_nb × n_cells) — retained for neighbour queries
    edges       :: EdgeList          # undirected edge list (primary flux data structure)
    sgs_tables  :: Vector{Any}       # SGSTable per cell (SGSFlow only) or empty
end
```

**Critical:** `velocity`, `vel_u`, `vel_v` are always zero — never computed. See Bug 36.

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
    flux      :: Vector{Float64}   # q (m²/s) at t-dt, signed: +ve = flow j→i
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

## SGS — Sub-Grid Sampling

`wse_from_volume(table, V)` converts stored volume to WSE via inverse interpolation
of the hypsometric volume curve. **As of Bug 52**, when `V > vol_curve[end]` (cell
overfull — water above the terrain ceiling), WSE is extrapolated linearly:

```julia
WSE = z_max + (V - vol_curve[end]) / cell_area
```

This preserves a positive driving head between overfull cells and allows the solver
to propagate volume across flat upstream basins. Prior to this fix, `wse_from_volume`
hard-clamped at `z_max`, stranding excess volume indefinitely.

---

## Source Types

### `--rainfall R` (mm/hr)
Applied uniformly to every cell every timestep. `dV = rate_m_s × dt × cell_area`.

### `--injection-point LAT,LON,RATE` (m³/s)
Fixed volumetric point source at nearest mesh cell. Rate in m³/s.

### `--rainpoint LAT,LON,RATE_MM_HR`
Localised rainfall at nearest mesh cell only. Rate in mm/hr, converted to m³/s
using the individual cell area: `rate_m3s = (mm_hr / 3_600_000) × cell_area`.
Unlike `--rainfall`, does not apply to other cells. Repeatable.
`_find_nearest_cell` warns if the requested point is >2 km from the nearest cell.

---

## Mesh Connectivity Check

`_check_mesh_connectivity(edges, n_cells)` runs after `_build_edge_list` using BFS.
Logs a warning if the mesh has >1 connected component (isolated cells that will never
receive flux). Logs component sizes and marks the largest as the expected source location.

**Key finding from validation:** Small AOIs at res 14 can produce disconnected meshes
where the two A5 sublattices form separate graph components. Use AOIs ≥5 km² at
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
| 37 | **`water_depth` not synced after source injection** | Added pre-physics-step sync: `water_depth[i] = volume[i] / cell_area[i]` |
| 38 | **NaN from `h_flow^(10/3)` underflow** | Added `h_flow = max(h_flow, 1e-6)` floor in `_bates_flux` |
| 39 | **Wrong cell areas from parquet** | Recompute from polygon for non-SGS meshes |
| 40 | **`_find_nearest_cell` silent out-of-mesh** | Added `@warn` when nearest cell is >2 km from requested point |
| 41 | **`launch.json` used `${file}`** | All Julia configs now use `${workspaceFolder}/FloodModel.jl` |
| 42 | **Julia 1.12 rejects `using` inside function body** | Moved to conditional top-level `if DO_ANAL; using HDF5; end` |
| 43 | **Julia 1.12 rejects `Main.varname =` in function** | Replaced with module-level `const _step_debug_count = Ref{Int}(0)` |
| 44 | **Disconnected mesh components** | Added `_check_mesh_connectivity` (BFS); warns with component sizes |
| 45 | **Wrong sign convention comment** | Corrected inline comment in `step_standard!` |
| 46 | **Volume limiter creates mass — asymmetric flux clipping** — The old cell-level post-hoc limiter clipped the donor's loss without reducing the recipient's gain. On a real DEM (e.g. Carlisle, 24.6 m elevation range) step-1 fluxes vastly exceed 50% of cell volume and `domain_vol` grew to ~2× `input_vol` by t=2h. | Moved limiter into flux accumulation loop as per-edge donor cap `V/DONOR_EDGE_DIVISOR`. Same clipped value applied to both donor and recipient — mass conserved exactly. Added constant `N_SIDES = 5`, `DONOR_EDGE_DIVISOR = 10`. |
| 47 | **SGS edge sill lookup fails silently for ordering-mismatched adjacency** — `_build_edge_list` searched only `lo→hi` adjacency slot; ordering mismatch between `grid_disk_neighbours_batch()` and parquet `adj_0..adj_4` caused NaN sill fallback to standard-solver sill, losing SGS benefit. | Sill lookup now tries both `lo→hi` and `hi→lo` perspectives. |
| 48 | **SGS dry-cell drives spurious flux — two-stage fix** — `wse_from_volume(V=0)` returns `z_min`; a high-elevation dry cell (z_min=14.84m) above a low sill (10m) drove flux from the dry cell causing ping-pong oscillation. | Dry cell effective WSE = `max(z_sill, tbl.z_min)`, requiring source WSE to exceed both channel bed AND receiving basin floor. |
| 49 | **SGS oscillations: 100% donor drain + h_flow amplification** — (A) donor cap `V/N_SIDES=V/5` allowed all 5 edges to drain 100% of volume in one step; (B) `z_sill << z_min` made `h_flow >> actual depth`. | (A) Changed cap to `V/DONOR_EDGE_DIVISOR = V/10`; (B) capped `h_flow` at `max(depth_ci, depth_cj)` in `step_sgs!`. |
| 50 | **`_shared_edge` crashes multi-threaded SGS mesh build** — `Set{Tuple}` constructor triggered GC → PyCall finaliser on non-main thread → `EXCEPTION_ACCESS_VIOLATION`. | Replaced `Set` with O(n²) linear scan (max 36 iterations). Use `--threads 1` for mesh generation until confirmed safe. |
| 51 | **Synthetic DEM notch carved too far upstream** — `in_emb` used `\|x − x_emb\| < 3σ`, extending 900 m west of the embankment. All DEM pixels in the upstream notch-latitude band were overwritten to `notch_elev = 0.705 m`, giving those cells `z_min = z_max = 0.705 m`. `wse_from_volume` was permanently clamped at 0.705 m regardless of injected volume → zero WSE gradient → zero flux → T2 deadlock. Also: default notch width was 1,000 m, wider than a res-14 cell (~355 m diameter), defeating the sub-cell purpose of the test. | Changed `in_emb` to `x_emb − σ ≤ x ≤ x_emb + 3σ` (carve only the embankment body itself). Reduced default notch width from 1,000 m to 300 m (= σ). Fix in `test/synthetic_dem/generate_synthetic_dem.py`. |
| 52 | **`wse_from_volume` clamps at `z_max`, stranding excess volume** — `V >= vol_curve[end] && return z_max` hard-clamped WSE regardless of how much additional volume was present. Once all upstream cells hit their terrain ceiling, all returned the same `z_max` → zero gradient → zero inter-cell flux → volume permanently stranded. This is a physics-correctness bug: any scenario where an upstream basin fills above its terrain ceiling produces silent volume trapping and no downstream routing. The synthetic DEM test exposed it because the flat parabolic bowl cells have `z_max < notch_sill`. | Changed to linear extrapolation: `V >= vol_curve[end] && return z_max + (V - vol_curve[end]) / cell_area`. Treats the cell as a vertical-walled container once fully submerged — physically correct. Fix in `A5Grid.jl`. All downstream callers verified safe (water_depth, wetted_area_from_wse, h_flow_cap). |

---

## Validation Results

### Flat-terrain point-source test (standard flow, no DEM)

| Run | Mesh | Rate | Duration | Mass balance | Max wet cells | Notes |
|-----|------|------|----------|-------------|---------------|-------|
| [08] | 29 cells (disconnected) | 50 mm/hr | 1 hr | Perfect | 17/29 (59%) | 12 cells isolated |
| [08b] | 29 cells | 50 mm/hr | 2 hr | Perfect | 17/29 | Same isolation issue |
| [08c] | 61 cells (connected) | 200 mm/hr | 1 hr | Perfect | 17/61 (28%) | Needs ~1.7hr for ring 3 |
| [08d] | 61 cells (connected) | 1000 mm/hr | 1 hr | Perfect | 33/61 (54%) | Ring 4 reached |

Mass balance exact at all checkpoints (`domain_vol = rate × t` to <0.002%).

### Carlisle Validation (2026-06-05)

Standard flow, res 14 (144 cells), single rainpoint 1000 mm/hr, 20h:
- Mass balance: exact (MB error ~7e-9 m³)
- Max depth: 5.7m; 58/144 cells wet at t=20h
- Behaviour: correct ponding in local depression

SGS flow, res 14, same configuration (after Bug 47–49 fixes):
- Mass balance: exact; max depth 6.1m; 58/144 cells wet at t=20h
- Residual oscillations at 1000 mm/hr acceptable (extreme rate for a single res-14 cell)

### SGS Synthetic DEM Validation (2026-06-07) ✓ ALL PASSING

Mesh: 88 cells, res 14, 4 km × 2 km domain. Parabolic bowl upstream, Gaussian embankment
(crest 1.505 m, σ = 300 m), notch (sill 0.705 m, width 300 m). Single rainpoint injection.

| Test | Description | Result |
|------|-------------|--------|
| T0 | Meshes load with cells on both sides of embankment | **PASS** (6/6) |
| T1 | No downstream flow before notch sill (30 min, 50 mm/hr) | **PASS** (1/1) |
| T2 | Downstream flow after notch sill exceeded (sustained injection) | **PASS** (2/2) |
| T3 | SGS routes more water downstream than standard at same injection | **PASS** (1/1) |
| T4 | Mass balance < 0.01% for both solvers over 600 steps | **PASS** (2/2) |

**T3 result (key scientific finding):** At 5,700 steps identical injection history, SGS
routed **902,354 m³ downstream** vs **164,029 m³ for standard** — 5.5× more. Standard
solver's cell-mean elevation at the embankment face is above the notch sill, requiring
upstream WSE to reach 1.089 m before flow. SGS detects the sub-cell channel and routes
at WSE = 0.708 m (just above the 0.705 m notch sill). Mass balance for both solvers:
error < 5×10⁻¹³%.

---

## Pending Work

### Immediate

1. **Bug 36: validate velocity computation** — `state.velocity`, `vel_u`, `vel_v` are
   computed by `_compute_velocity!` (flux-weighted vector sum) but have not been
   validated against a known flow field. Should verify with flat-terrain rainpoint test.

2. **NaN elevation cells in synthetic mesh (cells 36 and 46)** — pre-existing minor issue.
   Cells are on the domain boundary where DEM sampling finds <1 valid pixel. They are
   hydraulically inert (no flux on adjacent edges) and do not affect test results.
   Fix: extend DEM extent slightly beyond AOI, or use `--dem-strict` to flag at mesh-gen.

3. **Carlisle SGS at realistic injection rate** — run Carlisle SGS with 50 mm/hr uniform
   rainfall (not single-cell 1000 mm/hr) to confirm Bug 48/49 oscillation fixes hold at
   realistic rates. The `wse_from_volume` extrapolation fix (Bug 52) should also be
   exercised here.

4. **Bug 46 follow-up: two-pass proportional limiter as user option** — The current
   per-edge donor limit (`V/DONOR_EDGE_DIVISOR`) is conservative. A two-pass
   proportional limiter would scale all outgoing fluxes proportionally when net drain
   would exceed 50% of cell volume, preserving flux priority ordering. Expose as
   `--limiter-mode donor_per_edge|proportional`. Needs performance benchmarking before
   considering as default.

5. **LINZ DEM ingestion** — Kaiapoi domain at res 14/16 with 1m LiDAR.

### Important: SGS mesh generation thread safety (Bug 50)
Use `--threads 1` for mesh generation until the `A5Grid.jl` Bug 50 fix is confirmed
stable. The simulation loop (`step_standard!`, `step_sgs!`) is safe with any thread count.

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
- Thacker (1981) planar surface analytical benchmark (see Validation Strategy below)
- Circular dam-break benchmark (key A5 paper result — directional bias comparison)

### Phase 5 — Physics and performance
- Full over-relaxed non-orthogonality correction (Weller 2014)
- Riemann solver option
- GPU solver kernel
- Ensemble/calibration: Manning sweep, skill scores

---

## Validation Strategy (next steps)

### Option A — Thacker (1981) planar surface
Exact analytical solution for a planar water surface oscillating in a parabolic bowl.
Simpler variant: uniform slope, initially stationary free surface, released.
- Generate 0.1% east-west slope GeoTIFF covering the test AOI
- Run standard model; compare WSE profile to analytical solution
- Success criterion: flow eastward (correct direction), depth profiles within ~10% at t=300/600/900s

### Option B — Symmetric circular dam-break
Release a circular water column on flat terrain. Tests radial symmetry (no directional
bias) — the key scientific result for the A5 grid paper vs rectangular grids.
- Requires `--initial-depth DEPTH_M` flag (sets all cells to specified depth at t=0)
- Measure angular variance of the wet-cell front radius over time

### Option C — Mass-conservative equilibrium
Uniform rainfall on flat domain until uniform depth is reached. Already essentially
confirmed by point-source tests. Extend to verify no preferred-direction bias.

### Immediate practical step
Add `--initial-depth DEPTH_M` CLI flag to enable dam-break tests without separate mesh
preparation. Combined with a sloped DEM this is the cleanest Thacker-style setup.

---

## Python Environment

```
Python: conda env flooda5, Python 3.13
Packages: pya5, geopandas, pyarrow, shapely≥2.0, numpy, rasterio (for DEM work)
```

## Julia Packages

```julia
Pkg.add(["PyCall", "JSON3", "DataFrames", "CUDA", "BenchmarkTools",
         "Oxygen", "HTTP", "GLMakie", "HDF5", "Statistics", "ArchGDAL", "Arrow"])
```

## Quick Reference

```powershell
# Point-source test (flat terrain, no DEM)
julia --threads auto --project=. FloodModel.jl `
    --meshload test/flat_mesh_res14.parquet `
    --flow-model standard `
    --rainpoint -43.4043,172.6644,1000.0 `
    --sim-duration 3600 --dt-max 30 `
    --output test/out.h5 --output-interval 300

# SGS simulation (Carlisle res 14)
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh14_sgs.parquet `
    --rainpoint 54.908,-2.896,50.0 `
    --sim-duration 72000 --flow-model sgs --vis makie

# Generate SGS mesh (--threads 1 required until Bug 50 confirmed fixed)
julia --threads 1 --project=. FloodModel.jl `
    --meshgen test/carlisle/Carlisle_domain.geojson --meshres 14 `
    --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
    --meshout test/carlisle/carlisle_mesh14_sgs.parquet `
    --flow-model sgs --mesh-only

# SGS synthetic DEM validation (all tests)
julia --threads 1 --project=. test/synthetic_dem/test_sgs_synthetic.jl

# SGS unit test (5-cell chain)
julia --threads auto --project=. test_sgs_unit.jl

# Validate HDF5 output
julia --project=. test_flat_rainpoint.jl --analyse --file 1000

# Thread-sweep benchmark
.\run_benchmarks.ps1 -Out test/benchmark_results.csv
```
