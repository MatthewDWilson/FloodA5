# FloodA5 — Project State Summary

_Last updated: 2026-06-09 (Stage 2 complete: _manning_flux_ra kernel, all unit tests pass, T3 ratio 14.3x). Paste this into a new conversation to resume._

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

**Manning's n per edge (SGS):** `0.5 * (manning_n[ci] + manning_n[cj])` — arithmetic mean (standard, matches LISFLOOD-FP). Changed from `min(n_ci, n_cj)` in `sgs_flow_fixes`.
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
| 53 | **SGS Fix A — Froude limiter missing from `step_sgs!`** — `step_sgs!` called `_bates_flux` with no Froude cap, unlike `step_standard!` which uses `_bates_flux_limited`. On the pentagonal mesh the five independently-evolving `q_prev` values per cell have no directional damping; without a Froude cap, supercritical discharge can develop at SGS edges and drive oscillation. Fix: after the `_bates_flux` call in Phase A, compute `h_flow_eff = max(wse_ci_eff, wse_cj_eff) - z_sill_eff` and clamp `|q|` to `h_flow_eff × √(g × h_flow_eff) × FROUDE_LIMIT (0.8)`. `h_flow_eff` is already bounded by `max(depth_ci, depth_cj)` via the Bug 49 `z_sill_eff` correction, so the cap is physically sound. | Applied inline in `step_sgs!` Phase A after `_bates_flux` call. `FloodModel.jl`. |
| 54 | **SGS Fix C — inconsistent `q_stored` in `step_sgs!`** — The raw unlimited Bates `q` (`Q / width`) was stored to `edges.flux[e]` in Phase A, while Phase B independently capped the transferred volume via `DONOR_EDGE_DIVISOR`. The divergence between stored momentum (`q_prev`) and actual hydraulic transfer is the SGS analogue of the primary driver of standard-flow checkerboarding (Bug primary cause, `standard_flow_fixes`). Fix (primary): write `q_stored` (post-Froude) to `edges.flux[e]` in Phase A. Fix (full): if the Phase B donor cap further clips `ev`, back-propagate the clipped `q` to `edges.flux[e]` so `q_prev` next step always reflects what was transferred. Mirrors LISFLOOD-FP SGC behaviour where `QxSGold` is written within the flux kernel. | Applied in `step_sgs!` Phase A (primary) and Phase B (full). `FloodModel.jl`. |
| 55 | **SGS Manning's n — `min` vs arithmetic mean per edge** — `step_sgs!` used `min(manning_n[ci], manning_n[cj])` per edge. LISFLOOD-FP and standard shallow-water practice use the arithmetic mean `0.5 × (n_ci + n_cj)`. The minimum is marginally non-standard and slightly over-conductive (lower n = less friction = more flux). Not a stability bug, but a systematic discrepancy with the reference implementation. | Changed to `0.5 * (state.manning_n[ci] + state.manning_n[cj])` in `step_sgs!` Phase A. `FloodModel.jl`. |
| 56 | **`sqrt` DomainError in `step_sgs!` Fix A — `h_flow_eff` can be marginally negative** — floating-point rounding when `z_sill_eff` is computed as `max(wse_ci_eff, wse_cj_eff) - h_flow_cap` can yield `z_sill_eff` marginally larger than `max(wse_ci_eff, wse_cj_eff)`, giving `h_flow_eff < 0`. The Froude limiter then calls `sqrt(negative)` → `DomainError`. In this case `_bates_flux` already returned `0.0` so `q_raw = 0` and `q_max = 0` is the correct result. | Added `max(0.0, ...)` floor to `h_flow_eff` computation in Fix A block. `FloodModel.jl`. |
| 57 | **`Threads.@threads` in SGS edge sill loop crashes at res 18 — PyCall GC on non-main thread** — `build_sgs_tables!` used `Threads.@threads` for the edge sill DEM sampling loop (A5Grid.jl:2031). Inside that loop, `_crs_gis` + `ArchGDAL.createcoordtrans` use PyCall. At res 18 (29,902 cells) GC fires frequently on non-main threads, triggering `pydecref` → `PyObject_ClearWeakRefs` → `EXCEPTION_ACCESS_VIOLATION`. Same root cause as Bug 50 (`_shared_edge` thread crash). At res 14 (144 cells) GC fires rarely enough that the crash was not observed. | Changed `Threads.@threads for ci in 1:n` to `for ci in 1:n` (serial). The hypsometric curve build (Step 1) dominates wall time; the edge sill loop is fast serially. `A5Grid.jl`. |
| 58 | **SGS CFL driven to near-zero by extrapolated WSE above `z_max` — oscillation and performance collapse at res 18** — `_apply_dV_sgs!` set `water_depth[i] = wse - z_min` where `wse` is the Bug 52 linear extrapolation above `z_max`. For overfull cells this can produce `water_depth >> z_max - z_min` (e.g. 9m in a cell with 3m terrain range). `_cfl_dt` then uses this as `h_max`, forcing `dt ≈ 1.7s` for 30k cells, causing ~42k steps per simulation hour instead of ~2500. The small `dt` amplifies momentum exchange between the overfull cell and its neighbours, producing the observed domain-wide depth oscillation (`max_depth` alternating ~7.7m / ~9.2m every ~50 steps). Two fixes: (A) cap `water_depth` in `_apply_dV_sgs!` at terrain range `z_max - z_min` — extrapolated head drives WSE gradients correctly but the CFL depth should reflect the physical water column; (B) use the 99th-percentile wet-cell depth in `_cfl_dt` rather than the absolute maximum, so a single extreme cell cannot dominate the timestep for the whole domain. | `_apply_dV_sgs!`: cap `water_depth[i]` at `terrain_depth = max(0, z_max - z_min)`. `_cfl_dt`: collect wet depths, sort, use `wet_depths[round(0.99 × n)]` as `h_cfl`. `FloodModel.jl`. |
| 59 | **SGS checkerboard at res 18 — Froude cap insufficient at small dx/dt; Fix B missing from `step_sgs!`** — At res 18 (`dx ≈ 22m`, `dt ≈ 2.35s`), the Froude cap gives `q_max ≈ 57 m²/s` at h=8m, far larger than the volume-based cap `depth × width / (5×dt) ≈ 7.5 m²/s` used by `step_standard!`. The momentum term drives period-2 `max_depth` oscillation (~7.7m ↔ ~8.9m every 50 steps) at stable amplitude. `step_sgs!` was missing the equivalent of `_bates_flux_limited` Fix B (per-edge volume limiter). The guidance doc recommended omitting it for SGS because `water_depth` was poorly defined, but the Bug 58 terrain-depth cap now makes `state.water_depth` a bounded physical quantity, making Fix B safe to apply. | Added Fix B inline in `step_sgs!` Phase A after Froude clamp. Donor identified by sign of `q_stored` (q>0 → cj donor; q<0 → ci donor). `q_vol_max = depth_donor × width / (5 × dt)`. `q_stored` clamped accordingly before Fix C write-back. `FloodModel.jl`. |
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

### Carlisle SGS res 18 Validation (2026-06-08, Bug 59 Fix B) ✓ OSCILLATION ELIMINATED

Re-run of 1h 50mm/hr after adding Fix B (volume limiter) to `step_sgs!`.

- **`max_depth`:** Locked at 4.528m (terrain-depth cap) from step 250, **never moves**.
  Previously oscillated 7.6m ↔ 9.2m every 50 steps. Oscillation completely gone.
- **`dt`:** Smooth monotonic decline 10s → 4.3s. No erratic fluctuation.
- **`wet`:** 29,853 constant. **`mb_err`:** 0.0 throughout.
- **681 steps / 3,600s** — identical to pre-Fix B run; performance unchanged.
- Fix B does not bind in the early steps (low depth → high q_vol_max); only activates
  as depths build, suppressing the momentum-driven period-2 instability cleanly.

### Carlisle SGS res 18 Validation (2026-06-08, sgs_flow_fixes + Bug 58 fix) ✓ STABLE

SGS flow, res 18 (29,902 cells, 33 NaN boundary cells), uniform 50 mm/hr rainfall, 3,600s (1h),
`--dt-max 60`. Run after Bug 58 fix (CFL terrain-depth cap + 99th-percentile h_max).

- **Mass balance:** `mb_err = 0.0 m³` at every logged checkpoint. Perfect.
- **Stability:** `wet = 29,853` from step 1, constant throughout. No oscillation.
- **`dt` behaviour:** Monotonically declining 60s → 10s (step 1) → 4.3s (step 650) as domain
  fills. No erratic fluctuation. Final step 0.796s is partial step to reach `t = 3600.0s` exactly.
- **`max_depth` cap:** Locks at 4.528m from step 300 — this is `z_max - z_min` of the deepest
  cell (Bug 58 fix, terrain-depth cap). Cell is genuinely overfull; cap is working correctly.
- **681 steps for 3,600s** (avg dt ≈ 5.3s). Dt continues declining as domain fills with no
  outflow — inherent to closed-domain uniform-rainfall scenario. Would stabilise with open
  outflow BC or finite-duration rainfall event.
- **Performance:** ~3–4× slower per step than res 14 due to 200× more cells and SGS lookup.
  Step count ~680/hr vs ~2500/hr at res 14 (dt ~5s vs ~30s). Acceptable for current phase.

### Carlisle SGS Validation (2026-06-08, sgs_flow_fixes) ✓ PASSING

SGS flow, res 14 (144 cells, 2 NaN boundary cells), uniform 50 mm/hr rainfall, 72,000s (20h),
`--dt-max 60`. Run after Fixes 53–56.

- **Mass balance:** `mb_err = 0.0 m³` at every logged checkpoint across all 2,576 steps. Perfect.
- **Stability:** `wet = 143` from step 1, constant throughout. No cells flickering wet/dry.
- **No checkerboarding:** smooth `vol_sum` progression (~65,000 m³/10 steps); pre-fix SGS runs
  showed visible ping-pong in the Makie depth field. Absent here.
- **`dt` behaviour:** 60s → ~25–30s as depths build; settles 20–27s for second half.
  Occasional dips to ~20s correspond to transient CFL tightening at local depth spikes — correct.
- **`max_depth`:** 2.5m at t≈2,757s → peak ~13.1m at t≈64,095s; settles ~9–11m as
  water redistributes over the domain. Coherent closed-domain ponding behaviour.
- **SGS edge sills:** 2/318 edges NaN (no DEM data along those edges, same as pre-fix).

### Carlisle Validation (2026-06-05)

Standard flow, res 14 (144 cells), single rainpoint 1000 mm/hr, 20h:
- Mass balance: exact (MB error ~7e-9 m³)
- Max depth: 5.7m; 58/144 cells wet at t=20h
- Behaviour: correct ponding in local depression

SGS flow, res 14, same configuration (after Bug 47–49 fixes):
- Mass balance: exact; max depth 6.1m; 58/144 cells wet at t=20h
- Residual oscillations at 1000 mm/hr acceptable (extreme rate for a single res-14 cell)

### SGS Synthetic DEM Validation (2026-06-08, sgs_flow_fixes) ✓ ALL PASSING — regression confirmed

Re-run after Fixes 53–56 (Froude limiter, consistent q_stored, Manning's n mean, sqrt floor).
All results match the pre-fix baseline to within machine epsilon — Froude cap not binding
on this subcritical synthetic DEM; Fix C changes only stored momentum, not transferred volume.

| Test | Result | Key value |
|------|--------|-----------|
| T0 | **PASS** (6/6) | — |
| T1 | **PASS** (1/1) | Max upstream WSE = 0.512 m < notch sill 0.705 m |
| T2 | **PASS** (2/2) | Downstream = 902,353.7 m³ (was 902,353.9 m³) |
| T3 | **PASS** (1/1) | SGS 902,353.7 m³ vs Standard 164,028.5 m³ (5.5×) |
| T4 | **PASS** (2/2) | SGS error 4.84×10⁻¹³%; Standard 3.46×10⁻¹³% |

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

### sgs_ra_flux branch — R-A flux formulation for SGS

See `FloodA5_SGS_RA_Flux_Implementation.md` for full implementation plan.

**Stage 1 — Pre-processing (✅ complete, 2026-06-09):**
- ✅ `SGSTable` struct extended with `edge_area_curves`, `edge_perim_curves`
- ✅ `build_sgs_tables!` Step 2 extended to build edge flow area / wetted perimeter curves
- ✅ New parquet columns: `sgs_edge_area_curve`, `sgs_edge_perim_curve`
- ✅ New lookup functions: `flow_area_from_wse`, `wetted_perim_from_wse`, `hydraulic_radius_from_wse`
- ✅ `sgs_table()` updated with backward-compat fallback for legacy meshes
- ✅ Unit test T-EA1–T-EA6: 251 assertions, all pass
- ✅ `test_sgs_unit.jl` regression: all 49 tests pass
- ✅ Synthetic mesh regenerated with new parquet columns (3.5s build)
- ✅ T0–T4 regression: all pass, results identical to pre-Stage-1 baseline

**Stage 2 — Flux kernel (✅ complete, 2026-06-09):**
- ✅ `EdgeList.flux_Q` (m³/s, SGS R-A only) — new field, zeros for standard flow
- ✅ `_adjacency_slot` helper
- ✅ `_manning_flux_ra` kernel — LISFLOOD-FP SGC formulation with R^(4/3)·A denominator
- ✅ `step_sgs!` Phase A: replaced Bates+limiters with `_manning_flux_ra`
- ✅ Phase B: donor cap back-propagates to `flux_Q` (Fix C full)
- ✅ `_write_frame!`: writes `flux_Q` dataset when non-zero
- ✅ `test_sgs_unit.jl`: all 49 tests pass; water now reaches cell 5 (58m³ vs 0m³ with Bates+limiter)
- ✅ T0–T4 all pass; T3: SGS 360,419 m³ vs standard 25,153 m³ (14.3× ratio, up from 5.5×)
- ✅ T4 mass balance: 1.07×10⁻¹²% (machine epsilon)
- ✅ Carlisle res 14 legacy mesh: backward-compat fallback confirmed (no crash, no routing)
- ⏳ **Regenerate Carlisle res 14 SGS mesh** — must include new parquet columns
- ⏳ **Regenerate Carlisle res 18 SGS mesh** — must include new parquet columns
- ⏳ **Validate Carlisle res 14 50mm/hr** — confirm no oscillation, mb_err=0
- ⏳ **Validate Carlisle res 18 1hr** — key test: dt stable, no period-2 oscillation

**Stage 3 — Cleanup (⏳ after validation):**
- New `EdgeList.flux_Q` array (m³/s, SGS R-A only)
- New `_manning_flux_ra` function
- `_adjacency_slot` helper
- `step_sgs!` Phase A: replace `_bates_flux` + Bates limiters with `_manning_flux_ra`
- Unit test: `test/test_manning_ra_flux.jl` (T-RA1 through T-RA6)
- Regenerate Carlisle res 14 + res 18 meshes and validate

**Stage 3 — Cleanup:**
- Remove `sgs_flow_fixes` Bug 53/59 limiters from `step_sgs!` (superseded by R-A)
- Update `DATA_FORMATS.md`, `HYDRAULICS.md`
- Full validation suite across res 14 and res 18

### Important: mesh generation thread safety (Bugs 50, 57)
Use `--threads 1` for ALL mesh generation. Simulation loop is safe with any thread count.

### Other pending
- **Bug 36: velocity computation validation** — `_compute_velocity!` exists but untested.
- **LINZ DEM ingestion** — Kaiapoi domain res 14/16.
- **Open outflow BC** (Phase 2) — prevents dt declining on closed-domain rainfall tests.

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
