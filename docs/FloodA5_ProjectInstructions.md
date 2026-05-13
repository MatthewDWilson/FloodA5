# FloodA5 — Project Instructions & Implementation Pathway

_This document is intended as Claude project instructions for development conversations on the FloodA5 codebase. Last reviewed: 2026-05-11 against commit `4005084` (Julia 1.12 update)._

---

## 1. Project Overview

FloodA5 is a 2D surface-water flood model implemented in Julia, built on the **A5 pentagonal Discrete Global Grid System (DGGS)**. The A5 grid tessellates the sphere with equal-area pentagons — each cell has exactly five edge-sharing neighbours, eliminating the directional bias of rectangular grids while retaining simple, uniform topology.

**Key properties of the A5 grid:**
- Five neighbours per cell (uniform connectivity, no directional bias)
- Hierarchical resolution: each cell has a strict parent–child relationship at adjacent levels
- Cell boundaries and centres are geodetic (lon/lat); all geometry uses haversine distances
- No native Julia A5 library — Python (`pya5`) is called via subprocess through `a5_bridge.py`
- Cells at adjacent resolutions belong to two interleaved **sublattices** (e.g. `8e…` and `9e…` hex prefixes at res 14); both must be included in any mesh

**Development environment:**
- Julia 1.12, multi-threaded (`--threads auto`)
- Python 3.11 (`pya5`, `geopandas`, `pyarrow`, `shapely≥2.0`, `numpy`)
- Platform: Windows with OneDrive sync (`F:\OneDrive - University of Canterbury\Julia\FloodA5\`)
- Optional: CUDA.jl for GPU-accelerated point-in-polygon sampling

---

## 2. File Structure

```
FloodA5/
├── FloodModel.jl           # Entry point, CLI, flow model physics, simulation loop
├── A5Grid.jl               # Julia module: mesh generation, cell query API, SGS tables
├── VisualisationServer.jl  # Julia module: HTTP + WebSocket server (Oxygen.jl)
├── MakieVisualiser.jl      # Julia module: GLMakie native desktop viewer
├── FloodViewer.jl          # (auxiliary viewer utilities)
├── a5_bridge.py            # Python bridge: pya5 calls, GeoParquet I/O
├── a5_mesh_diagnostic.py   # Mesh diagnostic utilities
├── setup.jl                # Package setup / environment bootstrap
├── benchmark_pip.jl        # Point-in-polygon benchmark
├── test_a5grid.jl          # A5Grid unit tests
├── test_edge_geometry.jl   # Edge geometry unit tests
├── viz/
│   ├── index.html          # CesiumJS web viewer (single-file, self-contained)
│   ├── config.json         # Runtime config — gitignored, create from example
│   └── config.example.json # Committed template
├── PROJECT_STATE.md        # Cumulative bug log and pending work (kept updated)
├── README.md               # User-facing documentation
└── abstract.md             # Academic abstract (conference/journal submission)
```

---

## 3. Architecture

```
Julia (FloodModel.jl)
  ├── A5Grid.jl
  │     ├── Phase 1 — PIP sampling  [GPU (CUDA) or CPU threads]
  │     └── Phase 2 — subprocess → a5_bridge.py
  │                     ├── pya5: fill_polygon + uncompact (both sublattices)
  │                     ├── grid_disk neighbours + uncompact at target resolution
  │                     ├── coordinate normalisation (NumPy)
  │                     └── geopandas → GeoParquet (EPSG:4326)
  │
  ├── VisualisationServer.jl  (--vis cesium)
  │     Endpoints:
  │       GET /viz/{file}           static files
  │       GET /mesh                 GeoJSON + ordered cell_ids
  │       GET /frames/count         {count, vars}
  │       GET /frames/{idx}         {t, vars} metadata
  │       GET /frames/{idx}/{var}   raw Float32 binary, n_cells × 4 bytes
  │       GET /status               diagnostics
  │       WS  /live                 JSON notifications
  │
  └── MakieVisualiser.jl  (--vis makie)
        GLMakie Figure: poly!, Colorbar, variable Menu, diagnostics sidebar
```

### State variables (primary)

| Variable | Type | Description |
|----------|------|-------------|
| `volume` | `Vector{Float64}` | Stored water volume per cell (m³) — **primary state** |
| `edge_flux` (→ `edges.flux`) | `Vector{Float64}` | Unit discharge q (m²/s) per edge, signed cell_i→cell_j |

Derived at each step: `water_depth`, `velocity` (not yet computed — see §6), `saturation`.

### EdgeList (current undirected edge implementation)

Each edge stored exactly once with `cell_i < cell_j` (canonical lower-index ordering):

```julia
struct EdgeList
    n_edges   :: Int
    cell_i    :: Vector{Int}
    cell_j    :: Vector{Int}
    width     :: Vector{Float64}    # shared edge length (m)
    L         :: Vector{Float64}    # centre-to-centre haversine distance (m)
    cos_theta :: Vector{Float64}    # non-orthogonality correction (1.0 = orthogonal)
    sill      :: Vector{Float64}    # sill elevation (m)
    flux      :: Vector{Float64}    # q (m²/s) at t-dt, signed cell_i → cell_j
end
```

Sign convention: `flux > 0` → flow from `cell_i` to `cell_j`; `flux < 0` → flow from `cell_j` to `cell_i`. Callers apply `dV[ci] += Q*dt; dV[cj] -= Q*dt`.

---

## 4. Physics — Bates et al. (2010) Inertial Formulation

Both `step_standard!` and `step_sgs!` use `_bates_flux`:

```
q^t = [ q^{t-dt} - g · h_flow · dt · dWSE / L_eff ]
      / [ 1 + g · h_flow · dt · n² · |q^{t-dt}| / h_flow^(10/3) ]

Q^t = q^t · width
```

Where:
- `h_flow = max(WSE_i, WSE_j) - z_sill`  (depth above sill on the higher side)
- `dWSE = WSE_i - WSE_j`  (positive when i is higher — flow i→j gives negative q)
- `L_eff = L × cos θ`  (non-orthogonality correction; L is centre-to-centre haversine)
- `z_sill` = `max(elev_i, elev_j)` for standard flow; SGS pre-computed edge minimum for SGS

**Sign convention (critical):** When `WSE_i > WSE_j`, `dWSE > 0`, making `q_new < 0` and `Q < 0`. The caller then applies `dV[ci] += Q*dt` (ci loses volume) and `dV[cj] -= Q*dt` (cj gains volume) — physically correct.

> ⚠️ **Known documentation inconsistency:** The comment in `step_standard!` (line ~861) says "Q > 0 means flow from cell_i to cell_j" which contradicts the `_bates_flux` docstring and the actual numerics. The docstring and numerics are correct; the inline comment is wrong and should be fixed.

**CFL timestep:**
```
dt ≤ CFL × dx² / (2D)    where D ≈ (1/n) × h^(5/3) × √S_ref
```
`cfl = 0.5`, `S_ref ≈ 0.001` as representative slope, `dx = √(cell_area)`.

---

## 5. SGS (Sub-Grid Sampling)

The SGS approach replaces a single mean bed elevation per cell with a **hypsometric curve** — a pre-computed lookup table relating water surface elevation (WSE) to stored volume and wetted plan area.

**Pre-processing (`build_sgs_tables!`):**
1. Distribute 512 Halton quasi-random points within each cell polygon (point-in-polygon test)
2. Sample 1m LiDAR DEM at each point (bilinear interpolation via ArchGDAL)
3. Bin elevations into 100 quantile-spaced levels → cumulative volume and area curves
4. Store as Arrow list columns in the GeoParquet file (`sgs_elev_bins`, `sgs_vol_curve`, `sgs_area_curve`, `sgs_edge_sills`, etc.)

**At runtime:**
- Volume → WSE via inverse interpolation of `vol_curve` (`wse_from_volume`)
- WSE → wetted area via `area_curve` (`wetted_area_from_wse`)
- Sill elevation = minimum DEM elevation along the shared boundary (pre-computed per edge)

**Key benefit:** Partial cell wetting; sub-cell channels and ditches are resolved without finer grid resolution.

---

## 6. Known Issues & Current Hydraulics Bugs

### 6.1 `velocity` is never computed — always zero

`state.velocity` is initialised to `zeros` in `initialise_flow_model` and is never updated in `step_standard!`, `step_sgs!`, or anywhere in the simulation loop. Velocity is displayed in both visualisers and written to HDF5, but will always be zero.

**Fix needed:** After `_apply_dV_*`, compute a per-cell velocity estimate from edge fluxes. A reasonable approach is the flux-weighted average:

```julia
# After edge flux loop, accumulate |Q| × width per cell
# velocity[i] ≈ Σ|Q_e| / (Σwidth_e × depth_i)   (discharge-weighted mean speed)
```

or simply the maximum unit discharge across adjacent edges as a scalar magnitude.

### 6.2 Inline comment contradicts sign convention

In `step_standard!` (~line 861), the comment reads:
> `"Sign convention: Q > 0 means flow from cell_i to cell_j"`

This is **wrong**. When `WSE_i > WSE_j` (i is higher), `dWSE > 0` → `q_new < 0` → `Q < 0`, meaning flow goes from i to j. The `_bates_flux` docstring correctly states `Q > 0 → flow from cell_j to cell_i`. The inline comment should be corrected to avoid future confusion.

### 6.3 Standard flow sill = `max(elev_i, elev_j)` — physically defensible but restrictive

For standard (non-SGS) flow, the edge sill is set to `max(elev_i, elev_j)`. This means `h_flow > 0` only when the WSE on the higher side exceeds the higher bed elevation — a conservative choice that prevents spurious cross-ridge flow. It is consistent with the Bates 2010 literature but means water cannot flow downslope into a lower-elevation cell until it rises above the bed of the sending cell. Consider documenting this explicitly or adding a `--sill-method min|max` option.

### 6.4 Non-orthogonality correction: hard floor at `cos θ = 0.1`

The `_bates_flux` function clamps `cos_theta` to a minimum of 0.1 (`cos(84°)`), preventing `L_eff` from collapsing to zero for highly skewed cell pairs. This is a pragmatic guard, not a physically rigorous treatment. The correct fix (deferred to Phase 5) is the over-relaxed non-orthogonality decomposition of Weller (2014).

### 6.5 CFL estimate uses a fixed representative slope

`_cfl_dt` uses `S_ref = 0.001` (`√S ≈ 0.032`) as a fixed representative slope for all cells. On steep terrain, actual slopes may be much larger, making this estimate non-conservative. A per-cell or per-edge slope from the WSE gradient would be more robust but requires edge data in the CFL loop.

### 6.6 Pending: Projection correction for A5 non-orthogonality (slope bias)

A5 pentagon edges are generally not perpendicular to the centre-to-centre vector. The current correction `L_eff = L × cos θ` projects the distance onto the face normal but does not correct the flux width or the momentum term. A full Weller (2014) over-relaxed decomposition is the proper fix, deferred to Phase 5.

---

## 7. Cumulative Bug History (summary)

All 35 bugs resolved to date are documented in `PROJECT_STATE.md`. The most significant resolved issues affecting hydraulic correctness were:

| # | Issue | Resolution |
|---|-------|-----------|
| 24 | `j==0 && break` silently skipped valid neighbours | Changed to `continue` |
| 26 | `_edge_cos_theta` returned NaN → zeroed all flux | Returns 1.0 (orthogonal fallback) instead |
| 28 | Double-counted edge fluxes | Added `j<=i && continue` (later superseded by EdgeList) |
| 29 | Inverted continuity sign — water flowed uphill | Changed `dV[i] -= vol` to `dV[i] += vol` |
| 31 | `j<=i` guard order-dependent — cells skipped entirely | Replaced with undirected EdgeList (Option D) |
| 32 | Cell ID padding mismatch — zero edges built | Normalised all IDs to 16-char zero-padded hex |
| 34 | `grid_disk` returned wrong resolution (compact) | Added `uncompact(disk, resolution)` |
| 35 | **A5 dual-sublattice** — mesh missing half its cells | `a5_bridge.py` now includes neighbours of primary cells (both sublattices) |

---

## 8. CLI Reference

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

**Defaults:** `--flow-model sgs`, `--sim-duration 3600`, `--dt-max 60`, `--manning-n 0.03`, `--dem-method mean`, `--dem-samples 256`, `--vis-port 8080`, `--output-interval 60`, rainfall `0`.

---

## 9. HDF5 Output Structure

```
/mesh/   cell_ids, elevations, center_lons, center_lats, [static vars]
/frames/ /000001/  t, water_depth, volume, saturation, velocity
         /000002/  ...
```

Datasets: `chunk=(min(n_cells, 4096),)`, `deflate=4`. Read in Python with `h5py`, `xarray`; in Julia with `HDF5.jl`.

---

## 10. Visualisation

### CesiumJS (`--vis cesium`)

Binary wire protocol: one `Float32` array per variable per frame fetched on demand.

**Variable definitions (client-side colormaps):**
```javascript
const VAR_DEFS = {
    depth:      { colormap: 'turbo',  fixedMax: null },
    saturation: { colormap: 'blues',  fixedMax: 1.0  },
    volume:     { colormap: 'turbo',  fixedMax: null },
    velocity:   { colormap: 'plasma', fixedMax: null },
};
```

**Outstanding Cesium issues:**
- Flash on frame updates (in progress)
- `simcomplete` signal not always stopping live updates

### GLMakie (`--vis makie`)

Native desktop window. 1200×760 figure; `poly!` with Observable depth, colorbar, variable Menu dropdown, monospace diagnostics sidebar. No browser or Ion token required.

---

## 11. Python Bridge (`a5_bridge.py`)

**Key pya5 patterns:**
- `fill_polygon + uncompact` returns only one sublattice → must collect `grid_disk(cell, 1)` neighbours and add those within AOI to capture both sublattices
- All cell IDs normalised to 16-char zero-padded hex via `_to_hex` throughout Julia; `u64_to_hex` from pya5 may produce shorter strings
- `grid_disk` returns compact mixed-resolution results → always follow with `uncompact(disk, resolution)`
- Normalise all longitudes: `((lon + 180) % 360) - 180`

---

## 12. Development Phases & Roadmap

### Phase 1 — Core model (current, mostly complete)

- [x] A5 mesh generation (both sublattices, correct adjacency)
- [x] DEM ingestion (ArchGDAL, mean/centroid sampling, Halton points)
- [x] SGS hypsometric tables (pre-computed, stored in parquet)
- [x] Bates (2010) inertial formulation, undirected EdgeList
- [x] Non-orthogonality correction (cos θ factor, with hard floor)
- [x] Adaptive CFL timestep
- [x] Uniform rainfall, point injection sources
- [x] HDF5 output
- [x] CesiumJS + GLMakie visualisation
- [ ] **Velocity computation** (currently always zero — see §6.1)
- [ ] Fix inline sign convention comment (§6.2)
- [ ] Validate against analytical test cases (Thacker, circular dam-break)
- [ ] Ingest LINZ 1m LiDAR for Christchurch domain

### Phase 2 — Boundary conditions

- [ ] `TimeSeries` struct: tiered input (CSV, WaterML 2.0, CF-NetCDF)
- [ ] Upstream hydrograph: time-varying inflow as volume source
- [ ] Open/tidal BC: fixed WSE at boundary cells, ghost-cell enforced
- [ ] Downstream open boundary (free outflow / critical depth)

### Phase 3 — Multi-resolution AMR

- [ ] GPU data layout review (Strategy 1 preferred: solver on GPU, AMR bookkeeping on CPU)
- [ ] `MultiResMesh` data model: multi-level cell store, parent/child A5 IDs
- [ ] Coarse/fine flux interface: ghost-cell projection + buffer zone
- [ ] Static multi-resolution mesh (user-specified refinement regions)
- [ ] Refinement triggers: wet/dry front, WSE gradient, momentum threshold
- [ ] Coarsening criterion: all-5-children stable for N timesteps
- [ ] Dynamic AMR: runtime refine/coarsen with CFL-aware dt update

### Phase 4 — Output and validation

- [ ] NetCDF output: CF-1.8 + UGRID unstructured mesh convention
- [ ] MR-aware visualisation: show resolution level per cell in Cesium and Makie
- [ ] Validation suite: Thacker planar surface, circular dam-break benchmarks
- [ ] Comparison against observed Christchurch flood event

### Phase 5 — Physics and performance upgrades

- [ ] Full over-relaxed non-orthogonality correction (Weller 2014) replacing current cos θ factor
- [ ] LSQ gradient reconstruction for slope term
- [ ] Riemann solver option (swap `_bates_flux` kernel; ghost-cell interface compatible)
- [ ] GPU solver kernel for `step_standard!` / `step_sgs!`
- [ ] Ensemble/calibration: Manning sweep, skill scores (NSE, KGE, F2)

---

## 13. Immediate Next Steps (Priority Order)

1. **Implement velocity computation** in `step_standard!` / `step_sgs!` — edge flux magnitudes → per-cell scalar velocity. This unblocks a meaningful velocity display and HDF5 output.
2. **Fix the inline sign convention comment** in `step_standard!` (~line 861).
3. **Validation run** with Thacker planar-surface test (analytical solution available) to confirm mass conservation and correct flow direction.
4. **LINZ DEM ingestion** for the Christchurch domain at resolution 14 (~2953 cells).
5. **Phase 2 boundary conditions** — upstream hydrograph, open outflow.

---

## 14. Key References

- Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). A simple inertial formulation of the shallow water equations for efficient two-dimensional flood inundation modelling. _Journal of Hydrology_ 387(1–2), 33–45. https://doi.org/10.1016/j.jhydrol.2010.03.027
- Weller, H. (2014). Non-orthogonal version of the arbitrary polygonal C-grid and a new diamond grid. _Geoscientific Model Development_ 7, 779–797.
- A5 DGGS: https://a5geo.org

---

## 15. Suggested Additional Files for Project Context

The following files would improve continuity across separate development conversations if added to the project (e.g. as Claude project knowledge):

| File | Purpose |
|------|---------|
| `PROJECT_STATE.md` | Already exists — paste into new conversations to resume context. Keep updated after each session. |
| `HYDRAULICS.md` | A dedicated technical note on the Bates (2010) implementation, sign convention, cos θ correction, CFL derivation, and the velocity-computation gap. This document would be the reference for any physics-focused conversation. |
| `TESTING.md` | Test plan: Thacker analytical solution, circular dam-break, mass conservation checks, dry-cell stability. Include expected outputs and acceptance criteria. |
| `DATA_FORMATS.md` | Specification of the GeoParquet schema (all `sgs_*` columns), HDF5 layout, and the binary WebSocket wire format. Useful for any visualisation or post-processing conversation. |
| `A5_QUIRKS.md` | Consolidated reference for A5-specific gotchas: dual sublattice, ID padding, `uncompact` requirement after `grid_disk`, longitude normalisation, and the Phase 3 multi-resolution implications. |
| `example_aoi.geojson` | Already exists — include in project so mesh generation examples can be tested without a separate data file. |
