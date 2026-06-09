![Static Badge](https://img.shields.io/badge/Status:-In_development-red)


# A5 Pentagon Flood Model

A flood modelling application in Julia using the [A5 DGGS](https://a5geo.org)
— a pentagonal Discrete Global Grid System offering equal-area cells with
minimal shape distortion. Two visualisation backends are available: a
**CesiumJS** web viewer and a native **GLMakie** desktop window.

---

## File Structure

```
FloodA5/
├── FloodModel.jl          ← Entry point, flow model logic, CLI
├── mesh/
│   ├── A5Grid.jl              ← Julia module: mesh generation & cell query API
│   ├── a5_bridge.py           ← Python bridge: pya5 calls, GeoParquet I/O
│   └── a5_mesh_diagnostic.py
├── surfacewater/
│   └── flow2d.jl              ← Pure physics kernels (Bates, Manning R-A, CFL)
├── visualisation/
│   ├── MakieVisualiser.jl     ← Julia module: GLMakie native desktop viewer
│   ├── VisualisationServer.jl ← Julia module: HTTP + WebSocket server (Oxygen.jl)
│   ├── FloodViewer.jl         ← Standalone post-processing viewer
│   └── cesium/                ← CesiumJS web viewer static files
├── test/                      ← All tests and benchmarks
│   ├── index.html         ← CesiumJS web viewer
│   ├── config.json        ← Runtime config (gitignored — create from example)
│   └── config.example.json ← Committed config template
└── data/
    ├── christchurch_aoi.geojson
examples/
    └── example_aoi.geojson
```

---

## Architecture

```
Julia (FloodModel.jl)
  ├── mesh/A5Grid.jl
  │     ├── Phase 1 — PIP sampling  [Julia: GPU (CUDA) or CPU threads]
  │     └── Phase 2 — subprocess → mesh/a5_bridge.py
  │                     ├── pya5: fill_polygon + uncompact
  │                     ├── cell_to_boundary × N cells
  │                     ├── coordinate normalisation (NumPy)
  │                     └── geopandas → GeoParquet (EPSG:4326)
  │
  ├── visualisation/VisualisationServer.jl  (--vis cesium)
  │     ├── HTTP  GET  /viz/{file}           ← static files (index.html, config.json)
  │     ├── HTTP  GET  /mesh                 ← GeoJSON + ordered cell_ids
  │     ├── HTTP  GET  /frames/count         ← total frame count
  │     ├── HTTP  GET  /frames/{idx}         ← frame metadata (t, var names)
  │     ├── HTTP  GET  /frames/{idx}/{var}   ← binary Float32 array (per-variable)
  │     ├── HTTP  GET  /status               ← server diagnostics
  │     └── WS    /live                      ← live notifications
  │
  └── visualisation/MakieVisualiser.jl  (--vis makie)
        └── GLMakie window: polygons, colorbar, variable selector, diagnostics sidebar
```

### GPU / CPU selection

| Condition | Backend |
|-----------|---------|
| CUDA.jl available + points ≥ 500,000 | GPU (CUDA kernel) |
| Otherwise | CPU (multi-threaded) |

---

## Prerequisites

### Python

```bash
pip install pya5 geopandas pyarrow shapely numpy
```

Verify: `python mesh/a5_bridge.py check`

### Julia packages

```julia
using Pkg
Pkg.add(["PyCall", "JSON3", "DataFrames", "CUDA", "BenchmarkTools",
         "Oxygen", "HTTP", "GLMakie", "HDF5"])
```

`CUDA`, `BenchmarkTools`, and `GLMakie` are optional.

### Cesium Ion token (for `--vis cesium` with Google 3D tiles)

```bash
cp visualisation/cesium/config.example.json visualisation/cesium/config.json
# Edit visualisation/cesium/config.json and set your token from ion.cesium.com/tokens
```

Add `visualisation/cesium/config.json` to `.gitignore`. The viewer validates the token against
the Ion API at startup and falls back gracefully if it is missing.

---

## Usage

### Command line

All arguments are named flags — position-independent.

```
julia [--threads auto] FloodModel.jl  --meshgen <aoi.geojson>  --meshres <N>
                                      [--meshout <file>]
                                      [--dem <file.tif>]  [--dem-strict]
                                      [--flow-model sgs|standard]
                                      [--rainfall R]  [--sim-duration S]
                                      [--output <file.h5>]
                                      [--vis [mode]]  [--vis-port PORT]

julia [--threads auto] FloodModel.jl  --meshload <mesh.parquet>
                                      [same options as above]

julia FloodModel.jl --help | -h
```

#### Mesh flags

| Flag | Description |
|------|-------------|
| `--meshgen FILE` | GeoJSON AOI → triggers mesh generation |
| `--meshres N` | A5 resolution level (required with `--meshgen`) |
| `--meshout FILE` | Save mesh (`.parquet` or `.geojson`); omit to keep in memory |
| `--meshload FILE` | Load saved GeoParquet mesh, skip generation |

#### DEM flags

| Flag | Description |
|------|-------------|
| `--dem FILE` | GeoTIFF elevation; baked into mesh parquet for reuse |
| `--dem-strict` | Error on out-of-bounds cells (default: NaN + warning) |
| `--dem-method mean\|centroid` | Sampling method (default: `mean`) |
| `--dem-samples N` | Halton points for mean sampling (default: 256) |

#### Flow model flags

| Flag | Default | Description |
|------|---------|-------------|
| `--flow-model sgs\|standard` | `sgs` | SGS = sub-grid hypsometric; standard = mean elevation |
| `--manning-n N` | `0.03` | Global Manning's roughness |
| `--friction FILE` | — | GeoTIFF friction raster (overrides `--manning-n`) |
| `--rainfall R` | `0` | Uniform rainfall rate (mm/hr) |
| `--sim-duration S` | `3600` | Simulation duration (seconds) |
| `--dt-max S` | `60` | Maximum adaptive timestep (seconds) |

#### Output flags

| Flag | Description |
|------|-------------|
| `--output FILE` | Write HDF5 time-series output (`.h5` / `.hdf5`) |
| `--output-interval S` | Seconds between snapshots (default: 60) |

#### Visualisation flags

| Flag | Description |
|------|-------------|
| `--vis [MODE]` | Enable visualisation (off by default). Omitting `MODE` → `cesium` |
| `--vis-port PORT` | Cesium HTTP port (default: 8080; ignored for `--vis makie`) |

#### Examples

```bash
# Generate mesh, sample DEM, SGS flow, 1hr simulation, Cesium viewer
julia --threads auto FloodModel.jl \
    --meshgen christchurch_aoi.geojson --meshres 14 --meshout mesh_sgs.parquet \
    --dem linz_dem.tif --rainfall 30 --sim-duration 3600 --output sim.h5 --vis

# Load saved mesh, run with Makie viewer
julia --threads auto FloodModel.jl \
    --meshload mesh_sgs.parquet --rainfall 10 --sim-duration 7200 --vis makie

# Quick test — standard flow, no DEM
julia --threads auto FloodModel.jl \
    --meshload mesh_sgs.parquet --flow-model standard --sim-duration 1800

# No visualisation (fastest)
julia --threads auto FloodModel.jl --meshload mesh_sgs.parquet
```

---

## Visualisation

### CesiumJS viewer (`--vis cesium`)

Open `http://localhost:8080` in a browser while the model runs.

**Variable selector** — dropdown in the left panel switches between:
- Depth (m) · Saturation (0–1) · Volume (m³) · Velocity (m/s)

Each variable fetches binary `Float32` data from `GET /frames/{idx}/{varname}`
on demand — only one variable is fetched at a time, scaling to 1M+ cells.

**Controls:**
- Colorbar with auto-scale toggle and manual ceiling slider
- 3D extrusion with adjustable vertical exaggeration
- Cell outlines and zero-value cell visibility toggles
- Basemap selector: **Google 3D Photorealistic Tiles** / **OSM** / **None**
- Live mode (follows simulation) / Replay mode (scrub timeline)
- Live button disabled and status set to "complete" when simulation ends

**Wire protocol:**

| Endpoint / message | Description |
|--------------------|-------------|
| WS `mesh` | GeoJSON + `cell_order` array (sent on connect) |
| WS `framecount` | Historical frame count + variable list |
| WS `newframe` | `{idx, t, vars}` — lightweight notification only |
| WS `simcomplete` | `{frames}` — simulation finished |
| `GET /frames/{idx}/{var}` | Raw `Float32` binary, n_cells × 4 bytes |

### GLMakie viewer (`--vis makie`)

Native desktop window — no browser or Ion token required.

- Pentagon polygons coloured by selected variable
- Dropdown selector: Depth · Saturation · Volume · Velocity
- Colorbar with auto-scaling
- Sim-time axis title (live `@lift`)
- Monospace diagnostics sidebar: frame, wet cells, max value, mesh info, clock

---

## HDF5 Output

```
/mesh/
    cell_ids, elevations, center_lons, center_lats, [other static vars]
/frames/
    /000001/  t, water_depth, volume, saturation, velocity
    /000002/  ...
```

Datasets are chunked (`min(n_cells, 4096)`) and gzip-compressed (level 4).
Read in Python with `xarray`, `h5py`, or `pandas`; in Julia with `HDF5.jl`.

---

## Resolution Guide

| Level | Approx. cell area | Typical use case |
|-------|-------------------|------------------|
| 8 | ~250 km² | Large catchment |
| 10 | ~50 km² | Medium catchment |
| 12 | ~10 km² | Small catchment |
| 14 | ~2 km² | Urban / detailed |
| 17 | ~0.1 km² | High-resolution |

---

## Flow Model

### Solver options

**`--flow-model sgs`** (default) — diffusive wave with sub-grid hypsometric
storage tables. WSE is derived from stored volume via a hypsometric curve
built from DEM samples. Allows partial wetting of cells; accurately captures
channels and ditches narrower than a cell.

**`--flow-model standard`** — diffusive wave on mean cell elevation. Faster,
suitable for initial runs and validation without a DEM.

Both solvers use the **Bates et al. (2010) inertial formulation**:

```
q^t = [ q^{t-dt} - g·h_flow·dt·(dWSE/L) ]
      / [ 1 + g·h_flow·dt·n²·|q^{t-dt}| / h_flow^(10/3) ]
```

where `L` is the **centre-to-centre haversine distance** between adjacent cells
(stored in `FlowState.centre_dist`), and `q^{t-dt}` is the per-edge unit
discharge from the previous timestep (stored in `FlowState.edge_flux`).

### Flow model pending work

1. Projection correction for non-orthogonal A5 edges (slope bias)
2. DEM ingestion: LINZ 1m LiDAR for Christchurch
3. Boundary conditions: point sources, upstream inflow
4. Output: NetCDF in addition to HDF5

---

## Python Bridge

```bash
python mesh/a5_bridge.py check                                    # verify environment
python mesh/a5_bridge.py mesh_for_aoi aoi.geojson 14 out.parquet geoparquet
```

**Key pya5 quirks:**
- Cell IDs are Python `int` — use `u64_to_hex` for serialisation
- Normalise all longitudes: `((lon + 180) % 360) - 180`
- `grid_disk` → always follow with `uncompact()` for uniform meshes

---

## Performance

Resolution 14, Kaiapoi / Christchurch AOI (~200 cells in test runs):

| Stage | Time |
|-------|------|
| Mesh generation (res 14, ~3K cells) | ~75s |
| `--meshload` (skip generation) | <1s |
| SGS initialisation | ~5s |
| Simulation step | depends on dt, n_cells |

Use `--meshload` on repeat runs to eliminate the ~75s generation cost.

---

## Licence

Apache 2.0 (matching the A5 library licence)
