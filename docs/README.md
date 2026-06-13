![Status: In Development](https://img.shields.io/badge/Status-In_Development-green)
![Julia 1.12+](https://img.shields.io/badge/Julia-1.12%2B-blue)
![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue)

# FloodA5

A 2D surface-water flood model built on the [A5 pentagonal Discrete Global Grid System
(DGGS)](https://a5geo.org). FloodA5 discretises the landscape into equal-area pentagons,
routes water between them using the Bates et al. (2010) inertial shallow-water
formulation, and resolves sub-cell terrain variability through pre-computed hypsometric
storage tables.

Key characteristics:

- **Uniform five-connectivity** — every interior cell has exactly five edge-sharing
  neighbours, eliminating the directional bias of rectangular grids.
- **Sub-Grid Sampling (SGS)** — hypsometric volume curves built from LiDAR allow
  partial wetting of cells and accurate routing through channels narrower than a cell.
- **Geographic coordinates throughout** — mesh generation, geometry, and output all
  use WGS 84 (EPSG:4326); no projected CRS selection required.
- **LISFLOOD-FP compatible inputs** — reads `.bci`/`.bdy` boundary condition files
  for cross-portability with existing LISFLOOD-FP workflows.
- **Two visualisation backends** — a CesiumJS web viewer and a native GLMakie desktop
  window, both updating live during simulation.

---

## Documentation

| Document | Contents |
|---|---|
| [docs/METHODS.md](docs/METHODS.md) | Physics: A5 grid, SGS, Bates formulation, R-A flux, boundary conditions |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | Full CLI reference, workflow walkthrough, worked examples |
| [docs/DATA_FORMATS.md](docs/DATA_FORMATS.md) | GeoParquet schema, HDF5 layout, `.bci`/`.bdy` format |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | Module map, extension points, developer orientation |
| [docs/TEST_CASES.md](docs/TEST_CASES.md) | Square domain and Carlisle 2005 flood test cases |

---

## Quick Start

### 1. Prerequisites

**Julia 1.12 or later** — download from [julialang.org](https://julialang.org/downloads/).

**Python 3.10 or later** with `pya5` and supporting packages — installed automatically
by the setup script (see below).

Optional: an NVIDIA CUDA GPU accelerates the point-in-polygon DEM sampling step at
high resolutions.

### 2. Install

```bash
git clone https://github.com/MatthewDWilson/FloodA5.git
cd FloodA5
julia setup.jl
```

`setup.jl` activates the Julia project environment, installs all dependencies,
installs Python packages into PyCall's bundled Python, runs a smoke test of the
Python bridge, and checks for GPU availability. See [docs/USER_GUIDE.md](docs/USER_GUIDE.md)
for platform-specific notes.

### 3. Run

```bash
# Generate a mesh, sample a DEM, and build SGS tables (save for reuse)
julia --threads auto FloodModel.jl \
    --meshgen examples/example_aoi.geojson --meshres 14 \
    --dem your_dem.tif --meshout mesh.parquet --mesh-only

# Run a 1-hour simulation with uniform 30 mm/hr rainfall
julia --threads auto FloodModel.jl \
    --meshload mesh.parquet \
    --rainfall 30 --sim-duration 3600 \
    --output sim.h5 --vis makie
```

---

## File Structure

```
FloodA5/
├── FloodModel.jl              # Entry point, CLI, simulation loop
├── setup.jl                   # One-time installation script
├── mesh/
│   ├── A5Grid.jl              # Julia module: mesh generation, cell API, SGS tables
│   ├── a5_bridge.py           # Python bridge: pya5 calls, GeoParquet I/O
│   └── a5_mesh_diagnostic.py  # Mesh diagnostic utilities
├── surfacewater/
│   └── flow2d.jl              # Physics kernels (Bates, Manning R-A, CFL)
├── boundaryinputs/
│   ├── sources.jl             # Water source types (rainfall, injection, hydrograph)
│   ├── boundary_conditions.jl # Ghost-edge BC types and ghost edge geometry
│   └── timeseries_io.jl       # .bdy / CSV hydrograph readers
├── visualisation/
│   ├── MakieVisualiser.jl     # Native desktop viewer (GLMakie)
│   ├── VisualisationServer.jl # HTTP + WebSocket server for CesiumJS viewer
│   ├── view_h5_output.jl      # Post-processing viewer for HDF5 output
│   ├── visualise_mesh.jl      # Publication-quality mesh and SGS figures
│   └── cesium/                # CesiumJS static files
├── test/
│   ├── carlisle/              # Carlisle 2005 flood domain and boundary files
│   ├── square/                # Square domain test case
│   ├── synthetic_dem/         # Synthetic DEM SGS validation suite
│   └── *.jl                   # Unit and integration tests
├── examples/
│   └── example_aoi.geojson    # Example area of interest (centred at 0°N, 0°E)
└── docs/                      # Documentation
```

---

## Resolution Guide

A5 cell size scales by approximately 4× with each resolution level. Choose a level
appropriate to your domain and available DEM resolution.

| Level | Approx. cell area | Cell spacing | Typical use |
|-------|-------------------|--------------|-------------|
| 10 | ~50 km² | ~7 km | Large catchment |
| 12 | ~2 km² | ~1.4 km | Medium catchment |
| 14 | ~13 ha | ~355 m | Urban / detailed |
| 16 | ~8,000 m² | ~89 m | Small-scale urban |
| 18 | ~500 m² | ~22 m | Fine urban / channel |
| 20 | ~31 m² | ~5.6 m | Very high resolution |

SGS provides the most benefit at resolution levels 14–16, where individual cells
typically span heterogeneous terrain (channels, roads, embankments). At resolution 18
and finer, cells are small enough that the standard solver and SGS produce similar
results.

---

## Solvers

**`--flow-model sgs`** (default) — Sub-Grid Sampling. WSE is derived from stored
volume via a pre-computed hypsometric curve. The edge flux uses a Manning hydraulic
radius (R-A) formulation rather than the wide-channel approximation, giving
physically correct friction scaling in confined channels. Requires a DEM.

**`--flow-model standard`** — Diffusive wave on mean cell elevation. Faster and
useful for initial exploration or validation without a DEM.

Both solvers use the Bates et al. (2010) inertial formulation with:
- Wave-speed CFL timestep (Courant number 0.7)
- Froude limiter (Fr ≤ 0.8)
- Q-centred spatial momentum smoothing (θ = 0.9)
- Consistent momentum state (q_prev matches actual transferred flux)

See [docs/METHODS.md](docs/METHODS.md) for full details.

---

## Visualisation

FloodA5 includes two visualisation backends:

**GLMakie** (`--vis makie`) — a native desktop window that updates live during
simulation. No browser, no token, no configuration required. This is the
recommended option.

**CesiumJS** (`--vis cesium` or `--vis`) — a browser-based 3D viewer with
Google Photorealistic Tiles support.

> **Note:** The CesiumJS viewer is currently paused in development and may have
> rough edges. `--vis makie` is recommended for all current use.

By default, domain-edge cells use a **zero-gradient (transmissive) outflow**
condition — water exits freely without reflection. This can be overridden:

- `--closed-boundaries` — no outflow (closed walls; useful for mass-balance testing)
- `--bc-file` — GeoJSON file assigning BC types per segment
- `--inflow-bci` — LISFLOOD-FP `.bci` file for upstream hydrographs and BCs

See [docs/USER_GUIDE.md](docs/USER_GUIDE.md) for full details and
[docs/TEST_CASES.md](docs/TEST_CASES.md) for a worked example using Carlisle `.bci`/`.bdy` files.

---

## Output

Simulation output is written to HDF5 (`.h5`), with per-frame datasets for water
depth, stored volume, saturation fraction, and velocity. Read with Python (`h5py`,
`xarray`) or Julia (`HDF5.jl`).

Post-processing tools:

```bash
# Interactive viewer: map snapshot + cell time series
julia --project=. visualisation/view_h5_output.jl sim.h5 --cell <cell_id>

# Publication figures: mesh wire-frame, SGS cross-section, hypsometric curve
julia --project=. visualisation/visualise_mesh.jl --mesh mesh.parquet <cell_id>
```

See [docs/DATA_FORMATS.md](docs/DATA_FORMATS.md) for the full HDF5 schema.

---

## Testing

```bash
# SGS unit tests (5-cell chain, hypsometric lookups)
julia --threads auto --project=. test/test_sgs_unit.jl

# SGS synthetic DEM validation (T0–T4: mesh, routing, mass balance)
julia --threads 1 --project=. test/synthetic_dem/test_sgs_synthetic.jl

# Edge geometry tests
julia --threads auto --project=. test/test_edge_geometry.jl

# Inflow and boundary condition tests
julia --threads auto --project=. test/test_inflow_point.jl
julia --threads auto --project=. test/test_open_boundary.jl
```

> **Note:** mesh generation must use `--threads 1` due to a thread-safety
> constraint in the Python bridge. Simulation runs are safe with any thread count.

---

## Platform Notes

FloodA5 is developed on Windows. Linux and macOS should work but have received
limited testing — please open an issue if you encounter platform-specific problems.

The `.CondaPkg/` directory in the repository is a Windows-only artefact from
the CondaPkg.jl package manager. It can be ignored on Linux and macOS, where
`setup.jl` installs Python packages via PyCall's bundled Conda directly.

---

## Licence

Apache 2.0 — matching the A5 DGGS library licence.

---

## Citation

If you use FloodA5in research, please cite the accompanying paper (in preparation).
Details will be added here on publication.

---

## References

- Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). A simple inertial formulation
  of the shallow water equations for efficient two-dimensional flood inundation
  modelling. *Journal of Hydrology* 387(1–2), 33–45.
- Neal, J.C. et al. (2012). How much physical complexity is needed to model flood
  inundation? *Hydrological Processes* 26(15), 2264–2282.
- A5 DGGS: [a5geo.org](https://a5geo.org)
