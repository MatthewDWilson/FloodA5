# FloodA5 — User Guide

This guide covers installation, the standard modelling workflow, a complete CLI
reference, and post-processing options. For the underlying physics see
[METHODS.md](METHODS.md); for file format details see [DATA_FORMATS.md](DATA_FORMATS.md).

---

## 1. Installation

### 1.1 Prerequisites

**Julia 1.12 or later.** Download from [julialang.org](https://julialang.org/downloads/).
Install the binary for your platform; no compilation is required.

**Python 3.10 or later.** The Python packages required by the A5 library are
installed automatically by the setup script (see §1.2). If you manage Python
environments manually, the required packages are:

```
pya5  geopandas  pyarrow  shapely  numpy
```

**CUDA (optional).** An NVIDIA CUDA GPU accelerates the point-in-polygon DEM
sampling step. FloodA5 falls back to multi-threaded CPU automatically if no
functional GPU is detected.

### 1.2 Setup

Clone the repository and run the setup script from the project root:

```bash
git clone https://github.com/<your-org>/FloodA5.git
cd FloodA5
julia setup.jl
```

The setup script:
1. Activates the Julia project environment and installs all Julia packages.
2. Installs Python packages into PyCall's bundled Conda Python.
3. Runs a smoke test of the Python bridge (`mesh/a5_bridge.py`).
4. Checks for a functional CUDA GPU and reports whether GPU acceleration is available.

If any step fails, the script reports the error and continues. Failures in the
optional GPU step do not prevent the model from running.

### 1.3 Platform notes

FloodA5 is developed and tested on Windows. Linux and macOS are expected to work
but have received limited testing.

**Linux / macOS:** The `.CondaPkg/` directory in the repository is a Windows artefact
from the CondaPkg.jl package manager. It can be safely ignored. The setup script
uses PyCall's bundled Conda on all platforms, so no separate Conda installation is
needed.

**Linux:** If `pya5` fails to install via pip (some distributions require a Rust
toolchain), try:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# then re-run julia setup.jl
```

**macOS Apple Silicon:** CUDA is not available. All PIP sampling runs on CPU threads,
which is slower but fully functional.

### 1.4 Verify the installation

```bash
julia --threads auto FloodModel.jl --help
```

This should print the full CLI reference. If it errors, check that `setup.jl`
completed without fatal errors and that `python mesh/a5_bridge.py check` prints
`pya5 OK`.

---

## 2. Workflow Overview

A typical FloodA5 run follows three stages:

```
1. MESH GENERATION    AOI GeoJSON + resolution  →  .parquet mesh file
2. DEM SAMPLING       GeoTIFF DEM               →  elevation + SGS tables in mesh
3. SIMULATION         sources + BCs + duration  →  .h5 output + live viewer
```

Stages 1 and 2 are the most time-consuming and produce a reusable `.parquet` file.
Stage 3 can be re-run with different sources, roughness, or duration without
regenerating the mesh.

```bash
# Stage 1+2: generate mesh, sample DEM, build SGS tables (save to disk)
julia --threads auto FloodModel.jl \
    --meshgen my_aoi.geojson --meshres 14 \
    --dem my_dem.tif \
    --meshout my_mesh.parquet --mesh-only

# Stage 3: run simulation, reload mesh from disk
julia --threads auto FloodModel.jl \
    --meshload my_mesh.parquet \
    --rainfall 30 --sim-duration 7200 \
    --output my_sim.h5 --vis
```

> **Thread count:** always use `--threads 1` for mesh generation (stage 1+2),
> due to a thread-safety constraint in the Python bridge. Simulation runs
> (stage 3) are safe with `--threads auto`.

---

## 3. Preparing Your Inputs

### 3.1 Area of interest (AOI)

The AOI must be a GeoJSON file containing a single polygon in EPSG:4326
(decimal degrees, longitude before latitude). Multi-polygon AOIs are not
currently supported; use the outer ring only.

```json
{
  "type": "Feature",
  "properties": {},
  "geometry": {
    "type": "Polygon",
    "coordinates": [[
      [172.55, -43.60],
      [172.75, -43.60],
      [172.75, -43.45],
      [172.55, -43.45],
      [172.55, -43.60]
    ]]
  }
}
```

The polygon should be closed (first and last coordinate pair identical). Coordinates
are `[longitude, latitude]` — GeoJSON convention.

An example AOI centred at 0°N, 0°E is provided at `examples/example_aoi.geojson`.

### 3.2 DEM

FloodA5 accepts any GDAL-readable raster format (GeoTIFF recommended). The DEM can
be in any projected or geographic CRS — FloodA5 reprojects sample coordinates into
the DEM's native CRS automatically using ArchGDAL.

For SGS, a DEM with a resolution of 1–5 m gives the best results at resolution
levels 14–18. Coarser DEMs (10–30 m) are adequate for the standard solver or at
resolution levels 10–12.

Cells whose centres fall outside the DEM extent receive `NaN` elevation. These cells
are hydraulically inert (no flux). Use `--dem-strict` to treat out-of-bounds cells
as an error if this is not acceptable.

### 3.3 Friction raster

An optional per-cell Manning's roughness raster can be supplied via `--friction`.
The raster must be GDAL-readable. Cell values are sampled at cell centres. Cells
where the raster is `NaN` or outside the raster extent use the global `--manning-n`
value (default 0.03 s/m^(1/3)).

### 3.4 Hydrograph files

Time-varying inflows are specified as two-column CSV files (`t_s, Q_m3s`) or
LISFLOOD-FP `.bdy` files. See §5.4 and [DATA_FORMATS.md](DATA_FORMATS.md) for format details.

---

## 4. CLI Reference

All arguments are named flags; position does not matter.

### 4.1 Mesh options

| Flag | Description |
|---|---|
| `--meshgen FILE` | GeoJSON AOI file. Triggers mesh generation via the Python bridge. |
| `--meshres N` | A5 resolution level. Required with `--meshgen`. |
| `--meshout FILE` | Save generated mesh to `.parquet` (recommended) or `.geojson`. Omit to keep in memory only. |
| `--meshload FILE` | Load a previously saved GeoParquet mesh, skipping generation. |
| `--mesh-only` | Generate (or load) and save the mesh, then exit without simulating. Implied when no water source is provided. |

### 4.2 DEM options

| Flag | Default | Description |
|---|---|---|
| `--dem FILE` | — | GeoTIFF elevation file. Samples DEM and bakes elevation into the mesh parquet. On `--meshload`, errors if elevation is already in the parquet to prevent silent overwrites. |
| `--dem-strict` | off | Error if any sample point falls outside the DEM extent. Default: assign NaN and continue with a warning. |
| `--dem-method mean\|centroid` | `mean` | Sampling method. `mean`: arithmetic mean of Halton-distributed points (better physics). `centroid`: single sample at cell centre (fast). |
| `--dem-samples N` | 256 | Halton candidate points per cell for `mean` sampling. |
| `--dem-seed N` | 0 | Halton sequence offset. Use different values for independent Monte Carlo uncertainty estimates. |

### 4.3 Flow model options

| Flag | Default | Description |
|---|---|---|
| `--flow-model sgs\|standard` | `sgs` | Flow routing method. `sgs`: sub-grid hypsometric (R-A kernel). `standard`: diffusive wave on mean cell elevation. |
| `--sgs-bins N` | 100 | Elevation bins for SGS hypsometric curves. |
| `--sgs-samples N` | 512 | Halton points per cell for SGS cross-section pre-processing. |
| `--manning-n N` | 0.03 | Global Manning's roughness coefficient (s/m^(1/3)). |
| `--friction FILE` | — | GeoTIFF friction raster. Per-cell Manning's n sampled at cell centres; overrides `--manning-n` where the raster is finite. |
| `--sim-duration S` | 3600 | Simulation duration in seconds. |
| `--dt-max S` | 60 | Maximum adaptive timestep in seconds. Reduce to 10–30 s for high-resolution (level 18+) or steep-terrain runs. |

### 4.4 Water source options

All source flags are repeatable. Multiple sources are summed; sources targeting
the same cell are summed before flux routing.

| Flag | Description |
|---|---|
| `--rainfall R` | Uniform rainfall rate in mm/hr applied to every cell. |
| `--rainpoint LAT,LON,RATE` | Localised rainfall at the nearest cell (mm/hr). Unlike `--rainfall`, applies to one cell only. Useful for point-source validation. |
| `--injection-point LAT,LON,RATE` | Constant volumetric injection at the nearest cell. RATE in m³/s. |
| `--inflow-point LAT,LON,FILE[,LABEL]` | Time-varying hydrograph at the nearest cell. FILE is a two-column CSV (`t_s, Q_m3s`) or a LISFLOOD-FP `.bdy` file. LABEL selects a named series from a multi-series `.bdy` file; defaults to the filename stem. |
| `--inflow-bci FILE` | Load inflow configuration from a LISFLOOD-FP `.bci` file. `QVAR` entries reference a companion `.bdy` time-series file (same directory, same stem as the `.bci` by default). `QFIX` entries create constant injection points. `FREE` entries set open outflow at that boundary cell. |
| `--inflow-bdy FILE` | Explicit path to the `.bdy` file for `QVAR` entries in `--inflow-bci`. Use when the `.bdy` and `.bci` files have different names or are in different directories. |

### 4.5 Boundary condition options

| Flag | Default | Description |
|---|---|---|
| `--bc-file FILE` | — | GeoJSON file assigning BC types to boundary segments. Each feature must have a `"bc_type"` property: `Closed`, `ZeroGradient`, or `Critical`. Boundary cells within 1.5× the cell diameter of each feature are assigned the specified type. |
| `--closed-boundaries` | off | Treat all domain-edge cells as closed walls (no outflow). Overrides `--bc-file`. Useful for mass-balance benchmarking. |
| `--bc-epsg CODE` | WGS 84 | Treat coordinates in `--inflow-bci` files as projected (easting/northing) in the given EPSG code, converting to longitude/latitude internally. Example: `--bc-epsg 27700` for British National Grid. |

### 4.6 Output options

| Flag | Default | Description |
|---|---|---|
| `--output FILE` | off | Write HDF5 simulation output (`.h5` / `.hdf5`). |
| `--output-interval S` | 60 | Seconds of simulation time between output snapshots. |

### 4.7 Visualisation options

| Flag | Default | Description |
|---|---|---|
| `--vis [MODE]` | off | Enable live visualisation. Omitting MODE defaults to `cesium`. Available modes: `cesium`, `makie`. |
| `--vis-port PORT` | 8080 | HTTP port for the Cesium server. Ignored for `--vis makie`. |

---

## 5. Worked Examples

### 5.1 Minimal run — standard solver, no DEM

The fastest way to check that the installation works:

```bash
# Generate a small flat mesh (no DEM — all elevations zero)
julia --threads 1 FloodModel.jl \
    --meshgen examples/example_aoi.geojson --meshres 14 \
    --meshout test_mesh.parquet --flow-model standard --mesh-only

# Inject water at the mesh centre and watch it spread
julia --threads auto FloodModel.jl \
    --meshload test_mesh.parquet \
    --flow-model standard \
    --rainpoint 0.0,0.0,50.0 \
    --sim-duration 3600 \
    --vis makie
```

### 5.2 SGS run with a DEM

```bash
# Stage 1+2: mesh generation with DEM and SGS tables (run once; takes a few minutes)
julia --threads 1 FloodModel.jl \
    --meshgen my_aoi.geojson --meshres 14 \
    --dem my_dem.tif \
    --meshout my_mesh_sgs.parquet --mesh-only

# Stage 3: simulation (fast to re-run with different parameters)
julia --threads auto FloodModel.jl \
    --meshload my_mesh_sgs.parquet \
    --rainfall 30 \
    --sim-duration 7200 \
    --output my_sim.h5 --output-interval 300 \
    --vis
```

### 5.3 Upstream hydrograph via LISFLOOD-FP files

```bash
julia --threads auto FloodModel.jl \
    --meshload my_mesh_sgs.parquet \
    --inflow-bci my_domain.bci \
    --sim-duration 432000 \
    --dt-max 30 \
    --output event_sgs.h5 --output-interval 3600
```

If the `.bci` file uses projected coordinates (e.g. British National Grid):

```bash
    --inflow-bci my_domain.bci --bc-epsg 27700
```

### 5.4 Two-column CSV hydrograph

```bash
julia --threads auto FloodModel.jl \
    --meshload my_mesh_sgs.parquet \
    --inflow-point -43.386,172.648,river_inflow.csv,Upstream \
    --sim-duration 86400 \
    --output daily_sim.h5
```

`river_inflow.csv` format (header optional, auto-detected):

```
t_s,Q_m3s
0,0.0
3600,5.2
7200,18.1
10800,42.7
```

### 5.5 Closed-boundary mass balance check

Use `--closed-boundaries` to verify that no volume is lost during a simulation —
useful when validating a new domain or DEM.

```bash
julia --threads auto FloodModel.jl \
    --meshload my_mesh_sgs.parquet \
    --rainfall 50 --sim-duration 3600 \
    --closed-boundaries \
    --output closed_test.h5
```

Check the `mb_err` column in the log: it should be < 0.01% of domain volume
throughout.

### 5.6 Custom friction raster

```bash
julia --threads auto FloodModel.jl \
    --meshload my_mesh_sgs.parquet \
    --friction land_cover_manning.tif \
    --rainfall 30 --sim-duration 7200 \
    --output friction_test.h5
```

The friction raster should contain Manning's n values in the same CRS as the DEM.
Cells outside the raster extent use the global `--manning-n` value.

---

## 6. Visualisation

### 6.1 GLMakie viewer (`--vis makie`) — recommended

A native desktop window — no browser, no Cesium Ion token, no configuration
required. This is the recommended visualisation option.

- Pentagon polygons coloured by selected variable
- Dropdown: depth · saturation · volume · velocity
- Colorbar with auto-scaling
- Diagnostics sidebar: step, time, wet cells, max depth, mass balance error

```bash
julia --threads auto FloodModel.jl \
    --meshload my_mesh_sgs.parquet \
    --rainfall 30 --sim-duration 3600 \
    --vis makie
```

### 6.2 CesiumJS viewer (`--vis cesium` or `--vis`)

> **Note:** The CesiumJS viewer is currently paused in development. It is
> functional for basic use but may have rough edges. `--vis makie` is
> recommended for all current work.

A browser-based 3D viewer that updates live during simulation. Open
`http://localhost:8080` while the model is running.

**Controls:**
- Variable selector (dropdown) — depth, saturation, volume, velocity
- Colorbar with auto-scale toggle and manual ceiling slider
- 3D extrusion with adjustable vertical exaggeration
- Cell outline and zero-value cell visibility toggles
- Basemap selector: Google 3D Photorealistic Tiles / OSM / None
- Live mode (follows simulation) / Replay mode (scrub timeline)

**Google 3D Photorealistic Tiles** require a free Cesium Ion token. To set one up:

```bash
cp visualisation/cesium/config.example.json visualisation/cesium/config.json
# Edit config.json and set cesium_ion_token from ion.cesium.com/tokens
```

Add `visualisation/cesium/config.json` to `.gitignore` to keep your token out of
version control. The viewer falls back to OpenStreetMap tiles if the token is absent.

### 6.3 Post-processing viewer

After a simulation, view stored output interactively:

```bash
# Map at final frame, time series for a specific cell
julia --project=. visualisation/view_h5_output.jl sim.h5 \
    --cell 6345ddac19800000

# Map snapshot at t=7200s, volume variable
julia --project=. visualisation/view_h5_output.jl sim.h5 \
    --time 7200 --var volume

# Find a cell ID from coordinates
julia -e 'using .A5Grid; println(lonlat_to_cell(-2.895, 54.909, 18))'
```

### 6.4 Publication figures

```bash
# Wire-frame mesh, SGS cross-section, hypsometric curve for a cell
julia --project=. visualisation/visualise_mesh.jl \
    --mesh my_mesh_sgs.parquet 6345f2518e800000

# Save SVG + PNG to a directory without opening windows
julia --project=. visualisation/visualise_mesh.jl \
    --mesh my_mesh_sgs.parquet 6345f2518e800000 \
    --outdir figures/ --headless
```

---

## 7. Reading Output in Python

```python
import h5py
import numpy as np
import xarray as xr

with h5py.File("sim.h5") as f:
    # Static mesh
    cell_ids = f["mesh/cell_ids"][:].astype(str)
    lons     = f["mesh/center_lons"][:]
    lats     = f["mesh/center_lats"][:]

    # All frames
    frames = sorted(f["frames"].keys())
    times  = np.array([f[f"frames/{fr}/t"][()] for fr in frames])
    depths = np.stack([f[f"frames/{fr}/water_depth"][:] for fr in frames])

# As xarray Dataset
ds = xr.Dataset(
    {"water_depth": (["time", "cell"], depths)},
    coords={"time": times, "lon": ("cell", lons), "lat": ("cell", lats)}
)

# Maximum depth in the final frame
print(depths[-1].max())
```

See [DATA_FORMATS.md](DATA_FORMATS.md) for the complete HDF5 schema and GeoParquet
column definitions.

---

## 8. Performance Notes

**Mesh generation** is the slowest stage and only needs to be run once per
domain/resolution/DEM combination. Always save the mesh with `--meshout` and
use `--meshload` for subsequent runs.

**SGS pre-processing** at resolution 14 with the default 512 Halton samples takes
~5 s for a 2,900-cell domain. At resolution 18 (~30,000 cells) this scales to
several minutes. GPU acceleration (CUDA) reduces this roughly 10× for large domains.

**Simulation timestep** is set by the wave-speed CFL criterion and `--dt-max`.
A smaller `--dt-max` gives smoother output and better stability at the cost of
more steps. Typical values:

| Resolution | Recommended `--dt-max` |
|---|---|
| 14 | 60 s |
| 16 | 30 s |
| 18 | 10 s |
| 20 | 5 s |

**Multi-threading** helps most in the SGS table lookups and edge flux loops.
Use `julia --threads auto` to exploit all available cores during simulation. Use
`julia --threads 1` for mesh generation.

---

## 9. Troubleshooting

**`pya5` not found after setup.** Check that `setup.jl` completed the Python
package installation step without errors. Manually verify:

```bash
julia -e "using PyCall; pyimport(\"a5\")"
```

If this fails, try installing pya5 manually into PyCall's Python:

```julia
using PyCall
run(`$(PyCall.python) -m pip install pya5`)
```

**Zero edges in mesh.** The mesh was likely generated with a single A5 sublattice
(a known issue with older `pya5` versions). Regenerate the mesh with the current
`a5_bridge.py`. Check `mesh_summary(mesh)` — a well-formed mesh should have
approximately `2.5 × n_cells` edges.

**Simulation exits immediately with no output.** No water source was provided.
Add `--rainfall`, `--rainpoint`, `--injection-point`, or `--inflow-bci` to the run.

**`DomainError` or NaN in output.** Usually caused by very large initial fluxes
on a steep DEM. Reduce `--dt-max` (try `--dt-max 5`) and check that the DEM
covers the full AOI extent (`--dem-strict` will report cells outside the DEM).

**Mesh generation hangs.** Use `--threads 1` for mesh generation. The Python
bridge is not thread-safe during the mesh build phase.
