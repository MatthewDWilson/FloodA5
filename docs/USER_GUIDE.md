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

FloodA5 works with pya5 0.7.3 and 1.x. pya5 1.0 changed the cell-ID encoding, so
a mesh written with 0.7.x is refused under 1.x until you migrate it, using the
Python that has pya5 1.x:

```bash
python mesh/migrate_mesh_ids.py test/carlisle/carlisle_mesh18_standard.parquet test/carlisle/carlisle_mesh18_standard_a5v1.parquet
```

The migrated mesh gives bit-identical results. See `A5_QUIRKS.md` §11.

To upgrade an existing installation, `julia setup.jl --upgrade` installs pya5 1.x
into PyCall's Python (re-running `setup.jl` without it does not upgrade an
installed pya5) and lists meshes under `test/` that still use 0.x IDs;
`julia setup.jl --migrate-meshes` migrates them in place, keeping each original
as `*_a5v0.parquet`. `--pya5 "pya5==0.7.3"` rolls back. FloodA5 uses PyCall's
Python (by default Julia's own Conda.jl Python, `~/.julia/conda/3`), not any
other conda environment; `julia setup.jl --python PATH` switches it.

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

> **Thread count for mesh generation:** the recommended default is
> `--threads 1` for mesh generation (stages 1+2, i.e. any run using
> `--meshgen`). This is a conservative recommendation, not an absolute
> requirement: two specific crashes previously traced to threaded calls
> into the Python/GDAL bridge (`PyCall`'s garbage-collector finaliser firing
> on a non-main thread) have been fixed, and `--threads auto` mesh
> generation — including DEM sampling — has since run successfully.
> However, DEM sampling still makes coordinate-transform calls into the
> same Python/GDAL bridge from inside a threaded loop, the general pattern
> behind those crashes, and the failure mode was intermittent and
> mesh-size-dependent rather than immediate, so a clean run is encouraging
> but not a guarantee at every mesh size or resolution. **If mesh
> generation crashes with a low-level fault (typically reported as
> `EXCEPTION_ACCESS_VIOLATION` or similar) rather than a normal Julia
> error, retry with `--threads 1`** — this has always resolved it.
> Simulation runs (stage 3, `--meshload`) are unaffected either way and are
> always safe with `--threads auto`.

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
| `--meshgen FILE` | GeoJSON AOI file. Triggers A5 pentagon mesh generation via the Python bridge. |
| `--meshres N` | A5 resolution level. Required with `--meshgen`. |
| `--meshgen-cartesian FILE` | GeoJSON AOI file. Triggers Cartesian (4-connected, rectangular) mesh generation — pure Julia, no Python bridge (`mesh/CartesianGrid.jl`, `cartesian-formulation` branch). Covers the AOI's **bounding box**, not its polygon shape. Requires `--dx`. Mutually exclusive with `--meshgen`/`--meshload`. |
| `--dx METRES` | Cartesian cell size. Required with `--meshgen-cartesian`. |
| `--cartesian-lat-ref DEG` | Reference latitude for the `--dx` metres→degrees conversion (default: the AOI bounding box's own mid-latitude). Cells are exactly square in true metres only at this latitude, so a domain with large north–south extent shows ordinary equirectangular distortion away from it. Only meaningful with `--meshgen-cartesian`. |
| `--meshgen-dggs FILE` | GeoJSON AOI file. Generates a mesh of the discrete global grid named by `--dggs`: every cell that intersects the AOI (`mesh/DGGSGrid.jl` over `mesh/dggs_families.py`; `docs/MESHES.md`). Requires `--dggs` and one of `--meshres` or `--dx`. Use `--threads 1` for mesh generation. |
| `--dggs FAMILY` | `h3`, `igeo7`, `isea3h`, `healpix`, `rhealpix` (libraries `h3` and DGGAL, installed by `setup.jl`), or `a5` (the same as `--meshgen`). With `--meshres N` the grid's level N; with `--dx M` the level whose mean cell area is nearest M² — the levels are discrete (×3 to ×9 in area), so the cells are only roughly M across: the log line and the mesh file's metadata give the level. `python mesh/dggs_families.py FAMILY --aoi FILE --list-levels` lists the cell size of every level over an AOI. |
| `--meshout FILE` | Save generated mesh (any family) to `.parquet` (recommended) or `.geojson` (`--meshgen-dggs`: `.parquet` only). Omit to keep in memory only. |
| `--meshload FILE` | Load a previously saved GeoParquet mesh (any family — A5, Cartesian, the planar families of `mesh/mesh_families.py` and the DGGS meshes are interchangeable from this point in the pipeline onward), skipping generation. The mesh's metadata (family, DGGS level, …) is kept when the mesh is re-saved and is written to the output (`/mesh` attributes). |
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
| `--flow-model standard\|swe\|sgs` | `standard` (since 0.5.0; `sgs` before) | Hydraulic engine. `standard`: local-inertial on the mean cell elevation with the `--scheme` corrections. `swe`: the full shallow-water equations on the mean cell elevation (HLLC). It keeps inertia and advection: choose it where they matter (dam breaks, rapidly varying or supercritical flow, obstacles) or as the reference for the local-inertial engines (`HYDRAULICS.md` §15.10). It typically takes 5–10× the time steps, and the local-inertial scheme options do not apply to it. `sgs`: sub-grid hypsometric local-inertial (R-A kernel), opt-in. It needs a mesh built with `--flow-model sgs`, and its behaviour at high stage is under investigation. |
| `--swe-order 1\|2` | 2 | `swe` only: first order (Euler) or second order (MUSCL, Barth–Jespersen limiter, SSP-RK2). |
| `--backend cpu\|gpu\|ka-cpu` | `cpu` | (WP5) Where the time steps run. `cpu`: CPU threads. `gpu`: an NVIDIA GPU through CUDA.jl (full project environment); every flow model, with edge momentum and `--cfl-mode local` (the defaults); the state stays on the GPU and is copied back for output and log lines. Results agree with `cpu` to round-off (`ARCHITECTURE.md` §6a). `ka-cpu`: the GPU code path on the CPU, for testing. |
| `--precision float64\|float32` | `float64` | (0.6.0) Working precision of the time steps. `float32` needs `--backend gpu` (or `ka-cpu`) and `--flow-model standard` or `swe`: fluxes and momentum in single precision, cell volumes and the mass totals in double, water levels relative to a local datum — so mass is conserved as in a `float64` run, and depths agree with it to about 10⁻⁶ m over short runs (long runs drift apart as any two `float64` runs with different round-off do; `ARCHITECTURE.md` §6a). About twice the `float64` speed on GPUs whose double-precision arithmetic is slow. Output files stay Float64. Not available for `sgs`. |
| `--swe-viscosity C_E` | 0 | `swe` only: eddy viscosity ν_t = C_E·h·u* (u* the friction velocity; typical river values C_E ≈ 0.1–0.6). |
| `--sgs-bins N` | 100 | Elevation bins for SGS hypsometric curves. |
| `--sgs-samples N` | 512 | Halton points per cell for SGS cross-section pre-processing. |
| `--manning-n N` | 0.03 | Global Manning's roughness coefficient (s/m^(1/3)). |
| `--friction FILE` | — | GeoTIFF friction raster. Per-cell Manning's n sampled at cell centres; overrides `--manning-n` where the raster is finite. |
| `--sim-duration S` | 3600 | Simulation duration in seconds. |
| `--dt-max S` | 60 | Maximum adaptive timestep in seconds. Reduce to 10–30 s for high-resolution (level 18+) or steep-terrain runs. |

### 4.3a Directional-bias correction options (standard flow only)

**These flags are designed for, and validated only with, `--flow-model standard`.**
Under `--flow-model sgs`, `--face-flux-method diamond`
and `--momentum-model cell` are ignored (a startup warning is printed).
`--gradient-correction on`, `--gradient-correction-alpha` and
`--q-centre-theta` **do** affect `step_sgs!` (it applies the legacy WLSQ + V̂
correction), but that SGS path is not yet validated — a warning is printed.
(Earlier versions of this guide, and the startup warning, wrongly described
gradient correction as a no-op under SGS; corrected in WP0, 2026-10.)

All defaults below are the **legacy, uncorrected Bates (2010) behaviour** —
nothing in this table changes the default simulation output. They exist for
A/B comparison while the directional-bias correction work
(`flow-direction-fixes` and `directional-bias-reformulation` branches) is
validated on more test domains (see §8 of `HYDRAULICS.md` for the current
recommendation and open items).

| Flag | Default | Description |
|---|---|---|
| `--q-centre-mode legacy\|oriented\|consistent` | `legacy` (`consistent` under `--scheme recommended` and `v0.2`, standard flow; no effect at θ = 1) | (WP2c) How Q-centring uses its collinear partner faces. `legacy`: as stored — on A5 about two partners in three are oriented against the edge, which makes Q-centring a dt-dependent resistance (lower conveyance; `HYDRAULICS.md` §15.4). `oriented`: signs corrected, but pentagon partners are not parallel (deflects the flow). `consistent`: signs corrected and each partner mapped onto the edge normal — exact for uniform flow. All three are identical on Cartesian grids. |
| `--courant C` | `0.7` | Courant number of the adaptive timestep, in (0, 1]. |
| `--cell-order spatial\|file` | `spatial` | (WP4) Order of the cells in memory during the simulation. `spatial` sorts them along a Morton curve (1.4–1.7× faster steps); `file` keeps the mesh file's order. Results agree to round-off (but see `HYDRAULICS.md` §15.5 for the standard solver); HDF5 output is in the file's order either way. Ignored with legacy Q-centring, which depends on the order. |
| `--q-centre-theta N` | `0.9` (`1` under `--scheme recommended`) | Spatial momentum-smoothing parameter θ for the Q-centred scheme (checkerboard suppression), range `[0.0, 1.0]`. `0.9` is the LISFLOOD-FP standard (light smoothing); `1.0` disables the scheme (pure Bates semi-implicit momentum — may reveal period-2 checkerboarding on fine/irregular meshes). |
| `--gradient-correction on\|off` | `off` | Enables the WLSQ non-orthogonal gradient correction in place of the legacy `cos θ`-only scaling. When `on`, the driving-head term is constructed by `--face-flux-method` below. Must be given an explicit `on`/`off` value (a bare flag is an error, to avoid accidentally swallowing the next argument). |
| `--gradient-correction-alpha N` | `1.0` | Only applies when `--gradient-correction on --face-flux-method legacy`. Scales the tangential (non-orthogonal) correction term. `1.0` = full over-relaxed correction; `0.0` = orthogonal-only correction (direct WSE-difference term alone, no tangential term). **`0.0` is the empirically more stable interim setting** found during real-mesh testing (see `HYDRAULICS.md` §8) — `1.0` gives a larger correction but has been observed to overshoot into a mirrored bias under sustained ponding. Has no effect with `--face-flux-method diamond`. |
| `--face-flux-method legacy\|diamond\|diamond-wetdry` | `legacy` | Only applies when `--gradient-correction on`. Selects how the corrected driving head is constructed. `legacy`: cell-centred WLSQ gradient, averaged onto the face, then skewness-corrected. `diamond`: face-local reconstruction built directly from the two adjacent cell centres and the edge's own two shared vertices, with no cell-averaged-then-shared gradient. Proven exact for linear fields; falls back to `legacy` automatically on any individual edge whose diamond geometry is degenerate (rare, confined to AOI-boundary vertices — see `HYDRAULICS.md` §7). |
| `--momentum-model edge\|cell` | `edge` | Selects how previous-timestep momentum (`q_prev`) is represented. `edge`: one independent scalar per cell face (5 on A5, 4 on the Cartesian backend — see `--meshgen-cartesian` above) — the original Bates representation, derived for a rectangular grid where this is exact by construction. `cell`: a single 2D discharge vector per cell, reconstructed each step from all of that cell's face fluxes by WLSQ (Perot-style), then projected onto each face normal. Empirically the stronger of the two correction levers tested to date on A5 (see `HYDRAULICS.md` §7.3); as of this writing it is also the untested next step in diagnosing a separate, mesh-independent directional bias found on the Cartesian backend — see `HYDRAULICS.md` §14. |

| `--friction-model component\|vector` | `component` | Discharge magnitude in the friction term and the Froude limiter (WP1). `component`: the face-normal discharge only (Bates 2010 / LISFLOOD-FP ACC), which makes face conductance depend on face orientation — a directional bias that grid refinement cannot remove. `vector`: the full magnitude √(q_n² + q_t²), with q_t from the reconstructed cell vectors; rotationally consistent. Standard and SGS solvers. See `HYDRAULICS.md` §15. |
| `--momentum-model blend` | — | (WP2) θ·q_face + (1−θ)·(projection of the averaged cell vector), θ from `--q-centre-theta`; θ = 1 is edge momentum without smoothing, θ = 0 is `cell`. **Checkerboards on the Carlisle DEM at θ = 0.9** (as does `cell`); experimental. |
| `--limiter legacy\|positivity\|hybrid` | `legacy` | (WP2) `positivity`: one positivity-preserving scale factor per donor cell in place of the per-edge volume limiter and donor cap — **not stable on its own** on real terrain. `hybrid`: keeps the per-edge volume limiter and replaces only the per-edge donor cap; as smooth as `legacy`. |
| `--cfl-mode legacy\|local` | `legacy` | (WP2) `local`: dt = 0.7·min over wet cells of ℓ/√(g h), ℓ = 4·area/perimeter, instead of √(min area) with the 99th-percentile depth. About 25% more steps on Carlisle; removes the early rain-on-slopes noise. |
| `--boundary-limiter legacy\|positivity` | `legacy` | (D1b) `positivity`: open-boundary outflow is limited only by what the boundary cell holds, not by the 1/Δt volume limiter and donor cap. |
| `--open-boundary normal\|legacy` | `legacy` | (0.3.0) Open (zero-gradient) boundaries. `normal`: a free outflow at normal depth — the ghost cell continues the bed slope from the inward neighbour at the boundary cell's depth. `legacy`: before 0.3.0, the water-level slope from the most outward neighbour (along the edge), which held water back at outlets (`HYDRAULICS.md` §15.11). |
| `--eta-filter C_S` | `0` | (D4, standard flow) Conservative free-surface high-pass filter: water exchanged across wet faces in proportion to the jump in the residual of a local planar fit, at a rate C_s·√(gh)/ℓ that does not depend on Δt (`HYDRAULICS.md` §15.8). |
| `--momentum-damping C_M` | `0` | (D6′, standard flow; needs the diamond face gradient, edge momentum and vector friction) Exactly divergence-free damping of the circulating part of the face discharges, an eddy viscosity ν = C_m·V·ℓ at each mesh vertex (`HYDRAULICS.md` §15.9). |
| `--momentum-damping-scale wave\|flow` | `wave` | V = √(gh) (`wave`) or the local speed \|u\| (`flow`). |

`--face-flux-method diamond-wetdry` (WP2) is `diamond` with the two-point face
gradient wherever an edge's vertex stencil includes a dry cell.

**Defaults (0.3.0): `--scheme recommended`.** Any of
`--gradient-correction`, `--face-flux-method`, `--friction-model`,
`--limiter`, `--cfl-mode`, `--momentum-model`, `--q-centre-mode`,
`--q-centre-theta`, `--boundary-limiter`, `--eta-filter`, `--momentum-damping`
and `--momentum-damping-scale` that you do not give takes the scheme's value
for the flow model; explicit flags always win. `--scheme v0.2` restores the
0.2.x default and `--scheme legacy` the pre-2026-10 behaviour, both exactly.

| Flow model | `--scheme recommended` (default) | `--scheme v0.2` |
|---|---|---|
| standard | `--gradient-correction on --face-flux-method diamond-wetdry --friction-model vector --limiter positivity --boundary-limiter positivity --cfl-mode local --q-centre-theta 1 --eta-filter 0.05 --open-boundary normal` | `--gradient-correction on --face-flux-method diamond-wetdry --friction-model vector --limiter hybrid --cfl-mode local --q-centre-mode consistent` |
| sgs | `--friction-model vector --cfl-mode local --q-centre-theta 1 --open-boundary normal` | `--friction-model vector --cfl-mode local --q-centre-theta 1` |

`--momentum-damping` is off in every scheme. See `HYDRAULICS.md` §15.11 for
the evidence. The defaults listed in the
table above for the individual flags are their `--scheme legacy` values.

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
| `--output FILE` | off | Write HDF5 simulation output (`.h5` / `.hdf5`). Since 0.4.0 (format 2) the file also records the version, git commit, command line, resolved options, input files and the mesh, so the run can be reproduced with `--rerun` and replayed or analysed without the mesh file ([DATA_FORMATS.md](DATA_FORMATS.md) §2). |
| `--output-interval S` | 60 | Seconds of simulation time between output snapshots. |
| `--embed-mesh full\|geometry` | `full` | `full`: store a copy of the mesh file (GeoParquet) in the output. `geometry`: store only the cell polygons, neighbours and the mesh file's checksum — for sweeps on one large SGS mesh, whose tables would otherwise dominate every output (Carlisle res-18 SGS: 179 MB). `--rerun` then needs the original mesh file, unchanged. |
| `--rerun FILE.h5` | — | Run again exactly the simulation recorded in `FILE.h5`, from its embedded mesh and input files (output `FILE_rerun.h5` unless `--output` is given). Any other flag given with it replaces the recorded one, e.g. `--rerun a.h5 --backend gpu`. With the same code and backend the frames are bit-identical; a warning is printed if the code differs from the recorded commit. |

### 4.7 Visualisation options

| Flag | Default | Description |
|---|---|---|
| `--vis [MODE]` | off | Enable live visualisation. Modes: `makie` (GLMakie window; also a bare `--vis`, since 0.8.3), `cesium` and `giro3d` (web globe viewers, §6.2). |
| `--vis-port PORT` | 8080 | HTTP port for the globe viewer's server. Ignored for `--vis makie`. |

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

### 6.1 GLMakie viewer (`--vis makie` or `--vis`)

A native desktop window — no browser, no Cesium ion token, no configuration
required. A bare `--vis` gives this viewer (since 0.8.3; before, it gave the
CesiumJS one).

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

**Replay of an output file.** The same window replays any run from its HDF5
output — whether or not it was visualised live, and while it is still
running (frames written so far):

```bash
julia --project=. --threads auto visualisation/replay.jl sim.h5
julia --project=. --threads auto visualisation/replay.jl sim.h5 --var velocity --from 3600 --fps 8
julia --project=. --threads auto FloodModel.jl --replay sim.h5          # the same
julia --project=. visualisation/replay.jl sim.h5 --from 7200 --png t2h.png   # one frame to PNG, no window
```

A frame slider scrubs to any output time, a speed slider sets the playback
rate (output frames per second) and Pause/Play drives playback; the display
menu and flow arrows work on the shown frame. The time-series strip shows
the whole run with a dotted marker at the shown frame:

- *Domain* volume and wet cells: read from the file (every frame).
- *Added*: rebuilt from the run record (rainfall, point sources and the
  embedded hydrograph / `.bci` / `.bdy` files, integrated as the model's
  mass balance) and checked against the recorded end-of-run total; not
  shown if it cannot be rebuilt or does not match, and never for format-1
  files.
- *Removed*: the file stores outflow only at the end of the run. Closed
  boundaries remove nothing; for open boundaries the strip shows
  added − domain (= removed + mass-balance error, labelled so) and the
  recorded end-of-run removed volume and error as ◆ markers.

Options: `--fps N`, `--var depth|saturation|volume|velocity|ring`, `--from T`
(s), `--paused`, `--png FILE`, `--mesh FILE.parquet` (format-1 files only:
they carry no mesh geometry). In Julia: `rd = FloodA5Replay.open_replay("sim.h5")`
(`visualisation/replay_data.jl`), then `MakieVisualiser.replay(rd)`.

### 6.2 Globe viewers: CesiumJS (`--vis cesium`) and Giro3D (`--vis giro3d`)

A web viewer on a 3D globe: the water surface is drawn at its height, so
with Google Photorealistic 3D Tiles flooded buildings stand out of the water.
It follows a running simulation live, or replays an output file. There are
two, with the same controls and behaviour on two rendering engines —
**CesiumJS 1.147** and **Giro3D 2.0.7** (three.js) — while we compare them
(docs/GLOBE_VIEWERS.md). The server always offers both:
`/viz/cesium/index.html` and `/viz/giro3d/index.html`; `--vis` picks the one
`/` opens. Open the printed URL (default `http://localhost:8080`) in a
browser with WebGL 2 (any current Chrome, Edge or Firefox).

```bash
# live: the viewer follows the run; after the run it keeps serving until Ctrl-C
julia --project=. --threads auto FloodModel.jl --meshload mesh.parquet \
    --rainfall 30 --sim-duration 3600 --output run.h5 --vis cesium [--vis-port 8091]
#   … --vis giro3d                                       (the Giro3D viewer)

# replay an output file (format 2: the mesh is in the file)
julia --project=. visualisation/web_replay.jl run.h5 [--viewer giro3d] [--port 8091] [--offset 53.2] [--autoplay] [--open]
julia --project=. FloodModel.jl --replay run.h5 --vis cesium|giro3d [--vis-port 8091]   # the same

# no run: open output files from this computer in the page (read in the browser)
julia --project=. visualisation/web_replay.jl [--viewer giro3d] [--port 8091] [--open]
```

A replay opens at the first frame and waits for *Play*; `--autoplay` (or
`?autoplay=1` in the URL) starts playing when the page opens.

**What is drawn.** One surface for the whole mesh, each cell flat at its
water level: `datum + depth`, where the datum is the cell elevation (standard,
SWE) or the cell's lowest sub-grid elevation (SGS). Cells shallower than the
*wet if depth ≥* threshold (default 5 cm) are not drawn. The surface is translucent and
depth-tested, so the photorealistic tiles (or terrain) hide it where the
ground or a building is higher — wet cells partly under buildings show as
water between them. *Depth × exaggeration* draws `bed + depth × N` instead
(for flat basemaps). *Bed surface* draws the mesh elevations, coloured by
height: the default on the flat basemaps, where there is no terrain. It has
its own colours: click its colour bar for a palette (default *terrain*),
a range (*auto*: 1st–99th percentiles of the bed elevations, or *fixed*)
and an opacity (below 1 the bed is translucent).
*Cell borders* draws each cell's outline as a white hairline (off by
default; faded out where cells are smaller than about 8 pixels). Hover over
a cell for its depth, level, bed and id.

When zoomed out, 3D tiles and terrain are drawn from coarse levels of
detail that stand metres above the ground they simplify; the water is
therefore drawn slightly nearer the camera than it is — up to 1.5 % of its
distance, from 800 m away (nothing closer, so close views keep exact
occlusion) — so that shallow water is not hidden in large triangles.

**Colours.** Colour by water depth, water level, speed or (SGS) wetted
fraction. The range is, by default, fixed for the whole run: the 2nd to
98th percentiles over the wet cells of up to 24 frames spread over the run
(the server's `/api/ranges`; a live run's are updated as it goes). *frame*
uses the 1st–99th percentiles of the frame shown, *fixed* your own limits.
Click the colour bar to choose a palette: water (blue), brown flood water
(a narrow range around the colour of open brown floodwater in aerial
photographs), viridis, cividis, magma, turbo, blues, greys. The choice is
remembered per variable in the browser.

**Hazard variables** (since 0.8.9; runs with depth and speed). The variable
menu also offers, for each frame: *Depth × velocity* (m²/s); the Defra/EA
*Hazard rating* HR = d (v + 0.5) + DF (FD2321, with a debris factor DF by
depth, velocity and land use); the *Hazard class (ARR H1–H6)* of the
Australian combined hazard curves; and the *Hazard degree* of Chen et al.
(2019): speed over the critical speed at which a person of a given mass and
height is toppled (0–1; 1 = swept off their feet). The viewer's server
computes them from each frame's depth and speed with the same code as
`FloodA5Analysis` (`analysis/hazard.jl`: `FA.hazard_rating`,
`FA.hazard_degree`, …). *Hazard parameters* (in the Layers tab, kept in the
browser) set the land use for the debris factor — one for the whole run:
*any* (FD2320's precautionary 0.5 / 1, the default), *pasture/arable*,
*woodland*, *urban* (FD2321 Table 3.1) or *no debris*; a land-use map is a
later development — the person (adult 75.1 kg, 1.77 m; child 25.5 kg,
1.26 m; or your own) and the calibration of Chen et al. (θ = 75°,
α = 3.064, β = 0.303 recommended), and the depth from which a cell counts as
flooded (default 0.1 m). The *About* tab gives the definitions, classes,
references and limitations (Chen et al.: toppling only).

**Classes.** Water depth, D·V and the hazard rating can be drawn in classes
(*Classes*; the hazard rating is classed by default, the ARR class always):
depth bands 0.15 / 0.3 / 0.6 / 1.2 m; D·V 0.3 / 0.6 / 1.0 / 4.0 m²/s (the
ARR curves' D·V limits); HR 0.75 / 1.25 / 2.0 (FD2320: low, danger for
some, for most, for all). The legend then lists the classes; the breaks can
be edited (*Reset* restores them) and are remembered in the browser.

**Statistics over the run** (since 0.8.9). *Maximum depth*, *maximum speed*,
the maximum of each hazard variable, the *time of first flooding* and the
*duration of flooding* (depth ≥ the flooded threshold) are maps over the
whole run: while one is shown the time slider and playback are off, and the
water is drawn at each cell's maximum depth (cells never flooded are left
out). They are computed from the output frames (a live run: from the frames
the server keeps, recomputed as frames arrive), so at frame resolution —
peaks between frames are missed. In Julia: `FA.statistics_maps(run; params)`
(`FA.run_statistics` over arrays). With SGS the depth is that of the
cell's lowest part, so depth-based hazards describe the deepest part of a
cell.

**Run and About.** The panel's *Run* tab shows what the run records: the
FloodA5 version and commit, flow model, backend, the command line (*copy*
puts it on the clipboard), the mesh, frames, input files, all options and
the end-of-run summary (a replay: from the file; a live run: the record its
output file gets). *About* describes each variable and the colour range,
vertical offset, wet threshold and depth bias. The descriptions are in
`io/variables.toml`, which the static figures' captions (§6.5) also use.

**Time series at cells.** Right-click a cell on the map (a click — dragging
with the right button still turns or zooms the view) and choose *Add cell N
to plot*. A plot opens at the bottom of the window (⇆ docks it on the
right instead); each further cell adds a line, up to eight, with a numbered
marker of the same colour on the map. The plot has its own variable (depth,
water level — drawn only where the cell is wet —, speed, wetted fraction,
or a hazard variable, with the hazard parameters; not the statistics)
and shows the whole run, with a dashed line at the frame on the map; click
the plot to go to that time. ✕ in the legend removes a line; *CSV*
downloads the lines (time in seconds, one column per cell). In a live run
the lines grow as frames arrive (the frames the viewer keeps: thinned in
long runs). The plots use Chart.js, loaded from jsDelivr when the first
plot opens; for offline use point `FLOODA5_CHARTJS_URL` at a local copy of
`chart.umd.min.js` (from the `chart.js` package, served by any local web server).

**Files from this computer** (since 0.8.10). *Open file…* (top bar), or a
file dropped anywhere on the page, opens a FloodA5 output file (format 2:
FloodA5 0.4 or later) from your computer in this page. The browser reads it
itself — HDF5 compiled to WebAssembly (h5wasm 0.10.3) in a background worker,
reading only the parts needed (a frame, a cell's series) — and nothing is
uploaded, so a viewer hosted elsewhere never sees the data. Everything works
as in a replay, with the same results (the server's code ported to
JavaScript; `test/test_vis_local.jl` checks them against each other): the
variables, hazard variables, statistics over the run, colour ranges, the
time-series plot and the *Run* tab. The file replaces the run in this page
only (other pages, and the server, are not affected; reload the page to go
back to the server's run). `web_replay.jl` without a file starts a server
with no run that opens with a prompt to choose a file, and the viewer pages
also work **without a server** — e.g. the `visualisation/` directory on a
static web host (the settings then come from
`viewer-common/defaults.json`, the *About* text from
`viewer-common/about.json`; there is no live view and no Cesium ion token,
so the basemaps are OpenStreetMap and Giro3D's Terrarium terrain).
Opening the 400 MB, 30 000-cell Carlisle run takes 3–8 s (h5wasm, about 4 MB, is
loaded from jsDelivr the first time; for offline use copy `hdf5_hl.js` and
`hdf5_util.js` from the h5wasm package's `dist/esm/` to
`visualisation/viewer-common/vendor/` (git-ignored) and set
`FLOODA5_H5WASM_URL=/viz/viewer-common/vendor/hdf5_hl.js`);
a statistic over its 112 frames about half a second. Format-1 files (≤ 0.3.x)
have no cell polygons: replay them with `web_replay.jl --mesh`.

**Vertical offset.** Model heights are in the DEM's datum (usually mean sea
level: a geoid); Google tiles and Cesium World Terrain use heights above the
WGS84 ellipsoid. The *vertical offset* is added to every model height to
convert: there, the geoid height — between about −107 and +86 m worldwide
(look it up for your area, e.g. with PROJ and the national geoid grid).
Each basemap has its own offset, since their heights need not share a datum:
at Carlisle the Google tiles need +52.2 m, Cesium World Terrain about
+52.9 m, Giro3D's Terrarium terrain (above mean sea level) about 0 m. The
viewer uses, in order:

1. `--offset` (replay) or `vertical_offset` in `config.json`, for a basemap
   with ellipsoidal heights (Google tiles, Cesium World Terrain);
2. the value saved in this browser for this domain and basemap (set by hand
   or estimated before);
3. on a flat basemap (OpenStreetMap, none): minus the 2nd percentile of the
   bed elevations, so the model's lowest ground sits on the map;
4. otherwise an estimate from the basemap, made when the basemap is first
   shown and then saved.

*Estimate* redoes the estimate: it compares the cell elevations with the
basemap's heights at about 300 cells (Google tiles: the 30th percentile,
since buildings and trees only raise them; terrain: the median) and reports
the spread. The panel says where the current value came from.

**Live frames.** The server keeps one frame per second of wall-clock time
(`FLOODA5_VIS_INTERVAL`, seconds) within 1 GB (`FLOODA5_VIS_BUDGET_MB`): when
the budget is reached it drops every second frame, so the whole run stays
viewable at a coarser time step. The end state is always sent. The HDF5
output keeps every output interval; replay it for the full record. A page
opened after the run ended, or reloaded, shows the complete run.
`FLOODA5_VIS_LINGER=S` stops serving S seconds after the run instead of
waiting for Ctrl-C (scripted runs).

**Controls.** Timeline slider; first / previous / play / next / last
(keys: space, ←, →); playback speed in frames per second or simulation time
per second; *LIVE* follows the newest frame while the run goes on. A replay,
or a run that has finished, opens at the first frame.

**Basemaps and the Cesium ion token.** Both viewers: *Google
Photorealistic 3D*, *OpenStreetMap (flat)*, *None*. Terrain: *Cesium World
Terrain + OSM* (CesiumJS; needs the token) or *Terrarium terrain + OSM*
(Giro3D; AWS open elevation tiles, no key). *Google Photorealistic 3D* and
*Cesium World Terrain* need a (free) Cesium ion access token
(ion.cesium.com → Access Tokens). Put it in a **keys file** (since 0.8.5):
the same on every OS, read by every new run, never committed:

```bash
# per user, shared by every clone and worktree (recommended)
mkdir ~/.flooda5            # PowerShell: mkdir $HOME\.flooda5
cp keys.example.json ~/.flooda5/keys.json
# …or in the repository root:  cp keys.example.json keys.json   (git-ignored)
# then edit the file: "FLOODA5_CESIUM_ION_TOKEN": "eyJ…"
```

Keys are looked up in this order (`io/FloodA5Keys.jl`): the environment
variable of the same name (for CI and servers), the file named by
`FLOODA5_KEYS`, `keys.json` in the repository root, `~/.flooda5/keys.json`.
The viewer's server logs where its token came from (never the value), or
that there is none. Environment variables work too, but on Windows a value
set with `setx` reaches only processes started afterwards from Explorer — a
VS Code that was already running (even after closing its windows, if a
`Code.exe` stayed alive) keeps the old environment. The token is passed to
the browser, so keep it scoped (below).

Scope the token to `assets:read` and restrict its allowed URLs to
`http://localhost:*` in the ion dashboard. A Google Map Tiles API key
(`FLOODA5_GOOGLE_TILES_KEY`) is the other route to the Google tiles.
Without either, the viewer uses OpenStreetMap on a flat globe. Google's
terms: the Google Maps logo (CesiumJS) or the words *Google Maps* (Giro3D)
and the data attributions are shown; Google tiles must not be combined with
Bing imagery (the viewers use none).
Optional settings go in `visualisation/viewer-common/config.json`
(git-ignored; see `config.example.json` and docs/DATA_FORMATS.md §7). Both
viewers load their engine from a CDN (CesiumJS from cesium.com; Giro3D,
three.js and OpenLayers as ES modules from esm.sh, pinned in
`visualisation/giro3d/index.html`). For offline use of the CesiumJS viewer,
point `FLOODA5_CESIUM_URL` at a local CesiumJS build (`…/Build/Cesium/`).

The server listens on 127.0.0.1 only and answers only requests addressed to
localhost.

### 6.3 Post-processing viewer

After a simulation, view stored output interactively (for the full viewer,
replay the run: §6.1). Format-2 files (0.4.0+) carry the mesh, so
`--animate` and `visualisation/FloodViewer.jl sim.h5` need no mesh file;
format-1 files need it (`--animate --mesh mesh.parquet`,
`FloodViewer.jl mesh.parquet sim.h5`).

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

### 6.5 Static figures from output files (CairoMakie)

`visualisation/FloodA5Plots.jl` (`FP`) draws maps, comparisons and time
series of finished runs for notebooks and papers. It needs only the output
files (format 2 carries the mesh; for ≤ 0.3.x files load runs with
`parquet = …`) and the analysis module — not `FloodModel.jl` or Python.

```julia
include("visualisation/FloodA5Plots.jl"); const FP = FloodA5Plots; const FA = FP.FA
using CairoMakie
li  = FA.load_run("test/case_runs/carlisle/runs/a5_li.h5")
swe = FA.load_run("test/case_runs/carlisle/runs/a5_swe.h5")

FP.fieldmap(li; t = 4 * 3600)                    # depth map at 4 h (also var = :speed, :wse, …)
FP.fieldmaps([li, swe], [2, 4, 6] .* 3600)       # runs × times, one colour range
FP.difference_map(swe, li; t = 6 * 3600)         # Δh = li − swe on cells wet in either
FP.extent_map(swe, li; t = 6 * 3600)             # hits / misses / false alarms, CSI
FP.storage_plot([li, swe])                       # storage against time
FP.mass_balance_plot(li)                         # Σ volume, dS/dt, summary check
FP.timeseries([li, swe], [(lon = -2.925, lat = 54.897, label = "floodplain")])
save("maps.png", FP.fieldmaps([li, swe], [7200.0]); px_per_unit = 2)

# any per-cell values, on a window of cells; or a mesh drawn by its outlines
res = FA.wse_residual(li, FA.nframes(li))       # η − mean of neighbours (a checkerboard measure)
FP.cellmap(li, res.r; cells = FP.cells_within(li, -2.912, 54.898, 500), origin = (-2.912, 54.898))
FP.cellmap(li.mesh; outline = true, cells = FP.cells_within(li, -2.912, 54.898, 150))
```

- Each function returning a `Figure` has a `…!(fig[i, j], …)` form for
  composing panels; map panels return a `MapPanel` whose `coords` and
  `range` can be passed to the next panel, and `FP.mark_points!(panel,
  points; labels)` adds gauges or labels. `FP.activate!()` applies the house
  style to figures you lay out yourself.
- Maps draw the mesh polygons: `coords = :lonlat` (degrees, drawn to scale)
  or `:local` (metres about `origin = :centre`, `:source` or `(lon, lat)`).
  `fieldmap` takes a vector of per-cell values as `var`; `cellmap` draws any
  values (diverging colours by default) on a subset of cells
  (`FP.cells_within(run, lon, lat, halfwidth_m)`), or the cells alone.
  Cells with depth ≤ `dry_threshold` (0.01 m) show the `basemap`
  (`:elevation` greys, `:grey`, `:none`). Colour ranges default to the
  rounded 98th percentile of the wet values (symmetric for differences).
- Times are matched to each run's nearest frame (within 60 % of the frame
  interval); titles give the time as hh:mm.
- Line plots style runs by mesh and solver when that identifies them (A5
  blue, Cartesian orange; local-inertial solid, SWE dashed, SGS dotted),
  otherwise by series colour (`style = :series`).
- Captions (since 0.8.7): under each figure, one line per run (label, flow
  model, mesh, FloodA5 version and commit, from the file's run record) and
  the variable's definition (`io/variables.toml`, as the globe viewers'
  *About* tab). `caption = false` leaves them out of one figure (for
  export or publication); `FP.CAPTIONS[] = false` everywhere.
  `FP.figure_caption(runs, var)` gives the text alone.
- Speed: ~1 s per Carlisle-size map (30k cells) after the first call
  (~20 s of compilation); ~6 s for 5 × 10⁵ cells.

---

## 7. Reading Output

In Julia, read output through `io/FloodA5Output.jl` (it reads both the
current format and files from FloodA5 ≤ 0.3.x):

```julia
include("io/FloodA5Output.jl"); const FO = FloodA5Output
t     = FO.frame_times("sim.h5")
h     = FO.read_var("sim.h5", "water_depth")      # n_cells × n_frames
h_end = FO.read_frame("sim.h5", "water_depth")    # last frame
FO.run_info("sim.h5")["argv"]                     # the command that produced it
```

For analysis, use `analysis/FloodA5Analysis.jl`. It loads a run together with
its mesh geometry and computes the project's metrics, in pure Julia (no model,
no Python):

```julia
include("analysis/FloodA5Analysis.jl"); const FA = FloodA5Analysis
ref = FA.load_run("swe.h5"); run = FA.load_run("li.h5")
FA.compare(ref, run)                      # storage ratio, rms Δh, CSI, bias at the last common time
k = FA.frame_at(run, 3600.0)              # frame nearest t = 1 h
FA.wse(run, k)                            # water surface (SGS-aware)
# files from FloodA5 ≤ 0.3.x: pass the mesh for areas, polygons and neighbours
old = FA.load_run("old.h5"; parquet = "mesh.parquet")
```

In Python:

```python
import h5py
import numpy as np
import xarray as xr

with h5py.File("sim.h5") as f:
    cell_ids = f["mesh/cell_ids"][:].astype(str)
    lons     = f["mesh/center_lons"][:]
    lats     = f["mesh/center_lats"][:]
    times    = f["results/t"][:]
    depths   = f["results/water_depth"][:]        # (n_frames, n_cells)

ds = xr.Dataset(
    {"water_depth": (["time", "cell"], depths)},
    coords={"time": times, "lon": ("cell", lons), "lat": ("cell", lats)}
)
print(depths[-1].max())   # maximum depth in the final frame
```

See [DATA_FORMATS.md](DATA_FORMATS.md) §2 for the complete HDF5 schema (run
record, cell polygons, embedded mesh, format 1) and §1 for the GeoParquet
columns.

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
Use `julia --threads auto` to exploit all available cores during simulation.
For mesh generation, `--threads 1` is the recommended default — see the
callout in §2 for the full explanation and the fallback if you use
`--threads auto` and hit a crash.

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

**Mesh generation crashes with a low-level fault (e.g.
`EXCEPTION_ACCESS_VIOLATION`), rather than a normal Julia error.** This is a
known interaction between Julia's garbage collector and the Python/GDAL
bridge (`PyCall`) when a coordinate-transform call happens on a non-main
thread. Retry with `--threads 1` — see the callout in §2 for what's
actually going on and why `--threads auto` sometimes works fine anyway.
