![Status: In Development](https://img.shields.io/badge/Status-In_Development-green)
![Status: In Development](https://img.shields.io/badge/Status-Experimental-red)
![Version 0.2.0](https://img.shields.io/badge/Version-0.2.0-lightgrey)
![Julia 1.12+](https://img.shields.io/badge/Julia-1.12%2B-blue)
![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue)

# FloodA5

A 2D surface-water flood model built on the [A5 pentagonal Discrete Global Grid System
(DGGS)](https://a5geo.org). FloodA5 discretises the landscape into equal-area pentagons,
routes water between them using the Bates et al. (2010) inertial shallow-water
formulation, and resolves sub-cell terrain variability through pre-computed hypsometric
storage tables.

Note: this model is still in an experimental phase and should not be used in production. 

Key characteristics:

- **Uniform five-connectivity** — every interior cell has exactly five edge-sharing
  neighbours, avoiding the axis-aligned preference of a rectangular grid in
  principle. A5 faces are not perpendicular to the line between cell centres
  (non-orthogonality ~16–38°), and the original Bates (2010) scheme, derived for an
  orthogonal grid, is directionally biased on them. The bias has two causes — the
  face-normal water-surface gradient and the face-normal-only friction term — and
  **the default configuration corrects both** (a "diamond" face gradient and vector
  friction), cutting the planar-slope north/south asymmetry from 0.87 to about 0.05.
  See [Solvers](#solvers) and [docs/HYDRAULICS.md](docs/HYDRAULICS.md) §15.
- **Sub-Grid Sampling (SGS)** — hypsometric volume curves built from LiDAR allow
  partial wetting of cells and accurate routing through channels narrower than a cell.
- **Other mesh types through the same solver** — a rectangular Cartesian backend
  (`--meshgen-cartesian`, pure Julia; given a raster it builds a grid aligned
  cell-for-cell with it in its own projected CRS, for comparisons with LISFLOOD-FP),
  and hexagonal, rotated-Cartesian, Voronoi and imported polygon meshes
  (`mesh/mesh_families.py`). The solver and the SGS tables work for any number of
  neighbours per cell. The Cartesian grid also serves as an orthogonal reference that
  separates grid-geometry effects from solver effects.
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
| [docs/HYDRAULICS.md](docs/HYDRAULICS.md) | Flux-kernel derivations, stability limiters, the directional-bias work, and (§15) the evidence behind the current defaults, including what was tried and rejected |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | Full CLI reference, workflow walkthrough, worked examples |
| [docs/DATA_FORMATS.md](docs/DATA_FORMATS.md) | GeoParquet schema, HDF5 layout, `.bci`/`.bdy` format |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | Module map, extension points, developer orientation |
| [docs/A5_QUIRKS.md](docs/A5_QUIRKS.md) | A5-specific gotchas: dual sublattice, ID padding, mesh-generation thread safety |
| [docs/TEST_CASES.md](docs/TEST_CASES.md) | Square domain, planar-slope, and Carlisle 2005 flood test cases |
| [CHANGELOG.md](CHANGELOG.md) | Version history |
| [VERSIONING.md](VERSIONING.md) | This project's versioning policy |

**Not included in this public release:** the project's session-by-session
development history — implementation plans, debugging handovers, and
detailed experimental records, including negative results and reversed
decisions — is retained internally to support continued development and is
available on request. The documents above summarise the current state of
the model; they don't duplicate that history, but do point to it where it's
relevant (e.g. `HYDRAULICS.md`'s directional-bias sections).

**Current defaults (since October 2026):** `--scheme recommended`, i.e. the
directional-bias corrections, vector friction, a local CFL timestep and
consistent Q-centring — see [Solvers](#solvers). `--scheme legacy` reproduces
the original Bates (2010) configuration exactly. Results from earlier versions
differ: the October 2026 changes (listed in [CHANGELOG.md](CHANGELOG.md) and
marked ⚠ there) include two bugs that affected all earlier SGS and A5 results.

---

## Quick Start

### 1. Prerequisites

**Julia 1.12 or later** — download from [julialang.org](https://julialang.org/downloads/).

**Python 3.10 or later** with `pya5` and supporting packages — installed automatically
by the setup script (see below). Both pya5 0.7.x and 1.x are supported; pya5 1.0
changed the A5 cell-ID encoding, so meshes record their ID scheme and existing
meshes are migrated with `julia setup.jl --upgrade` then
`julia setup.jl --migrate-meshes` (see [docs/USER_GUIDE.md §1](docs/USER_GUIDE.md)
and [docs/A5_QUIRKS.md §11](docs/A5_QUIRKS.md)).

Optional: an NVIDIA CUDA GPU accelerates the point-in-polygon sampling step of
mesh generation at high resolutions. CUDA.jl is loaded only when that step needs
it; `FLOODA5_NO_CUDA=1` disables it.

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
│   ├── CartesianGrid.jl       # Cartesian (rectangular) mesh backend, pure Julia
│   ├── DiamondFlux.jl         # Diamond face-normal gradient reconstruction
│   ├── a5_bridge.py           # Python bridge: pya5 calls (0.x and 1.x), GeoParquet I/O
│   ├── migrate_mesh_ids.py    # pya5 0.x → 1.x cell-ID migration for existing meshes
│   ├── repair_sgs_slot_order.py # Repair SGS tables built before the 2026-10 slot-order fix
│   ├── mesh_families.py       # Hexagonal, rotated-Cartesian, Voronoi and imported polygon meshes
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
│   ├── planar_embankment/     # Planar-slope (flow-direction) test case
│   ├── square/                # Square domain test case
│   ├── synthetic_dem/         # Synthetic DEM SGS validation suite
│   └── *.jl                   # Unit and integration tests, diagnostics, benchmarks
├── ci/                        # CI environment and test runner (run_tests.sh)
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
| 10 | ~32.4 km² | ~5.7 km | Large catchment |
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

**`--flow-model standard`** — the same local-inertial scheme on the mean cell
elevation (no sub-grid tables). Faster, and useful for exploration or for runs
without a DEM.

**`--flow-model swe`** — a full shallow-water reference solver (HLLC, hydrostatic
reconstruction, MUSCL + SSP-RK2 with `--swe-order 2`, the default) on the same meshes,
beds, sources and boundaries. About ten times the cost of `standard`; it is there to
benchmark the local-inertial schemes against (`docs/HYDRAULICS.md` §15.10), not for
production runs.

Both solvers use the Bates et al. (2010) local-inertial formulation with a
Froude limiter (Fr ≤ 0.8), a wave-speed CFL timestep (`--courant`, default 0.7)
and a momentum state that matches the flux actually transferred.

**Default numerical options: `--scheme recommended`.** Any option you don't set
takes the value for your flow model; explicit flags always win, and
`--scheme legacy` restores the original configuration.

| Option | `standard` | `sgs` | `--scheme legacy` |
|---|---|---|---|
| Face gradient | diamond, two-point next to dry cells (`--gradient-correction on --face-flux-method diamond-wetdry`) | two-point | two-point |
| Friction | vector, \|q\| (`--friction-model vector`) | vector | face-normal only |
| Wet/dry limiter | hybrid (`--limiter hybrid`) | per-edge (legacy) | per-edge |
| Timestep | per-cell CFL (`--cfl-mode local`) | per-cell CFL | smallest cell, 99th-percentile depth |
| Momentum smoothing | consistent Q-centring, θ = 0.9 (`--q-centre-mode consistent`) | none (`--q-centre-theta 1`) | Q-centring θ = 0.9 as originally written |

Why these:
- **Directional bias** has two causes: the face-normal gradient on
  non-orthogonal faces, and the Bates/LISFLOOD-FP friction term, which uses
  only the face-normal discharge (this one also biases Cartesian grids, and
  refinement cannot remove it). Fixing both cuts the planar-slope N/S
  asymmetry on A5 from 0.87 to about 0.05, decreasing with refinement.
- **Q-centring as originally written is wrong on A5:** about two thirds of
  the "collinear" partner faces it averages were added with the wrong sign,
  so it acted as an artificial, timestep-dependent resistance. On the
  Carlisle test about twice as much water stayed in the domain as on a
  Cartesian grid. The consistent form maps each partner face onto the
  edge's own direction; on Cartesian grids it is exactly LISFLOOD-FP's
  Q-centring.

> **Results from before October 2026 are affected** by that Q-centring defect
> (A5, both solvers) and, for SGS, by a table bug that attached edge sills and
> R–A curves to the wrong neighbour for most cells. Old SGS meshes are refused
> until rebuilt or repaired (`mesh/repair_sgs_slot_order.py`).

Known limitations: with the default options, the standard solver still shows
more grid-scale water-surface noise in deep channels than the legacy scheme
did (the old Q-centring suppressed it, at the cost of conveyance). Results on
the Carlisle test are also sensitive to the Courant number. Both are open; see
[docs/HYDRAULICS.md §15.3–15.4](docs/HYDRAULICS.md), which also records the
options that were tried and rejected. The CLI flags are documented in
[docs/USER_GUIDE.md §4.3a](docs/USER_GUIDE.md), and the physics in
[docs/METHODS.md](docs/METHODS.md).

**Performance.** Each time step is multi-threaded (`julia --threads auto`).
Carlisle res-18 (30,000 cells, river inflows plus 10 mm/h rain, 8 threads)
runs at roughly 500 simulated hours per wall-clock hour with the standard
solver and 125 with SGS; mostly dry domains run much faster.

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

The unit and regression set that CI runs on every push (GitHub Actions, Linux,
pya5 0.7.3 and 1.0.0) is listed in `ci/run_tests.sh`. It covers the solver
options and their defaults, Q-centring, vector friction, limiters, SGS on
general meshes, the mesh families, pya5 compatibility, open boundaries and the
performance fast paths. Set `FLOODA5_NO_VIS=1` for test and batch runs to skip
loading the visualisers.

```bash
FLOODA5_NO_VIS=1 bash ci/run_tests.sh .                    # Git Bash / Linux / macOS
julia --project=. -t 1 test/check_determinism.jl           # compare FINGERPRINT lines with -t 4
julia --project=. -t 8 test/benchmark_baseline.jl          # step benchmark
```

Other individual tests and benchmarks:

```bash
# SGS unit tests (5-cell chain, hypsometric lookups)
julia --threads auto --project=. test/test_sgs_unit.jl

# SGS synthetic DEM validation (T0–T4: mesh, routing, mass balance)
julia --threads 1 --project=. test/synthetic_dem/test_sgs_synthetic.jl

# Edge geometry and non-orthogonal correction unit tests
julia --threads auto --project=. test/test_edge_geometry.jl
julia --threads auto --project=. test/test_noc_correction.jl

# Point-spread directional-bias benchmark (gradient correction validation)
julia --project=. test/test_point_spread.jl \
    --baseline  test/square/square_baseline.h5 \
    --corrected test/square/square_corrected.h5 \
    --frame 50 --source-lon 0.0 --source-lat 0.0

# Inflow and boundary condition tests
julia --threads auto --project=. test/test_inflow_point.jl
julia --threads auto --project=. test/test_open_boundary.jl

# Cell-vector momentum and diamond face-flux unit tests (directional-bias-reformulation)
julia --threads auto --project=. test/test_cell_momentum.jl
julia --threads auto --project=. test/test_mirror_symmetry.jl
julia --threads auto --project=. test/test_analytical_gradient_threeway.jl
```

> **Note:** `--threads 1` is the recommended default for mesh generation
> (`--meshgen`), due to a known interaction between Julia's garbage
> collector and the Python/GDAL bridge used for DEM sampling and mesh
> pre-processing. `--threads auto` mesh generation has run successfully in
> practice, but the failure mode this guards against is intermittent and
> mesh-size-dependent, so a clean run isn't a guarantee at every scale — if
> mesh generation ever crashes with a low-level fault (not a normal Julia
> error), retry with `--threads 1`. Simulation runs (`--meshload`) are
> unaffected either way and always safe with any thread count. See
> [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) §6 for the full explanation.

---

## Development Roadmap

FloodA5 is under active development. The list below reflects the areas of
future work identified through the project's testing and development history
so far — it is not a commitment or a timeline, and items may be reordered,
combined, or dropped as understanding improves.

**Done in October 2026** (details in [CHANGELOG.md](CHANGELOG.md)):
the directional bias explained and corrected (diamond face gradient + vector
friction), and made the default; a consistent, mesh-general Q-centring; SGS on
meshes with any number of neighbours; generic polygon mesh families; and a first
round of CPU performance work (3–4× faster on Carlisle, 10× for the SGS
point-source benchmark).

**Physics and correctness**
- Damping of grid-scale water-surface noise in deep channels that does not
  cost conveyance (the remaining weakness of the standard solver's defaults),
  or a symmetric form of the diamond face gradient
- A study of the timestep (Courant number) sensitivity seen on the Carlisle
  test, on both A5 and Cartesian meshes
- Validating the gradient correction for the SGS solver (vector friction and the
  local CFL already apply to it)
- Making the open-boundary (ghost-edge) flux consistent with the interior
  options (it still uses the original kernel and face-normal friction)
- Additional SGS methodology development, including further validation of
  the hydraulic-radius (R-A) flux kernel and extending sub-grid sampling to
  more general channel/culvert geometries
- An independent, minimal Riemann-solver (e.g. HLL) reference implementation
  on the same mesh, to benchmark the reduced-complexity Bates/Manning
  formulation against full shallow-water physics where the two should agree
- A like-for-like benchmark against LISFLOOD-FP on a matched Cartesian mesh
  (raster-aligned grids are available with `--meshgen-cartesian dem.tif`)
- Further mesh families: hexagonal DGGS (H3, IGEO7) and HEALPix
- Analytical benchmarks: Thacker (1981) planar-surface oscillation and a
  circular dam break (also the key demonstration of A5's rotational symmetry)

**Performance**
- Restructuring as a Julia package with precompilation (start-up and
  first-step compilation now dominate short runs), flat SGS tables, and
  non-allocating threading
- An active-cell list so that dry parts of the domain cost nothing, and
  spatially ordered cells for memory locality
- A GPU backend for the time-step kernels, from the same code

**Visualisation**
- Continued development of the CesiumJS (or a similar) web-based viewer —
  current known issues include a frame-update flash and inconsistent
  `simcomplete` signalling for live updates (see
  [docs/METHODS.md](docs/METHODS.md))

**Longer-term**
- Multi-resolution adaptive mesh refinement (Phase 3 in the project's
  internal roadmap) — several current data structures (`EdgeList` ordering,
  `adj_matrix` retained alongside the edge list, per-cell area storage) were
  deliberately designed with this in mind, but no refinement logic exists yet
- Additional boundary condition types: fixed/time-varying water-surface
  elevation (tidal) boundaries — the `.bci` `HFIX`/`HVAR` codes are already
  parsed but not yet implemented (see
  [docs/DATA_FORMATS.md](docs/DATA_FORMATS.md)) — plus spatially distributed
  rainfall, infiltration, and evaporation sources
  and additional hydrograph input formats (WaterML 2.0, CF-NetCDF)
- NetCDF/UGRID output alongside the existing HDF5 format
- Ensemble and calibration tooling (Manning's n sweeps, skill scores such as
  NSE, KGE, F2 against observed flood outlines)
- LiDAR DEM ingestion and validation for a real-world New Zealand catchment,
  as a companion to the Carlisle (UK) validation case

This project also retains a detailed internal development history —
including several investigations that did not lead anywhere directly useful,
and are recorded as such — which is available on request for anyone
continuing work in these areas.

---

## Acknowledgements

FloodA5's development has been a close collaboration between the project's
lead developer and Claude (Anthropic), used throughout as a coding and
research assistant — including implementation, debugging, code review, test
design, and the technical documentation in this repository. In the interests
of transparency, this is noted here explicitly rather than left unstated.
Scientific and engineering judgement, direction, and validation of results
against real-world data and domain knowledge remain the responsibility of
the project's human author(s).

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

If you use FloodA5 in research, please cite the accompanying paper (in preparation).
Details will be added here on publication.

---

## References

- Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). A simple inertial formulation
  of the shallow water equations for efficient two-dimensional flood inundation
  modelling. *Journal of Hydrology* 387(1–2), 33–45.
- Neal, J.C. et al. (2012). How much physical complexity is needed to model flood
  inundation? *Hydrological Processes* 26(15), 2264–2282.
- A5 DGGS: [a5geo.org](https://a5geo.org)
