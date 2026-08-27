![Status: In Development](https://img.shields.io/badge/Status-In_Development-green)
![Status: In Development](https://img.shields.io/badge/Status-Experimental-red)
![Version 0.1.0](https://img.shields.io/badge/Version-0.1.0-lightgrey)
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
  neighbours, eliminating the axis-aligned bias of a rectangular grid in principle.
  In practice, A5 pentagon edges are not perpendicular to the cell-centre-to-cell-centre
  vector (non-orthogonality angle ~16–38°), and the default Bates (2010) formulation —
  derived for an orthogonal Cartesian grid — shows a measurable residual directional
  bias as a result. **Three opt-in, experimental corrections are implemented and
  under active validation** (a non-orthogonal gradient correction, a face-local
  "diamond" reconstruction, and a cell-vector momentum representation); none is yet
  the default, and none fully eliminates the bias on its own. See
  [docs/HYDRAULICS.md](docs/HYDRAULICS.md) §7–§8 for the full picture, current
  recommendation, and open questions.
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
| [docs/HYDRAULICS.md](docs/HYDRAULICS.md) | Full flux-kernel derivations, stability limiters, and the current state of the directional-bias correction work (opt-in flags, what's validated, what isn't) |
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

**Current defaults:** the model runs with the original, uncorrected Bates
(2010) formulation unless you explicitly opt in to the directional-bias
corrections above (`--gradient-correction`, `--face-flux-method`,
`--momentum-model` — standard-flow solver only, see
[docs/USER_GUIDE.md](docs/USER_GUIDE.md) §4.3a). This is deliberate: the
corrections have real, measured benefit but are not yet validated on a
real DEM and at least one combination is known to make the bias worse
under some conditions. See [docs/HYDRAULICS.md](docs/HYDRAULICS.md) §8.

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

**`--flow-model standard`** — Diffusive wave on mean cell elevation. Faster and
useful for initial exploration or validation without a DEM.

Both solvers use the Bates et al. (2010) inertial formulation with:
- Wave-speed CFL timestep (Courant number 0.7)
- Froude limiter (Fr ≤ 0.8) and a per-edge volume limiter
- Q-centred spatial momentum smoothing (θ = 0.9, standard flow)
- Consistent momentum state (`q_prev` matches actual transferred flux, both solvers)

**Directional-bias correction (standard flow only, opt-in, off by default):**
three independent corrections for the flow-direction bias caused by A5's
non-orthogonal edge geometry are implemented and under active validation —
a WLSQ gradient correction, a face-local "diamond" reconstruction, and a
cell-vector momentum representation. None is enabled by default and none is
used by the SGS solver. See
[docs/HYDRAULICS.md §7–§8](docs/HYDRAULICS.md) for the full picture,
current best-evidenced experimental configuration, and open items, and
[docs/USER_GUIDE.md §4.3a](docs/USER_GUIDE.md) for the CLI flags.

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

**Physics and correctness**
- Further testing and correction of the flow-direction bias described in
  [docs/HYDRAULICS.md](docs/HYDRAULICS.md) §7–§8, including investigation of
  how the Froude/volume stability limiters interact with the directional-bias
  corrections (currently the leading suspect for the bias that remains even
  with correction enabled)
- Extending the directional-bias corrections (or an SGS-appropriate
  equivalent) to the SGS solver, which none of them currently touch
- Making the open-boundary (ghost-edge) flux calculation consistent with
  whichever interior flux method is selected, rather than always using the
  uncorrected kernel at domain edges
- Additional SGS methodology development, including further validation of
  the hydraulic-radius (R-A) flux kernel and extending sub-grid sampling to
  more general channel/culvert geometries
- An independent, minimal Riemann-solver (e.g. HLL) reference implementation
  on the same A5 mesh, to benchmark the reduced-complexity Bates/Manning
  formulation against full shallow-water physics on cases where the two are
  expected to agree closely
- Analytical benchmarks not yet run against the corrected flux kernels:
  Thacker (1981) planar-surface oscillation, circular dam-break (also the
  key demonstration case for A5's rotational-symmetry advantage over a
  rectangular grid)

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
- GPU-accelerated flux solver kernel (GPU acceleration currently applies only
  to SGS mesh pre-processing, not the simulation loop)
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
