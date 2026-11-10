# FloodA5 — Architecture

This document is for developers who want to understand how FloodA5 is structured,
where to find things, and how to extend the model. For the physics see
[METHODS.md](METHODS.md); for running the model see [USER_GUIDE.md](USER_GUIDE.md).

---

## 1. Overview

FloodA5 is a Julia application with a Python subprocess for grid operations. The
split is deliberate: Julia provides multithreaded performance for the simulation
loop; Python provides access to the `pya5` A5 DGGS library, which has no Julia
binding.

```
FloodModel.jl          Entry point, CLI, simulation loop, FlowState
mesh/A5Grid.jl         Julia module: mesh struct, DEM sampling, SGS tables
mesh/CartesianGrid.jl  Cartesian (4-connected) mesh backend — pure Julia,
                       no Python bridge (see mesh/CartesianGrid.jl
                       and HYDRAULICS.md §14)
mesh/DiamondFlux.jl    Diamond face-flux reconstruction (opt-in; see HYDRAULICS.md §7.2)
mesh/a5_bridge.py      Python subprocess: pya5 calls (0.x and 1.x), GeoParquet I/O
mesh/migrate_mesh_ids.py   pya5 0.x → 1.x cell-ID migration (A5_QUIRKS.md §11)
mesh/mesh_families.py  Generic polygon meshes: hex, rotated Cartesian, Voronoi, import
surfacewater/flow2d.jl Pure physics kernels (Bates, Manning R-A, CFL)
surfacewater/kernels.jl  Per-edge time-step kernels (StdOpts/SGSOpts type parameters)
surfacewater/reference_steps.jl  Frozen reference step functions (bit-equality tests)
boundaryinputs/        Source types, BC types, hydrograph I/O
visualisation/         CesiumJS and GLMakie backends, post-processing tools
```

All Julia modules are included (not imported) by `FloodModel.jl` using `include()`.
There are no registered Julia packages — the entire application lives in the
repository and is run directly with `julia FloodModel.jl`.

---

## 2. Module Responsibilities

### `FloodModel.jl`

The entry point and main application file. Responsibilities:

- **Type definitions** — `FlowState`, `EdgeList`, `InjectionPoint`, `RainPoint`,
  `AbstractSource`, `GhostEdge`, `BCType`, and the `FlowMethod` hierarchy
  (`StandardFlow`, `SGSFlow`). These are defined here rather than in their
  respective modules because `FlowState` references them and Julia requires
  types to be defined before use.
- **CLI parsing** — `_pop_flag` / `_pop_bool` helpers and the `main()` function.
  The `print_help()` function is the authoritative CLI reference.
- **Mesh initialisation** — `initialise_flow_model(mesh, method)` builds the
  `EdgeList`, assigns SGS tables, detects boundary cells, builds ghost edges,
  pre-computes the WLSQ gradient and cell-momentum weight matrices (always,
  regardless of whether the directional-bias corrections are enabled — see
  `HYDRAULICS.md` §7), and — only when `--face-flux-method diamond` is
  selected — builds the diamond flux table (`mesh/DiamondFlux.jl`), which is
  otherwise left as `nothing` since it is more expensive to construct than
  the always-built WLSQ weights.
- **Simulation loop** — `run_simulation!(state, method, sources, ...)` drives
  the timestep loop. Each outer-loop step:
  1. Applies all water sources (`apply_source!` for every `AbstractSource`)
  2. Syncs `water_depth` from `volume` (standard solver only — SGS derives
     WSE from the hypsometric table directly inside `step_sgs!`)
  3. Calls `step_standard!` or `step_sgs!`, which implement their own
     internal phases (see below)
  4. Logs progress and writes output frames

  The step functions have their own internal phase lettering, which is
  *not* the same as the outer loop's steps above:

  | Phase | `step_standard!` | `step_sgs!` |
  |---|---|---|
  | A | Parallel edge flux computation | Parallel edge flux computation |
  | B | Serial `dV` scatter | Serial `dV` scatter, per-edge donor limiter |
  | C | Parallel cell volume update | Parallel cell volume update |
  | D | Open-boundary (ghost-edge) outflow | Open-boundary (ghost-edge) outflow |
  | E | Velocity computation | Velocity computation |
  | F | Cell-vector reconstruction `_reconstruct_qvec!` (`--momentum-model cell\|blend` or `--friction-model vector`) | Same, for `--friction-model vector` only |

**Include order matters.** `flow2d.jl` is included after type definitions but
before the step functions, because the step functions call kernels defined there.
`boundaryinputs/` files are included after `InjectionPoint` and `RainPoint` are
defined (sources dispatch on these types). The include order is:

```julia
include("mesh/A5Grid.jl")
include("visualisation/VisualisationServer.jl")
include("visualisation/MakieVisualiser.jl")
# ... type definitions ...
include("surfacewater/flow2d.jl")
include("boundaryinputs/boundary_conditions.jl")
include("boundaryinputs/sources.jl")
include("boundaryinputs/timeseries_io.jl")
```

### `mesh/A5Grid.jl`

Julia module (`module A5Grid ... end`). Provides the mesh API used by
`FloodModel.jl`. Key responsibilities:

- **`A5Mesh` and `A5Cell` structs** — the in-memory mesh representation.
  `A5Mesh.cells` is a `Vector{A5Cell}`; static variables (elevation, friction)
  live in `A5Mesh.static_vars :: Dict{String, Vector{Float64}}`. Array variables
  (SGS curves) live in `A5Mesh.array_vars :: Dict{String, Matrix{Float64}}`.
- **Mesh generation** — `mesh_for_aoi(geojson, resolution)` delegates to
  `a5_bridge.py` via subprocess. Python handles all `pya5` calls and writes the
  GeoParquet file; Julia reads it back.
- **DEM sampling** — `sample_dem_mean!` and `sample_dem_centroid!` use ArchGDAL
  for GDAL-based raster access and CRS reprojection. Halton points are generated
  in `_halton(base, n, offset)`.
- **SGS pre-processing** — `build_sgs_tables!` builds hypsometric curves for
  cells and edges. The cell PIP step runs on GPU (CUDA) when available.
- **SGS runtime lookups** — `wse_from_volume`, `wetted_area_from_wse`,
  `flow_area_from_wse`, `wetted_perim_from_wse`, `hydraulic_radius_from_wse`.
  All are `@inline` for zero-overhead calls in the inner loop.
- **Cell query API** — `lonlat_to_cell`, `cell_to_lonlat`, `cell_to_boundary`,
  `cell_to_children`, `cell_to_parent`, `get_resolution`. These call `pya5` via
  PyCall with a `ReentrantLock` for thread safety.

### `mesh/a5_bridge.py`

Python subprocess called by `A5Grid.jl` for mesh generation. It:

1. Calls `pya5.fill_polygon` + `pya5.uncompact` (pya5 0.x) or
   `pya5.polygon_to_cells` + `pya5.uncompact` (pya5 1.x) to enumerate cells at
   the target resolution, covering both A5 sublattices by adding `grid_disk`
   neighbours. The cell-ID scheme is recorded in the output file's metadata
   (`A5_QUIRKS.md` §11).
2. Computes boundary polygons via `pya5.cell_to_boundary` for each cell.
3. Normalises all cell IDs to 16-character zero-padded hex and all longitudes to
   [-180, 180].
4. Writes the mesh as GeoParquet using geopandas and pyarrow.

The bridge is invoked as `python mesh/a5_bridge.py mesh_for_aoi <aoi_wkt> <res> <out>`.
A smoke-test mode is available: `python mesh/a5_bridge.py check`.

**Why subprocess and not PyCall?** PyCall can call `pya5` for single-cell queries
(with a GIL lock), but `pya5`'s batch operations interact poorly with PyCall's
garbage collector in multithreaded Julia contexts. The subprocess boundary provides
complete isolation.

### `mesh/DiamondFlux.jl`

Implements the diamond face-flux reconstruction (`--face-flux-method
diamond`; see `HYDRAULICS.md` §7.2). Deliberately does **not** add fields to
`EdgeList` — all diamond-specific geometry lives in a separate
`DiamondFluxTable`, keyed 1:1 with `EdgeList` by edge index, referenced from
`FlowState.diamond_table` (`nothing` unless `--face-flux-method diamond` is
selected). Key exports: `VertexRecord`, `eta_at_vertex`, `build_diamond_flux_table`.
When `--face-flux-method legacy` (the default) is selected, this file's
functions are loaded but never called — zero behavioural change.

### `surfacewater/flow2d.jl`

Pure physics kernels — no `FlowState`, no mesh, no I/O. All functions are
`@inline` and operate on scalar `Float64` arguments only.

Key functions:

| Function | Description |
|---|---|
| `_bates_flux` | Bates eq. 9 inertial formulation. Base flux kernel. |
| `_bates_flux_limited` | `_bates_flux` + Froude limiter + volume limiter + consistent `q_stored`. Used by `step_standard!`. |
| `_manning_flux_ra` | Manning R-A formulation (`A`, `R` based). Used by `step_sgs!`. |
| `_ghost_wse` | Ghost-cell WSE for each BC type. |
| `_bates_ghost_flux` | Outflow flux across a ghost edge (standard solver). |
| `_manning_ghost_flux` | Outflow flux across a ghost edge (SGS solver). |
| `_cfl_dt` | Wave-speed CFL timestep. |
| `_q_centred` | Q-centred (spatially smoothed) unit discharge. |
| `_bates_flux_corrected` / `_bates_flux_limited_corrected` | WLSQ or diamond-corrected driving head, standard flow. See `HYDRAULICS.md` §7. |

Constants: `_G = 9.81`, `HFLOW_THRESHOLD = 0.001`, `FROUDE_LIMIT = 0.8`.
`N_SIDES = 5` and `DONOR_EDGE_DIVISOR = 10` (`= 2 × N_SIDES`) still exist in
`FloodModel.jl`, but only as **A5-specific legacy defaults** kept for the
existing A5 test fixtures — production code no longer reads either
directly. The live, mesh-derived equivalents are `state.max_neighbours`
(5 on A5, 4 on the Cartesian backend, derived from the loaded mesh's own
adjacency) and `_donor_edge_divisor(state) = 2.0 * state.max_neighbours`,
and the Fix B volume limiter in `_bates_flux_limited`/
`_bates_flux_limited_corrected` takes a `vol_limit_divisor` keyword
(default `5.0`, i.e. unchanged A5 behaviour for any caller that doesn't
pass it) which `step_standard!` sets to `Float64(state.max_neighbours)`.
See the `FlowState.max_neighbours` comment in the struct listing below, and
the project's internal development history (available on request), for the full reasoning. `Q_CENTRE_THETA_DEFAULT = 0.9` is
only the *default* — the live value is `state.q_centre_theta`, set via
`--q-centre-theta` and threaded through explicitly (it was promoted from a
plain constant to a `FlowState` field so it could be varied per run for
the directional-bias correction work — see `HYDRAULICS.md` §7).

This file is designed to be includable independently (e.g. by a GPU kernel
wrapper or a 3D coupling plugin) without pulling in the full application.

### `boundaryinputs/`

Three files, included after type definitions in `FloodModel.jl`:

**`sources.jl`** — `AbstractSource` abstract type and `InflowPoint` concrete type.
Dispatch methods `apply_source!(state, src, t, dt)` and `cumulative_volume(src, t)`
for all source types. `InjectionPoint` and `RainPoint` are defined in `FloodModel.jl`
as `AbstractSource` subtypes; `InflowPoint` is defined here.

**`boundary_conditions.jl`** — `BCType` enum, `GhostEdge` struct, `BoundarySegment`
struct, `_build_ghost_edges(mesh, edges, ...)`, `_ghost_wse`, `_bates_ghost_flux`,
`_manning_ghost_flux`, `load_bc_file` (GeoJSON reader), `apply_bci_free_entries!`.

**`timeseries_io.jl`** — `AbstractTimeSeriesReader`, `LisfloodBDYReader`,
`TwoColumnCSVReader`, `read_timeseries` dispatch, `parse_bci_file`.

### `io/FloodA5Output.jl`

Module `FloodA5Output` (aliased `FO` in `FloodModel.jl`): the HDF5 output
format, writer and reader in one place (0.4.0; layout in `DATA_FORMATS.md`
§2). It depends only on HDF5, JSON3, SHA and Dates, so analysis scripts and
notebooks can load it without `FloodModel.jl`.

- **Writer** — `create_output` (run record, mesh, empty `/results`),
  `append_frame!`, `write_summary!`, `provenance` (git, environment).
  `FloodModel.jl` calls them from `_write_mesh_metadata!` / `_write_frame!`,
  with arrays already in the mesh file's cell order; `_run_mesh_source` and
  `_run_inputs` decide which mesh file and input files a run records.
- **Reader** — `n_frames`, `frame_times`, `frame_at`, `frame_vars`,
  `read_frame`, `read_var`, `read_series`, `cell_table`, `cell_area`,
  `cell_polygons`, `cell_neighbours`, `run_info`; each takes a path or an
  open `HDF5.File` and reads format 2 and format 1 (≤ 0.3.x) files. Scripts
  should read output through these, not by dataset path.
- **Reproduction** — `extract_mesh`, `rerun_argv` (used by `--rerun`). The
  CLI flags `rerun_argv` must drop or rewrite are listed in the module
  (`_DROP_FLAGS`, `_BOOL_FLAGS`): a new mesh-generation or file-valued flag
  in `main()` needs an entry there.

### `analysis/` — `FloodA5Analysis` (`FA`)

Loading runs and computing metrics, for scripts and notebooks. Pure Julia:
it does not include `FloodModel.jl` and needs no Python (loads in ~3 s). It
reads output only through `FO`. Format-2 files carry the mesh geometry; for
format-1 files, pass the mesh's GeoParquet (`load_run(p; parquet = …)`).

- `runs.jl`: `Run` and `MeshGeom` (cells in the output file's order:
  ids, centres, elevation, area, polygons, neighbours, scalar statics),
  `load_run`, `frame_at`, `matched_frames`, `field`/`depth`/`speed`/`wse`
  (SGS-aware), `storage`, `local_xy`, `source_cell`, `check_frame`,
  `print_frames`.
- `parquet_mesh.jl`: `load_mesh_parquet` for format-1 runs, via ArchGDAL's
  Parquet driver (Parquet2 cannot read the list columns), and
  `polygon_area_m2`. The result matches what the model uses: row order, closed
  rings, area `sgs_cell_area` or else polygon area, neighbour slot order,
  z = 0 without a DEM. GDAL's Arrow reader deadlocks in Julia unless
  `GDAL_NUM_THREADS=1`, which the loader sets while it reads.
- `compare.jl`: run vs reference (storage, rms Δh, CSI, bias, deep-water
  speed) and twin-run divergence with its e-folding time.
- `shape.jl`: point-spread and planar-slope metrics (harmonics, sector
  ratios, hull circularity, front radius by bearing, N/S split, deflection,
  flow direction).
- `noise.jl`: checkerboard and velocity-noise metrics.
- `channel.jl`: straight-channel series metrics.

The analysis scripts in `test/` (`wp1_metrics.jl`, `diag_vs_reference.jl`,
`diag_divergence.jl`, `wp2_checkerboard.jl`, `channel_bench/analyse_channel.jl`,
`test_point_spread.jl`, `test_octagonal_anisotropy.jl`, `test_planar_slope.jl`,
`test_planar_symmetry.jl`) are thin command-line wrappers over `FA`. Their
output is identical to the pre-module versions. Each keeps its own historical
projection constants (Earth radius, metres per degree, origin snapping),
passed to `FA` explicitly.

### `visualisation/`

Four scripts and one static directory:

| File | Description |
|---|---|
| `MakieVisualiser.jl` | GLMakie native desktop window. Module `MakieVisualiser`. |
| `VisualisationServer.jl` | Oxygen.jl HTTP + WebSocket server. Module `VisualisationServer`. |
| `view_h5_output.jl` | Standalone post-processing viewer for HDF5 output. |
| `visualise_mesh.jl` | Publication-quality mesh and SGS figures (wire-frame, cross-section, hypsometric curve). |
| `cesium/` | CesiumJS static files (`index.html`, `config.example.json`). |

---

## 3. Key Data Structures

### `FlowState`

The central mutable struct, passed through the entire simulation loop. The
fields below are grouped by purpose; see the struct definition in
`FloodModel.jl` (it carries extensive inline comments explaining the
non-obvious ones, particularly around the directional-bias fields) for the
authoritative, fully-commented version.

```julia
mutable struct FlowState
    # ── Core per-cell state ──────────────────────────────────────────────
    cell_ids    :: Vector{String}
    water_depth :: Vector{Float64}   # m — diagnostic; derived from volume each step
    volume      :: Vector{Float64}   # m³ — PRIMARY STATE VARIABLE
    velocity    :: Vector{Float64}   # m/s scalar magnitude
    vel_u       :: Vector{Float64}   # m/s eastward velocity component
    vel_v       :: Vector{Float64}   # m/s northward velocity component
    elevation   :: Vector{Float64}
    manning_n   :: Vector{Float64}
    cell_area   :: Vector{Float64}
    cell_lons   :: Vector{Float64}
    cell_lats   :: Vector{Float64}
    adjacency   :: Dict{String, Vector{String}}
    adj_matrix  :: Matrix{Int}       # (max_nb × n_cells)
    edges       :: EdgeList
    sgs_tables  :: Vector{Any}       # Vector{SGSTable} for SGS; empty for standard

    # ── Open boundary / ghost-edge state ─────────────────────────────────
    boundary_mask :: BitVector
    ghost_edges   :: Vector{Any}     # Vector{GhostEdge}
    ghost_cell_bc :: Vector{Any}     # Vector{BCType}
    vol_removed   :: Float64         # cumulative ghost-edge outflow (m³)

    # ── Directional-bias correction state (see HYDRAULICS.md §7) ────────
    # Always allocated and populated regardless of whether corrections are
    # enabled — the cost of building them is small relative to mesh
    # generation, and it keeps FlowState's shape independent of CLI flags.
    grad_wse                  :: Matrix{Float64}   # (2 × n_cells) per-cell WSE gradient
    wlsq_weights              :: Matrix{Float64}   # (2×max_neighbours × n_cells) pre-computed WLSQ projection
    q_centre_theta            :: Float64           # --q-centre-theta
    gradient_correction       :: Bool              # --gradient-correction
    gradient_correction_alpha :: Float64           # --gradient-correction-alpha (legacy method only)
    qvec_u, qvec_v            :: Vector{Float64}   # per-cell 2D discharge vector (--momentum-model cell)
    cell_edge_index           :: Matrix{Int}       # (slot, cell) → edge index
    mom_weights               :: Matrix{Float64}   # (2×max_neighbours × n_cells) WLSQ projection for qvec
    momentum_model            :: Symbol            # :edge (default) or :cell
    face_flux_method          :: Symbol            # :legacy (default) or :diamond
    diamond_table             :: Any               # DiamondFluxTable, or nothing if :legacy

    # ── Mesh-backend generalisation (cartesian-formulation branch) ──────
    # max_neighbours: the mesh's own neighbour-count ceiling — 5 for A5,
    # 4 for the Cartesian backend (mesh/CartesianGrid.jl), derived once in
    # initialise_flow_model from the loaded mesh's actual adjacency, NOT a
    # hardcoded constant. Sizes adj_matrix/wlsq_weights/mom_weights/
    # cell_edge_index above (hence the "2×max_neighbours" comments) and
    # sets the ×2-safety-margin donor/volume-limiter divisors (see
    # `_donor_edge_divisor` and `vol_limit_divisor` in flow2d.jl). A cell
    # with fewer actual neighbours than this ceiling (a boundary cell on
    # either backend) simply leaves the extra matrix slots zero-padded —
    # the ACTUAL per-cell neighbour count is always derivable on demand via
    # `count(!iszero, view(adj_matrix, :, i))`, so no separate per-cell
    # array is stored. `N_SIDES = 5` still exists as a module constant, but
    # is now used only by A5-specific test fixtures — production code reads
    # this field instead.
    max_neighbours :: Int

    # ── WP1/WP2 (2026-10; see HYDRAULICS.md §15) ─────────────────────────
    # Fields 33–36. Backward-compatible 31-, 32- and 33-argument outer
    # constructors default them to the pre-WP1 behaviour.
    friction_model :: Symbol           # :component (default) or :vector
    limiter        :: Symbol           # :legacy (default), :positivity or :hybrid
    cfl_mode       :: Symbol           # :legacy (default) or :local
    cell_length    :: Vector{Float64}  # ℓ_i = 4A_i/P_i, for cfl_mode :local
end
```

`momentum_model` may also be `:blend`, and `face_flux_method` `:diamond_wetdry`.

`A5Mesh` (in `mesh/A5Grid.jl`) has an `id_scheme` field since WP-M1: `"a5-v0"`
(pya5 0.x cell IDs), `"a5-v1"` (pya5 1.x), `"cartesian"` or `"polygon"`. It is
read from and written to the GeoParquet metadata, and `load_mesh_geoparquet`
refuses an A5 mesh whose scheme differs from the installed pya5's.

### `EdgeList`

Undirected edge list — each edge stored once, with `cell_i < cell_j`:

```julia
struct EdgeList
    n_edges     :: Int
    cell_i      :: Vector{Int}       # lower-index cell
    cell_j      :: Vector{Int}       # higher-index cell
    width       :: Vector{Float64}   # shared edge length (m)
    L           :: Vector{Float64}   # centre-to-centre haversine distance (m)
    cos_theta   :: Vector{Float64}   # oriented cos θ = d̂·n̂ ∈ [0,1] (n̂ oriented so d̂·n̂ ≥ 0)
    sill        :: Vector{Float64}   # sill elevation (m)
    flux        :: Vector{Float64}   # q (m²/s) at t-dt — standard flow, SGS Bates kernel
    flux_Q      :: Vector{Float64}   # Q (m³/s) at t-dt — SGS R-A kernel only, zero otherwise
    collinear_i :: Vector{Int}       # Q-centred scheme: most-collinear edge on the ci side
    collinear_j :: Vector{Int}       # Q-centred scheme: most-collinear edge on the cj side
    skew_x      :: Vector{Float64}   # V̂ₓ — tangential correction vector component (see below)
    skew_y      :: Vector{Float64}   # V̂ᵧ
    dx_m        :: Vector{Float64}   # cj−ci displacement, local metres, x
    dy_m        :: Vector{Float64}   # cj−ci displacement, local metres, y
    nf_x        :: Vector{Float64}   # oriented unit face normal, x
    nf_y        :: Vector{Float64}   # oriented unit face normal, y
end
```

**On `skew_x`/`skew_y`:** despite the field names (kept unchanged through a
revision to minimise diff size), these hold the *directional* correction
vector `V̂ = n̂ − c·d̂` (dimensionless, `|V̂| = sin θ`), **not** a positional
offset in metres. An earlier formulation genuinely did store a positional
offset here and was found to be geometrically incorrect — see the
`_edge_geometry` docstring in `mesh/A5Grid.jl` for the full account if
you're modifying anything that touches this field. This is exactly the
kind of thing worth reading before touching this code, not after.

**Sign convention:** `flux > 0` means flow from `cell_j` to `cell_i` (j is
higher). When `WSE_i > WSE_j`: `dWSE > 0` → `q_new < 0` → `Q < 0` →
`dV[i] += Q·dt` (i loses volume) — correct.

### `SGSTable`

Pre-computed hypsometric curves for one cell:

```julia
struct SGSTable
    elev_bins        :: Vector{Float64}   # WSE knots (n_bins)
    vol_curve        :: Vector{Float64}   # cumulative volume (m³)
    area_curve       :: Vector{Float64}   # wetted plan area (m²)
    cell_area        :: Float64
    z_min            :: Float64
    z_max            :: Float64
    edge_area_curves :: Matrix{Float64}   # (n_bins × 5) flow area per slot (m²)
    edge_perim_curves:: Matrix{Float64}   # (n_bins × 5) wetted perimeter per slot (m)
end
```

---

## 4. Extension Points

### Adding a new water source type

1. Define a new struct `<: AbstractSource` in `boundaryinputs/sources.jl`.
2. Implement `apply_source!(state, src::MySource, t, dt)` — add volume to
   `state.volume[src.cell_index]`.
3. Implement `cumulative_volume(src::MySource, t)` — return the cumulative volume
   added from t=0 to t (for mass balance logging).
4. Add a CLI flag in `FloodModel.jl`'s `main()` and construct the source there.

The simulation loop calls `apply_source!` for every element of `all_sources ::
Vector{AbstractSource}` each step, so no further changes are needed in the loop.

### Adding a new boundary condition type

1. Add a new variant to the `BCType` enum in `boundaryinputs/boundary_conditions.jl`.
2. Add a branch in `_ghost_wse(wse_ci, sill, bc, ...)` to return the appropriate
   ghost-cell WSE.
3. Add the string name to `load_bc_file` and `apply_bci_free_entries!` for CLI/file
   recognition.
4. The flux kernel (`_bates_ghost_flux` / `_manning_ghost_flux`) already passes
   `wse_ghost` from `_ghost_wse` — no change needed there unless the new BC type
   requires a fundamentally different flux computation.

`FixedWSE` (tidal BC) and `FixedQ` stubs are already present in `BCType`; their
`_ghost_wse` branches return `wse_series_val` and are awaiting implementation of the
time-series lookup.

### Adding a new physics kernel

`surfacewater/flow2d.jl` is designed to hold all scalar physics. New kernels
(infiltration, evaporation, a Riemann solver option) should be added here as
`@inline` functions taking only `Float64` scalar arguments. The step functions
in `FloodModel.jl` call these in their inner loops.

A companion file (e.g. `surfacewater/infiltration.jl`) is preferable for larger
additions, following the same pattern: pure scalar kernels, no `FlowState`
dependency, included by `FloodModel.jl`.

A pass over cells or edges that uses such a kernel must run on every compute
backend: write it as a kernel body launched by `foreach_index!` (§6a).

### Adding a new time-series format

1. Define a new struct `<: AbstractTimeSeriesReader` in `boundaryinputs/timeseries_io.jl`.
2. Implement `read_timeseries(reader::MyReader) →
   Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}`, where the dict maps
   series name → `(t_seconds, Q_m3s)`.
3. The `InflowPoint` constructor and `parse_bci_file` call `read_timeseries` through
   dispatch — no other changes are needed.

Planned future readers: `WaterML2Reader`, `CFNetCDFGaugeReader`, `HydroMLReader`.

### Adding a new directional-bias correction

The `--gradient-correction` / `--face-flux-method` / `--momentum-model`
machinery (see `HYDRAULICS.md` §7) is deliberately structured so a new
correction can be added without touching the others:

1. A new driving-head construction (like the diamond method) is a new
   function taking pre-computed geometry and cell/edge state and returning
   `dWSE_n`, dispatched on a new `face_flux_method` symbol in
   `step_standard!` Phase A. It does not need a new flux kernel — B1's
   proof (see `FloodA5_PhaseB_B1_FaceFluxEquation.md` in the project's
   internal development history) established that any exact gradient
   source plugs into the existing `_bates_flux_limited_corrected`
   unchanged.
2. A new momentum representation (like the cell-vector model) is a new
   `momentum_model` symbol, a new per-cell or per-edge state array on
   `FlowState`, and a dispatch branch computing `q_prev_eff` in
   `step_standard!` Phase A plus a reconstruction step in Phase F.
3. **Either kind of correction currently only reaches `step_standard!`.**
   Extending either to `step_sgs!` is open work — see `HYDRAULICS.md` §9.4
   for the specific known risk (interaction with the SGS dry-cell WSE
   clamp) that needs resolving first.

### Adding a new visualisation variable

1. Add the variable to `FlowState` and compute it in the step function.
2. Add the dataset to `_write_frame!` in `FloodModel.jl`.
3. Add the variable to `VAR_DEFS` in `visualisation/cesium/index.html`:

```javascript
const VAR_DEFS = {
    depth:      { colormap: "turbo",  fixedMax: null },
    saturation: { colormap: "blues",  fixedMax: 1.0  },
    // add new variable here
    myvar:      { colormap: "plasma", fixedMax: null },
};
```

Unknown variables sent by the server are handled by a fallback in the CesiumJS
viewer, so new variables will appear in the dropdown automatically. Explicit
entries in `VAR_DEFS` allow you to set a specific colormap and fixed maximum.

---

## 5. Python Bridge Details

### A5 cell ID conventions

- pya5 returns cell IDs as Python `int`. The bridge converts to hex via
  `hex(cell_id).lstrip("0x").zfill(16)` — always 16 characters, zero-padded.
- Julia normalises all IDs to 16-char hex on load: `_norm_id(id) = A5Grid._to_hex(parse(UInt64, id, base=16))`.
- Never compare cell IDs without normalising both sides to 16-char hex.

### Dual sublattice

`pya5.fill_polygon + uncompact` returns cells from only one A5 sublattice. The
bridge adds neighbours of these primary cells (via `grid_disk(cell, 1)`) and
includes those whose centres fall inside the AOI — capturing both sublattices.
Any mesh generated without this step will have zero edges and produce no flux.

### `grid_disk` returns compact results

`pya5.grid_disk(cell, 1)` returns a compact result — a mix of cells at different
resolution levels. Always follow with `pya5.uncompact(disk, target_resolution)`
(pya5 0.x) or `pya5.uncompact(disk)` and drop compaction markers (pya5 1.x) — see
`A5_QUIRKS.md` §11 and `A5Grid.grid_disk_neighbours`.

---

## 6. Thread Safety

- **Mesh generation** — `--threads 1` is the recommended default. The known
  hazard: `ArchGDAL.createcoordtrans` (a PyCall-backed GDAL call) is invoked
  from inside `Threads.@threads` loops in both `sample_dem_mean!` (standard
  DEM sampling) and `build_sgs_tables!`'s Step 1 (hypsometric curve build).
  Garbage collection firing on a non-main thread during one of these calls
  can trigger a `PyCall` finaliser crash (`EXCEPTION_ACCESS_VIOLATION` or
  similar). Two specific instances of this were found and fixed by making
  the affected loop serial (`_shared_edge`'s internal implementation, and
  `build_sgs_tables!`'s Step 2 edge-sill loop) — both are now safe under any
  thread count — but the two loops named above still carry the same general
  hazard and have not been made serial. In practice, `--threads auto` mesh
  generation (including DEM sampling) has run successfully; the failure mode
  was originally observed to depend on mesh size and GC timing rather than
  occurring every time, so a clean run is not a guarantee. If mesh generation
  ever crashes with a low-level fault rather than a Julia error, retry with
  `--threads 1`.
- **Simulation loop** — safe with any thread count, and (since WP0, 2026-10)
  deterministic across thread counts. Since WP4 stage 2 every phase is either
  a per-edge kernel (`surfacewater/kernels.jl`) that writes only its own
  edge's slots, or a per-cell gather over a cell→edge CSR in edge order (the
  volume update, positivity totals, velocity, cell vectors) — no shared
  writes, no atomics, the same summation order as a serial loop. Q-centring
  reads a start-of-step snapshot of the momentum array. Loops use Polyester
  `@batch per=core` (non-allocating; `per=core` avoids spinning on every
  logical core). Kernel helpers must accept `AbstractVector`, since
  Polyester passes captured arrays as `PtrArray`s.
  `test/check_determinism.jl` verifies identical results for 1 vs N threads,
  and `test/test_wp4_stage2.jl` checks the steps against the frozen
  `step_*_reference!` copies.
- **Active edges** (WP4 stage 3) — Phase A visits only the edges that can
  change: those of cells holding water, last step's set (so an edge that
  has just dried is visited once more and zeroed by the kernel itself), and
  edges that can carry flux while dry (none on normally built meshes). A full
  pass runs on the first step, when elevations/sills/areas change, and when
  more than half the cells hold water. The volume gather still covers every
  edge, so mass balance is exactly that of a full pass. See `_active_edges!`
  and `test/test_wp4_stage2.jl` R4.
- **Cell order** — `run_flood_model` permutes the mesh into a spatial
  (Morton) order before building the flow model (`--cell-order`); output is
  mapped back to the file's order (`SimOutput.perm`). Edge orientation
  (`cell_i < cell_j`) therefore follows the internal order, which is why the
  order is not changed under legacy (orientation-dependent) Q-centring.
- **Ghost-edge Phase E** — serial (ghost edges are owned by single cells).
- **Cell queries** (`lonlat_to_cell` etc.) — thread-safe via `ReentrantLock(_pycall_lock)`.
- **Compute backends** (WP5) — see §6a.

## 6a. Compute Backends (CPU threads and GPUs)

Every solver (`standard`, `sgs`, `swe`) runs on the CPU or on an NVIDIA GPU
from one code base (`--backend cpu|gpu`; `run_flood_model(...; backend =
:cpu | :gpu | :ka_cpu)`).

**Kernel bodies and launchers** (`surfacewater/backend.jl`). Each per-step
pass is an `@inline` body `body(i, args...)` for one cell, face, boundary
face or mesh vertex, writing only the slots that index owns. `args` is a
`NamedTuple` of arrays and scalars. `foreach_index!(backend, n, body,
args...)` launches it:

| Backend | What runs | Used for |
|---|---|---|
| `CPUThreads()` | Polyester `@batch per=core` | the default; bit-identical to the loops it replaced |
| `KADevice(CUDABackend())` | a KernelAbstractions kernel on the GPU | `--backend gpu` |
| `KADevice(KernelAbstractions.CPU())` | the device code path on the CPU | `--backend ka-cpu`; tests without a GPU |

`on_host_order!` runs a short loop whose order matters (point sources, open
boundaries) on one device thread. Reductions (`_minimum`, `_maximum`,
`_dot_sum`) keep the serial order on host arrays and use the device's own
reductions on GPU arrays. CUDA.jl is loaded only for `--backend gpu`
(`_cuda_backend()`), and the GPU run is entered through `Base.invokelatest`
so that it sees CUDA.jl's methods; CPU runs and CI need neither CUDA nor a
GPU.

**Data on the device.** `_swe_args`, `_std_device_args` and
`_sgs_device_args` gather everything a solver's passes touch into one
NamedTuple; `to_device` adapts it (CuArrays on a GPU; the same Arrays on the
CPU backends). Structures a GPU cannot hold are converted first: the edge
list and diamond table become NamedTuples of arrays (`_edge_view`,
`_dtab_view`), `BitVector`s become `Vector{Bool}`, ghost edges become
struct-of-arrays, and the per-cell SGS tables are packed into padded
arrays, read through `_tab(tables, i)` (the `SGSTable` on the CPU, a
NamedTuple of views on a device), so the A5Grid lookups run unchanged.
Never put the same array into a device NamedTuple under two names: `Adapt`
would make two device copies. Alias after the transfer instead (see
`_sgs_step!`).

**The run loop on a device** (`surfacewater/device_run.jl`). `DeviceRun{K}`
(K = `:swe`, `:std`, `:sgs`) holds the device arrays for the whole run. Each
step computes Δt, rainfall, point sources, the depth sync and the step on
the device. The host `FlowState` is brought up to date
(`_device_sync!`) only where it is read: output frames, the visualiser, the
every-50-steps log line, the debug lines and the end of the run.

**Supported on a device:** the options of `--scheme recommended` and `v0.2` (and
`legacy` with `--cfl-mode local`): edge momentum, the local CFL, any
Q-centring mode and face gradient; η filter, momentum damping,
positivity at boundaries; all source types; open, closed and critical
boundaries. Not on a device (refused with an error): `--momentum-model
cell|blend`, `--cfl-mode legacy` (a sort), `--limiter-stats`.

**Order-dependent passes.** Point sources are added on one device thread,
in source order. Open-boundary outflow (Phase D) on the CPU walks the ghost
edges in order and reads a neighbour's volume after earlier ghost edges have
changed it, so the order is part of the result. The device runs it in
dependency levels (`_ghost_schedule`): each edge one level after the latest
earlier edge that touched its cell or its interior neighbour, the levels in
order, the edges of a level in parallel, and the outflow volumes summed in
edge order — which reproduces the CPU loop exactly. With
`--open-boundary normal` (inward neighbours) there are a handful of levels
(4 on a 360 × 360 grid); with `legacy` the neighbours run along the boundary
and the levels approach one per edge.

**Throughput** (`test/benchmark_backend.jl`: all cells wet, Cartesian n × n,
million cell-steps per second; Xeon E-2278G, 8 cores / 16 threads, against
an RTX 4000 SFF Ada, whose FP64 rate is 1/64 of FP32):

| cells | standard CPU | standard GPU | SWE CPU | SWE GPU |
|---|---|---|---|---|
| 65 536 | 10.0 | 7.2 | 3.1 | 4.8 |
| 262 144 | 9.1 | 26.0 | 2.4 | 22.0 |
| 1 048 576 | 8.4 | 93.2 | 2.7 | 29.4 |

The GPU overtakes the CPU at about 10⁵ cells; below that, launch latency
(about 20 launches per standard step with the default options) dominates. Dry cells cost nothing
on the CPU (active-edge set, WP4) but are still visited on the GPU, so a
mostly dry domain moves the crossover upward.

**Precision and agreement.** Float64 throughout. The device path run on the
CPU is bit-identical to `CPUThreads` for all three solvers
(`test/test_backend.jl`, in CI). On the GPU, single steps agree to round-off
(10⁻¹⁵ relative); over a run the local-inertial solvers amplify round-off
at wet/dry thresholds (depth rms 10⁻⁸–10⁻⁶ m on a rained plane over 1 h,
shrinking with time), as they do between thread counts with a different
cell order (HYDRAULICS.md §15.5). Carlisle res 18, 6 h (`--cell-order
file`): on the Cartesian 20 m mesh the GPU run matches the CPU to 3·10⁻¹¹ m
in every frame; on A5 (standard and SGS) the two agree to 10⁻⁹ m for 30 min,
then diverge locally by up to 1–2 m — exactly as a CPU run in Morton cell
order diverges from one in file order, in the same cells and at the same
times — while total storage stays within 3·10⁻⁴ (CPU reordering: 5·10⁻⁴).
At this size (30 000 cells) the GPU is slower than 16 CPU threads for
standard flow (153 vs 104 s) and faster for SGS (127 vs 281 s). FP64 throughput varies a lot between GPUs
(1/64 of FP32 on most workstation and consumer cards); a Float32 option is
future work.

**Adding a pass** (or a phase to an existing solver): write the body in the
solver's `*_device.jl` (or `kernels.jl`), make the CPU function call it
through `foreach_index!(CPUThreads(), ...)`, add its arrays to the solver's
device-args builder, and call it from the device driver in the same place as
in the CPU step. Bodies must not allocate, call `@info`/`error` with
interpolated strings, read global `Ref`s, or dispatch on runtime values (use
`Val` or type parameters). `test/test_backend.jl` must stay bit-identical
on `KADevice(CPU())`, and the frozen-reference test must still pass.

---

## 7. Testing

Tests live in `test/`. Run individual test files with:

```bash
julia --threads auto --project=. test/test_sgs_unit.jl
```

Key test files:

| File | What it tests |
|---|---|
| `test_a5grid.jl` | A5Grid cell API |
| `test_edge_geometry.jl` | `_edge_geometry` (cos θ, tangential correction vector, face normal) |
| `test_sgs_unit.jl` | SGS hypsometric lookups, 5-cell chain routing |
| `test_sgs_edge_ra_tables.jl` | SGS edge flow area / wetted perimeter curves |
| `test_inflow_point.jl` | `InflowPoint`, `_interp_hydrograph`, `.bdy` reader, `.bci` parser |
| `test_open_boundary.jl` | Ghost edge geometry, BC types, outflow flux, mass balance |
| `test_output_format.jl` | HDF5 output format 2: run record, embedded mesh and inputs, geometry, `--rerun` bit-identical (O1–O4), format-1 reader (O5) |
| `bench_output_format.jl` | Format 1 vs format 2: write time, size, read times on an existing run (not in CI) |
| `test_analysis.jl` | `analysis/` (FA): one format-2 fixture run, then `test/analysis/test_*.jl` — core loading, parquet mesh = format-2 mesh, compare/divergence, shape, noise, channel metrics on synthetic fields with known answers (`FA_TEST_ONLY=shape,…` runs a subset) |
| `test_flat_rainpoint.jl` | Flat-terrain point-source integration test |
| `test_noc_correction.jl` | WLSQ gradient weights and reconstruction correctness (T-NOC1–5) |
| `test_point_spread.jl` | Point-spread directional-bias benchmark (convex hull, Polsby–Popper circularity) |
| `test_planar_symmetry.jl` | Planar-slope north/south volume-asymmetry benchmark — the primary directional-bias acceptance test |
| `test_mirror_symmetry.jl` | Synthetic mirror-symmetric geometry — isolates the driving-head formula from real-mesh noise |
| `test_cell_momentum.jl` | Cell-vector momentum reconstruction (`--momentum-model cell`) |
| `test_analytical_gradient_threeway.jl` | Gradient-reconstruction exactness, raw/limited flux comparison, `:edge`-vs-`:cell` convergence |
| `synthetic_dem/test_sgs_synthetic.jl` | Full SGS validation: mesh, routing, sub-cell notch, mass balance (T0–T4) |

The synthetic DEM test (`test_sgs_synthetic.jl`) is the most comprehensive
validation: it uses a purpose-built DEM with a Gaussian embankment containing
a sub-cell notch, and confirms that the SGS solver routes 5–14× more water through
the notch than the standard solver at the same injection history.

Several one-off audit/diagnostic scripts (`audit_vertex_valence.jl`,
`audit_diamond_gradient.jl`, `audit_centre_vs_centroid.jl`,
`audit_vertex_winding.jl`, `diagnose_skew_bias.jl`, and others) were written
during the directional-bias investigation to check specific geometric
properties against real mesh data rather than as ongoing regression tests.
They're kept in `test/` since they're cheap to re-run and useful if the
mesh-generation pipeline changes, but aren't part of the standard test
suite above.

Three further diagnostic scripts were added on the `cartesian-formulation`
branch, specifically for validating the Cartesian mesh backend
(`mesh/CartesianGrid.jl`) and investigating a mesh-independent directional
bias found while doing so — see `HYDRAULICS.md` §14 for the full story,
including a resolution-refinement test whose result (the bias got *worse*,
not better, at finer resolution) is still an open question:

| File | What it tests |
|---|---|
| `audit_cartesian_geometry.jl` | Adjacency correctness (count/symmetry/direction/distance, from `mesh.adjacency` directly) and edge-geometry orthogonality (`cos_theta`/`skew` via the real `_build_edge_list`) for any generated mesh — confirmed exact on the Cartesian backend. |
| `test_octagonal_anisotropy.jl` | Cardinal-vs-diagonal directional bias, via sectors centred exactly on the 8 compass bearings (fixes `test_point_spread.jl`'s sectors, which are centred *between* cardinal and diagonal bearings) plus a fixed-radius exact-bearing sample. The tool that found the open anisotropy result above. |
| `check_boundary_proximity.jl` | Reports how close a frame's wetted region is to each of the four domain edges (N/S/E/W separately) — a sanity check for whether a closed-boundary run's results might be contaminated by the front reaching a wall. |
