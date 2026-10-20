# Changelog

All notable changes to FloodA5 are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning follows
[Semantic Versioning](https://semver.org) — see `VERSIONING.md` for the
project's specific policy and the reasoning behind starting at 0.1.0.

## Historical development (pre-0.1.0)

FloodA5 was developed over many sessions before formal versioning began.
This is a brief prose summary of that trajectory, not a version-by-version
log — no version numbers are retrofitted to this period (see
`VERSIONING.md`). Full detail is retained in the project's internal
development history.

Development proceeded roughly through these phases: core A5 mesh
generation and the Bates (2010) inertial flow solver on a flat/synthetic
domain; sub-grid sampling (SGS) using hypsometric volume curves derived
from LiDAR, added to resolve sub-cell channel features the standard
solver's cell-mean elevation cannot see; a stability investigation
(prompted by checkerboard oscillation observed on the Carlisle domain)
that identified and fixed an inconsistent momentum-storage pattern and
added Froude/volume flux limiters, informed by a line-by-line comparison
against LISFLOOD-FP and CAESAR-Lisflood reference implementations; an
alternative hydraulic-radius (R-A) flux formulation for SGS edges; dynamic
fluvial inflows and open/closed boundary conditions with LISFLOOD-FP
`.bci`/`.bdy` format support; and — the most extensive single line of
work — a multi-session investigation into a systematic directional flow
bias on the A5 pentagonal mesh, producing three independent, opt-in
correction candidates (a WLSQ gradient correction, a face-local "diamond"
reconstruction, and a cell-vector momentum representation), none of which
fully resolves the bias on its own.

## [Unreleased]

Merged to `main` on 2026-10-20 via `wp0-correctness`, which is based on
`cartesian-formulation`. No version number assigned yet; under
`VERSIONING.md`'s policy the new CLI flags make this a MINOR bump (likely
0.2.0), to be decided by the maintainer.

### WP3 (part 1): mesh families and raster-aligned Cartesian grids (`wp3-mesh-families`)

Work package WP3 from `FloodA5_Review_and_ForwardPlan_2026-10.md` (§6.1, §9),
the parts that do not touch the solver: SGS generalisation (per-edge-side
curves) is left until WP1/WP2 are merged.

#### Added
- `mesh/mesh_families.py`: FloodA5 GeoParquet meshes of other families —
  regular hexagons, rotated Cartesian, Voronoi (jittered, optionally Lloyd-
  relaxed to centroidal; generator or centroid centres) — and an importer for
  any polygon GeoParquet/GeoJSON/GeoPackage (reprojected to EPSG:4326).
  Edge-based adjacency; ID scheme `polygon`.
- `generate_cartesian_mesh_raster` (`mesh/CartesianGrid.jl`) and
  `--meshgen-cartesian FILE.tif|.asc`: a Cartesian grid aligned cell-for-cell
  with a raster in its own (projected) CRS, with the raster as the bed
  elevation — for LISFLOOD-FP comparisons. Note: FloodA5 measures cell areas
  on its spherical metric, 0.38% below the British National Grid pixel area
  near Carlisle (uniform).
- `test/test_mesh_families.jl` (in CI): Tier 0 on every family — closure
  Σ w·n̂ = 0, valence, orthogonality where guaranteed, lake at rest under the
  legacy, WLSQ, diamond and cell-momentum options, exact mass balance, import
  round trip, raster-aligned grid.

#### Fixed
- `a5_bridge.mesh_to_geoparquet` padded polygons with fewer vertices than the
  longest with (0, 0) vertices (a spike to the origin) and overwrote an
  unclosed ring's last vertex. Only meshes with mixed vertex counts were
  affected (not A5 or Cartesian meshes).
- Domain-boundary detection (`_build_ghost_edges`): a cell is on the boundary
  if its polygon has more sides than neighbours. It was "fewer neighbours
  than the mesh maximum", which on a mixed-valence mesh gave a max-valence
  edge cell no ghost edge. Identical on A5 and Cartesian meshes (ghost
  edges unchanged on Carlisle res-18 and planar res-16).
- Degenerate-face guards in `A5Grid._edge_geometry`: 1 mm instead of 1 m.
  A 1 m face is normal on a Voronoi mesh (and on A5 at resolution ≥ 22); such
  faces got a zero normal. A5 edge geometry is bit-identical on the Carlisle
  and planar meshes.

### WP-M1: pya5 1.x support (`maint-pya5-1x`)

Work package WP-M1 from `FloodA5_Review_and_ForwardPlan_2026-10.md` (§9).
pya5 0.7.3 behaviour is unchanged; pya5 1.x is now supported.

#### Added
- One code path per pya5 major version in `mesh/a5_bridge.py` (`A5_MAJOR`,
  `_uncompact_to`, `_grid_disk_ring1`) and `mesh/A5Grid.jl`
  (`pya5_version()`, `a5_id_scheme()`, version-aware `grid_disk_neighbours`).
  Under 1.x the AOI seed set comes from `polygon_to_cells` (overlapping
  containment): Carlisle res-18 mesh generation 1 m 02 s vs 13 m 13 s.
- Cell-ID scheme recorded in every GeoParquet mesh (schema metadata key
  `flooda5`: `a5_id_scheme`, `pya5_version`) and in the HDF5 `/mesh`
  `id_scheme` attribute; new `A5Mesh.id_scheme` field (7th; the 6-argument
  constructor defaults it to the installed pya5's scheme; Cartesian meshes
  use `"cartesian"`).
- ⚠ `load_mesh_geoparquet` refuses a mesh whose A5 ID scheme differs from the
  installed pya5's (`check_id_scheme`; override with
  `FLOODA5_ALLOW_ID_SCHEME_MISMATCH=1`). Existing meshes are `a5-v0`, so under
  pya5 1.x they must be migrated first.
- `mesh/migrate_mesh_ids.py`: rewrites `cell_id`/`neighbours` with
  `a5.migrate`, verifying every centre. Verified on Carlisle res-18: a
  simulation on the migrated mesh is bit-identical to the original, and a
  mesh generated afresh under 1.x has the same cells and adjacency graph
  (centroids within 1e-13°, areas within 6e-10 relative).
- `A5Grid.migrate_cell_id`; `test/test_pya5_compat.jl` (in CI) with the fixture
  `test/fixtures/a5_id_schemes.csv`; CI matrix over pya5 0.7.3 and 1.0.0
  (pya5 removed from `ci/python-requirements.txt`); CI also runs on `maint-*`
  branches.

### WP0: correctness and test infrastructure (`wp0-correctness`)

Correctness and test-infrastructure work package WP0 from
`FloodA5_Review_and_ForwardPlan_2026-10.md` (§4.4, §9). **Some fixes change
numerical results** (marked ⚠); re-baseline any stored comparison outputs.

#### Fixed
- ⚠ **Data race in Q-centred momentum** (`step_standard!`, `step_sgs!`): the
  threaded edge loop read neighbouring edges' momentum from the same array it
  was overwriting. Q-centring now reads a start-of-step snapshot, giving a
  proper Jacobi update that is deterministic across thread counts. Affects
  every run with `--momentum-model edge` and `--q-centre-theta` < 1 (the
  defaults), standard and SGS.
- ⚠ **SGS CFL depth cap never applied** (`_cfl_dt`): `method == SGSFlow`
  compared an instance with a type (always false); now `method isa SGSFlow`.
  SGS runs with overfull channel cells will take larger, correctly capped
  timesteps.
- **`_adjacency_slot`** searched only 5 slots by default and returned slot 1
  silently when a neighbour was not found; it now searches the mesh's own
  neighbour ceiling, `step_sgs!` passes it explicitly, and
  `initialise_flow_model` verifies SGS slot consistency at start-up (and
  errors clearly if a mesh has more neighbours than the SGS tables support).
- **Diamond degenerate vertex** (`eta_at_vertex`) returned 0.0 m instead of
  NaN, defeating the documented per-edge fallback to the legacy gradient.
- **O(n²) allocation in SGS initialisation** (`A5Grid.sgs_table`): the
  `get(dict, key, fill(NaN, n))` default was evaluated eagerly for every cell.
- ⚠ (diagnostic output only) **Cell velocity** (`_compute_velocity!`) is now the
  depth-averaged velocity in m/s from the WLSQ cell-vector reconstruction.
  The previous scatter used degree-space directions with no cos(lat) factor
  and had units of s⁻¹. HDF5 `velocity`/`vel_u`/`vel_v` magnitudes change.
- Startup warning and docs wrongly said `--gradient-correction` has no effect
  under `--flow-model sgs`; it does (legacy WLSQ + V̂, unvalidated). New
  warning for `--face-flux-method diamond` under SGS.

#### Added
- `FLOODA5_NO_VIS=1` environment variable: stubs the visualisation modules so
  FloodModel.jl loads without GLMakie/Oxygen (headless, CI, AI-assisted runs).
- `test/test_wp0_regressions.jl`: regression tests for the fixes above, plus a
  lake-at-rest (well-balancing) test across all flux/momentum options.
- `test/check_determinism.jl`: thread-count determinism fingerprint.
- Continuous integration: `.github/workflows/ci.yml`, minimal `ci/Project.toml`
  environment, `ci/run_tests.sh`, `ci/python-requirements.txt`.

### Cartesian formulation (`cartesian-formulation`)

Full account: `HYDRAULICS.md` §14 and the project's internal development
history (available on request).

#### Added
- **Cartesian mesh backend** (`mesh/CartesianGrid.jl`): a rectangular,
  4-connected grid generator, pure Julia (no Python bridge for geometry
  generation). New CLI flags `--meshgen-cartesian FILE`, `--dx METRES`,
  `--cartesian-lat-ref DEG`. Covers the AOI's bounding box (not clipped to
  its polygon shape). Returns the existing `A5Mesh`/`A5Cell` types, so
  everything downstream of mesh acquisition is unchanged and
  backend-agnostic.
- Diagnostic scripts: `test/audit_cartesian_geometry.jl`,
  `test/test_octagonal_anisotropy.jl`, `test/check_boundary_proximity.jl`.

#### Changed
- **`N_SIDES` generalisation.** A new `FlowState.max_neighbours` field,
  derived from the loaded mesh's own adjacency, replaces hardcoded
  neighbour-count assumptions in matrix sizing (`adj_matrix`,
  `wlsq_weights`, `mom_weights`, `cell_edge_index`), the donor-edge safety
  divisor (`_donor_edge_divisor(state)`), the Fix B volume limiter (new
  `vol_limit_divisor` keyword on `_bates_flux_limited`/
  `_bates_flux_limited_corrected`, default `5.0` = unchanged A5 behaviour),
  and the ghost-edge boundary threshold. **No behavioural change on A5** —
  confirmed by the full existing regression suite. `N_SIDES`/
  `DONOR_EDGE_DIVISOR` remain as A5-only legacy constants for existing test
  fixtures. A backward-compatible 31-argument `FlowState` constructor keeps
  every existing call site working unmodified.
- **Cell-ID normalisation hardened.** Six previously-duplicated, unsafe
  inline `A5Grid._to_hex(parse(UInt64, id, base=16))` calls (which threw on
  any non-hex ID) now route through `_norm_cell_id` / `_to_hex_safe`, which
  fall back to the raw string — so a non-A5 backend's own ID scheme no
  longer crashes mesh loading.

#### Known issues / open
- A real, small, **resolution-independent** cardinal-vs-diagonal
  directional bias was found on the (provably exactly orthogonal)
  Cartesian backend; it got *worse* under grid refinement, ruling out
  ordinary truncation error. Leading hypothesis is the per-edge momentum
  representation; the confirming test (`--momentum-model cell`) has not
  been run. See `HYDRAULICS.md` §14.
- SGS (`--flow-model sgs`) is **not** generalised: `build_sgs_tables!` still
  hardcodes 5 neighbour slots, so a Cartesian mesh is unsupported with SGS.
  There is currently no guard preventing the combination.
- Projected-CRS support does not exist: the Cartesian generator works in
  decimal degrees via a single reference latitude, so exact raster
  alignment against a LISFLOOD-FP grid in its own projected CRS is not yet
  possible.

## [0.1.0] — 2026-08-26

First formally versioned release. Merges the `directional-bias-reformulation`
branch to `main` and brings the documentation set up to date with the
current state of the model.

### Added
- Diamond face-flux reconstruction (`--face-flux-method diamond`) — a
  face-local, polygon-native alternative to the cell-averaged WLSQ gradient
  correction, proven exact for linear WSE fields.
- Cell-vector momentum representation (`--momentum-model cell`) — a
  Perot-style reconstruction replacing five independent per-edge momentum
  scalars with one coherent 2D vector per cell.
- `--gradient-correction-alpha` — scaling factor for the tangential
  correction term in the legacy (WLSQ) gradient-correction method.
- A startup warning when `--momentum-model cell`, `--gradient-correction
  on`, `--face-flux-method`, or `--gradient-correction-alpha` are combined
  with `--flow-model sgs` (the default flow model), since none of them
  currently have any effect there.
- `--version`/`-v` CLI flag; version now printed in the startup banner and
  `--help` header.
- `VERSION`, `VERSIONING.md`, and this `CHANGELOG.md`.

### Fixed
- `gradient_correction`'s keyword default in `initialise_flow_model` and
  `run_flood_model` was `true`, contradicting the CLI's own default of
  `false`. CLI runs were unaffected (the CLI always passes an explicit
  value), but several test files calling `initialise_flow_model` directly
  without specifying this argument were silently using the corrected
  kernel rather than the intended uncorrected baseline. Both defaults
  corrected to `false`.
- `--help` was missing documentation for the (fully working)
  `--momentum-model` and `--gradient-correction-alpha` flags.
- The resolution-to-cell-area guidance shown in `--help` was incorrect by
  roughly 16× at resolution 14, and disagreed with the (correct, internally
  consistent) table in `docs/A5_QUIRKS.md`. Corrected to match.

### Documentation
- `docs/HYDRAULICS.md` rewritten in full — the previous version predated
  the stability-limiter fixes, the SGS R-A flux path, and all of the
  directional-bias correction work.
- `docs/METHODS.md` §5.4, `README.md`: corrected statements that the
  non-orthogonal gradient correction was enabled by default and had
  resolved the directional bias — neither is accurate; the correction is
  opt-in and the bias is reduced, not eliminated, even when enabled.
- `docs/USER_GUIDE.md`: added a complete CLI reference section for the
  five directional-bias-correction flags, previously undocumented.
- `docs/DATA_FORMATS.md`: corrected the documented mesh adjacency schema
  (a single variable-length `neighbours` column, not fixed `adj_0`..`adj_4`
  columns), added missing `vel_u`/`vel_v` HDF5 datasets, corrected
  `flux_Q`'s presence condition, fixed an off-by-one in the SGS elevation
  bins array length, and documented that `.bci` `HFIX`/`HVAR` entries are
  parsed but not yet implemented.
- `docs/ARCHITECTURE.md`: `FlowState`/`EdgeList` struct references updated
  to match the current code (both were missing roughly twenty fields added
  since the doc was last touched); added the `mesh/DiamondFlux.jl` module;
  corrected the simulation step-function phase descriptions; updated the
  test-file inventory; clarified the mesh-generation thread-safety guidance
  (see below).
- `docs/A5_QUIRKS.md`: corrected references to three functions
  (`_edge_cos_theta`, `_build_adjacency_grid_disk`) that no longer exist
  under those names, having been consolidated into `_edge_geometry` and
  `_build_adjacency_shared_vertices`/`_build_adjacency_matrix!`.
- `docs/TEST_CASES.md`: corrected the Square Domain test description to
  match the actual working test script (resolution and workflow had
  drifted from what's implemented) and removed references to two files
  that don't exist in the repository; added full documentation of the
  `test/planar_embankment/` test case, previously undocumented despite
  being the primary acceptance test for the directional-bias correction
  work.
- Mesh-generation thread-safety guidance (`--threads 1`) reworded across
  all docs to be evidence-based: two specific historical crash causes are
  confirmed fixed, but the same general hazard (a GDAL/PyCall coordinate
  transform called from inside a threaded loop) remains present in two
  other code paths (standard-flow DEM sampling and SGS hypsometric-curve
  building), so `--threads 1` remains the recommended default with an
  explicit fallback note rather than being removed.
- Added a "Development Roadmap" section and an "Acknowledgements" section
  to `README.md`; removed direct references to internal-only planning and
  handover documents from the public documentation set, replaced with
  general pointers to the project's internal development history
  (retained, not published, available on request).

[0.1.0]: https://github.com/ (tag pending)
