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
