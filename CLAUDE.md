# CLAUDE.md — working notes for AI-assisted sessions on FloodA5

This file is read automatically by Claude Code (CLI or VS Code extension) when
it is started in this repository. Keep it short and current.

## What this is
FloodA5: a 2-D local-inertial (Bates et al. 2010) flood model in Julia on the
A5 pentagonal DGGS (plus a Cartesian backend). The core is a flat script,
`FloodModel.jl` (not a package/module), which `include`s `mesh/`,
`surfacewater/`, `boundaryinputs/` and `visualisation/`.

Physics plan of record: `projectinstructions/FloodA5_Review_and_ForwardPlan_2026-10.md`
(work packages WP0–WP7). Branch per work package (CI runs on `main`, `wp*` and
`maint-*` pushes); every merged PR bumps the version and is tagged
(`VERSIONING.md`). Status (2026-11-26): 0.8.14 — the default engine is
`--flow-model standard` (SGS opt-in; `swe` a full option). Before: WP0–WP5, the SWE engine, the
Δt-consistent 0.3.0 scheme, CPU/GPU backends, the channel series, output
format 2 (0.4.0). The **consolidation** phase
is complete: reproducible self-describing HDF5 output (0.4.0, below), one
analysis module (0.4.1, `analysis/`), a test-case registry and runner
(0.4.2, `cases/`), HDF5 replay in the GLMakie viewer and static plots
(0.4.3, `visualisation/replay.jl`, `FloodA5Plots`), Jupyter notebooks
(0.4.4, `notebooks/`; commit them stripped: `notebooks/nbtools.jl strip`).
0.6.0: `--precision float32` on the device backends (standard, swe) and the
benchmark suite (`benchmarks/`, `notebooks/05_performance.ipynb`).
0.7.0: SGS conveyance rework (one table per face, connectivity, positivity,
Courant 0.5; `docs/HYDRAULICS.md` §9, §15.13); 5 m Carlisle references
(`carlisle_peak_fine`, GPU). 0.7.1: DGGS mesh families (Phase A).
0.8.0: DEM sampled at pixel centres (was half a pixel off; every case
re-recorded, old values in `cases/expected/v0.7.1/`), SGS critical-flow
limit (`--sgs-critical`), the 5 m SWE run is the Carlisle benchmark
(HYDRAULICS §15.14). SGS = the elevation representation; hydraulic
formulations on it are to be pluggable. 0.8.1: DGGS mesh families
(Phase B): `--meshgen-dggs`, mesh metadata kept into the output, the
benchmarks on eleven mesh families (`docs/TEST_CASES.md` §8). 0.8.2:
notebook 06 (mesh families; §8 = the SGS checkerboard per mesh — re-run it
with `RERUN_SGS = true` after SGS changes: SGS checkerboards on Cartesian
too, smooth only on A5, ISEA3H, jittered Voronoi). 0.8.3: the web globe
viewer rebuilt (`--vis cesium`, `visualisation/cesium_replay.jl`; USER_GUIDE
§6.2, DATA_FORMATS §6; Oxygen dropped, server on HTTP.jl; a bare `--vis` is
Makie). 0.8.4: the same viewer on Giro3D (`--vis giro3d`); both share
`visualisation/viewer-common/` (app, UI, protocol, surface GLSL; engines
are thin adapters); replay is `visualisation/web_replay.jl`. Choosing
between them: docs/GLOBE_VIEWERS.md. 0.8.5: access keys (Cesium ion token,
…) through `io/FloodA5Keys.jl` — environment first, then git-ignored
`keys.json` / `~/.flooda5/keys.json` (template `keys.example.json`); never
print a key's value (use `key_source`). 0.8.6: viewer fixes (offset per
basemap, run colour ranges `/api/ranges`, palettes, cell borders, depth
bias for coarse 3D tiles); the viewer plan (decisions 2026-11-25, Chen et
al. 2019 equations in appendix B):
`projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md`. 0.8.7 (plan
step 1): Run/About tabs (`/api/run`, `/api/about`), bed colours, brown
flood water palette; variable descriptions in one place, `io/variables.toml`
(`FO.variable_docs`), the display record `FO.run_record` — both also give
the static plots' captions (`caption = false` / `FP.CAPTIONS[]`).
0.8.8 (step 2): time series at cells (right-click → Add to plot; Chart.js,
`viewer-common/js/plot.js`; `/api/series`).
0.8.9 (plan step 3): hazard variables (D·V, Defra HR, ARR H1–H6, Chen et al.
hazard degree) and statistics over the run (maxima, first flooding,
duration) in `analysis/hazard.jl` — pure, included by both FA and the
viewers' server (exports only in FloodA5Analysis.jl); `/api/stats`,
hazard parameters in the query; classed palettes (`viewer-common/js/classes.js`).
0.8.10 (plan step 4): local files read in the browser (h5wasm in a worker,
`viewer-common/js/local/`: a JavaScript port of the replay server's /api and
of `analysis/hazard.jl` — change both together: `test/test_vis_local.jl`
checks them against each other, fixture from `test/vis_local_fixture.jl
--write`); the pages also run without a server (`defaults.json`, `about.json`
— regenerate it after editing `io/variables.toml`: `viewer-common/tools/make_about.jl`).
0.8.12: viewer fixes (parameters per layer, time of maximum depth). 0.8.13:
`/statistics` in the output, gathered by the model (`--stats-interval`,
`--flood-depth`) with `analysis/hazard.jl`'s `StatsAccumulator` — the same
code as the viewers' statistics; keep `run_statistics` a loop over it.
0.8.14: smooth water/bed surfaces (`geometry.js` `cornerIndex` /
`cornerValues`, corner texture; CesiumJS rebuilds its model with a corner
attribute once).
0.8.11 (plan step 5, palettes): cmocean and Crameri maps as fixed tables
(`viewer-common/colormaps.json`, from `tools/make_palettes.jl`), user
palettes (OKLab stops, editor, JSON export/import; `js/palettes.js`,
`js/palette-editor.js`) and a group's shared `palettes.json` (`/api/palettes`).
Next and open findings: `projectinstructions/PROJECT_STATE.md`
"Pending Work → Current plan". Physics results: `docs/HYDRAULICS.md` §15
(§15.11 the 0.3.0 defaults, §15.12 channel series, §15.13–15.14 SGS). SGS meshes
built before 0.8.0 (table layout < 4) are refused until rebuilt. Edge orientation (`cell_i < cell_j`)
is not geometric on A5: never combine stored edge fluxes of different edges
without `_orient_collinear!`-style signs.

## Running things
- Always run from the repo root with `--project=.` (full stack) or
  `--project=ci` (minimal headless stack, fast to load).
- Set `FLOODA5_NO_VIS=1` for anything that doesn't need a visualiser — it
  stubs GLMakie and cuts load time substantially. The globe viewers' server
  (`--vis cesium|giro3d`, plain HTTP.jl) still works; viewer JS tests: `cd
  visualisation/viewer-common && node --test test/*.test.js`, browser check
  (both engines): `visualisation/viewer-common/tools/smoke.mjs`.
- Unit/regression suite (same as CI):
  `bash ci/run_tests.sh .`  (Git Bash) — or run single files:
  `julia --project=. test/test_wp0_regressions.jl`
- Thread determinism: `julia --project=. -t 1 test/check_determinism.jl` vs `-t 4`;
  compare the `FINGERPRINT` lines.
- Step benchmark: `julia --project=. -t 8 test/benchmark_baseline.jl`; CPU/GPU,
  Float64/Float32 and whole runs, per machine: `julia --project=. benchmarks/run.jl`
  (results in `benchmarks/results/<machine>.json`, committed).
- Kernel bodies (`surfacewater/kernels.jl`, `*_device.jl`, `swe_reference.jl`)
  are generic in the float type: write literals inside `@fp T …`, never form
  z + h except through the level helpers (`_edge_levels`, `_swe_face_z`, …),
  and keep volume sums in `eltype(volume)` (`surfacewater/precision.jl`).
  Float64 results must stay bit-identical; Float32: `test/test_precision.jl`.
  Performance changes must keep results bit-identical unless stated: compare
  determinism fingerprints and Carlisle outputs frame by frame against `main`
  (use `--cell-order file`; spatial order changes round-off, and the standard
  solver amplifies it — HYDRAULICS.md §15.5). `surfacewater/reference_steps.jl`
  is frozen: the live steps must match it (`test/test_wp4_stage2.jl`).
  Phase A skips edges via the active set (`_active_edges!`): any new state an
  edge kernel reads must be covered by the skip condition there (an edge may
  be skipped only if visiting it would change nothing).
- Full simulations: `julia --project=. --threads auto FloodModel.jl --help`.
  Mesh *generation* (`--meshgen`) should use `--threads 1` (PyCall/GDAL hazard).
- Large data lives in `test/carlisle`, `test/planar_embankment`, `test/square`
  (git-ignored: `*.h5`, `*.parquet`, `*.tif`). Don't commit outputs; experiment
  runs go in `test/wp1_runs/` (git-ignored). Run matrices: `test/run_jobs.sh`.
- Test cases: `julia --project=. cases/run.jl list|show|run|check CASE`
  (`docs/TEST_CASES.md`). Cases generate their own inputs in
  `test/case_runs/<case>/`; never share a derived input (DEM, mesh) between
  cases. After a code change that can move results, run the affected cases
  with `--fresh` and `check`; re-record expected values (`baseline`) only for
  a deliberate change, and say so in the CHANGELOG. Old run scripts are in
  `test/archive/` and rely on pre-0.3.0 defaults.
- Python: FloodA5 uses PyCall's Python (Julia's Conda.jl Python,
  `~/.julia/conda/3`), not other conda envs; it needs conda's `Library\bin` on
  PATH, so call it from Julia rather than straight from Git Bash. pya5 is 1.x
  locally; meshes carry an ID scheme (`a5-v0`/`a5-v1`) and old ones are refused —
  migrate with `julia setup.jl --migrate-meshes`. To test against another pya5
  without touching that env, put a venv's `site-packages` first on `PYTHONPATH`.
- Long runs / parallel work: use `git worktree` so code edits can't change jobs
  that have not started yet. Case runs made in one worktree are reused from
  another (the runner maps repository paths); the runner does not see code
  changes — rerun with `--fresh` (or a notebook's rerun flag) after them.
  On Windows, stopping a background shell can leave its `julia.exe`/`xargs`
  children running — check and stop them explicitly.

## Julia gotchas that have bitten this project
1. `@printf` needs a literal format string — no `"a" * "b"`; put dynamic text
   in a separate `println`.
2. Never write `\"` inside a `$(...)` interpolation — it is a ParseError for
   the whole file.
3. `FlowState` / `EdgeList` positional constructors are called from many test
   fixtures. When adding a field: append it last and add a backward-compatible
   outer constructor (see `FlowState.max_neighbours`). Check every call site's
   argument count.
4. Don't read `N_SIDES` in production code — use `state.max_neighbours`.
5. Threaded loops (`Threads.@threads`) must only write per-iteration slots and
   must not read arrays that the same loop writes (see the WP0 Q-centring race).
   A closure must not capture a local that is assigned only in one branch
   (it gets boxed and every access becomes type-unstable — bind it first).
   Step loops use Polyester `@batch per=core` (never plain `@batch`: it spins
   on all logical cores); functions called inside must take `AbstractVector`.
6. Use `x isa SomeType`, not `x == SomeType`, for type tests.
7. Non-A5 polygon meshes have mixed valence: don't infer "boundary cell" from
   neighbour count vs the maximum (use unshared polygon sides), and don't use
   fixed 1 m geometric guards.

## Output files (format 2, since 0.4.0)
- The HDF5 format lives in `io/FloodA5Output.jl` (`FO` in `FloodModel.jl`):
  read output only through `FO.read_frame` / `read_var` / `read_series` /
  `frame_times` / `run_info` (they also read ≤ 0.3.x files), never by
  dataset path. Layout: `docs/DATA_FORMATS.md` §2.
- Each output records its run (argv, resolved options, git, inputs, the mesh)
  and `--rerun out.h5` must reproduce it bit for bit
  (`test/test_output_format.jl`). A new CLI flag that reads a file must be
  added to `_run_inputs` (FloodModel.jl); a new mesh-generation flag to
  `FO._DROP_FLAGS`; a new boolean flag to `FO._BOOL_FLAGS`.
- Experiment provenance: the command is in the file (`FO.run_info(p)["argv"]`),
  so a run's `.h5` is enough to report it.
- Analysis goes in `analysis/` (`FloodA5Analysis`, `FA`; tests in
  `test/analysis/`, run by `test/test_analysis.jl`): pure Julia, never
  include `FloodModel.jl` there. New metrics belong in the module, as pure
  functions of arrays plus a `Run`-level wrapper; scripts in `test/` should
  stay thin CLI wrappers over it.

## Conventions
- Physics changes: keep the old behaviour behind a flag until validated; update
  `CHANGELOG.md` `[Unreleased]`, and mark result-changing fixes with ⚠.
- Every behaviour change needs a test in `test/` that CI runs (`ci/run_tests.sh`).
- Report numbers with the run command that produced them; compare at matched
  simulation times.
