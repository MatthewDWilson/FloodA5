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
(`VERSIONING.md`). Status (2026-11-10): 0.3.2 released (WP0–WP5, the SWE
reference `--flow-model swe`, the Δt-consistent 0.3.0 default, CPU/GPU
backends, the channel series). Current phase: **consolidation** — reproducible
self-describing HDF5 output (done: 0.4.0, below), one analysis module, a test-case
registry and runner, HDF5 replay in the GLMakie viewer, then Jupyter notebooks
(planar → channel → Carlisle); see `projectinstructions/PROJECT_STATE.md`
"Pending Work → Current plan". Physics results: `docs/HYDRAULICS.md` §15
(§15.11 the 0.3.0 defaults, §15.12 channel series). SGS meshes built before
2026-10-21 are refused until rebuilt or
repaired (`mesh/repair_sgs_slot_order.py`). Edge orientation (`cell_i < cell_j`)
is not geometric on A5: never combine stored edge fluxes of different edges
without `_orient_collinear!`-style signs.

## Running things
- Always run from the repo root with `--project=.` (full stack) or
  `--project=ci` (minimal headless stack, fast to load).
- Set `FLOODA5_NO_VIS=1` for anything that doesn't need a visualiser — it
  stubs GLMakie/Oxygen and cuts load time substantially.
- Unit/regression suite (same as CI):
  `bash ci/run_tests.sh .`  (Git Bash) — or run single files:
  `julia --project=. test/test_wp0_regressions.jl`
- Thread determinism: `julia --project=. -t 1 test/check_determinism.jl` vs `-t 4`;
  compare the `FINGERPRINT` lines.
- Step benchmark: `julia --project=. -t 8 test/benchmark_baseline.jl`.
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
- Python: FloodA5 uses PyCall's Python (Julia's Conda.jl Python,
  `~/.julia/conda/3`), not other conda envs; it needs conda's `Library\bin` on
  PATH, so call it from Julia rather than straight from Git Bash. pya5 is 1.x
  locally; meshes carry an ID scheme (`a5-v0`/`a5-v1`) and old ones are refused —
  migrate with `julia setup.jl --migrate-meshes`. To test against another pya5
  without touching that env, put a venv's `site-packages` first on `PYTHONPATH`.
- Long runs / parallel work: use `git worktree` so code edits can't change jobs
  that have not started yet. On Windows, stopping a background shell can leave
  its `julia.exe`/`xargs` children running — check and stop them explicitly.

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

## Conventions
- Physics changes: keep the old behaviour behind a flag until validated; update
  `CHANGELOG.md` `[Unreleased]`, and mark result-changing fixes with ⚠.
- Every behaviour change needs a test in `test/` that CI runs (`ci/run_tests.sh`).
- Report numbers with the run command that produced them; compare at matched
  simulation times.
