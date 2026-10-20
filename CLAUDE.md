# CLAUDE.md — working notes for AI-assisted sessions on FloodA5

This file is read automatically by Claude Code (CLI or VS Code extension) when
it is started in this repository. Keep it short and current.

## What this is
FloodA5: a 2-D local-inertial (Bates et al. 2010) flood model in Julia on the
A5 pentagonal DGGS (plus a Cartesian backend). The core is a flat script,
`FloodModel.jl` (not a package/module), which `include`s `mesh/`,
`surfacewater/`, `boundaryinputs/` and `visualisation/`.

Current plan of record: `projectinstructions/FloodA5_Review_and_ForwardPlan_2026-10.md`
(work packages WP0–WP7). Branch per work package (`wp0-correctness`, `wp1-…`).

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
- Full simulations: `julia --project=. --threads auto FloodModel.jl --help`.
  Mesh *generation* (`--meshgen`) should use `--threads 1` (PyCall/GDAL hazard).
- Large data lives in `test/carlisle`, `test/planar_embankment`, `test/square`
  (git-ignored: `*.h5`, `*.parquet`, `*.tif`). Don't commit outputs.

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
6. Use `x isa SomeType`, not `x == SomeType`, for type tests.

## Conventions
- Physics changes: keep the old behaviour behind a flag until validated; update
  `CHANGELOG.md` `[Unreleased]`, and mark result-changing fixes with ⚠.
- Every behaviour change needs a test in `test/` that CI runs (`ci/run_tests.sh`).
- Report numbers with the run command that produced them; compare at matched
  simulation times.
