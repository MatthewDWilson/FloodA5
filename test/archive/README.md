# test/archive — superseded run scripts

These scripts drove the experiments of the directional-bias work (Phases A–C,
WP1–WP2), the FOSS4G 2026 runs and the early Carlisle tests. They are kept
for the record: the results they produced are in `docs/HYDRAULICS.md` (§7–8,
§15) and `projectinstructions/`.

**They do not reproduce their recorded numbers on the current code.** None
of them passes `--scheme`, and since 0.3.0 the default (`--scheme
recommended`) changes the limiter, Q-centring, friction and free-surface
filter they ran with. To repeat one of those experiments, add `--scheme
legacy` (0.1.0 behaviour) or `--scheme v0.2` to its commands. Several also
have their own bugs, listed in the 2026-11 inventory (wrong `--meshout` or
`--meshres` in the FOSS4G Carlisle blocks; meshgen without `--flow-model
standard` and with `--threads auto` in `carlisle_tests.ps1`; the last
comparison in `alpha_sweep.ps1` reads the wrong file).

The paths inside the scripts are relative to the repository root, as they were
before the move; run them from there (`pwsh test/archive/<script>.ps1`).

The current, runnable cases are in the test-case registry (`cases/`,
`docs/TEST_CASES.md`):

| Archived script | Registry case |
|---|---|
| `run_planar_slope_tests.ps1`, `benchmark_planar_slope_tests.ps1`, `run_diamond_dt_sweep.ps1`, `alpha_sweep.ps1`, `run_task4_convergence_sweep.ps1`, `run_flow_direction_tests.ps1`, `run_phaseC_level2_tests.ps1` | `planar_slope` (the current-scheme benchmark; the parameter sweeps themselves are not re-registered) |
| `run_wp1_matrix.sh` | `point_spread`, `planar_slope` |
| `run_wp2_carlisle.sh`, `carlisle/carlisle_tests.ps1`, `run_FOSS4G2026_simulations.ps1` | `carlisle` |
| `diag_d2_summary.sh` | — (summarises the damping runs; `test/diag_divergence.jl` and `analysis/compare.jl` do the analysis) |

Still in `test/`: `run_jobs.sh` (the generic matrix runner the channel
series used) and `run_benchmarks.{ps1,sh}` (thread-scaling performance).
