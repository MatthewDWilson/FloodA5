# FloodA5 notebooks

Worked, reproducible accounts of the test cases, in this order:

| Notebook | Cases (`cases/run.jl`) | What it covers |
|---|---|---|
| `01_planar.ipynb` | `planar_slope`, `planar_embankment`, `planar_wall` | Planar-slope symmetry, the sub-grid (SGS) embankment notch, and a wall with no gap (does a barrier inside a cell hold water?): DEM generation (Python), A5 and Cartesian mesh generation, runs, metrics, figures |
| `02_channel.ipynb` | `channel` (and `channel_series`) | The straight-channel representation series: A5 vs Cartesian, SWE vs local-inertial |
| `03_carlisle.ipynb` | `carlisle_peak`, `carlisle` | Carlisle against the full shallow-water reference: the 6 h before the Eden peak (river inflow only), and the recorded first-6 h test window (with its 10 mm/h test rain, and without) |
| `04_sgs_strip.ipynb` | `sgs_strip` | The sub-grid (SGS) model on a synthetic strip with a sub-cell channel, against analytic normal flow, storage and arrival times and a 2.5 m reference: the storage-measurement question, why SGS held water back until 0.6.0 (edge flow areas frozen above each cell's highest DEM sample; legacy donor cap), what 0.7.0 changed (face tables, 2-D conveyance, positivity, time step, flow-path geometry), the result against the 0.6.0 record, and 0.8.0 (DEM at pixel centres, critical-flow limit) |
| `05_performance.ipynb` | — (`benchmarks/`) | Performance on this machine and the others with results: step throughput on CPU threads, GPU Float64 and GPU Float32, thread scaling, Carlisle wall times, and the accuracy of GPU and Float32 runs against the CPU Float64 run. Runs the benchmarks the machine has no results for and saves them to `benchmarks/results/` |

Each notebook describes every simulation and shows the exact commands (from
the case registry, so text and commands cannot drift apart). It runs the
cases through the registry, which reuses completed runs, then analyses them
with `analysis/` (`FA`) and draws them with `visualisation/FloodA5Plots.jl` (`FP`).

## Kernel (one-time setup)

The notebooks use the kernel **"Julia FloodA5 1.12"** (`julia-flooda5`): Julia
started with `--project=<repo>` and `FLOODA5_NO_VIS=1`. IJulia lives in your
global Julia environment, not in the FloodA5 project:

```julia
using Pkg; Pkg.activate(); Pkg.add("IJulia")          # global environment (@v1.12)
using IJulia
IJulia.installkernel("Julia FloodA5", "--project=D:/FloodA5";   # your checkout
                     specname = "julia-flooda5", env = Dict("FLOODA5_NO_VIS" => "1"))
```

In VS Code (Jupyter extension), open a notebook and pick "Julia FloodA5 1.12"
as the kernel.

To run notebooks headless (`nbtools.jl execute`), Jupyter itself is needed.
Install it in its own conda environment, separate from the one PyCall/pya5
uses:

```bash
julia --project=. -e 'using Conda; Conda.add(["jupyter", "nbconvert", "nbclient"], :jupyter)'
```

## Outputs are not committed

Clear outputs before committing; CI fails on a notebook with outputs
(`test/test_notebooks.jl`):

```bash
julia --project=. notebooks/nbtools.jl strip          # all notebooks, in place
julia --project=. notebooks/nbtools.jl check          # what CI checks
julia --project=. notebooks/nbtools.jl execute notebooks/01_planar.ipynb   # → test/notebook_runs/
```

`strip` writes the file exactly as Jupyter does (sorted keys, one-space
indent), so it changes nothing but the outputs.

## Running time

The first run of a notebook makes its case runs: planar about 35 min (three cases, 12 runs and 12 meshes),
channel about 15 min, SGS strip about 40 min, Carlisle about 45 min (`carlisle`) plus 1.5 h (`carlisle_peak`:
the event to the Eden peak, two SWE runs of ~85 min) on 16 cores. Later runs reuse
them (`test/case_runs/`, git-ignored) and take a minute or two.

`05_performance.ipynb` is different: its first run on a machine runs the
benchmark suite (about an hour and a half; use an idle machine and the full
project, `--project=<repo>`, for the GPU), and its results file,
`benchmarks/results/<machine>.json`, is committed, so later runs — on any
machine — show every machine's results (`benchmarks/README.md`).
