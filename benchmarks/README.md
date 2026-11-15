# FloodA5 performance benchmarks

A benchmark suite that runs on any machine and keeps its results, one small
file per machine, so that hardware can be compared. The figures and tables
are in `notebooks/05_performance.ipynb`, which also runs whatever the
machine it is opened on has no results for.

```bash
julia --project=. benchmarks/run.jl                 # run what this machine has no results for
julia --project=. benchmarks/run.jl --force         # run everything again
julia --project=. benchmarks/run.jl --sections step,threads
julia --project=. benchmarks/run.jl --no-large      # without the 5 m Carlisle mesh
julia --project=. benchmarks/run.jl --label NAME    # name the machine (default: host name)
julia --project=. benchmarks/run.jl --quick         # 15-minute smoke test; not a benchmark
julia --project=. benchmarks/run.jl --list          # machines with results
```

Use the full project (`--project=.`): the GPU runs need CUDA.jl. Use a
machine that is doing nothing else. The full suite takes about two hours on a
16-thread workstation with a GPU. Without a usable NVIDIA GPU the
CPU parts run and the GPU columns stay empty.

## Sections

| Section | What runs | Measures |
|---|---|---|
| `step` | `test/benchmark_backend.jl`: square Cartesian grids 256² – 2048², all cells wet, at least 20 steps and 10⁸ cell-steps after 2 warm-up steps; standard flow (default scheme) and SWE (order 2); CPU threads, GPU Float64, GPU Float32 | million cell-steps per second; depth difference of each device run from the CPU run after the same steps |
| `threads` | the CPU step on 1024² at 1, 2, 4, … threads, one process each | thread scaling |
| `carlisle` | `FloodModel.jl` runs with the `carlisle` case's forcing (three inflows + 10 mm/h rain): A5 res 18 and Cartesian 20 m for 6 h, Cartesian 5 m (~580 000 cells) for 1 h; standard and SWE; CPU, GPU Float64, GPU Float32 (5 m SWE: GPU only) | wall time of the simulation; at the end, against the first configuration of the mesh and model: storage ratio, extent CSI, depth differences; conservation error |

The 20 m and A5 Carlisle meshes are the `carlisle` case's
(`julia --project=. cases/run.jl prep carlisle`; needs `test/carlisle/`).
The 5 m mesh is built by the suite from
`test/carlisle/Carlisle_LiDAR_5m_mean.tif`. Missing meshes are skipped with
a message. Run outputs and logs go to `test/case_runs/perf/` (git-ignored).

## Results files

`benchmarks/results/<machine id>.json` — commit them. The machine id is the
host name and the operating system (`myhost-windows`, `myhost-linux`), or the
`--label`. Running a section again replaces that section in the file.

```
schema            1
machine           id, label, hostname, os, os_version, cpu, cores, threads, memory_gb,
                  gpu, gpu_memory_gb, cuda_runtime, gpu_driver, julia
sections.<name>   date, flooda5_version, git_commit, git_dirty, julia, elapsed_s, params, rows
```

Rows of `step` and `threads`: `solver` (standard | swe), `n` (grid side),
`cells`, `config` (cpu | gpu64 | gpu32), `threads`, `steps`, `mcs` (million
cell-steps per second), `setup_s` (set-up and the two warm-up steps: on a
GPU, kernel compilation), `max_dh`, `rms_dh` (m, against the cpu row of the
same grid; null for cpu), `mass_rel` (relative change of volume stored plus
volume removed).

Rows of `carlisle`: `mesh` (a5_res18 | cart20 | cart5), `model`, `config`,
`threads`, `cells`, `sim_s`, `wall_s`, `first5_s` (set-up and the first five
steps), `steps`, `speed` (simulated seconds per wall second), `storage_m3`,
`conservation_rel`, `command`; and, for every run but the first of its mesh
and model, `ref_config`, `t_compare_s`, `storage_ratio`, `csi` and `bias`
(extent deeper than 0.1 m), `rms_dh_wet`, `median_dh`, `max_dh` (m),
`frac_dh_gt_1cm`.

## Older scripts

`test/benchmark_baseline.jl` (the CPU step, with allocations) and
`test/profile_backend.jl` (per-phase GPU times) are still current.
`test/benchmark_sim.jl` with `test/run_benchmarks.ps1` / `.sh` (a thread
sweep from before 0.3) needs a `benchmark_config.json` that is no longer in
the repository; the `threads` section here replaces it.
