# benchmarks/BenchSuite.jl
# ========================
# Performance benchmark suite (module `BenchSuite`, `BS`; 0.6.0): runs the
# benchmarks on this machine and keeps the results, one small JSON file per
# machine in `benchmarks/results/` (tracked in git), so that machines can be
# compared. Front ends: `benchmarks/run.jl` (command line) and
# `notebooks/05_performance.ipynb` (runs what is missing, draws the figures).
#
# Sections (each a list of rows, see benchmarks/README.md):
#   step      step throughput on square Cartesian grids, all cells wet
#             (test/benchmark_backend.jl): CPU threads, GPU Float64, GPU
#             Float32; standard flow and SWE; with the depth difference of
#             each device run from the CPU run after the same steps
#   threads   the CPU step at 1, 2, 4, … threads (one process per count)
#   carlisle  the Carlisle 6 h test runs (cases/registry/carlisle.jl) on A5
#             res 18 and Cartesian 20 m, from the command line: wall time
#             and, against the CPU Float64 run, storage, extent and depth
#             differences at 6 h; and a 5 m Cartesian mesh of the same domain
#             (~580 000 cells, the first hour) where a GPU is worthwhile
#
# Every measurement is a separate Julia process (thread counts are fixed at
# start-up; a failed or out-of-memory run does not take the suite with it).
# A machine is identified by host name and operating system (`machine()`;
# `label` or ENV["FLOODA5_BENCH_LABEL"] overrides the name). Running a
# section again replaces that section of this machine's file.

module BenchSuite

using JSON3, Dates, Printf

if isdefined(Main, :FloodA5Analysis)
    const FA = Main.FloodA5Analysis
else
    include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
    const FA = FloodA5Analysis
end
const FO = FA.FO

const REPO = rstrip(normpath(joinpath(@__DIR__, "..")), ['/', '\\'])
const RESULTS_DIR = joinpath(@__DIR__, "results")
const SCHEMA = 1
const SECTIONS = (:step, :threads, :carlisle)

"Project the runs use: the active one, or ENV[\"FLOODA5_PROJECT\"] (the GPU needs the full project)."
project() = get(ENV, "FLOODA5_PROJECT", (ap = Base.active_project(); ap === nothing ? REPO : dirname(ap)))
"Where the Carlisle benchmark runs and the 5 m mesh go (git-ignored)."
work_dir() = joinpath(get(ENV, "FLOODA5_CASES_DIR", joinpath(REPO, "test", "case_runs")), "perf")
carlisle_dir() = joinpath(get(ENV, "FLOODA5_CASES_DIR", joinpath(REPO, "test", "case_runs")), "carlisle")

# ───────────────────────────── this machine ─────────────────────────────

_slug(s) = strip(replace(lowercase(String(s)), r"[^a-z0-9]+" => "-"), '-')

function _julia(args::Vector{String}; threads = 1)
    env = copy(ENV); env["FLOODA5_NO_VIS"] = "1"
    return Cmd(setenv(`$(Base.julia_cmd()) --project=$(project()) --startup-file=no --threads $threads $args`, env);
               dir = REPO)
end

function _physical_cores()
    try
        if Sys.iswindows()
            out = readchomp(`powershell -NoProfile -Command "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum"`)
            return parse(Int, strip(out))
        elseif Sys.islinux()
            ids = Set{Tuple{String,String}}(); phys = ""
            for l in eachline("/proc/cpuinfo")
                startswith(l, "physical id") && (phys = strip(split(l, ":")[2]))
                startswith(l, "core id") && push!(ids, (phys, strip(split(l, ":")[2])))
            end
            return isempty(ids) ? 0 : length(ids)
        end
    catch
    end
    return 0
end

const _GPU = Ref{Any}(missing)
"""
    gpu() → NamedTuple or nothing

The GPU that `--backend gpu` would use (name, memory, CUDA runtime and driver
versions), asked of CUDA.jl in a separate process with the runs' project;
`nothing` if CUDA.jl is not in that project or finds no functional GPU.
"""
function gpu()
    _GPU[] === missing || return _GPU[]
    code = """
        using CUDA
        CUDA.functional() || exit(3)
        d = CUDA.device()
        println("GPU|", CUDA.name(d), "|", CUDA.totalmem(d), "|", CUDA.runtime_version(), "|",
                something(tryparse(VersionNumber, string(CUDA.driver_version())), CUDA.driver_version()))
        """
    info = nothing
    try
        out = read(pipeline(_julia(["-e", code]); stderr = devnull), String)
        for l in eachline(IOBuffer(out))
            p = split(l, "|")
            length(p) == 5 && p[1] == "GPU" || continue
            info = (name = String(p[2]), memory_gb = round(parse(Float64, p[3]) / 2^30; digits = 1),
                    cuda_runtime = String(p[4]), driver = String(p[5]))
        end
    catch
    end
    return _GPU[] = info
end

"""
    machine(; label = "") → Dict

This machine: `id` (from the host name and operating system, or from
`label` / ENV["FLOODA5_BENCH_LABEL"]), CPU, cores, memory, GPU, Julia.
"""
function machine(; label::AbstractString = "")
    label = label != "" ? label : get(ENV, "FLOODA5_BENCH_LABEL", "")
    os = Sys.iswindows() ? "windows" : Sys.islinux() ? "linux" : Sys.isapple() ? "macos" : lowercase(string(Sys.KERNEL))
    host = gethostname()
    g = gpu()
    return Dict{String,Any}(
        "id" => _slug(label != "" ? label : host * "-" * os), "label" => label != "" ? label : host,
        "hostname" => host, "os" => os, "os_version" => string(Sys.KERNEL, " ", Sys.MACHINE),
        "cpu" => strip(Sys.cpu_info()[1].model), "threads" => Sys.CPU_THREADS, "cores" => _physical_cores(),
        "memory_gb" => round(Sys.total_memory() / 2^30; digits = 1),
        "gpu" => g === nothing ? "" : g.name, "gpu_memory_gb" => g === nothing ? 0.0 : g.memory_gb,
        "cuda_runtime" => g === nothing ? "" : g.cuda_runtime, "gpu_driver" => g === nothing ? "" : g.driver,
        "julia" => string(VERSION))
end

# ───────────────────────────── results files ─────────────────────────────

results_path(id) = joinpath(RESULTS_DIR, "$id.json")

"""
All machines' results (`benchmarks/results/*.json`), as Dicts, by machine id.
The smoke-test files of `run_suite(quick = true)` (`*-quick.json`) are left
out unless `quick = true`.
"""
function load_results(; quick::Bool = false)
    isdir(RESULTS_DIR) || return Dict{String,Any}[]
    out = Dict{String,Any}[]
    for f in sort(readdir(RESULTS_DIR))
        endswith(f, ".json") || continue
        (!quick && endswith(f, "-quick.json")) && continue
        r = JSON3.read(read(joinpath(RESULTS_DIR, f), String), Dict{String,Any})
        get(r, "schema", 0) == SCHEMA || (@warn "$f: results schema $(get(r, "schema", "?")), expected $SCHEMA — skipped"; continue)
        push!(out, r)
    end
    return out
end

"This machine's results (an empty record if it has none yet)."
function load_mine(; label = "")
    m = machine(; label)
    p = results_path(m["id"])
    isfile(p) && return JSON3.read(read(p, String), Dict{String,Any})
    return Dict{String,Any}("schema" => SCHEMA, "machine" => m, "sections" => Dict{String,Any}())
end

function _save(r)
    mkpath(RESULTS_DIR)
    open(results_path(r["machine"]["id"]), "w") do io
        JSON3.pretty(io, r); println(io)
    end
    return results_path(r["machine"]["id"])
end

function _provenance()
    git(args...) = try readchomp(Cmd(`git $args`; dir = REPO)) catch; "" end
    return Dict{String,Any}("date" => Dates.format(now(), "yyyy-mm-dd HH:MM"),
                            "flooda5_version" => strip(read(joinpath(REPO, "VERSION"), String)),
                            "git_commit" => git("rev-parse", "--short", "HEAD"),
                            "git_dirty" => !isempty(git("status", "--porcelain", "--untracked-files=no")),
                            "julia" => string(VERSION))
end

"Rows of a section of a results record (empty if it has not been run)."
rows(r, section) = (s = get(r["sections"], String(section), nothing); s === nothing ? Any[] : s["rows"])

# ───────────────────────────── running things ─────────────────────────────

_log(io, a...) = (println(io, a...); flush(io))

# run a Julia process, its output to `logfile`; → success
function _spawn(args; threads, logfile)
    mkpath(dirname(logfile))
    return open(logfile, "w") do log
        success(pipeline(_julia(args; threads); stdout = log, stderr = log))
    end
end

function _thread_counts(maxt)
    c = [1]
    while 2 * c[end] < maxt; push!(c, 2 * c[end]); end
    c[end] == maxt || push!(c, maxt)
    return c
end

# one process of test/benchmark_backend.jl → its JSON rows
function _backend_rows(io, solver, sizes, steps, configs, threads; min_cellsteps = 0)
    tmp = joinpath(work_dir(), "step_$(solver)_t$(threads).jsonl")
    mkpath(work_dir()); rm(tmp; force = true)
    logfile = replace(tmp, ".jsonl" => ".log")
    ok = _spawn(["test/benchmark_backend.jl", "--solver", solver, "--sizes", join(sizes, ","),
                 "--steps", string(steps), "--min-cellsteps", string(min_cellsteps),
                 "--configs", join(configs, ","), "--json", tmp];
                threads, logfile)
    out = Any[]
    isfile(tmp) && for l in eachline(tmp)
        isempty(strip(l)) || push!(out, JSON3.read(l, Dict{String,Any}))
    end
    ok || _log(io, "    ! the process failed after $(length(out)) row(s) — see $logfile")
    return out
end

function _step_table(io, rws)
    for r in rws
        _log(io, @sprintf("    %-8s %5d² %9d cells  %-7s %3d thr  %8.1f Mc·s/s", r["solver"], r["n"], r["cells"],
                          r["config"], r["threads"], r["mcs"]),
             r["max_dh"] === nothing ? "" : @sprintf("   max |Δh| %.1e m", r["max_dh"]))
    end
end

"""
    run_step(; solvers, sizes, steps, configs, threads, io) → section Dict

Step throughput (file header): one `test/benchmark_backend.jl` process per solver.
Each grid runs at least `steps` steps and at least `min_cellsteps` cell-steps
(20 steps of a small grid are too short to time on a GPU).
"""
function run_step(; solvers = ["standard", "swe"], sizes = [256, 512, 1024, 2048], steps = 20,
                  min_cellsteps = 100_000_000,
                  configs = gpu() === nothing ? ["cpu"] : ["cpu", "gpu64", "gpu32"],
                  threads = Sys.CPU_THREADS, io = stdout)
    out = Any[]
    for solver in solvers
        _log(io, "  step: $solver, grids $(join(sizes, ", ")), $(join(configs, ", ")), $threads threads")
        r = _backend_rows(io, solver, sizes, steps, configs, threads; min_cellsteps)
        _step_table(io, r); append!(out, r)
    end
    return merge(_provenance(), Dict{String,Any}("rows" => out,
                 "params" => Dict{String,Any}("sizes" => sizes, "steps" => steps, "min_cellsteps" => min_cellsteps, "configs" => configs, "threads" => threads)))
end

"""
    run_threads(; solvers, n, steps, counts, io) → section Dict

The CPU step on an n × n grid at each thread count (one process each).
"""
function run_threads(; solvers = ["standard", "swe"], n = 1024, steps = 20,
                     counts = _thread_counts(Sys.CPU_THREADS), io = stdout)
    out = Any[]
    for solver in solvers, t in counts
        _log(io, "  threads: $solver, $(n)², $t thread(s)")
        r = _backend_rows(io, solver, [n], steps, ["cpu"], t)
        _step_table(io, r); append!(out, r)
    end
    return merge(_provenance(), Dict{String,Any}("rows" => out,
                 "params" => Dict{String,Any}("n" => n, "steps" => steps, "counts" => counts)))
end

# ── Carlisle ────────────────────────────────────────────────────────────────

const _CARL = joinpath("test", "carlisle")
const _CONFIG_FLAGS = Dict("cpu" => String[], "gpu64" => ["--backend", "gpu"],
                           "gpu32" => ["--backend", "gpu", "--precision", "float32"],
                           "kacpu64" => ["--backend", "ka-cpu"],
                           "kacpu32" => ["--backend", "ka-cpu", "--precision", "float32"])

_carl_args(mesh, model, duration) =
    ["--meshload", mesh, "--flow-model", model, "--inflow-bci", joinpath(REPO, _CARL, "carlisle_wgs84.bci"),
     "--rainfall", "10", "--sim-duration", string(duration), "--dt-max", "10",
     "--output-interval", string(min(900, duration))]

"The command of a Carlisle benchmark run, as typed from the repository."
carlisle_command(mesh, model, config, duration, threads, out) =
    "julia --project=. --threads $threads FloodModel.jl " *
    join((_rel(a) for a in vcat(_carl_args(mesh, model, duration), _CONFIG_FLAGS[config], ["--output", out])), " ")
_rel(a) = startswith(a, REPO) ? join(splitpath(relpath(a, REPO)), "/") : a

# the meshes: (name, path, prepared?) — A5 res 18 and Cartesian 20 m are the `carlisle` case's
function _carl_meshes(io, large::Bool)
    lidar = joinpath(REPO, _CARL, "Carlisle_LiDAR_5m_mean.tif")
    out = Tuple{String,String}[]
    for (name, file) in (("a5_res18", "carlisle_a5_res18.parquet"), ("cart20", "carlisle_cart20.parquet"))
        p = joinpath(carlisle_dir(), file)
        isfile(p) ? push!(out, (name, p)) :
            _log(io, "  carlisle: $(relpath(p, REPO)) is missing — make it with " *
                     "`julia --project=. cases/run.jl prep carlisle` (needs test/carlisle/, see docs/TEST_CASES.md §3)")
    end
    if large
        p = joinpath(work_dir(), "carlisle_cart5.parquet")
        if !isfile(p) && isfile(lidar)
            _log(io, "  carlisle: building the 5 m Cartesian mesh ($(relpath(p, REPO)))")
            _spawn(["FloodModel.jl", "--meshgen-cartesian", lidar, "--meshout", p, "--mesh-only"];
                   threads = 1, logfile = joinpath(work_dir(), "prep_cart5.log")) ||
                _log(io, "    ! mesh generation failed — see $(joinpath(work_dir(), "prep_cart5.log"))")
        end
        isfile(p) ? push!(out, ("cart5", p)) :
            _log(io, "  carlisle: no 5 m mesh (needs $(relpath(lidar, REPO)))")
    end
    return out
end

# depth differences of `run` from `ref` at the last frame both reach
function _carl_compare(ref, run)
    a = FA.load_run(ref); b = FA.load_run(run)
    c = FA.compare(a, b)
    ha = FA.depth(a, c.k_ref); hb = FA.depth(b, c.k_run)
    ok = [i for i in eachindex(ha) if isfinite(ha[i]) && isfinite(hb[i])]
    d = abs.(hb[ok] .- ha[ok])
    return Dict{String,Any}("t_compare_s" => c.t, "storage_ratio" => c.storage_ratio, "csi" => c.csi,
                            "bias" => c.bias, "rms_dh_wet" => c.rms, "max_dh" => maximum(d),
                            "median_dh" => sort(d)[cld(length(d), 2)],
                            "frac_dh_gt_1cm" => count(>(0.01), d) / length(d))
end

"""
    run_carlisle(; models, configs, large, duration, large_duration, cpu_threads, io) → section Dict

The Carlisle runs (file header). `configs` are run in order; each run after
the first of its mesh and model is compared with that first one (the CPU
Float64 run, or GPU Float64 where the CPU is left out).
"""
function run_carlisle(; models = ["standard", "swe"],
                      configs = gpu() === nothing ? ["cpu"] : ["cpu", "gpu64", "gpu32"],
                      large = true, duration = 21600, large_duration = 3600,
                      cpu_threads = Sys.CPU_THREADS, gpu_threads = min(4, Sys.CPU_THREADS), io = stdout)
    out = Any[]
    for (mname, mesh) in _carl_meshes(io, large), model in models
        big = mname == "cart5"
        dur = big ? large_duration : duration
        # the 5 m mesh in SWE on CPU threads takes hours: devices only
        cfgs = big && model == "swe" ? filter(!=("cpu"), configs) : configs
        ref = nothing
        for config in cfgs
            threads = config == "cpu" ? cpu_threads : gpu_threads
            h5 = joinpath(work_dir(), "$(mname)_$(model)_$(config).h5")
            rm(h5; force = true)
            _log(io, "  carlisle: $mname, $model, $config\n    ", carlisle_command(mesh, model, config, dur, threads, h5))
            ok = _spawn(vcat(["FloodModel.jl"], _carl_args(mesh, model, dur), _CONFIG_FLAGS[config], ["--output", h5]);
                        threads, logfile = replace(h5, ".h5" => ".log"))
            s = ok && isfile(h5) ? get(FO.run_info(h5), "summary", Dict()) : Dict()
            if !get(s, "completed", false)
                _log(io, "    ! the run did not finish — see $(replace(h5, ".h5" => ".log"))")
                continue
            end
            cells = length(FO.read_frame(h5, "volume"))
            wall = s["wall_simulation_s"]; first5 = get(s, "wall_first5_steps_s", 0.0)
            row = Dict{String,Any}(
                "mesh" => mname, "model" => model, "config" => config, "threads" => threads, "cells" => cells,
                "sim_s" => dur, "wall_s" => wall, "first5_s" => first5, "steps" => s["n_steps"],
                "speed" => dur / wall,                                   # simulated seconds per wall second
                "storage_m3" => s["vol_domain_m3"],
                "conservation_rel" => abs(s["conservation_error_m3"]) / max(s["vol_applied_m3"], eps()),
                "command" => carlisle_command(mesh, model, config, dur, threads, h5))
            if ref === nothing
                ref = (config, h5)
            else
                row["ref_config"] = ref[1]
                merge!(row, _carl_compare(ref[2], h5))
            end
            _log(io, @sprintf("    %d cells, %d steps, wall %.1f s (first 5 steps %.1f s), %.0f× real time",
                              cells, row["steps"], wall, first5, row["speed"]),
                 haskey(row, "max_dh") ? @sprintf("; vs %s: storage ×%.5f, CSI %.4f, max |Δh| %.2e m",
                                                  ref[1], row["storage_ratio"], row["csi"], row["max_dh"]) : "")
            push!(out, row)
        end
    end
    return merge(_provenance(), Dict{String,Any}("rows" => out,
                 "params" => Dict{String,Any}("duration" => duration, "large_duration" => large_duration,
                                              "configs" => configs, "models" => models)))
end

# ───────────────────────────── the suite ─────────────────────────────

"""
    run_suite(; sections = [:step, :threads, :carlisle], force = false, quick = false,
              label = "", large = true, io = stdout) → this machine's results (Dict)

Run the sections this machine has no results for (`force = true`: all of
`sections` again) and save them to `benchmarks/results/<machine id>.json`.
`quick = true` is a smoke test of the whole chain (small grids, 15 simulated
minutes of Carlisle, no 5 m mesh) — its numbers are not benchmarks, and it
writes to `<machine id>-quick.json`.
"""
function run_suite(; sections = collect(SECTIONS), force::Bool = false, quick::Bool = false,
                   label::AbstractString = "", large::Bool = true, io = stdout)
    all(in(SECTIONS), Symbol.(sections)) || throw(ArgumentError("sections: $(join(SECTIONS, ", "))"))
    quick && (label = (label != "" ? label : gethostname() * "-" * machine()["os"]) * "-quick")
    r = load_mine(; label)
    r["machine"] = machine(; label)
    m = r["machine"]
    _log(io, "machine $(m["id"]): $(m["cpu"]), $(m["threads"]) threads; ",
         m["gpu"] == "" ? "no GPU for --backend gpu (CUDA.jl in $(project())?)" : "$(m["gpu"]) ($(m["gpu_memory_gb"]) GB)")
    for sec in Symbol.(sections)
        if !force && !quick && haskey(r["sections"], String(sec))
            _log(io, "$sec: results of $(r["sections"][String(sec)]["date"]) kept (force = true runs it again)")
            continue
        end
        _log(io, "$sec: running")
        t0 = time()
        r["sections"][String(sec)] =
            sec === :step ? (quick ? run_step(; sizes = [128, 256], steps = 5, min_cellsteps = 0, io) : run_step(; io)) :
            sec === :threads ? (quick ? run_threads(; n = 256, steps = 5, counts = unique([1, Sys.CPU_THREADS]), io) :
                                        run_threads(; io)) :
            (quick ? run_carlisle(; models = ["standard"], large = false, duration = 900, io) : run_carlisle(; large, io))
        r["sections"][String(sec)]["elapsed_s"] = round(time() - t0; digits = 1)
        _log(io, "$sec: done in $(round(Int, time() - t0)) s → ", relpath(_save(r), REPO))
    end
    return r
end

end # module
