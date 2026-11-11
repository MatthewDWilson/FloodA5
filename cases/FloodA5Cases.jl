# cases/FloodA5Cases.jl
# =====================
# Test-case registry and runner (module `FloodA5Cases`, `FC`).
#
# A case is data: the inputs it prepares (DEM generators, mesh generation),
# the FloodModel.jl runs it makes (one argument list per job), the metrics it
# computes from them with the analysis module (`FA`), the physical checks
# those metrics must pass, and the expected values recorded for regression
# (`cases/expected/<case>.toml`). The command-line front end is
# `cases/run.jl`; the cases themselves are in `cases/registry/*.jl`.
#
# Paths in arguments may use the placeholders `{case}` (the case's directory,
# `<root>/<case>/`), `{runs}` (`{case}/runs`) and `{repo}` (the repository).
# The root is `test/case_runs/` (git-ignored), or `$FLOODA5_CASES_DIR`.
#
# Runs are separate `julia FloodModel.jl …` processes, `--output
# {runs}/<job>.h5` appended, log in `{runs}/<job>.log`. A run is reused if
# its output records the same command line and completed (format 2's run
# record); otherwise it is run again (`--fresh` forces it, e.g. after a code
# change). Preparation steps of one checkout must not run concurrently with
# mesh generation elsewhere in the same checkout before 0.4.2 (fixed temp
# files in mesh/A5Grid.jl).
#
# Registry files (`cases/registry/*.jl`, included in name order) share this
# module's namespace: prefix every helper with the case's abbreviation
# (`_pt_` point spread, `_ps_` planar slope, `_pe_` planar embankment,
# `_ch_` channel, `_carl_` Carlisle); `_planar_common.jl` holds shared ones.

module FloodA5Cases

using Printf, Dates, TOML

if isdefined(Main, :FloodA5Analysis)
    const FA = Main.FloodA5Analysis
else
    include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
    const FA = FloodA5Analysis
end
const FO = FA.FO

export Case, Prep, Job, Check, register!, cases, case, case_dir, run_case, prepare,
       run_jobs, analyse, check, write_baseline, show_case, job_output

const REPO = rstrip(normpath(joinpath(@__DIR__, "..")), ['/', '\\'])
case_root() = get(ENV, "FLOODA5_CASES_DIR", joinpath(REPO, "test", "case_runs"))
"Project the runs use: the active one (CI: `ci/`), or `ENV[\"FLOODA5_PROJECT\"]`."
const PROJECT = Ref{String}("")
function __init__()
    ap = Base.active_project()
    PROJECT[] = get(ENV, "FLOODA5_PROJECT", ap === nothing ? REPO : dirname(ap))
end

# ─────────────────────────────── definitions ───────────────────────────────

"""
    Prep(name, kind, args; outputs, threads = 1)

A preparation step. `kind`:
- `:floodmodel`: `julia FloodModel.jl args…` (e.g. mesh generation; one thread,
  the PyCall/GDAL rule for `--meshgen`)
- `:julia`: `julia args…` (a script and its arguments)
- `:python`: FloodA5's Python (PyCall's, with conda's `Library/bin` on PATH)
  running `args…`
- `:write`: `args = [path, text]` — write a small input file (a hydrograph
  CSV, an AOI) given inline in the case definition
It is skipped when all `outputs` (paths, placeholders allowed) exist.
"""
struct Prep
    name::String
    kind::Symbol
    args::Vector{String}
    outputs::Vector{String}
    threads::Int
end
Prep(name, kind, args; outputs, threads = 1) =
    Prep(String(name), kind, String.(args), String.(outputs), threads)

"""
    Job(name, args; threads = 4, tags...)

One FloodModel.jl run. `args` without `--output` (the runner adds it). `tags`
(e.g. `mesh = "a5", res = 16`) label the job for analysis and selection.
"""
struct Job
    name::String
    args::Vector{String}
    threads::Int
    tags::Dict{Symbol,Any}
end
Job(name, args; threads = 4, tags...) =
    Job(String(name), String.(args), threads, Dict{Symbol,Any}(tags))

"""
    Check(name, test; detail = "")

A physical check: `test(metrics) -> Bool`, where `metrics[job][metric]` are
the case's metrics. `detail(metrics) -> String` explains the outcome.
"""
struct Check
    name::String
    test::Function
    detail::Function
end
Check(name, test; detail = m -> "") = Check(String(name), test, detail)

"""
    Case(; name, title, doc, tier, description, prep, jobs, analyse, checks,
         rtol, atol)

- `tier`: `:smoke` (seconds; in CI), `:standard` (minutes), `:full` (hours).
- `doc`: where the case is described (e.g. "docs/TEST_CASES.md §5b").
- `analyse(ctx) -> Dict{job => Dict{metric => value}}`; `ctx` has `dir`,
  `runs` (job name → `FA.Run`), `jobs` (job name → `Job`).
- `rtol`/`atol`: regression tolerance per metric (`"*"` = default for all).
"""
Base.@kwdef struct Case
    name::String
    title::String
    doc::String = ""
    tier::Symbol = :standard
    description::String = ""
    prep::Vector{Prep} = Prep[]
    jobs::Vector{Job} = Job[]
    analyse::Function = ctx -> Dict{String,Dict{String,Float64}}()
    checks::Vector{Check} = Check[]
    rtol::Dict{String,Float64} = Dict("*" => 1e-6)
    atol::Dict{String,Float64} = Dict("*" => 1e-12)
end

const REGISTRY = Dict{String,Case}()
const ORDER = String[]

"Add a case to the registry (a later definition with the same name replaces it)."
function register!(c::Case)
    c.tier in (:smoke, :standard, :full) || error("$(c.name): tier must be :smoke, :standard or :full")
    names = [j.name for j in c.jobs]
    allunique(names) || error("$(c.name): job names are not unique")
    haskey(REGISTRY, c.name) || push!(ORDER, c.name)
    REGISTRY[c.name] = c
    return c
end

"Registered cases, in definition order (optionally of one tier)."
cases(; tier = nothing) = [REGISTRY[n] for n in ORDER if tier === nothing || REGISTRY[n].tier == tier]
case(name::AbstractString) = haskey(REGISTRY, name) ? REGISTRY[name] :
    error("no case \"$name\"; known: $(join(ORDER, ", "))")

case_dir(c::Case) = joinpath(case_root(), c.name)
function _subst(s::AbstractString, c::Case)
    occursin('{', s) || return String(s)
    t = replace(s, "{case}" => case_dir(c), "{runs}" => joinpath(case_dir(c), "runs"), "{repo}" => REPO)
    # one separator style: compound values ("lat,lon,{case}/q.csv") keep their commas
    return join((occursin('{', p) ? normpath(q) : q for (p, q) in zip(split(s, ','), split(t, ','))), ",")
end
_subst(v::Vector{String}, c::Case) = [_subst(s, c) for s in v]
job_output(c::Case, j::Job) = joinpath(case_dir(c), "runs", j.name * ".h5")
job_log(c::Case, j::Job) = joinpath(case_dir(c), "runs", j.name * ".log")
"The full FloodModel.jl argument list of a job, as recorded in its output."
job_argv(c::Case, j::Job) = vcat(_subst(j.args, c), ["--output", job_output(c, j)])

# ───────────────────────────────── commands ────────────────────────────────

const _PYTHON = Ref{String}("")
function _python()
    isempty(_PYTHON[]) || return _PYTHON[]
    out = readchomp(`$(Base.julia_cmd()) --project=$(PROJECT[]) --startup-file=no -e "using PyCall; print(PyCall.python)"`)
    _PYTHON[] = out
end

function _env(extra_path = String[])
    env = copy(ENV)
    env["FLOODA5_NO_VIS"] = "1"
    if !isempty(extra_path)
        sep = Sys.iswindows() ? ";" : ":"
        env["PATH"] = join(vcat(extra_path, get(ENV, "PATH", "")), sep)
    end
    return env
end

function _cmd(c::Case, p::Prep)
    args = _subst(p.args, c)
    jl = Base.julia_cmd()
    if p.kind === :floodmodel
        return setenv(`$jl --project=$(PROJECT[]) --threads $(p.threads) $(joinpath(REPO, "FloodModel.jl")) $args`, _env()), REPO
    elseif p.kind === :julia
        return setenv(`$jl --project=$(PROJECT[]) --threads $(p.threads) $args`, _env()), REPO
    elseif p.kind === :python
        py = _python(); d = dirname(py)
        return setenv(`$py $args`, _env([d, joinpath(d, "Library", "bin"), joinpath(d, "Scripts")])), REPO
    end
    error("$(c.name)/$(p.name): unknown step kind $(p.kind)")
end

_cmd(c::Case, j::Job; threads = j.threads) = setenv(`$(Base.julia_cmd()) --project=$(PROJECT[]) --threads $threads $(joinpath(REPO, "FloodModel.jl")) $(job_argv(c, j))`, _env())

_q(a) = occursin(' ', a) ? "\"$a\"" : a
"How a step is shown (the command as a user would type it, from the repository)."
function _show_prep(c::Case, p::Prep)
    a = join(_q.(_subst(p.args, c)), " ")
    p.kind === :floodmodel && return "julia --project=. --threads $(p.threads) FloodModel.jl $a"
    p.kind === :julia && return "julia --project=. --threads $(p.threads) $a"
    p.kind === :write && return "write $(_subst(p.args[1], c)):\n        " *
                                replace(rstrip(p.args[2]), "\n" => "\n        ")
    # FloodA5's Python is PyCall's (Julia's Conda env), with conda's Library/bin on PATH
    return "julia --project=. -e 'using PyCall; d = dirname(PyCall.python); " *
           "withenv(\"PATH\" => join([d, joinpath(d, \"Library\", \"bin\"), ENV[\"PATH\"]], Sys.iswindows() ? \";\" : \":\")) do; " *
           "run(`\$(PyCall.python) $a`); end'"
end
_show_job(c::Case, j::Job, threads = j.threads) =
    "julia --project=. --threads $threads FloodModel.jl " * join(_q.(job_argv(c, j)), " ")

# ─────────────────────────────── preparation ───────────────────────────────

"""
    prepare(c; force = false, dry_run = false, io = stdout)

Run the case's preparation steps whose outputs are missing (all with `force`).
"""
function prepare(c::Case; force::Bool = false, dry_run::Bool = false, io::IO = stdout)
    mkpath(joinpath(case_dir(c), "runs"))
    for p in c.prep
        outs = _subst(p.outputs, c)
        if !force && !isempty(outs) && all(isfile, outs)
            println(io, "  prep $(p.name): up to date")
            continue
        end
        if p.kind === :write                     # args = [path, text]
            println(io, "  prep $(p.name): write ", _subst(p.args[1], c))
            dry_run && continue
            mkpath(dirname(_subst(p.args[1], c)))
            write(_subst(p.args[1], c), p.args[2])
            continue
        end
        cmd, dir = _cmd(c, p)
        println(io, "  prep $(p.name): ", _show_prep(c, p))
        dry_run && continue
        log = joinpath(case_dir(c), "prep_$(p.name).log")
        t0 = time()
        ok = success(pipeline(Cmd(cmd; dir), stdout = log, stderr = log))
        ok || error("$(c.name): preparation step $(p.name) failed — see $log")
        missing_out = filter(!isfile, outs)
        isempty(missing_out) || error("$(c.name)/$(p.name) did not produce $(missing_out)")
        @printf(io, "  prep %s: done (%.0f s)\n", p.name, time() - t0)
    end
end

# ────────────────────────────────── runs ───────────────────────────────────

"""
    job_status(c, j) -> :done | :missing | :failed | :changed

`:done` if the output completed and records exactly this job's command line.
"""
function job_status(c::Case, j::Job)
    p = job_output(c, j)
    isfile(p) || return :missing
    try
        FO.format_version(p) >= 2 || return :changed
        ri = FO.run_info(p)
        ri["argv"] == job_argv(c, j) || return :changed
        get(get(ri, "summary", Dict()), "completed", false) == true || return :failed
        return :done
    catch
        return :failed
    end
end

"""
    run_jobs(c; select = nothing, parallel = 1, threads = nothing,
             fresh = false, dry_run = false, io = stdout) -> Bool

Run the jobs that are not `:done` (all with `fresh`), `parallel` at a time.
`select`: job names or a predicate on `Job`. `threads` overrides each job's
thread count. Returns false if any run failed.
"""
function run_jobs(c::Case; select = nothing, parallel::Integer = 1, threads = nothing,
                  fresh::Bool = false, dry_run::Bool = false, io::IO = stdout)
    jobs = _selected(c, select)
    todo = Job[]
    for j in jobs
        st = fresh ? :fresh : job_status(c, j)
        st === :done ? println(io, "  run $(j.name): up to date") : push!(todo, j)
    end
    isempty(todo) && return true
    mkpath(joinpath(case_dir(c), "runs"))
    ok = Ref(true)
    lk = ReentrantLock()
    sem = Base.Semaphore(max(parallel, 1))
    @sync for j in todo
        cmd = _cmd(c, j; threads = something(threads, j.threads))
        if dry_run
            println(io, "  run $(j.name): ", _show_job(c, j, something(threads, j.threads)))
            continue
        end
        @async begin
            Base.acquire(sem)
            try
                t0 = time()
                log = job_log(c, j)
                good = success(pipeline(Cmd(cmd; dir = REPO), stdout = log, stderr = log))
                good &= job_status(c, j) === :done
                lock(lk) do
                    good || (ok[] = false)
                    @printf(io, "  run %s: %s (%.0f s)\n", j.name, good ? "done" : "FAILED — see $log", time() - t0)
                end
            finally
                Base.release(sem)
            end
        end
    end
    return ok[]
end

function _selected(c::Case, select)
    select === nothing && return c.jobs
    if select isa Function
        return filter(select, c.jobs)
    end
    names = String.(collect(select))
    bad = setdiff(names, [j.name for j in c.jobs])
    isempty(bad) || error("$(c.name): no job(s) $(join(bad, ", "))")
    return filter(j -> j.name in names, c.jobs)
end

# ───────────────────────────── analysis, checks ────────────────────────────

const Metrics = Dict{String,Dict{String,Float64}}

"""
    analyse(c; io = stdout) -> Metrics

Load the case's completed runs and compute its metrics; also writes
`{case}/metrics.toml`. Every job also gets, from its summary:
- `mass_balance_rel`: |volume added by the stepping − stored − removed| /
  volume added — conservation alone (0.4.2+ summaries; NaN before)
- `input_error_rel`: (volume added by the stepping − exact time integral of
  the inputs) / exact — the sources' first-order time integration
"""
function analyse(c::Case; io::IO = stdout)
    runs = Dict{String,Any}()
    for j in c.jobs
        job_status(c, j) === :done || error("$(c.name): run $(j.name) is not complete — run the case first")
        runs[j.name] = FA.load_run(job_output(c, j))
    end
    ctx = (; dir = case_dir(c), runs, jobs = Dict(j.name => j for j in c.jobs), case = c)
    m = Metrics()
    for j in c.jobs
        s = get(runs[j.name].info, "summary", Dict())
        num(k) = Float64(get(s, k, NaN))
        applied, added = num("vol_applied_m3"), num("vol_added_m3")
        m[j.name] = Dict("mass_balance_rel" => abs(num("conservation_error_m3")) / max(applied, eps()),
                         "input_error_rel" => (applied - added) / max(added, eps()))
    end
    for (job, d) in c.analyse(ctx)
        merge!(get!(m, String(job), Dict{String,Float64}()), Dict(String(k) => Float64(v) for (k, v) in d))
    end
    open(joinpath(case_dir(c), "metrics.toml"), "w") do f
        TOML.print(f, _toml_metrics(m); sorted = true)
    end
    return m
end

_toml_metrics(m::Metrics) = Dict(j => Dict(k => (isfinite(v) ? v : string(v)) for (k, v) in d) for (j, d) in m)
_from_toml(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)

expected_path(c::Case) = joinpath(REPO, "cases", "expected", c.name * ".toml")

_tol(d::Dict{String,Float64}, metric) = get(d, metric, get(d, "*", 0.0))

"""
    check(c, m = analyse(c); io = stdout) -> Bool

The case's physical checks, then every metric against the expected values
in `cases/expected/<case>.toml` (|Δ| ≤ atol + rtol·|expected|). Prints a
table; returns true if everything passes. A case without expected values
passes on its physical checks alone (with a note).
"""
function check(c::Case, m::Metrics = analyse(c); io::IO = stdout)
    ok = true
    for ch in c.checks
        pass = try
            ch.test(m) === true
        catch err
            println(io, "  check $(ch.name): ERROR ", sprint(showerror, err)); false
        end
        ok &= pass
        det = try ch.detail(m) catch; "" end
        println(io, "  ", pass ? "PASS" : "FAIL", "  ", ch.name, isempty(det) ? "" : "  — " * det)
    end
    ep = expected_path(c)
    if !isfile(ep)
        println(io, "  (no expected values: cases/expected/$(c.name).toml — `run.jl baseline $(c.name)`)")
        return ok
    end
    exp = TOML.parsefile(ep)
    meta = get(exp, "meta", Dict())
    nfail = 0; n = 0
    for (job, d) in exp
        job == "meta" && continue
        for (metric, ev) in d
            n += 1
            e = _from_toml(ev)
            v = get(get(m, job, Dict{String,Float64}()), metric, NaN)
            tol = _tol(c.atol, metric) + _tol(c.rtol, metric) * abs(e)
            pass = (isnan(e) && isnan(v)) || (isfinite(e) && abs(v - e) <= tol) || (isinf(e) && v == e)
            if !pass
                nfail += 1
                @printf(io, "  FAIL  %s / %s = %.10g, expected %.10g (tol %.3g)\n", job, metric, v, e, tol)
            end
        end
    end
    println(io, "  ", nfail == 0 ? "PASS" : "FAIL", "  $(n - nfail)/$n metrics match expected values ",
            "(recorded with FloodA5 $(get(meta, "flooda5_version", "?")), $(get(meta, "git_commit", "?")[1:min(end, 8)]))")
    return ok && nfail == 0
end

"""
    write_baseline(c, m = analyse(c); io = stdout)

Record the current metrics as the case's expected values, with the version,
commit and platform that produced them.
"""
function write_baseline(c::Case, m::Metrics = analyse(c); io::IO = stdout)
    r = FA.load_run(job_output(c, first(c.jobs)))
    meta = Dict("flooda5_version" => get(r.info, "flooda5_version", "?"),
                "git_commit" => get(r.info, "git_commit", "?"),
                "git_dirty" => get(r.info, "git_dirty", false),
                "recorded_utc" => Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SSZ"),
                "platform" => string(Sys.KERNEL, " ", Sys.ARCH))
    out = merge(Dict{String,Any}("meta" => meta), _toml_metrics(m))
    mkpath(dirname(expected_path(c)))
    open(expected_path(c), "w") do f
        println(f, "# Expected metrics for case \"$(c.name)\" — written by `julia --project=. cases/run.jl baseline $(c.name)`.")
        println(f, "# Tolerances are in cases/registry (Case.rtol / Case.atol).")
        TOML.print(f, out; sorted = true)
    end
    println(io, "  wrote ", relpath(expected_path(c), REPO))
end

# ──────────────────────────────── front end ────────────────────────────────

"""
    run_case(c; select, parallel, threads, fresh, dry_run, io) -> Bool

Prepare, run, analyse and check a case.
"""
function run_case(c::Case; select = nothing, parallel::Integer = 1, threads = nothing,
                  fresh::Bool = false, dry_run::Bool = false, io::IO = stdout)
    println(io, "== $(c.name): $(c.title)")
    prepare(c; dry_run, io)
    ok = run_jobs(c; select, parallel, threads, fresh, dry_run, io)
    dry_run && return true
    ok || (println(io, "  some runs failed: not analysed"); return false)
    select === nothing || return true          # partial: analyse once all jobs are done
    return check(c, analyse(c; io); io)
end

"Print a case: description, preparation commands and runs."
function show_case(c::Case; io::IO = stdout)
    println(io, "$(c.name) — $(c.title)   [tier: $(c.tier)]")
    isempty(c.doc) || println(io, "  documented in ", c.doc)
    isempty(c.description) || println(io, "\n", strip(c.description), "\n")
    println(io, "Directory: ", case_dir(c))
    for p in c.prep
        println(io, "prep $(p.name):\n    ", _show_prep(c, p))
    end
    for j in c.jobs
        println(io, "run $(j.name)", isempty(j.tags) ? "" : "  " * join(("$k=$v" for (k, v) in j.tags), " "),
                ":\n    julia --project=. --threads $(j.threads) FloodModel.jl ", join(job_argv(c, j), " "))
    end
    isempty(c.checks) || println(io, "checks: ", join((ch.name for ch in c.checks), "; "))
end

# ─────────────────────────────── the registry ──────────────────────────────

for f in sort(readdir(joinpath(@__DIR__, "registry")))
    endswith(f, ".jl") && include(joinpath(@__DIR__, "registry", f))
end

# ───────────────── public names for notebooks and scripts ─────────────────
# The notebooks (notebooks/*.ipynb) use these; the underscore helpers in the
# registry files may be renamed freely, these may not (test/test_cases.jl).

"FloodA5's Python (PyCall's), as the `:python` preparation steps run it."
python_command() = _python()
"Environment for running FloodA5's Python: conda's `Library/bin` etc. on PATH."
python_env(extra_path = String[]) = _env(extra_path)
"The planar DEM generator's parameters (`planar_embankment_params.json` in `dir`)."
planar_params(dir) = _planar_params(dir)
"(lat, lon) of the planar generator's injection point for the given domain."
planar_point(; kw...) = _planar_point(; kw...)
"Half-width (m) of the band about the source row excluded from the N/S split (planar_slope)."
const PLANAR_SLOPE_BAND = _PS_BAND
"Depth (m) above which a cell counts as wet in the planar_embankment metrics."
const EMBANKMENT_WET_DEPTH = _PE_WET
"x-range (m) of a cell polygon in the planar_embankment frame."
embankment_xrange(mesh, p) = _pe_xrange(mesh, p)
"(added, removed, domain, error, rel) of a run from its log (the model's step lines)."
log_mass_balance(logpath) = _case_log_mass_balance(logpath)
"The channel series table (rms/spread of Δη over orientations) from metric rows."
channel_series_table(rows) = _ch_series(rows)

export python_command, python_env, planar_params, planar_point, PLANAR_SLOPE_BAND,
       EMBANKMENT_WET_DEPTH, embankment_xrange, log_mass_balance, channel_series_table

end # module
