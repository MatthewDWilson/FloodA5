# cases/run.jl — command-line front end of the test-case registry.
#
#   julia --project=. cases/run.jl list [--tier smoke|standard|full]
#   julia --project=. cases/run.jl show CASE
#   julia --project=. cases/run.jl prep CASE… [--force] [--dry-run]   # inputs only
#   julia --project=. cases/run.jl run CASE… [--jobs a,b] [-P N] [-T N] [--fresh] [--dry-run]
#   julia --project=. cases/run.jl check CASE…        # analyse existing runs, check
#   julia --project=. cases/run.jl baseline CASE…     # record current metrics as expected
#
# CASE may also be a tier (smoke, standard, full) for every case in it.
# -P: runs in parallel (default 1); -T: threads per run (default: per job).
# Outputs: test/case_runs/<case>/ (or $FLOODA5_CASES_DIR/<case>/).
# Exit code 1 if a run or a check fails.

using Printf
include(joinpath(@__DIR__, "FloodA5Cases.jl"))
const FC = FloodA5Cases

function _opt(args, flag, default)
    i = findfirst(==(flag), args)
    i === nothing && return default
    v = args[i+1]; deleteat!(args, i:i+1)
    return v
end
function _flag(args, flag)
    i = findfirst(==(flag), args)
    i === nothing && return false
    deleteat!(args, i); return true
end

function _cases(names)
    out = FC.Case[]
    for n in names
        if n in ("smoke", "standard", "full")
            append!(out, FC.cases(; tier = Symbol(n)))
        else
            push!(out, FC.case(n))
        end
    end
    isempty(out) && error("no case given")
    return out
end

function main(argv)
    args = String.(copy(argv))
    isempty(args) && (println(read(@__FILE__, String)[1:findfirst("\n\ninclude", read(@__FILE__, String))[1]]); return 0)
    cmd = popfirst!(args)
    if cmd == "list"
        tier = _opt(args, "--tier", nothing)
        for c in FC.cases(; tier = tier === nothing ? nothing : Symbol(tier))
            st = isempty(c.jobs) ? "" : string(count(j -> FC.job_status(c, j) === :done, c.jobs), "/", length(c.jobs), " runs done")
            @printf("%-24s %-9s %-16s %s\n", c.name, c.tier, st, c.title)
        end
        return 0
    elseif cmd == "show"
        foreach(c -> (FC.show_case(c); println()), _cases(args))
        return 0
    elseif cmd == "run"
        jobs = _opt(args, "--jobs", nothing)
        P = parse(Int, _opt(args, "-P", "1"))
        T = _opt(args, "-T", nothing)
        fresh = _flag(args, "--fresh"); dry = _flag(args, "--dry-run")
        ok = true
        for c in _cases(args)
            ok &= FC.run_case(c; select = jobs === nothing ? nothing : split(jobs, ","),
                              parallel = P, threads = T === nothing ? nothing : parse(Int, T),
                              fresh, dry_run = dry)
        end
        return ok ? 0 : 1
    elseif cmd == "prep"
        force = _flag(args, "--force"); dry = _flag(args, "--dry-run")
        for c in _cases(args)
            println("== $(c.name)")
            FC.prepare(c; force, dry_run = dry)
        end
        return 0
    elseif cmd == "check"
        ok = true
        for c in _cases(args)
            println("== $(c.name)")
            ok &= FC.check(c, FC.analyse(c))
        end
        return ok ? 0 : 1
    elseif cmd == "baseline"
        for c in _cases(args)
            println("== $(c.name)")
            FC.write_baseline(c, FC.analyse(c))
        end
        return 0
    end
    error("unknown command \"$cmd\" (list, show, prep, run, check, baseline)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
