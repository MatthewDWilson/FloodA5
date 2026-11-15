# benchmarks/run.jl — run the performance suite on this machine (benchmarks/BenchSuite.jl).
#
#   julia --project=. benchmarks/run.jl                       # the sections this machine has no results for
#   julia --project=. benchmarks/run.jl --force               # all of them again
#   julia --project=. benchmarks/run.jl --sections step,threads
#   julia --project=. benchmarks/run.jl --quick               # smoke test: small, not a benchmark
#   julia --project=. benchmarks/run.jl --label NAME          # name this machine (default: host name + OS)
#   julia --project=. benchmarks/run.jl --no-large            # leave out the 5 m Carlisle mesh
#   julia --project=. benchmarks/run.jl --list                # the machines with results
#
# Results: benchmarks/results/<machine>.json (commit it). Use a machine that is
# otherwise idle. The full suite takes about two hours on a 16-thread
# workstation with a GPU; figures: notebooks/05_performance.ipynb.

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "BenchSuite.jl"))
const BS = BenchSuite

function main(args)
    opt(name, default) = (i = findfirst(==(name), args); i === nothing ? default : args[i + 1])
    if "--help" in args || "-h" in args
        foreach(l -> startswith(l, "#") ? println(lstrip(l, ['#', ' '])) : nothing,
                Iterators.takewhile(l -> startswith(l, "#"), eachline(@__FILE__)))
        return
    end
    if "--list" in args
        for r in BS.load_results()
            m = r["machine"]
            println(rpad(m["id"], 28), m["cpu"], ", ", m["threads"], " threads; ",
                    m["gpu"] == "" ? "no GPU" : m["gpu"], "; sections: ", join(sort(collect(keys(r["sections"]))), ", "))
        end
        return
    end
    BS.run_suite(; sections = Symbol.(split(opt("--sections", join(BS.SECTIONS, ",")), ",")),
                 force = "--force" in args, quick = "--quick" in args,
                 label = opt("--label", ""), large = !("--no-large" in args))
    return
end

main(ARGS)
