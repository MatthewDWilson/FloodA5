#!/usr/bin/env julia
# test_bdy_load.jl
# Quick diagnostic: load a .bdy file and print parsed series summary.
# Usage:
#   julia --project=. test/test_bdy_load.jl test/carlisle/carlisle.bdy

include("../FloodModel.jl")

function call_read_timeseries(path::String)
    println("\nLoading: $path\n" * "="^60)
    reader = LisfloodBDYReader(path)
    series = read_timeseries(reader)
    println("\nParsed $(length(series)) series:\n")
    for (name, (t_s, Q)) in sort(collect(series))
        n      = length(t_s)
        t_min  = n > 0 ? t_s[1]   : NaN
        t_max  = n > 0 ? t_s[end] : NaN
        Q_min  = n > 0 ? minimum(Q) : NaN
        Q_max  = n > 0 ? maximum(Q) : NaN
        t_hr_min = t_min / 3600
        t_hr_max = t_max / 3600
        println("  Series: '$name'")
        println("    Steps : $n")
        println("    Time  : $(round(t_hr_min, digits=3)) – $(round(t_hr_max, digits=3)) hours " *
                "($(round(t_min, digits=1)) – $(round(t_max, digits=1)) s)")
        println("    Q     : $(round(Q_min, sigdigits=4)) – $(round(Q_max, sigdigits=4)) m³/s")
        println()
    end
end

path = length(ARGS) >= 1 ? ARGS[1] : "test/carlisle/carlisle.bdy"
call_read_timeseries(path)
