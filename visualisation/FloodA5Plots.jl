# visualisation/FloodA5Plots.jl
# =============================
# Static figures of FloodA5 runs with CairoMakie, for notebooks and papers.
#
# Works from output files alone: it loads analysis/FloodA5Analysis.jl (FA),
# not FloodModel.jl, and needs no Python. Without CairoMakie (the CI
# environment) the module still loads and its Makie-free logic works
# (colour ranges, frame selection, extent classes, map coordinates, sources:
# plots/core.jl); the plotting functions are then not defined.
#
#   include("visualisation/FloodA5Plots.jl"); const FP = FloodA5Plots
#   const FA = FP.FA
#   li  = FA.load_run("test/case_runs/carlisle/runs/a5_li.h5")
#   swe = FA.load_run("test/case_runs/carlisle/runs/a5_swe.h5")
#   FP.fieldmap(li; t = 4 * 3600)                     # one map → Figure
#   FP.fieldmaps([li, swe], [2, 4, 6] .* 3600)        # runs × times grid
#   FP.difference_map(swe, li; t = 6 * 3600)          # Δh, diverging
#   FP.extent_map(swe, li; t = 6 * 3600)              # hits/misses/false alarms, CSI
#   FP.storage_plot([li, swe]); FP.mass_balance_plot(li)
#   FP.timeseries([li, swe], [(-2.93, 54.895)])
#   save("fig.png", fig; px_per_unit = 2)
#
# Every `name(...)` returning a Figure has a `name!(gridposition, ...)` that
# draws into a figure you lay out yourself (map functions return a
# `MapPanel`, which `mark_points!` uses to place points and labels).
#
# Style: `house_theme()` (also used by test/channel_bench/plot_channel.jl) —
# INK/INK2/MUTED text, SURF background, SERIES colours in fixed order, A5
# blue / Cartesian orange (MESH_STYLE). FLOODA5_PLOTS_NO_MAKIE=1 skips
# CairoMakie even where it is installed.
#
# Files: plots/core.jl (Makie-free), plots/makie.jl (figures).

module FloodA5Plots

using Statistics, Printf

if isdefined(Main, :FloodA5Analysis)
    const FA = Main.FloodA5Analysis
else
    include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
    const FA = FloodA5Analysis
end

include(joinpath("plots", "core.jl"))

"True if CairoMakie was loaded and the plotting functions are defined."
const HAS_MAKIE = get(ENV, "FLOODA5_PLOTS_NO_MAKIE", "") != "1" &&
                  Base.find_package("CairoMakie") !== nothing

if HAS_MAKIE
    using CairoMakie
    include(joinpath("plots", "makie.jl"))
end

end # module
