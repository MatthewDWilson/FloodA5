# test/channel_bench/plot_channel.jl
# ----------------------------------
# Channel-representation series: plots from analyse_channel.jl's metrics.csv.
#   channel_vs_resolution.png  stage error and floodplain spill against cell
#                              size, per mesh (median over orientations and
#                              offsets; band = their range), SWE and
#                              local-inertial side by side
#   channel_vs_angle.png       stage error against channel orientation, per
#                              resolution and mesh (SWE and local-inertial)
#
#   julia --project=. test/channel_bench/plot_channel.jl RUNDIR
# (needs CairoMakie: the full project environment)

using CairoMakie, Statistics, Printf

const RUNDIR = ARGS[1]
const INK, INK2, MUTED, GRID, SURF = colorant"#0b0b0b", colorant"#52514e", colorant"#898781", colorant"#e1e0d9", colorant"#fcfcfb"
const MESHES = [("a5", "A5", colorant"#2a78d6", :circle), ("cart", "Cartesian", colorant"#eb6834", :rect)]
const SOLVERS = [("swe", "Full shallow-water (SWE)"), ("li", "Local-inertial (0.3.0 default scheme)")]

lines_ = readlines(joinpath(RUNDIR, "metrics.csv"))
hdr = split(lines_[1], ","); col(k) = findfirst(==(k), hdr)
rows = [split(l, ",") for l in lines_[2:end]]
num(r, k) = parse(Float64, r[col(k)])
get_(mesh, solver) = [r for r in rows if r[col("mesh")] == mesh && r[col("solver")] == solver]

set_theme!(Theme(fontsize = 13, textcolor = INK, backgroundcolor = SURF,
                 Axis = (backgroundcolor = SURF, xgridcolor = GRID, ygridcolor = GRID,
                         leftspinecolor = MUTED, bottomspinecolor = MUTED, topspinevisible = false,
                         rightspinevisible = false, xtickcolor = MUTED, ytickcolor = MUTED,
                         xticklabelcolor = INK2, yticklabelcolor = INK2, titlealign = :left,
                         titlefont = :bold)))

# ── metric against cell size ────────────────────────────────────────────────
fig = Figure(size = (980, 640))
for (c, (sv, stitle)) in enumerate(SOLVERS), (r, (metric, ylab)) in enumerate((("deta", "Stage error Δη (m)"), ("spill", "Floodplain spill width (m)")))
    ax = Axis(fig[r, c]; xscale = log2, xlabel = r == 2 ? "Cell size (m, √mean area)" : "",
              ylabel = c == 1 ? ylab : "", title = r == 1 ? stitle : "",
              xticks = ([5.5, 11.1, 22.2, 44.4], ["5.5", "11", "22", "44"]), xreversed = true)
    hlines!(ax, [0.0]; color = MUTED, linestyle = :dash, linewidth = 1)
    for (m, label, colr, mk) in MESHES
        rs = get_(m, sv); isempty(rs) && continue
        cells = sort(unique(num.(rs, "cell_m")))
        vals = [filter(isfinite, [num(x, metric) for x in rs if num(x, "cell_m") == cs]) for cs in cells]
        ok = .!isempty.(vals); cells = cells[ok]; vals = vals[ok]
        med = median.(vals); lo = minimum.(vals); hi = maximum.(vals)
        band!(ax, cells, lo, hi; color = (colr, 0.18))
        lines!(ax, cells, med; color = colr, linewidth = 2)
        scatter!(ax, cells, med; color = colr, marker = mk, markersize = 9, strokecolor = SURF, strokewidth = 1.5)
        up = m == "a5"                        # A5 labelled above its line, Cartesian below
        text!(ax, cells[1], med[1]; text = label, color = INK2, align = (:right, up ? :bottom : :top),
              offset = (-6, up ? 6 : -6), fontsize = 12)
    end
end
Legend(fig[3, 1:2], [[LineElement(color = colr, linewidth = 2), MarkerElement(color = colr, marker = mk, markersize = 9)]
                     for (_, _, colr, mk) in MESHES], [l for (_, l, _, _) in MESHES];
       orientation = :horizontal, framevisible = false, labelcolor = INK2)
Label(fig[0, 1:2], "Straight 20 m channel: error against resolution (median over 7 orientations × 2 offsets; band = range)";
      halign = :left, font = :bold, color = INK)
save(joinpath(RUNDIR, "channel_vs_resolution.png"), fig; px_per_unit = 2)

# ── stage error against orientation, one row per solver ─────────────────────
cellsz = sort(unique(num.(rows, "cell_m")))
fig2 = Figure(size = (240 * length(cellsz) + 60, 300 * length(SOLVERS) + 50))
for (row, (sv, stitle)) in enumerate(SOLVERS)
    rs_all = [r for r in rows if r[col("solver")] == sv]
    Label(fig2[2row - 1, 1:length(cellsz)], stitle; halign = :left, font = :bold, color = INK)
    for (k, cs) in enumerate(cellsz)
        ax = Axis(fig2[2row, k]; xlabel = row == length(SOLVERS) ? "Channel orientation (°)" : "",
                  ylabel = k == 1 ? "Stage error Δη (m)" : "", title = @sprintf("%.0f m cells", cs), xticks = 0:30:90)
        hlines!(ax, [0.0]; color = MUTED, linestyle = :dash, linewidth = 1)
        for (m, label, colr, mk) in MESHES
            rs = [r for r in rs_all if r[col("mesh")] == m && num(r, "cell_m") == cs]
            isempty(rs) && continue
            angs = sort(unique(num.(rs, "angle")))
            md = [mean(num(x, "deta") for x in rs if num(x, "angle") == a) for a in angs]
            lines!(ax, angs, md; color = colr, linewidth = 2)
            scatter!(ax, angs, md; color = colr, marker = mk, markersize = 9, strokecolor = SURF, strokewidth = 1.5)
        end
        ylims!(ax, -0.4, 1.5)
    end
end
Legend(fig2[2length(SOLVERS) + 1, 1:length(cellsz)], [[LineElement(color = colr, linewidth = 2), MarkerElement(color = colr, marker = mk, markersize = 9)]
                                   for (_, _, colr, mk) in MESHES], [l for (_, l, _, _) in MESHES];
       orientation = :horizontal, framevisible = false, labelcolor = INK2)
save(joinpath(RUNDIR, "channel_vs_angle.png"), fig2; px_per_unit = 2)
println("plots written to $RUNDIR")
