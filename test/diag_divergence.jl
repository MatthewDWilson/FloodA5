# test/diag_divergence.jl
# -----------------------
# Damping diagnostics (D1, 2026-10): how two runs that should agree drift
# apart, frame by frame — for twin runs (a perturbation of ~1e-10, e.g. the
# rainfall rate, or `--cell-order spatial` vs `file`) and for timestep
# convergence (the same case at different `--courant`).
#
#   FLOODA5_NO_VIS=1 julia --project=. test/diag_divergence.jl A.h5 B.h5 [label]
#
# Both files must cover the same mesh in the same (file) cell order and write
# frames at the same times. Per frame, over cells wet (> 5 cm) in either run:
#   rms / max depth difference, cells differing by > 1 cm, rms speed
#   difference, and the relative difference in total stored volume.
# For twin runs, `growth` is the e-folding time of the rms depth difference
# while it is still small (between 1e-9 m and 1e-3 m), from a log-linear fit.

using HDF5, Printf, Statistics

a, b = ARGS[1], ARGS[2]
label = length(ARGS) >= 3 ? ARGS[3] : basename(b)

rows = NamedTuple[]
h5open(a) do A; h5open(b) do B
    ida = read(A["mesh/cell_ids"]); idb = read(B["mesh/cell_ids"])
    ida == idb || error("the two files have different cell order or meshes")
    ka = sort(collect(keys(A["frames"])), by = x -> parse(Int, x))
    kb = Set(keys(B["frames"]))
    for k in ka
        k in kb || continue
        ga = A["frames"][k]; gb = B["frames"][k]
        t  = read(ga["t"])
        ha = read(ga["water_depth"]); hb = read(gb["water_depth"])
        ua = read(ga["velocity"]);    ub = read(gb["velocity"])
        va = read(ga["volume"]);      vb = read(gb["volume"])
        wet = findall(i -> ha[i] > 0.05 || hb[i] > 0.05, eachindex(ha))
        dh = isempty(wet) ? [0.0] : abs.(ha[wet] .- hb[wet])
        du = isempty(wet) ? [0.0] : abs.(ua[wet] .- ub[wet])
        Va = sum(filter(isfinite, va)); Vb = sum(filter(isfinite, vb))
        push!(rows, (t = t, nwet = length(wet), rms = sqrt(mean(dh .^ 2)), mx = maximum(dh),
                     n1cm = count(>(0.01), dh), urms = sqrt(mean(du .^ 2)),
                     dV = Va > 0 ? (Vb - Va) / Va : 0.0))
    end
end end

println("== ", label)
@printf("%8s %7s %11s %11s %7s %11s %11s\n", "t (s)", "wet", "rms dh (m)", "max dh (m)", ">1cm", "rms du", "dV/V")
for r in rows
    @printf("%8.0f %7d %11.3e %11.3e %7d %11.3e %11.3e\n", r.t, r.nwet, r.rms, r.mx, r.n1cm, r.urms, r.dV)
end
sel = [r for r in rows if 1e-9 < r.rms < 1e-3]
if length(sel) >= 2
    x = [r.t for r in sel]; y = log.([r.rms for r in sel])
    slope = sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
    @printf("growth: e-folding time %.0f s (%d frames in the 1e-9…1e-3 m range)\n", 1 / slope, length(sel))
else
    println("growth: not measurable (rms outside 1e-9…1e-3 m in all but ", length(sel), " frame(s))")
end
