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
# A thin wrapper over analysis/FloodA5Analysis.jl (`FA.divergence`,
# `FA.growth_time`).

include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
using Printf

a, b = ARGS[1], ARGS[2]
label = length(ARGS) >= 3 ? ARGS[3] : basename(b)

rows = FA.divergence(FA.load_run(a), FA.load_run(b))

println("== ", label)
@printf("%8s %7s %11s %11s %7s %11s %11s\n", "t (s)", "wet", "rms dh (m)", "max dh (m)", ">1cm", "rms du", "dV/V")
for r in rows
    @printf("%8.0f %7d %11.3e %11.3e %7d %11.3e %11.3e\n", r.t, r.nwet, r.rms, r.mx, r.n1cm, r.urms, r.dV)
end
g = FA.growth_time(rows)
if g.n >= 2
    @printf("growth: e-folding time %.0f s (%d frames in the 1e-9…1e-3 m range)\n", g.tau, g.n)
else
    println("growth: not measurable (rms outside 1e-9…1e-3 m in all but ", g.n, " frame(s))")
end
