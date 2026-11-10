# test/wp1_metrics.jl
# -------------------
# WP1 analysis of FloodA5 HDF5 output with the review's §3.4 metrics
# (FloodA5_Review_and_ForwardPlan_2026-10.md). Analysis only: a thin CLI over
# analysis/shape.jl (FloodA5Analysis, `FA.spread_metrics`, `FA.planar_metrics`);
# reads cell centres and the water_depth / volume frames (output formats 1 and 2).
#
# Point spread (flat bed, central source):
#   julia --project=. test/wp1_metrics.jl spread --t 3600 --lon 0 --lat 0 FILE.h5 ...
#   A4      Σ h r² cos 4φ / Σ h r²      (< 0: diagonal elongation)
#   A6      the same with cos 6φ       (for hexagonal meshes)
#   R(h>x)  cardinal/diagonal ratio of the 8-sector maximum radius at depth > x
#
# Planar slope (north/south asymmetry about the source latitude):
#   julia --project=. test/wp1_metrics.jl planar --lat 51.0001 [--t T] FILE.h5 ...
#   asym    (V_N − V_S)/(V_N + V_S), cells within ±50 m of the source
#           latitude excluded (as test/test_planar_symmetry.jl)
#   A2      Σ h r² sin 2φ / Σ h r² about the source — signed lateral skew of
#           the plume relative to the downslope (east) axis
#   Reports the frame nearest --t (default: last), plus max|asym| over all
#   frames and the time it occurs.
#
# The frame is chosen as the one nearest the requested time, and its time is
# printed: compare runs at matched simulation times (review §3.3).
#
# Coordinates are about the CENTRE OF THE SOURCE CELL — the mesh cell nearest
# (lon0, lat0), which is where FloodA5 injects (`snap = true`). Using the
# nominal source point instead adds a fixed offset (the cell centre is up to
# half a cell away) that reads as a large deflection while the plume is small.
# This script's Earth radius is 6_371_008.8 m (FA's default is R_EARTH =
# 6_371_000.0 m); it is passed explicitly so the printed numbers are unchanged.

using Printf
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis

const _R_EARTH = 6_371_008.8

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])

function spread(paths; t_target, lon0, lat0)
    @printf("%-48s %8s %7s %8s %8s %9s %9s %9s\n", "file", "t (s)", "n_wet",
            "A4", "A6", "R(>1e-4)", "R(>1cm)", "R(>5cm)")
    for p in paths
        r = FA.load_run(p)
        m = FA.spread_metrics(r, FA.frame_at(r, t_target); lon0, lat0, R = _R_EARTH, snap = true)
        @printf("%-48s %8.0f %7d %8.4f %8.4f %9.4f %9.4f %9.4f\n", basename(p),
                m.t, m.n_wet, m.A4, m.A6, m.ratio...)
    end
end

function planar(paths; t_target, lon0, lat0, t_sweep = nothing)
    if t_sweep === nothing
        @printf("%-52s %8s %8s %8s %8s   %9s %8s\n", "file", "t (s)", "asym", "A2",
                "defl°", "max|asym|", "(t>=3600)")
    else
        @printf("%-46s %s\n", "file", join([@sprintf("  t=%-6.0f asym/defl°", t) for t in t_sweep]))
    end
    pm(r, k) = FA.planar_metrics(r, k; lon0, lat0, R = _R_EARTH, snap = true)
    for p in paths
        r = FA.load_run(p)
        if t_sweep !== nothing
            cols = [(m = pm(r, FA.frame_at(r, t)); @sprintf("  %+7.3f / %+6.2f  ", m.asym, m.deflection))
                    for t in t_sweep]
            @printf("%-46s %s\n", basename(p), join(cols))
            continue
        end
        late   = findall(>=(3600.0), FA.times(r))      # skip the first hour (few wet cells)
        series = [pm(r, j).asym for j in late]
        m = pm(r, t_target === nothing ? FA.nframes(r) : FA.frame_at(r, t_target))
        imax = argmax(abs.(series))
        @printf("%-52s %8.0f %+8.4f %+8.4f %+8.2f   %9.4f %8.0f\n", basename(p),
                m.t, m.asym, m.A2, m.deflection, abs(series[imax]), r.t[late[imax]])
    end
end

mode  = ARGS[1]
paths = filter(a -> endswith(a, ".h5"), ARGS)
t_str = _arg("--t", nothing)
lon0  = parse(Float64, _arg("--lon", "0.0"))
lat0  = parse(Float64, _arg("--lat", "0.0"))
if mode == "spread"
    spread(paths; t_target = t_str === nothing ? Inf : parse(Float64, t_str), lon0, lat0)
elseif mode == "planar"
    sw = _arg("--sweep", nothing)                  # e.g. --sweep 7200,36000,72000
    planar(paths; t_target = t_str === nothing ? nothing : parse(Float64, t_str), lon0, lat0,
           t_sweep = sw === nothing ? nothing : parse.(Float64, split(sw, ",")))
else
    error("mode must be 'spread' or 'planar'")
end
