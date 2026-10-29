# test/diag_vs_reference.jl
# -------------------------
# R1 (2026-10): compare runs with a reference run on the same mesh (e.g. a
# local-inertial run against the `--flow-model swe` reference), frame by
# frame at matched simulation times, and at the final frame:
#   storage      total volume in the domain (Σ depth·area), and its ratio to
#                the reference
#   rms Δh       rms depth difference over cells wet (> 1 cm) in either run
#   CSI          flood-extent agreement, critical success index at a depth
#                threshold (default 0.1 m): hits / (hits + misses + false alarms)
#   bias         (wet in run) / (wet in reference) at the threshold
#   speed ≥ 3 m  mean speed in cells over 3 m deep, run / reference
#
#   FLOODA5_NO_VIS=1 julia --project=. test/diag_vs_reference.jl MESH.parquet REF.h5 RUN.h5 ... [--threshold 0.1] [--frame k]
#
# --frame k compares frame k (default: the last common frame); the frame
# times of each pair must agree within 60 s.

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using HDF5, Logging, Printf, Statistics

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const TH = parse(Float64, _arg("--threshold", "0.1"))
const FRAME = parse(Int, _arg("--frame", "0"))
files = filter(a -> endswith(a, ".h5"), ARGS)
mesh = with_logger(SimpleLogger(stderr, Logging.Error)) do; load_mesh_geoparquet(ARGS[1]); end
area = [A5Grid._polygon_area_m2(c.boundary) for c in mesh.cells]

function frames(p)
    h5open(p) do f
        ks = sort(collect(keys(f["frames"])))
        [(t = read(f["frames"][k]["t"]), h = read(f["frames"][k]["water_depth"]),
          u = read(f["frames"][k]["vel_u"]), v = read(f["frames"][k]["vel_v"])) for k in ks]
    end
end
function compare(R, X)
    h1 = R.h; h2 = X.h
    wet = [i for i in eachindex(h1) if isfinite(h1[i]) && isfinite(h2[i]) && (h1[i] > 0.01 || h2[i] > 0.01)]
    rms = sqrt(mean(abs2, h2[wet] .- h1[wet]))
    a = h1 .> TH; b = h2 .> TH
    hits = count(a .& b); miss = count(a .& .!b); fa = count(.!a .& b)
    deep = [i for i in eachindex(h1) if h1[i] >= 3.0]
    sp(F) = mean(hypot(F.u[i], F.v[i]) for i in deep)
    st(F) = sum(F.h[i] * area[i] for i in eachindex(F.h) if isfinite(F.h[i]))
    return (rms = rms, csi = hits / max(hits + miss + fa, 1), bias = count(b) / max(count(a), 1),
            stor = st(X) / 1e6, sratio = st(X) / st(R), speed = isempty(deep) ? NaN : sp(X) / sp(R))
end

ref = frames(files[1])
kr = FRAME > 0 ? FRAME : length(ref)
@printf("reference %s: storage %.3f Mm³ at t = %.0f s (frame %d); threshold %.2f m\n", basename(files[1]),
        sum(ref[kr].h[i] * area[i] for i in eachindex(area) if isfinite(ref[kr].h[i])) / 1e6, ref[kr].t, kr, TH)
@printf("%-34s %9s %8s %9s %7s %7s %10s\n", "run", "storage", "ratio", "rms Δh", "CSI", "bias", "speed≥3m")
for p in files[2:end]
    X = frames(p)
    k = FRAME > 0 ? FRAME : min(length(X), length(ref))
    abs(X[k].t - ref[k].t) < 60 || @warn "frame times differ: $(X[k].t) vs $(ref[k].t)"
    c = compare(ref[k], X[k])
    @printf("%-34s %9.3f %8.3f %8.3fm %7.3f %7.3f %10.3f\n", basename(p), c.stor, c.sratio, c.rms, c.csi, c.bias, c.speed)
end
