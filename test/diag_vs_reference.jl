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
#   FLOODA5_NO_VIS=1 julia --project=. test/diag_vs_reference.jl [MESH.parquet] REF.h5 RUN.h5 ... [--threshold 0.1] [--frame k]
#
# --frame k compares frame k (default: the last common frame); the frame
# times of each pair must agree within 60 s. MESH.parquet gives the cell
# areas of format-1 runs (FloodA5 ≤ 0.3.x); format-2 runs carry them, so it
# may be omitted. A thin wrapper over analysis/FloodA5Analysis.jl
# (`FA.compare_fields`, `FA.compare`).

include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
using Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const TH = parse(Float64, _arg("--threshold", "0.1"))
const FRAME = parse(Int, _arg("--frame", "0"))
files = filter(a -> endswith(a, ".h5"), ARGS)
const MESH = (i = findfirst(a -> endswith(a, ".parquet"), ARGS); i === nothing ? nothing : ARGS[i])

# Format-1 runs on a re-labelled copy of the mesh (e.g. pya5 0.x ids in the
# output, 1.x in the parquet) are matched to the parquet by row, as before.
function _load(p)
    MESH === nothing && return FA.load_run(p)
    try
        return FA.load_run(p; parquet = MESH)
    catch e
        (e isa ErrorException && occursin("re-labelled", e.msg)) || rethrow()
        @warn "$(basename(p)): cell ids are not those of $(basename(MESH)); matching cells by row"
        return FA.load_run(p; parquet = MESH, match = :row)
    end
end

ref = _load(files[1])
isempty(ref.mesh.area) && error("$(files[1]) is a format-1 file: pass the mesh's parquet first")
kr = FRAME > 0 ? FRAME : FA.nframes(ref)
@printf("reference %s: storage %.3f Mm³ at t = %.0f s (frame %d); threshold %.2f m\n", basename(files[1]),
        FA.storage(ref, kr) / 1e6, ref.t[kr], kr, TH)
@printf("%-34s %9s %8s %9s %7s %7s %10s\n", "run", "storage", "ratio", "rms Δh", "CSI", "bias", "speed≥3m")
for p in files[2:end]
    X = _load(p)
    k = FRAME > 0 ? FRAME : min(FA.nframes(X), FA.nframes(ref))
    abs(X.t[k] - ref.t[k]) < 60 || @warn "frame times differ: $(X.t[k]) vs $(ref.t[k])"
    c = FA.compare(ref, X; k_ref = k, k_run = k, threshold = TH)
    @printf("%-34s %9.3f %8.3f %8.3fm %7.3f %7.3f %10.3f\n", basename(p), c.storage / 1e6, c.storage_ratio,
            c.rms, c.csi, c.bias, c.speed_ratio)
end
