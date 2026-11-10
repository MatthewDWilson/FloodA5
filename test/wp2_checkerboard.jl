# test/wp2_checkerboard.jl
# ------------------------
# WP2 acceptance diagnostic: checkerboarding and noise in FloodA5 HDF5 output.
# Analysis only; a thin wrapper over analysis/FloodA5Analysis.jl (noise.jl).
#
#   FLOODA5_NO_VIS=1 julia --project=. test/wp2_checkerboard.jl [MESH.parquet] [--every N] [--series] [--sgs] RUN.h5 ...
#
# MESH: the runs' mesh, needed for format-1 output (FloodA5 ≤ 0.3.x); format-2
# files carry their geometry and do not use it. --sgs: SGS output, where
# water_depth is measured from the cell's lowest DEM point, so WSE =
# sgs_z_min + depth instead of the mean elevation + depth (format-2 runs
# record their flow model; --sgs is implied for SGS ones). If a format-1 run's
# cell IDs do not match MESH (e.g. a run on the pya5 0.x copy of a migrated
# mesh), cells are matched by row order.
#
# Over "fully wet interior" cells (depth > 5 cm, and every neighbour too):
#   r_rms    rms of r_i = WSE_i − mean(WSE of neighbours), mm. Smooth fields
#            give O(dx²·curvature); a checkerboard gives |r| ~ its amplitude.
#   extrema  fraction of those cells that are strict local WSE extrema with
#            |r_i| > 2 mm. ~0 for a smooth field, large for a checkerboard
#            (1 for a perfect two-colour one on a square grid).
#   dv_rms   rms of |v_i − mean(v of neighbours)|, m/s (velocity noise).
# Reported for the last frame and as the maximum over all frames (every Nth).

include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
using Printf

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])

mesh_path = (endswith(ARGS[1], ".h5") || startswith(ARGS[1], "--")) ? nothing : ARGS[1]
every     = parse(Int, _arg("--every", "1"))
paths     = filter(a -> endswith(a, ".h5"), ARGS)

"Load a run; format 1 takes its geometry from MESH, by cell id or else by row order."
function load(p)
    FA.FO.format_version(p) >= 2 && return FA.load_run(p)
    mesh_path === nothing && error("$p is a format-1 file: give the mesh's parquet first")
    try
        return FA.load_run(p; parquet = mesh_path, match = :id)
    catch err
        occursin("cell ids", sprint(showerror, err)) || rethrow()
        return FA.load_run(p; parquet = mesh_path, match = :row)
    end
end

@printf("%-44s %7s | %6s %8s %8s %8s | %8s %8s %8s %7s\n", "run", "t_end", "n",
        "r_rms", "extrema", "dv_rms", "max r", "max ext", "max dv", "at t")
for p in paths
    r = load(p)
    series = FA.noise_series(r; every, sgs = "--sgs" in ARGS || r.flow_model == "sgs")
    if "--series" in ARGS
        for s in series
            @printf("  %-40s t=%7.0f n=%6d r_rms=%7.2fmm extrema=%.4f dv_rms=%.4f\n",
                    basename(p), s.t, s.m.n, s.m.r_rms, s.m.extrema, s.m.dv_rms)
        end
    end
    sm = FA.noise_summary(series)
    last = sm.last
    @printf("%-44s %7.0f | %6d %7.2fmm %8.4f %8.4f | %6.2fmm %8.4f %8.4f %7.0f\n",
            basename(p), last.t, last.m.n, last.m.r_rms, last.m.extrema, last.m.dv_rms,
            sm.max_r_rms, sm.max_extrema, sm.max_dv_rms, sm.t_max_r_rms)
end
