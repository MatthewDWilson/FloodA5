# test/wp2_checkerboard.jl
# ------------------------
# WP2 acceptance diagnostic: checkerboarding and noise in FloodA5 HDF5 output
# (standard flow: WSE = elevation + depth). Analysis only.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/wp2_checkerboard.jl MESH.parquet [--every N] [--series] RUN.h5 ...
#
# Over "fully wet interior" cells (depth > 5 cm, and every neighbour too):
#   r_rms    rms of r_i = WSE_i − mean(WSE of neighbours), mm. Smooth fields
#            give O(dx²·curvature); a checkerboard gives |r| ~ its amplitude.
#   extrema  fraction of those cells that are strict local WSE extrema with
#            |r_i| > 2 mm. ~0 for a smooth field, approaching 0.5 for a
#            checkerboard.
#   dv_rms   rms of |v_i − mean(v of neighbours)|, m/s (velocity noise).
# Reported for the last frame and as the maximum over all frames (every Nth).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using HDF5, Printf, Statistics

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])

function neighbour_lists(mesh_path, h5_ids)
    mesh = load_mesh_geoparquet(mesh_path)
    idx  = Dict(_norm_cell_id(id) => k for (k, id) in enumerate(h5_ids))
    nbrs = [Int[] for _ in h5_ids]
    for (id, nb) in mesh.adjacency
        i = get(idx, _norm_cell_id(id), 0)
        i == 0 && continue
        for n in nb
            j = get(idx, _norm_cell_id(n), 0)
            j > 0 && push!(nbrs[i], j)
        end
    end
    return nbrs
end

function frame_metrics(z, h, u, v, nbrs; hmin = 0.05)
    r = Float64[]; dv = Float64[]; nextr = 0
    for i in eachindex(h)
        (h[i] > hmin && isfinite(z[i]) && length(nbrs[i]) >= 3) || continue
        all(j -> h[j] > hmin && isfinite(z[j]), nbrs[i]) || continue
        wi  = z[i] + h[i]
        wn  = [z[j] + h[j] for j in nbrs[i]]
        ri  = wi - mean(wn)
        push!(r, ri)
        (abs(ri) > 2e-3 && (all(<(wi), wn) || all(>(wi), wn))) && (nextr += 1)
        push!(dv, hypot(u[i] - mean(u[nbrs[i]]), v[i] - mean(v[nbrs[i]])))
    end
    isempty(r) && return (n = 0, r_rms = NaN, extrema = NaN, dv_rms = NaN)
    return (n = length(r), r_rms = 1e3 * sqrt(mean(r .^ 2)), extrema = nextr / length(r),
            dv_rms = sqrt(mean(dv .^ 2)))
end

mesh_path = ARGS[1]
every     = parse(Int, _arg("--every", "1"))
paths     = filter(a -> endswith(a, ".h5"), ARGS)
nbrs      = nothing

@printf("%-44s %7s | %6s %8s %8s %8s | %8s %8s %8s %7s\n", "run", "t_end", "n",
        "r_rms", "extrema", "dv_rms", "max r", "max ext", "max dv", "at t")
for p in paths
    h5open(p) do f
        global nbrs
        nbrs === nothing && (nbrs = neighbour_lists(mesh_path, read(f["mesh/cell_ids"])))
        z  = read(f["mesh/elevations"])
        ks = sort(collect(keys(f["frames"])))
        series = []
        for (n, k) in enumerate(ks)
            (n % every == 0 || n == length(ks)) || continue
            g = f["frames"][k]
            m = frame_metrics(z, read(g["water_depth"]), read(g["vel_u"]), read(g["vel_v"]), nbrs)
            push!(series, (t = read(g["t"]), m = m))
        end
        if "--series" in ARGS
            for s in series
                @printf("  %-40s t=%7.0f n=%6d r_rms=%7.2fmm extrema=%.4f dv_rms=%.4f
",
                        basename(p), s.t, s.m.n, s.m.r_rms, s.m.extrema, s.m.dv_rms)
            end
        end
        last = series[end]
        valid = filter(s -> s.m.n > 0, series)
        imax = isempty(valid) ? 0 : argmax([s.m.r_rms for s in valid])
        @printf("%-44s %7.0f | %6d %7.2fmm %8.4f %8.4f | %6.2fmm %8.4f %8.4f %7.0f\n",
                basename(p), last.t, last.m.n, last.m.r_rms, last.m.extrema, last.m.dv_rms,
                imax > 0 ? valid[imax].m.r_rms : NaN,
                isempty(valid) ? NaN : maximum(s.m.extrema for s in valid),
                isempty(valid) ? NaN : maximum(s.m.dv_rms for s in valid),
                imax > 0 ? valid[imax].t : NaN)
    end
end
