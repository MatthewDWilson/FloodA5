# test/channel_bench/analyse_channel.jl
# -------------------------------------
# Channel-representation series: metrics from the runs listed in a
# make_channel_meshes.jl manifest. For every case and solver, from the last
# output frame, over the measurement reach s ∈ [−0.3 L, 0.05 L] (300 m below
# the inflow, 500 m above the steep reach):
#   deta     mean water-surface error (m) in the cells that overlap the
#            channel, against the 2-D normal-flow surface Z0 − S·s + h_n.
#            Lost conveyance shows as a positive error.
#   spill    wet (> 5 cm) floodplain area per metre of reach (m) in cells
#            wholly outside the channel; 0 for in-bank flow
#   unsteady max |Δη| over the channel cells between the last two frames (m)
#   outflow  discharge leaving the domain over the last frame interval,
#            as a fraction of the inflow (1 at steady state)
# Writes RUNDIR/metrics.csv.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/channel_bench/analyse_channel.jl MESHDIR RUNDIR [TAG]

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
using HDF5, Logging, Printf, Statistics

const MESHDIR, RUNDIR = ARGS[1], ARGS[2]
const TAG = length(ARGS) >= 3 ? ARGS[3] : "chan"
const L, W, D, S, SX, Z0, Q, NMAN = 2000.0, 20.0, 2.0, 0.001, 0.005, 10.0, 20.0, 0.03
const LON0, LAT0 = 0.0, 51.0
const R = A5Grid._EARTH_R
const HN = (Q * NMAN / (W * sqrt(S)))^(3 / 5)
xy(lon, lat) = (deg2rad(lon - LON0) * R * cosd(LAT0), deg2rad(lat - LAT0) * R)

function lastframes(p)
    h5open(p) do f
        ks = sort(collect(keys(f["frames"])), by = x -> parse(Int, x))
        g(k) = (t = read(f["frames"][k]["t"]), h = read(f["frames"][k]["water_depth"]))
        g(ks[end - 1]), g(ks[end])
    end
end
"(time, cumulative outflow) of the run log's step line nearest to time t."
function logvol(p, t)
    best = (Inf, NaN, NaN)
    for l in eachline(p)
        m = match(r"t=([0-9.]+)s.*vol_removed=([-0-9.e]+)m³", l)
        m === nothing && continue
        tt = parse(Float64, m[1]); abs(tt - t) < best[1] &&
            (best = (abs(tt - t), tt, parse(Float64, m[2])))
    end
    return best[2], best[3]
end

# the GeoParquet reader works through a temporary directory that Windows can
# briefly hold locked (antivirus, indexer): retry a failed load
function _load_mesh(p)
    for attempt in 1:5
        try
            return with_logger(() -> load_mesh_geoparquet(p), SimpleLogger(stderr, Logging.Error))
        catch err
            attempt == 5 && rethrow(); sleep(3)
        end
    end
end

rows = [split(l, ",") for l in readlines(joinpath(MESHDIR, "manifest.csv"))[2:end]]
out = open(joinpath(RUNDIR, "metrics.csv"), "w")
println(out, "case,mesh,res,cell_m,angle,offset_m,solver,deta,spill,unsteady,outflow")
for r in rows
    case, mname, res, cell, α, off = r[1], r[2], parse(Int, r[3]), parse(Float64, r[4]), parse(Float64, r[5]), parse(Float64, r[6])
    mesh = _load_mesh(String(r[10]))
    z = mesh.static_vars["elevation"]; area = [A5Grid._polygon_area_m2(c.boundary) for c in mesh.cells]
    sn = map(mesh.cells) do c
        x, y = xy(c.center_lon, c.center_lat)
        (x * cosd(α) + y * sind(α), -x * sind(α) + y * cosd(α) - off)
    end
    zone = [-0.3L <= s <= 0.05L for (s, _) in sn]
    chan = [zone[i] && abs(sn[i][2]) <= W / 2 + 0.5cell for i in eachindex(sn)]
    flood = [zone[i] && abs(sn[i][2]) > W / 2 + cell for i in eachindex(sn)]
    for solver in ("swe", "li", "limd")     # limd: li + --momentum-damping 0.05, if run
        p = joinpath(RUNDIR, "$(TAG)_$(case)_$(solver).h5")
        isfile(p) || (solver == "limd" || @warn "missing $p"; continue)
        nfr = h5open(f -> haskey(f, "frames") ? length(keys(f["frames"])) : 0, p)
        nfr >= 2 || (@warn "no output frames in $p (failed run?)"; continue)
        F0, F1 = lastframes(p)
        wetc = [i for i in findall(chan) if F1.h[i] > 0.01]
        deta = isempty(wetc) ? NaN : mean(F1.h[i] + z[i] - (Z0 - S * sn[i][1] + HN) for i in wetc)
        spill = sum(area[i] for i in findall(flood) if F1.h[i] > 0.05; init = 0.0) / (0.35L)
        unst = isempty(wetc) ? NaN : maximum(abs(F1.h[i] - F0.h[i]) for i in wetc)
        lp = replace(p, ".h5" => ".log")
        of = NaN
        if isfile(lp)
            (t0, v0), (t1, v1) = logvol(lp, F0.t), logvol(lp, F1.t)
            t1 > t0 && (of = (v1 - v0) / (Q * (t1 - t0)))
        end
        @printf(out, "%s,%s,%d,%.3f,%.1f,%.3f,%s,%.5f,%.3f,%.5f,%.4f\n", case, mname, res, cell, α, off, solver, deta, spill, unst, of)
        @printf("%-18s %-3s Δη %+7.3f m  spill %7.1f m  unsteady %.4f  outflow %.3f\n", case, solver, deta, spill, unst, of)
    end
end
close(out)
