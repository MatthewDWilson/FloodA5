# test/channel_bench/make_channel_meshes.jl
# -----------------------------------------
# Channel-representation series (2026-10): a straight rectangular channel,
# narrower or wider than a cell, cut into a sloping floodplain at a range of
# orientations, on A5 and on Cartesian grids of the same mean cell area, over
# a range of resolutions. Writes one mesh per (mesh type, resolution, angle,
# offset) with the bed sampled from the analytic surface by the same Halton
# cell mean as `sample_dem_mean!`, the inflow file, a manifest (CSV) and a
# job list for test/run_jobs.sh-style runners.
#
#   FLOODA5_NO_VIS=1 julia --project=. --threads 1 test/channel_bench/make_channel_meshes.jl OUTDIR [--res 17,18,19]
#
# Geometry (local metric frame at the domain centre; s along the channel, n
# across it, both in metres):
#   bed z = Z0 − S·s − S2·max(s − SB, 0) + (|n'| > W/2 ? D + SX·(|n'| − W/2) : 0),
#   n' = n − off: a valley (banks D above the channel bed, floodplain rising
#   across the valley at SX, so water that leaves the channel drains back
#   towards it) that steepens to S + S2 below s = SB, a hydraulic control that
#   isolates the reach from the outflow boundary
#   domain: square, side L; every edge open (outflow only)
#   inflow Q at (s, n) = (−0.45 L, off), constant
# With n = 0.03 the 2-D normal depth in the channel (wide-channel Manning,
# R = h) is h_n = (Q n / (W √S))^(3/5) ≈ 0.97 m: in bank (D = 2 m).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
using Logging, Printf, Statistics

const OUT = ARGS[1]
_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const RES    = parse.(Int, split(_arg("--res", "17,18,19,20"), ","))
const ANGLES = parse.(Float64, split(_arg("--angles", "0,15,30,45,60,75,90"), ","))
const OFFS   = parse.(Float64, split(_arg("--offsets", "0,0.5"), ","))   # fraction of the cell size
const L, W, D, S, SX, Z0, Q, NMAN = 2000.0, 20.0, 2.0, 0.001, 0.005, 10.0, 20.0, 0.03
const S2, SB = 0.02, 0.3L
const LON0, LAT0 = 0.0, 51.0
const R = A5Grid._EARTH_R
_q(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))

dlon(x) = rad2deg(x / (R * cosd(LAT0))); dlat(y) = rad2deg(y / R)
xy(lon, lat) = (deg2rad(lon - LON0) * R * cosd(LAT0), deg2rad(lat - LAT0) * R)
bed(x, y, α, off) = (s = x * cosd(α) + y * sind(α); a = abs(-x * sind(α) + y * cosd(α) - off);
                     Z0 - S * s - S2 * max(s - SB, 0.0) + (a > W / 2 ? D + SX * (a - W / 2) : 0.0))

"Halton cell mean of `bed`, as sample_dem_mean! does for a raster."
function sample_bed!(mesh, α, off; n = 256)
    h2 = A5Grid._halton(2, n, 0); h3 = A5Grid._halton(3, n, 0)
    z = Vector{Float64}(undef, length(mesh.cells))
    for (ci, c) in enumerate(mesh.cells)
        b = c.boundary
        lo0 = minimum(Float64(v[1]) for v in b); lo1 = maximum(Float64(v[1]) for v in b)
        la0 = minimum(Float64(v[2]) for v in b); la1 = maximum(Float64(v[2]) for v in b)
        acc = 0.0; k = 0
        for j in 1:n
            lo = lo0 + h2[j] * (lo1 - lo0); la = la0 + h3[j] * (la1 - la0)
            A5Grid._point_in_ring(lo, la, b) || continue
            acc += bed(xy(lo, la)..., α, off); k += 1
        end
        z[ci] = k > 0 ? acc / k : bed(xy(c.center_lon, c.center_lat)..., α, off)
    end
    mesh.static_vars["elevation"] = z
    return mesh
end

mkpath(OUT)
qfile = "$OUT/q20.csv"; write(qfile, "t_s,Q_m3s\n0,$Q\n1000000,$Q\n")
aoi = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon","coordinates":[[[$(dlon(-L/2)),$(LAT0 + dlat(-L/2))],[$(dlon(L/2)),$(LAT0 + dlat(-L/2))],[$(dlon(L/2)),$(LAT0 + dlat(L/2))],[$(dlon(-L/2)),$(LAT0 + dlat(L/2))],[$(dlon(-L/2)),$(LAT0 + dlat(-L/2))]]]}}"""
manifest = open("$OUT/manifest.csv", "w")
println(manifest, "case,mesh,res,cell_m,angle,offset_m,n_cells,inflow_lat,inflow_lon,mesh_file")
jobs = open("$OUT/jobs.txt", "w")
const LI = "--flow-model standard"                      # the default scheme (0.3.0)
const RUN = "--sim-duration 10800 --dt-max 10 --output-interval 900 --manning-n $NMAN"
for res in RES
    a5 = _q(() -> A5Grid.mesh_for_aoi(aoi, res))
    cell = sqrt(mean(A5Grid._polygon_area_m2(c.boundary) for c in a5.cells))
    cart = _q(() -> generate_cartesian_mesh_bounds(dlon(-L/2), LAT0 + dlat(-L/2), dlon(L/2), LAT0 + dlat(L/2), cell; lat_ref = LAT0))
    for (mname, m) in (("a5", a5), ("cart", cart)), α in ANGLES, of in OFFS
        off = of * cell
        case = @sprintf("%s_r%d_a%02d_o%02d", mname, res, round(Int, α), round(Int, 100of))
        f = "$OUT/$case.parquet"
        if !isfile(f)                        # resumable: existing meshes are kept
            sample_bed!(m, α, off)
            for attempt in 1:5               # the save stages through a fixed temp file,
                try                          # which Windows can briefly hold locked
                    _q(() -> save_mesh_geoparquet(m, f)); break
                catch err
                    attempt == 5 && rethrow(); @warn "save of $case failed ($err); retrying"; sleep(5)
                end
            end
        end
        s0 = -0.45L; x = s0 * cosd(α) - off * sind(α); y = s0 * sind(α) + off * cosd(α)
        lat = LAT0 + dlat(y); lon = LON0 + dlon(x)
        @printf(manifest, "%s,%s,%d,%.3f,%.1f,%.3f,%d,%.9f,%.9f,%s\n", case, mname, res, cell, α, off, length(m.cells), lat, lon, f)
        src = "--meshload $f --inflow-point $(round(lat, digits = 9)),$(round(lon, digits = 9)),$qfile $RUN"
        println(jobs, "chan_$(case)_swe|$src --flow-model swe")
        println(jobs, "chan_$(case)_li|$src $LI")
        @printf("%s: %d cells, cell %.1f m\n", case, length(m.cells), cell)
    end
end
close(manifest); close(jobs)
