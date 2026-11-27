# cases/scripts/sgs_pool_dem.jl — the terrain of case `sgs_pool`: a closed,
# rough basin on a 5 m grid in which every pixel drains to one sink.
#
#   julia --project=. cases/scripts/sgs_pool_dem.jl OUTDIR
#
# Writes OUTDIR/sgs_pool_dem5.tif (100 × 80 pixels of 5 m, British National
# Grid, north-west corner at 338495 E, 557760 N — not on a multiple of 25 m,
# like the Carlisle DEM), OUTDIR/sgs_pool_params.json and the inflow
# hydrograph OUTDIR/inflow.csv.
#
# Terrain (x east, y north from the south-west corner, m): a bowl 2·(u² + v²)
# (u, v from −1 to 1 across the 500 × 400 m domain) on a 10 m datum, a tilt
# of 0.2% eastwards, ridges and hollows on 25–80 m wavelengths (±0.35 m),
# pixel noise (±0.15 m), and a sinuous channel one pixel wide and 1 m deep
# (y = 200 + 60·sin(2πx/250)) — the scales of ground a 25 m cell has to hold
# inside it. A 0.3 m hole at the inflow pixel (x = 252.5, y = 202.5) makes it
# the lowest point. Then every closed depression is filled to its spill level
# towards that sink (4-connected, as the models connect cells): each pixel
# gets the lowest, over the paths from it to the sink, of the highest ground
# on the path. Water at rest then stands at one level η over exactly the
# pixels below η, so the volume of a pool at η is Σ_p a_p·max(η − z_p, 0)
# with nothing disconnected — the exact answer the case compares with.
#
# The design volume V0 fills the pixels below the 45th percentile of the
# ground; the params file records it, the level η0 it gives (pixels of
# 25 m²), and the inflow: Q = V0 / (T_fill − 60 s), ramped over 60 s at
# both ends of T_fill = 2 h.

using ArchGDAL, Random, JSON3

length(ARGS) == 1 || error("usage: julia --project=. cases/scripts/sgs_pool_dem.jl OUTDIR")
const OUT = ARGS[1]
const P = 5.0
const NX, NY = 100, 80
const X0, Y0 = 338495.0, 557760.0             # north-west corner (EPSG:27700)
const SINK = (50, 39)                         # 0-based (column, row from the north): x = 252.5, y = 202.5
const T_FILL = 7200.0
const INFLOW_LATLON = (54.909001, -2.956895638)   # the sink pixel's centre in WGS84 (the case's argv)

xc(c) = (c + 0.5) * P                         # m from the west edge
yc(r) = (NY - r - 0.5) * P                     # m from the south edge
function ground(c, r, rng_phase, noise)
    x, y = xc(c), yc(r)
    u = (x - 250.0) / 250.0; v = (y - 200.0) / 200.0
    z = 10.0 + 2.0 * (u^2 + v^2) + 0.002 * x
    for (k, (lx, ly, a)) in enumerate(((80.0, 55.0, 0.15), (45.0, -70.0, 0.10), (25.0, 30.0, 0.10)))
        z += a * sin(2π * (x / lx + y / ly) + rng_phase[k])
    end
    z += noise[c+1, r+1]
    abs(y - (200.0 + 60.0 * sin(2π * x / 250.0))) < P / 2 && (z -= 1.0)
    return z
end

rng = MersenneTwister(20261127)
phase = 2π .* rand(rng, 3)
noise = 0.3 .* rand(rng, NX, NY) .- 0.15
z = [ground(c, r, phase, noise) for c in 0:(NX-1), r in 0:(NY-1)]        # (nx × ny), ArchGDAL's layout
z[SINK[1]+1, SINK[2]+1] = minimum(z) - 0.3

# fill every depression to its spill level towards the sink: f(p) = min over
# paths p → sink of the highest ground on the path (relaxation to the fixed
# point f(p) = max(z(p), min over 4-neighbours f(q)))
f = fill(Inf, NX, NY); f[SINK[1]+1, SINK[2]+1] = z[SINK[1]+1, SINK[2]+1]
changed = true; sweeps = 0
while changed
    global changed = false; global sweeps += 1
    for pass in 1:2
        cs = pass == 1 ? (1:NX) : (NX:-1:1); rs = pass == 1 ? (1:NY) : (NY:-1:1)
        for j in rs, i in cs
            (i, j) == (SINK[1] + 1, SINK[2] + 1) && continue
            m = Inf
            i > 1 && (m = min(m, f[i-1, j])); i < NX && (m = min(m, f[i+1, j]))
            j > 1 && (m = min(m, f[i, j-1])); j < NY && (m = min(m, f[i, j+1]))
            v = max(z[i, j], m)
            if v < f[i, j]
                f[i, j] = v; changed = true
            end
        end
    end
end
n_raised = count(f .> z .+ 1e-12)
z = f

zs = sort(vec(z))
η0 = zs[round(Int, 0.45 * length(zs))]
V0 = sum(P^2 * max(η0 - v, 0.0) for v in z)
Q = V0 / (T_FILL - 60.0)

mkpath(OUT)
dem = joinpath(OUT, "sgs_pool_dem5.tif")
ArchGDAL.create(dem; driver = ArchGDAL.getdriver("GTiff"), width = NX, height = NY,
                nbands = 1, dtype = Float64) do ds
    ArchGDAL.setgeotransform!(ds, [X0, P, 0.0, Y0, 0.0, -P])
    ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(27700)))
    ArchGDAL.write!(ds, z, 1)
end
# the inflow point (WGS84) back in the raster: it must fall in the sink pixel
px = [INFLOW_LATLON[2]]; py = [INFLOW_LATLON[1]]; pz = [0.0]
function crs_gis(epsg)                       # x = easting/longitude (as A5Grid._crs_gis)
    c = ArchGDAL.importEPSG(epsg)
    ArchGDAL.GDAL.osrsetaxismappingstrategy(c, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    return c
end
src = crs_gis(4326); dst = crs_gis(27700)
ArchGDAL.createcoordtrans(src, dst) do t
    ArchGDAL.transform!(px, py, pz, t)
end
inflow_pixel = (floor(Int, (px[1] - X0) / P), floor(Int, (Y0 - py[1]) / P))
open(joinpath(OUT, "sgs_pool_params.json"), "w") do io
    JSON3.write(io, Dict(
        "pixel_m" => P, "nx" => NX, "ny" => NY, "x0" => X0, "y0" => Y0, "epsg" => 27700,
        "sink_col" => SINK[1], "sink_row" => SINK[2],
        "inflow_lat" => INFLOW_LATLON[1], "inflow_lon" => INFLOW_LATLON[2],
        "inflow_pixel_ok" => inflow_pixel == SINK,
        "inflow_offset_m" => hypot(px[1] - (X0 + (SINK[1] + 0.5) * P), py[1] - (Y0 - (SINK[2] + 0.5) * P)),
        "fill_sweeps" => sweeps, "n_pixels_raised" => n_raised,
        "z_min" => zs[1], "z_max" => zs[end], "eta0" => η0, "V0_m3" => V0,
        "Q_m3s" => Q, "T_fill_s" => T_FILL))
end
open(joinpath(OUT, "inflow.csv"), "w") do io
    print(io, "t_s,Q_m3s
0,0
60,$Q
$(T_FILL - 60),$Q
$T_FILL,0
1000000,0
")
end
println("$dem: $NX x $NY pixels; η0 = $(round(η0; digits = 3)) m, V0 = $(round(V0; digits = 1)) m³, " *
        "Q = $(round(Q; digits = 4)) m³/s; $n_raised pixels raised by the fill ($sweeps sweeps); " *
        "inflow in the sink pixel: $(inflow_pixel == SINK)")
