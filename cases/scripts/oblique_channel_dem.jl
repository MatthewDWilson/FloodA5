# cases/scripts/oblique_channel_dem.jl — a straight channel at an angle to the
# grid, on 5 m pixels (case `oblique_channel`).
#
#   julia --project=. cases/scripts/oblique_channel_dem.jl OUTDIR ANGLE_DEG [P=5 W=15 D=2 Q=25 B=50]
#
# Writes OUTDIR/dem5.tif (EPSG:4326, pixels of 5 m at 51°N, the raster's
# origin on a whole 25 m block), OUTDIR/params.json and OUTDIR/inflow.csv.
#
# Geometry (s along the channel from its upstream end, pointing ANGLE_DEG
# south of east; n across it, m): a rectangular channel W = 15 m wide (three
# pixels) and D = 2 m deep, its bed falling at S = 0.002 to s = XB = 1300 m
# and at 0.02 beyond; floodplains at the bank top out to |n| = 50 m; walls
# (+10 m) beyond |n| = 50 m and upstream of s = 0. Past s = L = 1400 m the
# ground keeps falling at 0.02 in every direction to the raster's edges (open
# boundaries): the flow passes critical depth at the break, so no boundary
# formulation reaches the measuring reach s = 400–1000 m. Manning n = 0.035.
#
# Analytic answer (no simulation): uniform flow, no lateral exchange, the 2-D
# models' friction on the bed only: Q = W·h^(5/3)·S^(1/2)/n, so the normal
# depth is h_n = (Q·n / (W·√S))^(3/5) — 1.173 m for Q = 25 m³/s (in bank).
# Pixels are channel where their centre is within W/2 of the axis (the axis
# runs through a pixel centre); the channel pixels' area per metre of reach,
# W_r, depends on the angle (15.0 m at 0°, 15.71 m at 26.565°, 17.62 m at
# 45°), so `h_normal` is the normal depth for W_r (`h_normal_design` for W).

using ArchGDAL, JSON3

length(ARGS) >= 2 || error("usage: julia --project=. cases/scripts/oblique_channel_dem.jl OUTDIR ANGLE_DEG [P=5 W=15 D=2 Q=25 B=50]")
const OUT = ARGS[1]
const θ = parse(Float64, ARGS[2])
const OPT = Dict(String(first(split(a, "="))) => parse(Float64, last(split(a, "="))) for a in ARGS[3:end])
const P = get(OPT, "P", 5.0)                        # pixel (m); 25/P must be a whole number
const W, D, B = get(OPT, "W", 15.0), get(OPT, "D", 2.0), get(OPT, "B", 50.0)
const S, S2, XB, L, WALL = 0.002, 0.02, 1300.0, 1400.0, 10.0
const Z0 = 30.0
const Q, NMAN = get(OPT, "Q", 25.0), 0.035
const LAT0, LON0, R = 51.0, 0.0, 6_371_000.0
const S_IN = 20.0                                   # inflow point on the axis (m from the upstream end)
const OX, OY = P / 2, P / 2                         # the axis origin: a pixel centre (centres lie at P/2 + P·k)

d = (cosd(θ), -sind(θ)); nrm = (sind(θ), cosd(θ))   # along (east, north), across (to the left)
zaxis(s) = s <= XB ? Z0 - S * s : Z0 - S * XB - S2 * (s - XB)
function ground(x, y)                               # x east, y north (m); the axis starts at (OX, OY)
    x -= OX; y -= OY
    s = x * d[1] + y * d[2]; n = x * nrm[1] + y * nrm[2]
    s > L && return zaxis(s)                        # the outlet: falling ground to the edges
    (s < 0 || abs(n) > B) && return zaxis(max(s, 0.0)) + D + WALL
    abs(n) < W / 2 && return zaxis(s)
    return zaxis(s) + D
end
# raster box: the strip from s = −30 to L + 150 m, |n| ≤ B + 20 m, snapped to whole 25 m blocks
corners = [(OX + s * d[1] + n * nrm[1], OY + s * d[2] + n * nrm[2]) for s in (-30.0, L + 150.0), n in (-(B + 20.0), B + 20.0)]
xmin = 25 * floor(minimum(first, corners) / 25); xmax = 25 * ceil(maximum(first, corners) / 25)
ymin = 25 * floor(minimum(last, corners) / 25); ymax = 25 * ceil(maximum(last, corners) / 25)
nx = round(Int, (xmax - xmin) / P); ny = round(Int, (ymax - ymin) / P)
z = [ground(xmin + (c + 0.5) * P, ymax - (r + 0.5) * P) for c in 0:(nx-1), r in 0:(ny-1)]
mlon = R * cosd(LAT0) * π / 180; mlat = R * π / 180
lon(x) = LON0 + x / mlon; lat(y) = LAT0 + y / mlat
mkpath(OUT)
ArchGDAL.create(joinpath(OUT, "dem5.tif"); driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                nbands = 1, dtype = Float64) do ds
    ArchGDAL.setgeotransform!(ds, [lon(xmin), P / mlon, 0.0, lat(ymax), 0.0, -P / mlat])
    ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
    ArchGDAL.write!(ds, z, 1)
end
hn = (Q * NMAN / (W * sqrt(S)))^(3 / 5)
xi, yi = OX + S_IN * d[1], OY + S_IN * d[2]
# the rasterised channel's width: channel pixels per metre of the reach s = 400–1000 m
nch = count(c -> begin
        x = xmin + (c[1] + 0.5) * P - OX; y = ymax - (c[2] + 0.5) * P - OY
        sv = x * d[1] + y * d[2]; nv = x * nrm[1] + y * nrm[2]
        400 <= sv < 1000 && abs(nv) < W / 2
    end, Iterators.product(0:(nx-1), 0:(ny-1)))
W_raster = nch * P^2 / 600
open(joinpath(OUT, "params.json"), "w") do io
    JSON3.write(io, Dict("angle_deg" => θ, "pixel_m" => P, "W" => W, "D" => D, "S" => S, "S2" => S2, "XB" => XB,
                         "L" => L, "B" => B, "Z0" => Z0, "Q" => Q, "n" => NMAN, "h_normal_design" => hn, "h_normal" => (Q * NMAN / (W_raster * sqrt(S)))^(3 / 5),
                         "lat0" => lat(OY), "lon0" => lon(OX), "W_raster" => W_raster, "m_per_deg_lon" => mlon, "m_per_deg_lat" => mlat,
                         "inflow_lat" => lat(yi), "inflow_lon" => lon(xi), "nx" => nx, "ny" => ny,
                         "xmin" => xmin, "ymax" => ymax))
end
open(joinpath(OUT, "inflow.csv"), "w") do io
    print(io, "t_s,Q_m3s\n0,0\n600,$Q\n1000000,$Q\n")
end
println("$(joinpath(OUT, "dem5.tif")): $nx x $ny pixels, angle $(θ)°; normal depth $(round(hn; digits = 3)) m; rasterised width $(round(W_raster; digits = 2)) m (normal depth $(round((Q * NMAN / (W_raster * sqrt(S)))^(3 / 5); digits = 3)) m); " *
        "inflow at $(round(lat(yi); digits = 9)),$(round(lon(xi); digits = 9))")
