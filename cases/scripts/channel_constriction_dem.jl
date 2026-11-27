# cases/scripts/channel_constriction_dem.jl — a straight channel with an abrupt
# contraction and expansion, on 5 m pixels (case `channel_constriction`).
#
#   julia --project=. cases/scripts/channel_constriction_dem.jl OUTDIR [ANGLE_DEG]
#
# Writes OUTDIR/dem5.tif (EPSG:4326, 5 m at 51°N, the origin on a whole 25 m
# block), OUTDIR/params.json and OUTDIR/inflow.csv.
#
# Geometry (s along the channel, ANGLE_DEG south of east, default 0: aligned
# with the grid; n across it, m): a rectangular channel W1 = 65 m wide
# (13 pixels) between vertical walls (+10 m), its bed falling at S = 0.0005;
# at s = SC … SC + LC = 2000 … 2025 m it narrows abruptly to W2 = 35 m
# (7 pixels; both sides step in by 15 m, a bridge opening). Manning n = 0.03;
# Q = 300 m³/s (normal depth 2.99 m in the wide channel, Froude 0.29; about
# 0.65 in the opening). The control: a broad-crested sill at s = XS … XS + LS =
# 3500 … 3550 m, its crest p above the bed set so that critical flow over it
# needs the normal flow's energy, p = h_n + v_n²/2g − 3/2·h_c (1.16 m): the
# channel above it runs at normal depth (no drawdown to choke the opening);
# below it the bed falls at 0.02 and past s = L = 3800 m to the raster's
# edges (open boundaries).
#
# The case's answer is 1-D (cases/scripts/channel_constriction_1d.jl): the
# level profile upstream of the opening for each set of equations, started
# from the run's own level downstream of it.

using ArchGDAL, JSON3

(1 <= length(ARGS) <= 2) || error("usage: julia --project=. cases/scripts/channel_constriction_dem.jl OUTDIR [ANGLE_DEG]")
const OUT = ARGS[1]
const θ = length(ARGS) == 2 ? parse(Float64, ARGS[2]) : 0.0
const P = 5.0
const W1, W2, SC, LC = 65.0, 35.0, 2000.0, 25.0
const S, S2, XS, LS, L, WALL, Z0 = 0.0005, 0.02, 3500.0, 50.0, 3800.0, 10.0, 40.0
const Q, NMAN = 300.0, 0.03
const HN = (Q * NMAN / (W1 * sqrt(S)))^(3 / 5)
const HC = ((Q / W1)^2 / 9.81)^(1 / 3)
const PSILL = HN + (Q / (W1 * HN))^2 / (2 * 9.81) - 1.5 * HC
const LAT0, LON0, R = 51.0, 0.0, 6_371_000.0
const S_IN = 20.0
const OX, OY = P / 2, P / 2                          # the axis through a pixel centre

d = (cosd(θ), -sind(θ)); nrm = (sind(θ), cosd(θ))
zaxis(s) = s < XS ? Z0 - S * s : s < XS + LS ? Z0 - S * XS + PSILL : Z0 - S * XS - S2 * (s - XS - LS)
width(s) = SC <= s < SC + LC ? W2 : W1
function ground(x, y)
    x -= OX; y -= OY
    s = x * d[1] + y * d[2]; n = x * nrm[1] + y * nrm[2]
    s > L && return zaxis(s)
    (s < 0 || abs(n) >= width(s) / 2) && return zaxis(max(s, 0.0)) + WALL
    return zaxis(s)
end
corners = [(OX + s * d[1] + n * nrm[1], OY + s * d[2] + n * nrm[2]) for s in (-30.0, L + 150.0), n in (-W1 / 2 - 20, W1 / 2 + 20)]
xmin = 25 * floor(minimum(first, corners) / 25); xmax = 25 * ceil(maximum(first, corners) / 25)
ymin = 25 * floor(minimum(last, corners) / 25); ymax = 25 * ceil(maximum(last, corners) / 25)
nx = round(Int, (xmax - xmin) / P); ny = round(Int, (ymax - ymin) / P)
z = [ground(xmin + (c + 0.5) * P, ymax - (r + 0.5) * P) for c in 0:(nx-1), r in 0:(ny-1)]
# the rasterised widths: channel pixels per metre in the wide reach and in the opening
function wr(s0, s1)
    c = count(Iterators.product(0:(nx-1), 0:(ny-1))) do (ci, ri)
        x = xmin + (ci + 0.5) * P - OX; y = ymax - (ri + 0.5) * P - OY
        s = x * d[1] + y * d[2]; n = x * nrm[1] + y * nrm[2]
        s0 <= s < s1 && abs(n) < width(s) / 2
    end
    return c * P^2 / (s1 - s0)
end
W1r = wr(500.0, 1900.0); W2r = wr(SC, SC + LC)
mlon = R * cosd(LAT0) * π / 180; mlat = R * π / 180
lon(x) = LON0 + x / mlon; lat(y) = LAT0 + y / mlat
mkpath(OUT)
ArchGDAL.create(joinpath(OUT, "dem5.tif"); driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                nbands = 1, dtype = Float64) do ds
    ArchGDAL.setgeotransform!(ds, [lon(xmin), P / mlon, 0.0, lat(ymax), 0.0, -P / mlat])
    ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
    ArchGDAL.write!(ds, z, 1)
end
xi, yi = OX + S_IN * d[1], OY + S_IN * d[2]
open(joinpath(OUT, "params.json"), "w") do io
    JSON3.write(io, Dict("angle_deg" => θ, "pixel_m" => P, "W1" => W1, "W2" => W2, "W1_raster" => W1r, "W2_raster" => W2r,
                         "SC" => SC, "LC" => LC, "S" => S, "S2" => S2, "XS" => XS, "LS" => LS, "p_sill" => PSILL, "h_normal" => HN, "L" => L, "Z0" => Z0,
                         "Q" => Q, "n" => NMAN, "lat0" => lat(OY), "lon0" => lon(OX),
                         "m_per_deg_lon" => mlon, "m_per_deg_lat" => mlat,
                         "inflow_lat" => lat(yi), "inflow_lon" => lon(xi), "nx" => nx, "ny" => ny))
end
open(joinpath(OUT, "inflow.csv"), "w") do io
    print(io, "t_s,Q_m3s\n0,0\n900,$Q\n1000000,$Q\n")
end
println("$(joinpath(OUT, "dem5.tif")): $nx x $ny pixels, angle $(θ)°; rasterised widths $(round(W1r; digits = 2)) / " *
        "$(round(W2r; digits = 2)) m; inflow at $(round(lat(yi); digits = 9)),$(round(lon(xi); digits = 9))")
