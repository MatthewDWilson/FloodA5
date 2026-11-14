# cases/scripts/sgs_strip_dem.jl — the SGS strip's DEM and inputs (case
# `sgs_strip`, cases/registry/sgs_strip.jl; docs/TEST_CASES.md §2c).
#
#   julia --project=. cases/scripts/sgs_strip_dem.jl OUTDIR [--length L] [--width B]
#       [--channel-width W] [--channel-depth D] [--slope S] [--pixel P]
#       [--wall H] [--channel-offset Y] [--cross-slope SX] [--z0 Z0]
#       [--break-x XB] [--break-slope S2]
#
# A flume of width B along x (east), sloping at S, with a rectangular channel
# of width W and depth D cut along y = Y (default the centre line), in a local
# metric frame about (LON0, LAT0) = (0°, 51°):
#     x = rad(lon − LON0)·R·cos(LAT0),  y = rad(lat − LAT0)·R,  R = 6 371 000 m
# (the Cartesian mesh generator's own equirectangular frame with
# `--cartesian-lat-ref 51`, so mesh rows run exactly along the channel).
#
# Bed:  z_f(x) = Z0 − S·x − S2·max(x − XB, 0) (+ SX·(|y − Y| − W/2) outside the channel)
#       channel |y − Y| < W/2:  z_f − D
#       walls  |y| > B/2 (sides) and x < 0 (upstream): z_f + H
#
# The raster stores the bed at its pixel *corners* (A5Grid._bilinear treats a
# pixel's value as sitting at its integer pixel coordinate), on nodes at
# (k + ½)·P in x and y, so the DEM the model samples is exactly the analytic
# surface: planes, with each vertical wall replaced by a ramp of width P
# centred on it (channel walls at |y − Y| = W/2 ± P/2; side walls between
# |y| = B/2 + P/2 and B/2 + 3P/2, outside a Cartesian mesh of width B; the
# upstream wall between x = −3P/2 and −P/2). `sgs_strip_params.json` records
# the geometry; the case's analytic answers are computed from it.
#
# Outputs in OUTDIR: sgs_strip_dem.tif (EPSG:4326, Float32), sgs_strip_params.json,
# sgs_strip_aoi.geojson (the flume, x ∈ [0, L], |y| ≤ B/2, shrunk by 1 cm so the
# Cartesian generator's ceil() gives exactly L/dx × B/dx cells),
# sgs_strip_bc.geojson (`--bc-file`: a Closed line round the flume, then a
# ZeroGradient line across the downstream end x = L; features match boundary
# cells whose centres lie within 1.5 cell diameters of the line, later ones
# overriding, so the last ~1.5 diameters of the strip are open on all sides).

using ArchGDAL, Printf

const OUT = ARGS[1]
_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : parse(Float64, ARGS[i + 1]))
const L  = _arg("--length", 3000.0)
const B  = _arg("--width", 60.0)
const W  = _arg("--channel-width", 5.0)
const D  = _arg("--channel-depth", 1.0)
const S  = _arg("--slope", 0.001)
const P  = _arg("--pixel", 0.5)
const H  = _arg("--wall", 5.0)
const YC = _arg("--channel-offset", 0.0)
const SX = _arg("--cross-slope", 0.0)
const Z0 = _arg("--z0", 10.0)
const XB = _arg("--break-x", Inf)      # the bed steepens to S + S2 below XB: a supercritical
const S2 = _arg("--break-slope", 0.0)  # control that isolates the strip from its outflow boundary
const LON0, LAT0, R = 0.0, 51.0, 6_371_000.0
const MX, MY = 60.0, 90.0                   # DEM extent: x ∈ [−MX, L + MX], |y| ≤ MY

lon(x) = LON0 + rad2deg(x / (R * cosd(LAT0)))
lat(y) = LAT0 + rad2deg(y / R)

include(joinpath(@__DIR__, "sgs_strip_geometry.jl"))
const PAR = (; L, B, W, D, S, P, H, YC, SX, Z0, XB, S2)
bed(x, y) = _ss_bed(PAR, x, y)               # walls, channel, floodplain (one definition)

mkpath(OUT)
xs = collect((-MX + P / 2):P:(L + MX))                     # node x, west → east
ys = collect(reverse((-MY + P / 2):P:MY))                  # node y, north → south
nx, ny = length(xs), length(ys)
z = Float32[bed(x, y) for x in xs, y in ys]                # (nx × ny), ArchGDAL's layout
dlon = rad2deg(P / (R * cosd(LAT0))); dlat = rad2deg(P / R)
dem = joinpath(OUT, "sgs_strip_dem.tif")
ArchGDAL.create(dem; driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                nbands = 1, dtype = Float32) do ds
    # node (c, r) at integer pixel coordinates (c, r): ulx/uly = the first node
    ArchGDAL.setgeotransform!(ds, [lon(xs[1]), dlon, 0.0, lat(ys[1]), 0.0, -dlat])
    ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
    ArchGDAL.write!(ds, z, 1)
end

ring(x0, y0, x1, y1) = "[[$(lon(x0)),$(lat(y0))],[$(lon(x1)),$(lat(y0))],[$(lon(x1)),$(lat(y1))],[$(lon(x0)),$(lat(y1))],[$(lon(x0)),$(lat(y0))]]"
ε = 0.01
write(joinpath(OUT, "sgs_strip_aoi.geojson"),
      """{"type":"Feature","properties":{},"geometry":{"type":"Polygon","coordinates":[$(ring(ε, -B/2 + ε, L - ε, B/2 - ε))]}}""")
write(joinpath(OUT, "sgs_strip_bc.geojson"), """
{"type":"FeatureCollection","features":[
 {"type":"Feature","properties":{"bc_type":"Closed","label":"flume walls"},
  "geometry":{"type":"LineString","coordinates":$(ring(0.0, -B/2, L, B/2))}},
 {"type":"Feature","properties":{"bc_type":"ZeroGradient","label":"downstream end"},
  "geometry":{"type":"LineString","coordinates":[[$(lon(L)),$(lat(-MY))],[$(lon(L)),$(lat(MY))]]}}]}
""")
write(joinpath(OUT, "sgs_strip_params.json"), """
{"length_m": $L, "width_m": $B, "channel_width_m": $W, "channel_depth_m": $D, "slope": $S,
 "pixel_m": $P, "wall_m": $H, "channel_offset_m": $YC, "cross_slope": $SX, "z0_m": $Z0, "break_x_m": $(isfinite(XB) ? XB : 1e30), "break_slope": $S2,
 "lon0": $LON0, "lat0": $LAT0, "earth_radius_m": $R, "dem_nx": $nx, "dem_ny": $ny,
 "dem_x0_m": $(xs[1]), "dem_y0_m": $(ys[1])}
""")
@printf("%s: %d × %d nodes, %.2f m; flume %.0f × %.0f m, channel %.1f × %.1f m at y = %.1f, S = %g\n",
        dem, nx, ny, P, L, B, W, D, YC, S)
