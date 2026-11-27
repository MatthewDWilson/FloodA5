# cases/scripts/block_dem.jl — exact k × k block means of a raster.
#
#   julia --project=. cases/scripts/block_dem.jl IN.tif OUT.tif K
#
# Pixel (c, r) of OUT is the mean of IN's pixels K·c … K·c+K−1 by
# K·r … K·r+K−1 (0-based), counted from IN's own origin: OUT's pixels are whole
# blocks of IN's, so a Cartesian mesh from OUT (`--meshgen-cartesian OUT.tif`)
# lines up with IN's pixels exactly and its SGS tables can be built from them
# (`--dem IN.tif --sgs-sampling pixels`). Pixels of IN past the last whole
# block are dropped. A block with nodata pixels takes the mean of the valid
# ones (nodata if none). Unlike `average_dem.jl` (GDAL `-tap`, which snaps to
# multiples of DX and so pads the edge blocks when IN's origin is not on
# one), nothing is resampled. Written as Float64, in IN's CRS.

using ArchGDAL

length(ARGS) == 3 || error("usage: julia --project=. cases/scripts/block_dem.jl IN.tif OUT.tif K")
src, dest = ARGS[1], ARGS[2]
K = parse(Int, ARGS[3]); K >= 1 || error("K must be a positive integer")
mkpath(dirname(abspath(dest)))

gt, proj, z, nd = ArchGDAL.read(src) do ds
    b = ArchGDAL.getband(ds, 1)
    raw = ArchGDAL.getnodatavalue(b)
    ArchGDAL.getgeotransform(ds), ArchGDAL.getproj(ds), Float64.(ArchGDAL.read(ds, 1)),
    raw === nothing ? NaN : Float64(raw)
end
(gt[3] == 0.0 && gt[5] == 0.0) || error("block_dem: rotated rasters are not supported")
nx, ny = size(z)                                   # ArchGDAL layout: (columns, rows)
mx, my = nx ÷ K, ny ÷ K
(mx >= 1 && my >= 1) || error("block_dem: the raster ($nx × $ny) is smaller than one $K × $K block")
valid(v) = isfinite(v) && !(isfinite(nd) && v == nd)
const ND = -9999.0
out = fill(ND, mx, my)
for r in 0:(my - 1), c in 0:(mx - 1)
    s = 0.0; n = 0
    for rr in (K * r + 1):(K * r + K), cc in (K * c + 1):(K * c + K)
        v = z[cc, rr]; valid(v) || continue
        s += v; n += 1
    end
    n > 0 && (out[c + 1, r + 1] = s / n)
end
ArchGDAL.create(dest; driver = ArchGDAL.getdriver("GTiff"), width = mx, height = my,
                nbands = 1, dtype = Float64) do ds
    ArchGDAL.setgeotransform!(ds, [gt[1], K * gt[2], 0.0, gt[4], 0.0, K * gt[6]])
    isempty(proj) || ArchGDAL.setproj!(ds, proj)
    ArchGDAL.setnodatavalue!(ArchGDAL.getband(ds, 1), ND)
    ArchGDAL.write!(ds, out, 1)
end
println(dest, ": ", mx, " x ", my, " blocks of ", K, " x ", K, " pixels (dropped ",
        nx - K * mx, " columns, ", ny - K * my, " rows)")
