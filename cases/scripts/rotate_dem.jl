# cases/scripts/rotate_dem.jl — a DEM resampled onto a grid rotated against its own
# (a reference run whose grid meets the channels at other angles).
#
#   julia --project=. cases/scripts/rotate_dem.jl IN.tif OUT.tif ANGLE_DEG [PIXEL_M]
#
# The target CRS is an oblique Mercator through the raster's centre whose
# grid axes lie along and across a centre line at azimuth ANGLE_DEG (Hotine,
# alpha = ANGLE_DEG, gamma = 0: no rectification; IN's ellipsoid), so OUT is north-up in its own
# CRS and FloodModel's raster-aligned Cartesian mesh (`--meshgen-cartesian
# OUT.tif`) runs on it. Pixels PIXEL_M (default: IN's), bilinear resampling,
# -9999 outside IN. Bilinear resampling smooths features one or two pixels
# wide (walls, banks): compare runs on OUT with each other, not with runs on IN.
# A mesh from OUT keeps its nodata cells (inert): drop them with
# cases/scripts/drop_nodata_cells.jl so that the domain's edge is IN's edge.

using ArchGDAL

(3 <= length(ARGS) <= 4) || error("usage: julia --project=. cases/scripts/rotate_dem.jl IN.tif OUT.tif ANGLE_DEG [PIXEL_M]")
src, dest = ARGS[1], ARGS[2]
θ = parse(Float64, ARGS[3])
ArchGDAL.read(src) do ds
    gt = ArchGDAL.getgeotransform(ds); nx, ny = ArchGDAL.width(ds), ArchGDAL.height(ds)
    px = length(ARGS) == 4 ? ARGS[4] : string(abs(gt[2]))
    # the raster's centre in lon/lat
    sref = ArchGDAL.importWKT(ArchGDAL.getproj(ds))
    ArchGDAL.GDAL.osrsetaxismappingstrategy(sref, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    geo = ArchGDAL.importEPSG(4326)
    ArchGDAL.GDAL.osrsetaxismappingstrategy(geo, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    x = [gt[1] + gt[2] * nx / 2]; y = [gt[4] + gt[6] * ny / 2]; z = [0.0]
    ArchGDAL.createcoordtrans(t -> ArchGDAL.transform!(x, y, z, t), sref, geo)
    # the target: IN's own geographic CRS (datum and all) with an oblique Mercator
    tgt = ArchGDAL.importEPSG(4326)
    ArchGDAL.GDAL.osrcopygeogcsfrom(tgt.ptr, sref.ptr)
    ArchGDAL.GDAL.osrsethom(tgt.ptr, y[1], x[1], θ, 0.0, 1.0, 0.0, 0.0)
    ArchGDAL.GDAL.osrsetaxismappingstrategy(tgt.ptr, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    omerc = ArchGDAL.toWKT(tgt)
    # check: IN's centre to lon/lat directly and through the target CRS
    xa = [gt[1] + gt[2] * nx / 2]; ya = [gt[4] + gt[6] * ny / 2]; za = [0.0]
    ArchGDAL.createcoordtrans(t -> ArchGDAL.transform!(xa, ya, za, t), sref, tgt)
    ArchGDAL.createcoordtrans(t -> ArchGDAL.transform!(xa, ya, za, t), tgt, geo)
    println("target CRS: Hotine oblique Mercator, azimuth $(θ)°, on IN's datum; round trip of the centre ",
            round(hypot((xa[1] - x[1]) * 111320 * cosd(y[1]), (ya[1] - y[1]) * 111320); sigdigits = 2), " m")
    ArchGDAL.gdalwarp([ds], ["-t_srs", omerc, "-tr", px, px, "-r", "bilinear", "-dstnodata", "-9999",
                             "-ot", "Float64"]; dest = dest) do out
        println(dest, ": ", ArchGDAL.width(out), " x ", ArchGDAL.height(out), " pixels of ", px, " m, turned ", θ, "°")
    end
end
