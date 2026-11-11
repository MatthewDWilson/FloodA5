# cases/scripts/average_dem.jl — average a DEM to a coarser raster-aligned grid.
#
#   julia --project=. cases/scripts/average_dem.jl IN.tif OUT.tif DX
#
# GDAL warp with `-tr DX DX -r average -tap` (cell means of the input pixels,
# grid snapped to multiples of DX in the input CRS). This is how the Carlisle
# 10/20/40 m Cartesian DEMs of HYDRAULICS.md §15.10–15.11 were made; a
# Cartesian mesh then comes from `FloodModel.jl --meshgen-cartesian OUT.tif`.

using ArchGDAL

length(ARGS) == 3 || error("usage: julia --project=. cases/scripts/average_dem.jl IN.tif OUT.tif DX")
src, dest, dx = ARGS[1], ARGS[2], ARGS[3]
parse(Float64, dx) > 0 || error("DX must be positive")
mkpath(dirname(abspath(dest)))
ArchGDAL.read(src) do ds
    ArchGDAL.gdalwarp([ds], ["-tr", dx, dx, "-r", "average", "-tap"]; dest = dest) do out
        println(dest, ": ", ArchGDAL.width(out), " x ", ArchGDAL.height(out), " cells of ", dx, " m")
    end
end
