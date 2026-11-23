# test/test_dem_registration.jl — DEM sampling at pixel centres (0.8.0)
#
# A GeoTIFF pixel's value is its value at the pixel's centre (GDAL's
# pixel-is-area convention; the geotransform is corner-based). Until 0.7.1
# `A5Grid._bilinear` used the corner-based pixel coordinates of
# `_affine_to_pixel` as centre-based ones, so every sample was taken half a
# pixel east and south of its point.
#
#   R1  a plane written at pixel centres is sampled exactly (bilinear is exact
#       for a plane) at arbitrary interior points
#   R2  within half a pixel of the raster's edge the value is the nearest
#       centre line's (constant across the edge); outside the raster NaN
#   R3  nodata: a point in a nodata pixel is NaN; a point next to one takes
#       the value of the pixel containing it
#   R4  sample_dem_centroid! / sample_dem_mean! on a small mesh: the plane at
#       the cell centres (the mean to the Halton sampling error), and the
#       mesh records `dem_registration`
#
# Run:  julia --project=. test/test_dem_registration.jl
include(joinpath(@__DIR__, "..", "mesh", "A5Grid.jl"))
using .A5Grid
using Test, ArchGDAL, Logging

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

const DIR = mktempdir()
const PX = 1e-4                         # pixel size (degrees)
const LON0, LAT0 = -0.003, 0.002        # north-west corner of the raster
const NX, NY = 60, 40
plane(lon, lat) = 3.0 + 2000.0 * lon + 5000.0 * lat
centre(c, r) = (LON0 + (c + 0.5) * PX, LAT0 - (r + 0.5) * PX)   # 0-based pixel (c, r)

function write_dem(path, z; nodata = nothing)
    ArchGDAL.create(path; driver = ArchGDAL.getdriver("GTiff"), width = NX, height = NY,
                    nbands = 1, dtype = Float64) do ds
        ArchGDAL.setgeotransform!(ds, [LON0, PX, 0.0, LAT0, 0.0, -PX])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
        nodata === nothing || ArchGDAL.setnodatavalue!(ArchGDAL.getband(ds, 1), nodata)
        ArchGDAL.write!(ds, z, 1)
    end
    return path
end

const ZP = [plane(centre(c, r)...) for c in 0:(NX-1), r in 0:(NY-1)]   # (nx × ny), ArchGDAL's layout
const DEM = write_dem(joinpath(DIR, "plane.tif"), ZP)
const R = A5Grid._load_dem_raster(DEM)
sample(lon, lat) = A5Grid._bilinear(R.band, A5Grid._affine_to_pixel(lon, lat, R.gt)..., R.nx, R.ny, R.nodata_val)

@testset "R1 plane sampled exactly at interior points" begin
    err = 0.0
    for k in 1:500
        lon = LON0 + (0.5 + (NX - 1) * ((k * 0.6180339887) % 1)) * PX
        lat = LAT0 - (0.5 + (NY - 1) * ((k * 0.7548776662) % 1)) * PX
        err = max(err, abs(sample(lon, lat) - plane(lon, lat)))
    end
    @test err < 1e-9
    # the pixel's own value at its centre
    @test sample(centre(7, 11)...) ≈ ZP[8, 12] atol = 1e-12
    # (until 0.7.1 the error was half a pixel times the gradient: 0.35 m here)
    @test abs(sample(centre(7, 11)...) - plane(centre(7, 11)...)) < 1e-12
end

@testset "R2 edges" begin
    # a quarter pixel inside the west edge: the first centre column's value
    lat = centre(0, 20)[2]
    @test sample(LON0 + 0.25PX, lat) ≈ plane(centre(0, 20)[1], lat) atol = 1e-9
    # a quarter pixel inside the south-east corner
    @test sample(LON0 + (NX - 0.25) * PX, LAT0 - (NY - 0.25) * PX) ≈ ZP[NX, NY] atol = 1e-9
    @test isnan(sample(LON0 - 0.1PX, lat))
    @test isnan(sample(centre(5, 0)[1], LAT0 + 0.1PX))
end

@testset "R3 nodata" begin
    z = copy(ZP); z[11, 11] = -9999.0                       # pixel (10, 10)
    r = A5Grid._load_dem_raster(write_dem(joinpath(DIR, "nodata.tif"), z; nodata = -9999.0))
    s(lon, lat) = A5Grid._bilinear(r.band, A5Grid._affine_to_pixel(lon, lat, r.gt)..., r.nx, r.ny, r.nodata_val)
    lon, lat = centre(10, 10)
    @test isnan(s(lon, lat))
    @test isnan(s(lon + 0.4PX, lat - 0.4PX))
    # in pixel (11, 10), a quarter pixel from the nodata one: that pixel's value
    @test s(centre(11, 10)[1] - 0.25PX, lat) == ZP[12, 11]
    # far from it: the plane
    @test s(centre(30, 30)...) ≈ plane(centre(30, 30)...) atol = 1e-9
end

@testset "R4 mesh sampling and provenance" begin
    "A square cell of side 5 pixels centred on (lon, lat)."
    function cell(id, lon, lat)
        h = 2.5PX
        A5Cell(id, 0, lon, lat, [[lon - h, lat - h], [lon + h, lat - h], [lon + h, lat + h],
                                 [lon - h, lat + h], [lon - h, lat - h]])
    end
    pts = [(-0.0021, 0.0012), (-0.00137, -0.00093), (0.0011, 0.0004)]
    mk() = A5Mesh(0, "", [cell("c$k", p...) for (k, p) in enumerate(pts)])
    m = mk()
    _quiet(() -> sample_dem_centroid!(m, DEM))
    @test all(abs(m.static_vars["elevation"][k] - plane(pts[k]...)) < 1e-9 for k in eachindex(pts))
    @test m.static_vars["dem_registration"] == fill(1.0, length(pts))
    @test m.cells[2].elevation == m.static_vars["elevation"][2]
    m2 = mk()
    _quiet(() -> sample_dem_mean!(m2, DEM))
    # the plane's mean over the cell is its centre value; Halton points' mean
    # is off the centre by ≲ 1% of the cell (0.035 m here; the old offset: 0.35 m)
    @test all(abs(m2.static_vars["elevation"][k] - plane(pts[k]...)) < 0.05 for k in eachindex(pts))
    @test haskey(m2.static_vars, "dem_registration")
end
