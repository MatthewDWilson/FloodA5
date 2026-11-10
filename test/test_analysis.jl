# test/test_analysis.jl — analysis/FloodA5Analysis.jl (FA)
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_analysis.jl
#
# Builds one small format-2 fixture run (Cartesian, 25 m, a point inflow on a
# flat bed) and its GeoParquet mesh, then runs test/analysis/test_*.jl.
# Fixtures for the test files: FIX_H5 (format 2), FIX_MESH (its parquet),
# FIX_DIR (a temp dir), FIX_SRC = (lon, lat) of the inflow point.
# FA_TEST_ONLY=shape,compare runs only test/analysis/test_shape.jl, test_compare.jl.

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis
using Test, HDF5, Logging

const FIX_DIR = mktempdir()
const FIX_SRC = (0.0012, 0.0014)        # inside a cell, not on a corner
const FIX_MESH = joinpath(FIX_DIR, "m.parquet")
const FIX_H5 = joinpath(FIX_DIR, "fix.h5")
let aoi = joinpath(FIX_DIR, "aoi.geojson"), hyd = joinpath(FIX_DIR, "hydro.csv")
    write(aoi, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
      "geometry":{"type":"Polygon","coordinates":[[[0.0,0.0],[0.0027,0.0],[0.0027,0.0027],[0.0,0.0027],[0.0,0.0]]]}}]}""")
    write(hyd, "t_s,Q\n0,2\n600,2\n")
    args = ["--meshgen-cartesian", aoi, "--dx", "25", "--meshout", FIX_MESH,
            "--flow-model", "standard", "--closed-boundaries",
            "--inflow-point", "$(FIX_SRC[2]),$(FIX_SRC[1]),$hyd",   # LAT,LON,FILE
            "--sim-duration", "600", "--output", FIX_H5, "--output-interval", "60"]
    with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            main(copy(args))
        end
    end
end

const _ONLY = filter(!isempty, split(get(ENV, "FA_TEST_ONLY", ""), ","))
for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl") &&
                          (isempty(_ONLY) || f[6:end-3] in _ONLY),
                     readdir(joinpath(@__DIR__, "analysis"))))
    include(joinpath(@__DIR__, "analysis", f))
end
