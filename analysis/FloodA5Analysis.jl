# analysis/FloodA5Analysis.jl
# ===========================
# One module for loading FloodA5 runs and computing metrics from them.
#
# Pure Julia: it does not include FloodModel.jl and needs no Python, so
# analysis scripts and notebooks load in seconds. Output files are read only
# through io/FloodA5Output.jl (FO). Format-2 files (0.4.0+) carry the mesh
# geometry; for format-1 files pass the mesh's GeoParquet (`parquet = ...`).
#
#   include("analysis/FloodA5Analysis.jl"); const FA = FloodA5Analysis
#   r   = FA.load_run("test/wp1_runs/run.h5")
#   k   = FA.frame_at(r, 3600.0)
#   h   = FA.depth(r, k)
#   ref = FA.load_run("ref.h5"); FA.compare(ref, r)
#
# Files:
#   runs.jl         Run and MeshGeom, loading, frame selection, fields,
#                   local coordinates, storage
#   parquet_mesh.jl format-1 fallback: mesh geometry from GeoParquet
#   compare.jl      run vs reference (rms, CSI, storage, speed), twin-run
#                   divergence
#   shape.jl        point-spread and planar-slope shape/direction metrics
#   noise.jl        checkerboard / velocity-noise metrics
#   channel.jl      straight-channel series metrics
#
# Conventions: cell vectors are in the output file's cell order; frame k is
# the 1-based frame index; times in seconds; distances in metres.

module FloodA5Analysis

using HDF5, Statistics, Printf, LinearAlgebra

if isdefined(Main, :FloodA5Output)
    const FO = Main.FloodA5Output
else
    include(joinpath(@__DIR__, "..", "io", "FloodA5Output.jl"))
    const FO = FloodA5Output
end

include("runs.jl")
include("parquet_mesh.jl")
include("compare.jl")
include("shape.jl")
include("noise.jl")
include("channel.jl")

end # module
