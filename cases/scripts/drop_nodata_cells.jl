# cases/scripts/drop_nodata_cells.jl — remove the cells without an elevation from a mesh.
#
#   julia --project=. cases/scripts/drop_nodata_cells.jl IN.parquet OUT.parquet
#
# A raster-aligned Cartesian mesh (`--meshgen-cartesian dem.tif`) keeps the
# raster's nodata pixels as cells with a NaN elevation, which the model holds
# inert: their faces are walls. Dropped, the domain ends at the data, and the
# sides they shared become boundary sides (open unless --closed-boundaries).
# Neighbour lists are pruned to the cells kept; every other column is kept.

include(joinpath(@__DIR__, "..", "..", "mesh", "A5Grid.jl"))
using .A5Grid

length(ARGS) == 2 || error("usage: julia --project=. cases/scripts/drop_nodata_cells.jl IN.parquet OUT.parquet")
m = A5Grid.load_mesh_geoparquet(ARGS[1])
z = m.static_vars["elevation"]
keep = findall(isfinite, z)
ids = Set(m.cells[i].id for i in keep)
sv = Dict(k => v[keep] for (k, v) in m.static_vars)
av = Dict(k => v[:, keep] for (k, v) in m.array_vars)
adj = Dict(id => [nb for nb in nbs if nb in ids] for (id, nbs) in m.adjacency if id in ids)
out = A5Grid.A5Mesh(m.resolution, m.aoi_geojson, m.cells[keep], sv, av, adj, m.id_scheme, m.meta)
A5Grid.save_mesh_geoparquet(out, ARGS[2])
println(ARGS[2], ": kept ", length(keep), " of ", length(z), " cells")
