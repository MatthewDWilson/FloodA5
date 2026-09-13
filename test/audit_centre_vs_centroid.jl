#!/usr/bin/env julia
# test/audit_centre_vs_centroid.jl
#
# build_diamond_flux_table() (mesh/DiamondFlux.jl) uses raw
# c.center_lon/center_lat for cell positions. The B2/B3 validation
# scripts (audit_vertex_valence.jl, audit_diamond_gradient.jl) instead
# compute an explicit polygon centroid, per B1 Section 5's recommendation
# (Task 2 measured a ~1.3%-of-diameter offset between the two on an
# earlier mesh). This means the diamond reconstruction validated in the
# B2/B3 audits is not exactly the same operator running in production.
#
# This script quantifies that gap directly on a real mesh, and -- the
# important question -- checks whether it is a SYSTEMATIC (non-zero mean,
# not just non-zero scatter) directional bias, which would be a direct,
# independent contributor to any observed north/south or bearing-
# correlated asymmetry, separate from anything B1-B3 already addressed.
#
# Usage:
#   julia --project=. test/audit_centre_vs_centroid.jl <mesh.parquet>

using Statistics, Printf

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_path = length(ARGS) >= 1 ? ARGS[1] :
            joinpath("test", "square", "square_mesh_res14.parquet")

mesh  = load_mesh_geoparquet(mesh_path)
cells = mesh.cells
n     = length(cells)
println("Loaded $n cells from $mesh_path")

"""
    _area_centroid_local(boundary, lon0, lat0) -> (area, cx, cy)

Area-weighted polygon centroid (shoelace), computed in a LOCAL
equirectangular frame centred on (lon0, lat0). If (lon0,lat0) IS the
cell's own center_lon/center_lat, the returned (cx, cy) is exactly the
displacement vector (metres) from the A5 centre to the polygon centroid,
and `area` is the cell's true polygon area (used here for a diameter
figure, avoiding any dependency on how cell_area is stored in the mesh
struct -- computed independently instead).
"""
function _area_centroid_local(boundary::Vector{Vector{Float64}}, lon0::Float64, lat0::Float64)
    R = A5Grid._EARTH_R
    cl = cosd(lat0)
    to_xy(lon, lat) = (deg2rad(lon - lon0) * R * cl, deg2rad(lat - lat0) * R)
    pts = [to_xy(v[1], v[2]) for v in boundary]
    m = length(pts)
    Acc = 0.0; Cx = 0.0; Cy = 0.0
    j = m
    for i in 1:m
        xi, yi = pts[i]; xj, yj = pts[j]
        cross = xj * yi - xi * yj
        Acc += cross
        Cx  += (xj + xi) * cross
        Cy  += (yj + yi) * cross
        j = i
    end
    Acc *= 0.5
    abs(Acc) < 1e-9 && return (0.0, mean(first.(pts)), mean(last.(pts)))
    return (abs(Acc), Cx / (6Acc), Cy / (6Acc))
end

dx    = Vector{Float64}(undef, n)
dy    = Vector{Float64}(undef, n)
dmag  = Vector{Float64}(undef, n)
areas = Vector{Float64}(undef, n)
sublattice = Vector{Char}(undef, n)

for i in 1:n
    c = cells[i]
    area, cx, cy = _area_centroid_local(c.boundary, c.center_lon, c.center_lat)
    dx[i]   = cx
    dy[i]   = cy
    dmag[i] = sqrt(cx^2 + cy^2)
    areas[i] = area
    id_hex = A5Grid._to_hex(parse(UInt64, c.id, base=16))
    sublattice[i] = id_hex[1]
end

println("
--- Centre -> polygon-centroid displacement, whole mesh ---")
@printf("  n cells = %d\n", n)
@printf("  |displacement|: mean=%.4f m  median=%.4f m  p95=%.4f m  max=%.4f m\n",
        mean(dmag), median(dmag), quantile(dmag, 0.95), maximum(dmag))
@printf("  MEAN displacement vector: dx=%+.4f m  dy=%+.4f m\n", mean(dx), mean(dy))
println("  A non-zero MEAN here (not just scatter/std) is the important")
println("  number -- it means the centre-vs-centroid gap is a systematic")
println("  bias that will not average out, not random per-cell noise.")

mean_area = mean(areas)
mean_diam = 2 * sqrt(mean_area / pi)
@printf("
  Mean cell equivalent diameter (independently computed): %.2f m\n", mean_diam)
@printf("  Mean |displacement| as %% of diameter: %.3f%%\n", 100 * mean(dmag) / mean_diam)
println("  Compare against Task 2's original ~1.3% finding on an earlier mesh.")

println("
--- By sublattice (leading hex nibble) ---")
for sl in sort(unique(sublattice))
    idx = findall(==(sl), sublattice)
    @printf("  '%s'  n=%6d  mean dx=%+.4f  mean dy=%+.4f  mean|d|=%.4f\n",
            sl, length(idx), mean(dx[idx]), mean(dy[idx]), mean(dmag[idx]))
end
println("  If one sublattice's mean dx/dy is clearly non-zero while the")
println("  other's is near-zero (or they point in different directions),")
println("  that points at a sublattice-correlated construction difference")
println("  in how pya5 reports center_lon/center_lat -- directly relevant")
println("  to the still-open Task 1 question from the original Phase A work.")

println("
--- By latitude band ---")
lats = [c.center_lat for c in cells]
n_bands = 6
lo, hi = extrema(lats)
band_edges = range(lo, hi, length = n_bands + 1)
for b in 1:n_bands
    blo, bhi = band_edges[b], band_edges[b+1]
    idx = findall(i -> lats[i] >= blo && lats[i] < (b == n_bands ? bhi + 1e-9 : bhi), 1:n)
    isempty(idx) && continue
    @printf("  [%.4f, %.4f)  n=%6d  mean dx=%+.4f  mean dy=%+.4f\n",
            blo, bhi, length(idx), mean(dx[idx]), mean(dy[idx]))
end

bearing = mod.(atand.(dx, dy), 360.0)
is_north = (bearing .< 90.0) .| (bearing .> 270.0)
is_south = .!is_north
is_east  = (bearing .>= 0.0) .& (bearing .< 180.0)
is_west  = .!is_east
println("
--- Directional split of the displacement vectors themselves ---")
@printf("  north-pointing (n=%d): mean dy=%+.4f\n", count(is_north), mean(dy[is_north]))
@printf("  south-pointing (n=%d): mean dy=%+.4f\n", count(is_south), mean(dy[is_south]))
@printf("  east-pointing  (n=%d): mean dx=%+.4f\n", count(is_east),  mean(dx[is_east]))
@printf("  west-pointing  (n=%d): mean dx=%+.4f\n", count(is_west),  mean(dx[is_west]))
println("  This asks a different question than the latitude-band table above:")
println("  not 'where in the domain', but 'does the displacement vector itself")
println("  have a preferred compass direction, regardless of cell location'.")

println("
Done.")
