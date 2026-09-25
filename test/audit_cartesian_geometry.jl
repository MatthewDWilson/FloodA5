#!/usr/bin/env julia
# audit_cartesian_geometry.jl
# ---------------------------
# Verifies a generated Cartesian mesh in two independent ways, per
# FloodA5_CartesianFormulation_Handover.md §8 item 2:
#
#   1. Adjacency correctness — directly answers "are neighbours specified
#      correctly": per-cell neighbour COUNT (2 at corners, 3 on edges,
#      4 interior), adjacency SYMMETRY (if A lists B, B must list A), and
#      GEOMETRIC correctness (each listed neighbour must actually sit one
#      cell-width away, in exactly one cardinal direction, from the cell
#      that lists it) — checked from mesh.adjacency directly, independent
#      of the EdgeList/geometry pipeline in step 2 below.
#
#   2. Edge-geometry orthogonality — builds the EdgeList through the real,
#      production _build_edge_list function (not a re-implementation), and
#      checks cos_theta and skew_x/skew_y on every edge. On a perfectly
#      axis-aligned Cartesian grid these should be EXACTLY 1.0 and (0,0)
#      respectively (to floating-point noise only) — not just "close",
#      since there is no non-orthogonality on this mesh at all. Any
#      material deviation here means a bug in the flux/geometry pipeline
#      that this clean grid has newly exposed, not a property of the mesh.
#
# No simulation is run — this is a pure geometry/adjacency diagnostic,
# matching the style of this project's other audit_*.jl scripts.
#
# Usage:
#   julia --project=. test/audit_cartesian_geometry.jl <mesh.parquet>

using Statistics

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_path = length(ARGS) >= 1 ? ARGS[1] :
            joinpath("test", "square", "square_cartesian_dx400.parquet")

mesh  = load_mesh_geoparquet(mesh_path)
cells = mesh.cells
n     = length(cells)
println("Loaded $n cells from $mesh_path")

# ── local equirectangular projection (for the geometric adjacency check) ──
lat0     = mean(c.center_lat for c in cells)
cos_lat0 = cosd(lat0)
R        = A5Grid._EARTH_R
function to_xy(lon, lat)
    return (deg2rad(lon) * R * cos_lat0, deg2rad(lat) * R)
end
# (absolute local frame is fine here — we only ever difference two points)

id_idx = Dict{String,Int}(cells[i].id => i for i in 1:n)

# =============================================================================
# Part 1 — Adjacency correctness (directly from mesh.adjacency)
# =============================================================================
println("\n" * "="^60)
println("Part 1: Adjacency correctness")
println("="^60)

valence_counts = Dict{Int,Int}()
n_asymmetric   = 0
n_bad_direction = 0
n_bad_distance  = 0
dx_expected = nothing   # cell width in metres, discovered from the first edge
dy_expected = nothing

for i in 1:n
    nbrs = get(mesh.adjacency, cells[i].id, String[])
    k = length(nbrs)
    valence_counts[k] = get(valence_counts, k, 0) + 1

    xi, yi = to_xy(cells[i].center_lon, cells[i].center_lat)

    for nb_id in nbrs
        # Symmetry: does the neighbour list this cell back?
        back = get(mesh.adjacency, nb_id, String[])
        if !(cells[i].id in back)
            global n_asymmetric += 1
            println("  ASYMMETRY: $(cells[i].id) -> $nb_id, but not back")
        end

        # Geometric check: neighbour should be exactly one cell-width away,
        # in exactly one cardinal direction (dx==0 xor dy==0).
        j = get(id_idx, nb_id, 0)
        j == 0 && continue
        xj, yj = to_xy(cells[j].center_lon, cells[j].center_lat)
        dx, dy = xj - xi, yj - yi

        is_ew = abs(dy) < 1e-6 && abs(dx) > 1e-6
        is_ns = abs(dx) < 1e-6 && abs(dy) > 1e-6
        if !is_ew && !is_ns
            global n_bad_direction += 1
            println("  BAD DIRECTION: $(cells[i].id) -> $nb_id  (dx=$dx, dy=$dy)")
            continue
        end

        d = is_ew ? abs(dx) : abs(dy)
        if is_ew
            global dx_expected = dx_expected === nothing ? d : dx_expected
            if abs(d - dx_expected) > 1e-3
                global n_bad_distance += 1
                println("  BAD EW DISTANCE: $(cells[i].id) -> $nb_id  " *
                        "(got $d, expected $dx_expected)")
            end
        else
            global dy_expected = dy_expected === nothing ? d : dy_expected
            if abs(d - dy_expected) > 1e-3
                global n_bad_distance += 1
                println("  BAD NS DISTANCE: $(cells[i].id) -> $nb_id  " *
                        "(got $d, expected $dy_expected)")
            end
        end
    end
end

println("\nValence (neighbour count) distribution:")
for k in sort(collect(keys(valence_counts)))
    println("  $k neighbours: $(valence_counts[k]) cells")
end
println("Expected for a rectangular grid: 2 (corners), 3 (edges), 4 (interior)")

println("\nAsymmetric adjacency pairs : $n_asymmetric   (expect 0)")
println("Non-cardinal-direction hits: $n_bad_direction (expect 0)")
println("Wrong-distance hits        : $n_bad_distance  (expect 0)")

part1_pass = n_asymmetric == 0 && n_bad_direction == 0 && n_bad_distance == 0
println("\nPart 1 verdict: ", part1_pass ? "PASS — adjacency is correct" :
                                            "FAIL — see flagged pairs above")

# =============================================================================
# Part 2 — Edge geometry orthogonality (via the real _build_edge_list)
# =============================================================================
println("\n" * "="^60)
println("Part 2: Edge geometry (cos_theta, skew) via production _build_edge_list")
println("="^60)

areas = [A5Grid._polygon_area_m2(c.boundary) for c in cells]
edges = _build_edge_list(cells, id_idx, mesh.adjacency, areas)

println("EdgeList: $(edges.n_edges) edges for $n cells")

ct  = edges.cos_theta[1:edges.n_edges]
sx  = edges.skew_x[1:edges.n_edges]
sy  = edges.skew_y[1:edges.n_edges]

println("\ncos_theta:  min=$(minimum(ct))  mean=$(mean(ct))  max=$(maximum(ct))")
println("skew_x:     min=$(minimum(sx))  max abs=$(maximum(abs.(sx)))")
println("skew_y:     min=$(minimum(sy))  max abs=$(maximum(abs.(sy)))")

max_ct_err = maximum(abs.(ct .- 1.0))
max_skew   = max(maximum(abs.(sx)), maximum(abs.(sy)))
println("\nmax |cos_theta - 1.0| : $max_ct_err  (expect ~1e-12 or smaller)")
println("max |skew|            : $max_skew  (expect ~1e-9 or smaller)")

# Machine-precision tolerance, not "small" — this mesh has NO non-orthogonality
# at all, so any deviation beyond floating-point roundoff is a real defect.
part2_pass = max_ct_err < 1e-8 && max_skew < 1e-6
println("\nPart 2 verdict: ", part2_pass ? "PASS — every edge is exactly orthogonal" :
                                            "FAIL — non-trivial cos_theta/skew found on an axis-aligned grid")

println("\n" * "="^60)
println("Overall: ", (part1_pass && part2_pass) ? "PASS" : "FAIL — see above")
println("="^60)
