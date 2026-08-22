#!/usr/bin/env julia
# audit_vertex_valence.jl — Phase B, B2 empirical verification.
#
# Usage: julia --project=. test/audit_vertex_valence.jl <mesh.parquet>

using Statistics, LinearAlgebra

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_path = length(ARGS) >= 1 ? ARGS[1] :
            joinpath("test", "square", "square_mesh_res14.parquet")

mesh  = load_mesh_geoparquet(mesh_path)
cells = mesh.cells
n     = length(cells)
println("Loaded $n cells from $mesh_path")

# ── Build vertex -> list of (cell_index, local equirect position) ─────────
# Vertices are matched by rounding lon/lat to a small tolerance, since
# floating point boundary coordinates from adjacent cells should coincide
# but may differ in the last bit or two.
_key(lon, lat) = (round(lon, digits=9), round(lat, digits=9))

vertex_cells = Dict{Tuple{Float64,Float64}, Vector{Int}}()
for (i, c) in enumerate(cells)
    for pt in c.boundary
        k = _key(pt[1], pt[2])
        push!(get!(vertex_cells, k, Int[]), i)
    end
end
# Keep only vertices where >=2 DISTINCT cells meet (dedupe repeated cell
# hits from a cell visiting the same physical vertex more than once, which
# shouldn't happen for a simple polygon boundary, but guard anyway).
for (k, v) in vertex_cells
    vertex_cells[k] = unique(v)
end
real_vertices = [(k,v) for (k,v) in vertex_cells if length(v) >= 2]
println("Distinct vertices with >=2 cells: $(length(real_vertices))")

# ── Task 1: valence distribution ───────────────────────────────────────────
valences = [length(v) for (_,v) in real_vertices]
valence_counts = Dict{Int,Int}()
for k in valences
    valence_counts[k] = get(valence_counts, k, 0) + 1
end
println("\nVertex valence distribution:")
for k in sort(collect(keys(valence_counts)))
    frac = valence_counts[k] / length(valences)
    println("  valence $k : $(valence_counts[k])  ($(round(100*frac,digits=1))%)")
end
println("Euler's-formula prediction (large-mesh limit, all-pentagon tiling):")
println("  ~66.7% valence 3, ~33.3% valence 4, average 3.333")
println("  Actual average: $(round(mean(valences), digits=3))")

# ── Local equirectangular helper ────────────────────────────────────────────
function to_xy(lon, lat, lon0, lat0)
    R = A5Grid._EARTH_R
    cl = cosd(lat0)
    return (deg2rad(lon-lon0)*R*cl, deg2rad(lat-lat0)*R)
end

# ── Task 2 & 3: conditioning + linear exactness on REAL vertex stencils ────
"""
    _vertex_weights(dx, w) -> weights (length k), cond(M)

WLS constant-term weight vector, per §2.1 of the B2 derivation.
"""
function _vertex_weights(dx::Matrix{Float64}, w::Vector{Float64})
    k = size(dx, 1)
    Wd = Diagonal(w)
    Amat = hcat(ones(k), dx)           # k x 3
    M = Amat' * Wd * Amat              # 3x3
    R = Amat' * Wd                     # 3 x k
    sol = M \ R                        # 3 x k
    return sol[1, :], cond(M)
end

conds = Float64[]
const_errs = Float64[]
lin_errs = Float64[]

# synthetic linear field over the WHOLE mesh, in a global local frame
lon0, lat0 = mean(c.center_lon for c in cells), mean(c.center_lat for c in cells)
A_true = 3.7
B_true = (0.0021, -0.0034)  # per-metre gradient components (x,y)

cell_xy = [to_xy(c.center_lon, c.center_lat, lon0, lat0) for c in cells]
cell_eta_const = fill(5.0, n)
cell_eta_lin = [A_true + B_true[1]*xy[1] + B_true[2]*xy[2] for xy in cell_xy]

n_skipped_degenerate = 0
for (vkey, cidx) in real_vertices
    k_valence = length(cidx)
    k_valence < 3 && continue   # k=2/1 fallback cases handled separately, not exactness-tested here

    vlon, vlat = vkey
    vx, vy = to_xy(vlon, vlat, lon0, lat0)

    dx = Matrix{Float64}(undef, k_valence, 2)
    w  = Vector{Float64}(undef, k_valence)
    for (row, ci) in enumerate(cidx)
        cx, cy = cell_xy[ci]
        dx[row,1] = cx - vx
        dx[row,2] = cy - vy
        d2 = dx[row,1]^2 + dx[row,2]^2
        w[row] = d2 < 1e-6 ? 1e6 : 1.0/d2
    end

    weights, cnum = _vertex_weights(dx, w)
    if !isfinite(cnum) || cnum > 1e14
        n_skipped_degenerate += 1
        continue
    end
    push!(conds, cnum)

    eta_v_const = dot(weights, cell_eta_const[cidx])
    push!(const_errs, abs(eta_v_const - 5.0))

    true_eta_v = A_true + B_true[1]*vx + B_true[2]*vy
    eta_v_lin = dot(weights, cell_eta_lin[cidx])
    push!(lin_errs, abs(eta_v_lin - true_eta_v))
end

println("\n--- Conditioning (real mesh vertices, k>=3) ---")
println("  n vertices tested: $(length(conds))  (skipped as degenerate: $n_skipped_degenerate)")
println("  condition number: min=$(round(minimum(conds),digits=1))  median=$(round(median(conds),digits=1))  " *
        "p95=$(round(quantile(conds,0.95),digits=1))  max=$(round(maximum(conds),digits=1))")

println("\n--- Constant preservation (real mesh vertices) ---")
println("  max error: $(maximum(const_errs))")

println("\n--- Linear exactness (real mesh vertices) ---")
println("  max error: $(maximum(lin_errs))")
println("  (compare against the synthetic 20,000-trial sweep in the B2 memo:")
println("   well-conditioned -> ~1e-14/1e-16; only ill-conditioned stencils")
println("   (cond > ~1e11) should show larger errors)")

println("\nDone.")