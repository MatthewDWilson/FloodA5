#!/usr/bin/env julia
# audit_diamond_gradient.jl — Phase B, B3 empirical verification.
#
# Usage: julia --project=. test/audit_diamond_gradient.jl <mesh.parquet>

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

_norm_id(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
ids    = [_norm_id(c.id) for c in cells]
id_idx = Dict{String,Int}(ids[i] => i for i in 1:n)

adj = if !isempty(mesh.adjacency)
    Dict{String,Vector{String}}(_norm_id(k) => [_norm_id(v) for v in vs]
                                 for (k, vs) in mesh.adjacency)
else
    _build_adjacency_shared_vertices(mesh)
end

# ── Local equirectangular frame, global to the whole mesh (adequate at
#    these cell scales — same convention used throughout this session) ──
lon0, lat0 = mean(c.center_lon for c in cells), mean(c.center_lat for c in cells)
function to_xy(lon, lat)
    R = A5Grid._EARTH_R
    cl = cosd(lat0)
    return (deg2rad(lon-lon0)*R*cl, deg2rad(lat-lat0)*R)
end

# ── Polygon centroid (B1 §5's recommendation, not raw center_lon/lat) ──
function centroid_xy(boundary::Vector{Vector{Float64}})
    m = length(boundary)
    pts = [to_xy(v[1], v[2]) for v in boundary]
    A = 0.0; Cx = 0.0; Cy = 0.0
    j = m
    for i in 1:m
        xi, yi = pts[i]; xj, yj = pts[j]
        cross = xj*yi - xi*yj
        A += cross; Cx += (xj+xi)*cross; Cy += (yj+yi)*cross
        j = i
    end
    A *= 0.5
    abs(A) < 1e-9 && return mean(first.(pts)), mean(last.(pts))
    return Cx/(6A), Cy/(6A)
end
cell_xy = [centroid_xy(c.boundary) for c in cells]

# ── B2 vertex reconstruction (re-derived here, matching the B2 memo) ──
_key(lon, lat) = (round(lon, digits=9), round(lat, digits=9))
vertex_cells = Dict{Tuple{Float64,Float64}, Vector{Int}}()
for (i, c) in enumerate(cells)
    for pt in c.boundary
        k = _key(pt[1], pt[2])
        push!(get!(vertex_cells, k, Int[]), i)
    end
end
for (k, v) in vertex_cells
    vertex_cells[k] = unique(v)
end

function vertex_weights_and_cond(dx::Matrix{Float64}, w::Vector{Float64})
    k = size(dx, 1)
    Wd = Diagonal(w)
    Amat = hcat(ones(k), dx)
    M = Amat' * Wd * Amat
    R = Amat' * Wd
    sol = M \ R
    return sol[1, :], cond(M)
end

"""Reconstruct eta_v at vertex key `vk`, for a given per-cell eta vector."""
function eta_at_vertex(vk, eta_cell::Vector{Float64})
    cidx = vertex_cells[vk]
    k = length(cidx)
    k < 3 && return NaN, Inf   # fallback cases not exercised here, see B2 memo
    vlon, vlat = vk
    vx, vy = to_xy(vlon, vlat)
    dx = Matrix{Float64}(undef, k, 2)
    w  = Vector{Float64}(undef, k)
    for (row, ci) in enumerate(cidx)
        cx, cy = cell_xy[ci]
        dx[row,1] = cx - vx; dx[row,2] = cy - vy
        d2 = dx[row,1]^2 + dx[row,2]^2
        w[row] = d2 < 1e-6 ? 1e6 : 1.0/d2
    end
    weights, cnum = vertex_weights_and_cond(dx, w)
    return dot(weights, eta_cell[cidx]), cnum
end

# ── Diamond gradient (this document, §2) ──
function diamond_gradient(P::Matrix{Float64}, eta::Vector{Float64})
    k = 4
    area2 = 0.0
    S = zeros(2)
    for a in 1:k
        b = a % k + 1
        xa, ya = P[a,1], P[a,2]
        xb, yb = P[b,1], P[b,2]
        area2 += xa*yb - xb*ya
        eta_avg = 0.5*(eta[a]+eta[b])
        dxe, dye = xb-xa, yb-ya
        S[1] += eta_avg*dye
        S[2] += eta_avg*(-dxe)
    end
    area = area2/2.0
    return S ./ area, area
end

# ── Synthetic linear field over the whole mesh ──
A_true = 4.1
B_true = (0.0018, -0.0027)
eta_cell = [A_true + B_true[1]*cell_xy[i][1] + B_true[2]*cell_xy[i][2] for i in 1:n]

# ── Walk real edges, build each diamond, compare methods ──
areas = Float64[]
diamond_errs = Float64[]
vertex_conds = Float64[]
n_degenerate_diamond = 0
n_vertex_fallback = 0

seen = Set{Tuple{Int,Int}}()
for i in 1:n
    for nb_id in get(adj, ids[i], String[])
        j = get(id_idx, nb_id, 0)
        j == 0 && continue
        lo, hi = minmax(i, j)
        (lo, hi) in seen && continue
        push!(seen, (lo, hi))

        bi, bj = cells[lo], cells[hi]
        edge = A5Grid._shared_edge(bi.boundary, bj.boundary)
        edge === nothing && continue
        elon1, elat1, elon2, elat2 = edge
        va_key = _key(elon1, elat1)
        vb_key = _key(elon2, elat2)
        (!haskey(vertex_cells, va_key) || !haskey(vertex_cells, vb_key)) && continue

        eta_va, cond_a = eta_at_vertex(va_key, eta_cell)
        eta_vb, cond_b = eta_at_vertex(vb_key, eta_cell)
        isfinite(cond_a) && push!(vertex_conds, cond_a)
        isfinite(cond_b) && push!(vertex_conds, cond_b)
        if !isfinite(eta_va) || !isfinite(eta_vb)
            global n_vertex_fallback += 1
            continue
        end

        xi_, yi_ = cell_xy[lo]
        xj_, yj_ = cell_xy[hi]
        vax, vay = to_xy(elon1, elat1)
        vbx, vby = to_xy(elon2, elat2)

        P = [xi_ yi_; vax vay; xj_ yj_; vbx vby]
        eta4 = [eta_cell[lo], eta_va, eta_cell[hi], eta_vb]

        grad, area = diamond_gradient(P, eta4)
        push!(areas, abs(area))
        if abs(area) < 1.0
            global n_degenerate_diamond += 1
            continue
        end

        true_grad = [B_true[1], B_true[2]]
        push!(diamond_errs, norm(grad - true_grad))
    end
end

println("\nEdges processed: $(length(seen))")
println("Vertex-reconstruction fallback (k<3) hit on: $n_vertex_fallback edges")
println("Degenerate diamond (|area|<1 m^2) on: $n_degenerate_diamond edges")

println("\n--- Diamond area distribution (m^2) ---")
println("  min=$(round(minimum(areas),digits=1))  median=$(round(median(areas),digits=1))  " *
        "p05=$(round(quantile(areas,0.05),digits=1))  max=$(round(maximum(areas),digits=1))")

println("\n--- Vertex-stencil conditioning (vertices feeding a real diamond) ---")
println("  n stencils: $(length(vertex_conds))")
println("  min=$(round(minimum(vertex_conds),digits=1))  median=$(round(median(vertex_conds),digits=1))  " *
        "p95=$(round(quantile(vertex_conds,0.95),digits=1))  max=$(round(maximum(vertex_conds),digits=1))")
println("  (compare against B2's whole-mesh check — this is the same quantity,")
println("   restricted to vertices that actually feed a real edge's diamond)")

println("\n--- End-to-end linear-exactness (B2 vertex recon -> B3 diamond gradient) ---")
println("  n edges tested: $(length(diamond_errs))")
println("  max |gradient error|: $(maximum(diamond_errs))")
println("  median |gradient error|: $(median(diamond_errs))")
println("  (compare against the synthetic sweep in this document: well-conditioned")
println("   configurations should be at or near machine precision)")

println("\nDone.")