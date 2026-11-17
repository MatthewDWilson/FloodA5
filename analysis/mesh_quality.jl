# analysis/mesh_quality.jl — geometric quality of a polygon mesh.
#
# Metrics that describe how a mesh family looks to a face-based finite-volume
# scheme, for comparing A5, Cartesian, the polygon families
# (mesh/mesh_families.py) and the DGGS families (mesh/dggs_families.py):
#
#   valence           neighbours per interior cell (all polygon sides shared)
#   area              mean, coefficient of variation, max/min
#   non-orthogonality angle between the face normal and the line joining the
#                     two cell centres (0° = orthogonal; the model's
#                     acos(cos_theta))
#   skewness          where that line crosses the face, as the distance from
#                     the face midpoint over the face length (0 = midpoint)
#   centre offset     declared cell centre to polygon centroid, over √area
#                     (Phase A Task 2's control-volume-centre check)
#   face orientation  azimuths of the face normals (mod 180°) and their
#                     width-weighted moments |Σ w e^{ikφ}| / Σ w for
#                     k = 2, 4, 6: a square lattice has m4 = 1 (four face
#                     directions), a regular hexagonal one m6 = 1 and
#                     m2 = m4 = 0; a direction-free mesh has all three → 0
#   vertex valence    cells meeting at each vertex whose cells are all interior
#
# Faces are projected about their own mid-latitude, as the model does
# (A5Grid._edge_geometry), so the non-orthogonality is the model's. Pure
# functions of arrays (`face_geometry`, `orientation_moments`) plus a
# MeshGeom wrapper (`mesh_quality`).

export face_geometry, orientation_moments, mesh_quality

"Index pairs (a, b) of the vertices `P[:, a]` and `Q[:, b]` that coincide (to `tol` degrees)."
function _shared_vertices(P::AbstractMatrix, Q::AbstractMatrix; tol = 1e-7)
    np_ = _nvert(P); nq = _nvert(Q)
    out = Tuple{Int,Int}[]
    for a in 1:np_, b in 1:nq
        abs(P[1, a] - Q[1, b]) <= tol && abs(P[2, a] - Q[2, b]) <= tol && push!(out, (a, b))
    end
    return out
end

"Number of distinct vertices of a (possibly closed) 2 × nv ring."
_nvert(P) = (nv = size(P, 2); nv > 1 && P[1, 1] == P[1, nv] && P[2, 1] == P[2, nv] ? nv - 1 : nv)

"""
    face_geometry(polygons, lon, lat, neighbours) -> NamedTuple

Per face (each neighbouring pair i < j once): `i`, `j`, `width` (m),
`nonorth_deg`, `skew` (crossing point's distance from the face midpoint /
width), `azimuth_deg` (face-normal direction mod 180°, from east,
anticlockwise), and per cell `n_faces` and `interior` (every polygon side is
a face). `polygons` are 2 × nv (lon, lat) rings, `lon`/`lat` the declared
cell centres, `neighbours` a max_nb × n matrix of 1-based indices (0 = none).
"""
function face_geometry(polygons::AbstractVector{<:AbstractMatrix}, lon::AbstractVector,
                       lat::AbstractVector, neighbours::AbstractMatrix{<:Integer})
    n = length(polygons)
    fi = Int[]; fj = Int[]; w = Float64[]; no = Float64[]; sk = Float64[]; az = Float64[]
    n_faces = zeros(Int, n)
    for i in 1:n, s in axes(neighbours, 1)
        j = Int(neighbours[s, i])
        (j == 0 || j <= i) && continue
        sv = _shared_vertices(polygons[i], polygons[j])
        length(sv) == 2 || continue                  # touching at a corner only
        P = polygons[i]
        (a1, _), (a2, _) = sv
        lat0 = (P[2, a1] + P[2, a2]) / 2
        k = cosd(lat0)
        xy(lo, la) = (deg2rad(lo) * R_EARTH * k, deg2rad(la) * R_EARTH)
        x1, y1 = xy(P[1, a1], P[2, a1]); x2, y2 = xy(P[1, a2], P[2, a2])
        cix, ciy = xy(lon[i], lat[i]); cjx, cjy = xy(lon[j], lat[j])
        tx, ty = x2 - x1, y2 - y1
        L = hypot(tx, ty)
        nx, ny = ty / L, -tx / L
        dx, dy = cjx - cix, cjy - ciy
        d = hypot(dx, dy)
        c = clamp(abs(dx * nx + dy * ny) / d, 0.0, 1.0)
        # crossing of ci + a·d with v1 + s·t:  a·d − s·t = v1 − ci
        det = dx * (-ty) - dy * (-tx)
        sfrac = abs(det) < eps() * d * L ? NaN :
                (dx * (y1 - ciy) - dy * (x1 - cix)) / det
        push!(fi, i); push!(fj, j); push!(w, L); push!(no, acosd(c))
        push!(sk, abs(sfrac - 0.5)); push!(az, mod(atand(ny, nx), 180.0))
        n_faces[i] += 1; n_faces[j] += 1
    end
    interior = [n_faces[i] == _nvert(polygons[i]) for i in 1:n]
    return (; i = fi, j = fj, width = w, nonorth_deg = no, skew = sk, azimuth_deg = az,
            n_faces, interior)
end

"""
    orientation_moments(azimuth_deg, width; k = (2, 4, 6)) -> Vector{Float64}

Width-weighted moments |Σ w e^{ikφ}| / Σ w of face-normal azimuths φ (any
k even, since φ is defined mod 180°). 1 when every face normal lies in a
lattice of k/2 directions spaced 360°/k apart (square: k = 4; hexagon:
k = 6); 0 for a direction-free distribution.
"""
function orientation_moments(azimuth_deg::AbstractVector, width::AbstractVector; k = (2, 4, 6))
    W = sum(width)
    return [abs(sum(width[f] * cis(kk * deg2rad(azimuth_deg[f])) for f in eachindex(width))) / W
            for kk in k]
end

"""Polygon centroid (lon, lat) by the shoelace formula on a local
equirectangular plane. Coordinates are taken relative to the first vertex:
in absolute degrees the cross products cancel catastrophically (a 25 m
hexagon at 55° N came out 0.9 m off its centre)."""
function _centroid(P::AbstractMatrix)
    nv = _nvert(P)
    k = cosd(sum(@view P[2, 1:nv]) / nv)
    o1, o2 = P[1, 1], P[2, 1]
    A = 0.0; cx = 0.0; cy = 0.0
    for a in 1:nv
        b = a == nv ? 1 : a + 1
        xa, ya = (P[1, a] - o1) * k, P[2, a] - o2; xb, yb = (P[1, b] - o1) * k, P[2, b] - o2
        cr = xa * yb - xb * ya
        A += cr; cx += (xa + xb) * cr; cy += (ya + yb) * cr
    end
    return cx / (3A) / k + o1, cy / (3A) + o2
end

_q(v, p) = isempty(v) ? NaN : quantile(filter(isfinite, v), p)

"""
    mesh_quality(g::MeshGeom) -> NamedTuple

Summary and per-face geometry of a mesh (see the file header): `n_cells`,
`n_interior`, `n_faces`, `valence` (Dict valence => count, interior cells),
`area_mean`, `area_cv`, `area_ratio`, `nonorth_median`/`_p95`/`_max` (°),
`skew_median`/`_max`, `centre_offset_max` (/√area, interior cells),
`orientation` ([m2, m4, m6]), `vertex_valence` (Dict), and `faces` (the
`face_geometry` arrays). Needs polygons and neighbours (format-2 output or
`load_mesh_parquet`).
"""
function mesh_quality(g::MeshGeom)
    isempty(g.polygons) && error("mesh_quality: the mesh has no polygons")
    F = face_geometry(g.polygons, g.lon, g.lat, g.neighbours)
    inter = findall(F.interior)
    val = Dict{Int,Int}()
    for i in inter; val[F.n_faces[i]] = get(val, F.n_faces[i], 0) + 1; end
    a = filter(x -> isfinite(x) && x > 0, g.area)
    off = Float64[]
    for i in inter
        cl, ct = _centroid(g.polygons[i])
        dx = deg2rad(g.lon[i] - cl) * R_EARTH * cosd(ct); dy = deg2rad(g.lat[i] - ct) * R_EARTH
        push!(off, hypot(dx, dy) / sqrt(g.area[i]))
    end
    # Vertex valence over vertices whose cells are all interior: at the domain
    # edge a cell can be missing from around a vertex of an interior cell.
    vkey(P, a) = (round(Int64, P[1, a] * 1e7), round(Int64, P[2, a] * 1e7))
    vcells = Dict{Tuple{Int64,Int64},Vector{Int}}()
    for (i, P) in enumerate(g.polygons), a in 1:_nvert(P)
        push!(get!(vcells, vkey(P, a), Int[]), i)
    end
    vval = Dict{Int,Int}()
    for cs in values(vcells)
        all(i -> F.interior[i], cs) || continue
        vval[length(cs)] = get(vval, length(cs), 0) + 1
    end
    return (; n_cells = length(g.polygons), n_interior = length(inter), n_faces = length(F.width),
            valence = val, area_mean = mean(a), area_cv = std(a) / mean(a),
            area_ratio = maximum(a) / minimum(a),
            nonorth_median = _q(F.nonorth_deg, 0.5), nonorth_p95 = _q(F.nonorth_deg, 0.95),
            nonorth_max = isempty(F.nonorth_deg) ? NaN : maximum(F.nonorth_deg),
            skew_median = _q(F.skew, 0.5), skew_max = _q(F.skew, 1.0),
            centre_offset_max = isempty(off) ? NaN : maximum(off),
            orientation = orientation_moments(F.azimuth_deg, F.width),
            vertex_valence = vval, faces = F)
end
