# Mesh quality (analysis/mesh_quality.jl) on meshes with known answers: the
# fixture's 25 m Cartesian grid, a regular hexagonal patch and a sheared
# (rhombic) lattice, all near the equator.

"MeshGeom from (lon, lat) centres and 2 × nv rings: neighbours by shared vertices."
function _geom(lon, lat, polys)
    closed = [hcat(P, P[:, 1]) for P in polys]
    nb = FA._neighbours_shared_vertices(closed)
    area = [FA.polygon_area_m2(P) for P in closed]
    return FA.MeshGeom(["c$i" for i in eachindex(lon)], lon, lat, zeros(length(lon)), area,
                       closed, nb, Dict{String,Vector{Float64}}(), "polygon")
end

const _M2DEG = 1 / (deg2rad(1.0) * FA.R_EARTH)

"Pointy-top regular hexagons of side `s` metres, axial radius `rad`, centred at (lon0, lat0)."
function _hex_patch(s = 30.0, rad = 2; lon0 = 0.0, lat0 = 0.0)
    lon = Float64[]; lat = Float64[]; polys = Matrix{Float64}[]
    for q in -rad:rad, r in -rad:rad
        abs(q + r) <= rad || continue
        cx = s * sqrt(3) * (q + r / 2); cy = 1.5s * r
        kx = _M2DEG / cosd(lat0)
        P = [round(lon0 + (cx + s * cosd(30 + 60k)) * kx; digits = 13) for _ in 1:1, k in 0:5]
        Q = [round(lat0 + (cy + s * sind(30 + 60k)) * _M2DEG; digits = 13) for _ in 1:1, k in 0:5]
        push!(lon, lon0 + cx * kx); push!(lat, lat0 + cy * _M2DEG); push!(polys, vcat(P, Q))
    end
    return _geom(lon, lat, polys)
end

"Rhombi (sheared squares, angle `θ` between the lattice vectors), side `s` metres."
function _rhombi(s = 30.0, θ = 60.0, n = 4)
    a = (s, 0.0); b = (s * cosd(θ), s * sind(θ))
    lon = Float64[]; lat = Float64[]; polys = Matrix{Float64}[]
    pt(i, j) = (round((i * a[1] + j * b[1]) * _M2DEG; digits = 13),
                round((i * a[2] + j * b[2]) * _M2DEG; digits = 13))
    for i in 0:n-1, j in 0:n-1
        vs = [pt(i, j), pt(i + 1, j), pt(i + 1, j + 1), pt(i, j + 1)]
        push!(polys, [v[k] for k in 1:2, v in vs])
        c = pt(i + 0.5, j + 0.5); push!(lon, c[1]); push!(lat, c[2])
    end
    return _geom(lon, lat, polys)
end

@testset "Q — mesh quality" begin
    @testset "Cartesian (fixture)" begin
        q = FA.mesh_quality(FA.load_mesh_parquet(FIX_MESH))
        @test collect(keys(q.valence)) == [4]
        @test q.nonorth_max < 1e-6 && q.skew_max < 1e-6
        @test q.orientation[1] < 1e-6 && q.orientation[2] > 1 - 1e-9   # m2 = 0, m4 = 1
        @test collect(keys(q.vertex_valence)) == [4]
        @test q.centre_offset_max < 1e-6
        @test q.area_cv < 1e-3
    end
    @testset "regular hexagons" begin
        q = FA.mesh_quality(_hex_patch())
        @test q.n_cells == 19 && q.n_interior == 7
        @test q.valence == Dict(6 => 7)
        @test q.n_faces == 42
        @test q.nonorth_max < 1e-4 && q.skew_max < 1e-6
        m2, m4, m6 = q.orientation
        @test m2 < 1e-6 && m4 < 1e-6 && m6 > 1 - 1e-6
        @test collect(keys(q.vertex_valence)) == [3]
        @test q.centre_offset_max < 1e-6
        @test q.area_cv < 1e-6
    end
    @testset "regular hexagons at 55° N (centroid precision)" begin
        # absolute coordinates of ~1e-4° cells: the centroid must not cancel
        q = FA.mesh_quality(_hex_patch(15.0; lon0 = -2.9, lat0 = 54.9))
        @test q.centre_offset_max < 1e-6
        @test q.orientation[3] > 1 - 1e-4
    end
    @testset "rhombi (60°)" begin
        q = FA.mesh_quality(_rhombi())
        @test q.valence == Dict(4 => 4)
        # the centre line of a 60° rhombus lattice is 30° off the face normal,
        # and crosses the face at its midpoint
        @test all(isapprox.(q.faces.nonorth_deg, 30.0; atol = 1e-6))
        @test q.skew_max < 1e-6
    end
end
