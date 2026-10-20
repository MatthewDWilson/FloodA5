# test/test_pya5_compat.jl
# ------------------------
# pya5 0.x / 1.x compatibility (WP-M1, 2026-10 — FloodA5_Review_and_ForwardPlan
# _2026-10.md §9). Runs against whichever pya5 is installed; CI runs it with
# both 0.7.3 and 1.x.
#
#   T1  pya5 version → A5 cell-ID scheme ("a5-v0" for 0.x, "a5-v1" for 1.x)
#   T2  Known cells (test/fixtures/a5_id_schemes.csv) resolve to their stored
#       centres and resolutions under the installed scheme; migrate_cell_id
#       maps every v0 ID to its v1 ID (1.x only)
#   T3  grid_disk_neighbours: 5 same-resolution, mutually adjacent,
#       edge-sharing neighbours (exercises the 1.x uncompact/marker path)
#   T4  GeoParquet round trip records the scheme; loading a mesh of the
#       other A5 scheme is refused (or warns with the override); Cartesian
#       meshes are not checked
#   T5  mesh/migrate_mesh_ids.py: v0 mesh → v1 mesh (1.x); refuses under 0.x
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_pya5_compat.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test

const V1 = startswith(A5Grid.pya5_version(), "1.")
const FIXTURE = let rows = NamedTuple[]
    for line in eachline(joinpath(@__DIR__, "fixtures", "a5_id_schemes.csv"))
        (startswith(line, "#") || startswith(line, "cell_id")) && continue
        f = split(line, ',')
        push!(rows, (v0 = String(f[1]), res = parse(Int, f[2]), lon = parse(Float64, f[3]),
                     lat = parse(Float64, f[4]), v1 = String(f[5])))
    end
    rows
end
_id(r) = V1 ? r.v1 : r.v0

function _fixture_mesh(rows; scheme = A5Grid.a5_id_scheme(), ids = _id.(rows))
    cells = [A5Cell(id, r.res, r.lon, r.lat, A5Grid.cell_to_boundary(_id(r)), NaN)
             for (id, r) in zip(ids, rows)]
    return A5Mesh(rows[1].res, "", cells, Dict{String,Vector{Float64}}(),
                  Dict{String,Matrix{Float64}}(), Dict{String,Vector{String}}(), scheme)
end

# ---------------------------------------------------------------------------
@testset "T1 — pya5 version and ID scheme (pya5 $(A5Grid.pya5_version()))" begin
    @test A5Grid.pya5_version() != "unknown"
    @test A5Grid.a5_id_scheme() == (V1 ? "a5-v1" : "a5-v0")
end

@testset "T2 — known cells resolve under the installed scheme" begin
    for r in FIXTURE
        lon, lat = A5Grid.cell_to_lonlat(_id(r))
        @test abs(lon - r.lon) < 1e-9 && abs(lat - r.lat) < 1e-9
        @test A5Grid.get_resolution(_id(r)) == r.res
        @test A5Grid.lonlat_to_cell(r.lon, r.lat, r.res) == _id(r)
    end
    if V1
        @test all(A5Grid.migrate_cell_id(r.v0) == r.v1 for r in FIXTURE)
        # The hazard the scheme check exists for: an old ID is still a valid
        # cell under 1.x, just a different one, far away.
        lon, lat = A5Grid.cell_to_lonlat(FIXTURE[1].v0)
        @test hypot(lon - FIXTURE[1].lon, lat - FIXTURE[1].lat) > 1.0
    else
        @test_throws ErrorException A5Grid.migrate_cell_id(FIXTURE[1].v0)
    end
end

@testset "T3 — grid_disk_neighbours" begin
    _vkey(v) = (round(v[1]; digits = 7), round(v[2]; digits = 7))
    for r in filter(r -> r.res == 18, FIXTURE)
        id  = _id(r)
        nbs = A5Grid.grid_disk_neighbours(id)
        @test length(nbs) == 5
        @test id ∉ nbs
        @test all(A5Grid.get_resolution(n) == 18 for n in nbs)
        @test all(id in A5Grid.grid_disk_neighbours(n) for n in nbs)
        vi = Set(_vkey.(A5Grid.cell_to_boundary(id)))
        @test all(length(intersect(vi, Set(_vkey.(A5Grid.cell_to_boundary(n))))) == 2 for n in nbs)
    end
end

@testset "T4 — GeoParquet ID-scheme metadata and the mismatch check" begin
    rows = filter(r -> r.res == 18, FIXTURE)
    dir  = mktempdir()
    p    = joinpath(dir, "mesh.parquet")
    save_mesh_geoparquet(_fixture_mesh(rows), p)
    m = load_mesh_geoparquet(p)
    @test m.id_scheme == A5Grid.a5_id_scheme()
    @test [c.id for c in m.cells] == _id.(rows)

    other = V1 ? "a5-v0" : "a5-v1"
    q = joinpath(dir, "other.parquet")
    save_mesh_geoparquet(_fixture_mesh(rows; scheme = other), q)
    @test_throws ErrorException load_mesh_geoparquet(q)
    m2 = withenv("FLOODA5_ALLOW_ID_SCHEME_MISMATCH" => "1") do
        load_mesh_geoparquet(q)
    end
    @test m2.id_scheme == other

    dlat = rad2deg(20.0 / A5Grid._EARTH_R)
    cart = generate_cartesian_mesh_bounds(0.0, 0.0, 2.5dlat, 2.5dlat, 20.0; lat_ref = 0.0)
    @test cart.id_scheme == "cartesian"
    c = joinpath(dir, "cart.parquet")
    save_mesh_geoparquet(cart, c)
    @test load_mesh_geoparquet(c).id_scheme == "cartesian"
end

@testset "T5 — mesh/migrate_mesh_ids.py" begin
    rows = filter(r -> r.res == 18, FIXTURE)
    dir  = mktempdir()
    src  = joinpath(dir, "v0.parquet"); dst = joinpath(dir, "v1.parquet")
    m0 = _fixture_mesh(rows; scheme = "a5-v0", ids = [r.v0 for r in rows])
    m0.adjacency[rows[1].v0] = [rows[2].v0, rows[3].v0]       # some neighbour IDs to map
    m0.adjacency[rows[2].v0] = [rows[1].v0]
    m0.adjacency[rows[3].v0] = [rows[1].v0]
    withenv("FLOODA5_ALLOW_ID_SCHEME_MISMATCH" => "1") do
        save_mesh_geoparquet(m0, src)
    end
    tool = joinpath(@__DIR__, "..", "mesh", "migrate_mesh_ids.py")
    cmd  = `$(A5Grid.PYTHON_CMD[]) $tool $src $dst`
    if V1
        run(pipeline(cmd; stdout = devnull))
        m1 = load_mesh_geoparquet(dst)
        @test m1.id_scheme == "a5-v1"
        @test [c.id for c in m1.cells] == [r.v1 for r in rows]
        @test sort(m1.adjacency[rows[1].v1]) == sort([rows[2].v1, rows[3].v1])
        @test [c.center_lon for c in m1.cells] == [c.center_lon for c in m0.cells]
        @test !success(pipeline(`$(A5Grid.PYTHON_CMD[]) $tool $dst $(dst * ".again")`;
                                stdout = devnull, stderr = devnull))   # already v1
    else
        @test !success(pipeline(cmd; stdout = devnull, stderr = devnull))  # needs pya5 >= 1
    end
end

println("\npya5 compatibility tests complete (pya5 $(A5Grid.pya5_version())).")
