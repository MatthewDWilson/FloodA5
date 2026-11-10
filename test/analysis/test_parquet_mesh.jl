# Parquet mesh: load_mesh_parquet (format-1 fallback) against the format-2 /mesh.

_nbsets(m) = [Set(filter(!=(0), m.neighbours[:, i])) for i in 1:length(m)]
_idsets(m) = [Set(m.ids[j] for j in m.neighbours[:, i] if j != 0) for i in 1:length(m)]

"A minimal format-1 file (FloodA5 ≤ 0.3.x layout) with the fixture's cells in order `perm`."
function _write_v1(path, m, perm, ids, h)
    h5open(path, "w") do f
        g = create_group(f, "mesh")
        g["cell_ids"] = ids; g["center_lons"] = m.lon[perm]
        g["center_lats"] = m.lat[perm]; g["elevations"] = m.elevation[perm]
        fr = create_group(f, "frames")
        for k in 1:2
            gk = create_group(fr, lpad(string(k), 6, '0'))
            gk["t"] = 300.0k; gk["water_depth"] = h[perm]
        end
    end
    return path
end

@testset "A-parquet — load_mesh_parquet matches the format-2 mesh" begin
    m2 = FA.mesh_geometry(FIX_H5)
    mp = FA.load_mesh_parquet(FIX_MESH)
    @test mp.ids == m2.ids                       # the output's cell order is the parquet row order
    @test mp.lon == m2.lon && mp.lat == m2.lat
    @test isequal(mp.elevation, m2.elevation)
    @test mp.polygons == m2.polygons             # closed rings, as stored
    @test all(P -> P[:, 1] == P[:, end], mp.polygons)
    @test _nbsets(mp) == _nbsets(m2)
    @test mp.area == m2.area                     # the model's state.cell_area, bit for bit
    @test mp.id_scheme == m2.id_scheme == "cartesian"
    @test !isempty(mp.static)
    for (k, v) in mp.static
        @test haskey(m2.static, k) && isequal(v, m2.static[k])
    end
    # polygon_area_m2 is A5Grid._polygon_area_m2 bit for bit: fixture squares, and an
    # A5-like pentagon at Carlisle's latitude as a closed and as an open ring
    pent = [-2.952306 -2.951899 -2.956015 -2.957829 -2.954771 -2.952306;
            54.885633 54.888368 54.887691 54.885181 54.883360 54.885633]
    for P in (mp.polygons[1:50]..., pent, pent[:, 1:5])
        @test FA.polygon_area_m2(P) === A5Grid._polygon_area_m2([P[:, k] for k in 1:size(P, 2)])
    end
    @test FA.polygon_area_m2(pent[:, 1:2]) == 0.0
end

@testset "A-parquet — mesh_geometry(format 1; parquet) = format 2" begin
    m2 = FA.mesh_geometry(FIX_H5)
    n = length(m2)
    h = FO.read_frame(FIX_H5, "water_depth", FO.n_frames(FIX_H5))
    perm = vcat(2:2:n, 1:2:n)                    # the format-1 file in another cell order
    p = _write_v1(joinpath(FIX_DIR, "v1.h5"), m2, perm, m2.ids[perm], h)
    @test isempty(FA.mesh_geometry(p).area)      # without the parquet: the cell table only
    m1 = FA.mesh_geometry(p; parquet = FIX_MESH)
    @test m1.ids == m2.ids[perm] && m1.lon == m2.lon[perm]
    @test m1.polygons == m2.polygons[perm]
    @test m1.area == m2.area[perm]
    @test !isempty(m1.static) && all(k -> isequal(m1.static[k], m2.static[k][perm]), keys(m1.static))
    @test _idsets(m1) == _idsets(m2)[perm]       # neighbours, compared by id
    r1 = FA.load_run(p; parquet = FIX_MESH)
    @test r1.format == 1 && FA.nframes(r1) == 2
    @test FA.storage(r1) ≈ FA.storage(FA.load_run(FIX_H5)) rtol = 1e-12
    # a re-labelled copy of the mesh (e.g. pya5 0.x vs 1.x ids): ids differ, rows match
    q = _write_v1(joinpath(FIX_DIR, "v1_relabelled.h5"), m2, 1:n, ["x" * id for id in m2.ids], h)
    @test_throws ErrorException FA.mesh_geometry(q; parquet = FIX_MESH)
    mr = FA.mesh_geometry(q; parquet = FIX_MESH, match = :row)
    @test mr.polygons == m2.polygons && mr.area == m2.area && _nbsets(mr) == _nbsets(m2)
    @test_throws ErrorException FA.mesh_geometry(q; parquet = FIX_MESH, match = :nope)
end
