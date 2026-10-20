# test/test_mesh_families.jl
# --------------------------
# Tier 0 on generic polygon mesh families (WP3, 2026-10 — FloodA5_Review_and_
# ForwardPlan_2026-10.md §6.1, §9 WP3). Meshes come from mesh/mesh_families.py:
# regular hexagons, rotated Cartesian, centroidal Voronoi (generator centres —
# orthogonal) and jittered Voronoi with centroid centres (non-orthogonal).
#
#   G   geometry: every interior cell closes (Σ w·n̂ = 0); widths and centre
#       distances positive; valence as expected; faces orthogonal to the
#       centre-to-centre line exactly where the family guarantees it
#   L   lake at rest: a still, non-flat-bed water surface stays still under
#       the legacy, WLSQ-corrected, diamond and cell-momentum options
#   M   mass conservation: closed point-source run, exact volume balance
#   I   import: exporting a generated mesh's polygons and re-importing them
#       gives the same cells and adjacency
#   R   raster-aligned projected Cartesian grid (generate_cartesian_mesh_raster)
#       on a GeoTIFF in British National Grid: one cell per pixel, centres on
#       the pixel lattice, elevation and orientation from the raster
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_mesh_families.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

const GEN  = joinpath(@__DIR__, "..", "mesh", "mesh_families.py")
const BBOX = "0.0,0.0,0.006,0.006"          # ≈ 670 m square at the equator
const DIR  = mktempdir()

A5Grid._find_python()
function _family(name, args...)
    out = joinpath(DIR, "$name.parquet")
    run(pipeline(Cmd([A5Grid.PYTHON_CMD[], GEN, args..., "--bbox", BBOX, "--out", out]);
                 stdout = devnull))
    return with_logger(SimpleLogger(stderr, Logging.Warn)) do
        load_mesh_geoparquet(out)
    end
end

const FAMILIES = [
    # name                 generator args                                   valence  orthogonal
    ("hex",              ("hex", "--dx", "30"),                                  6, true),
    ("rotcart30",        ("rotcart", "--dx", "40", "--angle", "30"),             4, true),
    ("voronoi_cvt",      ("voronoi", "--dx", "40", "--jitter", "0.3", "--lloyd", "5"), 0, true),
    ("voronoi_centroid", ("voronoi", "--dx", "40", "--jitter", "0.45", "--centre", "centroid"), 0, false),
]
const MESHES = Dict(name => _family(name, args...) for (name, args, _, _) in FAMILIES)

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))
_state(mesh; kw...) = _quiet() do
    st = initialise_flow_model(mesh, StandardFlow(); kw...)
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    st
end

# ---------------------------------------------------------------------------
@testset "G — geometry ($name)" for (name, _, valence, orthogonal) in FAMILIES
    mesh = MESHES[name]
    @test mesh.id_scheme == "polygon"
    st = _state(mesh)
    e  = st.edges
    @test all(isfinite, e.width) && all(>(0.0), e.width)
    @test all(isfinite, e.L)     && all(>(0.0), e.L)
    valence > 0 && @test st.max_neighbours == valence
    valence == 0 && @test st.max_neighbours >= 6
    # Closure: Σ_faces w n̂ (outward) = 0 on every interior cell
    n = length(st.cell_ids)
    sx = zeros(n); sy = zeros(n); sw = zeros(n)
    for k in 1:e.n_edges
        i, j = e.cell_i[k], e.cell_j[k]
        sx[i] += e.width[k] * e.nf_x[k]; sy[i] += e.width[k] * e.nf_y[k]
        sx[j] -= e.width[k] * e.nf_x[k]; sy[j] -= e.width[k] * e.nf_y[k]
        sw[i] += e.width[k]; sw[j] += e.width[k]
    end
    # Boundary cells (some polygon side not shared) carry ghost edges; every
    # max-valence domain-edge cell must be among them on a mixed-valence mesh.
    has_ghost = falses(n)
    for ge in st.ghost_edges; has_ghost[ge.cell_index] = true; end
    @test has_ghost == st.boundary_mask
    interior = findall(!, has_ghost)
    @test !isempty(interior)
    @test maximum(hypot(sx[i], sy[i]) / sw[i] for i in interior) < 1e-6
    if orthogonal
        @test minimum(e.cos_theta) > 0.99999
    else
        @test minimum(e.cos_theta) < 0.99      # genuinely non-orthogonal
    end
end

@testset "L — lake at rest ($name, $label)" for (name, _, _, _) in FAMILIES,
        (label, kw) in [("legacy", (;)),
                        ("WLSQ", (; gradient_correction = true)),
                        ("diamond", (; gradient_correction = true, face_flux_method = :diamond)),
                        ("cell momentum", (; momentum_model = :cell))]
    mesh = MESHES[name]
    n = length(mesh.cells)
    mesh.static_vars["elevation"] = [0.6 + 0.3 * sin(0.9 * i) + 0.2 * isodd(i) for i in 1:n]
    st = _state(mesh; kw...)
    for i in 1:n
        st.water_depth[i] = 2.0 - st.elevation[i]
        st.volume[i]      = st.water_depth[i] * st.cell_area[i]
    end
    _quiet() do
        for _ in 1:100; step_standard!(st, 1.0); end
    end
    wse = [st.elevation[i] + st.volume[i] / st.cell_area[i] for i in 1:n]
    @test maximum(abs, wse .- 2.0) < 1e-10
    @test maximum(abs, st.edges.flux) < 1e-10
    delete!(mesh.static_vars, "elevation")
end

@testset "M — mass conservation ($name)" for (name, _, _, _) in FAMILIES
    mesh = MESHES[name]
    n  = length(mesh.cells)
    st = _state(mesh; gradient_correction = true, face_flux_method = :diamond)
    src = argmin([hypot(c.center_lon - 0.003, c.center_lat - 0.003) for c in mesh.cells])
    added = 0.0
    _quiet() do
        for _ in 1:200
            st.volume[src] += 2.0; added += 2.0
            for i in 1:n; st.water_depth[i] = st.volume[i] / st.cell_area[i]; end
            step_standard!(st, 0.5)
        end
    end
    @test minimum(st.volume) >= 0.0
    @test isapprox(sum(st.volume), added; rtol = 1e-12)
    @test count(>(1e-3), st.water_depth) >= 4            # water spread beyond the source
end

@testset "I — import round trip (hex)" begin
    src  = MESHES["hex"]
    gj   = joinpath(DIR, "hex_polygons.geojson")
    feats = [Dict("type" => "Feature", "properties" => Dict(),
                  "geometry" => Dict("type" => "Polygon", "coordinates" => [c.boundary]))
             for c in src.cells]
    open(gj, "w") do io
        JSON3.write(io, Dict("type" => "FeatureCollection", "features" => feats))
    end
    out = joinpath(DIR, "hex_imported.parquet")
    run(pipeline(Cmd([A5Grid.PYTHON_CMD[], GEN, "import", "--input", gj, "--out", out]);
                 stdout = devnull))
    imp = _quiet(() -> load_mesh_geoparquet(out))
    @test length(imp.cells) == length(src.cells)
    # same polygons in the same order → same adjacency (IDs are row indices)
    @test all(sort(imp.adjacency[c.id]) == sort(src.adjacency[c.id]) for c in src.cells)
    @test maximum(abs(a.center_lon - b.center_lon) + abs(a.center_lat - b.center_lat)
                  for (a, b) in zip(imp.cells, src.cells)) < 1e-9
end

@testset "R — raster-aligned projected Cartesian grid (EPSG:27700)" begin
    # 6 × 4 pixels of 25 m near Carlisle in British National Grid; elevation
    # = 10·row_from_top + col so orientation can be checked.
    nx, ny, px = 6, 4, 25.0
    x0, y0 = 340_000.0, 556_000.0                    # north-west corner
    tif = joinpath(DIR, "bng.tif")
    vals = Float32[10 * i + j for j in 0:(nx-1), i in 0:(ny-1)]   # (nx, ny), i from top
    ArchGDAL.create(tif; driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                    nbands = 1, dtype = Float32) do ds
        ArchGDAL.setgeotransform!(ds, [x0, px, 0.0, y0, 0.0, -px])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(27700)))
        ArchGDAL.write!(ds, vals, 1)
    end
    mesh = _quiet(() -> generate_cartesian_mesh_raster(tif))
    @test length(mesh.cells) == nx * ny
    @test mesh.id_scheme == "cartesian"
    # Centres back in BNG are the pixel centres; row 0 of the mesh is the
    # southernmost raster row.
    lons = [c.center_lon for c in mesh.cells]; lats = [c.center_lat for c in mesh.cells]
    zs = zeros(length(lons))
    ArchGDAL.createcoordtrans(A5Grid._crs_gis(4326), A5Grid._crs_gis(27700)) do xf
        ArchGDAL.transform!(lons, lats, zs, xf)
    end
    # Absolute positions to 1 cm: PROJ's OSGB36↔WGS84 forward and inverse
    # transforms differ by ~1 mm, a constant offset. Relative positions (the
    # pixel lattice) to 1 µm.
    for rw in 0:(ny-1), j in 0:(nx-1)
        k = rw * nx + j + 1; i = ny - 1 - rw
        @test abs(lons[k] - (x0 + (j + 0.5) * px)) < 1e-2
        @test abs(lats[k] - (y0 - (i + 0.5) * px)) < 1e-2
        @test abs((lons[k] - lons[1]) - j * px) < 1e-6
        @test abs((lats[k] - lats[1]) - rw * px) < 1e-6
        @test mesh.static_vars["elevation"][k] == 10 * i + j
    end
    st = _state(mesh)
    @test st.max_neighbours == 4
    # FloodA5 measures areas on its R = 6371 km sphere (equirectangular per
    # cell), which near Carlisle is 0.38% below the BNG pixel area; uniform.
    @test all(a -> isapprox(a, px^2; rtol = 1e-2), st.cell_area)
    @test maximum(st.cell_area) / minimum(st.cell_area) < 1 + 1e-5
    @test length(st.ghost_edges) == 2 * (nx + ny)
    for i in 1:length(mesh.cells)
        st.water_depth[i] = 100.0 - st.elevation[i]
        st.volume[i]      = st.water_depth[i] * st.cell_area[i]
    end
    _quiet() do; for _ in 1:50; step_standard!(st, 1.0); end; end
    @test maximum(abs, st.edges.flux) < 1e-10
end

println("\nMesh-family Tier 0 tests complete.")
