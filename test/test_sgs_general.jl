# test/test_sgs_general.jl
# ------------------------
# SGS tables on a general polygon mesh, and the slot-order fix (WP3 part 2,
# 2026-10 — see build_sgs_tables!). Builds SGS tables on a hexagonal mesh
# (6 neighbours; mesh/mesh_families.py) over a small synthetic GeoTIFF.
#
#   S1  layout: 6 edge slots, sgs_layout = 2
#   S2  slot order: slot k of cell i is its k-th neighbour. Checked through a
#       property the old grid_disk-order layout broke: the sill a cell stores
#       for neighbour j equals the sill j stores for it (same edge), and a
#       fully wet edge's perimeter never exceeds that edge's length
#   S3  SGS initialises on the 6-neighbour mesh, keeps a lake at rest and
#       conserves mass in a point-source run
#   S4  tables without sgs_layout (pre-fix) are refused unless
#       FLOODA5_ALLOW_OLD_SGS_LAYOUT=1
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_sgs_general.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

const DIR = mktempdir()
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

# Hex mesh, ~30 m cells over ~330 m × 330 m at the equator
A5Grid._find_python()
const HEX = joinpath(DIR, "hex.parquet")
run(pipeline(Cmd([A5Grid.PYTHON_CMD[], joinpath(@__DIR__, "..", "mesh", "mesh_families.py"),
                  "hex", "--dx", "30", "--bbox", "0.0,0.0,0.003,0.003", "--out", HEX]); stdout = devnull))

# DEM (EPSG:4326, 2 m pixels): gentle slope + a ridge with a notch, so edge
# sills vary from edge to edge.
const DEM = joinpath(DIR, "dem.tif")
let px = rad2deg(2.0 / A5Grid._EARTH_R), nx = 200, ny = 200
    z = Float32[0.002 * (2.0 * i) + (abs(2.0 * i - 180.0) < 6.0 && abs(2.0 * j - 200.0) > 20.0 ? 1.0 : 0.0) +
                0.05 * sin(0.13 * i) * cos(0.07 * j) for i in 0:(nx-1), j in 0:(ny-1)]
    ArchGDAL.create(DEM; driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                    nbands = 1, dtype = Float32) do ds
        ArchGDAL.setgeotransform!(ds, [-0.0001, px, 0.0, 0.0031, 0.0, -px])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
        ArchGDAL.write!(ds, z, 1)
    end
end

const MESH = _quiet() do
    m = load_mesh_geoparquet(HEX)
    sample_dem!(m, DEM)
    build_sgs_tables!(m, DEM; n_bins = 40, n_samples = 256)
    m
end

@testset "S1 — layout" begin
    @test MESH.static_vars["sgs_layout"][1] == 3.0          # face tables (0.7.0)
    @test size(MESH.array_vars["sgs_edge_sills"], 1) == 6
    @test size(MESH.array_vars["sgs_face_zr"], 1) == 6
    for k in ("sgs_face_A", "sgs_face_P", "sgs_face_K")
        @test size(MESH.array_vars[k], 1) == A5Grid.SGS_FACE_LEVELS * 6
    end
    @test !haskey(MESH.array_vars, "sgs_edge_area_curve")   # the per-cell edge curves are gone
end

@testset "S2 — slot k is the k-th neighbour" begin
    idx   = Dict(c.id => i for (i, c) in enumerate(MESH.cells))
    sills = MESH.array_vars["sgs_edge_sills"]
    tbl(i) = A5Grid.sgs_table(MESH, i)
    n_pairs = 0; n_sym = 0; n_perim_ok = 0; n_perim = 0; n_same = 0
    for (i, c) in enumerate(MESH.cells)
        nbrs = MESH.adjacency[c.id]
        for (k, nb) in enumerate(nbrs)
            j = idx[nb]
            kj = findfirst(==(c.id), MESH.adjacency[nb])
            kj === nothing && continue
            n_pairs += 1
            isapprox(sills[k, i], sills[kj, j]; atol = 1e-6) && (n_sym += 1)
            edge = A5Grid._shared_edge(c.boundary, MESH.cells[j].boundary)
            W = A5Grid._haversine_m(edge...)
            ti, tj = tbl(i), tbl(j)
            P = ti.face_P[end, k]
            n_perim += 1
            P <= W * (1 + 1e-6) && (n_perim_ok += 1)
            # one table per face: both cells hold the same numbers
            (ti.face_z0[k] == tj.face_z0[kj] && ti.face_zr[k] == tj.face_zr[kj] &&
             ti.face_A[:, k] == tj.face_A[:, kj] && ti.face_K[:, k] == tj.face_K[:, kj]) && (n_same += 1)
        end
    end
    @test n_pairs > 100
    @test n_sym == n_pairs              # every edge's sill agrees from both sides
    @test n_same == n_pairs             # … and its face table is the same table (0.7.0)
    @test n_perim_ok == n_perim         # no curve longer than its own edge
end

function _sgs_state()
    st = _quiet(() -> initialise_flow_model(MESH, SGSFlow()))
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    return st
end

@testset "S3 — SGS runs on the 6-neighbour mesh" begin
    st = _sgs_state()
    @test st.max_neighbours == 6
    n  = length(st.cell_ids)
    # Lake at rest: flat WSE above every cell's terrain (volume above z_max is
    # the full-cell column; check it round-trips through wse_from_volume)
    for i in 1:n
        t = st.sgs_tables[i]
        st.volume[i] = t.vol_curve[end] + t.cell_area * (3.0 - t.elev_bins[end])
        @assert isapprox(wse_from_volume(t, st.volume[i]), 3.0; atol = 1e-9)
    end
    v0 = copy(st.volume)
    _quiet() do; for _ in 1:50; step_sgs!(st, 1.0); end; end
    @test maximum(abs, st.volume .- v0) < 1e-6 * maximum(v0)
    @test maximum(abs, st.edges.flux_Q) < 1e-8
    # Mass balance, point source
    st2 = _sgs_state()
    src = div(n, 2); added = 0.0
    _quiet() do
        for _ in 1:200
            st2.volume[src] += 5.0; added += 5.0
            step_sgs!(st2, 1.0)
        end
    end
    @test isapprox(sum(st2.volume), added; rtol = 1e-10)
    @test minimum(st2.volume) >= 0.0
end

@testset "S4 — SGS tables older than layout 3 are refused" begin
    old = deepcopy(MESH)
    delete!(old.static_vars, "sgs_layout")                   # pre-2026-10: wrong slot order
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old, SGSFlow()))
    old2 = deepcopy(MESH)
    old2.static_vars["sgs_layout"] .= 2.0                    # 0.2–0.6: per-cell edge curves
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old2, SGSFlow()))
    old3 = deepcopy(MESH)
    delete!(old3.array_vars, "sgs_face_A")
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old3, SGSFlow()))
    @test length(_quiet(() -> initialise_flow_model(MESH, SGSFlow())).sgs_tables) == length(MESH.cells)
end

@testset "S5 — boundary sides have face tables; open-boundary outflow uses them (0.7.0)" begin
    @test A5Grid.boundary_side_slots(5, Bool[true, false, true]) == [2, 4, 5]
    @test A5Grid.boundary_side_slots(4, Bool[true, true, true, true]) == Int[]
    st = _quiet(() -> initialise_flow_model(MESH, SGSFlow()))
    G = st.ghost_edges
    @test !isempty(G)
    with_tab = [g for g in G if g.slot > 0]
    @test length(with_tab) > 0.8 * length(G)        # the DEM covers (nearly) every boundary side
    idx = Dict(c.id => i for (i, c) in enumerate(MESH.cells))
    for g in with_tab
        t = st.sgs_tables[g.cell_index]
        nb = MESH.adjacency[MESH.cells[g.cell_index].id]
        @test g.slot > length(nb) || !haskey(idx, nb[g.slot])   # a slot without a mesh neighbour
        @test g.sill == t.face_z0[g.slot]                        # the side's own lowest point
        @test t.face_P[end, g.slot] <= g.width * (1 + 1e-6)
        @test t.face_P[end, g.slot] > 0.5 * g.width
    end
    # area and radius of a side: the table, not the rectangle above the cell's lowest ground
    g = with_tab[argmax([st.sgs_tables[x.cell_index].face_zr[x.slot] for x in with_tab])]
    t = st.sgs_tables[g.cell_index]; η = g.sill + 0.3
    A, R = _sgs_ghost_area_radius(t, g.slot, η, g.sill, g.width, true)
    a, p, k = A5Grid.face_apk(t, g.slot, η)
    @test A == a && R == _sgs_radius(a, p, k, true)
    @test A < g.width * 0.3                                      # less than the full-width rectangle
    Ar, Rr = _sgs_ghost_area_radius(t, 0, η, g.sill, g.width, true)
    @test Ar ≈ g.width * 0.3 && Rr ≈ Ar / (g.width + 0.6)
    @test _sgs_ghost_area_radius(t, g.slot, g.sill, g.sill, g.width, true) == (0.0, 0.0)
    # open boundaries drain a flooded mesh; volume leaves only through the ghost edges
    for i in eachindex(st.volume)
        tt = st.sgs_tables[i]
        st.volume[i] = isnan(tt.z_min) ? 0.0 : tt.vol_curve[end] + tt.cell_area * 0.5
    end
    v0 = sum(st.volume)
    for _ in 1:200; step_sgs!(st, 0.5); end
    @test sum(st.volume) < v0
    @test isapprox(sum(st.volume) + st.vol_removed, v0; rtol = 1e-10)
    @test minimum(st.volume) >= 0.0
end

@testset "S6 — SGS open boundary: the normal-depth slope is the local valley floor's (0.7.0)" begin
    st = _quiet(() -> initialise_flow_model(MESH, SGSFlow(); open_boundary = :normal))
    G = st.ghost_edges
    @test all(g -> isnan(g.slope) || g.slope >= 1e-4, G)
    @test count(g -> isfinite(g.slope), G) > 0.9 * length(G)
    # the DEM rises 0.2% to the east: the ground leading to the west edge falls
    # towards it at about that slope (the 1 m ridge on the way steepens the fit);
    # on the east edge it falls away from the boundary: the floor, 1e-4
    lon = [st.cell_lons[g.cell_index] for g in G]
    west = [g.slope for (g, l) in zip(G, lon) if l < 0.0004 && isfinite(g.slope)]
    east = [g.slope for (g, l) in zip(G, lon) if l > 0.0026 && isfinite(g.slope)]
    @test length(west) >= 5 && length(east) >= 5
    @test all(s -> 0.001 <= s <= 0.008, west)
    @test all(==(1e-4), east)
    # the sides of one boundary cell share its slope; standard flow has none
    for g in G, h in G
        g.cell_index == h.cell_index && @test isequal(g.slope, h.slope)
    end
    st2 = _quiet(() -> initialise_flow_model(MESH, StandardFlow(); open_boundary = :normal))
    @test all(g -> isnan(g.slope), st2.ghost_edges)
    # the outflow follows it: a flooded mesh drains through the west edge, not the east
    for i in eachindex(st.volume)
        t = st.sgs_tables[i]
        st.volume[i] = isnan(t.z_min) ? 0.0 : t.vol_curve[end] + t.cell_area * max(0.0, 3.5 - t.z_max)
    end
    v0 = sum(st.volume)
    for _ in 1:100; step_sgs!(st, 0.5); end
    @test st.vol_removed > 0.0
    @test isapprox(sum(st.volume) + st.vol_removed, v0; rtol = 1e-10)
end

println("\nGeneral-mesh SGS tests complete.")
