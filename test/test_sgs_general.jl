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
#   S4  tables older than layout 4 (0.8.0) are refused
#   S5  boundary sides have face tables (0.7.0)
#   S6  the open-boundary slope is the valley floor's (0.7.0)
#   S7  critical-flow limit (0.8.0): water released through the notch in the
#       ridge; no face passes more than its critical discharge at the
#       upstream level, which the local-inertial flux alone exceeds
#   S8  SGS on discrete-global-grid meshes (0.8.1): H3, IGEO7 and HEALPix
#       (mesh/dggs_families.py) over the same DEM — S1–S3's layout, slot
#       order, lake at rest and mass balance, and an open-boundary drain
#       (a family whose library is not installed is skipped, except in CI)
#   S9  the command line builds SGS tables for a loaded mesh without them
#       when --dem is given (0.8.1), and refuses without --dem
#   S10 sources with the flux update (--source-order post, the SGS default):
#       at steady state a point source's cell holds the same level whatever
#       Δt is; added before the fluxes, its level moves by ≈ Q·ΔΔt/area —
#       the coupling that let the adaptive Δt drive the SGS checkerboard
#       (HYDRAULICS §15.15). Mass balance either way; the defaults resolve
#       (sgs: post, standard: pre)
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
    @test MESH.static_vars["sgs_layout"][1] == 4.0          # face tables (0.7.0) + critical integral (0.8.0)
    @test size(MESH.array_vars["sgs_edge_sills"], 1) == 6
    @test size(MESH.array_vars["sgs_face_zr"], 1) == 6
    for k in ("sgs_face_A", "sgs_face_P", "sgs_face_K", "sgs_face_C")
        @test size(MESH.array_vars[k], 1) == A5Grid.SGS_FACE_LEVELS * 6
    end
    @test !haskey(MESH.array_vars, "sgs_edge_area_curve")   # the per-cell edge curves are gone
end

"""
Slot-order properties of a mesh's SGS tables: (pairs, sills agreeing from both
sides, face tables identical from both sides, perimeters within their edge, faces).
"""
function _slot_props(M)
    idx   = Dict(c.id => i for (i, c) in enumerate(M.cells))
    sills = M.array_vars["sgs_edge_sills"]
    tbl(i) = A5Grid.sgs_table(M, i)
    n_pairs = 0; n_sym = 0; n_perim_ok = 0; n_perim = 0; n_same = 0
    for (i, c) in enumerate(M.cells)
        nbrs = M.adjacency[c.id]
        for (k, nb) in enumerate(nbrs)
            j = idx[nb]
            kj = findfirst(==(c.id), M.adjacency[nb])
            kj === nothing && continue
            n_pairs += 1
            isapprox(sills[k, i], sills[kj, j]; atol = 1e-6) && (n_sym += 1)
            edge = A5Grid._shared_edge(c.boundary, M.cells[j].boundary)
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
    return (; n_pairs, n_sym, n_same, n_perim_ok, n_perim)
end

@testset "S2 — slot k is the k-th neighbour" begin
    p = _slot_props(MESH)
    @test p.n_pairs > 100
    @test p.n_sym == p.n_pairs              # every edge's sill agrees from both sides
    @test p.n_same == p.n_pairs             # … and its face table is the same table (0.7.0)
    @test p.n_perim_ok == p.n_perim         # no curve longer than its own edge
end

function _sgs_state(M = MESH)
    st = _quiet(() -> initialise_flow_model(M, SGSFlow()))
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    return st
end

"Lake at rest and point-source mass balance under SGS (S3, S8)."
function _sgs_rest_and_mass(M)
    st = _sgs_state(M)
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
    st2 = _sgs_state(M)
    src = div(n, 2); added = 0.0
    _quiet() do
        for _ in 1:200
            st2.volume[src] += 5.0; added += 5.0
            step_sgs!(st2, 1.0)
        end
    end
    @test isapprox(sum(st2.volume), added; rtol = 1e-10)
    @test minimum(st2.volume) >= 0.0
    @test count(>(0.0), st2.volume) >= 4
    return st
end

@testset "S3 — SGS runs on the 6-neighbour mesh" begin
    st = _sgs_rest_and_mass(MESH)
    @test st.max_neighbours == 6
end

@testset "S4 — SGS tables older than layout 4 are refused" begin
    old = deepcopy(MESH)
    delete!(old.static_vars, "sgs_layout")                   # pre-2026-10: wrong slot order
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old, SGSFlow()))
    old2 = deepcopy(MESH)
    old2.static_vars["sgs_layout"] .= 2.0                    # 0.2–0.6: per-cell edge curves
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old2, SGSFlow()))
    old3 = deepcopy(MESH)
    delete!(old3.array_vars, "sgs_face_A")
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old3, SGSFlow()))
    old4 = deepcopy(MESH)                                    # 0.7.x: half a pixel off, no face_C
    old4.static_vars["sgs_layout"] .= 3.0; delete!(old4.array_vars, "sgs_face_C")
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old4, SGSFlow()))
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

@testset "S7 — critical-flow limit (0.8.0)" begin
    "Volume of a cell at level η (its tables; the full cell area above them)."
    function _vol_at(t, η)
        b = t.elev_bins; v = t.vol_curve
        η <= b[1] && return 0.0
        η >= b[end] && return v[end] + t.cell_area * (η - b[end])
        k = searchsortedlast(b, η)
        return v[k] + (η - b[k]) / (b[k+1] - b[k]) * (v[k+1] - v[k])
    end
    # A pond west of the ridge (x ≈ 180 m, 1 m high, a 40 m gap at y ≈ 200 m)
    # at 1.2 m, the rest dry: the water falls through the gap. Returns the
    # largest |Q| ÷ critical discharge over the faces, each step, with the
    # upstream level of that step (an upper bound of the kernel's cap, which
    # also scales it by the flow-path factor ≤ 1).
    function run_pond(crit::Bool)
        old = SGS_CRITICAL[]; SGS_CRITICAL[] = crit
        try
            st = _sgs_state()
            for i in eachindex(st.volume)
                t = st.sgs_tables[i]
                st.volume[i] = isnan(t.z_min) || st.cell_lons[i] > rad2deg(170.0 / A5Grid._EARTH_R) ? 0.0 :
                               _vol_at(t, 1.2)
            end
            E = st.edges; worst = 0.0; nbind = 0
            for _ in 1:40
                η = [st.volume[i] > 0 ? wse_from_volume(st.sgs_tables[i], st.volume[i]) : NaN for i in eachindex(st.volume)]
                _quiet(() -> step_sgs!(st, 0.5))
                for e in 1:E.n_edges
                    Q = E.flux_Q[e]; Q == 0.0 && continue
                    ci, cj = E.cell_i[e], E.cell_j[e]
                    up, dn = Q > 0 ? (cj, ci) : (ci, cj)
                    isnan(η[up]) && continue
                    slot = findfirst(==(cj), view(st.adj_matrix, :, ci))
                    Qc = _sgs_critical_q(A5Grid.face_crit(st.sgs_tables[ci], slot, η[up]))
                    Qc > 0 || continue
                    r = abs(Q) / Qc
                    worst = max(worst, r); r > 0.999 && (nbind += 1)
                end
            end
            return worst, nbind, sum(st.volume) + st.vol_removed
        finally
            SGS_CRITICAL[] = old
        end
    end
    w_on, nb_on, _ = run_pond(true)
    w_off, _, _ = run_pond(false)
    @test w_on <= 1.0 + 1e-9               # never above critical
    @test nb_on > 0                        # … and the limit binds in this release
    @test w_off > 1.2                      # the local-inertial flux alone passes more
end

# ---------------------------------------------------------------------------
# S8 — discrete global grids (0.8.1). A box inside the DEM, ~30 m cells.
const DGGS_GEN = joinpath(@__DIR__, "..", "mesh", "dggs_families.py")
_py(code) = success(pipeline(Cmd([A5Grid.PYTHON_CMD[], "-c", code]); stdout = devnull, stderr = devnull))
const DGGS_HAVE = Dict(
    "h3"    => _py("import h3"),
    "dggal" => _py("import sys; sys.path.insert(0, r'$(dirname(DGGS_GEN))'); import dggs_families as d; d._dggal()"))
const DGGS_SGS = [("h3", "h3", 6), ("igeo7", "dggal", 6), ("healpix", "dggal", 4)]
if get(ENV, "FLOODA5_REQUIRE_DGGS", "") == "1"
    @testset "S8 — DGGS libraries present" begin
        @test DGGS_HAVE["h3"] && DGGS_HAVE["dggal"]
    end
end

@testset "S8 — SGS on DGGS meshes ($fam)" for (fam, lib, valence) in DGGS_SGS
    if !DGGS_HAVE[lib]
        @warn "test_sgs_general: $lib not importable — S8 $fam skipped"
        continue
    end
    out = joinpath(DIR, "$fam.parquet")
    run(pipeline(Cmd([A5Grid.PYTHON_CMD[], DGGS_GEN, fam, "--bbox=0.0004,0.0004,0.0026,0.0026",
                      "--target-dx", "30", "--out", out]); stdout = devnull))
    M = _quiet() do
        m = load_mesh_geoparquet(out)
        sample_dem!(m, DEM)
        build_sgs_tables!(m, DEM; n_bins = 40, n_samples = 256)
        m
    end
    @test M.id_scheme == fam
    @test M.static_vars["sgs_layout"][1] == 4.0
    @test size(M.array_vars["sgs_edge_sills"], 1) == valence
    @test size(M.array_vars["sgs_face_A"], 1) == A5Grid.SGS_FACE_LEVELS * valence
    @test all(isfinite, M.static_vars["sgs_z_min"])            # the DEM covers every cell
    p = _slot_props(M)
    @test p.n_pairs > 50
    @test p.n_sym == p.n_pairs && p.n_same == p.n_pairs && p.n_perim_ok == p.n_perim
    st = _sgs_rest_and_mass(M)
    @test st.max_neighbours == valence
    # open boundaries: a flooded mesh drains, and volume leaves only through them
    so = _quiet(() -> initialise_flow_model(M, SGSFlow()))
    for i in eachindex(so.volume)
        t = so.sgs_tables[i]
        so.volume[i] = t.vol_curve[end] + t.cell_area * 0.5
    end
    v0 = sum(so.volume)
    _quiet() do; for _ in 1:200; step_sgs!(so, 0.5); end; end
    @test sum(so.volume) < v0 && so.vol_removed > 0
    @test isapprox(sum(so.volume) + so.vol_removed, v0; rtol = 1e-10)
    @test minimum(so.volume) >= 0.0
end

@testset "S9 — --meshload … --dem --flow-model sgs builds the tables (0.8.1)" begin
    out = joinpath(DIR, "hex_sgs_cli.parquet")
    with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            main(["--meshload", HEX, "--dem", DEM, "--flow-model", "sgs", "--sgs-bins", "40",
                  "--sgs-samples", "256", "--meshout", out, "--mesh-only"])
        end
    end
    m = _quiet(() -> load_mesh_geoparquet(out))
    @test m.static_vars["sgs_layout"][1] == 4.0 && size(m.array_vars["sgs_edge_sills"], 1) == 6
    @test m.id_scheme == "polygon" && m.meta["mesh_family"] == "hex"
    # the same tables as built in-process (S1: same mesh, DEM, bins, samples)
    @test m.array_vars["sgs_vol_curve"] == MESH.array_vars["sgs_vol_curve"]
    @test m.array_vars["sgs_face_A"] == MESH.array_vars["sgs_face_A"]
    # without --dem the loaded mesh is still refused
    @test_throws ErrorException with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            main(["--meshload", HEX, "--flow-model", "sgs", "--rainfall", "10", "--sim-duration", "10"])
        end
    end
end

@testset "S10 — sources with the flux update (--source-order post)" begin
    # A point source on the slope below the ridge, draining to the open west
    # side; Δt fixed by dt_max (the CFL step is ~10 s), run to steady state.
    Q = 40.0; T = 2700.0
    function src_run(order, dtmax)
        st = _quiet(() -> initialise_flow_model(MESH, SGSFlow(); scheme_options(:recommended, :sgs)...))
        c = argmin(hypot.(st.cell_lons .- 0.0011, st.cell_lats .- 0.0015))
        src = InjectionPoint(c, st.cell_ids[c], Q, st.cell_lons[c], st.cell_lats[c])
        _quiet() do
            run_simulation!(st, MESH, T, dtmax, nothing, :none; method = SGSFlow(),
                            courant = 0.5, all_sources = AbstractSource[src], source_order = order)
        end
        t = st.sgs_tables[c]
        return (wse = t.z_min + st.water_depth[c], a = t.cell_area, c = c, st = st)
    end
    pre1 = src_run(:pre, 0.5);  pre2 = src_run(:pre, 1.0)
    post1 = src_run(:post, 0.5); post2 = src_run(:post, 1.0)
    d_pre = pre1.wse - pre2.wse; d_post = post1.wse - post2.wse
    @test isapprox(d_pre, Q * 0.5 / pre1.a; rtol = 0.1)        # pre: the level follows Δt (8.6 mm) …
    @test abs(d_post) < 0.05 * abs(d_pre)                       # … post: it does not (< 0.1 mm)
    @test post1.st.vol_removed > 0.8 * Q * T                    # (steady: nearly all of it drains)
    for r in (pre1, pre2, post1, post2)                         # mass balance
        @test isapprox(sum(r.st.volume) + r.st.vol_removed, Q * T; rtol = 1e-9)
        @test minimum(r.st.volume) >= 0.0
        @test r.st.water_depth[r.c] ≈ wse_from_volume(r.st.sgs_tables[r.c], r.st.volume[r.c]) -
                                      r.st.sgs_tables[r.c].z_min   # depth kept in step with the source
    end
    @test_throws ArgumentError run_simulation!(_quiet(() -> initialise_flow_model(MESH, SGSFlow())),
                                               MESH, 1.0, 1.0, nothing, :none; source_order = :middle)
    @test SGS_DEFAULT_SOURCE_ORDER === :post
end

println("\nGeneral-mesh SGS tests complete.")
