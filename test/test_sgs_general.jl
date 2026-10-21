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
    @test MESH.static_vars["sgs_layout"][1] == 2.0
    @test size(MESH.array_vars["sgs_edge_sills"], 1) == 6
    @test size(MESH.array_vars["sgs_edge_area_curve"], 1) == 40 * 6
end

@testset "S2 — slot k is the k-th neighbour" begin
    idx   = Dict(c.id => i for (i, c) in enumerate(MESH.cells))
    sills = MESH.array_vars["sgs_edge_sills"]
    tbl(i) = A5Grid.sgs_table(MESH, i)
    n_pairs = 0; n_sym = 0; n_perim_ok = 0; n_perim = 0
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
            P = tbl(i).edge_perim_curves[end, k]
            n_perim += 1
            P <= W * (1 + 1e-6) && (n_perim_ok += 1)
        end
    end
    @test n_pairs > 100
    @test n_sym == n_pairs              # every edge's sill agrees from both sides
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

@testset "S4 — pre-fix SGS layout is refused" begin
    old = deepcopy(MESH)
    delete!(old.static_vars, "sgs_layout")
    @test_throws ErrorException _quiet(() -> initialise_flow_model(old, SGSFlow()))
    st = withenv("FLOODA5_ALLOW_OLD_SGS_LAYOUT" => "1") do
        _quiet(() -> initialise_flow_model(old, SGSFlow()))
    end
    @test length(st.sgs_tables) == length(old.cells)
end

println("\nGeneral-mesh SGS tests complete.")
