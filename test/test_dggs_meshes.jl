# test/test_dggs_meshes.jl
# ------------------------
# Tier 0 on discrete-global-grid meshes (WP3 part 3, 2026-11 —
# docs/notes/MeshFamilies_DGGS_Plan_2026-11.md). Meshes come from
# mesh/dggs_families.py: H3, IGEO7 (ISEA7H_Z7), ISEA3H, HEALPix, rHEALPix,
# over a ~650 m box near Carlisle at a 40 m target cell size. A family whose
# library (h3, dggal) is not installed is skipped with a message.
#
#   P   generator: level choice (--list-levels, --target-dx), vertex snapping
#       across bucket boundaries, reproducible output
#   G   geometry: id scheme and level recorded, interior cells have the
#       family's valence, every interior cell closes (Σ w·n̂ = 0), faces
#       within 40° of orthogonal, equal-area families uniform to 1e-5
#   L   lake at rest under the legacy, WLSQ, diamond and cell-momentum options
#   M   mass conservation: closed point-source run, exact volume balance
#   X   HEALPix corners against healpix-geo (WGS84), an independent library
#       (skipped if healpix-geo is not installed)
#   C   command line (0.8.1): --meshgen-dggs with --dx and --meshres; the
#       mesh metadata survives a re-save and reaches the output's /mesh
#       attributes; --rerun reproduces the run; a5 through the same flag;
#       a mesh without neighbour lists is refused by build_sgs_tables!
#       unless it is A5; a file without an ID scheme is loaded with a warning
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_dggs_meshes.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging, Statistics

const GEN  = joinpath(@__DIR__, "..", "mesh", "dggs_families.py")
const BBOX = "-2.905,54.897,-2.895,54.903"     # ≈ 640 × 670 m near Carlisle
const CTR  = (-2.900, 54.900)
const DIR  = mktempdir()

A5Grid._find_python()
_py(code) = success(pipeline(Cmd([A5Grid.PYTHON_CMD[], "-c", code]); stdout = devnull, stderr = devnull))
_pyout(args...) = read(Cmd([A5Grid.PYTHON_CMD[], args...]), String)
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

const HAVE = Dict(
    "h3"    => _py("import h3"),
    "dggal" => _py("import sys; sys.path.insert(0, r'$(dirname(GEN))'); import dggs_families as d; d._dggal()"),
)

const FAMILIES = [
    # name        library  valence  equal-area
    ("h3",        "h3",    6,       false),
    ("igeo7",     "dggal", 6,       true),
    ("isea3h",    "dggal", 6,       true),
    ("healpix",   "dggal", 4,       true),
    ("rhealpix",  "dggal", 4,       true),
]
const ACTIVE = [f for f in FAMILIES if HAVE[f[2]]]
# CI installs the libraries and sets FLOODA5_REQUIRE_DGGS=1: there a library
# that does not import (or DGGAL's setup crashing) is a failure, not a skip.
if get(ENV, "FLOODA5_REQUIRE_DGGS", "") == "1"
    @testset "libraries present" begin
        @test HAVE["h3"]
        @test HAVE["dggal"]
    end
end
for (name, lib, _, _) in FAMILIES
    HAVE[lib] || @warn "test_dggs_meshes: $lib not importable — $name skipped"
end

function _generate(name, args...; out = joinpath(DIR, "$name.parquet"))
    s = _pyout(GEN, name, "--bbox=$BBOX", args..., "--out", out)
    js = JSON3.read(last(split(strip(s), '\n')))
    haskey(js, :error) && error("dggs_families $name: $(js[:error])")
    return js
end

const SUMMARY = Dict(name => _generate(name, "--target-dx", "40") for (name, _, _, _) in ACTIVE)
const MESHES  = Dict(name => _quiet(() -> load_mesh_geoparquet(joinpath(DIR, "$name.parquet")))
                     for (name, _, _, _) in ACTIVE)

_state(mesh; kw...) = _quiet() do
    st = initialise_flow_model(mesh, StandardFlow(); kw...)
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    st
end

# ---------------------------------------------------------------------------
@testset "P — generator" begin
    @testset "level choice ($name)" for (name, _, _, _) in ACTIVE
        js = JSON3.read(_pyout(GEN, name, "--bbox=$BBOX", "--list-levels"))
        lv = js[:levels]
        @test issorted([r[:level] for r in lv])
        @test issorted([r[:cell_area_m2] for r in lv]; rev = true)
        # --target-dx picked the level closest to 40² m² in log terms
        best = argmin(r -> abs(log(r[:cell_area_m2] / 40.0^2)), lv)
        @test SUMMARY[name][:level] == best[:level]
        @test SUMMARY[name][:adjacency_checked] == true
    end
    @testset "vertex snapping" begin
        code = """
import sys; sys.path.insert(0, r'$(dirname(GEN))')
import dggs_families as d
s = d.VertexSnapper(1e-9)
a = s((0.9999999999999e-9 * 5, 1.0)); b = s((1.0000000000001e-9 * 5, 1.0))   # straddle a bucket edge
c = s((5e-9 + 3e-9, 1.0))                                                   # a distinct vertex
assert a is b and a != c and s.merged == 1, (a, b, c, s.merged)
"""
        @test _py(code)
    end
    if !isempty(ACTIVE)
        @testset "reproducible ($(ACTIVE[1][1]))" begin
            name = ACTIVE[1][1]
            _generate(name, "--target-dx", "40"; out = joinpath(DIR, "again.parquet"))
            m2 = _quiet(() -> load_mesh_geoparquet(joinpath(DIR, "again.parquet")))
            m1 = MESHES[name]
            @test [c.id for c in m1.cells] == [c.id for c in m2.cells]
            @test all(a.boundary == b.boundary for (a, b) in zip(m1.cells, m2.cells))
        end
    end
end

@testset "G — geometry ($name)" for (name, _, valence, equal_area) in ACTIVE
    mesh = MESHES[name]
    @test mesh.id_scheme == name
    @test all(c -> c.resolution == SUMMARY[name][:level], mesh.cells)
    st = _state(mesh)
    e  = st.edges
    @test st.max_neighbours == valence
    @test all(isfinite, e.width) && all(>(0.0), e.width)
    @test all(isfinite, e.L)     && all(>(0.0), e.L)
    n = length(st.cell_ids)
    sx = zeros(n); sy = zeros(n); sw = zeros(n); deg = zeros(Int, n)
    for k in 1:e.n_edges
        i, j = e.cell_i[k], e.cell_j[k]
        sx[i] += e.width[k] * e.nf_x[k]; sy[i] += e.width[k] * e.nf_y[k]
        sx[j] -= e.width[k] * e.nf_x[k]; sy[j] -= e.width[k] * e.nf_y[k]
        sw[i] += e.width[k]; sw[j] += e.width[k]
        deg[i] += 1; deg[j] += 1
    end
    has_ghost = falses(n)
    for ge in st.ghost_edges; has_ghost[ge.cell_index] = true; end
    @test has_ghost == st.boundary_mask
    interior = findall(!, has_ghost)
    @test length(interior) > n ÷ 2
    @test all(deg[i] == valence for i in interior)         # no pentagons in this box
    # Each face is projected about its own mid-latitude (A5Grid._edge_geometry),
    # so away from the equator a cell closes to O(L/R · tan φ), L the cell size
    # (≈ 1e-5 for 40 m cells at 55° N; the equatorial families of
    # test_mesh_families.jl close to 1e-6).
    tol = 2 * sqrt(mean(st.cell_area)) / A5Grid._EARTH_R * tand(CTR[2])
    @test maximum(hypot(sx[i], sy[i]) / sw[i] for i in interior) < tol
    @test minimum(e.cos_theta) > cosd(40.0)
    equal_area && @test std(st.cell_area) / mean(st.cell_area) < 1e-5
end

@testset "L — lake at rest ($name, $label)" for (name, _, _, _) in ACTIVE,
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

@testset "M — mass conservation ($name)" for (name, _, _, _) in ACTIVE
    mesh = MESHES[name]
    n  = length(mesh.cells)
    st = _state(mesh; gradient_correction = true, face_flux_method = :diamond)
    src = argmin([hypot(c.center_lon - CTR[1], c.center_lat - CTR[2]) for c in mesh.cells])
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
    @test count(>(1e-3), st.water_depth) >= 4
end

# HEALPix from DGGAL against healpix-geo on the WGS84 authalic sphere: the
# same cells (by centre) must have the same corners.
if haskey(MESHES, "healpix") && _py("import healpix_geo")
    @testset "X — HEALPix corners against healpix-geo" begin
        cells = MESHES["healpix"].cells
        lev = SUMMARY["healpix"][:level]
        pts = joinpath(DIR, "hp_centres.txt")
        open(pts, "w") do io
            for c in cells; println(io, c.center_lon, " ", c.center_lat); end
        end
        code = """
import numpy as np, healpix_geo.nested as N
lon, lat = np.loadtxt(r'$pts', ndmin=2).T
ip = N.lonlat_to_healpix(lon, lat, $lev, ellipsoid='WGS84')
vlo, vla = N.vertices(ip, $lev, ellipsoid='WGS84')
vlo = np.where(vlo > 180, vlo - 360, vlo)
for a, b in zip(vlo, vla): print(' '.join(f'{float(x)!r} {float(y)!r}' for x, y in zip(a, b)))
"""
        rows = split(strip(_pyout("-c", code)), '\n')
        @test length(rows) == length(cells)
        worst = 0.0
        for (c, r) in zip(cells, rows)
            v = parse.(Float64, split(r))
            ref = [(v[2k-1], v[2k]) for k in 1:4]
            mine = [(p[1], p[2]) for p in c.boundary[1:end-1]]
            for p in ref
                worst = max(worst, minimum(abs(p[1] - q[1]) + abs(p[2] - q[2]) for q in mine))
            end
        end
        @test worst < 1e-9        # degrees (≈ 0.1 mm)
    end
else
    @info "test_dggs_meshes: healpix-geo (or the HEALPix mesh) not available — X skipped"
end

# ---------------------------------------------------------------------------
# C — the command line (0.8.1)
const AOI = joinpath(DIR, "aoi.geojson")
let (x0, y0, x1, y1) = parse.(Float64, split(BBOX, ','))
    write(AOI, """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
      "geometry":{"type":"Polygon","coordinates":[[[$x0,$y0],[$x1,$y0],[$x1,$y1],[$x0,$y1],[$x0,$y0]]]}}]}""")
end
_main(args) = with_logger(NullLogger()) do
    redirect_stdout(devnull) do
        main(copy(args))
    end
end
_same_results(a, b) = FO.frame_times(a) == FO.frame_times(b) &&
    all(v -> isequal(FO.read_var(a, v), FO.read_var(b, v)), FO.frame_vars(a))

@testset "C — command line" begin
    @test a5_level_for_dx(22.2) == 18 && a5_level_for_dx(44.4) == 17
    @test_throws ErrorException dggs_generate("s2", AOI; target_dx = 40.0)
    @test_throws ErrorException dggs_generate("h3", AOI)                  # neither level nor dx
    @test_throws ErrorException dggs_generate("h3", AOI; level = 12, target_dx = 40.0)
    @test "--meshgen-dggs" in FO._DROP_FLAGS && "--dggs" in FO._DROP_FLAGS

    if HAVE["h3"]
        mout = joinpath(DIR, "cli_h3.parquet"); h5 = joinpath(DIR, "cli_h3.h5")
        args = ["--meshgen-dggs", AOI, "--dggs", "h3", "--dx", "40", "--meshout", mout,
                "--flow-model", "standard", "--rainfall", "40", "--sim-duration", "300",
                "--output", h5, "--output-interval", "100"]
        _main(args)
        m = _quiet(() -> load_mesh_geoparquet(mout))
        @test m.id_scheme == "h3"
        @test m.meta["mesh_family"] == "h3" && m.meta["dggs_level"] == SUMMARY["h3"][:level]
        @test m.meta["target_dx_m"] == 40.0
        # the same cells as the generator run directly
        @test [c.id for c in m.cells] == [c.id for c in MESHES["h3"].cells]
        # the metadata reaches the output, and the run reruns from it
        mi = FO.mesh_info(h5)
        @test mi["id_scheme"] == "h3" && mi["mesh_family"] == "h3"
        @test mi["dggs_level"] == SUMMARY["h3"][:level] && mi["dggs_equal_area"] == 0
        @test mi["dggs_refinement_ratio"] == 7
        @test FO.run_info(h5)["argv"] == args
        h5b = joinpath(DIR, "cli_h3_rerun.h5")
        _main(["--rerun", h5, "--output", h5b])
        @test _same_results(h5, h5b)
        ri = FO.run_info(h5b)
        @test "--meshload" in ri["argv"] && !("--meshgen-dggs" in ri["argv"]) && !("--dggs" in ri["argv"])
        @test FO.mesh_info(h5b)["mesh_family"] == "h3"
        # a re-save (as after DEM sampling) keeps the metadata
        again = joinpath(DIR, "cli_h3_resaved.parquet")
        _quiet(() -> save_mesh_geoparquet(m, again))
        m2 = _quiet(() -> load_mesh_geoparquet(again))
        @test m2.meta == m.meta && m2.id_scheme == "h3"
        # no neighbour lists: SGS cannot fall back to pya5's grid_disk
        m3 = deepcopy(m); empty!(m3.adjacency)
        err = try A5Grid.build_sgs_tables!(m3, FileDEM("nonexistent.tif")); nothing catch e; e end
        @test err isa ErrorException && occursin("only A5 meshes", err.msg)
        # a file without an ID scheme: the old guess (a5-v0), with a warning —
        # here wrong, so with pya5 1.x the A5 scheme check refuses it as well
        bare = joinpath(DIR, "cli_h3_bare.parquet")
        @test _py("import pyarrow.parquet as pq; t = pq.read_table(r'$mout'); " *
                  "md = {k: v for k, v in (t.schema.metadata or {}).items() if k != b'flooda5'}; " *
                  "pq.write_table(t.replace_schema_metadata(md), r'$bare')")
        mb = withenv("FLOODA5_ALLOW_ID_SCHEME_MISMATCH" => "1") do
            @test_logs (:warn, r"records no cell-ID scheme") match_mode = :any load_mesh_geoparquet(bare)
        end
        @test mb.id_scheme == "a5-v0" && isempty(mb.meta)
    end
    if HAVE["dggal"]
        mout = joinpath(DIR, "cli_healpix.parquet")
        lev = SUMMARY["healpix"][:level]
        _main(["--meshgen-dggs", AOI, "--dggs", "healpix", "--meshres", string(lev),
               "--meshout", mout, "--mesh-only"])
        m = _quiet(() -> load_mesh_geoparquet(mout))
        @test m.id_scheme == "healpix" && m.meta["dggs_level"] == lev
        @test m.meta["dggs_equal_area"] == true && !haskey(m.meta, "target_dx_m")
        @test length(m.cells) == length(MESHES["healpix"].cells)
    end
end

println("\nDGGS mesh Tier 0 tests complete.")
