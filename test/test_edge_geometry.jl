"""
test_edge_geometry.jl
---------------------
Tests for Phase 1 physics corrections and the Option D EdgeList refactor:
  - A5Grid._edge_geometry (cos θ + skewness; supersedes _edge_cos_theta)
  - _build_edge_list (EdgeList construction)
  - step_standard! / step_sgs! (EdgeList-based flux loop)

Run with:
    julia test_edge_geometry.jl

No external mesh or Python bridge required.  All geometry is constructed
synthetically from first principles so the tests run offline and fast.

Test structure
--------------
  Group 1 — _edge_geometry (cos θ component): unit geometry
      1.1  Perfectly orthogonal pair      → cos θ = 1.0 exactly
      1.2  45° skew                       → cos θ = cos(45°) ≈ 0.7071
      1.3  Arbitrary known angle          → cos θ within tolerance
      1.4  Non-adjacent cells             → 1.0 fallback (no NaN suppression)
      1.5  Near-antipodal cells           → no crash, finite result
      1.6  Symmetry: swap i and j         → same cos θ
      1.7  Resolution invariance          → same cos θ at different scales

  Group 2 — _build_edge_list: EdgeList population
      2.1  Edge count correct (3 cells → 2 edges)
      2.2  All EdgeList geometry fields finite and valid
      2.3  Canonical ordering: cell_i < cell_j for all edges
      2.4  centre_dist > 0 for all edges
      2.5  Mixed-resolution pair → 1 finite edge (MR forward-compat.)

  Group 3 — _bates_flux integration and sign convention
      3.1  cos θ = 1 reproduces old behaviour exactly
      3.2  cos θ < 1 increases effective L and reduces flux magnitude
      3.3  No flow when h_flow ≤ 0 (unchanged)
      3.4  cos_theta floor: pathological skew (cos θ = 0.05) does not blow up
      3.5  Sign: Q < 0 when wse_i > wse_j (flow i→j, i loses volume)
      3.6  Sign: Q > 0 when wse_j > wse_i (flow j→i, i gains volume)

  Group 4 — EdgeList flux, sign convention, volume update correctness
      4.1  Single edge: deeper cell loses, shallower gains (EdgeList)
      4.2  NaN cos_theta edge skipped; valid edge still processed
      4.3  step_standard!: deeper cell loses volume
      4.4  step_standard!: global volume conserved
      4.5  step_sgs!: deeper cell loses volume
      4.6  step_sgs!: global volume conserved
"""

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "mesh"))
include(joinpath(dirname(@__DIR__), "mesh", "A5Grid.jl"))

# FloodModel.jl is a flat script (no module wrapper).  We include only the
# parts needed for testing — the physics functions and types — without pulling
# in the heavy VisualisationServer / MakieVisualiser / HDF5 dependencies that
# FloodModel.jl normally loads.  We replicate the minimal preamble here.
using .A5Grid
using .A5Grid: SGSTable, wse_from_volume, wetted_area_from_wse,
               sgs_table, build_sgs_tables!,
               grid_disk_neighbours, grid_disk_neighbours_batch,
               _polygon_area_m2, _shared_edge, _haversine_m
using JSON3
using Dates
using Statistics: mean
using HDF5
using Printf

# Stub out the vis modules so FloodModel.jl can be included without GLMakie
module VisualisationServer
    start(;kwargs...)  = nothing
    stop(_)            = nothing
    push_frame!(args...) = nothing
    set_mesh!(args...)   = nothing
    notify_complete!(_)  = nothing
    struct VisServer end
end

module MakieVisualiser
    start(args...; kwargs...) = nothing
    stop(_)                   = nothing
    push_frame!(args...)      = nothing
end

include(joinpath(dirname(@__DIR__), "FloodModel.jl"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

const PASS = "\e[32m✓\e[0m"
const FAIL = "\e[31m✗\e[0m"

n_pass = Ref(0)
n_fail = Ref(0)

function check(name::String, cond::Bool, detail::String = "")
    if cond
        println("  $PASS  $name")
        n_pass[] += 1
    else
        println("  $FAIL  $name$(isempty(detail) ? "" : "  ($detail)")")
        n_fail[] += 1
    end
end

# Build a synthetic A5Cell from raw geometry (no Python bridge needed).
# The boundary is a regular pentagon centred at (lon, lat) with a given
# half-width in degrees, rotated by `rot_deg`.
function _synthetic_pentagon(lon::Float64, lat::Float64,
                              half_width_deg::Float64,
                              rot_deg::Float64 = 0.0)::A5Cell
    verts = Vector{Vector{Float64}}()
    for k in 0:4
        θ = deg2rad(rot_deg + 90.0 + 72.0 * k)
        push!(verts, [lon + half_width_deg * cos(θ),
                      lat + half_width_deg * sin(θ)])
    end
    push!(verts, verts[1])   # close the ring
    # Use a dummy id; resolution 14 ≈ half_width_deg 0.01°
    A5Cell("test_$(round(Int, lon*100))_$(round(Int, lat*100))",
           14, lon, lat, verts, NaN)
end

# Build two adjacent pentagons sharing one edge, with a controlled skew
# angle between the centre-to-centre vector and the shared edge normal.
#
# Strategy: place cell i at origin; place cell j offset so that the
# centre-to-centre vector makes angle `skew_deg` with the y-axis (north).
# Then construct a shared edge perpendicular to north (i.e. horizontal),
# so that the face normal is north and the c-to-c vector is at `skew_deg`.
# Expected cos θ = cos(skew_deg).
function _make_pair_with_skew(skew_deg::Float64;
                               hw::Float64 = 0.01)
    lon0, lat0 = 172.0, -43.5   # Christchurch-ish

    # Cell i centred at origin
    cell_i = _synthetic_pentagon(lon0, lat0, hw)

    # Cell j offset along direction skew_deg from north
    θ = deg2rad(skew_deg)   # angle from north (y-axis)
    # Separation: 2 × hw so boundaries just touch along a shared edge
    sep = 2.0 * hw
    dlon = sep * sin(θ)
    dlat = sep * cos(θ)
    cell_j = _synthetic_pentagon(lon0 + dlon, lat0 + dlat, hw)

    # Inject a synthetic shared edge: a horizontal segment at lat0 + hw
    # (i.e. the top edge of cell i and the bottom edge of cell j).
    # We override bnd_i and bnd_j to guarantee these two cells share
    # exactly two vertices at that latitude.
    shared_lat = lat0 + hw
    v1 = [lon0 - hw * 0.6, shared_lat]
    v2 = [lon0 + hw * 0.6, shared_lat]

    bnd_i = [[lon0 - hw, lat0 - hw],
              [lon0 + hw, lat0 - hw],
              [lon0 + hw, lat0 + hw],
              v2, v1,
              [lon0 - hw, lat0 - hw]]   # closed ring

    bnd_j = [v1, v2,
              [lon0 + hw * 2 + dlon, lat0 + hw],
              [lon0 + dlon, lat0 + dlat + hw],
              [lon0 - hw * 2 + dlon, lat0 + hw],
              v1]   # closed ring

    cell_i2 = A5Cell(cell_i.id, 14, lon0,         lat0,         bnd_i, NaN)
    cell_j2 = A5Cell(cell_j.id, 14, lon0 + dlon,  lat0 + dlat,  bnd_j, NaN)
    return cell_i2, cell_j2
end

# ---------------------------------------------------------------------------
# Group 1 — _edge_geometry (cos θ component) unit geometry
# ---------------------------------------------------------------------------

println("\nGroup 1 — _edge_geometry (cos θ component) unit geometry")
println("─" ^ 50)

# 1.1  Perfectly orthogonal pair: c-to-c is due north, shared edge is
#      horizontal → face normal is due north → cos θ = 1.0
let
    hw   = 0.01
    lon0 = 172.0; lat0 = -43.5
    bnd_i = [[-hw+lon0, -hw+lat0], [hw+lon0, -hw+lat0],
              [hw+lon0,  hw+lat0], [-hw+lon0, hw+lat0],
              [-hw+lon0, -hw+lat0]]
    bnd_j = [[-hw+lon0,  hw+lat0], [hw+lon0, hw+lat0],
              [hw+lon0,  3hw+lat0], [-hw+lon0, 3hw+lat0],
              [-hw+lon0,  hw+lat0]]
    ct, _, _ = A5Grid._edge_geometry(bnd_i, bnd_j, lon0, lat0-0.005,
                                              lon0, lat0+hw+0.005)
    check("1.1  orthogonal pair → cos θ ≈ 1.0",
          isfinite(ct) && abs(ct - 1.0) < 0.01,
          @sprintf("got %.4f", ct))
end

# 1.2  45° skew — constructed precisely in local projection space.
# Shared edge is horizontal at lat0.  Centre i is at lat0-hw (directly south);
# centre j is at lat0+hw (directly north of edge) but offset eastward so that
# the c-to-c vector is at exactly 45° from north.
# The total dy in metres = 2*hw * m_per_deg_lat (i south of edge, j north of edge).
# For 45°: dx_m == dy_m → dlon_deg = 2*hw * m_per_deg_lat / m_per_deg_lon.
let
    lon0 = 172.0; lat0 = -43.5; hw = 0.01
    # 1° of latitude and longitude in metres at lat0
    m_per_deg_lat = A5Grid._haversine_m(lon0, lat0, lon0, lat0 + 1.0)
    m_per_deg_lon = A5Grid._haversine_m(lon0, lat0, lon0 + 1.0, lat0)
    # dlon needed so that dx_m == dy_m  →  cos θ = cos(45°)
    dlon_deg = 2.0 * hw * m_per_deg_lat / m_per_deg_lon

    # Shared edge: horizontal segment at lat0
    v1 = [lon0 - hw, lat0]
    v2 = [lon0 + hw, lat0]
    bnd_i = [[lon0 - hw, lat0 - 2hw], [lon0 + hw, lat0 - 2hw],
              [lon0 + hw, lat0],        v1,
              [lon0 - hw, lat0 - 2hw]]
    bnd_j = [v1, v2,
              [lon0 + hw + dlon_deg, lat0 + 2hw],
              [lon0 - hw + dlon_deg, lat0 + 2hw],
              v1]

    lon_i = lon0;             lat_i = lat0 - hw   # centre i: south of edge
    lon_j = lon0 + dlon_deg;  lat_j = lat0 + hw   # centre j: north-east of edge
    ct, _, _ = A5Grid._edge_geometry(bnd_i, bnd_j, lon_i, lat_i, lon_j, lat_j)
    expected = cosd(45.0)
    check("1.2  45° skew → cos θ ≈ 0.7071",
          isfinite(ct) && abs(ct - expected) < 0.02,
          @sprintf("got %.4f, expected %.4f", ct, expected))
end

# 1.3  Arbitrary angle: 20°
let
    cell_i, cell_j = _make_pair_with_skew(20.0)
    ct, _, _ = A5Grid._edge_geometry(cell_i.boundary, cell_j.boundary,
                                 cell_i.center_lon, cell_i.center_lat,
                                 cell_j.center_lon, cell_j.center_lat)
    expected = cosd(20.0)
    check("1.3  20° skew → cos θ ≈ 0.9397",
          isfinite(ct) && abs(ct - expected) < 0.05,
          @sprintf("got %.4f, expected %.4f", ct, expected))
end

# 1.4  Non-adjacent cells → 1.0 fallback (orthogonal assumption, not NaN)
#      NaN would propagate into the step-function NaN guard and suppress flux,
#      so non-adjacent or degenerate pairs fall back to cos θ = 1.0.
let
    hw   = 0.01
    lon0 = 172.0; lat0 = -43.5
    bnd_i = [[-hw+lon0, -hw+lat0], [hw+lon0, -hw+lat0],
              [hw+lon0,  hw+lat0], [-hw+lon0, hw+lat0],
              [-hw+lon0, -hw+lat0]]
    # bnd_j far away, no shared vertices
    bnd_j = [[-hw+lon0+1.0, -hw+lat0+1.0], [hw+lon0+1.0, -hw+lat0+1.0],
              [hw+lon0+1.0,  hw+lat0+1.0], [-hw+lon0+1.0, hw+lat0+1.0],
              [-hw+lon0+1.0, -hw+lat0+1.0]]
    ct, _, _ = A5Grid._edge_geometry(bnd_i, bnd_j, lon0, lat0, lon0+1.0, lat0+1.0)
    check("1.4  non-adjacent cells → 1.0 fallback (no NaN suppression)",
          ct == 1.0, @sprintf("got %s", string(ct)))
end

# 1.5  Near-antipodal centres (edge near equator, centres far apart) — no crash
let
    hw   = 0.01
    lon0 = 0.0; lat0 = 0.0
    bnd_i = [[-hw, -hw], [hw, -hw], [hw, hw], [-hw, hw], [-hw, -hw]]
    bnd_j = [[-hw, hw],  [hw, hw],  [hw, 3hw], [-hw, 3hw], [-hw, hw]]
    ct, _, _ = A5Grid._edge_geometry(bnd_i, bnd_j, lon0, lat0-10.0, lon0, lat0+10.0)
    check("1.5  large centre separation → finite result, no crash",
          isfinite(ct))   # must be finite (NaN fallback → 1.0)
end

# 1.6  Symmetry: swap i and j → same cos θ
let
    cell_i, cell_j = _make_pair_with_skew(30.0)
    ct_ij, _, _ = A5Grid._edge_geometry(cell_i.boundary, cell_j.boundary,
                                    cell_i.center_lon, cell_i.center_lat,
                                    cell_j.center_lon, cell_j.center_lat)
    ct_ji, _, _ = A5Grid._edge_geometry(cell_j.boundary, cell_i.boundary,
                                    cell_j.center_lon, cell_j.center_lat,
                                    cell_i.center_lon, cell_i.center_lat)
    check("1.6  symmetry: swap i↔j → same cos θ",
          isfinite(ct_ij) && isfinite(ct_ji) && abs(ct_ij - ct_ji) < 1e-10,
          @sprintf("ij=%.6f  ji=%.6f", ct_ij, ct_ji))
end

# 1.7  Resolution invariance: scaling the geometry should not change cos θ
let
    hw_coarse = 0.1    # ~res 10 scale
    hw_fine   = 0.001  # ~res 17 scale
    lon0 = 172.0; lat0 = -43.5

    function orthogonal_pair(hw)
        bnd_i = [[-hw+lon0, -hw+lat0], [hw+lon0, -hw+lat0],
                  [hw+lon0,  hw+lat0], [-hw+lon0,  hw+lat0],
                  [-hw+lon0, -hw+lat0]]
        bnd_j = [[-hw+lon0,  hw+lat0], [hw+lon0,  hw+lat0],
                  [hw+lon0,  3hw+lat0], [-hw+lon0, 3hw+lat0],
                  [-hw+lon0,  hw+lat0]]
        A5Grid._edge_geometry(bnd_i, bnd_j,
                                lon0, lat0-hw/2,
                                lon0, lat0+hw*1.5)[1]   # cos_theta component
    end

    ct_coarse = orthogonal_pair(hw_coarse)
    ct_fine   = orthogonal_pair(hw_fine)
    check("1.7  resolution invariance: coarse ≈ fine cos θ",
          isfinite(ct_coarse) && isfinite(ct_fine) &&
          abs(ct_coarse - ct_fine) < 0.02,
          @sprintf("coarse=%.4f  fine=%.4f", ct_coarse, ct_fine))
end

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Group 2 — _build_edge_list and EdgeList population
# ---------------------------------------------------------------------------

println("\nGroup 2 — _build_edge_list and EdgeList population")
println("─" ^ 50)

# Build a minimal 3-cell synthetic mesh: cell 0 (centre) adjacent to cell 1
# (north) and cell 2 (east), giving 2 undirected edges total.
#
# NOTE: cell IDs must be valid hex strings, because _build_edge_list
# normalises every id via _norm(id) = _to_hex(parse(UInt64, id, base=16))
# before doing any adjacency lookup (matching the real pya5 cell-ID
# convention used in production).  Using non-hex placeholder strings like
# "c0"/"mr_coarse" here previously caused every adjacency lookup to miss
# silently (ids[i] became a zero-padded 16-char hex string but `adj`/`id_idx`
# still used the original short keys), producing an empty EdgeList with no
# error — this was a pre-existing bug in the test fixture, not in
# _build_edge_list or _edge_geometry.  Fixed by using hex IDs throughout and
# normalising the dictionary keys exactly as _build_edge_list will.
function _three_cell_mesh()
    hw   = 0.01
    lon0 = 172.0; lat0 = -43.5

    bnd0 = [[-hw+lon0, -hw+lat0], [hw+lon0, -hw+lat0],
             [hw+lon0,  hw+lat0], [0.0+lon0, 1.5hw+lat0],
             [-hw+lon0,  hw+lat0], [-hw+lon0, -hw+lat0]]

    bnd1 = [[-hw+lon0,  hw+lat0], [hw+lon0,  hw+lat0],
             [hw+lon0, 3hw+lat0], [0.0+lon0, 3.5hw+lat0],
             [-hw+lon0, 3hw+lat0], [-hw+lon0, hw+lat0]]

    bnd2 = [[hw+lon0,  -hw+lat0], [3hw+lon0, -hw+lat0],
             [3hw+lon0,  hw+lat0], [hw+lon0,  hw+lat0],
             [hw+lon0, -hw+lat0]]

    # Valid hex IDs (raw, pre-normalisation form — _build_edge_list will
    # pad these to 16 chars internally via _norm).
    id0, id1, id2 = "0", "1", "2"

    c0 = A5Cell(id0, 14, lon0,       lat0,       bnd0, NaN)
    c1 = A5Cell(id1, 14, lon0,       lat0+2hw,   bnd1, NaN)
    c2 = A5Cell(id2, 14, lon0+2hw,   lat0,       bnd2, NaN)

    cells  = [c0, c1, c2]

    # id_idx and adj must be keyed by the *normalised* hex string, since
    # that is what _build_edge_list looks them up by (ids[i] = _norm(c.id)).
    norm(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
    n0, n1, n2 = norm(id0), norm(id1), norm(id2)

    id_idx = Dict(n0 => 1, n1 => 2, n2 => 3)
    adj    = Dict(n0 => [n1, n2], n1 => [n0], n2 => [n0])
    areas  = [_polygon_area_m2(c.boundary) for c in cells]
    elevs  = [0.0, 0.0, 0.0]
    return cells, id_idx, adj, areas, elevs
end

cells, id_idx, adj, areas, elevs = _three_cell_mesh()
edges_test = _build_edge_list(cells, id_idx, adj, areas, nothing, elevs)

# 2.1  Edge count is correct (3 cells, 2 shared edges)
check("2.1  edge count correct (3 cells → 2 edges)",
      edges_test.n_edges == 2,
      "got $(edges_test.n_edges)")

# 2.2  All edge geometry fields are finite
check("2.2  all edge widths finite",
      all(isfinite, edges_test.width),
      "NaN count=$(count(isnan, edges_test.width))")
check("2.2  all edge distances finite",
      all(isfinite, edges_test.L),
      "NaN count=$(count(isnan, edges_test.L))")
check("2.2  all cos_theta finite and in [0,1]",
      all(x -> isfinite(x) && 0.0 ≤ x ≤ 1.0, edges_test.cos_theta),
      "values=$(edges_test.cos_theta)")
check("2.2  all sill values finite",
      all(isfinite, edges_test.sill),
      "NaN count=$(count(isnan, edges_test.sill))")

# 2.3  Canonical ordering: cell_i < cell_j for every edge
check("2.3  canonical ordering: cell_i < cell_j for all edges",
      all(edges_test.cell_i[e] < edges_test.cell_j[e]
          for e in 1:edges_test.n_edges),
      "violations: $([(edges_test.cell_i[e], edges_test.cell_j[e])
                       for e in 1:edges_test.n_edges
                       if edges_test.cell_i[e] >= edges_test.cell_j[e]])")

# 2.4  centre_dist > 0 for all edges
check("2.4  centre_dist > 0 for all edges",
      all(x -> x > 0.0, edges_test.L),
      "min=$(minimum(edges_test.L))")

# 2.5  Mixed-resolution pair: _build_edge_list works across resolution levels
#      (forward-compatible with Phase 3 MR mesh).
let
    hw_coarse = 0.01; hw_fine = 0.002
    lon0 = 172.0; lat0 = -43.5

    bnd_c = [[-hw_coarse+lon0, -hw_coarse+lat0],
              [ hw_coarse+lon0, -hw_coarse+lat0],
              [ hw_coarse+lon0,  hw_coarse+lat0],
              [-hw_coarse+lon0,  hw_coarse+lat0],
              [-hw_coarse+lon0, -hw_coarse+lat0]]
    bnd_f = [[-hw_fine+lon0,  hw_coarse+lat0],
              [ hw_fine+lon0,  hw_coarse+lat0],
              [ hw_fine+lon0,  hw_coarse+hw_fine*2+lat0],
              [-hw_fine+lon0,  hw_coarse+hw_fine*2+lat0],
              [-hw_fine+lon0,  hw_coarse+lat0]]

    # Valid hex IDs — see note in _three_cell_mesh() above. The original
    # "mr_coarse"/"mr_fine" placeholder strings were not valid hex at all
    # (contain 'r', 's', etc.), so parse(UInt64, ..., base=16) would throw
    # outright. Use plain hex IDs and key id_idx/adj by their normalised form.
    id_c, id_f = "10", "11"
    c_c = A5Cell(id_c, 14, lon0, lat0,                  bnd_c, NaN)
    c_f = A5Cell(id_f, 15, lon0, lat0+hw_coarse+hw_fine, bnd_f, NaN)

    cells_mr = [c_c, c_f]
    norm(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
    n_c, n_f = norm(id_c), norm(id_f)

    id_idx_mr = Dict(n_c => 1, n_f => 2)
    adj_mr    = Dict(n_c => [n_f], n_f => [n_c])
    areas_mr  = [_polygon_area_m2(c.boundary) for c in cells_mr]
    elevs_mr  = [0.0, 0.0]

    el_mr = _build_edge_list(cells_mr, id_idx_mr, adj_mr, areas_mr,
                              nothing, elevs_mr)
    check("2.5  mixed-resolution pair → 1 finite edge (MR forward-compat.)",
          el_mr.n_edges == 1 &&
          isfinite(el_mr.cos_theta[1]) && 0.0 ≤ el_mr.cos_theta[1] ≤ 1.0,
          @sprintf("n_edges=%d  cos_theta=%.4f", el_mr.n_edges,
                   isnan(el_mr.cos_theta[1]) ? -1.0 : el_mr.cos_theta[1]))
end

# ---------------------------------------------------------------------------
# Group 3 — _bates_flux integration
# ---------------------------------------------------------------------------

println("\nGroup 3 — _bates_flux integration")
println("─" ^ 50)

# Pull _bates_flux from FloodModel module directly
const _G = 9.81

# 3.1  cos θ = 1 reproduces old behaviour exactly
#      (compare to the pre-correction formula with L_eff = L)
let
    q_prev = 0.0; wse_i = 2.0; wse_j = 1.0
    z_sill = 0.0; width = 100.0; L = 1500.0
    n_mann = 0.03; dt = 10.0
    cos1   = 1.0

    q_new = _bates_flux(q_prev, wse_i, wse_j, z_sill,
                                    width, L, cos1, n_mann, dt)

    # Manual calculation with L_eff = L (old formula)
    h_flow = max(wse_i, wse_j) - z_sill
    dWSE   = wse_i - wse_j
    num    = q_prev - _G * h_flow * dt * dWSE / L
    den    = 1.0 + _G * h_flow * dt * n_mann^2 * abs(q_prev) / h_flow^(10/3)
    q_ref  = (num / den) * width

    check("3.1  cos θ = 1.0 reproduces pre-correction formula",
          abs(q_new - q_ref) < 1e-12,
          @sprintf("new=%.6f  ref=%.6f", q_new, q_ref))
end

# 3.2  cos θ < 1 shortens L_eff and increases the slope, increasing flux magnitude
let
    q_prev = 0.0; wse_i = 2.0; wse_j = 1.0
    z_sill = 0.0; width = 100.0; L = 1500.0
    n_mann = 0.03; dt = 10.0

    Q_ortho = _bates_flux(q_prev, wse_i, wse_j, z_sill,
                                      width, L, 1.0, n_mann, dt)
    Q_skew  = _bates_flux(q_prev, wse_i, wse_j, z_sill,
                                      width, L, 0.7, n_mann, dt)

    check("3.2  cos θ < 1 increases |flux| relative to orthogonal",
          abs(Q_skew) > abs(Q_ortho),
          @sprintf("|Q_ortho|=%.4f  |Q_skew|=%.4f", abs(Q_ortho), abs(Q_skew)))
end

# 3.3  No flow when h_flow ≤ 0
let
    q_prev = 1.0; wse_i = 0.5; wse_j = 0.3
    z_sill = 1.0   # sill above both WSEs → h_flow = max(0.5,0.3) - 1.0 < 0
    width  = 100.0; L = 1500.0; n_mann = 0.03; dt = 10.0

    Q = _bates_flux(q_prev, wse_i, wse_j, z_sill,
                                width, L, 0.8, n_mann, dt)
    check("3.3  no flow when h_flow ≤ 0",
          Q == 0.0, @sprintf("Q = %.4f", Q))
end

# 3.4  Pathological skew (cos θ = 0.05) does not produce Inf/NaN
let
    q_prev = 0.5; wse_i = 3.0; wse_j = 1.0
    z_sill = 0.0; width = 100.0; L = 1500.0
    n_mann = 0.03; dt = 10.0

    Q = _bates_flux(q_prev, wse_i, wse_j, z_sill,
                                width, L, 0.05, n_mann, dt)
    check("3.4  pathological skew (cos θ = 0.05) → finite flux",
          isfinite(Q), @sprintf("Q = %s", string(Q)))
end

# 3.5  Sign convention: Q < 0 when first WSE argument > second.
#      In EdgeList context: pass (wse_ci, wse_cj); Q < 0 means ci is higher,
#      water flows toward cj.  dV[ci] += Q*dt < 0 → ci loses volume. Correct.
let
    q_prev = 0.0; wse_i = 3.0; wse_j = 1.0   # i is higher → flow i→j
    z_sill = 0.0; width = 100.0; L = 1500.0
    n_mann = 0.03; dt = 10.0

    Q = _bates_flux(q_prev, wse_i, wse_j, z_sill, width, L, 1.0, n_mann, dt)
    check("3.5  sign: Q < 0 when wse_i > wse_j (flow goes i→j)",
          Q < 0.0, @sprintf("Q = %.4f (expected < 0)", Q))
end

# 3.6  Reverse: Q > 0 when second WSE arg > first → ci gains volume.
let
    q_prev = 0.0; wse_i = 1.0; wse_j = 3.0   # j is higher → flow j→i
    z_sill = 0.0; width = 100.0; L = 1500.0
    n_mann = 0.03; dt = 10.0

    Q = _bates_flux(q_prev, wse_i, wse_j, z_sill, width, L, 1.0, n_mann, dt)
    check("3.6  sign: Q > 0 when wse_j > wse_i (flow goes j→i, i gains)",
          Q > 0.0, @sprintf("Q = %.4f (expected > 0)", Q))
end

# ---------------------------------------------------------------------------
# Group 4 — EdgeList flux, sign convention, volume update correctness
# ---------------------------------------------------------------------------
# Tests the Option D EdgeList refactor and all bugs fixed along the way:
#   Bug 24: j == 0 && break → continue (sparse adj_matrix slots)
#   Bug 28/31: double-counted / order-dependent loops → EdgeList solves both
#   Bug 29: inverted continuity sign → dV[ci] += Q*dt; dV[cj] -= Q*dt
#   Bug 30: per-edge volume limiter → cell-level post-hoc cap
#
# EdgeList sign convention (cell_i < cell_j always):
#   Q > 0 → flow from cell_j to cell_i  (j is higher, water toward lower i)
#   Q < 0 → flow from cell_i to cell_j  (i is higher, water toward lower j)
#   dV[ci] += Q*dt   (ci gains when Q > 0)
#   dV[cj] -= Q*dt   (cj loses when Q > 0)
# ---------------------------------------------------------------------------

println("\nGroup 4 — EdgeList flux, sign convention, volume update")
println("─" ^ 50)

# Helpers: build a minimal EdgeList and FlowState without a mesh or Python bridge.

function _make_edge_list(edges_spec::Vector{Tuple{Int,Int,Float64}};
                          width::Float64=200.0, L::Float64=1500.0,
                          cos_theta::Float64=1.0)
    ne  = length(edges_spec)
    ci  = [e[1] for e in edges_spec]
    cj  = [e[2] for e in edges_spec]
    sls = [e[3] for e in edges_spec]
    EdgeList(ne, ci, cj,
             fill(width, ne), fill(L, ne), fill(cos_theta, ne),
             sls,
             zeros(Float64, ne),   # flux
             zeros(Float64, ne),   # flux_Q
             zeros(Int, ne),       # collinear_i
             zeros(Int, ne),       # collinear_j
             zeros(Float64, ne),   # skew_x (orthogonal synthetic edges)
             zeros(Float64, ne),   # skew_y
             fill(L, ne),          # dx_m = L (edges run east)
             zeros(Float64, ne),   # dy_m = 0
             fill(1.0, ne),        # nf_x = 1 (east, orthogonal)
             zeros(Float64, ne))   # nf_y = 0
end

function _minimal_state(; n_cells::Int,
                           elev::Vector{Float64}=zeros(n_cells),
                           depth::Vector{Float64}=zeros(n_cells),
                           volume::Vector{Float64}=zeros(n_cells),
                           edges::EdgeList,
                           manning_n::Vector{Float64}=fill(0.03, n_cells),
                           cell_area::Vector{Float64}=fill(1.5e6, n_cells),
                           sgs_tables::Vector{Any}=Any[],
                           q_centre_theta::Float64=0.9,
                           gradient_correction::Bool=false,
                           gradient_correction_alpha::Float64=1.0,
                           momentum_model::Symbol=:edge,
                           face_flux_method::Symbol=:legacy)
    adj_matrix = zeros(Int, 5, n_cells)   # not used by step functions in these tests

    cell_edge_index = zeros(Int, N_SIDES, n_cells)
    _build_cell_edge_index!(cell_edge_index, adj_matrix, edges, n_cells)
    mom_weights = zeros(Float64, 10, n_cells)
    _build_mom_weights!(mom_weights, cell_edge_index, edges, n_cells)

    FlowState(
        ["c$i" for i in 1:n_cells],
        copy(depth), copy(volume),
        zeros(n_cells), zeros(n_cells), zeros(n_cells),   # velocity, vel_u, vel_v
        copy(elev),
        copy(manning_n), copy(cell_area),
        zeros(n_cells), zeros(n_cells),   # cell_lons, cell_lats (unused by these tests)
        Dict{String,Vector{String}}(),
        adj_matrix,
        edges,
        sgs_tables,
        falses(n_cells),         # boundary_mask — no open boundaries in these synthetic tests
        Any[],                   # ghost_edges — empty: Phase D is a no-op
        Any[],                   # ghost_cell_bc
        0.0,                     # vol_removed
        zeros(Float64, 2, n_cells),    # grad_wse — unused while gradient_correction=false
        zeros(Float64, 10, n_cells),   # wlsq_weights — unused while gradient_correction=false
        q_centre_theta,
        gradient_correction,
        gradient_correction_alpha,
        zeros(Float64, n_cells), zeros(Float64, n_cells),   # qvec_u, qvec_v
        cell_edge_index, mom_weights,
        momentum_model,
        face_flux_method, nothing,   # face_flux_method, diamond_table (Phase C)
    )
end

# 4.1  Single edge between cell 1 (deeper) and cell 2 (shallower).
#      With EdgeList, cell_i=1, cell_j=2 (1 < 2).
#      wse_ci > wse_cj → Q < 0 → dV[1] += Q*dt < 0 (1 loses) → depth[1] falls.
let
    edges = _make_edge_list([(1, 2, 0.0)])   # one edge, sill=0
    depth  = [2.0, 0.5]
    volume = depth .* 1.5e6
    state  = _minimal_state(n_cells=2, depth=depth, volume=volume, edges=edges)

    step_standard!(state, 30.0)

    check("4.1  single edge: deeper cell loses, shallower cell gains",
          state.water_depth[1] < depth[1] && state.water_depth[2] > depth[2],
          @sprintf("depth[1]=%.6f depth[2]=%.6f", state.water_depth[1], state.water_depth[2]))
end

# 4.2  NaN cos_theta edge is skipped; a second valid edge still processes.
let
    # Two edges between cells 1 and 2: first has NaN cos_theta, second is valid.
    # In practice only one edge per pair exists, but we test the NaN guard.
    ne = 2
    edges = EdgeList(ne, [1,1], [2,2],
                     fill(200.0, ne), fill(1500.0, ne),
                     [NaN, 1.0],     # slot 0: NaN → skipped; slot 1: valid
                     fill(0.0, ne),
                     zeros(Float64, ne),    # flux
                     zeros(Float64, ne),    # flux_Q
                     zeros(Int, ne),        # collinear_i
                     zeros(Int, ne),        # collinear_j
                     zeros(Float64, ne),    # skew_x
                     zeros(Float64, ne),    # skew_y
                     fill(1500.0, ne),      # dx_m = L (east)
                     zeros(Float64, ne),    # dy_m = 0
                     fill(1.0, ne),         # nf_x = 1 (east, orthogonal)
                     zeros(Float64, ne))    # nf_y = 0
    depth  = [2.0, 0.5]
    volume = depth .* 1.5e6
    state  = _minimal_state(n_cells=2, depth=depth, volume=volume, edges=edges)

    d2_before = state.water_depth[2]
    step_standard!(state, 30.0)

    check("4.2  NaN cos_theta edge skipped; valid edge still processed",
          state.water_depth[2] > d2_before,
          @sprintf("depth[2] before=%.6f after=%.6f", d2_before, state.water_depth[2]))
end

# 4.3  step_standard!: deeper cell loses volume (sign convention correct).
let
    edges = _make_edge_list([(1, 2, 0.0)])
    depth  = [2.0, 0.5];  volume = depth .* 1.5e6
    state  = _minimal_state(n_cells=2, depth=depth, volume=volume, edges=edges)
    v0 = copy(state.volume)
    step_standard!(state, 30.0)

    check("4.3  step_standard!: deeper cell loses volume",
          state.volume[1] < v0[1],
          @sprintf("vol[1] before=%.1f after=%.1f", v0[1], state.volume[1]))
end

# 4.4  step_standard!: global volume conserved (no limiter triggered).
#      Three cells in a chain: 1-2, 2-3.  With cell_i < cell_j enforced.
let
    edges = _make_edge_list([(1,2,0.0), (2,3,0.0)])
    depth  = [2.0, 1.0, 0.5];  volume = [3.0e6, 1.5e6, 0.75e6]
    state  = _minimal_state(n_cells=3, depth=depth, volume=volume, edges=edges)

    total_before = sum(state.volume)
    step_standard!(state, 10.0)
    total_after  = sum(state.volume)

    rel_err = abs(total_after - total_before) / max(total_before, 1.0)
    check("4.4  step_standard!: global volume conserved",
          rel_err < 1e-10,
          @sprintf("rel_err=%.2e", rel_err))
end

# 4.5  step_sgs!: deeper cell loses volume.
let
    area  = 1.5e6
    n_bins = 5
    elev_bins  = collect(range(0.0, 4.0, length=n_bins))
    vol_curve  = area .* elev_bins
    area_curve = fill(area, n_bins)

    # step_sgs! always uses the R-A kernel (_manning_flux_ra) — there is no
    # Bates fallback in the current implementation despite the SGSTable
    # docstring's "falls back to Bates if zero" note (that note describes
    # intended legacy-mesh behaviour, not what step_sgs! actually does).
    # zeros(n_bins, 5) for edge_area_curves would make _manning_flux_ra's
    # `A <= 1e-6 && return 0.0` dry-edge guard fire on every edge, silently
    # producing zero flux — test 4.5's assertion would then fail, correctly
    # flagging that no water moved, but for the wrong reason (a degenerate
    # fixture, not the physics under test).  Build a simple rectangular
    # channel cross-section instead: width 50m, so at WSE = elev_bins[k]
    # the flow area through any edge is 50 × elev_bins[k] (m²) and the
    # wetted perimeter is 50 + 2×elev_bins[k] (m) — a trapezoidal-ish
    # rectangular trench, consistent with the T-EA1 analytical test case
    # described in FloodA5_SGS_RA_Flux_Implementation.md §4.1.
    chan_width = 50.0
    edge_area_curve  = repeat(chan_width .* elev_bins, 1, 5)
    edge_perim_curve = repeat(chan_width .+ 2.0 .* elev_bins, 1, 5)

    tbl = SGSTable(elev_bins, vol_curve, area_curve, area, 0.0, 4.0,
                   edge_area_curve, edge_perim_curve)

    edges  = _make_edge_list([(1, 2, 0.0)])
    volume = [3.0e6, 0.75e6]
    state  = _minimal_state(n_cells=2, volume=volume, edges=edges,
                             sgs_tables=Any[tbl, tbl])
    v0 = copy(state.volume)
    step_sgs!(state, 30.0)

    check("4.5  step_sgs!: deeper cell loses volume",
          state.volume[1] < v0[1],
          @sprintf("vol[1] before=%.1f after=%.1f", v0[1], state.volume[1]))
end

# 4.6  step_sgs!: global volume conserved.
let
    area  = 1.5e6
    n_bins = 5
    elev_bins  = collect(range(0.0, 4.0, length=n_bins))
    vol_curve  = area .* elev_bins
    area_curve = fill(area, n_bins)

    # See note in test 4.5 above: step_sgs! always uses the R-A kernel, so
    # the edge curves must be non-degenerate for this test to exercise real
    # flux rather than a silent zero-flow no-op.
    chan_width = 50.0
    edge_area_curve  = repeat(chan_width .* elev_bins, 1, 5)
    edge_perim_curve = repeat(chan_width .+ 2.0 .* elev_bins, 1, 5)

    tbl = SGSTable(elev_bins, vol_curve, area_curve, area, 0.0, 4.0,
                   edge_area_curve, edge_perim_curve)

    edges  = _make_edge_list([(1,2,0.0), (2,3,0.0)])
    volume = [3.0e6, 1.5e6, 0.75e6]
    state  = _minimal_state(n_cells=3, volume=volume, edges=edges,
                             sgs_tables=Any[tbl, tbl, tbl])

    total_before = sum(state.volume)
    step_sgs!(state, 10.0)
    total_after  = sum(state.volume)

    rel_err = abs(total_after - total_before) / max(total_before, 1.0)
    check("4.6  step_sgs!: global volume conserved",
          rel_err < 1e-10,
          @sprintf("rel_err=%.2e", rel_err))
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

println()
println("─" ^ 50)
total = n_pass[] + n_fail[]
@printf("Results: %d / %d passed", n_pass[], total)
if n_fail[] == 0
    println("  — all tests passed ✓")
else
    println()
    println("  $(n_fail[]) test(s) FAILED")
    exit(1)
end
