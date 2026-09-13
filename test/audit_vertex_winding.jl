#!/usr/bin/env julia
# ============================================================================
# audit_vertex_winding.jl  (v2)
#
# Phase A, Tasks 1 and 2 of FloodA5_PhaseA_ImplementationScope.md.
#
# Reads a real mesh (parquet) and runs two independent, simulation-free
# geometry audits against the CURRENT mesh/geometry pipeline
# (mesh/A5Grid.jl, FloodModel.jl's _build_edge_list). No flux kernel, no
# FlowState physics, nothing here is modified from what production code
# already does — this script only *queries* it.
#
# ── v2 change note ──────────────────────────────────────────────────────
# v1 grouped cells by "leading hex nibble", per docs/A5_QUIRKS.md §1's one
# worked example (res 14, "8e…"/"9e…"). The first real run of this script
# found 51 genuine geometric (shared-vertex) edges among 27 cells that ALL
# shared the SAME leading nibble — which is impossible under
# A5_QUIRKS.md §1's own definition of sublattice ("every cell surrounded
# exclusively by the OTHER sublattice"). That proves leading-hex-nibble is
# NOT a reliable sublattice discriminator in general (it was a coincidence
# of that one example, not a stated pya5 API guarantee) — it is NOT
# recomputed as ground truth here.
#
# v2 instead determines the true two-sublattice split (if it exists) by
# 2-colouring the mesh's own adjacency GRAPH directly (BFS bipartite
# check) — this is a direct test of A5_QUIRKS.md §1's actual geometric
# claim, with no dependency on pya5 ID internals at all. Leading-nibble is
# still reported, as a secondary cross-check, but the graph colouring is
# the ground truth used throughout.
#
# ── Task 1 — Vertex-winding / sublattice audit ─────────────────────────────
# Does pya5's raw polygon boundary vertex ordering (A5Cell.boundary) have a
# winding convention (CW/CCW) that differs systematically between the two
# A5 sublattices (as determined by graph 2-colouring, see above)? And does
# _shared_edge's raw (pre-orientation-flip) vertex ordering correlate with
# sublattice membership or compass bearing, rather than being "unrelated to
# which cell is i and which cell is j" as _edge_geometry's docstring
# currently claims (asserted, not yet tested)?
#
# ── Task 2 — Control-volume-centre audit ───────────────────────────────────
# What does A5Cell.center_lon/center_lat actually represent? Is it the
# polygon centroid? Do DEM sampling and the WLSQ gradient stencil both
# anchor on the same point? (The latter is answered by direct code
# inspection, printed below — see the note in that section.)
#
# Usage:
#   julia --project=. test/audit_vertex_winding.jl <mesh.parquet> [sample_n]
#
#   e.g. julia --project=. test/audit_vertex_winding.jl \
#            test/square/square_mesh_res14.parquet
#
# sample_n (optional, default 100000 i.e. effectively "all") caps the number
# of cells used for the centroid check (Task 2), for speed on very large
# meshes. The winding and graph-colouring checks always use every cell —
# a partial sample would break the bipartiteness/adjacency analysis.
# ============================================================================

using Statistics

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_path = length(ARGS) >= 1 ? ARGS[1] :
            joinpath("test", "square", "square_mesh_res14.parquet")
sample_n  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100_000

println("=" ^ 78)
println("audit_vertex_winding.jl (v2) — Phase A Tasks 1 & 2")
println("Mesh: $mesh_path")
println("=" ^ 78)

mesh  = load_mesh_geoparquet(mesh_path)
cells = mesh.cells
n     = length(cells)
println("\nLoaded $n cells.")

_norm_id(id) = A5Grid._to_hex(parse(UInt64, id, base=16))
ids    = [_norm_id(c.id) for c in cells]
id_idx = Dict{String,Int}(ids[i] => i for i in 1:n)

# Recover adjacency exactly as initialise_flow_model does: prefer the
# mesh's own pre-computed adjacency dict (mesh.adjacency, stored in the
# parquet at mesh-build time — geometrically exact shared-vertex adjacency,
# see mesh/a5_bridge.py's _build_adjacency_shared_vertices), falling back
# to shared-vertex detection for older parquets that predate it.
adj = if !isempty(mesh.adjacency)
    Dict{String,Vector{String}}(_norm_id(k) => [_norm_id(v) for v in vs]
                                 for (k, vs) in mesh.adjacency)
else
    println("(no pre-computed adjacency in parquet — falling back to")
    println(" shared-vertex detection, this may take a little longer)")
    _build_adjacency_shared_vertices(mesh)
end

adj_by_idx = [Int[] for _ in 1:n]
for i in 1:n
    for nb in get(adj, ids[i], String[])
        j = get(id_idx, nb, 0)
        j != 0 && push!(adj_by_idx[i], j)
    end
end
n_edges_graph = sum(length.(adj_by_idx)) ÷ 2
println("Adjacency graph: $n cells, $n_edges_graph edges")

nibble_class(id) = id[1]   # leading hex nibble — docs/A5_QUIRKS.md §1 heuristic;
                           # kept ONLY as a secondary cross-check, see v2 note above.
nibbles = [nibble_class(i) for i in ids]
nibble_counts = Dict{Char,Int}()
for s in nibbles
    nibble_counts[s] = get(nibble_counts, s, 0) + 1
end
println("Leading hex nibble grouping (secondary cross-check only): $nibble_counts")

# ────────────────────────────────────────────────────────────────────────────
# GROUND-TRUTH SUBLATTICE CHECK — adjacency-graph bipartite 2-colouring
# ────────────────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 78)
println("Ground-truth sublattice check — adjacency-graph bipartite 2-colouring")
println("=" ^ 78)
println("docs/A5_QUIRKS.md §1 defines the two sublattices geometrically: \"every")
println("cell is surrounded exclusively by the OTHER sublattice\" — i.e. the")
println("adjacency graph should be BIPARTITE, with the two colour classes being")
println("the two sublattices. This is checked directly via BFS 2-colouring, with")
println("no dependency on pya5's ID encoding.")

colour = fill(0, n)   # 0 = unvisited, 1/2 = colour classes
is_bipartite = true
odd_cycle_example = nothing
for start in 1:n
    colour[start] != 0 && continue
    colour[start] = 1
    queue = [start]
    while !isempty(queue)
        u = popfirst!(queue)
        for v in adj_by_idx[u]
            if colour[v] == 0
                colour[v] = colour[u] == 1 ? 2 : 1
                push!(queue, v)
            elseif colour[v] == colour[u] && is_bipartite
                global is_bipartite = false
                global odd_cycle_example = (u, v)
            end
        end
    end
end

n_col1 = count(==(1), colour)
n_col2 = count(==(2), colour)
println("\nBipartite 2-colouring result: class 1 = $n_col1 cells, class 2 = $n_col2 cells")

if is_bipartite
    println("Graph IS bipartite — a genuine two-sublattice structure exists on this")
    println("mesh, consistent with docs/A5_QUIRKS.md §1's geometric definition.")
else
    ci, cj = odd_cycle_example
    println("⚠ Graph is NOT bipartite (found same-colour adjacent cells, e.g. " *
            "$(ids[ci])  ↔  $(ids[cj])).")
    println("  docs/A5_QUIRKS.md §1's 'two complementary sublattices, each cell")
    println("  surrounded exclusively by the other' claim does NOT hold as a strict")
    println("  global property of this mesh's adjacency graph — plausibly because")
    println("  boundary/AOI-edge cells (without a full 5-neighbour ring) break the")
    println("  clean bipartite pattern locally, or because the underlying A5 tiling")
    println("  has isolated odd-degree defect cells (cf. the pentagon 'pole' cells")
    println("  in geodesic dodecahedral tilings generally — worth checking against")
    println("  a larger, more interior-heavy mesh before concluding this is a bug).")
    println("  Sublattice-based comparisons below use the best-effort 2-colouring")
    println("  from the BFS (correct on the bipartite-respecting bulk of the mesh,")
    println("  approximate near whatever broke bipartiteness).")
end

# Cross-check the graph colouring against the leading-nibble grouping.
nibble_vals = sort(collect(keys(nibble_counts)))
if length(nibble_vals) == 1
    bipartite_desc = is_bipartite ? "a genuine bipartition" : "an approximate 2-class split"
    println("\nOnly ONE leading-nibble value present ('$(nibble_vals[1])'), yet the")
    println("graph colouring found $bipartite_desc of sizes $n_col1/$n_col2.")
    println("→ CONFIRMS leading hex nibble is NOT the sublattice discriminator at")
    println("  this resolution/AOI (a real discriminator cannot classify all cells")
    println("  into one bucket when the graph itself splits into two). The v1 script's")
    println("  premise (nibble = sublattice) was wrong to generalise from")
    println("  docs/A5_QUIRKS.md §1's single res-14 worked example. Recommend fixing")
    println("  docs/A5_QUIRKS.md §1 to either state the correct pya5-level")
    println("  discriminator (if one exists) or drop the nibble claim and point to")
    println("  this graph-colouring method instead.")
elseif length(nibble_vals) == 2
    n11 = count(k -> nibbles[k]==nibble_vals[1] && colour[k]==1, 1:n)
    n12 = count(k -> nibbles[k]==nibble_vals[1] && colour[k]==2, 1:n)
    println("\nCross-tab (nibble '$(nibble_vals[1])' vs graph colour): colour1=$n11  colour2=$n12")
    if n11 == 0 || n12 == 0
        println("→ Leading-nibble grouping AGREES EXACTLY with the graph bipartition on")
        println("  this mesh — the nibble heuristic happens to be valid here too.")
    else
        println("→ Leading-nibble grouping does NOT match the graph bipartition even")
        println("  though both have 2 classes — nibble is not a reliable discriminator.")
    end
else
    println("\n$(length(nibble_vals)) distinct leading-nibble values present — more than the")
    println("2 sublattices the graph colouring found. Nibble is clearly encoding")
    println("something else (e.g. resolution/parent structure), not sublattice,")
    println("on this mesh.")
end

# Ground-truth grouping used for the rest of this script.
sublattice_of = Dict{Int,Symbol}(i => (colour[i] == 1 ? :A : :B) for i in 1:n)

# ────────────────────────────────────────────────────────────────────────────
# TASK 1a — Signed polygon area (winding sign) by (graph-true) sublattice
# ────────────────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 78)
println("TASK 1a — Signed polygon area (winding) by sublattice")
println("=" ^ 78)
println("Convention: local equirectangular frame (east=+x, north=+y) centred on")
println("each polygon's own vertex-mean latitude, standard shoelace sum WITHOUT")
println("abs(). Positive = counter-clockwise (CCW) vertex order; negative = CW.")
println("Sublattice grouping used here is the GRAPH-DERIVED colouring above, not")
println("leading hex nibble.")

"""
    _signed_area_m2(boundary) → Float64

Same projection as A5Grid._polygon_area_m2 (equirectangular, centred on the
polygon's own vertex-mean latitude) but returns the SIGNED shoelace sum, not
abs(). Positive = CCW, negative = CW, in a standard east=+x/north=+y frame.
"""
function _signed_area_m2(boundary::Vector{Vector{Float64}})::Float64
    m = length(boundary)
    m < 3 && return 0.0
    lat0 = mean(v[2] for v in boundary)
    cos_lat = cosd(lat0)
    R = A5Grid._EARTH_R
    area = 0.0
    j = m
    for i in 1:m
        xi = deg2rad(boundary[i][1]) * R * cos_lat
        yi = deg2rad(boundary[i][2]) * R
        xj = deg2rad(boundary[j][1]) * R * cos_lat
        yj = deg2rad(boundary[j][2]) * R
        area += (xj + xi) * (yj - yi)
        j = i
    end
    return area / 2.0   # signed — no abs()
end

sample_idx = n <= sample_n ? collect(1:n) :
             collect(1:max(1, div(n, sample_n)):n)   # deterministic stride sample

signed_areas = [_signed_area_m2(cells[k].boundary) for k in sample_idx]
winding = [a > 0 ? :CCW : (a < 0 ? :CW : :DEGENERATE) for a in signed_areas]

by_sublat = Dict{Symbol, Vector{Symbol}}()
for (pos, k) in enumerate(sample_idx)
    s = sublattice_of[k]
    push!(get!(by_sublat, s, Symbol[]), winding[pos])
end

println("\nWinding distribution by (graph-true) sublattice (n = $(length(sample_idx)) cells):")
all_consistent_within = true
for (s, ws) in sort(collect(by_sublat), by = x -> x[1])
    ccw = count(==(:CCW), ws); cw = count(==(:CW), ws); deg = count(==(:DEGENERATE), ws)
    println("  sublattice $s (n=$(length(ws))): CCW=$ccw  CW=$cw  degenerate=$deg")
    if ccw > 0 && cw > 0
        global all_consistent_within = false
    end
end

println("\nVerdict (Task 1a):")
if !all_consistent_within
    println("  ⚠ INCONSISTENT WINDING WITHIN AT LEAST ONE SUBLATTICE — this would be")
    println("    a genuine pya5/bridge bug (mixed winding within a single sublattice is")
    println("    not expected under any hypothesis). Investigate mesh/a5_bridge.py's")
    println("    boundary vertex ordering immediately.")
elseif length(by_sublat) == 2
    keys_sorted = sort(collect(keys(by_sublat)))
    w1 = unique(filter(!=(:DEGENERATE), by_sublat[keys_sorted[1]]))
    w2 = unique(filter(!=(:DEGENERATE), by_sublat[keys_sorted[2]]))
    if length(w1) == 1 && length(w2) == 1 && w1[1] != w2[1]
        println("  Winding is CONSISTENT within each sublattice, and the two")
        println("  sublattices wind in OPPOSITE directions (sublattice $(keys_sorted[1])=$(w1[1]), " *
                "sublattice $(keys_sorted[2])=$(w2[1])).")
        println("  This is a genuine, systematic, sublattice-correlated asymmetry in")
        println("  raw vertex storage order — flagged as the leading candidate root")
        println("  cause worth checking against Task 1b's flip-rate correlation below.")
    elseif length(w1) == 1 && length(w2) == 1 && w1[1] == w2[1]
        println("  Winding is CONSISTENT within each sublattice, and BOTH sublattices")
        println("  wind the SAME direction ($(w1[1])). No raw winding asymmetry found")
        println("  at this (the whole-polygon) level — any remaining bias must come from")
        println("  a finer-grained effect (e.g. Task 1b's per-edge first-vertex choice,")
        println("  or genuine geometric chirality independent of storage order).")
    else
        println("  Mixed/inconclusive result — see per-sublattice counts above.")
    end
else
    println("  Only one sublattice group present after graph colouring (mesh may be")
    println("  fully disconnected or a single connected bipartite component collapsed")
    println("  unexpectedly) — cannot compare between sublattices. Investigate before")
    println("  trusting this result.")
end

# ────────────────────────────────────────────────────────────────────────────
# TASK 1b — _shared_edge raw vertex ordering / pre-flip normal correlation
# ────────────────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 78)
println("TASK 1b — _shared_edge raw ordering vs sublattice / bearing")
println("=" ^ 78)
println("Replicates _edge_geometry's steps 1-5 inline (same maths, same call")
println("order as production _build_edge_list: cell_i is ALWAYS the lower-ID")
println("cell of the pair — confirmed by direct inspection of the i-loop in")
println("_build_edge_list). Records, per edge: whether the RAW pre-orientation-")
println("flip candidate normal already pointed i→j (no flip needed) or had to")
println("be flipped, plus the compass bearing of the edge.")
println()
println("IMPORTANT — corrected null hypothesis (see v1 run notes): if every cell")
println("shares the SAME winding direction (as found in Task 1a for at least one")
println("sublattice, possibly both), the raw candidate normal is NOT randomly")
println("oriented relative to cell_i vs cell_j — it deterministically points either")
println("toward or away from cell_i's OWN interior for every one of cell_i's edges,")
println("because it is a fixed 90° rotation of an edge vector that always traverses")
println("cell_i's boundary in the SAME rotational sense. Since cell_j sits on the")
println("opposite side of the shared edge from cell_i's interior, this predicts a")
println("flip rate systematically FAR from 50% (close to 0% or close to 100%,")
println("uniformly) as the NORMAL, non-buggy outcome — not a red flag by itself.")
println("The meaningful check is whether the flip rate is UNIFORM across every")
println("sublattice pair and bearing octant (expected, harmless) or VARIES by")
println("heading/pairing (would indicate real, position-dependent bias).")

flip_needed = Bool[]
edge_bearing_deg = Float64[]
edge_sublat_i = Symbol[]
edge_sublat_j = Symbol[]

seen = Set{Tuple{Int,Int}}()
for i in 1:n
    for nb_id in get(adj, ids[i], String[])
        nb_idx = get(id_idx, nb_id, 0)
        nb_idx == 0 && continue
        lo, hi = minmax(i, nb_idx)
        (lo, hi) in seen && continue
        push!(seen, (lo, hi))

        bi, bj = cells[lo], cells[hi]   # cell_i = lo, cell_j = hi (matches _build_edge_list)
        edge = A5Grid._shared_edge(bi.boundary, bj.boundary)
        edge === nothing && continue
        elon1, elat1, elon2, elat2 = edge

        lat0 = (elat1 + elat2) * 0.5
        cos_lat0 = cosd(lat0)
        R = A5Grid._EARTH_R
        to_xy(lon, lat) = (deg2rad(lon) * R * cos_lat0, deg2rad(lat) * R)

        x1, y1 = to_xy(elon1, elat1)
        x2, y2 = to_xy(elon2, elat2)
        ex, ey = x2 - x1, y2 - y1
        e_len = sqrt(ex^2 + ey^2)
        e_len < 1.0 && continue

        nx, ny = -ey / e_len, ex / e_len   # raw candidate normal, pre-flip

        ci_x, ci_y = to_xy(bi.center_lon, bi.center_lat)
        cj_x, cj_y = to_xy(bj.center_lon, bj.center_lat)
        dx, dy = cj_x - ci_x, cj_y - ci_y
        d_len = sqrt(dx^2 + dy^2)
        d_len < 1.0 && continue
        dx_hat, dy_hat = dx / d_len, dy / d_len

        c_raw = dx_hat * nx + dy_hat * ny
        push!(flip_needed, c_raw < 0.0)

        bearing = mod(atand(dx, dy), 360.0)   # 0° = north, clockwise, i→j direction
        push!(edge_bearing_deg, bearing)
        push!(edge_sublat_i, sublattice_of[lo])
        push!(edge_sublat_j, sublattice_of[hi])
    end
end

n_edges_walked = length(flip_needed)
n_flip = count(flip_needed)
overall_flip_rate = n_flip / max(n_edges_walked, 1)
println("\nEdges walked: $n_edges_walked")
println("Flip needed (raw candidate normal pointed j→i, not i→j): " *
        "$n_flip / $n_edges_walked ($(round(100*overall_flip_rate, digits=1))%)")

println("\nFlip rate by sublattice PAIR (sublattice(cell_i) → sublattice(cell_j)):")
pair_keys = sort(unique(zip(edge_sublat_i, edge_sublat_j)))
pair_rates = Float64[]
for (si, sj) in pair_keys
    m = (edge_sublat_i .== si) .& (edge_sublat_j .== sj)
    cnt = count(m)
    cnt == 0 && continue
    fr = count(flip_needed[m]) / cnt
    push!(pair_rates, fr)
    println("  $si→$sj : n=$cnt   flip rate = $(round(100*fr, digits=1))%")
end

println("\nFlip rate by bearing octant (compass direction of i→j edge vector):")
octant_names = ["N","NE","E","SE","S","SW","W","NW"]
octant_rates = Float64[]
for k in 0:7
    lo_b, hi_b = k*45.0, (k+1)*45.0
    m_oct = (edge_bearing_deg .>= lo_b) .& (edge_bearing_deg .< hi_b)
    cnt_oct = count(m_oct)
    cnt_oct == 0 && continue
    fr_oct = count(flip_needed[m_oct]) / cnt_oct
    push!(octant_rates, fr_oct)
    println("  $(octant_names[k+1]) [$(lo_b)°,$(hi_b)°): n=$cnt_oct   flip rate = $(round(100*fr_oct, digits=1))%")
end

pair_spread   = isempty(pair_rates)   ? 0.0 : maximum(pair_rates) - minimum(pair_rates)
octant_spread = isempty(octant_rates) ? 0.0 : maximum(octant_rates) - minimum(octant_rates)

println("\nVerdict (Task 1b):")
println("  Overall flip rate: $(round(100*overall_flip_rate, digits=1))% " *
        "(expected to be far from 50% given consistent per-cell winding — see note above).")
println("  Spread across sublattice pairs: $(round(100*pair_spread, digits=1)) percentage points.")
println("  Spread across bearing octants:  $(round(100*octant_spread, digits=1)) percentage points.")
if pair_spread < 0.10 && octant_spread < 0.10
    println("  Both spreads are small (<10 points): flip rate is UNIFORM across sublattice")
    println("  pairing and heading. No evidence of a systematic vertex-ordering bug —")
    println("  _edge_geometry's orientation-flip is behaving exactly as its docstring")
    println("  claims (correcting a globally-consistent, but position-independent, raw")
    println("  winding artefact). Any remaining directional bias found elsewhere (e.g.")
    println("  the skew_x north/south split from the Pentagon Chirality handover) is NOT")
    println("  explained by this mechanism and should be attributed to something else")
    println("  (genuine A5 chirality in the |V̂| correction term itself, per that handover).")
else
    println("  ⚠ Non-trivial spread found across sublattice pairing and/or bearing —")
    println("  this IS a systematic asymmetry upstream of _edge_geometry (in pya5's")
    println("  vertex ordering, or a5_bridge.py's boundary construction) and should be")
    println("  investigated before Phase B, independent of the orientation-flip logic")
    println("  itself (which is mathematically sound — see")
    println("  FloodA5_GradientCorrection_PentagonChirality_Handover.md §8.3).")
end

# ────────────────────────────────────────────────────────────────────────────
# TASK 2 — Control-volume-centre audit
# ────────────────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 78)
println("TASK 2 — Control-volume-centre audit")
println("=" ^ 78)
println("Confirms: (a) what pya5's cell centre represents relative to the polygon")
println("centroid, and (b) whether all current consumers of \"cell centre\"")
println("(DEM sampling, WLSQ stencil) use the SAME point.")
println()
println("(b) by direct code inspection (not re-derived here, since this is a")
println("static-code fact, not a numeric one):")
println("  - sample_dem_centroid!/sample_dem_mean! (A5Grid.jl) use cell.center_lon/")
println("    cell.center_lat directly as the DEM sample location.")
println("  - _build_wlsq_weights! (FloodModel.jl) uses cell_lons[i]/cell_lats[i]")
println("    (== cell.center_lon/center_lat, copied into FlowState at init) as the")
println("    stencil anchor for every neighbour displacement vector.")
println("  → CONFIRMED CONSISTENT: both consumers anchor on the same point.")

"""
    _polygon_centroid_offset_m(boundary, clon, clat) → (dx_m, dy_m)

Area-weighted polygon centroid, computed in a local equirectangular frame
centred on (clon, clat) (the cell's OWN stored centre — not the polygon's
own vertex-mean latitude, so the returned (dx_m, dy_m) is directly the
offset of the true centroid FROM the stored centre, in metres, in an
east/north frame anchored at the stored centre).
"""
function _polygon_centroid_offset_m(boundary::Vector{Vector{Float64}},
                                     clon::Float64, clat::Float64)
    m = length(boundary)
    m < 3 && return (0.0, 0.0)
    cos_lat0 = cosd(clat)
    R = A5Grid._EARTH_R
    to_xy(lon, lat) = (deg2rad(lon - clon) * R * cos_lat0, deg2rad(lat - clat) * R)

    pts = [to_xy(v[1], v[2]) for v in boundary]
    A = 0.0; Cx = 0.0; Cy = 0.0
    j = m
    for i in 1:m
        xi, yi = pts[i]
        xj, yj = pts[j]
        cross = xj * yi - xi * yj
        A  += cross
        Cx += (xj + xi) * cross
        Cy += (yj + yi) * cross
        j = i
    end
    A *= 0.5
    abs(A) < 1e-9 && return (0.0, 0.0)
    Cx /= (6.0 * A)
    Cy /= (6.0 * A)
    return (Cx, Cy)
end

offsets = [_polygon_centroid_offset_m(cells[k].boundary, cells[k].center_lon, cells[k].center_lat)
           for k in sample_idx]
dists = [sqrt(o[1]^2 + o[2]^2) for o in offsets]

# Equivalent circle diameter, for scale context (2*sqrt(area/pi)).
diam = [2 * sqrt(A5Grid._polygon_area_m2(cells[k].boundary) / pi) for k in sample_idx]
pct_of_diam = [d / dd * 100 for (d, dd) in zip(dists, diam) if dd > 0]

println("\nDistance from stored cell.center_lon/lat to true polygon centroid")
println("(n = $(length(sample_idx)) cells sampled):")
println("  mean   = $(round(mean(dists), digits=4)) m")
println("  median = $(round(median(dists), digits=4)) m")
println("  max    = $(round(maximum(dists), digits=4)) m")
println("  mean equivalent cell diameter = $(round(mean(diam), digits=1)) m")
println("  mean offset as % of cell diameter = $(round(mean(pct_of_diam), digits=3))%")
println("  max  offset as % of cell diameter = $(round(maximum(pct_of_diam), digits=3))%")

println("\nVerdict (Task 2):")
if mean(pct_of_diam) < 0.5
    println("  Stored centre and true polygon centroid agree to well under 0.5% of")
    println("  cell diameter on average — close enough to treat as equivalent for the")
    println("  Phase B diamond-quadrilateral derivation (ci_centre/cj_centre can safely")
    println("  be taken as pya5's cell.center_lon/lat with no separate centroid")
    println("  computation needed).")
else
    println("  ⚠ Stored centre diverges from the true polygon centroid by a")
    println("  non-negligible fraction of cell diameter — Phase B's diamond")
    println("  construction should use an explicitly-computed centroid rather than")
    println("  assuming cell.center_lon/lat is equivalent to it.")
end

println("\n" * "=" ^ 78)
println("Done. See FloodA5_PhaseA_ImplementationScope.md's decision matrix for how")
println("to combine this script's results with Task 3 (correlation) and Task 4")
println("(convergence sweep) before deciding on Phase B.")
println("=" ^ 78)
