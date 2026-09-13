# mesh/DiamondFlux.jl
# ------------------
# Phase C, Level 0-2: production implementation of the diamond
# face-normal derivative (B3) and its vertex-reconstruction input (B2).
#
# See:
#   FloodA5_PhaseB_B1_FaceFluxEquation.md   -- proves dWSE_n = -L*(g_n) is
#                                              exact for ANY gradient source,
#                                              so this file needs no new flux
#                                              kernel: it plugs directly into
#                                              the existing
#                                              _bates_flux_limited_corrected.
#   FloodA5_PhaseB_B2_VertexReconstruction.md -- vertex reconstruction proof.
#   FloodA5_PhaseB_B3_DiamondGradient.md      -- diamond gradient proof.
#
# Selected via --face-flux-method diamond (default: legacy). When legacy,
# this file's functions are loaded but never called -- zero behavioural
# change to existing runs.
#
# Deliberately NOT touching EdgeList: all diamond geometry lives in the new
# DiamondFluxTable struct below, keyed 1:n_edges identically to EdgeList,
# rather than adding more fields to an already-large, widely-constructed
# struct. FlowState gains two new fields (face_flux_method, diamond_table)
# to hold this; see FloodModel.jl's FlowState struct and
# initialise_flow_model for the wiring.

using LinearAlgebra

const VERTEX_COND_DEGENERATE_THRESHOLD  = 1.0e11   # B2 §3.4 real-mesh finding
const DIAMOND_DEGENERATE_AREA_THRESHOLD = 1.0      # m^2; matches B3 audit-script guard

"""
    VertexRecord

One mesh vertex's B2 reconstruction data: η_v = dot(weights, η[cell_indices]).

  fallback ∈ (:none, :distance_avg, :single, :degenerate)
    :none          k>=3, well-conditioned WLS solve (the normal case)
    :distance_avg  k==2, underdetermined for a full 2D fit -- distance-
                   weighted average instead (NOT exact for a general
                   linear field, only along the line through the two cells)
    :single        k<=1 -- trivial passthrough, exact only for constant fields
    :degenerate    k>=3 but condition number too high -- weights zeroed,
                   same safe-fallback convention as _build_wlsq_weights!
"""
struct VertexRecord
    cell_indices :: Vector{Int}
    weights      :: Vector{Float64}
    cond_number  :: Float64
    fallback     :: Symbol
end

"""
    eta_at_vertex(rec, eta_cell) -> Float64

η_v = dot(weights, η[cell_indices]). NaN for an empty record (k=0; should
not occur for a genuine mesh vertex).
"""
@inline function eta_at_vertex(rec::VertexRecord, eta_cell::Vector{Float64})::Float64
    isempty(rec.cell_indices) && return NaN
    return dot(rec.weights, view(eta_cell, rec.cell_indices))
end

"""
    _vertex_weights(dx, w) -> (weights, cond)

WLS constant-term weight vector for a k>=3 vertex stencil.
FloodA5_PhaseB_B2_VertexReconstruction.md §2.1.
"""
function _vertex_weights(dx::Matrix{Float64}, w::Vector{Float64})::Tuple{Vector{Float64},Float64}
    k    = size(dx, 1)
    Wd   = Diagonal(w)
    Amat = hcat(ones(k), dx)      # k x 3   [1, dx, dy]
    M    = Amat' * Wd * Amat       # 3 x 3
    Rm   = Amat' * Wd              # 3 x k
    sol  = M \ Rm                  # 3 x k
    return sol[1, :], cond(M)
end

"""
    _build_vertex_record(cell_indices, dx_centred) -> VertexRecord

`dx_centred[row, :]` is the row-th contributing cell's position MINUS the
vertex position (already vertex-centred by the caller). Implements the
k>=3 / k==2 / k<=1 fallback ladder, B2 §4.
"""
function _build_vertex_record(cell_indices::Vector{Int},
                               dx_centred::Matrix{Float64})::VertexRecord
    k = length(cell_indices)

    if k <= 1
        w = k == 1 ? [1.0] : Float64[]
        return VertexRecord(cell_indices, w, NaN, :single)
    end

    if k == 2
        d2_1 = dx_centred[1,1]^2 + dx_centred[1,2]^2
        d2_2 = dx_centred[2,1]^2 + dx_centred[2,2]^2
        w1 = 1.0 / max(d2_1, 1e-6)
        w2 = 1.0 / max(d2_2, 1e-6)
        wsum = w1 + w2
        return VertexRecord(cell_indices, [w1/wsum, w2/wsum], NaN, :distance_avg)
    end

    ws = Vector{Float64}(undef, k)
    for row in 1:k
        d2 = dx_centred[row,1]^2 + dx_centred[row,2]^2
        ws[row] = 1.0 / max(d2, 1e-6)
    end
    weights, cnum = _vertex_weights(dx_centred, ws)

    if !isfinite(cnum) || cnum > VERTEX_COND_DEGENERATE_THRESHOLD
        return VertexRecord(cell_indices, zeros(k), cnum, :degenerate)
    end
    return VertexRecord(cell_indices, weights, cnum, :none)
end

"""
    _diamond_normal_coeffs(x_i, x_va, x_j, x_vb, n_hat) -> (c_i, c_va, c_j, c_vb, area)

Static geometry coefficients such that, for any η at the four diamond
vertices in order (x_i, x_va, x_j, x_vb):

    (∂η/∂n)_e = c_i*η_i + c_va*η_va + c_j*η_j + c_vb*η_vb

FloodA5_PhaseB_B3_DiamondGradient.md §2-3. `n_hat` must be the SAME
oriented face normal already stored as (edges.nf_x[e], edges.nf_y[e]) --
not recomputed independently (see EdgeList's own field docs on why two
orientation conventions must never be allowed to drift apart).

Provably invariant to (x_i,x_va,x_j,x_vb) traversal winding direction
(B3 §3.3) -- no orientation enforcement needed on the raw `_shared_edge`
vertex order.
"""
function _diamond_normal_coeffs(x_i::Tuple{Float64,Float64}, x_va::Tuple{Float64,Float64},
                                 x_j::Tuple{Float64,Float64}, x_vb::Tuple{Float64,Float64},
                                 n_hat::Tuple{Float64,Float64})
    P = (x_i, x_va, x_j, x_vb)
    nx, ny = n_hat

    area2 = 0.0
    for a in 1:4
        b = a % 4 + 1
        xa, ya = P[a]; xb, yb = P[b]
        area2 += xa*yb - xb*ya
    end
    area = area2 / 2.0

    coeffs = zeros(4)
    for a in 1:4
        b = a % 4 + 1
        xa, ya = P[a]; xb, yb = P[b]
        dxe, dye = xb - xa, yb - ya
        # rot90(dxe,dye) = (dye,-dxe); dotted with n_hat and split 0.5/0.5
        # onto the two endpoints (trapezoidal boundary average)
        contrib = 0.5 * (nx*dye - ny*dxe)
        coeffs[a] += contrib
        coeffs[b] += contrib
    end
    coeffs ./= area

    return coeffs[1], coeffs[2], coeffs[3], coeffs[4], area
end

"""
    DiamondFluxTable

Precomputed, mesh-build-time diamond geometry, one record per EdgeList
edge (same 1:n_edges indexing as `edges.cell_i`/`edges.cell_j`).

`va_key`/`vb_key` are (round(lon,9), round(lat,9)) lookup keys into
`vertex_table`, matching `A5Grid._shared_edge`'s own tolerance
convention exactly (do not change the rounding digits independently of
`_shared_edge`).

`valid[e] == false` means this edge's diamond could not be built (no
shared edge found, degenerate face normal, or degenerate diamond area --
all expected to be extremely rare: the B3 real-mesh audit found 1
degenerate edge out of 74,130 on the Carlisle mesh). The caller MUST
check `valid[e]` and fall back to the legacy WLSQ+skew `dWSE_n`
construction for any edge where it is false -- there is no dedicated
kernel fallback inside this file.
"""
struct DiamondFluxTable
    vertex_table :: Dict{Tuple{Float64,Float64}, VertexRecord}
    va_key       :: Vector{Tuple{Float64,Float64}}
    vb_key       :: Vector{Tuple{Float64,Float64}}
    c_i          :: Vector{Float64}
    c_va         :: Vector{Float64}
    c_j          :: Vector{Float64}
    c_vb         :: Vector{Float64}
    valid        :: BitVector
end

"""
    build_diamond_flux_table(cells, edges) -> DiamondFluxTable

Mesh-build-time construction, called once from `initialise_flow_model`
when `face_flux_method == :diamond`.

  cells  mesh.cells :: Vector{A5Cell}
  edges  the already-built EdgeList (reuses `edges.nf_x`/`nf_y` -- the
         production oriented face normal -- and `edges.cell_i`/`cell_j`;
         does not recompute orientation independently)

Uses a single global local-equirectangular frame for the whole mesh,
centred on the mean of all cell centres. This matches the convention
already validated to ~1e-8 m/m absolute precision on real meshes up to
Carlisle's full extent in the B2/B3 audit scripts (audit_vertex_valence.jl,
audit_diamond_gradient.jl) -- simpler than EdgeList's own per-edge moving
frame (`dx_m`/`dy_m`), and that per-edge choice was made for EdgeList's
different purpose (Q-centred / cell-momentum geometry), not because a
single global frame is inadequate for this one.
"""
function build_diamond_flux_table(cells::Vector{A5Cell}, edges::EdgeList)::DiamondFluxTable
    n = length(cells)
    lon0 = mean(c.center_lon for c in cells)
    lat0 = mean(c.center_lat for c in cells)
    cos_lat0 = cosd(lat0)
    R = A5Grid._EARTH_R
    to_xy(lon, lat) = (deg2rad(lon - lon0) * R * cos_lat0, deg2rad(lat - lat0) * R)
    _key(lon, lat) = (round(lon, digits=9), round(lat, digits=9))

    cell_xy = [to_xy(c.center_lon, c.center_lat) for c in cells]

    # ── vertex table (Layer 1, B2) ────────────────────────────────────
    vertex_cells = Dict{Tuple{Float64,Float64}, Vector{Int}}()
    for i in 1:n
        for pt in cells[i].boundary
            k = _key(pt[1], pt[2])
            push!(get!(vertex_cells, k, Int[]), i)
        end
    end
    for (k, v) in vertex_cells
        vertex_cells[k] = unique(v)
    end

    vertex_table = Dict{Tuple{Float64,Float64}, VertexRecord}()
    for (vkey, cidx) in vertex_cells
        vx, vy = to_xy(vkey[1], vkey[2])
        k = length(cidx)
        dx = Matrix{Float64}(undef, k, 2)
        for (row, ci) in enumerate(cidx)
            cx, cy = cell_xy[ci]
            dx[row, 1] = cx - vx
            dx[row, 2] = cy - vy
        end
        vertex_table[vkey] = _build_vertex_record(cidx, dx)
    end

    n_deg = count(r -> r.fallback == :degenerate,    values(vertex_table))
    n_k1  = count(r -> r.fallback == :single,         values(vertex_table))
    n_k2  = count(r -> r.fallback == :distance_avg,   values(vertex_table))
    @info "Diamond flux: vertex table built — $(length(vertex_table)) vertices " *
          "($n_k2 k=2 fallback, $n_k1 k<=1 fallback, $n_deg ill-conditioned)"

    # ── per-edge diamond records (Layer 2, B3) ────────────────────────
    ne = edges.n_edges
    va_keys = Vector{Tuple{Float64,Float64}}(undef, ne)
    vb_keys = Vector{Tuple{Float64,Float64}}(undef, ne)
    c_is  = zeros(Float64, ne)
    c_vas = zeros(Float64, ne)
    c_js  = zeros(Float64, ne)
    c_vbs = zeros(Float64, ne)
    valid = falses(ne)

    n_no_shared_edge    = 0
    n_degenerate_normal = 0
    n_degenerate_area   = 0

    for e in 1:ne
        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        shared = _shared_edge(cells[ci].boundary, cells[cj].boundary)
        if shared === nothing
            va_keys[e] = (0.0, 0.0)
            vb_keys[e] = (0.0, 0.0)
            n_no_shared_edge += 1
            continue
        end
        valon, valat, vblon, vblat = shared
        va_keys[e] = _key(valon, valat)
        vb_keys[e] = _key(vblon, vblat)

        nx, ny = edges.nf_x[e], edges.nf_y[e]
        if nx == 0.0 && ny == 0.0
            n_degenerate_normal += 1
            continue   # valid[e] stays false
        end

        x_i  = cell_xy[ci]
        x_j  = cell_xy[cj]
        x_va = to_xy(valon, valat)
        x_vb = to_xy(vblon, vblat)

        c_i, c_va, c_j, c_vb, area = _diamond_normal_coeffs(x_i, x_va, x_j, x_vb, (nx, ny))
        if abs(area) < DIAMOND_DEGENERATE_AREA_THRESHOLD
            n_degenerate_area += 1
            continue
        end

        c_is[e]  = c_i
        c_vas[e] = c_va
        c_js[e]  = c_j
        c_vbs[e] = c_vb
        valid[e] = true
    end

    n_valid = count(valid)
    @info "Diamond flux: $n_valid/$ne edges have a valid diamond record " *
          "($n_no_shared_edge no shared edge, $n_degenerate_normal degenerate normal, " *
          "$n_degenerate_area degenerate area — these fall back to the legacy " *
          "gradient-correction path automatically)"

    return DiamondFluxTable(vertex_table, va_keys, vb_keys, c_is, c_vas, c_js, c_vbs, valid)
end

"""
    diamond_dWSE_n(table, e, wse_ci, wse_cj, wse_all, L_e) -> Float64

Runtime evaluation of `dWSE_n` for edge `e`, for direct use as the
`dWSE_n` argument to the EXISTING `_bates_flux_corrected` /
`_bates_flux_limited_corrected` kernels — no new flux kernel is needed.
See FloodA5_PhaseB_B1_FaceFluxEquation.md §3.3: `dWSE_n = -L·g_n` is
exact for any exact gradient source, including this one.

`wse_all` is the full per-cell WSE array (already materialised by the
caller for the legacy gradient-correction path — reused here, not
recomputed).

Returns NaN if `table.valid[e]` is false, or if either vertex's
reconstruction hit the `:degenerate` (ill-conditioned) fallback — the
caller must check for NaN and use the legacy `dWSE_n` construction for
that one edge in that case.
"""
@inline function diamond_dWSE_n(table::DiamondFluxTable, e::Int,
                                 wse_ci::Float64, wse_cj::Float64,
                                 wse_all::Vector{Float64}, L_e::Float64)::Float64
    table.valid[e] || return NaN
    eta_va = eta_at_vertex(table.vertex_table[table.va_key[e]], wse_all)
    eta_vb = eta_at_vertex(table.vertex_table[table.vb_key[e]], wse_all)
    (isnan(eta_va) || isnan(eta_vb)) && return NaN
    dn = table.c_i[e]*wse_ci + table.c_va[e]*eta_va +
         table.c_j[e]*wse_cj + table.c_vb[e]*eta_vb
    return -L_e * dn
end
