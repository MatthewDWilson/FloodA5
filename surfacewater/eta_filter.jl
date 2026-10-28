# surfacewater/eta_filter.jl
# --------------------------
# D4 (deep-channel damping, 2026-10): conservative free-surface high-pass
# filter, `--eta-filter C_s` (standard flow; 0 = off, the default).
#
# After the volume update, water is exchanged across wet faces in proportion
# to the difference of the high-frequency WSE residual
#     r_i = η_i − a_i,
# where a_i is the value at cell i's centre of the weighted affine fit
# a + g·d to the WSE of the cell and its neighbours. Constants and planes give
# r = 0, so a lake at rest and a uniform slope are untouched; a checkerboard of
# amplitude δ on a square grid has r = ±1.6δ. Including the cell itself keeps
# the fit determined: a cell with only two wet neighbours (three points) is
# fitted exactly and gets r = 0, and r goes to 0 continuously as neighbours
# dry out — the neighbour-only fit would invent a residual there.
# The face exchange (m³/s, positive into cell_i, like edge_vol)
#     F_f = C_s · √(g h_f) · w_f · (r_j − r_i),   h_f = min(η_i, η_j) − z_sill > 0
# is a diffusivity D_f = C_s·c_f·L_f times w_f/L_f: its rate, C_s·c_f/ℓ, is in
# physical units, so the damping per unit time does not depend on Δt (unlike
# Q-centring's (1 − θ)/Δt). Applied to the water surface after the volume
# update, it raises rather than lowers the forward–backward stability limit
# (HYDRAULICS.md §15.8), as long as the damping per step stays below 1: on
# a square grid it damps the checkerboard by a factor 1 − 12.8·C_s·C per
# step at Courant number C, so C_s < 0.11 at C = 0.7.
#
# Weights in the fit: inverse squared distance × the face-wetness ramp
# clamp(h_f / ETA_FILTER_HW, 0, 1) for neighbours, so they join and leave the
# fit continuously and dry neighbours (bed above the water) never enter it;
# the cell itself gets the mean geometric weight of its neighbour slots. The
# exchange is pairwise (conservative); a donor cell's filter outflows are
# scaled together so that it cannot lose more than it holds.

const ETA_FILTER_HW = 0.05      # m: face depth over which a neighbour's fit weight ramps to 1

# Neighbour offsets per adjacency slot (m, local equirectangular frame as in
# _build_wlsq_weights!) and the filter work arrays, built on first use.
function _eta_filter_ws!(state::FlowState, ws::StepWorkspace)
    n = length(state.cell_ids); nb = size(state.adj_matrix, 1); ne = state.edges.n_edges
    size(ws.f_dx) == (nb, n) && length(ws.f_flux) == ne && return ws
    R = A5Grid._EARTH_R
    dx = zeros(nb, n); dy = zeros(nb, n)
    @inbounds for i in 1:n
        c = cosd(state.cell_lats[i])
        for s in 1:nb
            j = state.adj_matrix[s, i]
            j == 0 && continue
            dx[s, i] = deg2rad(state.cell_lons[j] - state.cell_lons[i]) * R * c
            dy[s, i] = deg2rad(state.cell_lats[j] - state.cell_lats[i]) * R
        end
    end
    ws.f_dx = dx; ws.f_dy = dy
    ws.f_eta = zeros(n); ws.f_r = zeros(n); ws.f_out = zeros(n); ws.f_flux = zeros(ne)
    return ws
end

# r_i: η_i minus the weighted affine fit to cell i and its wet neighbours,
# evaluated at cell i (fit on δ_j = η_j − η_i for precision; cell i is the
# point d = 0, δ = 0). A tiny Tikhonov term on the gradient keeps one-neighbour
# and collinear stencils defined.
@inline function _eta_residual(i::Int, η::AbstractVector, V::AbstractVector,
                               adj::AbstractMatrix, cei::AbstractMatrix, sill::AbstractVector,
                               dxm::AbstractMatrix, dym::AbstractMatrix)   # Polyester passes PtrArrays
    @inbounds begin
        V[i] > 0.0 || return 0.0
        ηi = η[i]
        W = 0.0; sx = 0.0; sy = 0.0; sd = 0.0; wg = 0.0; ng = 0
        sxx = 0.0; sxy = 0.0; syy = 0.0; sxd = 0.0; syd = 0.0
        for s in 1:size(adj, 1)
            j = adj[s, i]; j == 0 && continue
            e = cei[s, i]; e == 0 && continue
            hf = min(ηi, η[j]) - sill[e]
            hf > 0.0 || continue                      # also false for NaN
            dx = dxm[s, i]; dy = dym[s, i]; d2 = dx * dx + dy * dy
            d2 > 0.0 || continue
            wg += 1.0 / d2; ng += 1
            w = min(hf / ETA_FILTER_HW, 1.0) / d2
            δ = η[j] - ηi
            W += w; sx += w * dx; sy += w * dy; sd += w * δ
            sxx += w * dx * dx; sxy += w * dx * dy; syy += w * dy * dy
            sxd += w * dx * δ; syd += w * dy * δ
        end
        W > 0.0 || return 0.0
        W += wg / ng                                  # the cell itself (d = 0, δ = 0)
        mx = sx / W; my = sy / W; md = sd / W
        λ = 1e-10 * (sxx + syy) / W                  # regularises exact degeneracy only
        cxx = sxx / W - mx * mx + λ; cyy = syy / W - my * my + λ; cxy = sxy / W - mx * my
        cxd = sxd / W - mx * md;     cyd = syd / W - my * md
        det = cxx * cyy - cxy * cxy
        gx = (cyy * cxd - cxy * cyd) / det
        gy = (cxx * cyd - cxy * cxd) / det
        return -(md - gx * mx - gy * my)              # η_i − (η_i + fitted δ at d = 0)
    end
end

"""
    _eta_filter!(state, ws, dt)

Phase D2 of `step_standard!` when `state.eta_filter_cs > 0` (see the file
header). Threaded per cell / per edge with fixed-order CSR gathers, so the
result does not depend on the thread count.
"""
function _eta_filter!(state::FlowState, ws::StepWorkspace, dt::Float64)
    _eta_filter_ws!(state, ws)
    Cs = state.eta_filter_cs
    n = length(state.cell_ids); E = state.edges; ne = E.n_edges
    η = ws.f_eta; r = ws.f_r; F = ws.f_flux; sc = ws.f_out
    adj = state.adj_matrix; cei = state.cell_edge_index; dxm = ws.f_dx; dym = ws.f_dy
    z = state.elevation; V = state.volume; A = state.cell_area; depth = state.water_depth
    ci = E.cell_i; cj = E.cell_j; sill = E.sill; wid = E.width
    cptr, cedge = _cell_csr!(state, ws)

    @batch per=core for i in 1:n
        @inbounds η[i] = z[i] + V[i] / max(A[i], 1.0)
    end
    @batch per=core for i in 1:n
        @inbounds r[i] = _eta_residual(i, η, V, adj, cei, sill, dxm, dym)
    end
    @batch per=core for e in 1:ne
        @inbounds begin
            i = ci[e]; j = cj[e]
            hf = min(η[i], η[j]) - sill[e]
            F[e] = hf > 0.0 && A[i] >= 1.0 && A[j] >= 1.0 ?
                   Cs * sqrt(_G * hf) * wid[e] * (r[j] - r[i]) * dt : 0.0
        end
    end
    # positivity: all filter outflows of a cell together ≤ its volume
    @batch per=core for i in 1:n
        o = 0.0
        @inbounds for k in cptr[i]:cptr[i + 1] - 1
            e = cedge[k]; f = F[e]
            o += ci[e] == i ? max(-f, 0.0) : max(f, 0.0)
        end
        @inbounds sc[i] = o > V[i] ? V[i] / o : 1.0
    end
    @batch per=core for e in 1:ne
        @inbounds begin
            f = F[e]
            iszero(f) || (F[e] = f * (f > 0.0 ? sc[cj[e]] : sc[ci[e]]))   # f > 0: cj donates
        end
    end
    @batch per=core for i in 1:n
        s = 0.0
        @inbounds for k in cptr[i]:cptr[i + 1] - 1
            e = cedge[k]; f = F[e]
            iszero(f) && continue
            s = ci[e] == i ? s + f : s - f
        end
        @inbounds if !iszero(s)
            V[i] = max(0.0, V[i] + s)
            depth[i] = V[i] / A[i]
        end
    end
    return nothing
end
