# surfacewater/std_device.jl
# --------------------------
# WP5 stage 2 (2026-10): the standard (local-inertial) step as kernel bodies,
# for any compute backend. The CPU step (`step_standard!`) calls the same
# bodies for the per-cell and per-edge phases (WSE gradients, positivity,
# gather, η filter, Phase F, momentum damping; Phase A is `_face_flux_std!`
# in kernels.jl), so there is one implementation of each. The diamond vertex
# values (`_veta_k!`/`_vwet_k!`, as DiamondFlux.jl's diamond_vertex_eta!/
# _wet!) and Phase D (below) have device versions that replicate the CPU
# loops; test/test_backend.jl checks that the device path on the CPU is
# bit-identical to `step_standard!`. `_std_step!` is the device driver,
# phase for phase in `step_standard!`'s order. Two differences, neither of
# which changes results:
#   - the device visits every edge (the CPU's active-edge list visits only the
#     edges that can change, and is bit-identical to the full pass);
#   - open-boundary outflow (Phase D) runs over struct-of-arrays copies of
#     the ghost edges in dependency levels (`_ghost_schedule`): the CPU loop
#     reads a neighbour's volume after earlier ghost edges have changed it,
#     so the order is part of the result; edges in one level share no cell
#     with each other's writes, and the levels run in order, which
#     reproduces the CPU loop exactly. The outflow volumes are summed in
#     edge order by one thread, as the CPU adds them to state.vol_removed.
# Supported on a device: the options of `--scheme recommended` and `v0.2`
# (edge momentum, any Q-centring mode, diamond/diamond-wetdry or two-point
# face gradient, any limiter, local CFL, η filter, momentum damping).
#
# 0.6.0: the bodies are written for any float type (`@fp`, T taken from one
# of the body's arrays) and for both level formulations (`_lev(a)`:
# AbsLevels, the original arithmetic, or RelLevels, which `--precision
# float32` uses) — surfacewater/precision.jl. With Float64 and AbsLevels the
# code generated is the same as before.

# ── Kernel bodies shared with the CPU path ─────────────────────────────────

@inline function _wse_k!(i, a)
    @inbounds a.wse[i] = a.elevation[i] + a.volume[i] / max(a.cell_area[i], 1.0)
    return nothing
end

# a cell's level is usable / the level difference to the neighbour in slot s
# (RelLevels: depth difference plus the bed offset dzadj = z_j − z_i)
@inline _lev_ok(::AbsLevels, a, i) = @inbounds isfinite(a.wse[i])
@inline _lev_ok(::RelLevels, a, i) = @inbounds isfinite(a.elevation[i])
@inline _dlev(::AbsLevels, a, s, i, j) = @inbounds a.wse[j] - a.wse[i]
@inline _dlev(::RelLevels, a, s, i, j) =
    @inbounds (a.water_depth[j] - a.water_depth[i]) + a.dzadj[s, i]

# WLSQ WSE gradient (as _compute_wse_gradients!: non-finite cells and
# neighbours are skipped)
@inline function _grad_wse_k!(i, a)
    T = eltype(a.grad_wse); lev = _lev(a)
    @fp T @inbounds begin
        G = a.grad_wse
        if !_lev_ok(lev, a, i)
            G[1, i] = 0.0; G[2, i] = 0.0
            return nothing
        end
        gx = 0.0; gy = 0.0; nb = a.nbmax
        for s in 1:nb
            j = a.adj[s, i]
            j == 0 && continue
            _lev_ok(lev, a, j) || continue
            dWSE = _dlev(lev, a, s, i, j)
            gx += a.W[s, i] * dWSE
            gy += a.W[nb + s, i] * dWSE
        end
        G[1, i] = gx; G[2, i] = gy
    end
    return nothing
end

# diamond vertex WSE and stencil wetness (as diamond_vertex_eta!/_wet!).
# RelLevels: the vertex value is Σ w·h; the bed part is in zva/zvb (_diamond_dn)
@inline _vlevel(::AbsLevels, a, c) = @inbounds a.wse[c]
@inline _vlevel(::RelLevels, a, c) = @inbounds a.water_depth[c]
@inline function _veta_k!(v, a)
    T = eltype(a.veta); lev = _lev(a)
    @fp T @inbounds if a.vnan[v]
        a.veta[v] = NaN
    else
        s = 0.0
        for k in a.vptr[v]:a.vptr[v + 1] - 1
            s += a.vw[k] * _vlevel(lev, a, a.vcell[k])
        end
        a.veta[v] = s
    end
    return nothing
end
@inline function _vwet_k!(v, a)
    T = eltype(a.water_depth)
    w = true
    @fp T @inbounds for k in a.vptr[v]:a.vptr[v + 1] - 1
        a.water_depth[a.vcell[k]] > HFLOW_THRESHOLD || (w = false; break)
    end
    @inbounds a.vwet[v] = w
    return nothing
end

# positivity (as _positivity_scale_gather!): per-cell outflow totals, then
# one scale factor per donor applied to each of its outflows
@inline function _pos_out_k!(i, a)
    T = eltype(a.edge_vol)
    @fp T begin
    s = 0.0
    @inbounds for k in a.cptr[i]:a.cptr[i + 1] - 1
        e = a.cedge[k]; ev = a.edge_vol[e]
        if ev > 0.0
            a.cj[e] == i && (s += ev)
        elseif ev < 0.0
            a.ci[e] == i && (s -= ev)
        end
    end
    @inbounds a.out[i] = s
    end
    return nothing
end
@inline function _pos_scale_k!(e, a)
    T = eltype(a.edge_vol)
    @fp T @inbounds begin
        ev = a.edge_vol[e]
        iszero(ev) && return nothing
        d = ev > 0.0 ? a.cj[e] : a.ci[e]
        o = a.out[d]
        V = max(a.volume[d], 0.0)
        if o > V
            sc = V / o
            a.edge_vol[e] = ev * sc
            a.flux[e] *= sc
            a.stats && (a.lflag[e] |= LF_POS)               # D1b
        end
    end
    return nothing
end

# Phases B + C (as _gather_apply!, standard flow). The sum is taken in the
# volume's type (Float64 at any working precision): each edge volume enters
# its two cells as the same number, so the exchange conserves mass
@inline function _gather_std_k!(i, a)
    T = eltype(a.edge_vol)
    s = zero(eltype(a.volume))
    @inbounds for k in a.cptr[i]:a.cptr[i + 1] - 1
        e = a.cedge[k]; ev = a.edge_vol[e]
        iszero(ev) && continue
        s = a.ci[e] == i ? s + ev : s - ev
    end
    @fp T @inbounds begin
        A_i = a.cell_area[i]
        A_i < 1.0 && return nothing
        a.volume[i]      = max(0.0, a.volume[i] + s)
        a.water_depth[i] = a.volume[i] / A_i
    end
    return nothing
end

# Phase F (as _reconstruct_qvec!): cell discharge vectors by WLSQ projection
# (Val: the unused discharge array is compiled out — SGS has no `flux` field)
@inline function _qvec_k!(i, a, ::Val{use_flux_Q}) where {use_flux_Q}
    T = eltype(a.qvec_u)
    @fp T begin
    qu = 0.0; qv = 0.0; nb = a.maxnb
    @inbounds for s in 1:nb
        e = a.cei[s, i]
        e == 0 && continue
        q_e = use_flux_Q ? a.flux_Q[e] / max(a.width[e], 1e-9) : a.flux[e]
        obs = -q_e
        qu += a.mom_weights[s, i] * obs
        qv += a.mom_weights[nb + s, i] * obs
    end
    @inbounds a.qvec_u[i] = qu
    @inbounds a.qvec_v[i] = qv
    end
    return nothing
end

# ── η filter (D4) ───────────────────────────────────────────────────────────

@inline function _ef_eta_k!(i, a)                    # AbsLevels only
    @inbounds a.f_eta[i] = a.elevation[i] + a.volume[i] / max(a.cell_area[i], 1.0)
    return nothing
end
@inline function _ef_r_k!(i, a)
    @inbounds a.f_r[i] = _eta_residual(_lev(a), i, a)
    return nothing
end
# face depth above the sill, from the lower of the two levels
@inline _ef_hf(::AbsLevels, a, e, i, j) = @inbounds min(a.f_eta[i], a.f_eta[j]) - a.sill[e]
@inline _ef_hf(::RelLevels, a, e, i, j) =
    @inbounds min(a.zsi[e] + a.water_depth[i], a.zsj[e] + a.water_depth[j])
@inline function _ef_flux_k!(e, a, Cs, dt)
    T = eltype(a.f_flux)
    @fp T @inbounds begin
        i = a.ci[e]; j = a.cj[e]; A = a.cell_area
        hf = _ef_hf(_lev(a), a, e, i, j)
        a.f_flux[e] = hf > 0.0 && A[i] >= 1.0 && A[j] >= 1.0 ?
                      Cs * sqrt(_G * hf) * a.width[e] * (a.f_r[j] - a.f_r[i]) * dt : 0.0
    end
    return nothing
end
@inline function _ef_out_k!(i, a)
    T = eltype(a.f_flux)
    @fp T begin
    o = 0.0
    @inbounds for k in a.cptr[i]:a.cptr[i + 1] - 1
        e = a.cedge[k]; f = a.f_flux[e]
        o += a.ci[e] == i ? max(-f, 0.0) : max(f, 0.0)
    end
    @inbounds a.f_out[i] = o > a.volume[i] ? T(a.volume[i] / o) : 1.0
    end
    return nothing
end
@inline function _ef_scale_k!(e, a)
    T = eltype(a.f_flux)
    @fp T @inbounds begin
        f = a.f_flux[e]
        iszero(f) || (a.f_flux[e] = f * (f > 0.0 ? a.f_out[a.cj[e]] : a.f_out[a.ci[e]]))
    end
    return nothing
end
@inline function _ef_apply_k!(i, a)
    s = zero(eltype(a.volume))                       # summed in the volume's type, as _gather_std_k!
    @inbounds for k in a.cptr[i]:a.cptr[i + 1] - 1
        e = a.cedge[k]; f = a.f_flux[e]
        iszero(f) && continue
        s = a.ci[e] == i ? s + f : s - f
    end
    @inbounds if !iszero(s)
        a.volume[i] = max(zero(s), a.volume[i] + s)
        a.water_depth[i] = a.volume[i] / a.cell_area[i]
    end
    return nothing
end

# ── Momentum damping (D6′) ──────────────────────────────────────────────────

@inline function _md_face_k!(e, a)
    T = eltype(a.m_h)
    @fp T @inbounds begin
        i = a.ci[e]; j = a.cj[e]
        ηi, ηj, sill = _edge_levels(_lev(a), a, e, i, j)
        hmin = min(ηi, ηj) - sill
        if hmin > 0.0 && a.m_tail[e] != 0
            hf = max(ηi, ηj) - sill
            nx = a.nf_x[e]; ny = a.nf_y[e]; dx = a.dx_m[e] / a.L[e]; dy = a.dy_m[e] / a.L[e]
            qt = 0.5 * ((a.qvec_u[i] + a.qvec_u[j]) * -ny + (a.qvec_v[i] + a.qvec_v[j]) * nx)
            a.m_h[e] = hf
            a.m_u[e] = (-a.flux[e] * (nx * dx + ny * dy) + qt * (-ny * dx + nx * dy)) / hf
            a.m_w[e] = min(hmin / MOM_DAMP_HW, 1.0)
        else
            a.m_h[e] = 0.0; a.m_u[e] = 0.0; a.m_w[e] = 0.0
        end
    end
    return nothing
end
@inline function _md_vertex_k!(v, a, Cm, wave::Bool, dt)
    T = eltype(a.m_psi)
    @fp T @inbounds begin
        ps = 0.0
        if a.m_interior[v]
            Γ = 0.0; Av = 0.0; ω = 1.0; hs = 0.0; us = 0.0; m = 0
            for k in a.m_vptr[v]:a.m_vptr[v + 1] - 1
                e = a.m_vedge[k]
                ω = min(ω, a.m_w[e])
                Γ += a.m_vsign[k] * a.m_u[e] * a.L[e]
                Av += a.width[e] * a.L[e] / 4
                hs += a.m_h[e]; us += abs(a.m_u[e]); m += 1
            end
            if ω > 0.0 && Av > 0.0
                hv = hs / m
                Vref = wave ? sqrt(_G * hv) : us / m
                ν = Cm * Vref * sqrt(Av)
                ps = -ν * dt * hv * (Γ / Av) * ω
            end
        end
        a.m_psi[v] = ps
    end
    return nothing
end
@inline function _md_apply_k!(e, a)
    @inbounds begin
        t = a.m_tail[e]
        if t != 0
            dQ = a.m_psi[a.m_head[e]] - a.m_psi[t]
            iszero(dQ) || (a.flux[e] -= dQ / a.width[e])
        end
    end
    return nothing
end

# ── Device-only bodies ──────────────────────────────────────────────────────

# Phase D for ghost edge k (as _apply_ghost_fluxes_standard!), launched level
# by level (`_ghost_pass!`); the edge's outflow volume goes to g_dv[k]
# Levels at a boundary face: the cell's level and bed, its interior neighbour's
# level and bed (NaN without one), the sill. RelLevels: relative to the sill,
# from the offsets g_zc = z_cell − sill and g_znb = z_neighbour − sill.
@inline function _ghost_levels(::AbsLevels, a, k, ci, nb)
    @inbounds begin
        wse_ci = a.elevation[ci] + a.volume[ci] / max(a.cell_area[ci], 1.0)
        wse_nb = nb > 0 ? a.elevation[nb] + a.volume[nb] / max(a.cell_area[nb], 1.0) : NaN
        return wse_ci, a.elevation[ci], wse_nb, nb > 0 ? a.elevation[nb] : NaN, a.g_sill[k]
    end
end
@inline function _ghost_levels(::RelLevels, a, k, ci, nb)
    T = eltype(a.g_zc)
    @inbounds return (a.g_zc[k] + a.water_depth[ci], a.g_zc[k],
                      nb > 0 ? a.g_znb[k] + a.water_depth[nb] : convert(T, NaN), a.g_znb[k], zero(T))
end

@inline function _ghost_std_k!(k, a, dt_in)
    T = eltype(a.g_fprev); dt = convert(T, dt_in)
    @fp T @inbounds begin
        a.g_dv[k] = 0.0
        bc = a.g_bc[k]
        bc === Closed && return nothing
        ci = a.g_cell[k]
        a.cell_area[ci] < 1.0 && return nothing
        isnan(a.elevation[ci]) && return nothing
        nb = a.g_nb[k]
        wse_ci, z_ci, wse_nb, z_nb, sill = _ghost_levels(_lev(a), a, k, ci, nb)
        isnan(wse_ci) && return nothing
        a.obnormal && nb > 0 && (wse_nb = _open_nb_wse(wse_ci, z_ci, z_nb))
        wse_ghost = _ghost_wse(wse_ci, sill, bc, wse_nb, a.g_nbL[k], a.g_L[k])
        wse_ghost == -Inf && return nothing
        Q_out, q_stored = _bates_ghost_flux(a.g_fprev[k], wse_ci, wse_ghost, sill,
                                            a.g_width[k], a.g_L[k], a.manning_n[ci], dt,
                                            a.bpos ? 0.0 : a.water_depth[ci])
        dV = convert(eltype(a.volume), Q_out * dt)
        dV = min(dV, a.bpos ? a.volume[ci] : a.volume[ci] / a.donor_div)
        dV = max(dV, 0.0)
        a.volume[ci]      = max(0.0, a.volume[ci] - dV)
        a.water_depth[ci] = a.volume[ci] / max(a.cell_area[ci], 1.0)
        a.g_dv[k]         = dV
        a.g_fprev[k]      = q_stored
    end
    return nothing
end

"""
    _ghost_schedule(g_cell, g_nb, g_bc) → (order, lptr)

Dependency levels for Phase D on a device. Ghost edge k reads the volumes of
its cell and its interior neighbour and writes its cell's; it is placed one
level after the latest earlier edge that touched either cell, so edges in a
level are independent and running the levels in order reproduces the serial
CPU loop exactly. Closed faces do nothing and are left out. `order[lptr[l]:
lptr[l+1]-1]` are the edges of level l (in edge order).
"""
function _ghost_schedule(g_cell, g_nb, g_bc)
    last = Dict{Int,Int}(); lev = zeros(Int, length(g_cell))
    for k in eachindex(g_cell)
        g_bc[k] === Closed && continue
        c = g_cell[k]; nb = g_nb[k]
        l = get(last, c, 0)
        nb > 0 && (l = max(l, get(last, nb, 0)))
        lev[k] = l + 1
        last[c] = lev[k]
        nb > 0 && (last[nb] = max(get(last, nb, 0), lev[k]))
    end
    nlev = maximum(lev; init = 0)
    order = [k for l in 1:nlev for k in eachindex(lev) if lev[k] == l]
    lptr = ones(Int, nlev + 1)
    for l in 1:nlev; lptr[l + 1] = lptr[l] + count(==(l), lev); end
    return order, lptr
end

@inline function _ghost_lvl_k!(j, body::F, a, off, dt) where {F}
    @inbounds body(a.g_order[off + j], a, dt)
    return nothing
end
@inline function _ghost_sum_k!(_, a)
    r = a.removed[1]
    @inbounds for k in eachindex(a.g_dv)
        r += a.g_dv[k]
    end
    a.removed[1] = r
    return nothing
end

"Phase D on a device: the levels of `lptr` in order, then the ordered outflow sum."
function _ghost_pass!(be::ComputeBackend, body::F, a, lptr::Vector{Int}, dt) where {F}
    for l in 1:length(lptr) - 1
        foreach_index!(be, lptr[l + 1] - lptr[l], _ghost_lvl_k!, body, a, lptr[l] - 1, dt)
    end
    on_host_order!(be, 1, _ghost_sum_k!, a)
    return nothing
end

# local CFL rate √(g h)/ℓ per cell (as _cfl_dt_local, standard flow)
@inline function _cfl_rate_k!(i, a)
    T = eltype(a.cfl)
    @fp T @inbounds begin
        h = a.water_depth[i]
        r = 0.0
        if h > 1e-4
            ℓ = a.cell_length[i]
            ℓ > 0.0 && (r = sqrt(_G * h) / ℓ)
        end
        a.cfl[i] = r
    end
    return nothing
end

@inline function _depth_std_k!(i, a)
    @inbounds a.cell_area[i] >= one(eltype(a.cell_area)) &&
              (a.water_depth[i] = a.volume[i] / a.cell_area[i])
    return nothing
end

# ── Device arrays and driver ────────────────────────────────────────────────

"""
    _rel_offsets(state, T, use_diamond) → NamedTuple

The bed offsets `RelLevels` needs (surfacewater/precision.jl), differenced in
Float64: per edge the two cells' beds relative to the sill (`zsi`, `zsj`)
and, with a diamond table, the two vertices' interpolated beds (`zva`,
`zvb`); per adjacency slot the neighbour's bed relative to the cell's
(`dzadj`); per boundary face the cell's and its interior neighbour's bed
relative to the face's sill (`g_zc`, `g_znb`).
"""
function _rel_offsets(state::FlowState, T, use_diamond::Bool)
    E = state.edges; z = state.elevation; ne = E.n_edges; n = length(z)
    adj = state.adj_matrix; G = state.ghost_edges
    zsi = [z[E.cell_i[e]] - E.sill[e] for e in 1:ne]
    zsj = [z[E.cell_j[e]] - E.sill[e] for e in 1:ne]
    dzadj = [adj[s, i] > 0 ? z[adj[s, i]] - z[i] : 0.0 for s in 1:size(adj, 1), i in 1:n]
    zva = Float64[]; zvb = Float64[]
    if use_diamond
        zv = [T.vnan[v] ? NaN : sum(T.vw[k] * z[T.vcell[k]] for k in T.vptr[v]:T.vptr[v + 1] - 1; init = 0.0)
              for v in 1:length(T.vertices)]
        zva = [T.valid[e] ? zv[T.va_idx[e]] - E.sill[e] : NaN for e in 1:ne]
        zvb = [T.valid[e] ? zv[T.vb_idx[e]] - E.sill[e] : NaN for e in 1:ne]
    end
    return (zsi = zsi, zsj = zsj, dzadj = dzadj, zva = zva, zvb = zvb,
            g_zc = [z[g.cell_index] - g.sill for g in G],
            g_znb = [g.interior_nb_idx > 0 ? z[g.interior_nb_idx] - g.sill : NaN for g in G])
end

"""
    _std_device_args(state; lev = AbsLevels()) → NamedTuple

Everything `_std_step!` reads or writes, as host arrays (references) and
scalars; `to_device` copies it to the device. Builds the workspaces of the
active options first (η filter, momentum damping, diamond vertex values).
`lev = RelLevels()` adds the bed offsets (`_rel_offsets`).
"""
function _std_device_args(state::FlowState; lev = AbsLevels())
    ws = _workspace!(state); E = state.edges
    cptr, cedge = _cell_csr!(state, ws)
    T = state.diamond_table
    use_diamond = state.gradient_correction && T !== nothing &&
                  state.face_flux_method in (:diamond, :diamond_wetdry)
    nv = use_diamond ? length(T.vertices) : 0
    length(ws.veta) == nv || resize!(ws.veta, nv)
    length(ws.vwet) == nv || resize!(ws.vwet, nv)
    state.eta_filter_cs > 0.0 && _eta_filter_ws!(state, ws)
    state.mom_damping_cs > 0.0 && _mom_damp_ws!(state, ws)
    G = state.ghost_edges
    M = _std_arrays(state, ws, ws.flux_prev)
    dtab = use_diamond ? merge(M.dtab, (valid = Vector{Bool}(T.valid),)) : nothing
    rel = lev isa RelLevels ? _rel_offsets(state, T, use_diamond) : (;)
    return merge(M, rel, (lev = lev,
        dtab = dtab, stats = false, lflag = UInt8[], lhf = Float64[],
        flux_prev = ws.flux_prev, wse = ws.wse, out = ws.out,
        adj = state.adj_matrix, W = state.wlsq_weights, nbmax = size(state.adj_matrix, 1),
        cei = state.cell_edge_index, mom_weights = state.mom_weights, maxnb = state.max_neighbours,
        flux_Q = E.flux_Q, dx_m = E.dx_m, dy_m = E.dy_m, cptr = cptr, cedge = cedge,
        cell_length = length(state.cell_length) == length(state.cell_ids) ? state.cell_length :
                      sqrt.(state.cell_area),
        cfl = zeros(length(state.cell_ids)),
        vptr = use_diamond ? T.vptr : Int[1], vcell = use_diamond ? T.vcell : Int[],
        vw = use_diamond ? T.vw : Float64[], vnan = use_diamond ? T.vnan : Bool[],
        f_dx = ws.f_dx, f_dy = ws.f_dy, f_eta = ws.f_eta, f_r = ws.f_r, f_out = ws.f_out,
        f_flux = ws.f_flux,
        m_tail = ws.m_tail, m_head = ws.m_head, m_vptr = ws.m_vptr, m_vedge = ws.m_vedge,
        m_vsign = ws.m_vsign, m_interior = Vector{Bool}(ws.m_interior), m_psi = ws.m_psi,
        m_u = ws.m_u, m_w = ws.m_w, m_h = ws.m_h,
        g_cell = [g.cell_index for g in G], g_width = [g.width for g in G], g_L = [g.L for g in G],
        g_sill = [g.sill for g in G], g_fprev = [g.flux_prev for g in G],
        g_nb = [g.interior_nb_idx for g in G], g_nbL = [g.interior_nb_L for g in G],
        g_bc = BCType[b for b in state.ghost_cell_bc], removed = [state.vol_removed],
        g_dv = zeros(length(G)), g_order = _ghost_schedule([g.cell_index for g in G],
                                                            [g.interior_nb_idx for g in G],
                                                            state.ghost_cell_bc)[1],
        obnormal = state.open_boundary === :normal, bpos = BOUNDARY_POSITIVITY[],
        n = length(state.cell_ids), ne = E.n_edges, ng = length(G), nv = nv,
        nvm = length(ws.m_psi),
        opts = StdOpts(state), use_diamond = use_diamond,
        wetdry = state.face_flux_method === :diamond_wetdry,
        gc = state.gradient_correction, edge_mom = state.momentum_model === :edge,
        positivity = state.limiter !== :legacy,
        need_qvec = state.momentum_model !== :edge || state.friction_model === :vector ||
                    state.q_centre_mode === :consistent,
        eta_cs = state.eta_filter_cs, md_cs = state.mom_damping_cs,
        md_wave = state.mom_damping_scale === :wave))
end

"""
    _std_step!(be, a, dt, glptr)

One standard-flow step on the arrays `a` (`_std_device_args`, on any
backend), in `step_standard!`'s phase order; `glptr` are the Phase D level
bounds (`_ghost_schedule`).
"""
function _std_step!(be::ComputeBackend, a, dt64::Float64, glptr::Vector{Int})
    dt = convert(eltype(a.flux), dt64)
    abs_lev = _lev(a) isa AbsLevels                  # RelLevels needs no z + h arrays
    if a.gc
        abs_lev && foreach_index!(be, a.n, _wse_k!, a)
        foreach_index!(be, a.n, _grad_wse_k!, a)
    end
    if a.use_diamond
        foreach_index!(be, a.nv, _veta_k!, a)
        a.wetdry && foreach_index!(be, a.nv, _vwet_k!, a)
    end
    a.edge_mom && copyto!(a.flux_prev, a.flux)
    foreach_index!(be, a.ne, _face_flux_std!, a.opts, a, dt)
    if a.positivity
        foreach_index!(be, a.n, _pos_out_k!, a)
        foreach_index!(be, a.ne, _pos_scale_k!, a)
    end
    foreach_index!(be, a.n, _gather_std_k!, a)
    a.ng > 0 && _ghost_pass!(be, _ghost_std_k!, a, glptr, dt)
    if a.eta_cs > 0.0
        abs_lev && foreach_index!(be, a.n, _ef_eta_k!, a)
        foreach_index!(be, a.n, _ef_r_k!, a)
        foreach_index!(be, a.ne, _ef_flux_k!, a, a.eta_cs, dt)
        foreach_index!(be, a.n, _ef_out_k!, a)
        foreach_index!(be, a.ne, _ef_scale_k!, a)
        foreach_index!(be, a.n, _ef_apply_k!, a)
    end
    a.need_qvec && foreach_index!(be, a.n, _qvec_k!, a, Val(false))
    if a.md_cs > 0.0
        foreach_index!(be, a.ne, _md_face_k!, a)
        foreach_index!(be, a.nvm, _md_vertex_k!, a, a.md_cs, a.md_wave, dt)
        foreach_index!(be, a.ne, _md_apply_k!, a)
    end
    return nothing
end

"Local-CFL Δt on the arrays `a` (as `_cfl_dt_local`, standard flow)."
function _std_cfl(be::ComputeBackend, a, courant::Float64)
    foreach_index!(be, a.n, _cfl_rate_k!, a)
    r = _maximum(a.cfl)
    return r > 0.0 ? courant / r : 60.0
end
