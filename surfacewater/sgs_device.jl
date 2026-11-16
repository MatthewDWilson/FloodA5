# surfacewater/sgs_device.jl
# --------------------------
# WP5 stage 3 (2026-10): the SGS (sub-grid, R–A kernel) step on any compute
# backend. As std_device.jl for the standard step: kernel bodies shared with
# the CPU path where the CPU loops are per cell (WSE, gather), the device
# driver `_sgs_step!` in `step_sgs!`'s phase order, and Phase D in
# dependency levels (std_device.jl, `_ghost_schedule`).
#
# The per-cell hypsometric tables (`Vector{SGSTable}`, each holding its own
# arrays) cannot go to a device as they are. `_sgs_flat` packs them into
# padded arrays (bins × cells, bins × slots × cells), and `_tab(tables, i)`
# returns cell i's table — the SGSTable itself on the CPU, a NamedTuple of
# views with the same field names on a device — so the A5Grid lookups
# (wse_from_volume, face_apk) and the face kernel run unchanged.

_tab(T::AbstractVector, i) = @inbounds T[i]
@inline function _tab(T::NamedTuple, i)
    @inbounds nb = T.nbins[i]
    return (elev_bins = view(T.elev, 1:nb, i), vol_curve = view(T.vol, 1:nb, i),
            area_curve = view(T.area, 1:nb, i), cell_area = (@inbounds T.ca[i]),
            z_min = (@inbounds T.zmin[i]), z_max = (@inbounds T.zmax[i]),
            face_z0 = view(T.fz0, :, i), face_zr = view(T.fzr, :, i),
            face_A = view(T.fA, :, :, i), face_P = view(T.fP, :, :, i), face_K = view(T.fK, :, :, i),
            wet_cx = view(T.wcx, 1:nb, i), wet_cy = view(T.wcy, 1:nb, i))
end

"Pack the SGS tables into padded arrays (see the file header)."
function _sgs_flat(tables)
    n = length(tables)
    nb = maximum((length(t.elev_bins) for t in tables); init = 0)
    ns = maximum((length(t.face_z0) for t in tables); init = 0)
    nf = maximum((size(t.face_A, 1) for t in tables); init = 0)
    elev = zeros(nb, n); vol = zeros(nb, n); area = zeros(nb, n)
    wcx = zeros(nb, n); wcy = zeros(nb, n)
    fz0 = fill(NaN, ns, n); fzr = zeros(ns, n)
    fA = zeros(nf, ns, n); fP = zeros(nf, ns, n); fK = zeros(nf, ns, n)
    nbins = zeros(Int, n); ca = zeros(n); zmin = zeros(n); zmax = zeros(n)
    for (i, t) in enumerate(tables)
        k = length(t.elev_bins); nbins[i] = k
        elev[1:k, i] .= t.elev_bins; vol[1:k, i] .= t.vol_curve; area[1:k, i] .= t.area_curve
        wcx[1:k, i] .= t.wet_cx; wcy[1:k, i] .= t.wet_cy
        s = length(t.face_z0)
        size(t.face_A, 1) == nf || error("SGS face tables of different lengths cannot be packed")
        fz0[1:s, i] .= t.face_z0; fzr[1:s, i] .= t.face_zr
        fA[:, 1:s, i] .= t.face_A; fP[:, 1:s, i] .= t.face_P; fK[:, 1:s, i] .= t.face_K
        ca[i] = t.cell_area; zmin[i] = t.z_min; zmax[i] = t.z_max
    end
    return (elev = elev, vol = vol, area = area, wcx = wcx, wcy = wcy,
            fz0 = fz0, fzr = fzr, fA = fA, fP = fP, fK = fK,
            nbins = nbins,
            ca = ca, zmin = zmin, zmax = zmax)
end

# ── Kernel bodies ───────────────────────────────────────────────────────────

@inline function _sgs_wse_k!(i, a)
    @inbounds begin
        t = _tab(a.tables, i)
        w = wse_from_volume(t, a.volume[i])
        a.wse[i] = w
        a.wcx[i], a.wcy[i] = isnan(t.z_min) ? (0.0, 0.0) : wet_centroid(t, w)
        a.wfrac[i] = isnan(t.z_min) ? 1.0 : wetted_area_from_wse(t, w) / t.cell_area
    end
    return nothing
end

# Phases B + C (as _gather_apply!, SGS): no cell-area guard; depth from the table
@inline function _gather_sgs_k!(i, a)
    s = 0.0
    @inbounds for k in a.cptr[i]:a.cptr[i + 1] - 1
        e = a.cedge[k]; ev = a.edge_vol[e]
        iszero(ev) && continue
        s = a.ci[e] == i ? s + ev : s - ev
    end
    @inbounds begin
        a.volume[i] = max(0.0, a.volume[i] + s)
        t = _tab(a.tables, i)
        a.water_depth[i] = isnan(t.z_min) ? 0.0 :
                           max(0.0, wse_from_volume(t, a.volume[i]) - t.z_min)
    end
    return nothing
end

# Phase D for ghost edge k (as _apply_ghost_fluxes_sgs!), level by level
@inline function _ghost_sgs_k!(k, a, dt)
    @inbounds begin
        a.g_dv[k] = 0.0
        bc = a.g_bc[k]
        bc === Closed && return nothing
        ci = a.g_cell[k]
        tbl = _tab(a.tables, ci)
        isnan(tbl.z_min) && return nothing
        wse_ci = wse_from_volume(tbl, a.volume[ci])
        isnan(wse_ci) && return nothing
        nb = a.g_nb[k]
        wse_nb = nb > 0 && nb <= a.n ? wse_from_volume(_tab(a.tables, nb), a.volume[nb]) : NaN
        a.obnormal && nb > 0 && (wse_nb = _open_nb_wse(wse_ci, a.elevation[ci], a.elevation[nb]))
        wse_ghost = _ghost_wse(wse_ci, a.g_sill[k], bc, wse_nb, a.g_nbL[k], a.g_L[k])
        wse_ghost == -Inf && return nothing
        if a.g_isQ[k]
            A_ge, R_ge = _sgs_ghost_area_radius(tbl, a.g_slot[k], wse_ci, a.g_sill[k], a.g_width[k], a.k2d)
            Q_out, q_stored = _manning_ghost_flux(a.g_fprev[k], wse_ci, wse_ghost, a.g_sill[k],
                                                  A_ge, R_ge, a.g_L[k], a.manning_n[ci], dt)
        else
            Q_out, q_stored = _bates_ghost_flux(a.g_fprev[k], wse_ci, wse_ghost, a.g_sill[k],
                                                a.g_width[k], a.g_L[k], a.manning_n[ci], dt,
                                                a.bpos ? 0.0 : a.water_depth[ci])
        end
        dV = Q_out * dt
        dV = min(dV, a.bpos ? a.volume[ci] : a.volume[ci] / a.donor_div)
        dV = max(dV, 0.0)
        a.volume[ci] = max(0.0, a.volume[ci] - dV)
        a.water_depth[ci] = max(0.0, wse_from_volume(tbl, a.volume[ci]) - tbl.z_min)
        a.g_dv[k] = dV
        a.g_fprev[k] = q_stored
    end
    return nothing
end

# local CFL rate (as _cfl_dt_local, SGS: the hydraulic depth, _sgs_cfl_depth)
@inline function _cfl_rate_sgs_k!(i, a)
    @inbounds begin
        r = 0.0
        h = a.water_depth[i]
        if h > 1e-4
            t = _tab(a.tables, i)
            if !isnan(t.z_min)
                aw = max(wetted_area_from_wse(t, t.z_min + h), SGS_CFL_AREA_FLOOR * t.cell_area)
                r = _sgs_face_rate(i, a.cptr, a.cedge, a.face_A, a.L, aw)
                h = _sgs_cfl_depth(t, a.volume[i], h)
            end
            if h > 1e-4
                ℓ = a.cell_length[i]
                ℓ > 0.0 && (r = max(r, sqrt(_G * h) / ℓ))
            end
        end
        a.cfl[i] = r
    end
    return nothing
end

# ── Device arrays and driver ────────────────────────────────────────────────

"""
    _sgs_device_args(state) → NamedTuple

Everything `_sgs_step!` reads or writes (host references and scalars, with the
tables packed by `_sgs_flat`); `to_device` copies it to the device.
"""
function _sgs_device_args(state::FlowState)
    ws = _edge_slots!(state, _workspace!(state)); E = state.edges
    cptr, cedge = _cell_csr!(state, ws)
    G = state.ghost_edges
    M = _sgs_arrays(state, ws)
    return merge(M, (
        tables = _sgs_flat(state.sgs_tables), stats = false, lflag = UInt8[], lhf = Float64[],
        out = ws.out, cptr = cptr, cedge = cedge,
        adj = state.adj_matrix, W = state.wlsq_weights, nbmax = size(state.adj_matrix, 1),
        cei = state.cell_edge_index, mom_weights = state.mom_weights, maxnb = state.max_neighbours,
        water_depth = state.water_depth, cell_area = state.cell_area,
        cell_length = length(state.cell_length) == length(state.cell_ids) ? state.cell_length :
                      sqrt.(state.cell_area),
        cfl = zeros(length(state.cell_ids)),
        g_cell = [g.cell_index for g in G], g_width = [g.width for g in G], g_L = [g.L for g in G],
        g_sill = [g.sill for g in G], g_fprev = [g.flux_prev for g in G],
        g_slot = [g.slot for g in G], k2d = SGS_CONVEYANCE[] === :twod,
        g_isQ = [g.is_Q_flux for g in G],
        g_nb = [g.interior_nb_idx for g in G], g_nbL = [g.interior_nb_L for g in G],
        g_bc = BCType[b for b in state.ghost_cell_bc], removed = [state.vol_removed],
        g_dv = zeros(length(G)), g_order = _ghost_schedule([g.cell_index for g in G],
                                                            [g.interior_nb_idx for g in G],
                                                            state.ghost_cell_bc)[1],
        obnormal = state.open_boundary === :normal, bpos = BOUNDARY_POSITIVITY[],
        n = length(state.cell_ids), ne = E.n_edges, ng = length(G),
        opts = SGSOpts(state), gc = state.gradient_correction,
        positivity = state.limiter !== :legacy,
        need_qvec = state.friction_model === :vector || state.q_centre_mode === :consistent))
end

"One SGS step on the arrays `a` (any backend), in `step_sgs!`'s phase order."
function _sgs_step!(be::ComputeBackend, a, dt::Float64, glptr::Vector{Int})
    foreach_index!(be, a.n, _sgs_wse_k!, a)
    a.gc && foreach_index!(be, a.n, _grad_wse_k!, a)
    copyto!(a.flux_Q_prev, a.flux_Q)
    foreach_index!(be, a.ne, _face_flux_sgs!, a.opts, a, dt)
    if a.positivity                    # scales flux_Q (Fix C); the alias is made on the
        p = merge(a, (flux = a.flux_Q,))  # device arrays, so both names are one array
        foreach_index!(be, a.n, _pos_out_k!, p)
        foreach_index!(be, a.ne, _pos_scale_k!, p)
    end
    foreach_index!(be, a.n, _gather_sgs_k!, a)
    a.ng > 0 && _ghost_pass!(be, _ghost_sgs_k!, a, glptr, dt)
    a.need_qvec && foreach_index!(be, a.n, _qvec_k!, a, Val(true))
    return nothing
end

"Local-CFL Δt on the arrays `a` (as `_cfl_dt_local`, SGS)."
function _sgs_cfl(be::ComputeBackend, a, courant::Float64)
    foreach_index!(be, a.n, _cfl_rate_sgs_k!, a)
    r = _maximum(a.cfl)
    return r > 0.0 ? courant / r : 60.0
end
