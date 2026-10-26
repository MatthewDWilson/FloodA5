# surfacewater/kernels.jl
# -----------------------
# WP4 stage 2 (2026-10): the time-step phases as kernels — one function per
# edge or per cell, taking flat arrays (a NamedTuple) and the solver options as
# type parameters, so the option branches are resolved at compile time and the
# same bodies can later be launched on a GPU (WP5). Each kernel must give
# exactly the result of the frozen reference (surfacewater/reference_steps.jl);
# see test/test_wp4_stage2.jl.

"""
    StdOpts{GC, FF, MOM, QC, VF, POS, FIXB}

Standard-flow options as type parameters: gradient correction (Bool), face
gradient (:none, :diamond, :diamond_wetdry — :none when no diamond table is in
use), momentum model (:edge, :cell, :blend), Q-centring mode (:legacy,
:oriented, :consistent), vector friction (Bool), positivity limiting (Bool),
per-edge volume limiter Fix B (Bool).
"""
struct StdOpts{GC, FF, MOM, QC, VF, POS, FIXB} end

function StdOpts(state::FlowState)
    use_diamond = state.gradient_correction && state.diamond_table !== nothing &&
                  state.face_flux_method in (:diamond, :diamond_wetdry)
    ff = use_diamond ? state.face_flux_method : :none
    return StdOpts{state.gradient_correction, ff, state.momentum_model, state.q_centre_mode,
                   state.friction_model === :vector, state.limiter !== :legacy,
                   state.limiter !== :positivity}()
end

"""
    _std_arrays(state, ws, flux_prev) → NamedTuple

The arrays and scalars `_face_flux_std!` reads or writes, gathered once per step
(references, not copies).
"""
function _std_arrays(state::FlowState, ws::StepWorkspace, flux_prev)
    E = state.edges
    return (ci = E.cell_i, cj = E.cell_j, width = E.width, L = E.L, cos_theta = E.cos_theta,
            sill = E.sill, nf_x = E.nf_x, nf_y = E.nf_y, skew_x = E.skew_x, skew_y = E.skew_y,
            col_i = E.collinear_i, col_j = E.collinear_j, edges = E,
            elevation = state.elevation, volume = state.volume, cell_area = state.cell_area,
            water_depth = state.water_depth, manning_n = state.manning_n,
            qvec_u = state.qvec_u, qvec_v = state.qvec_v, grad_wse = state.grad_wse,
            flux = E.flux, flux_prev = flux_prev, edge_vol = ws.edge_vol,
            veta = ws.veta, vwet = ws.vwet, dtab = state.diamond_table,
            theta = state.q_centre_theta, alpha = state.gradient_correction_alpha,
            max_nb = Float64(state.max_neighbours), donor_div = _donor_edge_divisor(state),
            lflag = _lf_flag(), lhf = _lf_hf())               # D1b (--limiter-stats)
end

"""
    _face_flux_std!(e, opts, M, dt)

Phase A of `step_standard!` for edge `e`: face discharge (stored in M.flux[e],
Fix C) and signed edge volume after the per-edge donor cap (M.edge_vol[e]).
Transcribed from the stage-1 loop body; bit-identical to it.
"""
@inline function _face_flux_std!(e::Int, ::StdOpts{GC, FF, MOM, QC, VF, POS, FIXB}, M,
                                 dt::Float64) where {GC, FF, MOM, QC, VF, POS, FIXB}
    @inbounds begin
        ci = M.ci[e]; cj = M.cj[e]
        # Degenerate edge guard (NaN in geometry or elevation)
        if (isnan(M.width[e])        || isnan(M.L[e])            ||
            isnan(M.cos_theta[e])    || isnan(M.sill[e])         ||
            isnan(M.elevation[ci])   || isnan(M.elevation[cj]))
            M.edge_vol[e] = 0.0
            M.flux[e] = 0.0
            LIMITER_STATS[] && (M.lflag[e] = 0x00)
            return nothing
        end
        wse_ci = M.elevation[ci] + M.volume[ci] / max(M.cell_area[ci], 1.0)
        wse_cj = M.elevation[cj] + M.volume[cj] / max(M.cell_area[cj], 1.0)
        # Dry edge: both kernels return (0, 0)
        if max(wse_ci, wse_cj) - M.sill[e] <= HFLOW_THRESHOLD
            M.flux[e] = 0.0
            M.edge_vol[e] = 0.0
            LIMITER_STATS[] && (M.lflag[e] = 0x00)
            return nothing
        end
        depth_donor = !FIXB ? 0.0 :
                      wse_ci >= wse_cj ? M.water_depth[ci] : M.water_depth[cj]
        q_proj = MOM === :edge ? 0.0 :
                 -0.5 * ((M.qvec_u[ci] + M.qvec_u[cj]) * M.nf_x[e] +
                         (M.qvec_v[ci] + M.qvec_v[cj]) * M.nf_y[e])
        q_prev_eff = if MOM === :cell
            q_proj
        elseif MOM === :blend
            M.theta * M.flux_prev[e] + (1.0 - M.theta) * q_proj
        elseif QC === :consistent
            _q_centred_consistent(M.flux_prev, e, M.edges, M.qvec_u, M.qvec_v, M.theta)
        else
            _q_centred(M.flux_prev, e, M.col_i[e], M.col_j[e], M.theta)
        end
        q_t = VF ? 0.5 * ((M.qvec_v[ci] + M.qvec_v[cj]) * M.nf_x[e] -
                          (M.qvec_u[ci] + M.qvec_u[cj]) * M.nf_y[e]) : 0.0
        if GC
            dWSE_n = NaN
            if FF === :diamond
                dWSE_n = diamond_dWSE_n_v(M.dtab, e, wse_ci, wse_cj, M.veta, M.L[e])
            elseif FF === :diamond_wetdry
                dWSE_n = diamond_stencil_wet_v(M.dtab, e, M.vwet) ?
                         diamond_dWSE_n_v(M.dtab, e, wse_ci, wse_cj, M.veta, M.L[e]) :
                         M.cos_theta[e] * (wse_ci - wse_cj)
            end
            if isnan(dWSE_n)
                gx_f = 0.5 * (M.grad_wse[1, ci] + M.grad_wse[1, cj])
                gy_f = 0.5 * (M.grad_wse[2, ci] + M.grad_wse[2, cj])
                Vhat_dot = gx_f * M.skew_x[e] + gy_f * M.skew_y[e]
                dWSE_n   = M.cos_theta[e] * (wse_ci - wse_cj) - M.alpha * M.L[e] * Vhat_dot
            end
            h_flow = max(wse_ci, wse_cj) - M.sill[e]
            Q, q_stored = _bates_flux_limited_corrected(
                q_prev_eff, h_flow, dWSE_n, M.width[e], M.L[e],
                min(M.manning_n[ci], M.manning_n[cj]), dt, depth_donor;
                vol_limit_divisor = M.max_nb, q_t = q_t, vector_friction = VF)
        else
            Q, q_stored = _bates_flux_limited(
                q_prev_eff, wse_ci, wse_cj, M.sill[e], M.width[e], M.L[e], M.cos_theta[e],
                min(M.manning_n[ci], M.manning_n[cj]), dt, depth_donor;
                vol_limit_divisor = M.max_nb, q_t = q_t, vector_friction = VF)
        end
        M.flux[e] = q_stored                   # Fix C
        ev = Q * dt
        ev_raw = ev
        if !POS && !iszero(ev)                 # per-edge donor cap
            ev = ev > 0.0 ? min(ev,  M.volume[cj] / M.donor_div) :
                            max(ev, -M.volume[ci] / M.donor_div)
        end
        M.edge_vol[e] = ev
        if LIMITER_STATS[]                     # D1b: which limiters bound (no effect on results)
            hs = max(max(wse_ci, wse_cj) - M.sill[e], 1e-6)
            fl = LF_WET
            qmax = hs * sqrt(_G * hs) * FROUDE_LIMIT
            if VF
                hypot(q_stored, q_t) >= qmax * (1 - 1e-12) && (fl |= LF_FROUDE)
            else
                abs(q_stored) == qmax && (fl |= LF_FROUDE)
            end
            depth_donor > 0.0 &&
                abs(q_stored) == (depth_donor * M.width[e]) / (M.max_nb * dt) && (fl |= LF_FIXB)
            ev != ev_raw && (fl |= LF_DONOR)
            GC && FF === :diamond_wetdry && !diamond_stencil_wet_v(M.dtab, e, M.vwet) && (fl |= LF_WD)
            M.lflag[e] = fl; M.lhf[e] = hs
        end
    end
    return nothing
end

# Launch over all edges with the options as compile-time constants (function
# barrier: `opts` is concrete inside).
function _face_flux_std_all!(opts::StdOpts, M, ne::Int, dt::Float64)
    @batch per=core for e in 1:ne
        _face_flux_std!(e, opts, M, dt)
    end
    return nothing
end

# Active-edge list mode (WP4 stage 3): the same kernel over act[1:n_act].
function _face_flux_std_list!(opts::StdOpts, M, act::Vector{Int}, n_act::Int, dt::Float64)
    @batch per=core for j in 1:n_act
        _face_flux_std!(act[j], opts, M, dt)
    end
    return nothing
end

"""
    SGSOpts{GC, QCC, VF, POS}

SGS (R–A kernel) options as type parameters: gradient correction (Bool),
consistent Q-centring (Bool; otherwise `_q_centred` with whatever partner
signs the edge list holds), vector friction (Bool), positivity limiting (Bool).
"""
struct SGSOpts{GC, QCC, VF, POS} end

SGSOpts(state::FlowState) =
    SGSOpts{state.gradient_correction, state.q_centre_mode === :consistent,
            state.friction_model === :vector, state.limiter !== :legacy}()

function _sgs_arrays(state::FlowState, ws::StepWorkspace)
    E = state.edges
    return (ci = E.cell_i, cj = E.cell_j, width = E.width, L = E.L, cos_theta = E.cos_theta,
            sill = E.sill, nf_x = E.nf_x, nf_y = E.nf_y, skew_x = E.skew_x, skew_y = E.skew_y,
            col_i = E.collinear_i, col_j = E.collinear_j, edges = E,
            elevation = state.elevation, volume = state.volume, manning_n = state.manning_n,
            qvec_u = state.qvec_u, qvec_v = state.qvec_v, grad_wse = state.grad_wse,
            tables = state.sgs_tables, slot_i = ws.slot_i, slot_j = ws.slot_j, wse = ws.wse,
            flux_Q = E.flux_Q, flux_Q_prev = ws.flux_prev, edge_vol = ws.edge_vol,
            theta = state.q_centre_theta, alpha = state.gradient_correction_alpha,
            donor_div = _donor_edge_divisor(state),
            lflag = _lf_flag(), lhf = _lf_hf())               # D1b (--limiter-stats)
end

"""
    _face_flux_sgs!(e, opts, M, dt)

Phase A of `step_sgs!` for edge `e` (R–A kernel; Q in M.flux_Q[e], edge volume
after the donor cap with Fix C back-propagation in M.edge_vol[e]). Transcribed
from the stage-1 loop body; bit-identical to it.
"""
@inline function _face_flux_sgs!(e::Int, ::SGSOpts{GC, QCC, VF, POS}, M,
                                 dt::Float64) where {GC, QCC, VF, POS}
    @inbounds begin
        ci = M.ci[e]; cj = M.cj[e]
        if isnan(M.elevation[ci]) || isnan(M.elevation[cj])
            M.edge_vol[e] = 0.0                # flux_Q left as is (as before)
            LIMITER_STATS[] && (M.lflag[e] = 0x00)
            return nothing
        end
        ti = M.tables[ci]; tj = M.tables[cj]
        z_sill = M.sill[e]
        isnan(z_sill) && (z_sill = min(ti.z_min, tj.z_min))
        if isnan(M.width[e]) || isnan(M.L[e]) || isnan(M.cos_theta[e])
            M.edge_vol[e] = 0.0
            LIMITER_STATS[] && (M.lflag[e] = 0x00)
            return nothing
        end
        # Bug 48/49: a dry cell's effective WSE clears the sill and its own floor
        wse_ci_eff = M.volume[ci] > 0.0 ? M.wse[ci] : max(z_sill, ti.z_min)
        wse_cj_eff = M.volume[cj] > 0.0 ? M.wse[cj] : max(z_sill, tj.z_min)
        if M.volume[ci] <= 0.0 && M.volume[cj] <= 0.0
            M.edge_vol[e] = 0.0
            M.flux_Q[e] = 0.0
            LIMITER_STATS[] && (M.lflag[e] = 0x00)
            return nothing
        end
        wse_flow = max(wse_ci_eff, wse_cj_eff)
        A_i, R_i = A5Grid.edge_area_radius(ti, M.slot_i[e], wse_flow)
        A_j, R_j = A5Grid.edge_area_radius(tj, M.slot_j[e], wse_flow)
        A_edge = 0.5 * (A_i + A_j)
        R_edge = 0.5 * (R_i + R_j)
        Q_prev_eff = QCC ?
            _q_centred_consistent(M.flux_Q_prev, e, M.edges, M.qvec_u, M.qvec_v, M.theta;
                                  volumetric = true) :
            _q_centred(M.flux_Q_prev, e, M.col_i[e], M.col_j[e], M.theta)
        Q_t = VF ? 0.5 * ((M.qvec_v[ci] + M.qvec_v[cj]) * M.nf_x[e] -
                          (M.qvec_u[ci] + M.qvec_u[cj]) * M.nf_y[e]) * M.width[e] : 0.0
        Q_new = if GC
            gx_f = 0.5 * (M.grad_wse[1, ci] + M.grad_wse[1, cj])
            gy_f = 0.5 * (M.grad_wse[2, ci] + M.grad_wse[2, cj])
            Vhat_dot = gx_f * M.skew_x[e] + gy_f * M.skew_y[e]
            dWSE_n   = M.cos_theta[e] * (wse_ci_eff - wse_cj_eff) - M.alpha * M.L[e] * Vhat_dot
            h_flow = wse_flow - z_sill
            _manning_flux_ra_corrected(Q_prev_eff, h_flow, A_edge, R_edge, dWSE_n, M.L[e],
                                       0.5 * (M.manning_n[ci] + M.manning_n[cj]),
                                       dt; Q_t = Q_t, vector_friction = VF)
        else
            _manning_flux_ra(Q_prev_eff, wse_ci_eff, wse_cj_eff, z_sill, A_edge, R_edge,
                             M.L[e], M.cos_theta[e],
                             0.5 * (M.manning_n[ci] + M.manning_n[cj]),
                             dt; Q_t = Q_t, vector_friction = VF)
        end
        M.flux_Q[e] = Q_new
        ev = Q_new * dt
        ev_raw = ev
        if !POS && !iszero(ev)                 # donor cap, Fix C back-propagation
            if ev > 0.0
                ev_capped = min(ev, M.volume[cj] / M.donor_div)
                ev_capped < ev && (M.flux_Q[e] = ev_capped / dt)
            else
                ev_capped = max(ev, -M.volume[ci] / M.donor_div)
                ev_capped > ev && (M.flux_Q[e] = ev_capped / dt)
            end
            ev = ev_capped
        end
        M.edge_vol[e] = ev
        if LIMITER_STATS[]                     # D1b (no effect on results)
            hs = wse_flow - z_sill
            M.lflag[e] = (hs > HFLOW_THRESHOLD && A_edge > 1e-6 ? LF_WET : 0x00) |
                         (ev != ev_raw ? LF_DONOR : 0x00)
            M.lhf[e] = hs
        end
    end
    return nothing
end

function _face_flux_sgs_all!(opts::SGSOpts, M, ne::Int, dt::Float64)
    @batch per=core for e in 1:ne
        _face_flux_sgs!(e, opts, M, dt)
    end
    return nothing
end

function _face_flux_sgs_list!(opts::SGSOpts, M, act::Vector{Int}, n_act::Int, dt::Float64)
    @batch per=core for j in 1:n_act
        _face_flux_sgs!(act[j], opts, M, dt)
    end
    return nothing
end
