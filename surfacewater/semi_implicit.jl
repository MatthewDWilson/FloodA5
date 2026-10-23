# surfacewater/semi_implicit.jl
# -----------------------------
# D3 (damping work, 2026-10): research option — a semi-implicit (θ-method)
# free-surface force for the standard solver, `--free-surface semi-implicit
# --fs-alpha α`. Explicit (α = 0, the Bates scheme) is the default and is not
# touched by this file.
#
# Per face the explicit kernel is
#     q⁺ = [q_eff − g·h·Δt·D_f(ηⁿ)] / (1 + F),
# with D_f the face term of the water-surface gradient (per unit length):
#   kind 1  two-point   D = a·(η_i − η_j), a = 1/max(L·max(cosθ, 0.1), 1)
#                       (legacy kernel) or cosθ/max(L, 1) (corrected kernel,
#                       two-point fallback of diamond-wetdry)
#   kind 2  diamond     D = diamond dWSE_n(η)/max(L, 1)
#   kind 3  WLSQ        D = cosθ·(η_i − η_j)/max(L, 1) + r (r from the cell
#                       gradients, kept explicit; edges without a diamond record)
# and F the semi-implicit friction factor. The semi-implicit form evaluates
# D at (1−α)·ηⁿ + α·ηⁿ⁺¹ — the whole face gradient, diamond included, since
# treating only its symmetric two-point part implicitly with the rest explicit
# was found to *increase* the growth of the asymmetric diamond modes:
#     q⁺ = q* − β·D_f(ηⁿ⁺¹),  q* = [q_eff − g·h·Δt·(1−α)·D_f(ηⁿ) − g·h·Δt·r] / (1 + F),
#     β  = g·h·Δt·α / (1 + F).
# Continuity, A_i·ηⁿ⁺¹ = A_i·ηⁿ + Δt·Σ_f s_if·w_f·q⁺_f, gives
#     A_i·η_i + Δt·Σ_f s_if·w_f·β_f·D_f(η) = A_i·ηⁿ_i + Δt·Σ_f s_if·w_f·q*_f,
# symmetric positive definite for two-point faces (solved by CG) and
# non-symmetric with diamond faces (BiCGStab); both matrix-free and Jacobi
# preconditioned. The face discharges q⁺ then get the explicit kernels' Froude
# cap and volume limiter, and the volume update is the explicit path's (donor
# cap or positivity, per-cell gather): mass balance is exact whatever the
# solver tolerance.

# (SIWork is defined in FloodModel.jl, before StepWorkspace, which holds it.)

function _si_work!(w::SIWork, n::Int, ne::Int, nv::Int)
    for f in (:qs, :beta, :hf, :qt, :dd, :te, :aj)
        v = getfield(w, f); length(v) == ne || resize!(v, ne)
    end
    length(w.kind) == ne || resize!(w.kind, ne)
    for f in (:eta0, :b, :x, :r, :z, :p, :Ap, :dg, :rh, :vv, :yy, :ss, :tt)
        v = getfield(w, f); length(v) == n || resize!(v, n)
    end
    length(w.vtmp) == nv || resize!(w.vtmp, nv)
    return w
end

# face term D_f(x) for edge e (x a cell field; vt its diamond vertex values)
@inline function _si_D(e, W, M, x, vt)
    @inbounds begin
        k = W.kind[e]; ci = M.ci[e]; cj = M.cj[e]
        k == 2 && return diamond_dWSE_n_v(M.dtab, e, x[ci], x[cj], vt, M.L[e]) / max(M.L[e], 1.0)
        return W.aj[e] * (x[ci] - x[cj])               # kinds 1 and 3
    end
end

"""
    _face_coef_si!(e, opts, M, W, dt, α)

Phase A of the semi-implicit step for edge `e`: kind, q*, β, and what the
post-solve limiter needs (flow depth, q_t, donor depth). Mirrors
`_face_flux_std!`; `M.veta`/`M.vwet` hold the vertex values of ηⁿ.
"""
@inline function _face_coef_si!(e::Int, ::StdOpts{GC, FF, MOM, QC, VF, POS, FIXB}, M, W,
                                dt::Float64, α::Float64) where {GC, FF, MOM, QC, VF, POS, FIXB}
    @inbounds begin
        W.qs[e] = 0.0; W.beta[e] = 0.0; W.hf[e] = 0.0; W.kind[e] = Int8(0); W.aj[e] = 0.0
        ci = M.ci[e]; cj = M.cj[e]
        if (isnan(M.width[e]) || isnan(M.L[e]) || isnan(M.cos_theta[e]) || isnan(M.sill[e]) ||
            isnan(M.elevation[ci]) || isnan(M.elevation[cj]))
            return nothing
        end
        wse_ci = W.eta0[ci]; wse_cj = W.eta0[cj]
        h_raw = max(wse_ci, wse_cj) - M.sill[e]
        h_raw <= HFLOW_THRESHOLD && return nothing
        h = max(h_raw, 1e-6)
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
        jump = wse_ci - wse_cj
        Ls = max(M.L[e], 1.0)
        rem = 0.0
        if GC
            kind = Int8(3)
            if FF === :diamond || FF === :diamond_wetdry
                if FF === :diamond_wetdry && !diamond_stencil_wet_v(M.dtab, e, M.vwet)
                    kind = Int8(1)                               # two-point fallback
                elseif !isnan(diamond_dWSE_n_v(M.dtab, e, wse_ci, wse_cj, M.veta, M.L[e]))
                    kind = Int8(2)
                end
            end
            if kind == 3                                         # WLSQ: explicit correction
                gx_f = 0.5 * (M.grad_wse[1, ci] + M.grad_wse[1, cj])
                gy_f = 0.5 * (M.grad_wse[2, ci] + M.grad_wse[2, cj])
                rem = -M.alpha * M.L[e] * (gx_f * M.skew_x[e] + gy_f * M.skew_y[e]) / Ls
            end
            W.aj[e] = M.cos_theta[e] / Ls
        else
            kind = Int8(1)
            W.aj[e] = 1.0 / max(M.L[e] * max(M.cos_theta[e], 0.1), 1.0)
        end
        W.kind[e] = kind
        D0 = _si_D(e, W, M, W.eta0, M.veta)
        den = 1.0 + _G * h * dt * min(M.manning_n[ci], M.manning_n[cj])^2 *
                    _friction_q(q_prev_eff, q_t, VF) / h^(10.0 / 3.0)
        W.qs[e]   = (q_prev_eff - _G * h * dt * ((1.0 - α) * D0 + rem)) / den
        W.beta[e] = _G * h * dt * α / den
        W.hf[e] = h; W.qt[e] = q_t; W.dd[e] = depth_donor
    end
    return nothing
end

"""
    _face_post_si!(e, opts, M, W, dt)

After the solve: q⁺ = q* − β·D_f(ηⁿ⁺¹), Froude cap, volume limiter, stored
flux (Fix C) and edge volume after the donor cap — as `_face_flux_std!`.
"""
@inline function _face_post_si!(e::Int, ::StdOpts{GC, FF, MOM, QC, VF, POS, FIXB}, M, W,
                                dt::Float64) where {GC, FF, MOM, QC, VF, POS, FIXB}
    @inbounds begin
        h = W.hf[e]
        if h == 0.0
            M.flux[e] = 0.0; M.edge_vol[e] = 0.0
            return nothing
        end
        ci = M.ci[e]; cj = M.cj[e]
        q = W.qs[e] - W.beta[e] * _si_D(e, W, M, W.x, W.vtmp)
        q = _froude_cap(q, W.qt[e], h * sqrt(_G * h) * FROUDE_LIMIT, VF)
        if W.dd[e] > 0.0
            qv = (W.dd[e] * M.width[e]) / (M.max_nb * dt)
            q = clamp(q, -qv, qv)
        end
        M.flux[e] = q
        ev = q * M.width[e] * dt
        if !POS && !iszero(ev)
            ev = ev > 0.0 ? min(ev,  M.volume[cj] / M.donor_div) :
                            max(ev, -M.volume[ci] / M.donor_div)
        end
        M.edge_vol[e] = ev
    end
    return nothing
end

# y = (A + Δt·Σ s·w·β·D)(x), matrix-free (vertex pass if any diamond faces).
function _si_apply!(y, x, W::SIWork, M, cptr, cedge, area, dt, n, ne, use_d::Bool)
    use_d && diamond_vertex_eta!(W.vtmp, M.dtab, x)
    te = W.te; beta = W.beta; w = M.width; ci = M.ci
    @batch per=core for e in 1:ne
        @inbounds te[e] = beta[e] == 0.0 ? 0.0 : w[e] * beta[e] * _si_D(e, W, M, x, W.vtmp)
    end
    @batch per=core for i in 1:n
        s = 0.0
        @inbounds for k in cptr[i]:cptr[i + 1] - 1
            e = cedge[k]
            s = ci[e] == i ? s + te[e] : s - te[e]
        end
        @inbounds y[i] = area[i] * x[i] + dt * s
    end
    return y
end

_si_dot(a, b) = (s = 0.0; @inbounds for i in eachindex(a); s += a[i] * b[i]; end; s)

# Jacobi diagonal of the operator: A_i + Δt·Σ w·β·(∂D/∂x_i), using the
# two-point coefficient for every face (an approximation for diamond faces).
function _si_diag!(W, M, cptr, cedge, area, dt, n, ne)
    te = W.te
    @batch per=core for e in 1:ne
        @inbounds te[e] = M.width[e] * W.beta[e] * W.aj[e]
    end
    @batch per=core for i in 1:n
        s = 0.0
        @inbounds for k in cptr[i]:cptr[i + 1] - 1
            s += te[cedge[k]]
        end
        @inbounds W.dg[i] = area[i] + dt * s
    end
end

function _si_cg!(W, M, cptr, cedge, area, dt, n, ne; rtol, maxiter)
    x = W.x; r = W.r; z = W.z; p = W.p; Ap = W.Ap; dg = W.dg; b = W.b
    _si_apply!(Ap, x, W, M, cptr, cedge, area, dt, n, ne, false)
    @inbounds for i in 1:n; r[i] = b[i] - Ap[i]; z[i] = r[i] / dg[i]; p[i] = z[i]; end
    bn = sqrt(_si_dot(b, b)); bn == 0.0 && (bn = 1.0)
    rz = _si_dot(r, z); it = 0; rel = sqrt(_si_dot(r, r)) / bn
    while rel > rtol && it < maxiter
        _si_apply!(Ap, p, W, M, cptr, cedge, area, dt, n, ne, false)
        a = rz / _si_dot(p, Ap)
        @inbounds for i in 1:n; x[i] += a * p[i]; r[i] -= a * Ap[i]; z[i] = r[i] / dg[i]; end
        rz2 = _si_dot(r, z); bk = rz2 / rz; rz = rz2
        @inbounds for i in 1:n; p[i] = z[i] + bk * p[i]; end
        it += 1; rel = sqrt(_si_dot(r, r)) / bn
    end
    W.iters = it; W.resid = rel
end

function _si_bicgstab!(W, M, cptr, cedge, area, dt, n, ne; rtol, maxiter)
    x = W.x; r = W.r; rh = W.rh; p = W.p; v = W.vv; y = W.yy; s = W.ss; t = W.tt
    z = W.z; dg = W.dg; b = W.b
    _si_apply!(v, x, W, M, cptr, cedge, area, dt, n, ne, true)
    @inbounds for i in 1:n; r[i] = b[i] - v[i]; rh[i] = r[i]; p[i] = 0.0; v[i] = 0.0; end
    bn = sqrt(_si_dot(b, b)); bn == 0.0 && (bn = 1.0)
    ρ = 1.0; αk = 1.0; ω = 1.0; it = 0; rel = sqrt(_si_dot(r, r)) / bn
    while rel > rtol && it < maxiter
        ρn = _si_dot(rh, r)
        βk = (ρn / ρ) * (αk / ω); ρ = ρn
        @inbounds for i in 1:n; p[i] = r[i] + βk * (p[i] - ω * v[i]); y[i] = p[i] / dg[i]; end
        _si_apply!(v, y, W, M, cptr, cedge, area, dt, n, ne, true)
        αk = ρ / _si_dot(rh, v)
        @inbounds for i in 1:n; s[i] = r[i] - αk * v[i]; end
        if sqrt(_si_dot(s, s)) / bn <= rtol
            @inbounds for i in 1:n; x[i] += αk * y[i]; end
            it += 1; rel = sqrt(_si_dot(s, s)) / bn
            break
        end
        @inbounds for i in 1:n; z[i] = s[i] / dg[i]; end
        _si_apply!(t, z, W, M, cptr, cedge, area, dt, n, ne, true)
        ω = _si_dot(t, s) / _si_dot(t, t)
        @inbounds for i in 1:n; x[i] += αk * y[i] + ω * z[i]; r[i] = s[i] - ω * t[i]; end
        it += 1; rel = sqrt(_si_dot(r, r)) / bn
    end
    W.iters = it; W.resid = rel
end

"""
    _si_phase_a!(state, ws, opts, M, dt)

The semi-implicit replacement for Phase A of `step_standard!`: face
coefficients, the cell solve, and the post-solve face discharges and edge
volumes (written to `M.flux` and `M.edge_vol`, exactly where Phase A writes).
"""
function _si_phase_a!(state::FlowState, ws::StepWorkspace, opts::StdOpts, M, dt::Float64)
    n = length(state.cell_ids); ne = state.edges.n_edges
    nv = M.dtab === nothing ? 0 : length(M.dtab.vertices)
    W = _si_work!(ws.si, n, ne, nv)
    cptr, cedge = _cell_csr!(state, ws)
    α = state.fs_alpha
    area = state.cell_area
    @batch per=core for i in 1:n
        @inbounds W.eta0[i] = state.elevation[i] + state.volume[i] / max(area[i], 1.0)
    end
    @batch per=core for e in 1:ne
        _face_coef_si!(e, opts, M, W, dt, α)
    end
    # right-hand side: A·ηⁿ + Δt·Σ s·w·q*  (A as in the WSE definition, ≥ 1 m²)
    te = W.te; ci = M.ci
    @batch per=core for e in 1:ne
        @inbounds te[e] = M.width[e] * W.qs[e]
    end
    A1 = (a -> max(a, 1.0)).(area)
    @batch per=core for i in 1:n
        s = 0.0
        @inbounds for k in cptr[i]:cptr[i + 1] - 1
            e = cedge[k]
            s = ci[e] == i ? s + te[e] : s - te[e]
        end
        @inbounds W.b[i] = A1[i] * W.eta0[i] + dt * s
    end
    copyto!(W.x, W.eta0)
    _si_diag!(W, M, cptr, cedge, A1, dt, n, ne)
    use_d = any(==(Int8(2)), W.kind)
    if use_d
        _si_bicgstab!(W, M, cptr, cedge, A1, dt, n, ne; rtol = 1e-12, maxiter = 500)
        diamond_vertex_eta!(W.vtmp, M.dtab, W.x)        # vertex values of ηⁿ⁺¹
    else
        _si_cg!(W, M, cptr, cedge, A1, dt, n, ne; rtol = 1e-12, maxiter = 1000)
    end
    @batch per=core for e in 1:ne
        _face_post_si!(e, opts, M, W, dt)
    end
    return nothing
end
