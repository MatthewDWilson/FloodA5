# surfacewater/swe_reference.jl
# -----------------------------
# R1 (2026-10): full 2-D shallow-water reference solver on the polygon mesh,
# `--flow-model swe`. A benchmark for the reduced-complexity (local-inertial)
# solver, not a replacement for it: it resolves advection and has a
# dissipation set by the Riemann solver (wave speeds), not by the time step,
# so it is the "robust dissipation target" against which the local-inertial
# damping choices (HYDRAULICS.md §15.7–15.9) can be judged.
#
#   ∂h/∂t  + ∇·(h u)                       = sources
#   ∂hu/∂t + ∇·(h u⊗u + ½ g h² I)          = −g h ∇z − g n² |u| u / h^{1/3}
#                                            [+ ∇·(ν_t h ∇u)]
#
# Finite volumes on the cells (volume V = h·A and discharge per unit width
# q = h·u per cell; see SWEWork):
# - HLLC fluxes in each face's normal frame (Toro). HLLC rather than HLL:
#   HLL smears the tangential velocity across the contact wave, a numerical
#   viscosity on exactly the lateral momentum exchange this solver is meant
#   to measure.
# - Hydrostatic reconstruction (Audusse et al. 2004) at the face bed
#   z_f = max(sill, z_i, z_j): positivity-preserving with wet/dry fronts.
#   Each side's momentum flux is taken relative to its own cell's
#   hydrostatic pressure ½g·h_i²·n (a term whose sum over the cell's faces,
#   ½g·h_i²·Σ w·n, is zero for a closed planar polygon). On the sphere the
#   face lengths and normals of a cell do not close exactly, and without
#   this a lake at rest feels the closure error (10⁻⁴ m drift on A5); with
#   it the lake at rest balances face by face, whatever the geometry.
# - Order 1: first-order states, forward Euler. Order 2: MUSCL — WLSQ
#   gradients of η, u, v (the solver's wlsq_weights), Barth–Jespersen
#   limiter, first order in cells next to a dry one — with SSP-RK2 (Heun).
#   η (not h) is reconstructed, which keeps a lake at rest exact; but under
#   a thin film on a stepped bed η follows the bed, so its gradient is
#   additionally limited to keep every face depth within [0, 2·h_i]
#   (|∇η·r| ≤ h_i): deep water keeps second order, thin films fall back to
#   first order continuously.
# - Manning friction point-implicit in each stage.
# - Optional eddy viscosity ν_t = C_E·h·u* (u* = √g·n·|u|/h^{1/6}),
#   `--swe-viscosity C_E`, as a face-based momentum diffusion.
# - Boundaries: Closed → reflecting wall; ZeroGradient → transmissive for
#   outflow and a wall when the flow at the face points inwards (FloodA5's
#   open boundary is outflow-only, as in the standard solver, which clamps
#   the boundary volume change at ≥ 0; a transmissive ghost would otherwise
#   copy an inward velocity and draw water from an unlimited reservoir;
#   `swe_open_inflow = true` makes them fully transmissive, the textbook
#   zero-gradient condition, for benchmarks such as uniform channel flow);
#   FixedWSE/FixedQ, not yet implemented, fall back to it, as in the
#   standard solver; Critical → outflow at least at critical speed.
#   Outflow volume goes to state.vol_removed, so the mass balance closes.
# - Time step (`_cfl_dt`): Δt = C·min_i A_i / Σ_f w_f·s_f (s_f the faster of
#   the face's two cells), the per-cell sufficient condition for positivity
#   of the first-order scheme (C ≤ 1); `--courant` sets C. Any volume added
#   by clamping a negative cell is counted (SWEWork.clamp_vol) and reported.
# Data layout and the GPU-port rules: surfacewater/swe_types.jl.

const SWE_HDRY = 1e-6      # m: below this a cell is dry (no momentum)
const SWE_HVEL = 1e-3      # m: velocity desingularisation u = q·h / max(h², H²)

@inline _swe_vel(q, h) = q * h / max(h * h, SWE_HVEL * SWE_HVEL)

# ── Geometry ────────────────────────────────────────────────────────────────

function _swe_bc_code(bc::BCType)::Int8
    bc === Closed   && return SWE_BC_WALL
    bc === Critical && return SWE_BC_CRITICAL
    return SWE_BC_TRANSMISSIVE
end

"""
    _build_swe_work(mesh, state; order = 2, visc_ce = 0.0) → SWEWork

Face normals and midpoints from the cell polygons (the true geometric face
normal, which on A5 differs from the centre-to-centre direction), boundary
faces from `state.ghost_edges`, the cell → boundary-face CSR, and the work
arrays.
"""
function _build_swe_work(mesh::A5Mesh, state::FlowState; order::Int = 2, visc_ce::Float64 = 0.0,
                         open_inflow::Bool = false)
    order in (1, 2) || throw(ArgumentError("swe order must be 1 or 2 (got $order)"))
    (isfinite(visc_ce) && visc_ce >= 0.0) || throw(ArgumentError("swe viscosity C_E must be ≥ 0"))
    E = state.edges; ne = E.n_edges; n = length(state.cell_ids); R = A5Grid._EARTH_R
    lon = state.cell_lons; lat = state.cell_lats
    nx = zeros(ne); ny = zeros(ne); rxi = zeros(ne); ryi = zeros(ne); rxj = zeros(ne); ryj = zeros(ne)
    loc(lo, la, i) = (deg2rad(lo - lon[i]) * R * cosd(lat[i]), deg2rad(la - lat[i]) * R)
    nfallback = 0
    for e in 1:ne
        i = E.cell_i[e]; j = E.cell_j[e]
        dx = E.dx_m[e]; dy = E.dy_m[e]; dl = hypot(dx, dy)
        sh = A5Grid._shared_edge(mesh.cells[i].boundary, mesh.cells[j].boundary)
        degenerate = sh === nothing || !(dl > 0.0)
        if !degenerate
            lo1, la1, lo2, la2 = sh
            degenerate = !(hypot(deg2rad(lo2 - lo1) * R * cosd(la1), deg2rad(la2 - la1) * R) > 1e-3)
        end
        if degenerate                                  # no (or a vanishing) shared side
            nfallback += 1
            nx[e] = dl > 0.0 ? dx / dl : 1.0; ny[e] = dl > 0.0 ? dy / dl : 0.0
            mlo = 0.5 * (lon[i] + lon[j]); mla = 0.5 * (lat[i] + lat[j])
        else
            lo1, la1, lo2, la2 = sh
            mlo = 0.5 * (lo1 + lo2); mla = 0.5 * (la1 + la2)
            tx = deg2rad(lo2 - lo1) * R * cosd(mla); ty = deg2rad(la2 - la1) * R
            tl = hypot(tx, ty)
            a, b = ty / tl, -tx / tl
            a * dx + b * dy < 0.0 && (a = -a; b = -b)           # orient cell_i → cell_j
            nx[e] = a; ny[e] = b
        end
        rxi[e], ryi[e] = loc(mlo, mla, i)
        rxj[e], ryj[e] = loc(mlo, mla, j)
    end
    nfallback > 0 && @info "SWE: $nfallback faces without a usable shared polygon side — centre-line normal used"
    ng = length(state.ghost_edges)
    gcell = [g.cell_index for g in state.ghost_edges]
    gnx = [g.nx for g in state.ghost_edges]; gny = [g.ny for g in state.ghost_edges]
    gw = [g.width for g in state.ghost_edges]; gz = [g.sill for g in state.ghost_edges]
    grx = [g.rx for g in state.ghost_edges]; gry = [g.ry for g in state.ghost_edges]
    gbc = Int8[_swe_bc_code(b) for b in state.ghost_cell_bc]
    cnt = zeros(Int, n)
    for c in gcell; cnt[c] += 1; end
    gptr = Vector{Int}(undef, n + 1); gptr[1] = 1
    for i in 1:n; gptr[i + 1] = gptr[i] + cnt[i]; end
    gface = zeros(Int, ng); fill!(cnt, 0)
    for (k, c) in enumerate(gcell); gface[gptr[c] + cnt[c]] = k; cnt[c] += 1; end
    perim = zeros(n)
    for e in 1:ne; perim[E.cell_i[e]] += E.width[e]; perim[E.cell_j[e]] += E.width[e]; end
    for k in 1:ng; perim[gcell[k]] += gw[k]; end
    z(k) = zeros(k)
    return SWEWork(order, visc_ce, open_inflow, nx, ny, rxi, ryi, rxj, ryj,
                   gcell, gnx, gny, gw, gz, grx, gry, gbc, gptr, gface, perim,
                   z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n),
                   z(ne), z(ne), z(ne), z(ne), z(ne), z(ne), z(ne),
                   z(ng), z(ng), z(ng), z(n), z(n), z(n), z(n))
end

# ── HLLC flux in the face frame ─────────────────────────────────────────────

"""
    _swe_hllc(hL, unL, utL, hR, unR, utR) → (mass, normal momentum, tangential momentum)

HLLC flux per unit face length for the shallow-water equations in the face
frame (Toro): two-rarefaction wave-speed estimates with the dry-bed limits,
HLL for mass and normal momentum, the contact wave carrying the tangential
velocity upwind.
"""
@inline function _swe_hllc(hL::Float64, unL::Float64, utL::Float64,
                           hR::Float64, unR::Float64, utR::Float64)
    wetL = hL > SWE_HDRY; wetR = hR > SWE_HDRY
    (wetL || wetR) || return (0.0, 0.0, 0.0)
    cL = sqrt(_G * hL); cR = sqrt(_G * hR)
    if !wetL
        sL = unR - 2.0 * cR; sR = unR + cR
    elseif !wetR
        sL = unL - cL; sR = unL + 2.0 * cL
    else
        us = 0.5 * (unL + unR) + cL - cR
        cs = 0.5 * (cL + cR) + 0.25 * (unL - unR)
        sL = min(unL - cL, us - cs); sR = max(unR + cR, us + cs)
    end
    fL1 = hL * unL; fL2 = hL * unL * unL + 0.5 * _G * hL * hL
    fR1 = hR * unR; fR2 = hR * unR * unR + 0.5 * _G * hR * hR
    sL >= 0.0 && return (fL1, fL2, fL1 * utL)
    sR <= 0.0 && return (fR1, fR2, fR1 * utR)
    inv = 1.0 / (sR - sL)
    f1 = (sR * fL1 - sL * fR1 + sL * sR * (hR - hL)) * inv
    f2 = (sR * fL2 - sL * fR2 + sL * sR * (hR * unR - hL * unL)) * inv
    den = hR * (unR - sR) - hL * (unL - sL)
    ss = den != 0.0 ? (sL * hR * (unR - sR) - sR * hL * (unL - sL)) / den : 0.0
    return (f1, f2, f1 * (ss >= 0.0 ? utL : utR))
end

# ── Stage passes ────────────────────────────────────────────────────────────

# primitives per cell: h, η, desingularised u, v; dry or NaN-bed cells inert
function _swe_primitives!(sw::SWEWork, V, qx, qy, z, A)
    h = sw.h; u = sw.u; v = sw.v; eta = sw.eta
    @batch per=core for i in eachindex(h)
        @inbounds begin
            if isfinite(z[i]) && A[i] >= 1.0
                hi = max(V[i], 0.0) / A[i]
                h[i] = hi; eta[i] = z[i] + hi
                if hi > SWE_HDRY
                    u[i] = _swe_vel(qx[i], hi); v[i] = _swe_vel(qy[i], hi)
                else
                    u[i] = 0.0; v[i] = 0.0
                end
            else
                h[i] = 0.0; eta[i] = NaN; u[i] = 0.0; v[i] = 0.0
            end
        end
    end
    return nothing
end

# Barth–Jespersen factor for one variable: the reconstruction at every face
# midpoint of cell i stays within the range of cell i and its neighbours
@inline function _bj(φi, φmin, φmax, gx, gy, i, cptr, cedge, ci, rxi, ryi, rxj, ryj, gptr, gface, grx, gry)
    α = 1.0
    @inbounds begin
        for k in cptr[i]:cptr[i + 1] - 1
            e = cedge[k]
            rx = ci[e] == i ? rxi[e] : rxj[e]; ry = ci[e] == i ? ryi[e] : ryj[e]
            d = gx * rx + gy * ry
            d > 0.0 ? (α = min(α, (φmax - φi) / d)) : d < 0.0 ? (α = min(α, (φmin - φi) / d)) : nothing
        end
        for k in gptr[i]:gptr[i + 1] - 1
            f = gface[k]; d = gx * grx[f] + gy * gry[f]
            d > 0.0 ? (α = min(α, (φmax - φi) / d)) : d < 0.0 ? (α = min(α, (φmin - φi) / d)) : nothing
        end
    end
    return max(α, 0.0)
end

# limited WLSQ gradients of η, u, v (order 2); zero next to dry cells
function _swe_gradients!(sw::SWEWork, state::FlowState, cptr, cedge)
    h = sw.h; u = sw.u; v = sw.v; eta = sw.eta
    adj = state.adj_matrix; W = state.wlsq_weights; nb = size(adj, 1)
    ci = state.edges.cell_i
    rxi = sw.rxi; ryi = sw.ryi; rxj = sw.rxj; ryj = sw.ryj
    gptr = sw.gptr; gface = sw.gface; grx = sw.grx; gry = sw.gry
    dex = sw.dex; dey = sw.dey; dux = sw.dux; duy = sw.duy; dvx = sw.dvx; dvy = sw.dvy
    @batch per=core for i in eachindex(h)
        @inbounds begin
            ok = h[i] > SWE_HDRY
            ex = 0.0; ey = 0.0; ux = 0.0; uy = 0.0; vx = 0.0; vy = 0.0
            emin = eta[i]; emax = eta[i]; umin = u[i]; umax = u[i]; vmin = v[i]; vmax = v[i]
            if ok
                for s in 1:nb
                    j = adj[s, i]; j == 0 && continue
                    if !(h[j] > SWE_HDRY)                    # dry (or NaN-bed) neighbour
                        ok = false; break
                    end
                    wx = W[s, i]; wy = W[nb + s, i]
                    de = eta[j] - eta[i]; du = u[j] - u[i]; dv = v[j] - v[i]
                    ex += wx * de; ey += wy * de; ux += wx * du; uy += wy * du; vx += wx * dv; vy += wy * dv
                    emin = min(emin, eta[j]); emax = max(emax, eta[j])
                    umin = min(umin, u[j]); umax = max(umax, u[j])
                    vmin = min(vmin, v[j]); vmax = max(vmax, v[j])
                end
            end
            if ok
                ae = _bj(eta[i], emin, emax, ex, ey, i, cptr, cedge, ci, rxi, ryi, rxj, ryj, gptr, gface, grx, gry)
                # depth bound: face depth h_i + ae·∇η·r within [0, 2h_i]
                dm = 0.0
                for k in cptr[i]:cptr[i + 1] - 1
                    e = cedge[k]
                    rx = ci[e] == i ? rxi[e] : rxj[e]; ry = ci[e] == i ? ryi[e] : ryj[e]
                    dm = max(dm, abs(ex * rx + ey * ry))
                end
                for k in gptr[i]:gptr[i + 1] - 1
                    f = gface[k]; dm = max(dm, abs(ex * grx[f] + ey * gry[f]))
                end
                ae * dm > h[i] && (ae = h[i] / dm)
                au = _bj(u[i], umin, umax, ux, uy, i, cptr, cedge, ci, rxi, ryi, rxj, ryj, gptr, gface, grx, gry)
                av = _bj(v[i], vmin, vmax, vx, vy, i, cptr, cedge, ci, rxi, ryi, rxj, ryj, gptr, gface, grx, gry)
                dex[i] = ae * ex; dey[i] = ae * ey
                dux[i] = au * ux; duy[i] = au * uy
                dvx[i] = av * vx; dvy[i] = av * vy
            else
                dex[i] = 0.0; dey[i] = 0.0; dux[i] = 0.0; duy[i] = 0.0; dvx[i] = 0.0; dvy[i] = 0.0
            end
        end
    end
    return nothing
end

# interior faces: reconstructed states, hydrostatic reconstruction, HLLC,
# pressure corrections, optional viscous stress
function _swe_face_fluxes!(sw::SWEWork, state::FlowState)
    E = state.edges; z = state.elevation; nman = state.manning_n
    ci = E.cell_i; cj = E.cell_j; sill = E.sill; L = E.L
    h = sw.h; u = sw.u; v = sw.v; eta = sw.eta
    dex = sw.dex; dey = sw.dey; dux = sw.dux; duy = sw.duy; dvx = sw.dvx; dvy = sw.dvy
    nx = sw.nx; ny = sw.ny; rxi = sw.rxi; ryi = sw.ryi; rxj = sw.rxj; ryj = sw.ryj
    F1 = sw.F1; F2 = sw.F2; F3 = sw.F3; PL = sw.PL; PR = sw.PR; T2 = sw.T2; T3 = sw.T3
    ce = sw.visc_ce
    @batch per=core for e in 1:E.n_edges
        @inbounds begin
            i = ci[e]; j = cj[e]
            if !(h[i] > SWE_HDRY || h[j] > SWE_HDRY) || !isfinite(eta[i]) || !isfinite(eta[j])
                F1[e] = 0.0; F2[e] = 0.0; F3[e] = 0.0; PL[e] = 0.0; PR[e] = 0.0; T2[e] = 0.0; T3[e] = 0.0
            else
                ax = rxi[e]; ay = ryi[e]; bx = rxj[e]; by = ryj[e]
                ηL = eta[i] + dex[i] * ax + dey[i] * ay; ηR = eta[j] + dex[j] * bx + dey[j] * by
                uL = u[i] + dux[i] * ax + duy[i] * ay;   uR = u[j] + dux[j] * bx + duy[j] * by
                vL = v[i] + dvx[i] * ax + dvy[i] * ay;   vR = v[j] + dvx[j] * bx + dvy[j] * by
                hL = h[i] > SWE_HDRY ? max(ηL - z[i], 0.0) : 0.0
                hR = h[j] > SWE_HDRY ? max(ηR - z[j], 0.0) : 0.0
                zf = max(sill[e], z[i], z[j])
                hLs = max(hL + z[i] - zf, 0.0); hRs = max(hR + z[j] - zf, 0.0)
                n1 = nx[e]; n2 = ny[e]
                unL = uL * n1 + vL * n2; utL = -uL * n2 + vL * n1
                unR = uR * n1 + vR * n2; utR = -uR * n2 + vR * n1
                f1, fn, ft = _swe_hllc(hLs, unL, utL, hRs, unR, utR)
                F1[e] = f1; F2[e] = fn * n1 - ft * n2; F3[e] = fn * n2 + ft * n1
                PL[e] = 0.5 * _G * (hL * hL - hLs * hLs - h[i] * h[i])   # relative to cell i's pressure
                PR[e] = 0.5 * _G * (hR * hR - hRs * hRs - h[j] * h[j])
                if ce > 0.0 && h[i] > SWE_HDRY && h[j] > SWE_HDRY
                    hf = min(h[i], h[j])
                    sp = 0.5 * (hypot(u[i], v[i]) + hypot(u[j], v[j]))
                    ustar = sqrt(_G) * 0.5 * (nman[i] + nman[j]) * sp / hf^(1.0 / 6.0)
                    k = ce * hf * ustar * hf / L[e]                 # ν·h/L
                    T2[e] = k * (u[j] - u[i]); T3[e] = k * (v[j] - v[i])
                else
                    T2[e] = 0.0; T3[e] = 0.0
                end
            end
        end
    end
    return nothing
end

# boundary faces: the ghost state from the BC code, fluxes out of the cell
function _swe_boundary_fluxes!(sw::SWEWork, state::FlowState)
    z = state.elevation
    h = sw.h; u = sw.u; v = sw.v; eta = sw.eta
    dex = sw.dex; dey = sw.dey; dux = sw.dux; duy = sw.duy; dvx = sw.dvx; dvy = sw.dvy
    gcell = sw.gcell; gnx = sw.gnx; gny = sw.gny; gz = sw.gz; grx = sw.grx; gry = sw.gry; gbc = sw.gbc
    G1 = sw.G1; G2 = sw.G2; G3 = sw.G3; open_in = sw.open_inflow
    @batch per=core for k in eachindex(gcell)
        @inbounds begin
            c = gcell[k]
            if !(h[c] > SWE_HDRY) || !isfinite(eta[c])
                G1[k] = 0.0; G2[k] = 0.0; G3[k] = 0.0
            else
                ax = grx[k]; ay = gry[k]; n1 = gnx[k]; n2 = gny[k]
                ηL = eta[c] + dex[c] * ax + dey[c] * ay
                uL = u[c] + dux[c] * ax + duy[c] * ay; vL = v[c] + dvx[c] * ax + dvy[c] * ay
                hL = max(ηL - z[c], 0.0)
                zf = max(gz[k], z[c]); hs = max(hL + z[c] - zf, 0.0)
                unL = uL * n1 + vL * n2; utL = -uL * n2 + vL * n1
                b = gbc[k]
                if b == SWE_BC_WALL || (b == SWE_BC_TRANSMISSIVE && unL <= 0.0 && !open_in)
                    f1, fn, ft = _swe_hllc(hs, unL, utL, hs, -unL, utL); f1 = 0.0; ft = 0.0
                elseif b == SWE_BC_CRITICAL
                    f1, fn, ft = _swe_hllc(hs, unL, utL, hs, max(unL, sqrt(_G * hs)), utL)
                else
                    f1, fn, ft = _swe_hllc(hs, unL, utL, hs, unL, utL)
                end
                P = 0.5 * _G * (hL * hL - hs * hs - h[c] * h[c])
                G1[k] = f1; G2[k] = (fn + P) * n1 - ft * n2; G3[k] = (fn + P) * n2 + ft * n1
            end
        end
    end
    return nothing
end

# per-cell gather and update, point-implicit Manning friction
function _swe_update!(sw::SWEWork, state::FlowState, dt::Float64, cptr, cedge)
    E = state.edges; ci = E.cell_i; w = E.width
    V = state.volume; qx = state.qvec_u; qy = state.qvec_v; A = state.cell_area; z = state.elevation
    nman = state.manning_n
    F1 = sw.F1; F2 = sw.F2; F3 = sw.F3; PL = sw.PL; PR = sw.PR; T2 = sw.T2; T3 = sw.T3
    nx = sw.nx; ny = sw.ny
    G1 = sw.G1; G2 = sw.G2; G3 = sw.G3; gw = sw.gw; gptr = sw.gptr; gface = sw.gface
    clampv = sw.clamp_vol
    @batch per=core for i in eachindex(V)
        @inbounds begin
            (isfinite(z[i]) && A[i] >= 1.0) || continue
            dm = 0.0; dx = 0.0; dy = 0.0
            for k in cptr[i]:cptr[i + 1] - 1
                e = cedge[k]; we = w[e]
                if ci[e] == i                         # flux leaves i along n
                    dm -= we * F1[e]
                    dx -= we * (F2[e] + PL[e] * nx[e] - T2[e])
                    dy -= we * (F3[e] + PL[e] * ny[e] - T3[e])
                else                                  # enters j
                    dm += we * F1[e]
                    dx += we * (F2[e] + PR[e] * nx[e] - T2[e])
                    dy += we * (F3[e] + PR[e] * ny[e] - T3[e])
                end
            end
            for k in gptr[i]:gptr[i + 1] - 1
                f = gface[k]; wf = gw[f]
                dm -= wf * G1[f]; dx -= wf * G2[f]; dy -= wf * G3[f]
            end
            Vn = V[i] + dt * dm
            if Vn < 0.0
                clampv[i] -= Vn; Vn = 0.0
            end
            hn = Vn / A[i]
            qxn = qx[i] + dt * dx / A[i]; qyn = qy[i] + dt * dy / A[i]
            if hn > SWE_HDRY
                qm = hypot(qxn, qyn)
                fr = 1.0 + dt * _G * nman[i]^2 * qm / hn^(7.0 / 3.0)
                qxn /= fr; qyn /= fr
            else
                qxn = 0.0; qyn = 0.0
            end
            V[i] = Vn; qx[i] = qxn; qy[i] = qyn
        end
    end
    return nothing
end

# one forward-Euler stage U ← U + Δt·L(U); returns the boundary outflow volume
function _swe_stage!(state::FlowState, sw::SWEWork, dt::Float64, cptr, cedge)
    _swe_primitives!(sw, state.volume, state.qvec_u, state.qvec_v, state.elevation, state.cell_area)
    if sw.order == 2
        _swe_gradients!(sw, state, cptr, cedge)
    end
    _swe_face_fluxes!(sw, state)
    _swe_boundary_fluxes!(sw, state)
    out = 0.0
    @inbounds for k in eachindex(sw.G1); out += dt * sw.gw[k] * sw.G1[k]; end
    _swe_update!(sw, state, dt, cptr, cedge)
    return out
end

"""
    step_swe!(state, dt; courant = 0.7) → number of sub-steps

Advance the R1 shallow-water reference solver by `dt` (see the file header).
The run loop picks `dt` before the sources are applied, so a source can put
water into a dry cell that the step was never sized for (an inflow of
100 m³/s puts metres of water into one cell). The step therefore re-checks
the time-step limit on the current state and sub-steps as needed: `dt` is a
target, never a stability constraint. Order 1: one Euler stage per
sub-step. Order 2: SSP-RK2, U = ½(Uⁿ + (U⁽¹⁾ + Δt·L(U⁽¹⁾))).
"""
function step_swe!(state::FlowState, dt::Float64; courant::Float64 = 0.7)
    sw = state.swe
    sw === nothing && throw(ArgumentError("step_swe! needs a state initialised with SWEFlow()"))
    ws = _workspace!(state)
    cptr, cedge = _cell_csr!(state, ws)
    @inbounds for k in eachindex(sw.gbc); sw.gbc[k] = _swe_bc_code(state.ghost_cell_bc[k]); end
    left = dt; nsub = 0
    while left > 0.0
        h = min(left, _cfl_dt(state, SWEFlow(); courant = courant))
        h >= left * (1 - 1e-9) && (h = left)                # no sliver sub-step at the end
        _swe_substep!(state, sw, h, cptr, cedge)
        left -= h; nsub += 1
    end
    A = state.cell_area; V = state.volume; d = state.water_depth
    @batch per=core for i in eachindex(V)
        @inbounds (A[i] >= 1.0) && (d[i] = V[i] / A[i])
    end
    return nsub
end

function _swe_substep!(state::FlowState, sw::SWEWork, dt::Float64, cptr, cedge)
    if sw.order == 1
        state.vol_removed += _swe_stage!(state, sw, dt, cptr, cedge)
    else
        copyto!(sw.V0, state.volume); copyto!(sw.qx0, state.qvec_u); copyto!(sw.qy0, state.qvec_v)
        o1 = _swe_stage!(state, sw, dt, cptr, cedge)
        o2 = _swe_stage!(state, sw, dt, cptr, cedge)
        V = state.volume; qx = state.qvec_u; qy = state.qvec_v
        V0 = sw.V0; qx0 = sw.qx0; qy0 = sw.qy0
        @batch per=core for i in eachindex(V)
            @inbounds begin
                V[i] = 0.5 * (V0[i] + V[i]); qx[i] = 0.5 * (qx0[i] + qx[i]); qy[i] = 0.5 * (qy0[i] + qy[i])
            end
        end
        state.vol_removed += 0.5 * (o1 + o2)
    end
    return nothing
end

"""
    _swe_velocity!(state)

Cell velocity for output: u = q/h (desingularised), its magnitude.
"""
function _swe_velocity!(state::FlowState)
    @inbounds for i in eachindex(state.volume)
        h = state.water_depth[i]
        if h > SWE_HDRY
            state.vel_u[i] = _swe_vel(state.qvec_u[i], h); state.vel_v[i] = _swe_vel(state.qvec_v[i], h)
        else
            state.vel_u[i] = 0.0; state.vel_v[i] = 0.0
        end
        state.velocity[i] = hypot(state.vel_u[i], state.vel_v[i])
    end
    return nothing
end

"""
    _cfl_dt(state, ::SWEFlow; courant)

Δt = C·min over cells of A_i / Σ_f w_f·s_f, with s_f the larger of the two
adjacent cells' |u| + √(g h) (boundary faces: the cell's own): the
first-order positivity condition (C ≤ 1) — a shallow cell next to a deep
one can lose water at the deep cell's wave speed. Plus the explicit limit of
the optional eddy viscosity, 4·ν_i·P_i/ℓ_i with ℓ_i = 2A_i/P_i.
"""
function _cfl_dt(state::FlowState, ::SWEFlow; courant::Float64 = 0.7)::Float64
    sw = state.swe
    ws = _workspace!(state); cptr, cedge = _cell_csr!(state, ws)
    E = state.edges; ci = E.cell_i; cj = E.cell_j; w = E.width
    sp = sw.u                                   # scratch: per-cell signal speed
    @inbounds for i in eachindex(state.volume)
        A = state.cell_area[i]; h = A >= 1.0 && isfinite(state.elevation[i]) ? state.volume[i] / A : 0.0
        sp[i] = h > SWE_HDRY ? hypot(_swe_vel(state.qvec_u[i], h), _swe_vel(state.qvec_v[i], h)) + sqrt(_G * h) : 0.0
    end
    best = Inf
    @inbounds for i in eachindex(state.volume)
        A = state.cell_area[i]
        (A >= 1.0 && isfinite(state.elevation[i])) || continue
        rate = 0.0
        for k in cptr[i]:cptr[i + 1] - 1
            e = cedge[k]; rate += w[e] * max(sp[ci[e]], sp[cj[e]])
        end
        for k in sw.gptr[i]:sw.gptr[i + 1] - 1
            rate += sw.gw[sw.gface[k]] * sp[i]
        end
        if sw.visc_ce > 0.0
            h = state.volume[i] / A
            if h > SWE_HDRY
                u = hypot(_swe_vel(state.qvec_u[i], h), _swe_vel(state.qvec_v[i], h))
                ν = sw.visc_ce * h * sqrt(_G) * state.manning_n[i] * u / h^(1.0 / 6.0)
                rate += 4.0 * ν * sw.perim[i] / (2A / sw.perim[i])
            end
        end
        rate > 0.0 && (best = min(best, A / rate))
    end
    return isfinite(best) ? courant * best : 60.0
end
