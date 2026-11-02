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
#   z_f = max(z_iL, z_jR), with each side's face bed z_iL = η_iL − h_iL
#   from the reconstructed states (z_i at first order) and a centred source
#   ½g·(h_iL + h_i)·(z_iL − z_i) per face — the second-order scheme of
#   Audusse et al. (2004, §4). Positivity-preserving with wet/dry fronts;
#   an edge sill above both cell beds (an obstacle) also bounds z_f.
#   Reconstructing the face bed matters on thin films over a slope: with
#   cell-centre beds every face removes the bed step Δz from the
#   downstream depth, which on a rained 0.1% slope at 22 m cells made the
#   sheet flow 6–11 mm (10–50%) too deep (HYDRAULICS.md §15.10).
#   Each side's momentum flux is taken relative to its own cell's
#   hydrostatic pressure ½g·h_i²·n (a term whose sum over the cell's faces,
#   ½g·h_i²·Σ w·n, is zero for a closed planar polygon). On the sphere the
#   face lengths and normals of a cell do not close exactly, and without
#   this a lake at rest feels the closure error (10⁻⁴ m drift on A5); with
#   it the lake at rest balances face by face, whatever the geometry.
# - Order 1: first-order states, forward Euler. Order 2: MUSCL — WLSQ
#   gradients of η, h, u, v (the solver's wlsq_weights), Barth–Jespersen
#   limiter, first order in cells next to a dry one — with SSP-RK2 (Heun).
#   η and h are both reconstructed: η keeps a lake at rest exact, and the
#   face depth h_iL comes from h's own limited reconstruction, so it stays
#   within the neighbours' depths (non-negative) even where η follows a
#   stepped bed under a thin film.
# - Manning friction semi-implicit in each stage, with |q| from the stage's
#   start: q ← (q + Δt·F)/(1 + Δt·g·n²·|q₀|/h^{7/3}), whose fixed point is
#   exactly Manning's normal flow. (Taking |q| from the predicted q + Δt·F
#   moves the fixed point to q² + Δt·F·q = q_Manning²; on thin films, where
#   Δt·g·h·S is as large as q itself, that slowed sheet flow by up to 40%.)
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
                   z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n), z(n),
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

# ── Arrays ──────────────────────────────────────────────────────────────────

"""
    _swe_args(state) → NamedTuple

Everything the SWE passes read or write: host arrays (references, not
copies) and scalars. `to_device` adapts it to a GPU (WP5, surfacewater/backend.jl).
"""
function _swe_args(state::FlowState)
    sw = state.swe; E = state.edges
    ws = _workspace!(state); cptr, cedge = _cell_csr!(state, ws)
    return (V = state.volume, qx = state.qvec_u, qy = state.qvec_v, z = state.elevation,
            A = state.cell_area, nman = state.manning_n, d = state.water_depth,
            vu = state.vel_u, vv = state.vel_v, vel = state.velocity,
            ci = E.cell_i, cj = E.cell_j, w = E.width, L = E.L, sill = E.sill,
            adj = state.adj_matrix, W = state.wlsq_weights, nb = size(state.adj_matrix, 1),
            cptr = cptr, cedge = cedge,
            nx = sw.nx, ny = sw.ny, rxi = sw.rxi, ryi = sw.ryi, rxj = sw.rxj, ryj = sw.ryj,
            gcell = sw.gcell, gnx = sw.gnx, gny = sw.gny, gw = sw.gw, gz = sw.gz,
            grx = sw.grx, gry = sw.gry, gbc = sw.gbc, gptr = sw.gptr, gface = sw.gface,
            perim = sw.perim,
            h = sw.h, u = sw.u, v = sw.v, eta = sw.eta, dex = sw.dex, dey = sw.dey,
            dux = sw.dux, duy = sw.duy, dvx = sw.dvx, dvy = sw.dvy, dhx = sw.dhx, dhy = sw.dhy,
            F1 = sw.F1, F2 = sw.F2, F3 = sw.F3, PL = sw.PL, PR = sw.PR, T2 = sw.T2, T3 = sw.T3,
            G1 = sw.G1, G2 = sw.G2, G3 = sw.G3, V0 = sw.V0, qx0 = sw.qx0, qy0 = sw.qy0,
            clamp = sw.clamp_vol,
            order = sw.order, ce = sw.visc_ce, open_in = sw.open_inflow,
            n = length(state.volume), ne = E.n_edges, ng = length(sw.gcell))
end

# ── Stage passes (kernel bodies: one cell, face or boundary face each) ─────

# primitives per cell: h, η, desingularised u, v; dry or NaN-bed cells inert
@inline function _swe_prim_k!(i, a)
    @inbounds begin
        if isfinite(a.z[i]) && a.A[i] >= 1.0
            hi = max(a.V[i], 0.0) / a.A[i]
            a.h[i] = hi; a.eta[i] = a.z[i] + hi
            if hi > SWE_HDRY
                a.u[i] = _swe_vel(a.qx[i], hi); a.v[i] = _swe_vel(a.qy[i], hi)
            else
                a.u[i] = 0.0; a.v[i] = 0.0
            end
        else
            a.h[i] = 0.0; a.eta[i] = NaN; a.u[i] = 0.0; a.v[i] = 0.0
        end
    end
    return nothing
end

# Barth–Jespersen factor for one variable: the reconstruction at every face
# midpoint of cell i stays within the range of cell i and its neighbours
@inline function _bj(φi, φmin, φmax, gx, gy, i, a)
    α = 1.0
    @inbounds begin
        for k in a.cptr[i]:a.cptr[i + 1] - 1
            e = a.cedge[k]
            rx = a.ci[e] == i ? a.rxi[e] : a.rxj[e]; ry = a.ci[e] == i ? a.ryi[e] : a.ryj[e]
            d = gx * rx + gy * ry
            d > 0.0 ? (α = min(α, (φmax - φi) / d)) : d < 0.0 ? (α = min(α, (φmin - φi) / d)) : nothing
        end
        for k in a.gptr[i]:a.gptr[i + 1] - 1
            f = a.gface[k]; d = gx * a.grx[f] + gy * a.gry[f]
            d > 0.0 ? (α = min(α, (φmax - φi) / d)) : d < 0.0 ? (α = min(α, (φmin - φi) / d)) : nothing
        end
    end
    return max(α, 0.0)
end

# limited WLSQ gradients of η, h, u, v (order 2); zero next to dry cells
@inline function _swe_grad_k!(i, a)
    @inbounds begin
        h = a.h; u = a.u; v = a.v; eta = a.eta; adj = a.adj; W = a.W; nb = a.nb
        ok = h[i] > SWE_HDRY
        ex = 0.0; ey = 0.0; hx = 0.0; hy = 0.0; ux = 0.0; uy = 0.0; vx = 0.0; vy = 0.0
        emin = eta[i]; emax = eta[i]; hmin = h[i]; hmax = h[i]
        umin = u[i]; umax = u[i]; vmin = v[i]; vmax = v[i]
        if ok
            for s in 1:nb
                j = adj[s, i]; j == 0 && continue
                if !(h[j] > SWE_HDRY)                    # dry (or NaN-bed) neighbour
                    ok = false; break
                end
                wx = W[s, i]; wy = W[nb + s, i]
                de = eta[j] - eta[i]; dh = h[j] - h[i]; du = u[j] - u[i]; dv = v[j] - v[i]
                ex += wx * de; ey += wy * de; hx += wx * dh; hy += wy * dh
                ux += wx * du; uy += wy * du; vx += wx * dv; vy += wy * dv
                emin = min(emin, eta[j]); emax = max(emax, eta[j])
                hmin = min(hmin, h[j]); hmax = max(hmax, h[j])
                umin = min(umin, u[j]); umax = max(umax, u[j])
                vmin = min(vmin, v[j]); vmax = max(vmax, v[j])
            end
        end
        if ok
            ae = _bj(eta[i], emin, emax, ex, ey, i, a)
            ah = _bj(h[i], hmin, hmax, hx, hy, i, a)
            au = _bj(u[i], umin, umax, ux, uy, i, a)
            av = _bj(v[i], vmin, vmax, vx, vy, i, a)
            a.dex[i] = ae * ex; a.dey[i] = ae * ey
            a.dhx[i] = ah * hx; a.dhy[i] = ah * hy
            a.dux[i] = au * ux; a.duy[i] = au * uy
            a.dvx[i] = av * vx; a.dvy[i] = av * vy
        else
            a.dex[i] = 0.0; a.dey[i] = 0.0; a.dhx[i] = 0.0; a.dhy[i] = 0.0
            a.dux[i] = 0.0; a.duy[i] = 0.0; a.dvx[i] = 0.0; a.dvy[i] = 0.0
        end
    end
    return nothing
end

# interior faces: reconstructed states, hydrostatic reconstruction, HLLC,
# pressure corrections, optional viscous stress
@inline function _swe_face_k!(e, a)
    @inbounds begin
        h = a.h; u = a.u; v = a.v; z = a.z
        i = a.ci[e]; j = a.cj[e]
        if !(h[i] > SWE_HDRY || h[j] > SWE_HDRY) || !isfinite(a.eta[i]) || !isfinite(a.eta[j])
            a.F1[e] = 0.0; a.F2[e] = 0.0; a.F3[e] = 0.0; a.PL[e] = 0.0; a.PR[e] = 0.0
            a.T2[e] = 0.0; a.T3[e] = 0.0
        else
            ax = a.rxi[e]; ay = a.ryi[e]; bx = a.rxj[e]; by = a.ryj[e]
            # reconstructed depth and face bed: z_L = z_i + Δη − Δh (z_i at first order)
            dhL = a.dhx[i] * ax + a.dhy[i] * ay; dhR = a.dhx[j] * bx + a.dhy[j] * by
            dzL = a.dex[i] * ax + a.dey[i] * ay - dhL; dzR = a.dex[j] * bx + a.dey[j] * by - dhR
            uL = u[i] + a.dux[i] * ax + a.duy[i] * ay;   uR = u[j] + a.dux[j] * bx + a.duy[j] * by
            vL = v[i] + a.dvx[i] * ax + a.dvy[i] * ay;   vR = v[j] + a.dvx[j] * bx + a.dvy[j] * by
            hL = h[i] > SWE_HDRY ? max(h[i] + dhL, 0.0) : 0.0
            hR = h[j] > SWE_HDRY ? max(h[j] + dhR, 0.0) : 0.0
            zL = z[i] + dzL; zR = z[j] + dzR
            zf = max(zL, zR)
            a.sill[e] > max(z[i], z[j]) && (zf = max(zf, a.sill[e]))   # an obstacle on the face
            hLs = max(hL + zL - zf, 0.0); hRs = max(hR + zR - zf, 0.0)
            n1 = a.nx[e]; n2 = a.ny[e]
            unL = uL * n1 + vL * n2; utL = -uL * n2 + vL * n1
            unR = uR * n1 + vR * n2; utR = -uR * n2 + vR * n1
            f1, fn, ft = _swe_hllc(hLs, unL, utL, hRs, unR, utR)
            a.F1[e] = f1; a.F2[e] = fn * n1 - ft * n2; a.F3[e] = fn * n2 + ft * n1
            # relative to cell i's pressure, plus the centred bed source
            a.PL[e] = 0.5 * _G * (hL * hL - hLs * hLs - h[i] * h[i] + (hL + h[i]) * dzL)
            a.PR[e] = 0.5 * _G * (hR * hR - hRs * hRs - h[j] * h[j] + (hR + h[j]) * dzR)
            if a.ce > 0.0 && h[i] > SWE_HDRY && h[j] > SWE_HDRY
                hf = min(h[i], h[j])
                sp = 0.5 * (hypot(u[i], v[i]) + hypot(u[j], v[j]))
                ustar = sqrt(_G) * 0.5 * (a.nman[i] + a.nman[j]) * sp / hf^(1.0 / 6.0)
                k = a.ce * hf * ustar * hf / a.L[e]             # ν·h/L
                a.T2[e] = k * (u[j] - u[i]); a.T3[e] = k * (v[j] - v[i])
            else
                a.T2[e] = 0.0; a.T3[e] = 0.0
            end
        end
    end
    return nothing
end

# boundary faces: the ghost state from the BC code, fluxes out of the cell
@inline function _swe_bnd_k!(k, a)
    @inbounds begin
        c = a.gcell[k]; h = a.h; z = a.z
        if !(h[c] > SWE_HDRY) || !isfinite(a.eta[c])
            a.G1[k] = 0.0; a.G2[k] = 0.0; a.G3[k] = 0.0
        else
            ax = a.grx[k]; ay = a.gry[k]; n1 = a.gnx[k]; n2 = a.gny[k]
            dhL = a.dhx[c] * ax + a.dhy[c] * ay; dzL = a.dex[c] * ax + a.dey[c] * ay - dhL
            uL = a.u[c] + a.dux[c] * ax + a.duy[c] * ay; vL = a.v[c] + a.dvx[c] * ax + a.dvy[c] * ay
            hL = max(h[c] + dhL, 0.0); zL = z[c] + dzL
            zf = a.gz[k] > z[c] ? max(zL, a.gz[k]) : zL; hs = max(hL + zL - zf, 0.0)
            unL = uL * n1 + vL * n2; utL = -uL * n2 + vL * n1
            b = a.gbc[k]
            if b == SWE_BC_WALL || (b == SWE_BC_TRANSMISSIVE && unL <= 0.0 && !a.open_in)
                f1, fn, ft = _swe_hllc(hs, unL, utL, hs, -unL, utL); f1 = 0.0; ft = 0.0
            elseif b == SWE_BC_CRITICAL
                f1, fn, ft = _swe_hllc(hs, unL, utL, hs, max(unL, sqrt(_G * hs)), utL)
            else
                f1, fn, ft = _swe_hllc(hs, unL, utL, hs, unL, utL)
            end
            P = 0.5 * _G * (hL * hL - hs * hs - h[c] * h[c] + (hL + h[c]) * dzL)
            a.G1[k] = f1; a.G2[k] = (fn + P) * n1 - ft * n2; a.G3[k] = (fn + P) * n2 + ft * n1
        end
    end
    return nothing
end

# per-cell gather and update, semi-implicit Manning friction
@inline function _swe_update_k!(i, a, dt)
    @inbounds begin
        A = a.A; ci = a.ci; w = a.w
        (isfinite(a.z[i]) && A[i] >= 1.0) || return nothing
        dm = 0.0; dx = 0.0; dy = 0.0
        for k in a.cptr[i]:a.cptr[i + 1] - 1
            e = a.cedge[k]; we = w[e]
            if ci[e] == i                         # flux leaves i along n
                dm -= we * a.F1[e]
                dx -= we * (a.F2[e] + a.PL[e] * a.nx[e] - a.T2[e])
                dy -= we * (a.F3[e] + a.PL[e] * a.ny[e] - a.T3[e])
            else                                  # enters j
                dm += we * a.F1[e]
                dx += we * (a.F2[e] + a.PR[e] * a.nx[e] - a.T2[e])
                dy += we * (a.F3[e] + a.PR[e] * a.ny[e] - a.T3[e])
            end
        end
        for k in a.gptr[i]:a.gptr[i + 1] - 1
            f = a.gface[k]; wf = a.gw[f]
            dm -= wf * a.G1[f]; dx -= wf * a.G2[f]; dy -= wf * a.G3[f]
        end
        Vn = a.V[i] + dt * dm
        if Vn < 0.0
            a.clamp[i] -= Vn; Vn = 0.0
        end
        hn = Vn / A[i]
        qxn = a.qx[i] + dt * dx / A[i]; qyn = a.qy[i] + dt * dy / A[i]
        if hn > SWE_HDRY
            qm = hypot(a.qx[i], a.qy[i])          # |q| at the stage start (fixed point = Manning)
            fr = 1.0 + dt * _G * a.nman[i]^2 * qm / hn^(7.0 / 3.0)
            qxn /= fr; qyn /= fr
        else
            qxn = 0.0; qyn = 0.0
        end
        a.V[i] = Vn; a.qx[i] = qxn; a.qy[i] = qyn
    end
    return nothing
end

@inline function _swe_rk_avg_k!(i, a)
    @inbounds begin
        a.V[i] = 0.5 * (a.V0[i] + a.V[i]); a.qx[i] = 0.5 * (a.qx0[i] + a.qx[i])
        a.qy[i] = 0.5 * (a.qy0[i] + a.qy[i])
    end
    return nothing
end

@inline function _swe_depth_k!(i, a)
    @inbounds (a.A[i] >= 1.0) && (a.d[i] = a.V[i] / a.A[i])
    return nothing
end

@inline function _swe_vel_k!(i, a)
    @inbounds begin
        h = a.d[i]
        if h > SWE_HDRY
            a.vu[i] = _swe_vel(a.qx[i], h); a.vv[i] = _swe_vel(a.qy[i], h)
        else
            a.vu[i] = 0.0; a.vv[i] = 0.0
        end
        a.vel[i] = hypot(a.vu[i], a.vv[i])
    end
    return nothing
end

# CFL: per-cell signal speed into `u`, per-cell admissible step into `v`
# (both are stage scratch, recomputed from the state at the next stage)
@inline function _swe_speed_k!(i, a)
    @inbounds begin
        Ai = a.A[i]; h = Ai >= 1.0 && isfinite(a.z[i]) ? a.V[i] / Ai : 0.0
        a.u[i] = h > SWE_HDRY ? hypot(_swe_vel(a.qx[i], h), _swe_vel(a.qy[i], h)) + sqrt(_G * h) : 0.0
    end
    return nothing
end

@inline function _swe_dtcell_k!(i, a)
    @inbounds begin
        Ai = a.A[i]
        if !(Ai >= 1.0 && isfinite(a.z[i]))
            a.v[i] = Inf
            return nothing
        end
        sp = a.u; rate = 0.0
        for k in a.cptr[i]:a.cptr[i + 1] - 1
            e = a.cedge[k]; rate += a.w[e] * max(sp[a.ci[e]], sp[a.cj[e]])
        end
        for k in a.gptr[i]:a.gptr[i + 1] - 1
            rate += a.gw[a.gface[k]] * sp[i]
        end
        if a.ce > 0.0
            h = a.V[i] / Ai
            if h > SWE_HDRY
                u = hypot(_swe_vel(a.qx[i], h), _swe_vel(a.qy[i], h))
                ν = a.ce * h * sqrt(_G) * a.nman[i] * u / h^(1.0 / 6.0)
                rate += 4.0 * ν * a.perim[i] / (2Ai / a.perim[i])
            end
        end
        a.v[i] = rate > 0.0 ? Ai / rate : Inf
    end
    return nothing
end

# ── Drivers (any backend; `a` from _swe_args, or its device copy) ──────────

"Δt of `_cfl_dt(state, SWEFlow())` on the arrays `a`."
function _swe_cfl(be::ComputeBackend, a, courant::Float64)
    foreach_index!(be, a.n, _swe_speed_k!, a)
    foreach_index!(be, a.n, _swe_dtcell_k!, a)
    best = _minimum(a.v)
    return isfinite(best) ? courant * best : 60.0
end

# one forward-Euler stage U ← U + Δt·L(U); returns the boundary outflow volume
function _swe_stage!(be::ComputeBackend, a, dt::Float64)
    foreach_index!(be, a.n, _swe_prim_k!, a)
    a.order == 2 && foreach_index!(be, a.n, _swe_grad_k!, a)
    foreach_index!(be, a.ne, _swe_face_k!, a)
    foreach_index!(be, a.ng, _swe_bnd_k!, a)
    out = _dot_sum(dt, a.gw, a.G1)
    foreach_index!(be, a.n, _swe_update_k!, a, dt)
    return out
end

# one sub-step; returns the volume that left through the boundary
function _swe_substep!(be::ComputeBackend, a, dt::Float64)
    a.order == 1 && return _swe_stage!(be, a, dt)
    copyto!(a.V0, a.V); copyto!(a.qx0, a.qx); copyto!(a.qy0, a.qy)
    o1 = _swe_stage!(be, a, dt)
    o2 = _swe_stage!(be, a, dt)
    foreach_index!(be, a.n, _swe_rk_avg_k!, a)
    return 0.5 * (o1 + o2)
end

"""
    _swe_advance!(be, a, dt, courant, removed) → number of sub-steps

`step_swe!` on the arrays `a` (host or device): sub-steps to `dt`, adds each
sub-step's boundary outflow to `removed[]`, and syncs the depth.
"""
function _swe_advance!(be::ComputeBackend, a, dt::Float64, courant::Float64,
                       removed::Base.RefValue{Float64})
    left = dt; nsub = 0
    while left > 0.0
        h = min(left, _swe_cfl(be, a, courant))
        h >= left * (1 - 1e-9) && (h = left)                # no sliver sub-step at the end
        removed[] += _swe_substep!(be, a, h)
        left -= h; nsub += 1
    end
    foreach_index!(be, a.n, _swe_depth_k!, a)
    return nsub
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
    @inbounds for k in eachindex(sw.gbc); sw.gbc[k] = _swe_bc_code(state.ghost_cell_bc[k]); end
    removed = Ref(state.vol_removed)
    nsub = _swe_advance!(CPUThreads(), _swe_args(state), dt, courant, removed)
    state.vol_removed = removed[]
    return nsub
end

"""
    _swe_velocity!(state)

Cell velocity for output: u = q/h (desingularised), its magnitude.
"""
function _swe_velocity!(state::FlowState)
    a = _swe_args(state)
    foreach_index!(CPUThreads(), a.n, _swe_vel_k!, a)
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
_cfl_dt(state::FlowState, ::SWEFlow; courant::Float64 = 0.7)::Float64 =
    _swe_cfl(CPUThreads(), _swe_args(state), courant)
