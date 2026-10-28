# surfacewater/momentum_damping.jl
# --------------------------------
# D6′ (deep-channel damping, 2026-10): divergence-free damping of the
# circulating face discharges, `--momentum-damping C_m` (standard flow;
# 0 = off, the default). Companion of the η filter (eta_filter.jl): the
# filter damps the surface modes, this damps the modes the surface cannot see.
#
# Discharge patterns that circulate around a mesh vertex move no water, so
# no surface filter reaches them; Q-centring damped them, at a rate
# (1 − θ)/Δt and at a cost in conveyance (HYDRAULICS.md §15.8). Here, per
# interior mesh vertex v (cells all round it):
#     Γ_v = Σ_{faces e at v} s_{v,e}·u_e·L_e          circulation of the face velocity
#     ζ_v = Γ_v / A_v,  A_v = Σ_e w_e·L_e / 4          vorticity, dual area
#     ψ_v = −ν_v·Δt·h_v·ζ_v·ω_v                        streamfunction (m³/s)
# and every face's discharge changes by ΔQ_e = ψ_head − ψ_tail (its two end
# vertices, ordered so that the tangent head − tail is the normal turned
# +90°). Around any cell the ΔQ telescope to zero: the change is exactly
# divergence-free, so it moves no water and leaves the gravity-wave
# dynamics (and their stability limit) alone. It vanishes for a uniform
# velocity field (Γ = 0) and acts like an eddy viscosity ν on the rotational
# part of the flow, ∂u/∂t = −ν·curl ζ.
#
#   ν_v = C_m·V·ℓ_v,  ℓ_v = √A_v,   V = √(g h_v)  (:wave, a numerical floor)
#                                   V = |u|_v     (:flow, eddy-viscosity-like)
#
# The rate is in physical units: the damping per unit time does not depend
# on Δt. Explicit stability of the vortex diffusion needs 8·C_m·C < 2 on a
# square grid for :wave at Courant number C (C_m < 0.36 at C = 0.7); :flow is
# far from that limit.
#
# u_e is the velocity along the line between the two cell centres (the dual
# edge), u_e = (q_n·(n̂·d̂) + q_t·(t̂·d̂)) / h_flow: the face discharge q_n plus
# the face-tangential part q_t of the mean cell discharge vector (Phase F),
# with h_flow = max(η_i, η_j) − z_sill as in the Bates kernel. On a skewed
# mesh (A5) n̂ ∦ d̂, and the normal part alone would give a uniform flow a
# circulation; with both parts Γ = U·Σ d_e = 0 exactly for any uniform U.
# h_v is the mean of h_flow over the vertex's faces. Sign convention: the
# stored edges.flux is positive from cell_j to cell_i (the Bates update is
# q − g·h·Δt·(η_i − η_j)/L), so the discharge along n̂ (cell_i → cell_j) is
# q_n = −flux, while the cell vectors are physical (qvec·n̂ = −flux); the
# streamfunction change ΔQ along n̂ is stored as flux −= ΔQ / w. ω_v ramps from
# 0 to 1 as the shallowest face at v goes from 0 to MOM_DAMP_HW of depth
# (min(η_i, η_j) − z_sill), so a vertex next to a dry or barely wet face is
# left alone, continuously. Vertices on the domain edge keep ψ = 0, which
# keeps the change divergence-free for boundary cells too. The change is
# applied to the stored face discharges at the end of the step (after the
# Phase F reconstruction), so it acts through the next step's momentum (like
# Q-centring's smoothing did).

const MOM_DAMP_HW = 0.05       # m: face depth over which a vertex's weight ramps to 1

# Vertex topology from the diamond table's per-face vertex keys: oriented
# tail/head vertex per face, and a CSR list of each interior vertex's faces
# with the circulation sign (+1 if the vertex is the face's head).
function _mom_damp_ws!(state::FlowState, ws::StepWorkspace)
    E = state.edges; ne = E.n_edges; n = length(state.cell_ids)
    length(ws.m_tail) == ne && return ws
    T = state.diamond_table
    T === nothing && throw(ArgumentError("--momentum-damping needs the diamond table " *
                                         "(--gradient-correction on --face-flux-method diamond|diamond-wetdry)"))
    keyidx = Dict{Tuple{Float64,Float64},Int}()
    vid(k) = get!(keyidx, k, length(keyidx) + 1)
    tail = zeros(Int, ne); head = zeros(Int, ne)
    R = A5Grid._EARTH_R
    @inbounds for e in 1:ne
        ka = T.va_key[e]; kb = T.vb_key[e]
        (ka == (0.0, 0.0) || kb == (0.0, 0.0) || ka == kb) && continue
        nx = E.nf_x[e]; ny = E.nf_y[e]
        (nx == 0.0 && ny == 0.0) && continue
        c = cosd(0.5 * (ka[2] + kb[2]))
        tx = deg2rad(kb[1] - ka[1]) * R * c; ty = deg2rad(kb[2] - ka[2]) * R
        a = vid(ka); b = vid(kb)
        if nx * ty - ny * tx >= 0.0                # t = (b − a) is n turned +90°
            tail[e] = a; head[e] = b
        else
            tail[e] = b; head[e] = a
        end
    end
    nv = length(keyidx)
    cnt = zeros(Int, nv)
    for e in 1:ne
        tail[e] == 0 && continue
        cnt[tail[e]] += 1; cnt[head[e]] += 1
    end
    vptr = Vector{Int}(undef, nv + 1); vptr[1] = 1
    for v in 1:nv; vptr[v + 1] = vptr[v] + cnt[v]; end
    vedge = zeros(Int, vptr[end] - 1); vsign = zeros(Float64, vptr[end] - 1); fill!(cnt, 0)
    for e in 1:ne                                  # faces of each vertex in increasing edge order
        tail[e] == 0 && continue
        for (v, s) in ((tail[e], -1.0), (head[e], 1.0))
            k = vptr[v] + cnt[v]; vedge[k] = e; vsign[k] = s; cnt[v] += 1
        end
    end
    # interior: ≥ 3 faces and every cell round the vertex on exactly two of them
    interior = falses(nv)
    for v in 1:nv
        m = vptr[v + 1] - vptr[v]
        m >= 3 || continue
        seen = Dict{Int,Int}()
        for k in vptr[v]:vptr[v + 1] - 1
            e = vedge[k]
            seen[E.cell_i[e]] = get(seen, E.cell_i[e], 0) + 1
            seen[E.cell_j[e]] = get(seen, E.cell_j[e], 0) + 1
        end
        interior[v] = length(seen) == m && all(==(2), values(seen))
    end
    ws.m_tail = tail; ws.m_head = head; ws.m_vptr = vptr; ws.m_vedge = vedge; ws.m_vsign = vsign
    ws.m_interior = interior; ws.m_psi = zeros(nv); ws.m_u = zeros(ne); ws.m_w = zeros(ne)
    ws.m_h = zeros(ne)
    return ws
end

"""
    _momentum_damping!(state, ws, dt)

Phase D3 of `step_standard!` when `state.mom_damping_cs > 0` (see the file
header). Threaded per edge / per vertex, fixed-order sums: the result does
not depend on the thread count.
"""
function _momentum_damping!(state::FlowState, ws::StepWorkspace, dt::Float64)
    _mom_damp_ws!(state, ws)
    Cm = state.mom_damping_cs; wave = state.mom_damping_scale === :wave
    E = state.edges; ne = E.n_edges
    z = state.elevation; V = state.volume; A = state.cell_area
    ci = E.cell_i; cj = E.cell_j; sill = E.sill; L = E.L; wid = E.width; q = E.flux
    u = ws.m_u; wf = ws.m_w; hfl = ws.m_h; psi = ws.m_psi
    nfx = E.nf_x; nfy = E.nf_y; dxm = E.dx_m; dym = E.dy_m
    qu = state.qvec_u; qv = state.qvec_v
    tail = ws.m_tail; head = ws.m_head; vptr = ws.m_vptr; vedge = ws.m_vedge; vsign = ws.m_vsign
    interior = ws.m_interior; nv = length(psi)

    @batch per=core for e in 1:ne
        @inbounds begin
            i = ci[e]; j = cj[e]
            ηi = z[i] + V[i] / max(A[i], 1.0); ηj = z[j] + V[j] / max(A[j], 1.0)
            hmin = min(ηi, ηj) - sill[e]
            if hmin > 0.0 && tail[e] != 0          # also false for NaN
                hf = max(ηi, ηj) - sill[e]
                nx = nfx[e]; ny = nfy[e]; dx = dxm[e] / L[e]; dy = dym[e] / L[e]
                qt = 0.5 * ((qu[i] + qu[j]) * -ny + (qv[i] + qv[j]) * nx)   # t̂ = n̂ turned +90°
                hfl[e] = hf
                u[e] = (-q[e] * (nx * dx + ny * dy) + qt * (-ny * dx + nx * dy)) / hf
                wf[e] = min(hmin / MOM_DAMP_HW, 1.0)
            else
                hfl[e] = 0.0; u[e] = 0.0; wf[e] = 0.0
            end
        end
    end
    @batch per=core for v in 1:nv
        @inbounds begin
            ps = 0.0
            if interior[v]
                Γ = 0.0; Av = 0.0; ω = 1.0; hs = 0.0; us = 0.0; m = 0
                for k in vptr[v]:vptr[v + 1] - 1
                    e = vedge[k]
                    ω = min(ω, wf[e])
                    Γ += vsign[k] * u[e] * L[e]
                    Av += wid[e] * L[e] / 4
                    hs += hfl[e]; us += abs(u[e]); m += 1
                end
                if ω > 0.0 && Av > 0.0
                    hv = hs / m
                    Vref = wave ? sqrt(_G * hv) : us / m
                    ν = Cm * Vref * sqrt(Av)
                    ps = -ν * dt * hv * (Γ / Av) * ω
                end
            end
            psi[v] = ps
        end
    end
    @batch per=core for e in 1:ne
        @inbounds begin
            a = tail[e]
            if a != 0
                dQ = psi[head[e]] - psi[a]
                iszero(dQ) || (q[e] -= dQ / wid[e])
            end
        end
    end
    return nothing
end
