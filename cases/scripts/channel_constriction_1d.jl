# cases/scripts/channel_constriction_1d.jl — steady 1-D level profiles through the
# opening of `channel_constriction`, for each set of equations (included by
# cases/registry/channel_constriction.jl).
#
# Friction as the 2-D models have it (on the bed only, no walls): with
# q = Q / W(s) the discharge per metre, S_f = n² q² / h^(10/3). Marching
# upstream from a level η_d at s_d (downstream of the opening):
#
#   :li      local-inertial: no convective acceleration, so the steady
#            momentum balance is g h ∂η/∂s = −g n² q²/h^(7/3): dη/ds = −S_f,
#            and η is continuous through the transitions.
#   :energy  the full equations, no losses at the transitions: the energy head
#            H = η + q²/(2g h²) falls by S_f along the channel and is
#            continuous through them (the water surface drops into the opening
#            and recovers after it).
#   :borda   :energy, with the abrupt expansion's loss: the shallow-water
#            equations conserve momentum, so across a sudden widening the head
#            falls by the Borda–Carnot loss (v_o − v_w)²/2g (v_o, v_w the
#            velocities in the opening and in the wide channel just below it).
#            The abrupt contraction is taken as lossless.
#   :abrupt  :energy, with the step-backwater coefficients for abrupt
#            transitions (HEC-RAS: 0.6 of the change in velocity head at a
#            contraction, 0.8 at an expansion) — the losses of a real
#            sharp-cornered opening, flow separation and vena contracta
#            included; an upper bound for the shallow-water equations here.
#
# Subcritical roots throughout (an error if there is none).

const _CC_G = 9.81

_cc_sf(q, h, n) = n^2 * q^2 / h^(10 / 3)

"Bisection for the root of an increasing `f` on [lo, hi]."
function _cc_root(f, lo, hi)
    for _ in 1:200
        m = 0.5 * (lo + hi); f(m) > 0 ? (hi = m) : (lo = m)
    end
    return 0.5 * (lo + hi)
end

"Subcritical depth with specific energy `E` for discharge per metre `q`."
function _cc_depth(E, q)
    hc = (q^2 / _CC_G)^(1 / 3)
    Es(h) = h + q^2 / (2 * _CC_G * h^2)
    E < Es(hc) - 1e-9 && error("channel_constriction_1d: no subcritical depth (energy $E below critical $(Es(hc)))")
    return _cc_root(h -> Es(h) - E, hc, max(E, hc) + 1.0)
end

"""
    cc_profile(model, η_d, s_d, s_up; Q, n, W, zb, edges, ds = 0.5) -> (s, η)

The level from s_d (level η_d) upstream to s_up; `model` ∈ (:li, :energy, :abrupt,
:borda), `W(s)` the width, `zb(s)` the bed, `edges` the positions where the
width changes.
"""
function cc_profile(model::Symbol, η_d, s_d, s_up; Q, n, W, zb, edges, ds = 0.5)
    model in (:li, :energy, :borda, :abrupt) || throw(ArgumentError("model must be :li, :energy, :borda or :abrupt"))
    ss = [Float64(s_d)]; ηs = [Float64(η_d)]
    s = Float64(s_d); η = Float64(η_d)
    while s > s_up + 1e-9
        s1 = max(s - ds, s_up)
        for e in edges
            s1 < e < s && (s1 = e)                     # stop at a width change
        end
        sm = 0.5 * (s + s1); q = Q / W(sm); L = s - s1
        if model === :li
            ηm = η + _cc_sf(q, η - zb(s), n) * L / 2
            η += _cc_sf(q, ηm - zb(sm), n) * L
        else
            H = η + q^2 / (2 * _CC_G * (η - zb(s))^2)
            Hm = H + _cc_sf(q, η - zb(s), n) * L / 2
            H += _cc_sf(q, _cc_depth(Hm - zb(sm), q), n) * L
            η = zb(s1) + _cc_depth(H - zb(s1), q)
        end
        s = s1
        if model !== :li && s > s_up && any(e -> abs(s - e) < 1e-9, edges)
            Wd, Wu = W(s + 1e-6), W(s - 1e-6)           # widths below and above the change
            vd = Q / (Wd * (η - zb(s))); Ed = η - zb(s) + vd^2 / (2 * _CC_G)
            qu = Q / Wu; hc = (qu^2 / _CC_G)^(1 / 3)
            if model === :borda && Wu < Wd
                # the opening above a sudden expansion: E_u − (v_u − v_d)²/2g = E_d
                f(h) = h + (qu / h)^2 / (2 * _CC_G) - (qu / h - vd)^2 / (2 * _CC_G) - Ed
                η = zb(s) + _cc_root(f, hc, η - zb(s) + 5.0)
            elseif model === :abrupt
                # expansion (opening above): E_u = E_d + 0.8 (v_u² − v_d²)/2g;
                # contraction (wide above): E_u = E_d + 0.6 (v_d² − v_u²)/2g:
                # h + k v_u²/2g = rhs
                k, rhs = Wu < Wd ? (1 - 0.8, Ed - 0.8 * vd^2 / (2 * _CC_G)) : (1 + 0.6, Ed + 0.6 * vd^2 / (2 * _CC_G))
                lo = max(hc, (k * qu^2 / _CC_G)^(1 / 3))
                g(h) = h + k * (qu / h)^2 / (2 * _CC_G) - rhs
                g(lo) > 0 && error("channel_constriction_1d: no subcritical depth at the transition (:abrupt)")
                η = zb(s) + _cc_root(g, lo, η - zb(s) + 5.0)
            else
                η = zb(s) + _cc_depth(Ed, Q / Wu)
            end
        end
        push!(ss, s); push!(ηs, η)
    end
    return ss, ηs
end
