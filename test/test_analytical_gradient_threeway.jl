#!/usr/bin/env julia
# test/test_analytical_gradient_threeway.jl
#
# ChatGPT review, Priority/Test 1-3: the sharpest available diagnostic for
# separating gradient-RECONSTRUCTION error (diamond vs legacy) from
# everything downstream of it (Bates kernel, Froude/volume limiters).
#
# A manufactured LINEAR WSE field is assigned directly to cell values --
# no DEM, no real simulation, no wet/dry state -- so the analytically
# exact face-normal derivative is known exactly:
#     g_n_exact = Bx*n_x + By*n_y
# (n_x, n_y = the same oriented face normal already stored as
# edges.nf_x/nf_y -- not recomputed independently).
#
# Three stages, matching the review's Test 1/2/3:
#   Test 1: g_n itself (pure geometry/reconstruction accuracy), broken
#           down by non-orthogonality angle, bearing, and vertex-fallback
#           class (review Priority 7)
#   Test 2: raw (unlimited) Bates Q from each of the three g_n sources
#   Test 3: Froude+volume-LIMITED Q from each source, plus a proxy check
#           of how often the Froude cap actually engages for each
#
# Uses a flat bed (elevation=0 everywhere) and a "reservoir" WSE field
# (large constant offset A plus a small linear tilt B) so h_flow is large
# and roughly uniform, isolating the slope term's effect on the flux
# kernel without wet/dry-front complications. Full time-evolution effects
# (q_prev feedback compounding over many steps) are NOT tested here --
# this is a single-step (q_prev=0) kernel test. See the note at the end
# for how to extend this to a multi-step version if this doesn't resolve
# things on its own.
#
# Usage:
#   julia --project=. test/test_analytical_gradient_threeway.jl \
#       <mesh.parquet> [A_offset] [Bx] [By]
#
# A_offset: reservoir depth offset (m), default 5.0
# Bx, By:   WSE gradient components (m/m), default 0.001, -0.0015
#           (comparable magnitude to the 0.1% planar-slope tests already
#           used elsewhere in this project)

using Statistics, LinearAlgebra, Printf

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

mesh_path = length(ARGS) >= 1 ? ARGS[1] :
            joinpath("test", "square", "square_mesh_res14.parquet")
A_off = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 5.0
Bx    = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.001
By    = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : -0.0015

mesh = load_mesh_geoparquet(mesh_path)
println("Loaded $(length(mesh.cells)) cells from $mesh_path")
println("Manufactured field: WSE = $A_off + $Bx*x + $By*y  (local metres, flat bed z=0)")

state = initialise_flow_model(mesh, StandardFlow();
                               gradient_correction       = true,
                               gradient_correction_alpha = 1.0,
                               face_flux_method           = :diamond)

n     = length(state.cell_ids)
edges = state.edges
ne    = edges.n_edges

# ── manufactured field: flat bed, synthetic WSE, no real simulation ────
state.elevation .= 0.0
lon0 = mean(state.cell_lons); lat0 = mean(state.cell_lats)
cos_lat0 = cosd(lat0); R = A5Grid._EARTH_R
to_xy(lon, lat) = (deg2rad(lon - lon0) * R * cos_lat0, deg2rad(lat - lat0) * R)

wse_all = Vector{Float64}(undef, n)
for i in 1:n
    x, y = to_xy(state.cell_lons[i], state.cell_lats[i])
    wse_all[i] = A_off + Bx * x + By * y
end

_compute_wse_gradients!(state, wse_all)
alpha = state.gradient_correction_alpha

# ── Test 1: g_n comparison ──────────────────────────────────────────
g_exact     = fill(NaN, ne)
g_diamond   = fill(NaN, ne)
g_legacy    = fill(NaN, ne)
theta_deg   = fill(NaN, ne)
bearing_deg = fill(NaN, ne)
va_fb = fill(:na, ne)
vb_fb = fill(:na, ne)
is_diamond_valid = falses(ne)

for e in 1:ne
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]

    nx, ny = edges.nf_x[e], edges.nf_y[e]
    g_exact[e] = Bx * nx + By * ny

    gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
    gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
    vhat_dot = gx_f * edges.skew_x[e] + gy_f * edges.skew_y[e]
    dwse_legacy = edges.cos_theta[e] * (wse_ci - wse_cj) - alpha * edges.L[e] * vhat_dot
    g_legacy[e] = -dwse_legacy / max(edges.L[e], 1.0)

    dn_dia = diamond_dWSE_n(state.diamond_table, e, wse_ci, wse_cj, wse_all, edges.L[e])
    if !isnan(dn_dia)
        g_diamond[e] = -dn_dia / max(edges.L[e], 1.0)
        is_diamond_valid[e] = true
        va_fb[e] = state.diamond_table.vertex_table[state.diamond_table.va_key[e]].fallback
        vb_fb[e] = state.diamond_table.vertex_table[state.diamond_table.vb_key[e]].fallback
    end

    theta_deg[e] = acosd(clamp(edges.cos_theta[e], -1.0, 1.0))
    lat_mid = 0.5 * (state.cell_lats[ci] + state.cell_lats[cj])
    dlat = state.cell_lats[cj] - state.cell_lats[ci]
    dlon = (state.cell_lons[cj] - state.cell_lons[ci]) * cosd(lat_mid)
    bearing_deg[e] = mod(atand(dlon, dlat), 360.0)
end

valid_idx = findall(is_diamond_valid)
println("
$(length(valid_idx))/$ne edges have a valid diamond record.")

err_dia = g_diamond[valid_idx] .- g_exact[valid_idx]
err_leg = g_legacy[valid_idx]  .- g_exact[valid_idx]

println("
=== Test 1: face-normal derivative g_n, vs analytical exact ===")
@printf("  diamond error: mean=%.3e  std=%.3e  max|.|=%.3e\n",
        mean(err_dia), std(err_dia), maximum(abs.(err_dia)))
@printf("  legacy  error: mean=%.3e  std=%.3e  max|.|=%.3e\n",
        mean(err_leg), std(err_leg), maximum(abs.(err_leg)))
if maximum(abs.(err_leg)) > 1e-9
    @printf("  ratio max|diamond err| / max|legacy err| = %.6f\n",
            maximum(abs.(err_dia)) / maximum(abs.(err_leg)))
end
println("  If diamond error is ~machine precision and legacy error is not,")
println("  Layer 3 (gradient reconstruction) is confirmed solved on THIS")
println("  real mesh, and the residual bias lives downstream (Bates/limiter).")
println("  If diamond error is NOT ~machine precision, geometry itself --")
println("  this specific mesh's vertex reconstruction -- still needs work.")

println("
--- Diamond error vs non-orthogonality angle theta ---")
th = theta_deg[valid_idx]
n_bins = 5
lo_t, hi_t = extrema(th)
edges_t = range(lo_t, hi_t, length = n_bins + 1)
for b in 1:n_bins
    tlo, thi = edges_t[b], edges_t[b+1]
    idx = findall(i -> th[i] >= tlo && th[i] < (b == n_bins ? thi + 1e-9 : thi), 1:length(th))
    isempty(idx) && continue
    @printf("  theta in [%.1f, %.1f) deg  n=%5d  mean|err|=%.3e  max|err|=%.3e\n",
            tlo, thi, length(idx), mean(abs.(err_dia[idx])), maximum(abs.(err_dia[idx])))
end

bea = bearing_deg[valid_idx]
cos2b = cosd.(2 .* bea)
println("
--- Diamond error bearing correlation ---")
@printf("  corr(err_diamond, cos(2*bearing)) = %.4f\n", cor(err_dia, cos2b))

println("
--- Diamond error by vertex-fallback class (review Priority 7) ---")
va_v = va_fb[valid_idx]; vb_v = vb_fb[valid_idx]
both_none     = (va_v .== :none) .& (vb_v .== :none)
one_fallback  = (va_v .!= :none) .⊻ (vb_v .!= :none)
both_fallback = (va_v .!= :none) .& (vb_v .!= :none)
for (lbl, mask) in (("both vertices k>=3 well-conditioned", both_none),
                     ("one vertex fallback",                 one_fallback),
                     ("both vertices fallback",               both_fallback))
    nb = count(mask)
    nb == 0 && continue
    @printf("  %-38s n=%6d (%.1f%%)  mean|err|=%.3e  max|err|=%.3e\n",
            lbl, nb, 100 * nb / length(valid_idx), mean(abs.(err_dia[mask])), maximum(abs.(err_dia[mask])))
end
println("  This directly tests the review's point: 'linear exact' should be")
println("  qualified to 'linear exact for edges whose vertex reconstructions")
println("  are k>=3 and well-conditioned'. If error is concentrated in the")
println("  fallback rows above, that's the qualification made concrete.")

# ── Test 2/3: through the actual Bates kernel ────────────────────────
println("
=== Test 2/3: raw and limited Bates flux, analytic vs diamond vs legacy ===")
dt = 2.0
n_mann = 0.03

nv = length(valid_idx)
Q_raw_exact = fill(NaN, nv); Q_raw_dia = fill(NaN, nv); Q_raw_leg = fill(NaN, nv)
Q_lim_exact = fill(NaN, nv); Q_lim_dia = fill(NaN, nv); Q_lim_leg = fill(NaN, nv)
froude_bound_exact = falses(nv); froude_bound_dia = falses(nv); froude_bound_leg = falses(nv)

for (k, e) in enumerate(valid_idx)
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]
    # Flat-bed design (state.elevation was zeroed above): z_sill = 0.0
    # directly, NOT edges.sill[e] -- that field was precomputed from the
    # REAL DEM inside _build_edge_list before this script ever touches
    # state.elevation, so it still holds real terrain values (including
    # NaN wherever an edge touches one of the mesh's NaN-elevation cells).
    # Using it here silently broke the flat-bed assumption and, via NaN
    # propagation through mean()/std()/maximum(), poisoned every Test 2/3
    # aggregate even though only a few edges were directly affected.
    h_flow = max(wse_ci, wse_cj) - 0.0
    depth_donor = max(wse_ci, wse_cj)   # flat bed: WSE == depth exactly

    dWSE_exact = -g_exact[e]   * edges.L[e]
    dWSE_dia   = -g_diamond[e] * edges.L[e]
    dWSE_leg   = -g_legacy[e]  * edges.L[e]

    Q_raw_exact[k] = _bates_flux_corrected(0.0, h_flow, dWSE_exact, edges.width[e], edges.L[e], n_mann, dt)
    Q_raw_dia[k]   = _bates_flux_corrected(0.0, h_flow, dWSE_dia,   edges.width[e], edges.L[e], n_mann, dt)
    Q_raw_leg[k]   = _bates_flux_corrected(0.0, h_flow, dWSE_leg,   edges.width[e], edges.L[e], n_mann, dt)

    Ql_e, _ = _bates_flux_limited_corrected(0.0, h_flow, dWSE_exact, edges.width[e], edges.L[e], n_mann, dt, depth_donor)
    Ql_d, _ = _bates_flux_limited_corrected(0.0, h_flow, dWSE_dia,   edges.width[e], edges.L[e], n_mann, dt, depth_donor)
    Ql_l, _ = _bates_flux_limited_corrected(0.0, h_flow, dWSE_leg,   edges.width[e], edges.L[e], n_mann, dt, depth_donor)
    Q_lim_exact[k] = Ql_e; Q_lim_dia[k] = Ql_d; Q_lim_leg[k] = Ql_l

    q_max = h_flow > 0 ? h_flow * sqrt(_G * h_flow) * FROUDE_LIMIT : 0.0
    froude_bound_exact[k] = abs(Q_raw_exact[k] / max(edges.width[e], 1e-6)) > q_max
    froude_bound_dia[k]   = abs(Q_raw_dia[k]   / max(edges.width[e], 1e-6)) > q_max
    froude_bound_leg[k]   = abs(Q_raw_leg[k]   / max(edges.width[e], 1e-6)) > q_max
end

err_raw_dia = Q_raw_dia .- Q_raw_exact
err_raw_leg = Q_raw_leg .- Q_raw_exact
err_lim_dia = Q_lim_dia .- Q_lim_exact
err_lim_leg = Q_lim_leg .- Q_lim_exact

@printf("
  RAW (unlimited) Q:      diamond err mean=%.3e std=%.3e max=%.3e\n",
        mean(err_raw_dia), std(err_raw_dia), maximum(abs.(err_raw_dia)))
@printf("                          legacy  err mean=%.3e std=%.3e max=%.3e\n",
        mean(err_raw_leg), std(err_raw_leg), maximum(abs.(err_raw_leg)))

@printf("
  LIMITED (Froude+vol) Q: diamond err mean=%.3e std=%.3e max=%.3e\n",
        mean(err_lim_dia), std(err_lim_dia), maximum(abs.(err_lim_dia)))
@printf("                          legacy  err mean=%.3e std=%.3e max=%.3e\n",
        mean(err_lim_leg), std(err_lim_leg), maximum(abs.(err_lim_leg)))

d_exact_clip = Q_lim_exact .- Q_raw_exact
@printf("
  (reference) limiter's own effect on the EXACT source: mean|clip|=%.3e  max|clip|=%.3e\n",
        mean(abs.(d_exact_clip)), maximum(abs.(d_exact_clip)))
println("  If LIMITED errors are much larger, relative to RAW errors, than")
println("  the limiter's own reference clipping effect above, that's direct")
println("  evidence the limiter is AMPLIFYING the diamond/legacy discrepancy")
println("  rather than just passing it through unchanged.")

println("
--- Froude cap engagement rate (proxy: |q_raw| > q_max) ---")
@printf("  exact:   %d / %d edges (%.1f%%)\n", count(froude_bound_exact), nv, 100 * count(froude_bound_exact) / nv)
@printf("  diamond: %d / %d edges (%.1f%%)\n", count(froude_bound_dia),   nv, 100 * count(froude_bound_dia)   / nv)
@printf("  legacy:  %d / %d edges (%.1f%%)\n", count(froude_bound_leg),   nv, 100 * count(froude_bound_leg)   / nv)
println("  If diamond binds the Froude cap on a meaningfully different rate")
println("  or set of edges than legacy/exact, the limiter is a live")
println("  contributor to the residual bias, not a passive pass-through.")

println("
NOTE: this is a single-step (q_prev=0) kernel test -- it does not capture")
println("compounding effects from q_prev feedback over many timesteps. Test 4")
println("below extends this into exactly that multi-step version.")

# ── Test 4: momentum representation (:edge vs :cell) under repeated forcing ─
#
# Replicates the ACTUAL production formulas verbatim (read directly from
# FloodModel.jl, not reconstructed from memory):
#   :edge -- _q_centred(flux, e, collinear_i[e], collinear_j[e], theta),
#            flow2d.jl lines ~80-91
#   :cell -- Phase F reconstruction, step_standard! lines ~2395-2409:
#              qvec_i = Σ_s mom_weights[s,i]*(-flux[e_s]),  (u-row: 1:5)
#                       Σ_s mom_weights[5+s,i]*(-flux[e_s]) (v-row: 6:10)
#            then Phase A's q_prev_eff for :cell (line ~2242-2244):
#              q_prev_eff = -0.5*((qvec_u_ci+qvec_u_cj)*nf_x +
#                                  (qvec_v_ci+qvec_v_cj)*nf_y)
#
# FROZEN-COEFFICIENT design: the same diamond dWSE_n, h_flow, and donor
# depth from Test 2/3 are reused every virtual step -- volume/depth are
# NOT evolved. This deliberately isolates the momentum-storage + limiter
# ITERATION's own dynamics from real mass-conservation dynamics: given a
# FIXED driving force, does repeatedly applying Bates+limiter+momentum-
# storage converge to a stable flux field, and does :edge vs :cell storage
# reach a different fixed point, or show different oscillation behaviour
# getting there? That's a well-posed, standard "frozen-coefficient"
# stability question, and it's exactly what distinguishes the two momentum
# models -- a single q_prev=0 evaluation (Test 2/3 above) cannot, since
# both models start from identical, zero prior momentum.
#
# Uses full ne-sized arrays throughout (not just the valid-diamond subset)
# because _q_centred looks up COLLINEAR edges by global EdgeList index --
# those may not themselves be valid-diamond edges. On this mesh diamond
# coverage is 100% so this distinction doesn't change any number below,
# but indexing by the compact valid-only subset would silently break on
# any mesh with partial diamond coverage. Edges without a valid diamond
# record are simply never updated (flux stays 0, excluded from stats).

println("
=== Test 4: momentum representation (:edge vs :cell) under repeated forcing ===")

n_steps = 300
theta   = 0.9   # matches state.q_centre_theta's production default

h_flow_full       = fill(NaN, ne)
dWSE_full         = fill(NaN, ne)
depth_donor_full  = fill(NaN, ne)
for e in 1:ne
    is_diamond_valid[e] || continue
    ci, cj = edges.cell_i[e], edges.cell_j[e]
    wse_ci, wse_cj = wse_all[ci], wse_all[cj]
    h_flow_full[e]      = max(wse_ci, wse_cj) - 0.0   # flat bed, matches Test 2/3
    dWSE_full[e]        = -g_diamond[e] * edges.L[e]
    depth_donor_full[e] = max(wse_ci, wse_cj)
end

flux_edge = zeros(ne)
flux_cell = zeros(ne)
qvec_u    = zeros(n)
qvec_v    = zeros(n)

hist_edge = zeros(n_steps, ne)   # full history retained for convergence check
hist_cell = zeros(n_steps, ne)

n_mann4 = 0.03
dt4 = 2.0

for step in 1:n_steps
    global flux_edge, flux_cell, qvec_u, qvec_v

    # ── :edge trajectory (Jacobi update: read old flux_edge throughout) ──
    flux_edge_new = copy(flux_edge)
    for e in 1:ne
        is_diamond_valid[e] || continue
        q_prev_eff = _q_centred(flux_edge, e, edges.collinear_i[e], edges.collinear_j[e], theta)
        _, q_stored = _bates_flux_limited_corrected(q_prev_eff, h_flow_full[e], dWSE_full[e],
                                                      edges.width[e], edges.L[e], n_mann4, dt4,
                                                      depth_donor_full[e])
        flux_edge_new[e] = q_stored
    end
    flux_edge = flux_edge_new
    hist_edge[step, :] = flux_edge

    # ── :cell trajectory (Jacobi update via qvec from PREVIOUS step) ──
    flux_cell_new = copy(flux_cell)
    for e in 1:ne
        is_diamond_valid[e] || continue
        ci, cj = edges.cell_i[e], edges.cell_j[e]
        q_prev_eff = -0.5 * ((qvec_u[ci] + qvec_u[cj]) * edges.nf_x[e] +
                              (qvec_v[ci] + qvec_v[cj]) * edges.nf_y[e])
        _, q_stored = _bates_flux_limited_corrected(q_prev_eff, h_flow_full[e], dWSE_full[e],
                                                      edges.width[e], edges.L[e], n_mann4, dt4,
                                                      depth_donor_full[e])
        flux_cell_new[e] = q_stored
    end
    flux_cell = flux_cell_new
    hist_cell[step, :] = flux_cell

    # ── Phase F: reconstruct qvec from the just-updated flux_cell ──
    qvec_u_new = zeros(n)
    qvec_v_new = zeros(n)
    for i in 1:n
        qu = 0.0; qv = 0.0
        for s in 1:N_SIDES
            e = state.cell_edge_index[s, i]
            e == 0 && continue
            obs = -flux_cell[e]
            qu += state.mom_weights[s,     i] * obs
            qv += state.mom_weights[5 + s, i] * obs
        end
        qvec_u_new[i] = qu
        qvec_v_new[i] = qv
    end
    qvec_u = qvec_u_new
    qvec_v = qvec_v_new
end

live4 = findall(is_diamond_valid)
println("  $n_steps virtual steps, $(length(live4))/$ne edges updated (diamond-valid only)")

# ── convergence: std of the last 20 steps, per edge, per model ──
window = min(20, n_steps)
conv_edge = [std(hist_edge[end-window+1:end, e]) for e in live4]
conv_cell = [std(hist_cell[end-window+1:end, e]) for e in live4]
@printf("
  Convergence (std of last %d steps, should -> 0 if settled to a fixed point):\n", window)
@printf("    :edge  max=%.3e  mean=%.3e  (n edges with std>1e-6: %d)\n",
        maximum(conv_edge), mean(conv_edge), count(>(1e-6), conv_edge))
@printf("    :cell  max=%.3e  mean=%.3e  (n edges with std>1e-6: %d)\n",
        maximum(conv_cell), mean(conv_cell), count(>(1e-6), conv_cell))
println("  A large max here for :edge but small for :cell is direct evidence")
println("  cell-vector momentum suppresses persistent per-edge oscillation")
println("  that edge-scalar storage does not -- the checkerboard-suppression")
println("  mechanism this whole reconstruction effort is ultimately about.")

# ── final converged flux: edge vs cell, and vs the single-shot Test 2/3 answer ──
final_edge = flux_edge[live4]
final_cell = flux_cell[live4]
diff_ec = final_edge .- final_cell

println("
--- Final (step $n_steps) converged flux: :edge vs :cell ---")
@printf("  :edge final: mean=%.3e  std=%.3e  max|.|=%.3e\n", mean(final_edge), std(final_edge), maximum(abs.(final_edge)))
@printf("  :cell final: mean=%.3e  std=%.3e  max|.|=%.3e\n", mean(final_cell), std(final_cell), maximum(abs.(final_cell)))
@printf("  diff (edge-cell): mean=%.3e  std=%.3e  max|diff|=%.3e\n", mean(diff_ec), std(diff_ec), maximum(abs.(diff_ec)))

# how far did repeated feedback move each model from the naive single-shot answer?
# (live4 == valid_idx by construction -- both are findall(is_diamond_valid)
# against the same unmodified array -- so Q_lim_dia[k] already corresponds
# to live4[k] directly; no re-lookup needed)
q_singleshot_dia = Q_lim_dia ./ [edges.width[e] for e in live4]
shift_edge = final_edge .- q_singleshot_dia
shift_cell = final_cell .- q_singleshot_dia
println("
--- How far momentum feedback shifted the result vs a single q_prev=0 evaluation ---")
@printf("  :edge shift from single-shot: mean=%.3e  std=%.3e  max|.|=%.3e\n",
        mean(shift_edge), std(shift_edge), maximum(abs.(shift_edge)))
@printf("  :cell shift from single-shot: mean=%.3e  std=%.3e  max|.|=%.3e\n",
        mean(shift_cell), std(shift_cell), maximum(abs.(shift_cell)))
println("  Larger shift = more distortion introduced by repeated momentum")
println("  feedback under a fixed driving force -- more room for that")
println("  feedback to compound into an orientation-dependent bias over")
println("  many real timesteps.")

# ── bearing correlation of the converged state, each model ──
bea4 = bearing_deg[live4]
cos2b4 = cosd.(2 .* bea4)
println("
--- Bearing correlation of converged flux, each model ---")
@printf("  corr(:edge final, cos(2*bearing)) = %.4f\n", cor(final_edge, cos2b4))
@printf("  corr(:cell final, cos(2*bearing)) = %.4f\n", cor(final_cell, cos2b4))

println("
Done.")
