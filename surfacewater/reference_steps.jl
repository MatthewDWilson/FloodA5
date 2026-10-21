# surfacewater/reference_steps.jl
# --------------------------------
# WP4 stage 2 (2026-10): frozen copies of step_standard! and step_sgs! as they
# were at the end of WP4 stage 1 (main @ b2482c1), kept as the reference that
# the restructured step functions (and, later, the GPU kernels) must match bit
# for bit — see test/test_wp4_equivalence.jl. Do not optimise these; change
# them only together with a deliberate, documented change of the physics.

function step_standard_reference!(state::FlowState, dt::Float64; velocity::Bool = true)
    n     = length(state.cell_ids)
    edges = state.edges
    ne    = edges.n_edges

    # ── Debug: log NaN/depth state before flux loop (first 2 calls) ───────
    if _step_debug_count[] < 2
        _step_debug_count[] += 1
        k = _step_debug_count[]
        n_nan_vol   = count(isnan, state.volume)
        n_nan_elev  = count(isnan, state.elevation)
        n_nan_area  = count(isnan, state.cell_area)
        n_zero_area = count(a -> a < 1.0, state.cell_area)
        vol_sum     = sum(v for v in state.volume if isfinite(v); init=0.0)
        max_depth   = maximum(state.water_depth; init=0.0)
        @info "step_standard! call $k: NaN_vol=$n_nan_vol  NaN_elev=$n_nan_elev" *
              " NaN_area=$n_nan_area  zero_area=$n_zero_area" *
              " vol_sum=$(round(vol_sum,sigdigits=4))" *
              "  max_water_depth=$(round(max_depth,sigdigits=4))"
        if n_nan_vol > 0
            bad = findall(isnan, state.volume)
            @warn "  NaN volumes at indices: $(bad[1:min(5,end)])"
        end
        # Sample up to the first 3 cells plus the max-volume cell, for context
        # on small synthetic meshes (unit tests can have as few as 1-2 cells).
        # unique() avoids printing the same cell twice when argmax coincides
        # with one of the first three, or when n < 3.
        sample_idxs = unique(filter(i -> 1 <= i <= n, [1, 2, 3, argmax(state.volume)]))
        for i in sample_idxs
            wse = state.elevation[i] + state.volume[i] / max(state.cell_area[i], 1.0)
            @info "  cell[$i]: vol=$(round(state.volume[i],sigdigits=4))" *
                  "  depth=$(round(state.water_depth[i],sigdigits=4))" *
                  "  wse=$(round(wse,sigdigits=6))"
        end
    end

    # ── Phase A: parallel edge flux computation ─────────────────────────────────────
    # Each iteration reads cell state (read-only) and writes to unique edge
    # slots edge_vol[e] and edges.flux[e].  No two threads share a write
    # target, and (since WP0, 2026-10) Q-centring reads the read-only
    # snapshot `flux_prev`, not the live edges.flux being written here —
    # previously that read was a data race.
    #
    # _bates_flux_limited applies three stability fixes vs the earlier call:
    #   • Froude limiter      (CAESAR froude_limit=0.8: caps supercritical q)
    #   • Volume limiter      (CAESAR depth/5 cap: ≤ 1/5 donor depth per edge)
    #   • Consistent q_prev  (stores post-limiting q, not raw Bates q)
    # See _bates_flux_limited docstring for full rationale.
    #
    # ── WLSQ non-orthogonal gradient correction (flow-direction-fixes) ─────
    # When state.gradient_correction is true, _bates_flux_limited_corrected
    # is used instead of _bates_flux_limited, with the WSE gradient
    # reconstructed once (serially, below) before the threaded edge loop —
    # _compute_wse_gradients! needs the full per-cell WSE array, which the
    # legacy code path only ever computed inline per-edge, not as a
    # standalone vector. Materialising `wse` here costs O(n_cells), trivial
    # next to the edge loop itself, and is skipped entirely when
    # gradient_correction is true by default; use --gradient-correction off
    # to select the legacy uncorrected path for benchmarking.
    # See FloodA5_NonOrthogonal_Correction_Plan.md §5.2, §5.5.
    ws = _workspace!(state)
    # Bound unconditionally (WP4): a closure capturing a conditionally
    # assigned local boxes it, which made every Phase A access type-unstable.
    wse_all = ws.wse
    if state.gradient_correction
        @inbounds for i in 1:n
            wse_all[i] = state.elevation[i] + state.volume[i] / max(state.cell_area[i], 1.0)
        end
        _compute_wse_gradients!(state, wse_all)
    end
    # Diamond vertex values once per vertex (WP4; previously two Dict lookups
    # and a dot product per edge, twice for diamond-wetdry).
    dtab = state.diamond_table
    use_diamond = state.gradient_correction && dtab !== nothing &&
                  state.face_flux_method in (:diamond, :diamond_wetdry)
    veta = ws.veta; vwet = ws.vwet
    if use_diamond
        diamond_vertex_eta!(veta, dtab, wse_all)
        state.face_flux_method === :diamond_wetdry &&
            diamond_vertex_wet!(vwet, dtab, state.water_depth, HFLOW_THRESHOLD)
    end

    edge_vol = ws.edge_vol

    # ── WP0 (2026-10): read-only momentum snapshot for Q-centring ──────────
    # The threaded loop below overwrites edges.flux[e]; _q_centred reads the
    # flux of OTHER (collinear) edges. Reading the live array was a data race
    # (results depended on edge order and thread scheduling). Q-centring now
    # reads this start-of-step copy. Under :cell momentum, Phase A never reads
    # edges.flux, and under :blend (WP2) it reads only its own edge's flux
    # before writing it, so no copy is needed and the alias is harmless.
    flux_prev = state.momentum_model === :edge ? copyto!(ws.flux_prev, edges.flux) : edges.flux

    # WP2 (2026-10): with --limiter positivity the per-edge volume limiter
    # (Fix B) is disabled — the kernels skip it when depth_donor ≤ 0 — and
    # positivity is enforced per donor cell in Phase B instead. --limiter
    # hybrid keeps Fix B and replaces only the per-edge donor cap.
    positivity = state.limiter !== :legacy
    keep_fixB  = state.limiter !== :positivity

    # WP1 (2026-10): vector friction. qvec_u/qvec_v hold the cell discharge
    # vectors reconstructed by Phase F at the end of the previous step (zero
    # before the first step, which reduces to component friction).
    vector_friction = state.friction_model === :vector

    Threads.@threads for e in 1:ne
        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        # Degenerate edge guard (NaN in geometry or elevation).
        if (isnan(edges.width[e])      || isnan(edges.L[e])          ||
            isnan(edges.cos_theta[e])  || isnan(edges.sill[e])       ||
            isnan(state.elevation[ci]) || isnan(state.elevation[cj]))
            edge_vol[e] = 0.0
            edges.flux[e] = 0.0   # clear stale momentum on degenerate edges
            continue
        end

        wse_ci = state.elevation[ci] + state.volume[ci] / max(state.cell_area[ci], 1.0)
        wse_cj = state.elevation[cj] + state.volume[cj] / max(state.cell_area[cj], 1.0)

        # Dry edge (WP4): both flux kernels return (0, 0) when
        # max(wse) − sill ≤ HFLOW_THRESHOLD, whatever q_prev, so skip the
        # momentum/gradient work. Same test (NaN falls through) — bit-identical.
        if max(wse_ci, wse_cj) - edges.sill[e] <= HFLOW_THRESHOLD
            edges.flux[e] = 0.0
            edge_vol[e]   = 0.0
            continue
        end

        # Identify donor depth for the volume limiter (higher-WSE side donates);
        # 0.0 disables the per-edge volume limiter under --limiter positivity.
        depth_donor = !keep_fixB ? 0.0 :
                      wse_ci >= wse_cj ? state.water_depth[ci] : state.water_depth[cj]

        # ── q_prev: edge-scalar (legacy) or face projection of cell vector ──
        # When momentum_model == :cell, the face unit discharge predictor is
        # obtained by projecting the average of the two adjacent cells' vector
        # unit discharges onto the face normal (arithmetic face averaging;
        # sign negated because flux convention is opposite to n̂_f direction):
        #   q_prev_eff = −0.5·((qvec_u_ci + qvec_u_cj)·nf_x +
        #                      (qvec_v_ci + qvec_v_cj)·nf_y)
        # When momentum_model == :edge, the legacy per-edge scalar is used
        # with optional Q-centred smoothing (theta). Q-centred is disabled
        # for the cell-vector path (already spatially coupled via the vector
        # reconstruction; adding theta would double-smooth).
        #
        # When momentum_model == :blend (WP2, 2026-10), the face's own stored
        # discharge is blended with the face-normal projection of the averaged
        # cell vector, using θ = q_centre_theta:
        #   q_prev_eff = θ·q_f + (1−θ)·q_proj
        # This replaces Q-centring's "most collinear face" average, which on
        # pentagons mixes faces ~144° apart (non-parallel components) and by
        # itself deflects the steady flow (WP1 frozen-coefficient test). It is
        # mesh-general and unifies the two models in one parameter: θ = 1 is
        # pure edge momentum (no smoothing), θ = 0 is exactly :cell. On a
        # Cartesian grid q_proj = (q_w + 2q_f + q_e)/4, so it acts like
        # Q-centring there. Review §4.2.
        q_proj = state.momentum_model === :edge ? 0.0 :
                 -0.5 * ((state.qvec_u[ci] + state.qvec_u[cj]) * edges.nf_x[e] +
                         (state.qvec_v[ci] + state.qvec_v[cj]) * edges.nf_y[e])
        q_prev_eff = if state.momentum_model === :cell
            q_proj
        elseif state.momentum_model === :blend
            state.q_centre_theta * flux_prev[e] + (1.0 - state.q_centre_theta) * q_proj
        elseif state.q_centre_mode === :consistent
            _q_centred_consistent(flux_prev, e, edges, state.qvec_u, state.qvec_v,
                                  state.q_centre_theta)
        else
            _q_centred(flux_prev, e,
                       edges.collinear_i[e], edges.collinear_j[e],
                       state.q_centre_theta)
        end

        # Face-tangential discharge for vector friction: the average of the
        # two adjacent cell vectors projected onto the face tangent
        # t̂ = (−nf_y, nf_x). Only its magnitude enters (via q_n² + q_t²), so
        # the tangent's orientation does not matter.
        q_t = vector_friction ?
              0.5 * ((state.qvec_v[ci] + state.qvec_v[cj]) * edges.nf_x[e] -
                     (state.qvec_u[ci] + state.qvec_u[cj]) * edges.nf_y[e]) : 0.0

        if state.gradient_correction
            # ── Local-orthogonal + WLSQ non-orthogonal correction ───────────
            # (2026-08-18, standard-flow-momentum-instability session)
            #
            # SUPERSEDES the "purely gradient-driven" (2026-07-22) formula
            # below, which computed BOTH the dominant driving term and the
            # tangential correction from the WLSQ-reconstructed cell-centred
            # gradient (grad_wse), with no direct local WSE-difference anchor
            # anywhere in the calculation. That formula passed every
            # correctness/symmetry test available at the time (T-NOC1-5,
            # test_mirror_symmetry.jl MS1-4 — a hand-built, exactly
            # mirror-symmetric synthetic mesh confirms the formula is
            # per-step symmetric in isolation) but produced a severe,
            # dt-dependent (non-monotonic, V-shaped) north/south instability
            # on the real planar-slope mesh once combined with cell-vector
            # momentum (--momentum-model cell): two independent WLSQ
            # reconstructions over the same 5-neighbour stencil, both fed
            # back into the same flux equation every step, with nothing
            # tying either one to a locally-measured, self-correcting value.
            #
            # This is the standard OpenFOAM "orthogonal-corrected" /
            # "over-relaxed" decomposition: the dominant (orthogonal) part of
            # the face gradient is computed directly from the two adjacent
            # cell values — cheap, local, and inherently self-correcting,
            # since it reacts immediately to this specific edge's own state —
            # and only the smaller non-orthogonal remainder (the V̂ term,
            # typically |V̂| ~ 0.3-0.6 on this mesh, i.e. a minority
            # contribution) comes from the reconstructed gradient:
            #
            #   dWSE_n = c·(wse_ci - wse_cj) − alpha·L·(∇WSE_f · V̂)
            #
            # where c = edges.cos_theta[e] (oriented cosine, ≥ 0) and
            # V̂ = (edges.skew_x[e], edges.skew_y[e]). Recovers the legacy
            # kernel exactly at skew=0 (T-NOC4 still applies unchanged, since
            # the direct-difference term is bit-identical to _bates_flux's
            # dWSE = wse_i - wse_j). Should be re-validated with
            # test_mirror_symmetry.jl (still expected to pass — wse_ci-wse_cj
            # is trivially mirror-symmetric by the same construction) and
            # with the real-mesh planar-symmetry dt sweep before flipping
            # gradient_correction's default.
            #
            # ── Diamond face-flux method (Phase C, directional-bias- ─────
            # reformulation, see mesh/DiamondFlux.jl and
            # FloodA5_PhaseB_Complete_PhaseC_Handoff.md). When
            # face_flux_method == :diamond, dWSE_n is instead constructed
            # directly per edge from the B2/B3 diamond face-normal
            # derivative — no cell-averaged-then-projected gradient (gx_f/
            # gy_f below) is involved for that edge at all. This is NOT a
            # new flux kernel: FloodA5_PhaseB_B1_FaceFluxEquation.md §3.3
            # proves dWSE_n = -L·g_n is exact for ANY exact gradient
            # source, so the diamond gradient plugs directly into the same
            # _bates_flux_limited_corrected call below.
            #
            # Falls back to the legacy WLSQ+skew construction (unchanged,
            # below) whenever face_flux_method == :legacy, OR when
            # face_flux_method == :diamond but this specific edge's
            # diamond record is invalid (diamond_dWSE_n returns NaN — see
            # DiamondFluxTable's docstring for when this occurs; expected
            # to be extremely rare, ~1 edge in 74,000 on the Carlisle
            # mesh in the Phase B empirical verification).
            dWSE_n = NaN
            #
            # :diamond_wetdry (WP2, 2026-10): as :diamond, but where the
            # edge's vertex stencil includes a dry cell (whose "WSE" is just
            # its bed elevation) the face gradient is the two-point form
            # c·(wse_ci − wse_cj) instead — see diamond_stencil_wet. Removes
            # the spurious bank-side gradients seen on the Carlisle DEM.
            if use_diamond && state.face_flux_method === :diamond
                dWSE_n = diamond_dWSE_n_v(dtab, e, wse_ci, wse_cj, veta, edges.L[e])
            elseif use_diamond     # :diamond_wetdry
                dWSE_n = diamond_stencil_wet_v(dtab, e, vwet) ?
                         diamond_dWSE_n_v(dtab, e, wse_ci, wse_cj, veta, edges.L[e]) :
                         edges.cos_theta[e] * (wse_ci - wse_cj)
            end
            if isnan(dWSE_n)
                gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
                gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
                Vhat_dot = gx_f * edges.skew_x[e]  + gy_f * edges.skew_y[e]
                dWSE_n   = edges.cos_theta[e] * (wse_ci - wse_cj) -
                           state.gradient_correction_alpha * edges.L[e] * Vhat_dot
            end
            h_flow = max(wse_ci, wse_cj) - edges.sill[e]

            Q, q_stored = _bates_flux_limited_corrected(
                q_prev_eff, h_flow, dWSE_n,
                edges.width[e], edges.L[e],
                min(state.manning_n[ci], state.manning_n[cj]), dt,
                depth_donor;
                vol_limit_divisor = Float64(state.max_neighbours),
                q_t = q_t, vector_friction = vector_friction)
        else
            Q, q_stored = _bates_flux_limited(
                q_prev_eff, wse_ci, wse_cj, edges.sill[e],
                edges.width[e], edges.L[e], edges.cos_theta[e],
                min(state.manning_n[ci], state.manning_n[cj]), dt,
                depth_donor;
                vol_limit_divisor = Float64(state.max_neighbours),
                q_t = q_t, vector_friction = vector_friction)
        end

        edges.flux[e] = q_stored   # Fix C: post-limiting q → consistent q_prev
        edge_vol[e]   = Q * dt     # signed volume (m³)
    end

    # ── Phase B: serial dV scatter ─────────────────────────────────────────────
    # _bates_flux_limited already bounds the transfer via the Froude and volume
    # limiters.  The DONOR_EDGE_DIVISOR cap is retained as a last-resort
    # mass-conservation guard (should rarely trigger).
    # state.volume[] is still the start-of-step snapshot (Phase A did not
    # modify it), so the donor-limit comparison is consistent.
    #
    # WP2: under --limiter positivity the per-edge donor cap is replaced by
    # one scale factor per donor cell (all of its outflows together; written
    # back to edges.flux — Fix C). See _positivity_scale!.
    positivity && _positivity_scale!(edge_vol, edges.flux, state.volume,
                                     edges.cell_i, edges.cell_j; out = ws.out)
    dV = fill!(ws.dV, 0.0)
    @inbounds for e in 1:ne
        ev = edge_vol[e]
        iszero(ev) && continue

        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        if !positivity
            if ev > 0.0
                ev = min(ev,  state.volume[cj] / _donor_edge_divisor(state))   # cj is donor
            else
                ev = max(ev, -state.volume[ci] / _donor_edge_divisor(state))   # ci is donor
            end
        end

        dV[ci] += ev   # ci gains when ev > 0
        dV[cj] -= ev   # cj loses when ev > 0
    end

    # ── Phase C: parallel cell volume update ──────────────────────────────
    _apply_dV_standard!(state, dV)

    # ── Phase D: open boundary outflow ────────────────────────────────────
    # Apply ghost-edge flux for each non-Closed boundary edge.
    # Serial: ghost edges are per-cell; no write hazard.
    # vol_removed accumulates for mass balance accounting.
    if !isempty(state.ghost_edges)
        _apply_ghost_fluxes_standard!(state, dt)
    end

    # ── Phase E: velocity ─────────────────────────────────────────────────
    # (WP4: `velocity = false` skips it; run_simulation! then computes it only
    # for the steps that are written or pushed to the visualiser — the result
    # depends only on the end-of-step state, so the output is unchanged.)
    velocity && _compute_velocity!(state)

    # ── Phase F: cell-vector momentum reconstruction ───────────────────────
    # Active when --momentum-model cell is set (qvec is then the momentum
    # state feeding Phase A q_prev_eff on the next step — breaking the
    # per-edge directional-feedback mechanism identified in the alpha sweep,
    # 2026-08-14) or --friction-model vector (qvec supplies the face-
    # tangential discharge q_t for the friction term; WP1, 2026-10). Under
    # :edge momentum + :vector friction, edges.flux remains the momentum
    # state and qvec is used only for q_t. Under :blend (WP2) qvec supplies
    # the projected part of the blend. See _reconstruct_qvec!.
    if state.momentum_model !== :edge || state.friction_model === :vector ||
       state.q_centre_mode === :consistent
        _reconstruct_qvec!(state, false)
    end
end

function step_sgs_reference!(state::FlowState, dt::Float64; velocity::Bool = true)
    n     = length(state.cell_ids)
    edges = state.edges
    ne    = edges.n_edges

    # ── Step 0: WSE and wetted area from hypsometric lookup (serial) ──────
    # Must be serial: wse_from_volume / wetted_area_from_wse each read a
    # cell's SGS table — no shared writes, but the lookup is not @inline
    # and benefits from sequential cache access patterns.
    # Can be made @threads if profiling shows it is a bottleneck.
    # (WP4: the wetted area computed here was never used — removed.)
    ws  = _edge_slots!(state, _workspace!(state))
    wse = ws.wse
    Threads.@threads for i in 1:n                 # per-cell writes only (WP4)
        @inbounds wse[i] = wse_from_volume(state.sgs_tables[i], state.volume[i])
    end

    # ── WLSQ non-orthogonal gradient correction (flow-direction-fixes) ─────
    # Unlike step_standard!, `wse` here is already a full per-cell array
    # (Step 0 above) — no separate materialisation needed before calling
    # _compute_wse_gradients!. See FloodA5_NonOrthogonal_Correction_Plan.md
    # §5.2, §11 (SGS extension).
    if state.gradient_correction
        _compute_wse_gradients!(state, wse)
    end

    # ── Phase A: parallel edge flux computation ────────────────────────────
    # Reads: wse[] (just computed, read-only), state.elevation, manning_n,
    #        edges geometry.  Writes only to edge_vol[e] and edges.flux[e]
    #        (unique per thread — no data race).
    #
    # Dry-cell WSE correction (Bug 48):
    # wse_from_volume(V=0) returns tbl.z_min (the cell thalweg elevation).
    # For a high-elevation dry cell, z_min can be >> z_sill of an adjacent edge,
    # creating a spurious driving head that drives large flux from a dry cell.
    # The fix: for flux computation only, cap a dry cell's effective WSE at z_sill,
    # so it contributes zero head above the sill — matching the standard solver's
    # behaviour where a dry cell never drives flux.
    edge_vol = ws.edge_vol
    slot_ia  = ws.slot_i; slot_ja = ws.slot_j

    # WP0 (2026-10): read-only snapshot of flux_Q for Q-centring — see the
    # matching comment in step_standard! (same data race on edges.flux_Q).
    flux_Q_prev = copyto!(ws.flux_prev, edges.flux_Q)
    max_nb_sgs  = state.max_neighbours
    # WP1 (2026-10): vector friction — see the matching comment in step_standard!.
    vector_friction = state.friction_model === :vector

    Threads.@threads for e in 1:ne
        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        (isnan(state.elevation[ci]) || isnan(state.elevation[cj])) &&
            (edge_vol[e] = 0.0; continue)

        z_sill = edges.sill[e]
        if isnan(z_sill)
            z_sill = min(state.sgs_tables[ci].z_min, state.sgs_tables[cj].z_min)
        end

        if (isnan(edges.width[e]) || isnan(edges.L[e]) || isnan(edges.cos_theta[e]))
            edge_vol[e] = 0.0
            continue
        end

        # Bug 48 / Bug 49 refinement: effective WSE for dry cells.
        # A dry cell (volume=0) contributes zero pooling head.  Its effective WSE
        # must clear TWO physical barriers before it can receive water:
        #   1. z_sill   — the channel bed (lowest DEM along the shared edge arc)
        #   2. tbl.z_min — the basin floor (lowest terrain in the receiving cell)
        wse_ci_eff = state.volume[ci] > 0.0 ? wse[ci] :
                     max(z_sill, state.sgs_tables[ci].z_min)
        wse_cj_eff = state.volume[cj] > 0.0 ? wse[cj] :
                     max(z_sill, state.sgs_tables[cj].z_min)

        # If both cells are dry, skip
        if state.volume[ci] <= 0.0 && state.volume[cj] <= 0.0
            edge_vol[e] = 0.0
            edges.flux_Q[e] = 0.0
            continue
        end

        # ── R-A flux kernel ──────────────────────────────────────────────────
        wse_flow = max(wse_ci_eff, wse_cj_eff)

        # (WP4: slots precomputed once; one bracket search per side.)
        A_i, R_i = A5Grid.edge_area_radius(state.sgs_tables[ci], slot_ia[e], wse_flow)
        A_j, R_j = A5Grid.edge_area_radius(state.sgs_tables[cj], slot_ja[e], wse_flow)
        A_edge = 0.5 * (A_i + A_j)
        R_edge = 0.5 * (R_i + R_j)

        # Q-centred smoothing on flux_Q (volumetric momentum for R-A kernel).
        # Same collinear lookup as standard flow — the checkerboard mode exists
        # on the pentagonal mesh regardless of the flux kernel used.
        # θ = state.q_centre_theta (default 0.9; set via --q-centre-theta).
        Q_prev_eff = state.q_centre_mode === :consistent ?
            _q_centred_consistent(flux_Q_prev, e, edges, state.qvec_u, state.qvec_v,
                                  state.q_centre_theta; volumetric = true) :
            _q_centred(flux_Q_prev, e,
                       edges.collinear_i[e], edges.collinear_j[e],
                       state.q_centre_theta)

        # Face-tangential discharge for vector friction, as volumetric Q:
        # qvec is reconstructed from the face unit discharges flux_Q/width
        # (Phase F below), so the matching tangential Q is q_t·width.
        Q_t = vector_friction ?
              0.5 * ((state.qvec_v[ci] + state.qvec_v[cj]) * edges.nf_x[e] -
                     (state.qvec_u[ci] + state.qvec_u[cj]) * edges.nf_y[e]) * edges.width[e] : 0.0

        Q_new = if state.gradient_correction
            # ── Local-orthogonal + WLSQ non-orthogonal correction ───────────
            # (2026-08-18 — matches step_standard!'s corrected branch; see
            # its comment for the full rationale and the real-mesh instability
            # this supersedes.) Uses the *effective* (dry-cell corrected) WSE
            # values for both h_flow and the direct WSE-difference term, as
            # the legacy SGS kernel and the original 2026-06-24 formula both
            # did — only the tangential (V̂) correction comes from the WLSQ
            # gradient. This does NOT resolve the separate open question
            # noted below (dry-cell gradient contamination via gx_f/gy_f,
            # which are still built from the raw wse[] array) — that remains
            # the acceptance gate for this SGS path once real-mesh testing
            # resumes here.
            #
            # ⚠️  KNOWN OPEN QUESTION, not yet resolved by testing: a dry
            # cell's *raw* wse[] is tbl.z_min (wse_from_volume(V=0) returns
            # z_min — see the Bug 48 comment above), which is exactly the
            # value Bug 48's clamp exists to suppress in the direct
            # wse_ci_eff − wse_cj_eff term. Because _compute_wse_gradients!
            # is called once on the raw `wse` array for ALL cells (dry and
            # wet alike — see its docstring, "agnostic to which flow method
            # produced it"), a dry cell's gradient is built from its
            # raw/unclamped z_min, not its clamped effective WSE. This means
            # the V̂ correction term (gx_f, gy_f) could reintroduce a small
            # amount of the same spurious-head signal Bug 48 was designed to
            # eliminate, via a different path than the one Bug 48 patches.
            # MUST be checked against the synthetic DEM T0–T4 regression
            # suite (dry-cell-heavy by design — see SGS_VALIDATION_SUMMARY.md)
            # before this SGS corrected path is considered validated. If
            # T0/T1 (no spurious downstream flow before the notch sill)
            # regress with --gradient-correction on, this is the first place
            # to look.
            gx_f = 0.5 * (state.grad_wse[1, ci] + state.grad_wse[1, cj])
            gy_f = 0.5 * (state.grad_wse[2, ci] + state.grad_wse[2, cj])
            Vhat_dot = gx_f * edges.skew_x[e]  + gy_f * edges.skew_y[e]
            dWSE_n   = edges.cos_theta[e] * (wse_ci_eff - wse_cj_eff) -
                       state.gradient_correction_alpha * edges.L[e] * Vhat_dot
            h_flow = wse_flow - z_sill

            _manning_flux_ra_corrected(Q_prev_eff, h_flow, A_edge, R_edge,
                                        dWSE_n, edges.L[e],
                                        0.5 * (state.manning_n[ci] + state.manning_n[cj]),
                                        dt; Q_t = Q_t, vector_friction = vector_friction)
        else
            _manning_flux_ra(Q_prev_eff, wse_ci_eff, wse_cj_eff, z_sill,
                              A_edge, R_edge,
                              edges.L[e], edges.cos_theta[e],
                              0.5 * (state.manning_n[ci] + state.manning_n[cj]),
                              dt; Q_t = Q_t, vector_friction = vector_friction)
        end

        # Fix C: store Q (m³/s) as momentum state for next step.
        edges.flux_Q[e] = Q_new
        edge_vol[e]     = Q_new * dt
    end

    # ── Phase B: serial dV scatter with per-edge donor limiter ────────────
    # WP2: --limiter positivity|hybrid replaces the per-edge donor cap with
    # one scale factor per donor cell, written back to flux_Q (Fix C). (The
    # R-A kernel has no per-edge volume limiter, so the two are the same here.)
    positivity = state.limiter !== :legacy
    positivity && _positivity_scale!(edge_vol, edges.flux_Q, state.volume,
                                     edges.cell_i, edges.cell_j; out = ws.out)
    dV = fill!(ws.dV, 0.0)
    @inbounds for e in 1:ne
        ev = edge_vol[e]
        iszero(ev) && continue

        ci = edges.cell_i[e]
        cj = edges.cell_j[e]

        # Last-resort per-edge donor cap: V/DONOR_EDGE_DIVISOR = V/10.
        # The R-A form is self-stabilising so this should rarely bind.
        # Fix C (full): back-propagate clipped Q to flux_Q so q_prev next
        # step reflects what was actually transferred.
        if positivity
            # already scaled above
        elseif ev > 0.0
            ev_capped = min(ev, state.volume[cj] / _donor_edge_divisor(state))
            if ev_capped < ev
                edges.flux_Q[e] = ev_capped / dt
            end
            ev = ev_capped
        else
            ev_capped = max(ev, -state.volume[ci] / _donor_edge_divisor(state))
            if ev_capped > ev
                edges.flux_Q[e] = ev_capped / dt
            end
            ev = ev_capped
        end

        dV[ci] += ev
        dV[cj] -= ev
    end

    # ── Phase C: parallel cell volume update ──────────────────────────────
    _apply_dV_sgs!(state, dV)

    # ── Phase D: open boundary outflow ────────────────────────────────────
    if !isempty(state.ghost_edges)
        _apply_ghost_fluxes_sgs!(state, dt)
    end

    # ── Phase E: velocity ─────────────────────────────────────────────────
    # (WP4: `velocity = false` skips it; run_simulation! then computes it only
    # for the steps that are written or pushed to the visualiser — the result
    # depends only on the end-of-step state, so the output is unchanged.)
    velocity && _compute_velocity!(state)

    # ── Phase F: cell discharge vectors for vector friction (WP1) ─────────
    (vector_friction || state.q_centre_mode === :consistent) && _reconstruct_qvec!(state, true)
end
