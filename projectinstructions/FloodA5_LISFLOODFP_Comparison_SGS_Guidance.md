# FloodA5 — LISFLOOD-FP ACC Comparison & SGS Implementation Guidance

**Branch context:** `standard_flow_fixes` merged to `main`  
**Reference code:** LISFLOOD-FP (C++, `fp_acc.cpp`, `sgc.cpp`, `iterateq.cpp`)  
**Date:** June 2026

---

## 1. Purpose

This document records the findings of a line-by-line comparison between FloodA5 and the LISFLOOD-FP acceleration (ACC) formulation, covering both standard (floodplain) flow and the sub-grid channel (SGC) model. It is intended as complete, self-contained guidance for implementing the outstanding SGS fixes in a fresh conversation.

The standard flow fixes applied in `standard_flow_fixes` (wave-speed CFL, Froude limiter, volume limiter, consistent `q_prev`) are confirmed against LISFLOOD-FP below. The SGS section provides a precise implementation plan informed by the LISFLOOD-FP SGC solver, which is architecturally different from FloodA5's SGS but illuminates the key design choices.

---

## 2. LISFLOOD-FP ACC Architecture Summary

### 2.1 Call sequence (ACC mode, `acceleration == ON`)

```
IterateQ()
  └── FloodplainQ()
        ├── CalcFPQxAcc()   ← flux kernel, x-direction
        └── CalcFPQyAcc()   ← flux kernel, y-direction
  └── UpdateH()             ← depth update from flux divergence
  └── UpdateQs()            ← copy Qx/Qy → Qxold/Qyold (q_prev update)
  └── CalcT()               ← CFL timestep for next step
```

The key structural point: **`UpdateQs` is called after `UpdateH`**, meaning `q_prev` is updated with the flux values that were actually used in the depth update of that step. This is the mechanism that achieves consistent `q_prev` in LISFLOOD-FP — not a per-edge limiter, but a disciplined end-of-step copy.

### 2.2 SGC call sequence (SGC mode)

```
IterateQ()
  └── SGC_FloodplainQ()
        ├── CalcFPQxSGC()   ← combined SGC + FP flux kernel, x-direction
        └── CalcFPQySGC()   ← combined SGC + FP flux kernel, y-direction
  └── SGC_UpdateH()         ← volume-based depth update (hypsometric inverse)
  └── CalcT()               ← CFL using updated max depth
```

Note: in the LISFLOOD-FP SGC mode there is no separate `UpdateQs` step — both `Qxold`/`QxSGold` and their new values are written directly within `CalcFPQxSGC`/`CalcFPQySGC`.

---

## 3. Standard Flow (ACC): Line-by-Line Verification

### 3.1 Flux kernel — `CalcFPQxAcc` / `CalcFPQyAcc`

**LISFLOOD-FP:**
```cpp
hflow = getmax(z0+h0, z1+h1) - getmax(z0, z1);
hflow = getmax(hflow, 0);
hflow = getmin(hflow, Solverptr->MaxHflow);    // upper cap (rare)

Sf = -dh / Parptr->dx;                        // dh = y0 - y1 (signed WSE diff)

// Q-centred scheme (default):
Q = ((theta*q0 + 0.5*(1-theta)*(qup+qdown)) - g*dt*hflow*Sf)
  / (1 + g*dt*hflow*fn*fn*fabs(qvect) / pow(hflow, 10.0/3.0))
  * dx;

// Fallback to Bates semi-implicit if centred scheme reverses flow direction:
if (Q * dh < 0.0) {
    Q = (q0 - g*dt*hflow*Sf)
      / (1 + g*dt*hflow*fn*fn*fabs(qvect) / pow(hflow, 10.0/3.0))
      * dx;
}
```

**Key observations:**

**LISFLOOD-FP uses `10.0/3.0` as a floating-point literal** — this is correct and avoids the integer division bug found in CAESAR-Lisflood (`10/3 = 3`). FloodA5 also uses the correct form (`h_flow^(10.0/3.0)`). ✅

**LISFLOOD-FP uses a Q-centred (theta) scheme by default**, where `qup` and `qdown` are the upstream and downstream neighbour fluxes along the same direction. This introduces spatial averaging of the inertial term and provides additional damping of the checkerboard mode. When the centred scheme would reverse the flow direction (`Q*dh < 0`), it falls back to the standard Bates semi-implicit scheme. **FloodA5 does not implement the Q-centred scheme — it uses only the Bates semi-implicit form.** This is a valid simplification but means FloodA5 lacks the extra damping that the centred scheme provides. This is unlikely to be a problem given that the Froude limiter and consistent `q_prev` (now applied) address the checkerboard instability directly.

**LISFLOOD-FP uses `qvect` rather than `|q0|` in the denominator friction term when `fricSolver2D == ON`.** In 2D mode, `qvect = sqrt(qx² + qy_avg²)` — the magnitude of the full velocity vector rather than just the x-component. This more accurately represents friction on a 2D flow field. FloodA5 uses `|q_prev|` for the edge in question, which is the 1D approximation. On the A5 pentagonal grid the 2D extension would require averaging `q` from all five edges of the two adjacent cells, which is feasible but not currently implemented. This is a lower-priority enhancement.

**Manning's n:** LISFLOOD-FP uses `fn = 0.5 * (n[p0] + n[p1])` — arithmetic mean of the two adjacent cells. FloodA5 uses `min(n[ci], n[cj])`. The mean is the more standard choice and avoids the slight over-conductivity of using the minimum. **This is a small discrepancy worth noting for the CAESAR benchmark but is not expected to affect stability.**

**Dry-edge threshold:** LISFLOOD-FP tests `hflow > Solverptr->DepthThresh` before computing any flux. `DepthThresh` is a user-configurable parameter, typically 1e-3 m (0.001 m). This is equivalent to FloodA5's `HFLOW_THRESHOLD = 0.001 m` which was added in `standard_flow_fixes`. ✅

**No explicit Froude limiter in LISFLOOD-FP ACC.** There is no equivalent of the `FROUDE_LIMIT = 0.8` cap from CAESAR-Lisflood in the LISFLOOD-FP ACC kernel. The Q-centred scheme and the flow-reversal correction (`if (Q*dh < 0)`) provide equivalent stability by different means. FloodA5's Froude limiter is a valid and more direct approach to the same problem. **Retain the Froude limiter in FloodA5.**

**`UpdateQs` and consistent `q_prev`:** After `UpdateH`, LISFLOOD-FP calls `UpdateQs()` which copies `Qx → Qxold` and `Qy → Qyold` (dividing by `dx` to convert from m³/s to m²/s). This means `q_prev` for the next step is exactly the flux that was computed and used in the current step — not an uncapped value. This is the LISFLOOD-FP mechanism for Fix C. FloodA5's implementation of Fix C (writing `q_stored` — the post-limiting unit discharge — directly to `edges.flux[e]` in Phase A) is equivalent in effect. ✅

### 3.2 CFL timestep — `CalcT`

**LISFLOOD-FP:**
```cpp
NUMERIC_TYPE cfl = Solverptr->cfl;   // default 0.7
MH = CalcMaxH(Parptr, Arrptr);       // global maximum water depth
if (MH > Solverptr->DepthThresh) {
    locT = cfl * Parptr->dx / sqrt(g * MH);
    Solverptr->Tstep = getmin(Solverptr->Tstep, locT);
}
else {
    Solverptr->Tstep = Solverptr->InitTstep;
}
```

This is exactly the wave-speed Courant formula: `dt = C × dx / √(g × h_max)` with `C = 0.7`. FloodA5's updated `_cfl_dt` uses the same formula with `dx_min = sqrt(min(cell_area))` — appropriate for variable-size A5 cells. ✅

The commented-out per-cell loop in `CalcT` shows an older approach that LISFLOOD-FP considered and rejected, matching the decision to use the global maximum depth for computational efficiency.

### 3.3 Volume/depth update — `UpdateH`

**LISFLOOD-FP:**
```cpp
dV = Tstep * (Qx[i,j] - Qx[i+1,j] + Qy[i,j] - Qy[i,j+1]);
H[i,j] += dV / dA;
if (H[i,j] < 0) H[i,j] = 0;
```

Standard finite-volume divergence on a staggered grid. No per-edge donor limiter — the staggered grid and Q-centred scheme provide sufficient stability without one. FloodA5's Phase B donor limiter (`DONOR_EDGE_DIVISOR`) fills the equivalent role given that FloodA5 does not use a staggered grid or Q-centred scheme. ✅

### 3.4 Standard flow summary: all fixes confirmed

| Feature | LISFLOOD-FP | FloodA5 (post-fixes) | Status |
|---|---|---|---|
| Bates eq. 9 exponent `10/3` | `10.0/3.0` ✓ | `10.0/3.0` ✓ | ✅ Match |
| CFL formula | Wave-speed, C=0.7 | Wave-speed, C=0.7 | ✅ Match |
| Dry-edge threshold | `DepthThresh` ~0.001 m | `HFLOW_THRESHOLD = 0.001 m` | ✅ Match |
| `q_prev` consistency | `UpdateQs` end-of-step copy | `q_stored` written in Phase A | ✅ Equivalent |
| Froude limiter | No (uses Q-centred scheme instead) | Yes, Fr ≤ 0.8 | ✅ Alternative approach |
| Manning's n per edge | Mean of adjacent cells | Min of adjacent cells | ⚠️ Minor discrepancy |
| Q-centred scheme | Yes (with semi-implicit fallback) | No | ℹ️ Not needed given Froude limiter |
| 2D friction vector | Optional (`fricSolver2D`) | No (1D per-edge) | ℹ️ Future enhancement |

---

## 4. SGC (Sub-Grid Channel) Model: Architecture Comparison

The LISFLOOD-FP SGC model and FloodA5 SGS are conceptually similar — both resolve sub-cell hydraulic geometry — but differ structurally. Understanding the differences is essential for implementing the SGS fixes correctly.

### 4.1 Architectural differences

| Aspect | LISFLOOD-FP SGC | FloodA5 SGS |
|---|---|---|
| Sub-cell geometry | Parametric cross-section (rectangular, power-law, parabolic, trapezoidal) stored per cell | Hypsometric curve (elevation vs. volume/area) sampled from LiDAR DEM |
| Channel bed | `SGCz` (derived from hydraulic geometry or file) | Pre-computed from DEM samples; stored as `sgs_elev_bins` |
| WSE from volume | `CalcSGC_UpH()` — analytical inverse (power law etc.) | Inverse interpolation of `vol_curve` (`wse_from_volume`) |
| Flow area | `CalcSGC_A()` — parametric | Inferred from hypsometric curve |
| Hydraulic radius | `CalcSGC_R()` — parametric | Not computed; uses `h_flow` cap instead (Bug 49) |
| Floodplain fraction | Separate FP flux computed when `we < cell_width` | Single unified flux per edge |
| SGC flux type | Volumetric (m³/s), uses cross-sectional area A in numerator | Unit discharge (m²/s), uses `h_flow` |

### 4.2 LISFLOOD-FP SGC flux kernel (`CalcFPQxSGC`)

The SGC flux in LISFLOOD-FP is fundamentally different from the floodplain Bates formulation:

```cpp
// SGC flux (volumetric, m³/s):
// Uses cross-sectional flow area A and hydraulic radius R
Sf = -dh / (dx * m);          // meander-corrected gradient
CalcSGC_A(gr, hflow, bf, &A, &w, SGCptr);
R  = CalcSGC_R(gr, hflow, bf, w, width, A, SGCptr);

QxSGold[pq0] = (qc - g*A*dt*Sf)
             / (1 + dt*g*cn*fabs(qc) / (pow(R, 4.0/3.0) * A));
```

This is **not** the Bates eq. 9 form — it is derived directly from the inertial shallow-water equations using hydraulic radius R and cross-sectional area A rather than flow depth h. The denominator uses `R^(4/3)` (Manning's hydraulic radius form) rather than `h^(10/3)` (Bates wide-channel approximation). `cn` is Manning's `n²` (pre-squared).

The FP component in SGC cells is computed separately using standard Bates eq. 9 applied only to the overbank depth `h - bankfull_depth`, over width `dx - channel_width`. The two fluxes are then summed: `Q = Qc + Q_fp`.

**FloodA5 SGS uses a single unified `_bates_flux` call per edge** with a modified sill (`z_sill_eff`) and `h_flow` cap. This is simpler than LISFLOOD-FP's split approach but means the hydraulic radius correction and the parametric cross-section are not represented. The hypsometric approach compensates for this by capturing sub-cell geometry statistically from LiDAR rather than analytically from cross-section parameters.

### 4.3 `q_prev` in LISFLOOD-FP SGC

In `CalcFPQxSGC`, the SGC flux is written directly to `QxSGold[pq0]` and the FP flux to `Qxold[pq0]` within the same function call. These are both read as `q_prev` on the next call. There is no separate `UpdateQs` step in SGC mode — the update is in-place within the flux kernel.

This means the stored `q_prev` in LISFLOOD-FP SGC is always consistent with what was computed, because there is no separate limiting step that diverges from the stored value.

---

## 5. SGS Implementation Plan

This section provides complete, actionable guidance for implementing Fixes A and C in `step_sgs!`. Fix B (volume limiter) is handled by the existing `DONOR_EDGE_DIVISOR` cap as explained below.

### 5.1 Fix A — Froude Limiter (apply inline in `step_sgs!` Phase A)

**Rationale:** The Froude cap is a function of `h_flow`, `g`, and `q` only. In SGS, `h_flow` has already been modified by Bug 49 Fix B: it is capped at `max(depth_ci, depth_cj)` before being passed to `_bates_flux`, making it equivalent to the actual water depth at the edge. The Froude cap applied after the Bates call uses this same effective `h_flow`, so it is physically correct.

LISFLOOD-FP ACC does not have an explicit Froude cap — it achieves the same effect via the Q-centred scheme. LISFLOOD-FP SGC does not have one either; it instead uses the `A` and `R` formulation which naturally constrains unrealistic fluxes. Since FloodA5 SGS does not use either of those mechanisms, the explicit Froude cap is warranted.

**Implementation:** After the `_bates_flux` call in SGS Phase A, extract `q_new = Q / width`, apply the Froude cap, then recompute `Q`.

```julia
# After existing _bates_flux call in step_sgs! Phase A:
Q_raw = _bates_flux(edges.flux[e], wse_eff_i, wse_eff_j, z_sill_eff,
                    edges.width[e], edges.L[e], edges.cos_theta[e],
                    min(state.manning_n[ci], state.manning_n[cj]), dt)

q_raw = Q_raw / max(edges.width[e], 1e-6)

# Fix A: Froude limiter using h_flow_eff (already computed above for z_sill_eff)
q_max    = h_flow_eff * sqrt(_G * h_flow_eff) * FROUDE_LIMIT
q_stored = clamp(q_raw, -q_max, q_max)

Q = q_stored * edges.width[e]
```

Note: `h_flow_eff` is the capped `h_flow` already computed in Phase A for the `z_sill_eff` calculation (Bug 49 Fix B). It is already available at the point where the Froude limiter is applied.

### 5.2 Fix B — Volume Limiter (retain existing Phase B cap; no change needed)

**Rationale:** The CAESAR-Lisflood volume limiter uses `depth_donor × width / (5 × dt)`. In the SGS context, `water_depth` is the hypsometric depth (above `z_min`), not a meaningful measure of the volume available to flux at a given edge — the actual movable water depends on the shape of the hypsometric curve below the current WSE. A depth-based cap can over-restrict flow in cells where most of the volume is concentrated at low elevations (e.g. a deep narrow channel with a high bankfull level).

LISFLOOD-FP SGC does not apply a volume limiter inside the flux kernel. The LISFLOOD-FP SGC volume update uses the full computed flux and relies on the parametric `A` and `R` formulation to produce physically bounded discharges. Since FloodA5 SGS uses the Bates formulation rather than the R-A formulation, some volume control is still needed, but the existing `DONOR_EDGE_DIVISOR = 10` cap in Phase B (`volume[donor] / 10` per edge, guaranteeing ≤ 50% total drain across all 5 edges) is the appropriate mechanism. It operates on actual stored volume rather than depth and is therefore more appropriate for hypsometric cells.

**Action: no change to Phase B for SGS.**

### 5.3 Fix C — Consistent `q_stored` (apply inline in `step_sgs!` Phase A; highest priority)

**Rationale:** This is the structural analogue of what LISFLOOD-FP SGC achieves by writing `QxSGold` within the flux kernel itself. In LISFLOOD-FP the stored `q_prev` is always the value that was actually used, because storage and computation happen in the same function. In FloodA5 SGS the raw Bates `q` was stored before the Phase B limiter could modify it, creating the divergence.

After Fix A is applied, `q_stored` is the post-Froude-limited unit discharge. This should be written to `edges.flux[e]` before Phase B runs. The Phase B `DONOR_EDGE_DIVISOR` cap then operates on `edge_vol[e] = q_stored × width × dt`, not on an unlimited value, but if Phase B does clip the volume further (which it should rarely do after Fix A), the stored `q` will still be slightly inconsistent. This residual inconsistency is small and acceptable — the primary source of divergence (the unlimited Bates `q`) is eliminated.

For a fully consistent implementation matching LISFLOOD-FP SGC behaviour, the `DONOR_EDGE_DIVISOR` cap in Phase B could also update `edges.flux[e]` after clipping `ev`. This would be Fix C full implementation:

```julia
# Phase B (SGS) — after donor cap:
if ev > 0.0
    ev_capped = min(ev, state.volume[cj] / DONOR_EDGE_DIVISOR)
    if ev_capped < ev
        # Update stored q to reflect the further cap
        edges.flux[e] = ev_capped / (max(edges.width[e], 1e-6) * dt)
    end
    ev = ev_capped
else
    ev_capped = max(ev, -state.volume[ci] / DONOR_EDGE_DIVISOR)
    if ev_capped > ev
        edges.flux[e] = ev_capped / (max(edges.width[e], 1e-6) * dt)
    end
    ev = ev_capped
end
```

This full implementation is recommended but can be deferred. The primary fix — writing `q_stored` (post-Froude) rather than `q_raw` — should be implemented first and validated before adding the Phase B correction.

### 5.4 `_bates_flux` return value refactor

`_bates_flux` currently returns only `Q` (volumetric flux, m³/s). To implement Fix C cleanly, it should also return `q_new` (unit discharge, m²/s). The options are:

**Option A — Derive `q_new` from `Q` in the caller:**
```julia
Q_raw  = _bates_flux(...)
q_raw  = Q_raw / max(edges.width[e], 1e-6)
```
Simple, no function signature change. Slightly wasteful (multiply then divide by `width`), but negligible.

**Option B — Modify `_bates_flux` to return `(Q, q_new)`:**
```julia
@inline function _bates_flux(...)::Tuple{Float64, Float64}
    # ... existing logic ...
    return (q_new * width, q_new)
end
```
Cleaner long-term. Aligns `_bates_flux` with `_bates_flux_limited`. Requires updating the SGS call site and any other callers.

**Recommendation: Option A for the initial SGS fix to minimise diff size and risk. Migrate to Option B in a subsequent cleanup.**

### 5.5 Complete Phase A change for `step_sgs!`

The following shows the complete modification to the inner body of the Phase A thread loop in `step_sgs!`, replacing the existing `_bates_flux` call and `edges.flux[e]` write:

```julia
# ── existing SGS pre-processing (unchanged) ────────────────────────────────
wse_eff_i = ...   # dry-cell effective WSE (Bug 48 fix)
wse_eff_j = ...
z_sill_eff = ...  # edge sill capped at max(wse) - max(depth) (Bug 49 fix)
h_flow_eff = max(wse_eff_i, wse_eff_j) - z_sill_eff   # already computed

# ── Bates flux (unchanged call signature) ──────────────────────────────────
Q_raw = _bates_flux(edges.flux[e], wse_eff_i, wse_eff_j, z_sill_eff,
                    edges.width[e], edges.L[e], edges.cos_theta[e],
                    min(state.manning_n[ci], state.manning_n[cj]), dt)

# ── Fix A: Froude limiter ──────────────────────────────────────────────────
q_raw    = Q_raw / max(edges.width[e], 1e-6)
q_max    = h_flow_eff * sqrt(_G * h_flow_eff) * FROUDE_LIMIT
q_stored = clamp(q_raw, -q_max, q_max)

# ── Fix C: write post-limiting q as q_prev for next step ──────────────────
edges.flux[e] = q_stored          # was: edges.flux[e] = Q_raw / max(edges.width[e], 1e-6)
edge_vol[e]   = q_stored * edges.width[e] * dt
```

### 5.6 Dry-edge handling for SGS

When `h_flow_eff <= HFLOW_THRESHOLD`, `_bates_flux` already returns `0.0` (updated in `standard_flow_fixes`). In this case `Q_raw = 0`, `q_raw = 0`, `q_stored = 0`, and `edges.flux[e] = 0`. This correctly clears stale momentum on dry edges, equivalent to LISFLOOD-FP's explicit `QxSGold[pq0] = 0` when `hflow <= DepthThresh`.

No additional handling is required for the dry-edge case.

---

## 6. Additional Observations for Future Work

### 6.1 Manning's n — mean vs min

FloodA5 uses `min(manning_n[ci], manning_n[cj])` per edge. LISFLOOD-FP uses the arithmetic mean `0.5 * (n[p0] + n[p1])`. The mean is the standard choice in shallow-water modelling (it represents equal weighting of resistance from both sides of the edge). The minimum is slightly non-standard and makes the edge marginally more conductive. **Recommended: change to arithmetic mean for better alignment with LISFLOOD-FP and the literature. Low priority but worth noting before the formal benchmark.**

### 6.2 Q-centred scheme (future enhancement)

The LISFLOOD-FP Q-centred scheme averages `q_prev` with the upstream and downstream neighbours along the same flow direction:

```
q_eff = theta × q0 + 0.5 × (1 - theta) × (q_upstream + q_downstream)
```

with `theta ≈ 0.9` (default). This introduces spatial smoothing of the momentum field that further damps checkerboarding. It requires access to the adjacent edge fluxes along each edge direction, which on the A5 pentagonal mesh means identifying the two edges of `ci` and `cj` most collinear with the current edge — geometrically non-trivial but feasible using `cos_theta` as a proxy for edge alignment. This is a Phase 5 enhancement and is not needed given the current Froude limiter.

### 6.3 2D friction vector (future enhancement)

LISFLOOD-FP's `fricSolver2D` mode uses `qvect = sqrt(qx² + qy_avg²)` in the denominator, where `qy_avg` is averaged from the four surrounding y-direction edge fluxes. On the A5 mesh the equivalent would be averaging the four non-current edge fluxes at the shared cell face. This better represents friction on a 2D flow field and would be most beneficial in areas with strong cross-flow (e.g. flow around obstacles). Phase 5 enhancement.

### 6.4 LISFLOOD-FP SGC not directly portable to FloodA5

The LISFLOOD-FP SGC uses parametric cross-sections (rectangular, power-law, parabolic) with analytically computed area `A` and hydraulic radius `R`, and a `R^(4/3)` friction term in the denominator. FloodA5 SGS uses LiDAR-derived hypsometric curves and the Bates `h^(10/3)` form. These are different formulations of the same physical problem. **The LISFLOOD-FP SGC code should not be directly ported** — instead, the Fixes A and C should be applied to the existing FloodA5 SGS formulation as described in §5.

---

## 7. Action Summary for Next Session

Listed in priority order:

| # | Action | Location | Section |
|---|---|---|---|
| 1 | Apply Fix C (consistent `q_stored`) to `step_sgs!` Phase A | `FloodModel.jl` | §5.3, §5.5 |
| 2 | Apply Fix A (Froude limiter) to `step_sgs!` Phase A | `FloodModel.jl` | §5.1, §5.5 |
| 3 | Validate SGS on Carlisle domain — confirm oscillations eliminated | Test | — |
| 4 | Run SGS with realistic uniform rainfall (50 mm/hr) — confirm stability | Test | — |
| 5 | Change Manning's n per edge from `min` to mean (`0.5*(n_ci + n_cj)`) | `FloodModel.jl` | §6.1 |
| 6 | Refactor `_bates_flux` to return `(Q, q_new)` (Option B) | `FloodModel.jl` | §5.4 |
| 7 | Consider Phase B `edges.flux[e]` update after `DONOR_EDGE_DIVISOR` cap | `FloodModel.jl` | §5.3 |

Items 1–4 should be completed and validated before proceeding to 5–7.

---

## 8. Key Constants and Parameters (Reference)

| Constant | Value | Set in | Applies to |
|---|---|---|---|
| `HFLOW_THRESHOLD` | 0.001 m | `FloodModel.jl` | Both standard and SGS |
| `FROUDE_LIMIT` | 0.8 | `FloodModel.jl` | Standard flow (applied); SGS (pending) |
| `DONOR_EDGE_DIVISOR` | 10 (= 2 × N_SIDES) | `FloodModel.jl` | Both (Phase B last-resort cap) |
| `N_SIDES` | 5 | `FloodModel.jl` | A5 pentagon topology |
| CFL Courant number | 0.7 | `_cfl_dt` | Both |
| `DepthThresh` (LISFLOOD-FP equiv.) | `HFLOW_THRESHOLD` | — | Both |
| `theta` (Q-centred, LISFLOOD-FP) | 0.9 | Not in FloodA5 | Not applicable |

---

*FloodA5 — University of Canterbury | Prepared for SGS implementation session | June 2026*
