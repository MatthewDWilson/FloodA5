# FloodA5 — Standard Flow Stability Fixes: Findings, Changes & SGS Roadmap

**Branch:** `standard_flow_fixes` | **Date:** June 2026

---

## 1. Background & Motivation

Checkerboarding (spatial oscillation alternating between high and low water depths on adjacent cells) was observed in real-world test runs on the Carlisle domain using the standard flow solver. This is a classical instability in explicit inertial shallow-water models and had not been apparent on the flat-terrain validation meshes, because those runs lacked the topographic gradients and sustained fluxes needed to excite the instability mode.

The investigation was conducted in two stages: first, a step-by-step audit of the FloodA5 standard flow pipeline to identify structural risk factors; then a line-by-line comparison against the CAESAR-Lisflood 2.0 reference implementation, which uses the same Bates et al. (2010) inertial formulation on a regular Cartesian grid and does not exhibit checkerboarding.

---

## 2. Root Cause Analysis

Four causes were identified through the CAESAR-Lisflood comparison. Two were confirmed as primary drivers of the instability by the elimination of checkerboarding after fixes were applied.

### 2.1 Primary Cause: Inconsistent `q_prev`

The most significant structural difference between the two solvers concerned how unit discharge `q` is stored between timesteps. In the Bates (2010) formulation, `q_prev` appears in the numerator of eq. 9 as the inertial momentum term — its value at step `t` directly scales the flux computed at step `t + dt`.

In FloodA5, `edges.flux[e]` was being set to the unlimited Bates unit discharge `Q / width` after Phase A. However, the volume actually transferred was then separately capped in Phase B via the `DONOR_EDGE_DIVISOR` limiter. The result was a divergence between the stored momentum state and the actual hydraulic transfer: `q_prev` on the next step reflected a flux larger than what was physically moved. On a pentagonal mesh, where each of the five edges carries an independent and geometrically inconsistent `q_prev`, this overshoot drives a ping-pong between adjacent cells — the classical checkerboard pattern.

CAESAR-Lisflood avoids this entirely by applying its volume limiter directly to `q` before storage, ensuring the stored momentum is always consistent with the actual transfer.

### 2.2 Primary Cause: Missing Froude Limiter

CAESAR-Lisflood caps unit discharge at `q_max = h_flow × √(g × h_flow) × 0.8` (Froude number ≤ 0.8) after computing the Bates eq. 9 value. FloodA5 had no equivalent limiter. On a regular rectangular grid this is less critical because directional consistency damps supercritical oscillations; on a pentagonal mesh with five independently-evolving `q_prev` values per cell there is no such damping and the supercritical mode propagates freely.

### 2.3 Contributing Factor: Incorrect CFL Formula

The existing CFL timestep used a diffusive-wave stability criterion (`dt ≤ CFL × dx² / 2D`), which is derived for parabolic diffusion equations. The Bates (2010) formulation is an inertial (hyperbolic) equation and requires the wave-speed Courant criterion: `dt ≤ C × dx / √(g × h)`. These give different timestep bounds and the diffusive form was not guaranteeing stability in the inertial regime.

### 2.4 Contributing Factor: No `h_flow` Threshold

FloodA5 computed flux whenever `h_flow > 0`. CAESAR-Lisflood uses an `hflow_threshold = 0.001 m` below which edges are treated as dry and flux is set to zero. Without this threshold, near-dry edges carry stale `q_prev` momentum through negligible WSE gradients, introducing small spurious fluxes that compound the instability over many steps.

### 2.5 Incidental Finding: Integer Division Bug in CAESAR-Lisflood

During the comparison, a bug was identified in the CAESAR-Lisflood C# source. The denominator friction term uses `Math.Pow(hflow, (10 / 3))`. In C#, `10 / 3` is integer division, evaluating to `3` rather than `3.333...`. The correct exponent is `10.0 / 3.0`. At depths below 1 m this under-estimates friction, making the denominator smaller and the flux larger than the Bates formulation intends. FloodA5 already used the correct floating-point form. This bug should be corrected in the CAESAR-Lisflood codebase.

---

## 3. Changes Made (`standard_flow_fixes` branch)

All changes are confined to `FloodModel.jl`. The SGS path is unmodified in this branch.

### 3.1 Wave-Speed CFL (`_cfl_dt`)

The diffusive-wave stability formula was replaced with the wave-speed Courant criterion matching CAESAR-Lisflood:

```
# Before
dt = 0.5 × dx² / (2D)     where D ≈ (1/n) × h^(5/3) × √0.001

# After
dt = 0.7 × dx_min / √(g × h_max)
```

The Courant number of 0.7 matches the CAESAR-Lisflood default. The formula uses the minimum cell length scale (`dx = √area`) and the global maximum water depth, which is faster to compute than the previous per-cell loop and gives more appropriate timesteps across the full depth range.

### 3.2 New Constants

| Constant | Value | Source | Purpose |
|---|---|---|---|
| `HFLOW_THRESHOLD` | 0.001 m | CAESAR `hflow_threshold` | Minimum `h_flow` for a live edge |
| `FROUDE_LIMIT` | 0.8 | CAESAR `froude_limit` | Maximum subcritical Froude number |

### 3.3 Updated `_bates_flux`

The dry-edge return threshold was changed from `h_flow <= 0.0` to `h_flow <= HFLOW_THRESHOLD`. This function continues to be used by the SGS path, which receives the threshold update automatically.

### 3.4 New Function: `_bates_flux_limited`

A new inline function `_bates_flux_limited` was created for use exclusively by `step_standard!`. It extends `_bates_flux` with three linked fixes and returns `(Q, q_stored)` rather than just `Q`:

**Fix A — Froude limiter:** after computing Bates eq. 9, `q` is clamped to `±(h_flow × √(g × h_flow) × 0.8)`. Prevents supercritical discharge and suppresses the oscillation mode on the pentagonal mesh.

**Fix B — Volume limiter:** `q` is further clamped so that `|Q × dt| ≤ depth_donor × width / 5`. No more than ~20% of the donor cell's water can leave via one edge per step. Matches the CAESAR-Lisflood `depth/4` threshold / `depth/5` cap.

**Fix C — Consistent `q_stored`:** the post-limiting unit discharge is returned as `q_stored` and written to `edges.flux[e]`. This ensures the momentum state carried into the next step is consistent with what was actually transferred — resolving the primary divergence identified in §2.1.

### 3.5 Updated `step_standard!` Phase A

Phase A now calls `_bates_flux_limited` instead of `_bates_flux`. The donor depth (water depth of the higher-WSE cell) is identified with a single ternary before the call. `edges.flux[e]` is set to `q_stored` (the post-limiting value), not `Q / width`. The `DONOR_EDGE_DIVISOR` cap in Phase B is retained as a last-resort mass-conservation guard but should rarely bind now that the primary limiters are applied inside the flux kernel.

---

## 4. Validation Result

The patched standard flow solver was tested on the Carlisle domain — the scenario where checkerboarding was originally observed. The changes entirely eliminated the checkerboard instability. Mass balance remains exact.

This confirms that the `q_prev` inconsistency and the missing Froude limiter were the root causes, and that the pentagonal A5 topology is not itself a fundamental obstacle to stable inertial shallow-water modelling.

The next step is a formal benchmark comparison against CAESAR-Lisflood on the Carlisle test case. The main remaining structural difference between the two solvers is the `L_eff = L × cos θ` non-orthogonality correction, which CAESAR does not require on its orthogonal rectangular grid. Small systematic differences in flux magnitude on skewed A5 edges are expected and reflect the grid geometry rather than a bug.

---

## 5. SGS Implementation Plan

The SGS solver (`step_sgs!`) currently calls the original `_bates_flux`. It receives the `HFLOW_THRESHOLD` update automatically (since that was changed in `_bates_flux` itself), but none of Fixes A, B, or C. The SGS path has its own bespoke pre-processing — dry-cell effective WSE, `z_sill_eff` from the hypsometric edge sill, and Bug 49's `h_flow` cap — which means the fixes cannot be applied identically to the standard flow approach.

### 5.1 Fix A — Froude Limiter (apply; straightforward)

The Froude cap is purely a function of `h_flow`, `g`, and `q`. It can be applied to the SGS path after the `_bates_flux` call, inline in Phase A of `step_sgs!`. The SGS `h_flow` is already capped at `max(depth_ci, depth_cj)` (Bug 49 Fix B), so the effective `h_flow` seen by the Froude limiter will be bounded by the actual water depth — this is physically correct and no special handling is needed.

### 5.2 Fix B — Volume Limiter (reformulate for SGS)

The CAESAR volume limiter uses `depth_donor × width / (5 × dt)` as the per-edge cap. In the SGS context, `water_depth` is derived from the hypsometric inverse (depth above `z_min`) and is a poor proxy for how much water is available to flux at a given edge — the sub-cell geometry means the actual movable water may be much less than the cell-average depth suggests.

The existing `DONOR_EDGE_DIVISOR` cap in Phase B already limits the SGS transfer to `volume[donor] / 10` per edge, which is a volume-based criterion and more appropriate than a depth-based one for hypsometric cells. The recommendation is to **retain the Phase B cap as the primary volume limiter for SGS** rather than adding a depth-based per-edge cap inside the flux kernel, avoiding misrepresentation of sub-cell storage.

### 5.3 Fix C — Consistent `q_stored` (apply; highest priority)

This is the highest priority fix for the SGS path. The same divergence between stored `q_prev` and actual transfer that drove standard flow checkerboarding is present in `step_sgs!`. The residual SGS oscillations observed at Carlisle (Bugs 48–49) were partially addressed by the dry-cell WSE fix and `z_sill_eff`, but the `q_prev` inconsistency was never resolved. After Fix A is applied, the post-limiting `q` should be written back to `edges.flux[e]` rather than the raw Bates value.

### 5.4 Implementation Approach for SGS

Rather than creating a `_bates_flux_limited_sgs` variant, the SGS fixes should be applied **inline in Phase A of `step_sgs!`**, immediately after the `_bates_flux` call. This keeps the SGS-specific pre-processing (`wse_eff`, `z_sill_eff`, `h_flow` cap) visible and co-located with the stability fixes rather than hidden inside a shared kernel.

The SGS Phase A additions, after the `_bates_flux` call, would be:

```julia
Q, q_raw = _bates_flux(...)   # refactor to return q_new as well as Q

# Fix A: Froude limiter (uses h_flow_eff, already computed above for z_sill_eff)
q_max    = h_flow_eff * sqrt(_G * h_flow_eff) * FROUDE_LIMIT
q_stored = clamp(q_raw, -q_max, q_max)

# Fix C: write post-limiting q back (volume limiter stays in Phase B as-is)
edges.flux[e] = q_stored
edge_vol[e]   = q_stored * edges.width[e] * dt
```

Note that `_bates_flux` currently returns only `Q` (volumetric flux). To apply Fix C cleanly it will need to also return `q_new` (unit discharge), or the caller can derive it as `Q / width` before clamping. Returning both from `_bates_flux` is the cleaner long-term solution and brings it into alignment with `_bates_flux_limited`.

---

## 6. Fix Summary

| Fix | Standard Flow | SGS | Notes |
|---|---|---|---|
| Wave-speed CFL | ✅ Applied | ✅ Shared (`_cfl_dt`) | Both solvers use the same function |
| `HFLOW_THRESHOLD` (0.001 m) | ✅ Applied | ✅ Applied (via `_bates_flux`) | Changed in shared `_bates_flux` |
| Fix A: Froude limiter | ✅ Applied | ⏳ Pending | Apply inline in `step_sgs!` Phase A |
| Fix B: Volume limiter | ✅ Applied (depth-based, in kernel) | — Use Phase B cap | Volume-based cap more appropriate for SGS hypsometric cells |
| Fix C: Consistent `q_prev` | ✅ Applied | ⏳ Pending (highest priority) | Apply inline after Froude clamp in `step_sgs!` Phase A |

---

*FloodA5 — University of Canterbury | branch: `standard_flow_fixes` | June 2026*
