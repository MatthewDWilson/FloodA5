# FloodA5 — SGS R-A Flux Formulation: Implementation Plan

**Branch:** `sgs_ra_flux`  
**Supersedes:** `sgs_flow_fixes` (merged to `main`)  
**Date:** 2026-06-08  
**Status:** Stage 1 in progress

---

## 1. Motivation

The existing SGS flux kernel (`step_sgs!`) uses the Bates et al. (2010) wide-channel
approximation:

```
q = (q_prev - g·h_flow·dt·dWSE/L) / (1 + g·dt·n²·|q_prev| / h_flow^(7/3))
Q = q · width
```

This formulation treats the edge as a wide, shallow sheet flow. For sub-cell channel
geometry at high resolution (res 18, `dx ≈ 22 m`) it dramatically underestimates
friction resistance for concentrated channel flow, producing momentum-driven
checkerboard oscillation that limiters alone cannot cure.

LISFLOOD-FP's SGC solver uses the Manning hydraulic radius formulation instead:

```
Q = (Q_prev - g·A·dt·Sf) / (1 + g·dt·n²·|Q_prev| / (R^(4/3)·A))
```

where `A` (m²) is the cross-sectional flow area and `R = A/P` (m) is the hydraulic
radius. As the channel fills, `A` grows and `R` grows, making the friction term
`n²·|Q| / (R^(4/3)·A)` self-stabilising — large fluxes are naturally damped without
any external limiter. This is physically correct for channel flow at any resolution.

FloodA5 already samples the DEM along each shared edge to find the sill elevation
(`sgs_edge_sills`). The extension to the R-A formulation requires retaining the
**full edge elevation profile** — not just the minimum — and building per-edge
hypsometric curves for `A(wse)` and `P(wse)`. This is the same incremental
hypsometric algorithm already used for cell-level curves.

---

## 2. Architecture Overview

The change spans three layers:

```
A5Grid.jl                    FloodModel.jl               EdgeList
─────────────────────────    ───────────────────────     ──────────────────
build_sgs_tables! Step 2  →  SGSTable.edge_area_curves   edges.flux_Q (new)
  edge elevation profile  →  SGSTable.edge_perim_curves  (m³/s, SGS only)
  sgs_edge_area_curve     →  flow_area_from_wse()
  sgs_edge_perim_curve    →  wetted_perim_from_wse()
  (new parquet columns)   →  _manning_flux_ra()
                          →  step_sgs! Phase A (replaces _bates_flux)
```

---

## 3. Stage 1: Pre-processing (`A5Grid.jl`)

### 3.1 Extended edge profile sampling in `build_sgs_tables!`

The existing Step 2 (edge sill sampling) already samples the DEM along each shared
edge arc and records the minimum elevation. The change is to **retain the full sorted
profile** and build cross-sectional hypsometric curves at the same `n_bins` knot
points used for cell curves.

**Algorithm for each edge (ci → neighbour slot s):**

```julia
# Existing: sample edge arc, take minimum
edge_elevs = [_bilinear(...) for each arc sample point]
filter!(isfinite, edge_elevs)
sort!(edge_elevs)
edge_sills_mat[s, ci] = edge_elevs[1]   # unchanged

# New: build flow area and wetted perimeter curves from the sorted profile
W  = _haversine_m(lon1, lat1, lon2, lat2)   # edge length (m)
dx = W / length(edge_elevs)                  # width element per sample

# Use the same knot points as the cell hypsometric curve for this cell
bins = elev_bins_mat[:, ci]   # n_bins elevation knots

cum_A = 0.0;  cum_P = 0.0;  bin_ptr = 1;  n_wet = 0;  prev_wse = edge_elevs[1]

for k in 1:n_bins
    wse_k = bins[k]
    # depth increase for already-wet samples
    cum_A += n_wet * dx * (wse_k - prev_wse)
    # bring in new wet samples
    while bin_ptr <= length(edge_elevs) && edge_elevs[bin_ptr] <= wse_k
        cum_A += dx * (wse_k - edge_elevs[bin_ptr])
        cum_P += dx
        n_wet += 1
        bin_ptr += 1
    end
    edge_area_mat[k, s, ci] = cum_A   # flow area (m²)
    edge_perim_mat[k, s, ci] = cum_P  # wetted perimeter (m)
    prev_wse = wse_k
end
```

**Storage shape:** `(n_bins, 5, n_cells)` — flattened to `(n_bins×5,)` per row for
parquet (same pattern as other array columns).

### 3.2 New parquet columns

| Column | Array length per row | Description |
|--------|---------------------|-------------|
| `sgs_edge_area_curve` | n_bins × 5 | Cross-sectional flow area at each elevation knot, per neighbour slot (m²). Zero for unused slots. |
| `sgs_edge_perim_curve` | n_bins × 5 | Wetted perimeter at each elevation knot, per neighbour slot (m). Zero for unused slots. |

Both stored as flat Arrow list columns of length `n_bins × 5`, reshaped to
`(n_bins, 5)` on load (same as `sgs_edge_sills` extended to n_bins rows).

### 3.3 `SGSTable` struct additions

```julia
struct SGSTable
    elev_bins        :: Vector{Float64}     # WSE knots (n_bins)
    vol_curve        :: Vector{Float64}     # cumulative volume (m³)
    area_curve       :: Vector{Float64}     # wetted plan area (m²)
    cell_area        :: Float64
    z_min            :: Float64
    z_max            :: Float64
    # New fields for R-A flux formulation:
    edge_area_curves :: Matrix{Float64}     # (n_bins × 5) flow area per slot (m²)
    edge_perim_curves:: Matrix{Float64}     # (n_bins × 5) wetted perimeter per slot (m)
end
```

Backward compatibility: meshes without the new columns (pre-existing parquet files)
will load the old struct fields as before; `sgs_table()` checks for presence of
`sgs_edge_area_curve` and fills with zeros if absent, flagging a warning that
the R-A kernel will fall back to Bates.

### 3.4 New lookup functions

```julia
"""
    flow_area_from_wse(tbl, slot, wse) → Float64

Cross-sectional flow area (m²) at the edge in adjacency `slot` at water surface
elevation `wse`, by interpolating the pre-computed edge_area_curve.
"""
@inline function flow_area_from_wse(t::SGSTable, slot::Int, wse::Float64)::Float64
    col = view(t.edge_area_curves, :, slot)
    wse <= t.elev_bins[1]   && return 0.0
    wse >= t.elev_bins[end] && return col[end]
    lo, hi = 1, length(t.elev_bins)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        t.elev_bins[mid] <= wse ? lo = mid : hi = mid
    end
    frac = (wse - t.elev_bins[lo]) / (t.elev_bins[hi] - t.elev_bins[lo] + 1e-15)
    return col[lo] + frac * (col[hi] - col[lo])
end

"""
    wetted_perim_from_wse(tbl, slot, wse) → Float64

Wetted perimeter (m) at the edge in adjacency `slot` at WSE `wse`.
"""
@inline function wetted_perim_from_wse(t::SGSTable, slot::Int, wse::Float64)::Float64
    # identical structure to flow_area_from_wse, using edge_perim_curves
end

"""
    hydraulic_radius_from_wse(tbl, slot, wse) → Float64

Hydraulic radius R = A/P (m) at the edge in adjacency `slot` at WSE `wse`.
Returns 0.0 when the edge is dry (P = 0).
"""
@inline function hydraulic_radius_from_wse(t::SGSTable, slot::Int,
                                            wse::Float64)::Float64
    A = flow_area_from_wse(t, slot, wse)
    P = wetted_perim_from_wse(t, slot, wse)
    P < 1e-6 && return 0.0
    return A / P
end
```

### 3.5 `sgs_table()` function update

```julia
function sgs_table(mesh::A5Mesh, i::Int)::SGSTable
    # ... existing fields ...
    n_bins = size(mesh.array_vars["sgs_elev_bins"], 1)
    if haskey(mesh.array_vars, "sgs_edge_area_curve")
        # New mesh: reshape flat (n_bins×5,) stored columns to (n_bins, 5)
        edge_area  = reshape(mesh.array_vars["sgs_edge_area_curve"][:, i],  n_bins, 5)
        edge_perim = reshape(mesh.array_vars["sgs_edge_perim_curve"][:, i], n_bins, 5)
    else
        # Legacy mesh: fill with zeros, log warning once
        edge_area  = zeros(n_bins, 5)
        edge_perim = zeros(n_bins, 5)
    end
    SGSTable(...existing..., edge_area, edge_perim)
end
```

### 3.6 Parquet schema change and backward compatibility

New meshes gain two new array columns. Old meshes without these columns will:
- Load successfully (Arrow schema evolution — missing columns return `nothing`)
- Run with the Bates SGS kernel (fallback, same as current behaviour)
- Log a `@warn` at init time: "SGS edge hydraulic curves not present in mesh — R-A flux unavailable; using Bates approximation. Regenerate mesh with `build_sgs_tables!` to enable R-A flux."

This means **existing mesh parquet files must be regenerated** to use the R-A kernel.
Meshes for res 14 Carlisle and synthetic DEM will need a one-off rebuild.

---

## 4. Stage 1 Tests

### 4.1 Unit test: `test_sgs_edge_ra_tables.jl`

**T-EA1 — Rectangular cross-section (analytical):**
A synthetic edge elevation profile is a flat trench: all edge samples at elevation
`z_sill`, flanked by walls above. At WSE = `z_sill + d`:
- `A_expected = d × W` (width × depth)
- `P_expected = W + 2d` (approximate; exact for rectangular cross-section)
- `R_expected = A/P`

Test that `flow_area_from_wse` and `wetted_perim_from_wse` match analytical values
to within 1% for a trench 10m wide, 5m deep.

**T-EA2 — Zero flow area below sill:**
`flow_area_from_wse(tbl, slot, z_sill - 0.01) == 0.0`

**T-EA3 — Monotonically non-decreasing:**
Both `edge_area_curve` and `edge_perim_curve` are non-decreasing with WSE
for all valid synthetic profiles.

**T-EA4 — Hydraulic radius bounds:**
`0 ≤ R ≤ depth_max` for all WSE in range. `R` returns 0 when `P < 1e-6`.

**T-EA5 — Backward compatibility:**
`SGSTable` constructed without edge curve fields (zeros) does not crash
`flow_area_from_wse`; returns 0.0.

### 4.2 Regression: synthetic DEM T0–T4

Re-run `test/synthetic_dem/test_sgs_synthetic.jl` after adding the new fields.
All T0–T4 must still pass — the flux kernel is unchanged in Stage 1, only
pre-processing and table loading are touched.

The new parquet columns add ~2× to the edge table storage. The synthetic mesh
(88 cells, 186 edges) is small enough that this is inconsequential.
Carlisle res 14 mesh will need to be regenerated before Stage 2 testing.

---

## 5. Stage 2: Flux Kernel (`FloodModel.jl`)

### 5.1 New function `_manning_flux_ra`

```julia
"""
    _manning_flux_ra(Q_prev, wse_i, wse_j, z_sill, A, R,
                     L, cos_theta, n_mann, dt) → Float64

Manning R-A inertial flux kernel for SGS edges.  Implements the LISFLOOD-FP
SGC formulation:

    Q = (Q_prev - g·A·dt·Sf) / (1 + g·dt·n²·|Q_prev| / (R^(4/3)·A))

where:
  Q_prev  volumetric discharge at previous step (m³/s)  [NOT unit discharge]
  A       cross-sectional flow area (m²) at max(wse_i, wse_j)
  R       hydraulic radius A/P (m)
  Sf      = dWSE / L_eff  (positive when i > j → flow i→j → Q < 0)

Sign convention matches _bates_flux:
  Q < 0  →  flow from cell_i to cell_j  (i is higher)
  Q > 0  →  flow from cell_j to cell_i  (j is higher)

The R-A form is self-stabilising: as |Q| grows, the denominator grows as
n²|Q|/(R^(4/3)·A), providing physically correct friction scaling without
requiring an explicit Froude cap or volume limiter.

Returns 0.0 for dry edges (A ≤ 1e-6 m² or h_flow ≤ HFLOW_THRESHOLD).
"""
@inline function _manning_flux_ra(Q_prev    :: Float64,
                                   wse_i     :: Float64,
                                   wse_j     :: Float64,
                                   z_sill    :: Float64,
                                   A         :: Float64,
                                   R         :: Float64,
                                   L         :: Float64,
                                   cos_theta :: Float64,
                                   n_mann    :: Float64,
                                   dt        :: Float64)::Float64
    h_flow = max(wse_i, wse_j) - z_sill
    (h_flow <= HFLOW_THRESHOLD || A <= 1e-6) && return 0.0
    R = max(R, 1e-4)

    dWSE  = wse_i - wse_j
    ct    = max(cos_theta, 0.1)
    L_eff = max(L * ct, 1.0)

    numerator   = Q_prev - _G * A * dt * dWSE / L_eff
    denominator = 1.0 + _G * dt * n_mann^2 * abs(Q_prev) / (R^(4.0/3.0) * A)
    return numerator / denominator
end
```

### 5.2 `edges.flux_Q` — new volumetric flux array in `EdgeList`

```julia
struct EdgeList
    n_edges   :: Int
    cell_i    :: Vector{Int}
    cell_j    :: Vector{Int}
    width     :: Vector{Float64}
    L         :: Vector{Float64}
    cos_theta :: Vector{Float64}
    sill      :: Vector{Float64}
    flux      :: Vector{Float64}    # q (m²/s) — standard flow and SGS Bates
    flux_Q    :: Vector{Float64}    # Q (m³/s) — SGS R-A flux only; zeros otherwise
end
```

`flux_Q` is initialised to `zeros(ne)` at `_build_edge_list` time regardless of
flow model. For standard flow it is never written and stays zero. For SGS R-A it is
the primary momentum state; `flux` is left at zero and ignored.

The HDF5 output (`_write_frame!`) writes `flux_Q` as a new dataset
`/frames/{idx}/flux_Q` when non-zero (i.e. when using the R-A SGS kernel).

### 5.3 `step_sgs!` Phase A change

The adjacency slot for each edge must be looked up to retrieve the correct column
of `edge_area_curves` and `edge_perim_curves`. The slot is the position of `cj` in
`ci`'s adjacency list (and vice versa for the reverse lookup). This slot is already
used in the existing edge sill lookup (`_sgs_edge_sill`); the same mechanism applies.

```julia
# Phase A — replace existing _bates_flux block with:

wse_flow = max(wse_ci_eff, wse_cj_eff)

# Look up R-A geometry from the SGS tables
slot_i = _adjacency_slot(state, ci, cj)   # slot of cj in ci's adjacency
A_i = flow_area_from_wse(state.sgs_tables[ci], slot_i, wse_flow)
R_i = hydraulic_radius_from_wse(state.sgs_tables[ci], slot_i, wse_flow)

# Average A and R from both cell perspectives (symmetric treatment)
slot_j = _adjacency_slot(state, cj, ci)
A_j = flow_area_from_wse(state.sgs_tables[cj], slot_j, wse_flow)
R_j = hydraulic_radius_from_wse(state.sgs_tables[cj], slot_j, wse_flow)

A_edge = 0.5 * (A_i + A_j)
R_edge = 0.5 * (R_i + R_j)

Q_new = _manning_flux_ra(edges.flux_Q[e], wse_ci_eff, wse_cj_eff, z_sill_eff,
                          A_edge, R_edge, edges.L[e], edges.cos_theta[e],
                          0.5 * (state.manning_n[ci] + state.manning_n[cj]), dt)

# Fix C: write-back (Q in m³/s directly)
edges.flux_Q[e] = Q_new
edge_vol[e]     = Q_new * dt
```

Note: the `wse_eff` pre-processing (Bug 48 dry-cell fix) and `z_sill_eff` (Bug 49
h_flow cap reformulated for R-A) are still needed. The `h_flow_cap` becomes
`max(A_i, A_j) / max(edges.width[e], 1e-6)` — the effective depth for R-A is `A/W`,
not the raw WSE difference. The Froude and volume limiters (Bugs 53, 59) are removed.

### 5.4 `_adjacency_slot` helper

A small function is needed to find the slot index of `cj` in `ci`'s adjacency:

```julia
@inline function _adjacency_slot(state::FlowState, ci::Int, cj::Int)::Int
    for s in 1:5
        state.adj_matrix[s, ci] == cj && return s
    end
    return 1   # fallback; shouldn't occur for valid edges
end
```

This is O(5) — negligible per edge call.

### 5.5 Removal of limiters

The following code blocks introduced in `sgs_flow_fixes` are **removed** from
`step_sgs!` when the R-A kernel is active:
- Fix A (Froude limiter) — replaced by natural R-A self-stabilisation
- Fix B (volume limiter, Bug 59) — same reason
- The `h_flow_cap` / `z_sill_eff` reformulation is retained but now bounds `A`
  rather than `h_flow` directly

The Phase B `DONOR_EDGE_DIVISOR` cap is retained as a last-resort mass guard.
Fix C write-back is retained (now writing Q in m³/s to `flux_Q`).

---

## 6. Stage 2 Tests

### 6.1 Unit test: `test_manning_ra_flux.jl`

**T-RA1 — Zero gradient → Q = 0:**
`_manning_flux_ra(0.0, 5.0, 5.0, 4.0, 2.0, 0.5, 100.0, 1.0, 0.03, 30.0) ≈ 0.0`

**T-RA2 — Known Manning's equation (steady state, Q_prev = Q_new):**
For a rectangular channel 10m wide, 2m deep, n=0.03, Sf=0.001:
`Q_expected = (1/0.03) × A × R^(2/3) × √0.001`
At steady state `Q_prev = Q_new`, the R-A formula reduces to Manning's equation
exactly. Verify numerically by iterating until |Q_new - Q_prev| < 1e-6.

**T-RA3 — Sign convention:**
WSE_i > WSE_j → Q < 0 (flow i→j). WSE_j > WSE_i → Q > 0. ✓

**T-RA4 — Dry edge:**
`A = 0 → Q = 0.0`

**T-RA5 — Mass conservation:**
Closed 5-cell chain, single injection, 300 steps: `mb_err < 0.01%`

**T-RA6 — No oscillation at res 18 proxy:**
Small domain (10 cells), uniform rainfall, `dt = 2.5s`, 1000 steps.
`max(abs(diff(max_depth_series))) < 0.05m` (no period-2 oscillation).

### 6.2 Regression

- `test_sgs_unit.jl` — all pass
- `test/synthetic_dem/test_sgs_synthetic.jl` — T0–T4 all pass
- Carlisle res 14 50mm/hr 20h: `mb_err = 0.0`, `wet` stable, no visible oscillation

---

## 7. Stage 3: Cleanup and Documentation

- Remove `sgs_flow_fixes` Bugs 53 (Froude), 59 (Fix B volume) limiters from
  `step_sgs!` (they're inert after Stage 2, but should be explicitly removed with
  a comment explaining why)
- Update `HYDRAULICS.md` §11 two-solver architecture table
- Update `DATA_FORMATS.md` §1 array columns table
- Update `PROJECT_STATE.md` bug table and validation results
- Regenerate all SGS mesh parquet files (res 14 and res 18 Carlisle, synthetic DEM)

---

## 8. Stage 4 Tests (res 18 regression)

- Carlisle res 18 SGS 50mm/hr 1h: `dt` should stabilise rather than decline,
  `max_depth` stable, `mb_err = 0.0`
- Compare depth field to standard flow res 18 run at same rainfall — qualitative
  agreement expected with SGS showing more ponding in channels

---

## 9. Key Constants and Parameters

| Constant | Value | Notes |
|---|---|---|
| `HFLOW_THRESHOLD` | 0.001 m | Unchanged; dry-edge guard |
| `FROUDE_LIMIT` | 0.8 | Retained for standard flow; not used by R-A SGS |
| `DONOR_EDGE_DIVISOR` | 10 | Retained as last-resort Phase B mass guard |
| `n_bins` | 100 | Unchanged; edge curves use same knot points as cell curves |

---

## 10. Files Changed (summary)

| File | Stage | Changes |
|---|---|---|
| `A5Grid.jl` | 1 | `build_sgs_tables!` Step 2: build edge area/perim curves. New parquet columns. `SGSTable` struct: add `edge_area_curves`, `edge_perim_curves`. New functions: `flow_area_from_wse`, `wetted_perim_from_wse`, `hydraulic_radius_from_wse`. `sgs_table()`: populate new fields with backward-compat fallback. |
| `FloodModel.jl` | 2 | `EdgeList`: add `flux_Q`. `_build_edge_list`: initialise `flux_Q`. New `_manning_flux_ra` function. `_adjacency_slot` helper. `step_sgs!` Phase A: replace `_bates_flux` + limiters with `_manning_flux_ra`. `_write_frame!`: write `flux_Q`. |
| `DATA_FORMATS.md` | 3 | Document `sgs_edge_area_curve`, `sgs_edge_perim_curve`, `flux_Q` |
| `HYDRAULICS.md` | 3 | Update SGS section |
| `test/` | 1, 2 | New: `test_sgs_edge_ra_tables.jl`, `test_manning_ra_flux.jl` |

---

## 11. Note on Existing `sgs_flow_fixes` Limiters

The Froude limiter (Bug 53) and volume limiter (Bug 59) added in `sgs_flow_fixes`
were correct responses to the Bates formulation's instability, given the analysis
available at the time. They are superseded by the R-A formulation, which is
self-stabilising by construction. They should be removed from `step_sgs!` in
Stage 3 to keep the code clean, but are harmless if left in place during Stages 1–2.

---

*FloodA5 — University of Canterbury | Branch: sgs_ra_flux | 2026-06-08*
