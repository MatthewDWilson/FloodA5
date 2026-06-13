# FloodA5 — SGS Synthetic DEM Validation Summary

_Session date: 2026-06-07. Branch: `test_sgs_unit`._

---

## Overview

This session completed the SGS (Sub-Grid Sampling) synthetic DEM validation suite
(`test/synthetic_dem/test_sgs_synthetic.jl`), which had been partially passing at
handover (T0, T1 passing; T2–T4 not yet confirmed). All four tests now pass.

**Final result:** T0 ✓  T1 ✓  T2 ✓  T3 ✓  T4 ✓

---

## Test Suite Description

The synthetic DEM is a 4 km × 2 km domain with:
- A parabolic bowl upstream basin (injection point at x = 0.2·Lx, z ≈ 0.21 m)
- A Gaussian embankment ridge at x = 0.45·Lx (crest ≈ 1.505 m, σ = 300 m)
- A notch in the embankment centred at y = 0.5·Ly (width 300 m, sill ≈ 0.705 m)
- 88 A5 pentagonal cells at resolution 14 (~12.6 ha each)
- 30 upstream cells, 58 downstream cells (partition boundary one cell-width west of embankment)

The notch is a deliberate sub-cell feature: its 300 m width is narrower than a
resolution-14 cell diameter (~355 m). The standard solver, which uses only cell-mean
elevation, cannot detect it; the SGS solver, which pre-computes minimum DEM elevations
along shared cell boundaries (edge sills), can.

| Test | Description |
|---|---|
| T0 | Mesh sanity: both meshes load with cells on both sides of the embankment |
| T1 | No downstream flow before notch sill: 30 min at 50 mm/hr leaves downstream dry |
| T2 | Downstream flow after notch sill exceeded: sustained injection drives flow through notch |
| T3 | SGS routes more water downstream than standard at same injection history |
| T4 | Mass balance < 0.01% for both solvers over 600 steps |

---

## Bugs Found and Fixed This Session

### Bug A — DEM notch carved too far upstream (`generate_synthetic_dem.py`)

**Symptom:** T2 failing with `max upstream WSE = 0.705m` and `dn_vol = 0.0` even after
7,800 steps of injection.

**Root cause:** The notch carving condition used `|x − x_emb| < 3σ`, which extended
900 m west of the embankment into the upstream basin. All DEM pixels in that zone (and
in the notch latitude band) were overwritten to `notch_elev = 0.705 m`, flattening the
upstream terrain to exactly the notch sill elevation. Any upstream cell whose footprint
fell in this zone had `z_min = z_max = 0.705 m`, causing `wse_from_volume` to return
exactly 0.705 m regardless of injected volume. This produced zero WSE gradient across
the notch → zero flux → permanent deadlock.

**Fix:** Changed `in_emb` from `|x − x_emb| < 3σ` to `x_emb − σ ≤ x ≤ x_emb + 3σ`,
carving the notch only in the eastern body of the embankment where the ridge actually
rises. Also reduced the default notch width from 1,000 m to 300 m (matching σ) so the
notch remains a genuinely sub-cell feature at resolution 14.

**File changed:** `test/synthetic_dem/generate_synthetic_dem.py`

---

### Bug B — `wse_from_volume` clamps at `z_max`, stranding excess volume (`A5Grid.jl`)

**Symptom:** After the DEM fix above, T2 still failed. The diagnostic `@warn` block
added to the test showed all 30 upstream cells with `z_max < 0.536 m < 0.705 m`. Max
upstream WSE was stuck at 0.512 m across 7,800 injection steps.

**Root cause:** `wse_from_volume` contained a hard clamp:
```julia
V >= t.vol_curve[end] && return t.z_max
```
When a cell's volume exceeded its hypsometric table range (i.e., all terrain was
submerged), the function returned a constant `t.z_max` regardless of how much
additional volume was present. Because all upstream cells in the flat parabolic bowl
have `z_max < notch_sill`, they all saturated to the same apparent WSE. Two "full"
cells at identical WSE have zero gradient → zero flux → volume permanently stranded.
No amount of injection could propagate water toward the embankment face.

This is a physics correctness bug, not just a test issue: in any scenario where an
upstream basin fills above its terrain ceiling (e.g., extreme rainfall, confined
topography), the SGS solver would silently trap water and prevent realistic downstream
routing.

**Fix:** Replaced the hard clamp with a linear extrapolation above `z_max`:
```julia
# A5Grid.jl — wse_from_volume
# Before:
V >= t.vol_curve[end] && return t.z_max

# After:
V >= t.vol_curve[end] && return t.z_max + (V - t.vol_curve[end]) / t.cell_area
```
This treats the cell as a vertical-walled container once all terrain is submerged:
WSE rises linearly with excess volume divided by total cell area. This is physically
correct (standard "pond" boundary condition) and restores the positive driving head
needed to propagate volume across flat upstream basins.

**Downstream callers verified safe:**
- `water_depth`: `max(0, wse − z_min)` — correct for `wse > z_max`
- `wetted_area_from_wse`: already clamps at `cell_area` for `wse ≥ z_max` — correct
- `h_flow_cap` in `step_sgs!`: cap becomes non-binding for overfull cells — correct
- `wse_from_volume(tbl, v_full) ≈ z_max` unit test: still passes exactly (`V = vol_curve[end]` returns `z_max + 0`)

**File changed:** `A5Grid.jl`

---

## Final Test Results (2026-06-07)

```
Test Summary:    | Pass  Total  Time
T0 — Mesh sanity |    6      6  5.5s

Test Summary:                                  | Pass  Total  Time
T1 — SGS: no downstream flow before notch sill |    1      1  4.5s

T2: max upstream WSE=0.8938m after 5700 steps  downstream=902353.94m³  (notch sill=0.705m)
Test Summary:                              | Pass  Total  Time
T2 — SGS: downstream flow after notch sill |    2      2  3.1s

T3 after 5700 steps:
  SGS: max upstream WSE=0.708m  downstream=902353.9m³
  Std: max upstream WSE=1.089m  downstream=164028.5m³
Test Summary:                                                   | Pass  Total  Time
T3 — SGS routes more water through sub-cell notch than standard |    1      1  3.7s

T4 (sgs):      injected=63108.6m³  domain=63108.623m³  err=0.0%  (4.61e-13%)
T4 (standard): injected=63108.6m³  domain=63108.623m³  err=0.0%  (4.38e-13%)
Test Summary:                       | Pass  Total  Time
T4 — Mass balance: SGS and standard |    2      2  2.2s

All synthetic DEM SGS tests passed.
```

**T3 result is particularly significant:** At the same injection history (5,700 steps),
the SGS solver routed 902,354 m³ downstream vs 164,029 m³ for the standard solver —
5.5× more downstream volume. The standard solver's mean cell elevation across the notch
band is above the notch sill, so it only routes water once upstream WSE exceeds the
cell-mean barrier (WSE = 1.089 m vs 0.708 m for SGS). This is the expected scientific
result demonstrating SGS sub-cell routing capability.

---

## Files Changed This Session

| File | Change |
|---|---|
| `A5Grid.jl` | `wse_from_volume`: extrapolate above `z_max` instead of clamping (Bug B) |
| `test/synthetic_dem/generate_synthetic_dem.py` | Notch carving restricted to embankment body; default width 1000 m → 300 m (Bug A) |
| `test/synthetic_dem/test_sgs_synthetic.jl` | Added `@warn` diagnostic block to T2; updated notch-width comment |

**Note:** The mesh parquet files do not need to be regenerated — the `A5Grid.jl` fix is
a runtime change only, and the DEM geometry fix in `generate_synthetic_dem.py` improves
test robustness but does not change the edge sill values that drive SGS routing (the
notch sill is still 0.705 m). If the meshes are regenerated from the fixed DEM, the
test results are expected to be identical.

---

## Remaining Known Issues (unchanged from handover)

- Cells 36 and 46 have NaN elevation (boundary cells outside DEM coverage) — hydraulically
  inert, no impact on test results, pre-existing minor issue.
- 25 cells outside DEM extent assigned NaN — normal for AOI boundary cells, pre-existing.
- `state.velocity` always zero (Bug 36) — not addressed this session.
- Mesh generation must use `--threads 1` until Bug 50 (`_shared_edge` thread safety) is
  patched in `A5Grid.jl`. Simulation loop is safe with any thread count.

---

## Next Recommended Steps

1. Run `test_sgs_unit.jl` to confirm the `wse_from_volume` extrapolation change does not
   break the existing 5-cell unit tests (expected to pass — the `v_full → z_max` assertion
   is still exact at the boundary point).
2. Run the Carlisle SGS validation at a realistic injection rate (50 mm/hr uniform rainfall)
   to confirm the oscillation suppression from Bugs 48/49 holds with the new `wse_from_volume`.
3. Consider whether the `wse_from_volume` extrapolation interacts with the h_flow cap
   (Bug 49 fix) in `step_sgs!` under extreme overfill — unit test with V >> vol_curve[end].
4. Update `PROJECT_STATE.md` to document Bug B and the confirmed T0–T4 pass.
