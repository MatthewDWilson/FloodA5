# FloodA5 — Hydraulics Reference

_Technical reference for the Bates (2010) inertial formulation as implemented in FloodModel.jl. Updated 2026-05-12 after flat-terrain point-source validation._

---

## 1. Governing Equation

Both `step_standard!` and `step_sgs!` use the inertial shallow-water formulation of Bates, Horritt & Fewtrell (2010), equation 9:

```
q^t = [ q^{t-dt} - g · h_flow · dt · (dWSE / L_eff) ]
      / [ 1 + g · h_flow · dt · n² · |q^{t-dt}| / h_flow^(10/3) ]

Q^t = q^t · width
```

| Symbol | Units | Description |
|--------|-------|-------------|
| `q^{t-dt}` | m2/s | Unit discharge from previous timestep, signed i to j |
| `h_flow` | m | Flow depth at edge: `max(WSE_i, WSE_j) - z_sill`, floored at 1e-6 m |
| `dWSE` | m | `WSE_i - WSE_j`; positive when cell i is higher |
| `L_eff` | m | `centre_dist x cos_theta` |
| `width` | m | Shared edge length (m) |
| `n` | s/m^(1/3) | Manning's roughness |
| `g` | m/s2 | 9.81 |
| `Q^t` | m3/s | Volumetric flux |

---

## 2. Sign Convention

When `WSE_i > WSE_j` (cell i is higher, flow direction is i to j):
- `dWSE > 0` -> numerator decreases -> `q_new < 0` -> `Q < 0`
- `dV[ci] += Q*dt` (negative = ci loses) correct
- `dV[cj] -= Q*dt` (positive = cj gains) correct

| Condition | dWSE | Q | dV[ci] | dV[cj] |
|-----------|------|---|--------|--------|
| WSE_i > WSE_j (flow i to j) | + | - | loses | gains |
| WSE_j > WSE_i (flow j to i) | - | + | gains | loses |
| WSE_i = WSE_j (equilibrium) | 0 | ~0 | none | none |

EdgeList sign: `flux > 0` means flow from `cell_j` to `cell_i`. The `_bates_flux` docstring is correct. An old inline comment in `step_standard!` that said "Q>0 means flow from cell_i to cell_j" was wrong and has been corrected (Bug 45).

---

## 3. h_flow Floor (Bug 38)

```julia
h_flow = max(WSE_i, WSE_j) - z_sill
h_flow <= 0.0 && return 0.0
h_flow  = max(h_flow, 1e-6)   # prevents h^(10/3) -> 0 causing NaN denominator
```

Without the floor, at depths ~4e-4 m, `h_flow^(10/3) ~ 1.4e-13`. When `q_prev` is non-zero, the denominator term `g*h*dt*n2*|q_prev| / h^(10/3)` becomes `Inf`, giving `q_new = NaN` and propagating to NaN volumes by step ~50. The 1 micrometre floor is hydrologically negligible.

---

## 4. Sill Elevation

**Standard flow:** `z_sill = max(elev_i, elev_j)` -- conservative; water cannot flow until WSE exceeds the receiving bed elevation.

**SGS flow:** `z_sill = min DEM elevation along shared boundary` -- pre-computed in `build_sgs_tables!`, stored as `sgs_edge_sills`.

---

## 5. Non-Orthogonality Correction

```
L_eff = L x cos_theta
```

`cos_theta` computed in `_edge_cos_theta` (local equirectangular projection, dot product of centre-to-centre unit vector with edge face-normal). Clamped to min 0.10.

Typical values at res 14: min=0.791, mean=0.913, max=0.969 (from validation).

Full Weller (2014) over-relaxed decomposition deferred to Phase 5.

---

## 6. Water Depth Pre-Sync (Bug 37)

Sources write directly to `state.volume`. Before each physics step, `water_depth` must be synced so the flux loop and CFL see the correct WSE:

```julia
if !use_sgs
    for i in eachindex(state.cell_ids)
        state.cell_area[i] >= 1.0 &&
            (state.water_depth[i] = state.volume[i] / state.cell_area[i])
    end
end
```

Without this, the first step has `wet=0` despite non-zero volume, and CFL returns the default 60s fallback.

---

## 7. Velocity -- Not Yet Implemented (Bug 36)

`state.velocity` is always zero. Recommended fix after the flux loop:

```julia
fill!(state.velocity, 0.0)
for e in 1:edges.n_edges
    ci = edges.cell_i[e]; cj = edges.cell_j[e]
    absQ = abs(edges.flux[e]) * edges.width[e]
    state.velocity[ci] += absQ; state.velocity[cj] += absQ
end
depth_thresh = 1e-4
for i in eachindex(state.cell_ids)
    h = state.water_depth[i]
    state.velocity[i] = h > depth_thresh ?
        state.velocity[i] / (state.cell_area[i] * max(h, depth_thresh)) : 0.0
end
```

---

## 8. CFL Timestep

```julia
dt <= CFL * dx^2 / (2*D)
```

- `CFL = 0.5`, `dx = sqrt(cell_area)`, `D = (1/n) * h^(5/3) * sqrt(S_ref)`, `S_ref = 0.001`
- Floor: `dt = max(dt, 0.1)`; fallback if all cells dry: `dt = 60.0`

Validation finding: at 1000 mm/hr and 0.152m depth on res-14 mesh, CFL gives dt ~690,000s. `dt_max = 30s` is always binding at realistic depths.

---

## 9. Volume Limiter

Applied once per cell after all edge fluxes accumulate:

```julia
if dV[i] < -0.5 * state.volume[i]
    dV[i] = -0.5 * state.volume[i]
end
state.volume[i] = max(0.0, state.volume[i] + dV[i])
```

Cap at 50% drainage per step. Previous per-edge application (Bug 30) allowed 250% over-drainage.

---

## 10. Mass Balance

Confirmed exact in all validation runs (2026-05-12). `domain_vol = rate * t` to <0.002% at all checkpoints across four rainfall rates.

Debug timing note: `vol_sum` in debug logs reads volume after rain injection but before `t` increments -- it reflects one step ahead of `domain_vol` in the progress log. Both are correct; they measure different moments.

---

## 11. Two-Solver Architecture

| | `step_standard!` | `step_sgs!` |
|---|---|---|
| WSE from | `elevation + volume / cell_area` | Hypsometric inverse interpolation |
| Sill | `max(elev_i, elev_j)` | Pre-computed edge minimum DEM |
| Pre-sync needed | Yes (Bug 37) | No |
| Requires SGS tables | No | Yes |

---

## 12. Validated Behaviour (2026-05-12)

On a 61-cell flat mesh (res 14, no DEM, fully connected), single `--rainpoint` at 1000 mm/hr for 1 hour:
- Mass balance: exact at all checkpoints
- Flux direction: correct (outward from source)
- Ring cascade: rings 1-4 wet in sequence, timing consistent with analytical estimates
- No NaN, no negative volumes, no uphill flow
- CFL: dt_max always binding; adaptive CFL formula not active at these depths

---

## 13. Planned Physics Upgrades (Phase 5)

- Weller (2014) full non-orthogonality decomposition
- LSQ gradient reconstruction for slope term
- Riemann solver option
- GPU solver kernel
- Vector velocity with edge-normal decomposition

---

## 14. References

- Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). Journal of Hydrology 387(1-2), 33-45.
- Weller, H. (2014). Geoscientific Model Development 7, 779-797.
- Thacker, W.C. (1981). Journal of Fluid Mechanics 107, 499-508.
