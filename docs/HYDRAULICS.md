# FloodA5 — Hydraulics Reference

_Technical reference for the physics implemented in `FloodModel.jl` and
`surfacewater/flow2d.jl`. Supersedes the 2026-05-12 version of this document,
which predates the stability fixes (§6), the SGS R-A flux path (§9), and all
of the directional-bias correction work (§7–§8). This document reflects the
state of the `directional-bias-reformulation` branch at merge time, plus
targeted updates (§6, §11, new §14) for the `cartesian-formulation` branch's
mesh-backend generalisation (`N_SIDES` → mesh-derived `state.max_neighbours`,
A5 or the new Cartesian backend), including an open directional-bias finding
(§14) directly relevant to §8's own open items._

---

## 1. Governing Equation

Both `step_standard!` and `step_sgs!` use the inertial shallow-water
formulation of Bates, Horritt & Fewtrell (2010), eq. 9, as their default
(uncorrected) flux kernel:

```
q^t = [ q^{t-dt} - g · h_flow · dt · (dWSE / L_eff) ]
      / [ 1 + g · h_flow · dt · n² · |q^{t-dt}| / h_flow^(10/3) ]

Q^t = q^t · width
```

| Symbol | Units | Description |
|--------|-------|-------------|
| `q^{t-dt}` | m²/s | Unit discharge from previous timestep, signed i→j |
| `h_flow` | m | Flow depth at edge: `max(WSE_i, WSE_j) - z_sill`, floored at `HFLOW_THRESHOLD` |
| `dWSE` | m | `WSE_i - WSE_j`; positive when cell i is higher |
| `L_eff` | m | `L × cos θ` (legacy); superseded by the corrected forms in §7–§8 when enabled |
| `width` | m | Shared edge length |
| `n` | s·m⁻¹ᐟ³ | Manning's roughness — arithmetic mean of the two adjacent cells (`0.5·(n_i+n_j)`) for both solvers |
| `g` | m/s² | 9.81 |
| `Q^t` | m³/s | Volumetric flux |

The `10.0/3.0` exponent is deliberately a floating-point literal, not
integer division (`10/3` would evaluate to `3` in some languages, a
confirmed bug class in at least one reference implementation FloodA5 was
checked against).

The SGS solver (`step_sgs!`) uses the same equation with a different `z_sill`
source and a `h_flow` cap (§9) but is otherwise the same formulation, unless
the R-A flux kernel is active (§9.3).

---

## 2. Sign Convention

`EdgeList` stores each edge once with `cell_i < cell_j` (canonical
lower-index ordering). `flux > 0` means flow from `cell_j` to `cell_i`.

When `WSE_i > WSE_j` (i is higher):
```
dWSE > 0 → numerator decreases → q_new < 0 → Q < 0
dV[ci] += Q·dt   (negative — ci loses volume)   OK
dV[cj] -= Q·dt   (positive — cj gains volume)   OK
```

| Condition | dWSE | Q | dV[ci] | dV[cj] |
|---|---|---|---|---|
| WSE_i > WSE_j (flow i→j) | + | − | loses | gains |
| WSE_j > WSE_i (flow j→i) | − | + | gains | loses |
| WSE_i = WSE_j (equilibrium) | 0 | ~0 | none | none |

---

## 3. Dry-Edge Threshold and the `h_flow` Floor

```julia
h_flow = max(WSE_i, WSE_j) - z_sill
h_flow <= HFLOW_THRESHOLD && return 0.0   # HFLOW_THRESHOLD = 0.001 m
h_flow = max(h_flow, 1e-6)                # underflow guard, both solvers
```

`HFLOW_THRESHOLD` (0.001 m, matching the LISFLOOD-FP `DepthThresh`
convention) treats near-dry edges as fully dry rather than carrying stale
momentum through a negligible WSE gradient. The additional 1 µm floor
prevents `h_flow^(10/3) → 0` producing an `Inf` denominator and `NaN` flux —
hydrologically negligible, numerically necessary.

---

## 4. Sill Elevation

**Standard flow:** `z_sill = max(elev_i, elev_j)` — a conservative choice:
water cannot cross into a lower cell until the sending cell's WSE exceeds
the *higher* of the two bed elevations. Matches Bates (2010) literature
convention.

**SGS flow:** `z_sill = min(DEM elevation along the shared cell boundary)`,
pre-computed once in `build_sgs_tables!` and stored as `sgs_edge_sills`.
This is what allows the SGS solver to detect sub-cell channels (a notch or
ditch narrower than a cell) that the standard solver's cell-mean elevation
cannot see. See §9 for the full SGS formulation, including a `h_flow` cap
applied on top of this sill to suppress dry-cell spurious flux (a bug found
and fixed early in SGS development).

---

## 5. CFL Timestep

**Current formula (wave-speed Courant criterion, both solvers):**

```julia
dt = C * dx_min / sqrt(g * h_cfl)     # C = 0.7 (Courant number)
```

- `dx_min = sqrt(min(cell_area))` across the mesh
- `h_cfl` = the 99th-percentile wet-cell depth (not the absolute maximum —
  a single overfull SGS cell should not dictate the timestep for the whole
  domain)
- SGS solver additionally caps `water_depth` used for this calculation at
  the cell's terrain range (`z_max - z_min`), since `wse_from_volume` can
  legitimately extrapolate WSE far above `z_max` for an overfull cell
  (§9.2) — that extrapolated head is correct for driving flux but is not a
  physical water column depth for CFL purposes
- Floor: `dt = max(dt, 0.1)`; fallback if the whole domain is dry: `dt = dt_max`

This replaced an earlier diffusive-wave formula (`dt <= 0.5*dx^2/(2D)`),
which is the correct stability criterion for a parabolic diffusion equation
but not for Bates' inertial (hyperbolic) formulation, and was found not to
reliably bound the timestep at realistic depths.

---

## 6. Stability Limiters (Standard Flow)

Three linked fixes are applied in `_bates_flux_limited` (standard flow
only; see §9 for how SGS handles the equivalent problem):

**Fix A — Froude limiter.** After the Bates kernel is evaluated, `q` is
clamped to `±(h_flow · sqrt(g·h_flow) · FROUDE_LIMIT)`, `FROUDE_LIMIT = 0.8`.
Prevents supercritical discharge, which — on a mesh with five independently
evolving `q_prev` values per cell instead of two on a rectangular grid — has
no natural directional damping and can drive a checkerboard oscillation.

**Fix B — Volume limiter.** `q` is further clamped so that
`|Q·dt| <= depth_donor · width / vol_limit_divisor`, i.e. no more than
`1/vol_limit_divisor` of the donor cell's water can leave via any one edge
in a single step. `vol_limit_divisor` defaults to `5.0` (A5's own neighbour
count, i.e. ~20% per edge, unchanged original behaviour) but
`step_standard!` passes `Float64(state.max_neighbours)` explicitly — `4.0`
on the Cartesian backend (`mesh/CartesianGrid.jl`), so the ~20%-per-edge
reasoning reflects each mesh's own actual connectivity rather than an
A5-specific assumption baked into the default. A coarser, cell-level
`DONOR_EDGE_DIVISOR = 10` (A5) / `_donor_edge_divisor(state) =
2×max_neighbours` cap (§11) is retained on top of this as a last-resort
mass-conservation guard; it should rarely bind once Fix A/B are active.

**Fix C — Consistent `q_stored`.** The *post-limiting* unit discharge (not
the raw Bates output) is what gets written to `edges.flux[e]` for use as
`q_prev` on the next step. This was the single largest fix found during
this project's investigation of checkerboarding on the Carlisle domain: the
raw, unlimited Bates `q` was previously stored while the *volume actually
transferred* was independently capped downstream, so the stored momentum
state diverged from what had really moved — a divergence that compounds
every step and, on a 5-edge-per-cell pentagon with no directional damping,
produces a self-sustaining oscillation. Storing the post-limiting value
closes this loop.

All three fixes together eliminated checkerboarding on the Carlisle
domain (standard flow). Full root-cause investigation against a
CAESAR-Lisflood reference implementation is retained in the project's
internal development history (see the note at the end of this document).

---

## 7. Directional Bias — Problem and Corrections (opt-in, standard flow only)

A5 pentagon edges are generally not perpendicular to the centre-to-centre
vector **d**. The legacy `L_eff = L·cos θ` scaling (§1) corrects the
*magnitude* of the driving gradient for this but not its *direction* — the
gradient is still implicitly treated as acting along **d**, when the true
face-normal direction **n̂** differs from **d** by the non-orthogonality
angle θ (measured on real A5 meshes: mean ~23°, p95 ~37–38°, max ~38° —
comfortably inside the "corrected-scheme-safe" envelope of <=70° cited in
the finite-volume literature, but far from negligible). This produces a
measurable systematic bias: point-source flood fronts are elongated along a
preferred axis rather than circular, and flow on a planar slope deviates
from the analytically correct downslope direction.

**Three independent, opt-in corrections have been implemented and tested to
different degrees. None is the default.** All are controlled by CLI flags
documented in `USER_GUIDE.md` §4.3a. They are validated only for the
standard solver. **Correction (WP0, 2026-10):** earlier versions of this
document said none of them affect the SGS solver. In fact `step_sgs!` does
apply Correction 1 (WLSQ + V̂, with `--gradient-correction-alpha`) when
`--gradient-correction on` is set, and honours `--q-centre-theta`; only the
diamond method and cell-vector momentum are ignored under SGS. See §9.4.

### 7.1 Correction 1 — WLSQ gradient + skewness term (`--face-flux-method legacy`)

A 2D WSE gradient is reconstructed once per cell per timestep by weighted
least squares over that cell's (up to five) neighbours:

```
minimise sum_k w_k [(WSE_jk - WSE_i) - grad_WSE_i . (x_jk - x_i)]^2,   w_k = 1/|x_jk-x_i|^2
```

The projection matrix depends only on geometry and is pre-computed once at
init (`wlsq_weights`); the per-step cost is one matrix–vector multiply per
cell. At each edge, the corrected driving head is:

```
dWSE_n = c*(WSE_i - WSE_j) - alpha * L*(grad_WSE_f . V_hat)

  c        = d_hat . n_hat = cos(theta)      (face normal oriented so c >= 0)
  grad_WSE_f = 0.5*(grad_WSE_i + grad_WSE_j)  (face-averaged cell gradient)
  V_hat    = n_hat - c*d_hat                  (tangential remainder of the face
                                                normal; |V_hat| = sin(theta);
                                                stored per edge as skew_x/skew_y)
  alpha    = --gradient-correction-alpha      (default 1.0)
```

For an orthogonal edge, `V_hat = 0`, `c = 1`, and the formula reduces
exactly to the uncorrected `WSE_i - WSE_j` — confirmed backward-compatible.

**Real-mesh finding (multiple development sessions — see the project's
internal development history for the full evidence trail):** at `alpha=1`
(full correction) combined with per-edge momentum storage, this correction
can overshoot into a *mirrored* directional bias under sustained ponding
against a closed boundary, with a dt-sensitivity that is V-shaped (best at
an intermediate `dt`, not monotonically improving as `dt -> 0`) — evidence
of a genuine compounding interaction with the momentum/limiter machinery
(§6), not a simple discretisation-error story. `alpha=0` (orthogonal term
only, no tangential correction) is a smaller but stable, non-flipping
effect and is the better-supported setting of the two if this correction
path is used at all.

A real-mesh diagnostic (`test/diagnose_skew_bias.jl`) additionally found a
clean, statistically significant north/south asymmetry in the `V_hat`
x-component on a long, thin east-west test domain. Direct code inspection
found no logic bug in the orientation/geometry pipeline; the leading
hypothesis is that this reflects a genuine, inherent chirality of the A5
pentagon tiling rather than a software defect, though this has not been
proven and the upstream mesh-generation pipeline (`pya5` vertex ordering)
has not been fully audited — the project's internal development history
has the complete
evidence and open questions.

### 7.2 Correction 2 — Diamond face-flux reconstruction (`--face-flux-method diamond`)

Rather than sharing one cell-centred gradient across all five of a cell's
faces, this constructs a bespoke gradient for *each edge* directly from a
quadrilateral ("diamond") formed by the two adjacent cell centres and the
edge's own two shared vertices:

1. **Vertex reconstruction** (`mesh/DiamondFlux.jl`, `VertexRecord`): each
   mesh vertex's WSE is a fixed, pre-computed weighted least-squares
   combination of the surrounding cells' WSE (3 or 4 cells at ~2:1 ratio
   on real A5 meshes, confirmed by an Euler's-formula argument and direct
   measurement to within 0.2 percentage points on a 16,567-cell mesh).
   Proven exact for constant and linear fields when the stencil is
   non-degenerate (`k >= 3`); a `k = 2` fallback (distance-weighted average)
   is not exact for a general 2D field and is confined to AOI-boundary
   vertices (~1.2% of edges on the test mesh used for validation).
2. **Diamond gradient**: the discrete Green-Gauss (divergence-theorem)
   gradient over the four-vertex diamond — a general property of any simple
   polygon, exact for a linear field regardless of shape.
3. This composes with (1) to give a gradient exact for a linear field
   wherever the diamond is non-degenerate and both vertex reconstructions
   are well-conditioned — both conditions hold essentially universally on
   real A5 meshes (companion audit scripts:
   `test/audit_vertex_valence.jl`, `test/audit_diamond_gradient.jl`).
4. The resulting `dWSE_n` is algebraically **identical** to §7.1's formula
   for any exact gradient input (proven directly by substitution — every
   `c`, `L_n`, `L_t` term cancels) — meaning diamond and legacy differ only
   in the *accuracy* of the gradient each supplies, not in how the
   driving-head decomposition itself works. No new flux kernel was needed;
   diamond plugs directly into the existing `_bates_flux_limited_corrected`.
5. Falls back to the legacy formula automatically, per edge, if that
   edge's diamond record is degenerate.

Full derivation and proofs are retained in the project's internal
development history (available on request — not part of this public
release; see the note at the end of this document).

### 7.3 Correction 3 — Cell-vector momentum (`--momentum-model cell`)

Independent of §7.1/§7.2: replaces the five independent per-edge `q_prev`
scalars (one per pentagon face) with a single 2D discharge vector per cell,
`(qvec_u, qvec_v)`, reconstructed each step by weighted least squares from
that cell's five (post-limiter) face fluxes — a Perot-style reconstruction
for collocated unstructured meshes. At each face, `q_prev_eff` is the
symmetric projection of the two adjacent cells' vectors onto that face's
normal:

```
q_prev_eff = -0.5 * [(qvec_u_i + qvec_u_j)*n_ex + (qvec_v_i + qvec_v_j)*n_ey]
```

**Rationale:** Bates' momentum term (`q_prev` in eq. 9, §1) was derived for
a Cartesian cell with exactly two flux directions, where `q_prev` at a face
*is*, by construction, the same physical direction one step earlier. On an
A5 pentagon, five independently-evolving per-edge scalars have no
construction ensuring they represent one coherent 2D state — this is a
distinct problem from §7.1/§7.2's slope-direction issue, and was only
identified after the gradient corrections above had already been built and
found insufficient on their own.

**Empirically the stronger lever tested to date**, and the only one shown
to reach a stable fixed point under repeated forcing where the per-edge
scalar model does not (a controlled frozen-coefficient convergence
experiment — full detail in the project's internal development history).

**Only wired into `step_standard!`.** Selecting `--momentum-model cell`
under `--flow-model sgs` has no effect (a startup warning is printed).

---

## 8. Current Status and Recommendation

> **Superseded (2026-10):** the recommendation below predates WP1/WP2. Cell-vector
> momentum checkerboards on the Carlisle DEM; see §15 for the current evidence
> and recommended experimental configuration.

**None of §7's corrections are the default.** All three default to the
legacy, uncorrected Bates behaviour, unchanged from the model's original
implementation. This is a deliberate, repeatedly-reaffirmed decision across
several sessions of this project — each correction has real, measured
benefit, but none fully eliminates the bias, at least one combination is
known to make it worse under some conditions (§7.1 at `alpha=1`), and
validation to date is limited to synthetic planar-slope domains at a small
number of resolutions, not the real Carlisle DEM.

**Best-evidenced experimental configuration**, for anyone wanting to
compare against the default:

```
--flow-model standard --gradient-correction on --face-flux-method diamond --momentum-model cell
```

This combination gave the largest reduction in the north/south
volume-asymmetry benchmark of any tested configuration (down from a
baseline of ~0.87-0.92 to ~0.61), without the sign-flip instability seen
under §7.1 at `alpha=1`. It is still an incomplete fix, not a validated
production default.

**Open items, in priority order:**

1. **Limiter interaction is the leading suspect for the remaining bias.**
   Controlled testing found the Froude/volume limiter (§6) can
   substantially negate whatever directional correction the driving-force
   term supplies — plausibly why the provably-more-accurate diamond method
   did not clearly outperform the legacy `alpha=1` correction in full
   simulation, despite being strictly more accurate at the gradient level.
   No reformulation has been attempted yet. **Independent support for this
   suspicion, from a different grid entirely:** the `cartesian-formulation`
   branch found a real, small, *resolution-independent* cardinal-vs-diagonal
   directional bias on the new Cartesian mesh backend — a grid that is
   provably exactly orthogonal (`cos_theta≡1`, `skew≡0`), so A5's own bias
   mechanism (non-orthogonal face/cell-averaged-gradient misalignment)
   cannot be the cause there. The leading hypothesis for *that* finding is
   the per-edge momentum representation (`--momentum-model edge`'s
   axis-local Q-centred smoothing), untested as of this writing. If
   confirmed, it would be direct, independent, cross-grid evidence that
   momentum/limiter representation — not grid geometry — is the dominant
   remaining source of directional bias in this solver family generally.
   See §14.
2. **Boundary/ghost-edge treatment is architecturally inconsistent.**
   Interior edges use whichever correction is selected; ghost-cell boundary
   outflow (open, `ZeroGradient` boundaries) always uses the original
   uncorrected kernel regardless of `--face-flux-method`/`--momentum-model`.
   All validation to date has used `--closed-boundaries`, under which this
   is a non-issue (ghost flux is identically zero) — it has **not** been
   tested under the default open-boundary condition. This is exactly what
   a Carlisle run under this configuration will exercise for the first time.
3. Terrain-ponding (embankment) behaviour has not been separately
   characterised — only excluded from the clean benchmark via a
   terrain-feature-free test mesh.
4. No Carlisle (real DEM) validation, no analytical benchmark (Thacker,
   circular dam-break), no second-resolution convergence check.
5. No formal stability proof for the cell-vector momentum model composed
   with the limiters — the evidence is empirical, not analytical.
6. `k=2` vertex-reconstruction fallback (§7.2) is not linear-exact;
   confined to AOI-boundary edges, low priority unless a domain's active
   flow reaches the boundary.

---

## 9. SGS (Sub-Grid Sampling) Solver

### 9.1 Hypsometric tables

The SGS approach replaces a single mean bed elevation per cell with a
hypsometric curve — a pre-computed lookup table relating WSE to stored
volume and wetted plan area, built from Halton-sampled LiDAR DEM points
within each cell polygon (`build_sgs_tables!`). `wse_from_volume` performs
the runtime inverse lookup.

**Overfull-cell extrapolation:** when a cell's volume exceeds its
hypsometric table range (fully submerged), WSE is extrapolated linearly
above `z_max` (`WSE = z_max + (V - vol_curve[end]) / cell_area`) rather
than clamped. Clamping was an earlier bug that silently stranded volume in
any scenario where an upstream basin filled above its terrain ceiling — a
real physics-correctness issue, not just a test artefact.

### 9.2 SGS-specific stability handling

The SGS solver has its own dry-cell and h_flow-cap handling, independent of
§6's standard-flow limiters:

- **Dry-cell effective WSE:** a dry cell's raw `wse_from_volume(V=0)`
  returns `z_min`, which can be well above a low channel sill and drive
  spurious flux. The effective WSE used in the flux kernel is
  `max(z_sill, tbl.z_min)`.
- **`h_flow` cap:** capped at `max(depth_ci, depth_cj)` rather than the raw
  WSE difference, since an edge sill far below both cells' mean elevation
  (a genuine sub-cell channel) would otherwise give an `h_flow` far larger
  than the actual water column.
- **Froude and volume limiters** equivalent to §6's Fix A/B are applied
  inline in `step_sgs!` Phase A, using the capped `h_flow` above; the
  `DONOR_EDGE_DIVISOR` cap (§11) remains the primary volume-based
  last-resort guard, since cell-mean `water_depth` is a poor proxy for
  actually-available sub-cell volume.

### 9.3 R-A (hydraulic radius) flux kernel, `flux_Q`

An alternative SGS flux formulation using cross-sectional flow area `A` and
hydraulic radius `R = A/P` (Manning form, `Q ~ R^(4/3)*A` in the
denominator) rather than the wide-channel `h^(10/3)` approximation. Requires
extended per-edge hypsometric curves (`sgs_edge_area_curve`,
`sgs_edge_perim_curve` — new mesh parquet columns; **meshes generated before
this feature must be regenerated to use it**, older meshes fall back to the
Bates SGS kernel automatically with a warning). Self-stabilising by
construction (friction grows with `|Q|` without an explicit Froude cap), so
Fix A/B above are inert once this kernel is active. Selected automatically
whenever the extended parquet columns are present in the loaded mesh.

### 9.4 Interaction with §7/§8

**Only Correction 1 (WLSQ + V̂) reaches the SGS solver, and it is
unvalidated there**; the diamond method and cell-vector momentum do not
apply to SGS (corrected WP0, 2026-10 — this section previously said none
apply). This is a scope boundary, not an oversight: the
directional-bias-reformulation work was deliberately staged to validate on
the standard solver first. Extending it to `step_sgs!` — including
resolving how the corrected gradient interacts with the dry-cell WSE clamp
in §9.2, flagged as an open risk in prior sessions but not yet tested
against the synthetic-DEM regression suite — remains future work.

---

## 10. Velocity Computation

`state.velocity` (scalar magnitude), `state.vel_u`, `state.vel_v`
(directional components) are computed after the flux/volume update each
step from a flux-weighted average across each cell's edges, and written to
HDF5 output (`/frames/{idx}/velocity`, `vel_u`, `vel_v`). Per-cell velocity
direction was found to be a noisy diagnostic in practice — near the
wet/dry front, the locally dominant direction is the front-normal, not the
regional flow direction, which swamps the regional signal in an
instantaneous per-cell average. Integrated quantities (e.g. north/south
volume asymmetry, §7) are the more reliable diagnostic for directional
bias and are what the corrections in §7 are validated against.

---

## 11. Mass Balance and the Volume Limiter

Per-edge donor cap, applied during flux accumulation (not as a separate
post-hoc pass — an earlier post-hoc version created mass by clipping only
the donor's loss without reducing the recipient's gain):

```julia
edge_vol_capped = min(edge_vol, state.volume[donor] / _donor_edge_divisor(state))
# _donor_edge_divisor(state) = 2 * state.max_neighbours
# A5: 2*5=10 -> <=50% total drain across all 5 edges
# Cartesian (mesh/CartesianGrid.jl): 2*4=8 -> <=50% total drain across all 4 edges
# (DONOR_EDGE_DIVISOR = 10 still exists as the old A5-only constant, kept
#  for existing test fixtures — production code uses the line above.)
```

The same clipped value is applied to both donor and recipient, so mass is
conserved exactly regardless of whether this cap binds. With open
(`ZeroGradient`) boundaries active (the default BC), mass balance is
checked as `input_vol - domain_vol - vol_removed ~= 0`, where `vol_removed`
accumulates ghost-edge outflow (see `boundaryinputs/boundary_conditions.jl`).

---

## 12. Two-Solver Architecture Summary

| | `step_standard!` | `step_sgs!` |
|---|---|---|
| WSE from | `elevation + volume/cell_area` | Hypsometric inverse interpolation |
| Sill | `max(elev_i, elev_j)` | Pre-computed edge-minimum DEM (`sgs_edge_sills`) |
| Water-depth pre-sync needed | Yes | No |
| Requires SGS tables | No | Yes (`build_sgs_tables!` at mesh-build time) |
| Stability limiters | Fix A/B/C (§6) | SGS-specific dry-cell + h_flow-cap equivalents (§9.2), or self-stabilising R-A (§9.3) |
| `--gradient-correction` | **Yes** | Legacy WLSQ + V̂ only, **unvalidated** (§9.4) |
| `--face-flux-method`, `--momentum-model` | **Yes** | **No effect** (§7, §9.4) |
| Ghost-edge (open boundary) flux | Uncorrected kernel regardless of flags (§8 item 2) | Uncorrected kernel |

---

## 13. References

- Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). A simple inertial
  formulation of the shallow water equations for efficient two-dimensional
  flood inundation modelling. *Journal of Hydrology* 387(1-2), 33-45.
- Jasak, H. (1996). *Error analysis and estimation for the finite volume
  method with applications to fluid flows.* PhD thesis, Imperial College
  London.
- Moukalled, F., Mangani, L., Darwish, M. (2016). *The Finite Volume Method
  in Computational Fluid Dynamics.* Springer.
- Perot, B. (2000). Conservation properties of unstructured staggered mesh
  schemes. *Journal of Computational Physics* 159(1), 58-89.
- Weller, H. (2014). Non-orthogonal version of the arbitrary polygonal
  C-grid and a new diamond grid. *Geoscientific Model Development* 7,
  779-797.
- Thacker, W.C. (1981). Some exact solutions to the nonlinear shallow-water
  wave equations. *Journal of Fluid Mechanics* 107, 499-508.
- OpenFOAM User Guide, §4.5 Numerical Schemes (non-orthogonality limits).

---

## 14. Mesh-Backend Generalisation and the Cartesian Directional-Bias Finding

Added for the `cartesian-formulation` branch. Two independent pieces of
work, both touching this document's earlier sections:

**Generalisation.** Everything in §6 and §11 that used to read a hardcoded
`N_SIDES = 5` (A5's pentagon valence) now reads a mesh-derived
`state.max_neighbours` instead — 5 on A5, 4 on the new Cartesian backend
(`mesh/CartesianGrid.jl`, a plain rectangular grid, pure Julia, no Python
bridge), computed once from the loaded mesh's own adjacency rather than
assumed. `N_SIDES`/`DONOR_EDGE_DIVISOR` still exist as module constants in
`FloodModel.jl`, but only as legacy defaults for existing A5 test fixtures
— production code no longer reads either directly. Everything else in this
document (the governing equation, sign convention, CFL formula, dry-edge
threshold, sill elevation, SGS solver, velocity computation) is unchanged
and applies identically regardless of mesh backend; the face-flux/geometry
machinery (`_edge_geometry`, WLSQ, the diamond method, §7-§8) is likewise
backend-agnostic by construction — none of it assumes a specific neighbour
count. (The SGS solver is the exception: `build_sgs_tables!` still hardcodes
five neighbour slots, so SGS is A5-only for now.)

**The finding.** The Cartesian backend's own geometry was validated first
and is exact: `cos_theta≡1.0`, `skew≡0.0` on every edge, confirmed at two
mesh sizes with zero floating-point deviation. This positively rules out
A5's own bias mechanism (a cell-averaged gradient shared across several
non-orthogonal faces) as a possible cause of anything measured on this
grid. Despite that, a flat-terrain, centred point-injection test found a
small (~3% in front-radius terms, ~20% in a depth-at-fixed-radius metric)
but real cardinal-vs-diagonal directional preference — with the *diagonal*
directions reaching slightly further, the opposite of the classic
hop-counting prediction for a four-connected grid — and, critically, a
refinement test (dx=25m → dx=12.5m) found the deviation from isotropy
**increasing**, not decreasing (cardinal/diagonal front-radius ratio 0.9728
→ 0.9662), which rules out ordinary discretisation truncation error as the
explanation (that would be expected to shrink by roughly 4× for a halved
grid spacing). This is the same diagnostic signature used on the A5 side to
distinguish genuine structural bias from truncation error (A5's own bias
metric also grew, not shrank, with refinement, in the project's internal
Phase A diagnostics).

The leading hypothesis is the per-edge momentum representation
(`--momentum-model edge`, default): its Q-centred smoothing operates
axis-locally, so a cardinal-direction cell's inertia is fed along one
smoothed axis while a diagonal-direction cell's is fed by two independently
-smoothed axes converging — a structural, grid-spacing-independent
difference, unlike a geometric/truncation effect. The specific test that
would confirm or refute this — re-running the same flat-domain experiment
with `--momentum-model cell` — had not been run as of this document's last
update. Diagnostic tooling: `test/test_octagonal_anisotropy.jl`,
`test/audit_cartesian_geometry.jl`, `test/check_boundary_proximity.jl`; see
`TEST_CASES.md` §5. The full numbers and reasoning are in
the project's internal development history (available on request).

---

**Full detail and complete experimental history:** this project's internal
development history retains session-by-session records of the full
derivations, proofs, and negative results this document summarises —
including the dead ends and reversed decisions along the way. These
records are not part of this public release, but are retained internally
and available on request. They are kept rather than discarded precisely
because the negative results and dead ends are as load-bearing for future
development as the things that worked.

---

## 15. Friction, Momentum and Limiter Consistency (WP1/WP2, 2026-10)

From `FloodA5_Review_and_ForwardPlan_2026-10.md` §3–§4 and §9 (WP1, WP2). All
options are opt-in; the defaults are bit-identical to before.

### 15.1 Vector friction (`--friction-model vector`)

The Bates/LISFLOOD-FP friction term uses only the face-normal discharge |q_n|.
In the friction-dominated limit the face conductance then depends on face
orientation: q_f ∝ g_n/√|g_n| instead of g_n/√|∇η|. On a Cartesian grid this
deflects the flow by up to 12.8° and makes diagonal flow 2^¼ too conductive; it
is a property of the continuum limit, so refinement cannot remove it (Cozzolino
et al. 2021). Vector friction uses |q| = √(q_n² + q_t²), with q_t the
face-tangential part of the averaged adjacent cell vectors (reconstructed every
step by the Phase F WLSQ operator); the Froude limiter scales q_n by
min(1, q_max/|q|). Implemented in all Bates kernels and the SGS R–A kernels.

**Frozen-coefficient test through `step_standard!`** (uniform depth on a planar
WSE, volumes reset each step; `test/frozen_isotropy_lib.jl`). Cartesian: component
friction gives exactly the theoretical deflection atan(√tan φ) − φ (max 12.78°)
and conductance ratio 2^¼ = 1.1892; vector friction gives 0.000° and 1.0000
(`test_wp1_vector_friction.jl` V5, in CI). A5 planar meshes
(`test/wp1_frozen_isotropy_mesh.jl`), max mean deflection / conductance ratio:

| Gradient | Friction | Momentum | res 16 | res 18 |
|---|---|---|---|---|
| legacy two-point | component | edge θ=0.9 | 15.05° / 1.706 | 20.38° / 2.002 |
| legacy two-point | vector | edge θ=0.9 | 15.52° / 1.719 | 20.56° / 2.001 |
| diamond | component | edge θ=0.9 | 7.09° / 1.024 | 2.11° / 1.185 |
| diamond | vector | edge θ=0.9 | 6.65° / 1.025 | 2.27° / 1.181 |
| diamond | component | edge θ=1 | 6.58° / 1.111 | 6.58° / 1.111 |
| **diamond** | **vector** | **edge θ=1** | **0.00° / 1.000** | **0.00° / 1.000** |
| diamond | vector | cell | 0.17° / 1.001 | 0.04° / 1.000 |

Both a consistent face gradient (diamond) and vector friction are needed, and
Q-centring at θ = 0.9 — which averages a face with its "most collinear" faces,
~144° apart on pentagons — by itself deflects the steady flow.

**Full simulations.** Cartesian flat point spread (10 m³/s, closed, t = 3600 s;
`test/run_wp1_matrix.sh`, `test/wp1_metrics.jl spread`), four-fold anisotropy A4:

| dx | component, edge | vector, edge | vector, cell |
|---|---|---|---|
| 25 m | −0.063 | +0.047 | +0.092 |
| 12.5 m | −0.089 | +0.024 | +0.033 |
| 6.25 m | −0.101 | +0.012 | +0.015 |

Component friction converges to a non-zero bias; vector friction to isotropy at
first order. A5 planar slope (0.1 m³/s, closed, dt_max 2 s, t = 72 000 s;
`wp1_metrics.jl planar`), N/S volume asymmetry and plume deflection:

| Configuration | res 16 | res 18 |
|---|---|---|
| default (legacy, component, edge) | +0.68 / +15.3° | +0.866 / +14.0° |
| legacy gradient, vector | +0.65 / +14.8° | +0.874 / +14.4° |
| diamond, component, edge | −0.43 / −5.1° | −0.573 / −4.8° |
| diamond, component, cell | −0.09 / −1.1° | −0.277 / −2.3° |
| diamond, vector, edge | +0.154 / +1.9° | +0.078 / +1.3° |
| diamond, vector, cell | +0.371 / +6.2° | +0.311 / +4.0° |
| **diamond-wetdry, vector, edge, hybrid, local CFL** | +0.166 / +2.5° | **+0.050 / +0.5°** |

With vector friction the residual decreases with refinement; with component
friction it grows. Vector friction does nothing with the legacy gradient.

### 15.2 Momentum, limiters and CFL on real terrain (WP2)

Carlisle res-18, standard flow, three QVAR inflows + 10 mm/h rain, 12 h
(`test/run_jobs.sh`, `test/wp2_checkerboard.jl`). Over fully wet interior
cells: fraction that are strict local WSE extrema with |r| > 2 mm (checkerboard
indicator) and rms velocity noise, end of run:

| Change from the default | extrema | dv_rms (m/s) | |
|---|---|---|---|
| none (default) | 0.6% | 0.096 | smooth |
| `--friction-model vector` | 0.3% | 0.092 | smooth |
| `--cfl-mode local` | 0.5% | 0.096 | smooth; no early transient |
| `--limiter hybrid` | 0.3% | 0.089 | smooth |
| `--face-flux-method diamond` | 2.3% | 0.159 | bounded, not growing |
| `--limiter positivity` | 10.9% | 0.525 (max 5.4) | **unstable** |
| `--momentum-model blend` (θ=0.9) | 9.3% | 0.241 | **checkerboard** |
| `--momentum-model cell` | 28.6% | 0.211 | **checkerboard** |
| diamond-wetdry + vector + hybrid + local | 2.7% | 0.159 | bounded |

- The per-edge volume limiter does real stabilising work: positivity scaling
  alone only prevents negative volumes. `hybrid` keeps it and removes only the
  face-count-dependent donor cap.
- Cell-vector momentum (and blend, which mixes it in) checkerboards on real
  terrain; it had been tested only on the planar slope.
- The wet/dry diamond fallback did not reduce diamond's extra noise on
  Carlisle (dry-cell vertex contamination is not its main source), but
  improves the planar result at res-18.

**Carlisle res-20** (478k cells; same forcing, 3 h): the default configuration
begins to checkerboard at this resolution, while the recommended one stays
clean, at about 2.4× the steps (local CFL):

| res-20, t = 3 h | extrema | r_rms | dv_rms (m/s) | steps |
|---|---|---|---|---|
| default | 6.1% | 252 mm | 0.144 | 9 313 |
| diamond-wetdry + vector + hybrid + local | 1.6% | 11 mm | 0.073 | 22 561 |

**Recommended experimental configuration:** `--gradient-correction on
--face-flux-method diamond-wetdry --friction-model vector --limiter hybrid
--cfl-mode local` (edge momentum, θ = 0.9). It is stable on Carlisle res-18
(somewhat noisier than the default there) and markedly smoother than the default
at res-20. The review's θ-blend and positivity replacements are not usable as
proposed.

### 15.3 Defaults (2026-10-21) and the checks behind them

Since WP2b the CLI default is `--scheme recommended` (`scheme_options` in
`FloodModel.jl`; `--scheme legacy` restores the earlier behaviour):

| Flow model | Recommended options |
|---|---|
| standard | `--gradient-correction on --face-flux-method diamond-wetdry --friction-model vector --limiter hybrid --cfl-mode local --q-centre-mode consistent` |
| SGS | `--friction-model vector --cfl-mode local --q-centre-theta 1` (legacy limiter, no gradient correction) |

(Q-centring entries added by WP2c, §15.4; the checks below predate them.)

**SGS tables were wrong before 2026-10-21.** `build_sgs_tables!` stored edge
sills and R–A curves in `grid_disk` order while the solver reads them in
neighbour-list order; the two differ for 90–96% of A5 cells, so most edges used
another neighbour's sill and curves (phantom dams — on the synthetic notch DEM
the pond rose ~1 m above the embankment crest). Fixed in WP3 part 2; existing
meshes are refused until rebuilt or repaired (`mesh/repair_sgs_slot_order.py`).
SGS results obtained before then should be treated as affected.

**SGS choice** (Carlisle res-18, repaired tables, inflows + 10 mm/h, 12 h;
local WSE extrema / rms WSE residual at the end): default 1.4% / 188 mm;
vector 1.1% / 189 mm; **vector + local CFL 1.2% / 147 mm**; + hybrid 1.2% /
130 mm; + gradient correction 4.6% (noisier — the legacy WLSQ correction is
unvalidated under SGS, §9.4); hybrid alone spiked (velocity 4 m/s) — under SGS
it is positivity-only scaling, as the R–A kernel has no per-edge volume
limiter. Planar slope under 20 mm/h rain (open boundaries, 6 h): plume
deflection −2.3° (default) vs −0.2° (vector + local), extrema 0.11% vs 0.09%.

**Open boundaries, standard flow.** Carlisle (inflows + rain; free and
zero-gradient boundaries) is stable at res-18 and res-20 (§15.2). Planar slope
under 20 mm/h rain (0.4–0.5 Mm³ leaves through the boundaries): mass balance
exact; N/S asymmetry +0.032 / +3.5° (default) vs −0.004 / −0.5° (recommended);
local WSE extrema 4.4% vs 2.4%. Ghost (boundary) edges still use the component
kernel; the noise near the boundary is the same for both configurations.

**Known trade-off.** On the rained planar slope the recommended set has more
interior velocity noise (0.056 vs 0.015 m/s rms cell-to-neighbour difference).
Single-change runs attribute all of it to the diamond face gradient (diamond
alone 0.054; vector 0.014, hybrid 0.015, local CFL 0.015), consistent with
Carlisle (§15.2). The water surface is nevertheless smoother, and the diamond
gradient is what removes the A5 directional bias, so it stays in the default;
reducing its flux noise is an open item, alongside a stable mesh-general
replacement for Q-centring.

### 15.4 Q-centring on A5 (WP2c, 2026-10-22)

**The defect.** `_q_centred` adds the stored discharges of each edge's two
"most collinear" faces. Stored discharges are signed by edge orientation, and
`_build_edge_list` orients edges from the lower to the higher cell index —
geometric on a row-major Cartesian grid, arbitrary on A5. 67% (Carlisle
res-18) and 68% (planar res-18) of partners are reversed (Cartesian 0%), so on
A5 the partner term mostly *opposes* the edge's own flux. In the steady state
q(1−θ+F) = (1−θ)·q̄_partner + G (F = g·dt·n²|q|/h^{7/3}, ≈ 0.02 for 0.5 m sheet
flow at dt = 5 s) the partner term dominates whenever F ≪ 1−θ, which is the
usual case.

**Frozen-coefficient test** (planar res-16, h = 0.5 m, S = 10⁻⁴, diamond +
vector friction, 36 bearings; `test/frozen_isotropy_lib.jl`; friction balance
10.50):

| Q-centring | dt | conveyance | max mean deflection | cell rms |
|---|---|---|---|---|
| none (θ = 1) | 1, 5 s | 10.50 | 0.0° | 0.0° |
| legacy θ = 0.9 | 1 s | **1.67** | 12.7° | 75° |
| legacy θ = 0.9 | 5 s | **2.72** | 6.7° | 36° |
| oriented θ = 0.9 | 5 s | – | 57° | 30° |
| blend θ = 0.9 | 1, 5 s | 10.50 | 0.1° | 0.2° |
| consistent θ = 0.9 | 1 s / 5 s | 10.49 / 10.50 | 0.3° / 0.15° | 0.25° / 0.18° |

Oriented partners are worse: on a pentagon the partner face is ~36° off the
edge, and once the partner term dominates, the edge copies a rotated flux.
**Consistent Q-centring** takes from each oriented partner only the component
along its own normal m and the remainder from the cell vectors,
q_est = s·q_p·(m·n_e) − Ū·(n_e − (m·n_e)·m) — exact for uniform flow, and
exactly LISFLOOD-FP's Q-centring when m = n_e (Cartesian).

**Carlisle res-18**, standard flow, recommended scheme unless stated, inflows +
10 mm/h, 12 h (`test/run_jobs.sh`; `test/wp2_checkerboard.jl`): water stored in
the domain at 12 h (Mm³; inflow + rain ≈ 14.6) at Courant 0.7 / 0.35, and
local WSE extrema at Courant 0.7. Reference: raster-aligned 20 m Cartesian
grid on the same DEM (`--meshgen-cartesian`, LiDAR averaged to 20 m), where
Q-centring is correctly oriented.

| run | stored, C 0.7 / 0.35 | extrema | r_rms | dv_rms |
|---|---|---|---|---|
| Cartesian 20 m, legacy scheme θ 0.9 | 5.7 / 8.9 | – | – | – |
| Cartesian 20 m, recommended θ 0.9 | 6.1 / 7.5 | – | – | – |
| Cartesian 20 m, legacy scheme θ 1 | 4.9 / 3.5 | – | – | – |
| A5, legacy scheme, legacy QC | **11.3 / 12.8** | 0.6% | 43 mm | 0.10 |
| A5, legacy QC | **12.2 / 13.0** | 2.7% | 61 mm | 0.16 |
| A5, consistent QC | 8.9 / **7.7** | 6.6% | 41 mm | 0.22 |
| A5, θ = 1 | 6.9 / 5.4 | 15.5% | 85 mm | 0.39 |
| A5, blend θ 0.9 / 0.7 / 0.5 | – | 7.3 / 17.1 / 22.3% | – | – |
| A5, oriented QC | – | 24.2% | 54 mm | 0.38 |

A5 with legacy Q-centring holds about twice the water of the Cartesian
reference, and more as dt is reduced (its damping grows like 1/dt);
consistent Q-centring is within 2% of the Cartesian recommended run at
Courant 0.35. Even on Cartesian, storage changes by up to ~50% between Courant
0.7 and 0.35 (in opposite directions for θ = 0.9 and θ = 1): Carlisle results
are dt-sensitive in general, which needs its own study. On the rained planar
slope (§15.3) interior mean depth at 6 h is 0.0722 m (legacy QC), 0.0648 m
(consistent), 0.0601 m (θ = 1) against a dt-converged 0.059 m (θ = 1 at
Courant 0.5 and 0.35); consistent QC gives the best symmetry of all runs
(asym −0.0003) and 1.5% local extrema (legacy QC 2.4%, θ = 1 5.6%).

**Why θ = 1 checkerboards.** Linear analysis of the forward–backward update:
for the symmetric two-point operator all modes are neutral (|z| = 1) up to
κ = g·h·dt²·λ = 4, i.e. the critical Courant number is 0.707 on Cartesian and
0.85–0.88 on A5 (power iteration on the actual operator). The diamond
operator is not symmetric: on planar res-16, 190 of 1109 eigenvalues are
complex and the frictionless growth is 1.0063 per step at Courant 0.7 (two-
point 1.0000); mixing a fraction ω of the two-point jump reduces but does not
remove it (1.0018 at ω = 0.25). Bed friction outweighs this in shallow water
but not in the deep channel, where the checkerboard sits (17–24% of cells
deeper than 1 m). Even pure two-point fluxes with θ = 1 checkerboard on
Carlisle (8.6%): the neutral modes, once excited, are barely damped in deep
water. Some damping is therefore needed — the original purpose of
Q-centring — and on A5 it must not cost conveyance.

**Tried and rejected** (prototypes, not in the code): a two-level time filter
q_eff = (1−α)qⁿ + αqⁿ⁻¹ (damps all modes but lowers the stability limit to
κ ≤ 4(1−α); Carlisle 6–7% extrema); damping towards a running mean
(Carlisle 1.7% extrema, but planar depth +14–20% and doubled velocity noise);
consistent QC + running-mean damping (5.1%); smaller Courant numbers alone
(θ = 1: 10% at 0.35).

**SGS** (Carlisle res-18, SGS tables repaired, recommended SGS options, 12 h):
legacy QC 1.2% extrema, r_rms 147 mm, stored 11.3 Mm³; θ = 1 0.4%, 131 mm,
8.0 Mm³ (velocity noise 0.118 vs 0.057 m/s); consistent QC 5.5%, 131 mm,
10.3 Mm³; oriented 3.1%, 120 mm, 11.4 Mm³. The R–A kernel uses the
symmetric two-point gradient, and there θ = 1 is clean.

**Decision (2026-10-22).** `--scheme recommended` uses consistent
Q-centring (θ = 0.9) for standard flow and θ = 1 for SGS. Consistent
Q-centring fixes the conveyance defect and is the most accurate option on the
planar tests and against the Cartesian reference, at the cost of more local
WSE extrema on Carlisle (6.6% vs 2.7%) with a smoother surface overall (r_rms
41 vs 61 mm). Open: a damping of the deep-water grid modes that costs no
conveyance (or a symmetric diamond form), and the dt-sensitivity of the
Carlisle results seen on both meshes.

### 15.5 Sensitivity to round-off (WP4 stage 2, 2026-10-23)

Running the same Carlisle res-18 case (inflows + 10 mm/h, 2 h, default
options) with the cells in two different memory orders changes only the order
of floating-point summations. The two runs then differ by:

| | rms depth diff. (wet cells) | max depth diff. | cells > 1 cm | max speed diff. | total storage |
|---|---|---|---|---|---|
| SGS | 2×10⁻⁵ m | 2×10⁻³ m | 0 | 2×10⁻³ m/s | identical to 10⁻⁷ |
| standard | 1.3×10⁻² m | 0.29 m | 348 of 4,545 | 6.5 m/s | 0.02% |

SGS is robust. The standard solver amplifies round-off to centimetres in the
deep channel: the grid-scale water-surface noise there (§15.4) is not a fixed
pattern but a weakly unstable, chaotic one, bounded only by the limiters.
Integral quantities are unaffected, but point values in those cells (and in
any comparison between runs, branches or machines) are not reproducible below
the noise level. This is the same weakness as the deep-channel extrema and the
Courant sensitivity, and the strongest argument yet for a conveyance-neutral
damping of those modes (§15.4, open). Until then, compare standard-flow runs
through integral or smoothed quantities.


### 15.6 Damping diagnostics (2026-10-23)

Tools: `test/diag_step_spectrum.jl`, `test/diag_operator_spectrum.jl`,
`test/diag_divergence.jl`, `test/diag_d2_summary.sh`.

**Linear stability of the full step** (finite-difference Jacobian; lake at
rest, 3 m deep, Courant 0.7): on A5 the explicit diamond grows by 0.14% per
step without Q-centring, two-point is neutral and consistent Q-centring damps
every mode (rough gravity modes by ~3% per step). On a square grid,
Q-centring θ = 0.9 is **unstable** at Courant 0.7 (|μ| = 1.6; stable at 0.6):
on the checkerboard mode it acts as a damping factor 2θ − 1, which lowers the
forward–backward limit to κ ≤ 4 − 2β < 3.92.

**Twin runs** (Carlisle res-18, inflows + 10 mm/h, rainfall perturbed by
10⁻¹⁰, 6 h): rms depth difference 2.3 cm (A5 recommended), 4.6 cm (θ = 1),
5.9 cm (legacy), 11 cm (Cartesian legacy), 3 mm (SGS), but 2.5×10⁻⁶ m for
Cartesian recommended. Single-change runs identify discrete switches as the
main source: Cartesian recommended with the legacy CFL diverges by 45 cm (its
99th-percentile depth is discontinuous in the state); A5 recommended with the
two-point gradient by 6×10⁻⁷ m, with plain diamond (no wet/dry switch) by
1.5 cm through smooth exponential growth (e-fold ~25 min, the diamond
operator's asymmetry).

**Courant sensitivity** of storage between Courant 0.7 and 0.35 is
systematic in every configuration: +10% (A5 recommended), −20% (A5 θ = 1),
+4% (A5 legacy), +16% (Cartesian recommended), +20% (Cartesian legacy), −17%
(SGS); rms depth differences 8–38 cm. Its sign follows the momentum
treatment, consistent with Δt-dependent terms: Q-centring's damping rate
(1 − θ)/Δt, and the per-step volume limiter and donor cap.

**Semi-implicit free surface** (θ-method on the whole face gradient, research
branch `d3-semi-implicit`, not merged): linear analysis predicts strong
damping of grid-scale gravity modes (|μ| ≈ 0.90 per step), but on Carlisle it
changes neither the extrema, nor the twin divergence, nor the Courant trend
(storage identical to three decimals), at 3–6× the cost. The deep-channel
noise is therefore not linear gravity-wave dynamics but is re-injected by the
nonlinear switches and limiters.

### 15.7 Limiter activity and Δt-consistency (D1b, 2026-10-25)

Tools: `--limiter-stats` (per-edge flags recorded in the kernels; results
unchanged) and the experimental `--boundary-limiter positivity`. Runs:
Carlisle res-18, inflows + 10 mm/h, 6 h, `--dt-max 10 --output-interval 900`,
e.g.

```
julia --project=. --threads auto FloodModel.jl \
  --meshload test/carlisle/carlisle_mesh18_standard.parquet --flow-model standard \
  --inflow-bci test/carlisle/carlisle_wgs84.bci --rainfall 10 \
  --sim-duration 21600 --dt-max 10 --output-interval 900 [--courant 0.35] \
  [--limiter-stats] [--limiter positivity --q-centre-theta 1 --boundary-limiter positivity]
```

**Which limiters bind** (% of wet edge-steps, by flow depth over the sill):

| | Froude cap (< 0.3 m) | volume limiter (1–3 m / ≥ 3 m) | donor cap (all / ≥ 3 m) | wet/dry fallback | toggles (< 0.3 m / 1–3 m) |
|---|---|---|---|---|---|
| A5 recommended, C 0.7 | 33 | 2.3 / 1.3 | — | 11 | 20 / 1.2 |
| A5 recommended, C 0.35 | 33 | 0.026 / 0.045 | — | 11 | 19 / 0.57 |
| A5 recommended θ = 1, C 0.7 | 34 | 8.0 / 6.3 | — | 11 | 21 / 3.9 |
| Cartesian recommended, C 0.7 | 32 | 0.015 / 0.003 | — | 10 | 17 / 2.3 |
| SGS recommended, C 0.7 | — | — | 54 / 7.1 | — | 0.35 / 0.05 |

("Toggles": wet edges whose set of active limiters changed since the previous
step. For SGS, whose R–A kernel has a different limiter set, only the donor
cap is recorded.) The Froude cap binds on a third of shallow edges whatever Δt is. The
volume limiter (Fix B, q ≤ h·w/(N·Δt)) binds in *deep* water on A5, about
30–90 times more often at Courant 0.7 than at 0.35, and more without Q-centring.
SGS's donor cap (V/(2N) per step) binds on half of all wet edge-steps. All
three per-step limits scale with 1/Δt.

**Courant series with the Δt-dependent terms removed** (storage at 6 h, Mm³):

| configuration | C 0.7 | C 0.5 | C 0.35 | 0.7 → 0.35 |
|---|---|---|---|---|
| A5 recommended (hybrid limiter, Q-centring θ 0.9) | 4.902 | 5.060 | 5.375 | +9.6% |
| A5 recommended, θ = 1 | 4.383 | 3.835 | 3.488 | −20% |
| A5 `--limiter positivity` (θ 0.9) | 4.919 | 5.108 | 5.381 | +9.4% |
| A5 `--limiter positivity`, θ = 1 | 3.990 | 3.777 | 3.391 | −15% |
| … and `--boundary-limiter positivity` | **3.382** | **3.382** | **3.386** | **+0.1%** |
| SGS recommended | 4.334 | 3.952 | 3.602 | −17% |
| SGS `--limiter positivity` | 3.582 | | 3.326 | −7.1% |
| … and `--boundary-limiter positivity` | **3.136** | | **3.138** | **+0.06%** |
| Cartesian recommended | 3.804 | 4.082 | 4.426 | +16% |
| Cartesian positivity, boundary positivity, θ 0.9 | 3.524 | | 4.251 | +21% |
| Cartesian positivity, boundary positivity, θ = 1 | 2.751 | 2.550 | 2.480 | −10% (noise) |

- With Q-centring off (θ = 1), interior Fix B and the donor cap replaced by
  positivity scaling, and the boundary outflow freed of the same two limits,
  **A5 and SGS converge in Δt**. Storage agrees to 0.1%, outflow to 0.4%
  (2523 / 2528 / 2533 k m³), and the flood wave reaches the outflow boundary at
  the same time (arrival of 50 000 m³ outflow: 9938 / 9856 / 9877 s). The rms
  depth difference between Courant 0.7 and 0.35 falls from 39 cm to 4.7 cm,
  below the twin-run difference (6.5 cm).
- Each Δt-dependent term moves the result in its own direction. Q-centring's
  damping rate (1 − θ)/Δt grows as Δt shrinks and slows the wave (arrival
  13 332 → 14 829 s, recommended). The per-step limits bind less as Δt
  shrinks and speed it (θ = 1: 11 677 → 10 218 s). With interior limits
  removed, arrival converges but outflow did not (+32%) until the boundary
  limits were removed too.
- The **boundary limits were also a noise source**: interior positivity,
  θ = 1 has an end-of-run WSE residual r_rms of 343 mm (§15.2 metric) with the
  legacy boundary limits and 52 mm without (twin rms 36 cm → 6.5 cm). The
  outflow boundary of a deep channel is where the per-step volume limiter
  binds hardest.
- The converged configuration is still noisier than the recommended scheme
  (extrema 0.11 vs 0.06, dv_rms 0.65 vs 0.24 m/s), but now the same at every Δt
  (r_rms 52 / 49 mm at Courant 0.7 / 0.35).
- The square grid needs damping: Cartesian with no Q-centring and no limiters
  is chaotic (twin rms 56 cm, r_rms 0.87 m at Courant 0.7, close to its 1/√2
  limit). Its −10% trend is within that noise.
- In the recommended A5 scheme the boundary option changes little (outflow
  1009 → 1025 k m³ at Courant 0.7, none at 0.35): Q-centring dominates there.

**Consequence.** The Courant sensitivity in §15.6 is not a property of the
local-inertial scheme on A5. It comes from three Δt-scaled terms: Q-centring,
the per-edge volume limiter and the donor cap, in the interior and at the
boundary. The Δt-converged answer stores markedly less water at 6 h (3.38 Mm³)
than the recommended scheme at any Courant number (4.9–5.4 Mm³). A Δt-
consistent replacement must provide the damping that the A5 diamond operator
and the square grid need, with a coefficient in physical units rather than per
step. That is D4, the conservative η filter. Positivity scaling with free
boundary outflow is the Δt-consistent limiter baseline it would run on.
