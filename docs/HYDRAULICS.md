# FloodA5 — Hydraulics Reference

_Technical reference for the physics implemented in `FloodModel.jl` and
`surfacewater/flow2d.jl`. Supersedes the 2026-05-12 version of this document,
which predates the stability fixes (§6), the SGS R-A flux path (§9), and all
of the directional-bias correction work (§7–§8). This document reflects the
state of the `directional-bias-reformulation` branch at merge time._

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
`|Q·dt| <= depth_donor · width / N_SIDES` (`N_SIDES = 5`), i.e. no more than
~20% of the donor cell's water can leave via any one edge in a single step.
A coarser, cell-level `DONOR_EDGE_DIVISOR = 10` cap (§11) is retained on top
of this as a last-resort mass-conservation guard; it should rarely bind
once Fix A/B are active.

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
documented in `USER_GUIDE.md` §4.3a, and **none of them affect the SGS
solver** — `step_sgs!` always uses the legacy uncorrected kernel regardless
of these settings.

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
   No reformulation has been attempted yet.
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

**None of the directional-bias corrections in §7 apply to the SGS
solver.** This is a scope boundary, not an oversight: the
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
edge_vol_capped = min(edge_vol, state.volume[donor] / DONOR_EDGE_DIVISOR)
# DONOR_EDGE_DIVISOR = 10 = 2 x N_SIDES -> <=50% total drain across all 5 edges
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
| `--gradient-correction`, `--face-flux-method`, `--momentum-model` | **Yes** | **No effect** (§7, §9.4) |
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

**Full detail and complete experimental history:** this project's internal
development history retains session-by-session records of the full
derivations, proofs, and negative results this document summarises —
including the dead ends and reversed decisions along the way. These
records are not part of this public release, but are retained internally
and available on request. They are kept rather than discarded precisely
because the negative results and dead ends are as load-bearing for future
development as the things that worked.
