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
the R-A flux kernel is active (§9.2).

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

**SGS flow:** `z_sill` is the sill of the face's table (`sgs_edge_sills`):
the lowest point of the DEM profile along the shared side that is connected,
below the water, to the main basins of both cells (§9.1, §9.3). This is what
lets the SGS solver see a channel, notch or ditch narrower than a cell, which
a cell-mean elevation cannot — and, since 0.7.0, a wall narrower than a cell.

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
- SGS (since 0.7.0): the depth is the hydraulic depth V/A_wet, and a second
  limit follows from the faces' flow areas against the wetted area (§9.4).
  Until 0.6.0 the depth was capped at the cell's terrain range
  (`z_max - z_min`), which ignored all water above a cell's highest ground
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
diamond method and cell-vector momentum are ignored under SGS. See §9.5.

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

Rewritten for 0.7.0. The scheme up to 0.6.0 (per-cell edge curves, averaged
between the two cells; legacy donor cap; relief-capped time step) is
described in the 0.6.0 version of this file; its results are in
`cases/expected/v0.6.0/` and its defects in §15.13.

### 9.1 What a cell and a face know

`build_sgs_tables!` (`mesh/A5Grid.jl`) covers every cell with a regular
**working grid** in a local metric frame — spacing
min(DEM pixel, √(bounding box / `--sgs-samples`)), default 1600 points — and
samples the DEM (bilinear) at the grid points inside the polygon. Everything
below is built from that grid and from the DEM profile along each polygon
side (sampled at half a pixel), so a cell and its faces describe one surface.

**Cell: storage.** A hypsometric curve: stored volume V(η) and wetted area at
100 levels between the cell's lowest and highest ground (the extremes taken
over the grid and the side profiles). `wse_from_volume` inverts it; above
the highest ground the level rises with the full cell area. Storage is a
**level pool**: all ground below η holds water, connected or not.

**Cell: where the water is.** The centroid of the wetted area at each level
(`wet_centroid`): the lowest point of the main basin for a dry cell, the
polygon's centroid for a submerged one.

**Face: conveyance.** One table per face, held identically by both cells
(`face_table_from_profile`, `face_apk`). From the profile z_p of the face
(points of width dx) and an activation level a_p ≥ z_p per point (§9.3):

    A(η) = Σ dx·(η − z_p)⁺·[η ≥ a_p]          flow area
    P(η) = Σ dx·[η ≥ a_p]                      wetted width
    K(η) = Σ dx·((η − z_p)⁺)^(5/3)·[η ≥ a_p]   conveyance integral

at 64 levels from the face's sill (its lowest activation level) to above its
highest point, spaced quadratically (dense near the sill, looked up without a
search). Above the top level every point is wet: A grows with the full width
W, P = W, and K = W·(d + c)^(5/3), d the height above the top level and c
fixed by K there. Nothing stops growing with stage. Polygon sides on the
mesh boundary have the same table, used by the open-boundary outflow.

### 9.2 The face flux

The R–A form of the local-inertial momentum equation (LISFLOOD-FP sub-grid
channel form), for the discharge Q through the face:

    Q ← (Q − g·A·Δt·Δη/L) / (1 + g·Δt·n²·|Q| / (R^(4/3)·A))

with A and R evaluated at the higher of the two cells' levels.

**Conveyance form** (`--sgs-conveyance`, `_sgs_radius`). `2d` (default):
R = (K/A)^(3/2), so the steady flux is K·√S/n with K = Σ dx·h^(5/3) — every
strip of the face carries its own Manning flow. This is the uniform-flow
limit of the 2-D equations, i.e. what a grid that resolves the face
converges to. `lumped`: R = A/P, one hydraulic radius for the whole face;
overbank this gives the level of a cell-mean bed (strip: +0.04 m, +8%
storage) and is kept for the record.

**Flow-path geometry** (`--sgs-path`, `_sgs_path`). The head difference
Δη is between the water bodies of the two cells. For partly wet cells those
are around their wetted-area centroids, not the cell centres: a channel that
crosses the mesh obliquely has its water on a line u shorter than the line
through the centres, and faces cut it obliquely, wider than it is. For flow
along u through a face with normal n̂ the kernel uses

    A·|n̂·û|   and   L = |u|.

Fully wet cells keep the kernel's own geometry (A; L projected on the face
normal); the two are blended by the mean wetted fraction w of the two cells,
s = clamp(2w − 1, 0, 1), so cells less than half wet are treated as channels.
On the A5 strip in bank this takes the level from 0.054 m above the analytic
one to 0.009 m below it (§15.13).

**Critical-flow limit** (`--sgs-critical`, default on since 0.8.0;
`_sgs_critical_q`). The local-inertial equations have no advection term, so
nothing in them stops a face passing more than critical flow where the water
falls away from it — a notch, a weir, a free overfall, an outflow edge. Until
0.7.x SGS passed ≈ 9.8 m³/s through the 4 m notch of `planar_embankment` at
1.07 m of head, against ≈ 7.5 m³/s at critical flow. The kernel now caps

    |Q| ≤ √g·(2/3)^(3/2)·C(E)·a,     C(E) = Σ dx·(E − a_p)^(3/2)

with E the upstream cell's level (still water: the approach velocity head is
neglected), a_p the points' activation levels and a the flow-path factor
above: each strip of the face passes at most its own critical discharge with
the head E − a_p over it (the counterpart for critical flow of K = Σ h^(5/3)
for friction). For a rectangle of width W it is the broad-crested-weir
discharge 1.705·W·H^(3/2). C is tabulated per face like K (`sgs_face_C`,
layout 4). The open-boundary outflow of a side is capped the same way. On a
slope the head includes the drop from the upstream level to each point of
the face, so flow down a moderate slope is not held back; but the head is
the still-water level (no velocity head), so on a steep reach the cap can
bind on supercritical flow: on the 2.1% outlet reach of `sgs_strip` it holds
1–2% more of the strip's water (§15.14). The limit is
a choice of the hydraulics, not of the SGS geometry: `off` gives the 0.7.x
flux on the same tables.

**Dry cells.** A dry cell's level for the flux is max(sill, its lowest
ground): water must clear the face and reach the cell's floor.

### 9.3 Connectivity: activation levels

Storage is a level pool, so by itself a barrier *inside* a cell only removes
volume: a face on either side of it has its own profile at natural ground,
and water crosses the cell without meeting the barrier (the "hop"; a wall
narrower than a cell passed water 2.9 m below its crest in case
`planar_wall`).

`_cell_connectivity` finds, on the working grid, the cell's **main basin**:
the first region of ground below a rising level (4-connected) that reaches
two different sides of the cell — the lowest route through it — and for
every grid point the lowest level at which it is connected to that basin
through ground below that level. A face point's activation level is the
larger of its own elevation and that connection level in *both* cells, so it
conveys only when the water at it is connected to the water bodies of the
two cells it joins. Consequences:

- a barrier inside a cell closes the faces on the side cut off from the main
  basin, up to the level at which water gets over or round it; its position
  is moved by at most one cell;
- a gap through the barrier is found by the search and conveys at its own
  sill (case `planar_embankment`);
- a hollow cut by a face that leads nowhere does not convey.

`--sgs-connectivity off` (mesh generation) gives every point its own
elevation. `sgs_face_raise` in the mesh records by how much each face's sill
was raised.

**Limit.** One storage per cell is one basin per cell. Where a cell holds two
basins of comparable importance (a river and a parallel side channel behind
a levee), the second is merged in storage and closed in conveyance until the
levee level. 5–8% of Carlisle res-18 cells (5 m DEM; a lower bound) have a
second basin over 5% of the cell and deeper than 0.25 m (Phase 0 audit,
`docs/notes/SGS_Phase0_Audit_2026-11.md`). Splitting such cells into one
storage per basin is the next step.

### 9.4 Limiter and time step

- **Limiter:** positivity (one scale factor per donor cell, applied only when
  the cell's outflow would exceed its volume), in the interior and at open
  boundaries, under `--scheme recommended`. The legacy donor cap, V/(2N) per
  face per step, makes a channel cell carry Q only once it holds 2N·Q·Δt and
  so holds in-bank channels brim-full; `--limiter legacy` still selects it.
- **Time step** (`_cfl_dt_local`): two limits. (i) The hydraulic depth
  V/A_wet in the usual ℓ/√(g h) form (`_sgs_cfl_depth`). (ii) From the faces
  (`_sgs_face_rate`): a level change of a cell changes its net inflow in
  proportion to Σ A_f/L_f over its faces, so its level oscillates at
  ω = √(g·Σ(A_f/L_f)/A_wet) and the limit is Δt ≤ 2·courant/ω. For a fully
  wet square cell (ii) is (i). For a cell that holds only a sliver of a
  channel its faces convey, (ii) is shorter by √(cell length / sliver
  length); without it such cells oscillate and hold the channel up (0.5 m
  between neighbouring cells on the A5 strip). The wetted area is floored at
  2% of the cell so that a wetting front does not set the step.

### 9.5 Interaction with §7/§8

Only Correction 1 (WLSQ + V̂) reaches the SGS solver, and it is unvalidated
there (off under `--scheme recommended`); the diamond method and cell-vector
momentum do not apply. Vector friction and the local time step do.

---

## 10. Velocity Computation

`state.velocity` (scalar magnitude), `state.vel_u`, `state.vel_v`
(directional components) are computed after the flux/volume update each
step from a flux-weighted average across each cell's edges, and written to
HDF5 output (`/results/velocity`, `vel_u`, `vel_v`; `/frames/{idx}/…` before 0.4.0). Per-cell velocity
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
| Sill | `max(elev_i, elev_j)` | The face's lowest activation level (`sgs_edge_sills`; §9.1, §9.3) |
| Water-depth pre-sync needed | Yes | No |
| Requires SGS tables | No | Yes (`build_sgs_tables!` at mesh-build time) |
| Stability limiters | Fix A/B/C (§6) | Positivity scaling; dry-cell level; the face-area time-step limit (§9.2, §9.4) |
| `--gradient-correction` | **Yes** | Legacy WLSQ + V̂ only, **unvalidated** (§9.5) |
| `--face-flux-method`, `--momentum-model` | **Yes** | **No effect** (§7, §9.5) |
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
`test/archive/run_wp1_matrix.sh`, `test/wp1_metrics.jl spread`), four-fold anisotropy A4:

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
unvalidated under SGS, §9.5); hybrid alone spiked (velocity 4 m/s) — under SGS
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
`test/diag_divergence.jl`, `test/archive/diag_d2_summary.sh`.

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

### 15.8 Δt-consistent damping: the η filter (D4, 2026-10)

**Why the square grid needs damping.** Linearised and without friction, the
local-inertial step is a forward–backward integration of the wave equation:
q is updated from the old surface, then the surface from the new q. For a
mode whose free-surface operator eigenvalue is λ, with κ = g·h·Δt²·λ, the
amplification factors solve μ² − (2 − κ)μ + 1 = 0. For κ ≤ 4 both have
|μ| = 1 exactly. Every resolved mode is neutral and the scheme has **no
numerical dissipation**; it turns through a phase φ per step, with
cos φ = 1 − κ/2. Friction is the only sink. Manning friction is quadratic, so
it vanishes to first order around slow, deep water: in a 3 m channel at
0.5 m/s with n = 0.035, a perturbation decays with an e-folding time of
6–12 minutes (along or across the flow), against a time step of a few
seconds, and in still water not at all.

On a square grid the stiffest mode is the 2Δx checkerboard, λ = 8/Δx², so
κ = 8C² and the scheme's limit is C* = 1/√2 = 0.707 (`lam.jl` on the Carlisle
Cartesian mesh: 0.707). At the LISFLOOD-FP value C = 0.7 the checkerboard sits
at κ = 3.92, 98% of the limit, and flips sign almost every step (φ = 0.91π).
A5 has no checkerboard: its cell graph is not bipartite (where three cells
meet at a vertex they are mutually adjacent, so the graph has triangles, and
no field can alternate in sign across every face). Its stiffest mode lies lower (C* = 0.876 on Carlisle
res-18), so at C = 0.7 κ = 2.55 and φ = 0.59π.

Meanwhile the discrete switches inject grid-scale perturbations every step:
the Froude cap binds on a third of shallow edges, and 16–20% of shallow edges
change their limiter state from one step to the next (§15.7). With nothing to
remove it, that energy accumulates in the neutral grid-scale modes. On the
square grid near C* the checkerboard is also close to resonance with
anything that alternates step by step. Measured on Carlisle with no
Q-centring and no per-step limiters (§15.7): the Cartesian end-of-run WSE
residual is 0.87 m at Courant 0.7 and 0.21 m at 0.35 (A5: 52 mm at both),
and its twin-run divergence is 56 cm. This is the known LISFLOOD-FP
behaviour that de Almeida et al. (2012) addressed with their θ weighting of
the discharge, which Q-centring implements here.

**Why damping q is the wrong place, and η the right one.** Let a damping
remove a fraction β of a mode per step. Applied to the discharge, as
Q-centring does (q′ = (1 − β)q − …), the 2×2 amplification matrix has
det = 1 − β and trace = 2 − β − κ. The Jury conditions then give
κ < 4 − 2β: **the damping lowers the stability limit**. On the square grid at
Courant 0.7 (κ = 3.92), any β > 0.04 on the checkerboard is unstable. That is
the |μ| = 1.6 of Q-centring θ = 0.9 found in §15.6. Applied instead to the
water surface after the continuity update (η′ = (1 − β)·η⁺), det = 1 − β and
trace = 1 + (1 − β)(1 − κ), and the conditions become β < 1 and
κ < (4 − 2β)/(1 − β). **The limit rises.** Damping the surface stabilises
the very mode that damping the discharge destabilises, provided the damping
per step stays below 1.

**The filter** (`--eta-filter C_s`, `surfacewater/eta_filter.jl`; the
"conservative free-surface high-pass filter" of the damping recommendations).
After the volume update, each cell's high-frequency residual is computed:
r_i = η_i − a_i, where a_i is the cell-centre value of a weighted planar fit
to the cell and its wet neighbours. Water is then exchanged across every wet
face:

  F_f = C_s · √(g h_f) · w_f · (r_j − r_i),   h_f = min(η_i, η_j) − z_sill > 0.

- *Conservative*: the exchange is pairwise; a donor cell's filter outflows
  are scaled so that it cannot lose more than it holds.
- *Well-balanced and consistent*: constants and planes have r = 0, so a lake
  at rest and a uniform slope are untouched (`test/test_eta_filter.jl` E2,
  E6), and smooth surfaces barely move (r ∝ ℓ²·curvature). It acts like a
  fourth-order (∇⁴) operator, not a diffusion.
- *Δt-consistent*: the rate C_s·√(gh)/ℓ is in physical units. A field is
  damped by the same amount over the same time whatever Δt is (E4).
- *Wet/dry*: neighbours enter the fit with weight ramp clamp(h_f/5 cm, 0, 1),
  so a dry bank never enters it and nothing switches discontinuously; no
  water crosses a face whose sill is above either side's surface.
- *Strength*: on a square grid r = ±1.6δ for a checkerboard of amplitude δ,
  and it decays by 1 − 12.8·C_s·C per step (E3). β < 1 then requires
  C_s < 0.11 at Courant 0.7.

**Linear stability** (`test/diag_step_spectrum.jl`, lake at rest, 3 m,
Courant 0.7; positivity limiter, θ = 1):

| C_s | Cartesian max\|μ\| (growing modes) | A5 max\|μ\| (growing modes) |
|---|---|---|
| 0 | 1.000000 (1052 neutral) | 1.001420 (304) |
| 0.02 | 1.000000 (0) | 1.000000 (0) |
| 0.05 | 0.999999 (0) | 1.000000 (0) |
| 0.08 | 0.999998 (0) | 1.000000 (0) |
| 0.12 | 1.240578 (8) | 1.000000 (0) |

(Recommended scheme with Q-centring θ = 0.9 on the square grid: 1.601.) The
filter removes the diamond operator's slow growth on A5 already at
C_s = 0.02. On the square grid it is stable up to the predicted bound of
0.11. On a sloping, friction-dominated base flow (A5, 1‰), every C_s is
stable.

**Uniform-flow conveyance** (frozen-coefficient test, planar res-16, slope
10⁻³, n = 0.03, mean conductance relative to Manning / maximum mean
deflection over 12 bearings):

| | h = 1 m, Δt = 2 / 8 s | h = 3 m, Δt = 2 s | h = 3 m, Δt = 8 s |
|---|---|---|---|
| recommended (Q-centring θ 0.9, hybrid limiter) | 1.000 / 0.2° | 0.999 / 0.2° | **0.889 / 8.0°** |
| recommended, θ = 1 (hybrid limiter) | 1.000 / 0° | 1.000 / 0° | **0.927 / 4.0°** |
| positivity, θ = 1 (± η filter 0.05) | 1.000 / 0° | 1.000 / 0° | 1.000 / 0° |

In deep water at a realistic time step, the volume limiter
(q ≤ h·w/(N·Δt)) removes up to 11% of the Manning conveyance, depending on
direction. That is the steady-flow face of the Courant sensitivity in
§15.7. Without it the scheme is exact, and the filter does not change that:
the surface is planar, so r = 0.

**Carlisle** (res-18 A5 and 20 m Cartesian, inflows + 10 mm/h, 6 h;
`--limiter positivity --q-centre-theta 1 --boundary-limiter positivity
[--eta-filter C_s]`; runs in `test/wp1_runs/d4/`):

| | storage at 6 h, C 0.7 / 0.5 / 0.35 (Mm³) | twin rms | rms C 0.7 vs 0.35 | extrema | r_rms |
|---|---|---|---|---|---|
| A5 recommended | 4.902 / 5.060 / 5.375 | 2.3 cm | 17.5 cm | 6.0% | 31 mm |
| A5, no filter | 3.382 / 3.382 / 3.386 | 6.5 cm | 4.7 cm | 10.6% | 52 mm |
| A5, C_s = 0.02 | 3.382 / 3.388 / 3.388 | 3.2 cm | 4.2 cm | 7.4% | 51 mm |
| A5, C_s = 0.05 | 3.381 / 3.385 / 3.386 | 2.8 cm | 3.5 cm | 5.8% | 54 mm |
| Cartesian recommended | 3.804 / 4.082 / 4.426 | 2.5 µm | — | 1.8% | 43 mm |
| Cartesian, no filter | 2.751 / 2.550 / 2.480 | 56 cm | 46 cm | 21% | 872 mm |
| Cartesian, C_s = 0.02 | 2.559 / 2.487 / 2.378 | 26 cm | — | 5.7% | 299 mm |
| Cartesian, C_s = 0.05 | 2.447 / 2.356 / 2.354 | **2.3 µm** | 12.5 cm (0.5 vs 0.35: 4.4 cm) | 1.5% | 68 mm |

(Extrema and r_rms at Courant 0.7, end of run, §15.2 metric.)

- On A5 the filter keeps the Δt convergence of §15.7 (storage within 0.2%)
  and brings local extrema down to the recommended scheme's level. Twin-run
  and Courant divergence are both about 3 cm, against 17.5 cm Courant
  divergence for the recommended scheme.
- On the square grid, C_s = 0.05 turns a chaotic run into a reproducible one:
  twin divergence is at round-off (2.3 µm, as for the recommended Cartesian
  scheme) and extrema are lower than recommended. Storage converges between
  Courant 0.5 and 0.35. At 0.7, just under the square grid's limit of 0.707,
  it is 4% higher. C_s = 0.02 is too weak for the square grid.
- Mass balance is exact in every run, and the filter costs 15–22% in step
  throughput (a full pass over all cells; it could use the active-edge set).
- Deep-channel speed (cells ≥ 3 m deep): 1.34 m/s for A5 and 1.39 m/s for
  Cartesian (Δt-consistent, filtered), against 0.42 m/s for recommended A5.
  The 6 h storage difference (3.4 vs 4.9–5.4 Mm³) is mainly this extra
  resistance from the Δt-scaled terms.

**Planar slope, res-18.**
- *Point source* (closed boundaries, 20 h, Δt ≤ 2 s): N/S asymmetry +0.176 for
  the baseline, with or without the filter (unchanged to 10⁻⁴), against
  +0.193 for recommended.
- *Rained* (20 mm/h, open boundaries, 6 h, Courant 0.7 / 0.35): interior mean
  depth 0.0592 / 0.0596 m (baseline), 0.0603 / 0.0602 m (C_s = 0.05) and
  0.0649 / 0.0628 m (recommended), against the Δt-converged 0.059 m of §15.4.
  Asymmetry is below 0.01 for all. Local extrema 6.2% → 4.3% (C_s = 0.05) at
  Courant 0.7, against 1.3% for recommended; interior velocity noise 0.12 vs
  0.05 m/s.

**What the filter does not do.** The remaining noise is in the velocity
field. These are discharge patterns that circulate without moving the water
surface (one per independent cycle of the face graph, E − N ≈ 1.5 per cell on A5
with its 2.5 faces per cell, against 1 on a square grid), so a surface filter cannot see them. Q-centring damped them,
at the cost of Δt-dependence and conveyance. On Carlisle part of the A5
velocity "noise" is shear between fast channel cells and slow bank cells.
Relative to the mean speed it rises only from 0.78 to 0.89 in cells over 3 m
deep (Cartesian 0.19 → 0.30). On the shallow rained slope it is a genuine
2.4× increase over recommended.

**Status.** The filter is opt-in; defaults are unchanged. The Δt-consistent
candidate for standard flow, `--limiter positivity --q-centre-theta 1
--boundary-limiter positivity --eta-filter 0.05`, meets these criteria:
- Δt convergence;
- reproducibility at the recommended scheme's level (Cartesian at round-off);
- exact uniform-flow conveyance;
- Carlisle surface extrema at the recommended scheme's level.

It does not yet meet the noise criterion for velocity, or for surface extrema
on shallow sheet flow. The missing piece is a Δt-consistent damping of the
circulating discharge modes. That is face-momentum diffusion acting on the
edge discharges, with a coefficient in physical units and vanishing for
uniform velocity: item 4 of the damping recommendations.

### 15.9 Δt-consistent momentum damping (D6′, 2026-10)

The η filter (§15.8) cannot see discharge patterns that circulate without
moving the water surface. `--momentum-damping C_m`
(`surfacewater/momentum_damping.jl`) damps those, and only those:
- Γ_v is the circulation of the face velocities around each interior mesh
  vertex. The velocity used is the component along the centre-to-centre
  line, which includes the face-tangential part from the cell vectors, so a
  uniform flow has Γ = 0 exactly on a skewed mesh too.
- A streamfunction ψ_v = −ν·Δt·h_v·Γ_v/A_v is formed at each vertex.
- Each face's stored discharge changes by ψ_head − ψ_tail. Around every
  cell these changes telescope to zero, so the correction is **exactly
  divergence-free**: it moves no water, and it leaves the gravity-wave
  dynamics and their stability limit untouched. The full-step spectra are
  unchanged (`test/diag_step_spectrum.jl`).
- ν = C_m·V·ℓ_v, with ℓ_v = √A_v and V = √(gh) (`wave`, a numerical floor)
  or the local flow speed |u| (`flow`, eddy-viscosity-like). The rate is in
  physical units, so the result does not depend on Δt.
- Wetness weights ramp the correction off near dry faces, and vertices on
  the domain edge are left alone. It requires the diamond face gradient and
  vector friction (cell vectors every step), as in the recommended scheme.
- Tests: `test/test_momentum_damping.jl`, M1–M5 (exact divergence-freeness,
  uniform flow untouched on A5 and Cartesian, Δt-independent vortex decay,
  rest and mass in the full step).
- Uniform-flow conveyance (frozen test, h = 3 m, Δt = 8 s): 0.999 × Manning
  with 0.1° deflection on A5, and exact on Cartesian. The recommended scheme
  gives 0.889 and 8°.

**Carlisle**, A5 and Cartesian side by side. All runs use `--limiter
positivity --q-centre-theta 1 --boundary-limiter positivity --eta-filter
0.05` plus the momentum damping shown; runs are in `test/wp1_runs/d6/`.
Columns: storage at 6 h; twin-run and Courant 0.7 vs 0.35 rms depth
difference; end-of-run local extrema and velocity noise dv_rms (§15.2
metric, Courant 0.7); mean speed in cells over 3 m deep.

| mesh | momentum damping | storage C 0.7 / 0.35 (Mm³) | twin rms | Courant rms | extrema | dv_rms | speed ≥ 3 m |
|---|---|---|---|---|---|---|---|
| A5 | recommended scheme | 4.902 / 5.375 | 2.3 cm | 17.5 cm | 6.0% | 0.24 | 0.42 |
| A5 | none (η filter only) | 3.381 / 3.386 | 2.8 cm | 3.5 cm | 5.8% | 0.65 | 1.34 |
| A5 | wave 0.05 | 3.444 / 3.445 | 2.5 cm | 2.9 cm | 6.1% | 0.30 | 0.95 |
| A5 | wave 0.2 | 4.202 / 4.186 | 2.5 cm | 4.7 cm | 5.5% | 0.18 | 0.66 |
| A5 | flow 0.5 | 4.089 / 4.071 | 7.4 cm | 6.3 cm | 6.1% | 0.25 | 0.84 |
| Cartesian | recommended scheme | 3.804 / 4.426 | 2.5 µm | — | 1.8% | 0.11 | 0.67 |
| Cartesian | none (η filter only) | 2.447 / 2.354 | 2.3 µm | 12.5 cm | 1.5% | 0.25 | 1.39 |
| Cartesian | wave 0.05 | 2.682 / 2.623 | 6 µm | 11.5 cm | 1.4% | 0.16 | 1.10 |
| Cartesian | wave 0.2 | 2.896 / 2.938 | 8×10⁻¹¹ m | 4.8 cm | 2.3% | 0.14 | 1.00 |
| Cartesian | flow 0.5 | 2.739 / 2.750 | 3×10⁻¹⁰ m | 5.6 cm | 1.8% | 0.15 | 1.08 |

(The A5 and Cartesian meshes differ in resolution and in how the DEM is
sampled, so their storages are not directly comparable; compare changes
within a mesh.)

**Rained planar slope** (res-18, 20 mm/h, Courant 0.7), against the
recommended scheme (in brackets):
- interior mean depth 0.0605–0.0608 m (0.0649; Δt-converged 0.059);
- extrema 3.1% / 2.7% / 2.7% for wave 0.05 / wave 0.2 / flow 0.5 (1.3%);
- dv_rms 0.10 / 0.054 / 0.063 m/s (0.052);
- N/S asymmetry +0.011 / +0.028 / +0.025 (−0.005; η filter only +0.006).

**Findings.**
- The circulating modes are the velocity noise. Damping them removes it
  without touching the water-surface dynamics. At wave 0.2 it is at or below
  the recommended scheme on every mesh (A5 0.18 vs 0.24 m/s), and the
  runs stay Δt-converged (A5 storage within 0.4%, Cartesian within 2.2%).
  Cartesian runs become reproducible to 10⁻¹⁰ m.
- **It is an eddy viscosity, and its strength is a modelling choice with a
  first-order effect.** Lateral momentum exchange slows the channels (deep
  speed 1.34 → 0.66 m/s on A5) and raises 6 h storage, by 2% (wave 0.05) to
  24% (wave 0.2) on A5 and 11% to 25% on Cartesian (at Courant 0.35). In 3 m of water on these
  meshes ν ≈ 5–20 m²/s for C_m = 0.05–0.2, one to two orders of magnitude
  above typical depth-averaged turbulent viscosities. The recommended
  scheme's Q-centring adds the same kind of resistance implicitly, more
  strongly and Δt-dependently (deep speed 0.42 m/s).
- **Polygon dependence.** For the strong settings the relative effect is
  similar on both meshes (wave 0.2: +24% A5, +25% Cartesian; flow 0.5: +20%,
  +17%). For wave 0.05 it is not (+2%, +11%). On A5 the vertex operator adds
  a small directional bias on the planar slope (asymmetry up to 0.03). The
  vertex fans are irregular there (three or four cells per vertex), and the
  operator weights them by dual area only.
- Cost: 10–30% more step time on top of the filter.

**Consequence.** With the η filter and the momentum damping, a
Δt-consistent scheme is reproducible and convergent, and its noise matches
the recommended scheme. What it cannot decide is how much lateral momentum
exchange is right. The local-inertial equations have neither advection nor
a turbulence closure, so that amount is set by C_m, and the answer moves
6 h storage by up to a quarter. Settling it needs a reference that resolves
these terms: a full shallow-water (HLL/HLLC) solver on the same polygon
meshes, with a physical eddy viscosity.

### 15.10 A full shallow-water reference solver (R1, 2026-10; corrected in 0.2.6)

`--flow-model swe` (`surfacewater/swe_reference.jl`, `swe_types.jl`) solves
the full 2-D shallow-water equations, advection included, on the same meshes,
beds, sources and boundaries as the local-inertial solvers. It exists to be
the dissipation target that §15.9 asked for: a scheme whose momentum
transport and numerical dissipation come from the equations and a standard
Riemann solver, not from a tuned damping coefficient.

> **0.2.6 correction.** The 0.2.5 solver had two defects on thin films over
> a slope (rain-on-grid), both found by comparing with analytic solutions:
> the second-order hydrostatic reconstruction used cell-centre beds at the
> faces, and the semi-implicit friction took |q| from the predicted state.
> Together they slowed shallow sheet flow by 36–80% and made the reference
> too dissipative: 6 h Carlisle storage fell from 3.44 to 2.69 Mm³ (A5) and
> from 2.63 to 2.33 Mm³ (Cartesian) when they were fixed. The comparisons
> below are with the corrected solver; the 0.2.5 conclusions about the
> momentum damping (that C_m = 0.05 matched the reference) do not hold.

**Scheme.** A cell-centred finite-volume method on arbitrary polygons:
- HLLC flux (Toro) in each face's normal–tangential frame, with the
  two-rarefaction wave-speed estimate and the dry-bed limits.
- Hydrostatic reconstruction (Audusse et al. 2004) at the face bed
  z_f = max(z_iL, z_jR), with each side's face bed z_iL = η_iL − h_iL taken
  from the reconstructed states (z_i at first order) and the centred bed
  source ½g(h_iL + h_i)(z_iL − z_i) — the second-order scheme of Audusse et
  al. (2004, §4). An edge sill above both cell beds also bounds z_f. This
  keeps depths non-negative and a lake at rest exactly at rest. Each cell's
  pressure correction is taken relative to its own ½gh_i², so the balance
  holds even where the polygon's face normals don't close exactly (the A5
  cells are spherical).
- Order 2 (default; `--swe-order 1` for first order):
  - MUSCL with weighted-least-squares gradients of η, h, u, v and the
    Barth–Jespersen limiter (the face depth comes from h's own limited
    reconstruction, so it stays within the neighbours' depths);
  - first order next to dry cells;
  - SSP-RK2 in time.
- Manning friction semi-implicit in each stage, with |q| from the stage's
  start, so that its fixed point is exactly Manning's normal flow. An
  optional eddy viscosity ν_t = C_E·h·u* (`--swe-viscosity C_E`, default 0).
- Boundaries:
  - Closed: a reflecting wall;
  - Critical;
  - ZeroGradient (and FixedWSE/FixedQ until implemented): transmissive for
    outflow only, matching the standard solver. A fully transmissive variant
    is used only in benchmarks (`swe_open_inflow`).
- Δt = C·min_i A_i/Σ_f w_f·s_f, with the face wave speed s_f the faster of
  the two cells. `step_swe!` sub-steps internally when the run loop's Δt
  (chosen before the sources are applied) exceeds the stable step.

Volume is `state.volume` and the cell discharge (hu, hv) is `state.qvec_u/v`,
so rainfall, inflows and outflow accounting are the ones the other solvers
use. Every pass is a per-face kernel or a per-cell gather over a fixed CSR
index, with no atomics, Dicts or Symbols and with Int8 boundary codes, the
layout a GPU backend needs.

**Verification** (`test/test_swe_reference.jl`, A5 and Cartesian):
- Lake at rest over an island stays at rest to round-off (S1). Closed
  sloshing conserves mass to round-off and stays positive (S2).
- Stoker dam break: L1 depth error 0.0219 → 0.0131 (order 1) and
  0.0094 → 0.0042 (order 2) as the grid halves (S3). Ritter dry-bed front
  within 5% (S4).
- Open boundaries: storage + outflow = initial volume (S5).
- Uniform Manning flow at 1 m holds normal depth and velocity: A5 0.14% /
  0.4%, Cartesian 0.05% / 0.15% (S6).
- Thin sheet flow at 2 cm, where the bed step across a face exceeds the
  depth: Cartesian 0.04% / 0.1%, A5 within 8% (S10). Order 1 cannot do
  this (velocities 37–67% low): first-order hydrostatic reconstruction is
  not a rain-on-grid reference.
- Rained 0.1% plane, 20 mm/h, 6 h: the interior depth profile matches the
  kinematic steady state h = (r·x·n/√S)^0.6 within 1–3 mm on both meshes.
- The eddy viscosity removes shear energy (S8). HDF5 frames carry the SWE
  velocity (S9).

**Polygon consistency on a smooth problem** (`test/bench_swe_dambreak.jl`).
A circular dam break, 60 s, flat frictionless bed, on A5 res-16 and on a
Cartesian grid of the same cell area, each compared with a Cartesian grid
four times finer:
- Mean |h − fine|: A5 0.097 m and Cartesian 0.092 m at order 1;
  0.063 m and 0.056 m at order 2.
- Azimuthal spread (0 for an isotropic scheme): 0.034 and 0.020 m at order
  1; 0.045 and 0.036 m at order 2.

The two meshes are equally accurate to within about 10%. A5 is slightly less
isotropic.

**Carlisle** (6 h, 10 mm/h rain, `--dt-max 10`; runs in `test/wp1_runs/r1b/`).
The reference converges in its own numerical parameters:

| | A5 res-18 | Cartesian 20 m |
|---|---|---|
| storage, order 2, C 0.7 (Mm³) | 2.690 | 2.331 |
| C 0.35 | 2.690 (rms Δh 0.2 cm) | 2.331 (1.6 cm) |
| twin (rain × (1 + 10⁻¹⁰)), rms Δh | 3×10⁻⁴ m | — |
| `--swe-viscosity 0.3` | 2.690 (0.2 cm) | 2.332 (0.1 cm) |
| order 1 | 3.776 | 2.911 |
| arrival (s; 50 000 m³ out) | 6 688 | 5 736 |

Clamped volume is 0 and the mass-balance error a few m³ in every run.
Second- against first-order is a 25–40% storage difference: first order is
not a reference here. A CPU run costs about 8× the steps of the
local-inertial schemes, because the step is bounded by |u| + √(gh) summed
over the faces; on the GPU (`--backend gpu`) the 10 m Cartesian reference
(146 000 cells) took 36 minutes.

**Local-inertial schemes against the reference, same mesh** (frame 24,
t ≈ 20 700 s; `test/diag_vs_reference.jl`).
- Columns: storage ratio; rms depth difference over wet cells; CSI and bias
  of the extent at 0.1 m; mean speed ratio in cells over 3 m deep.
- Every Δt-consistent run uses `--limiter positivity --q-centre-theta 1
  --boundary-limiter positivity --eta-filter 0.05` plus the momentum damping
  shown. "0.2 scheme" is the 0.2.x default.

| mesh | run | storage | rms Δh | CSI | bias | speed ≥ 3 m |
|---|---|---|---|---|---|---|
| A5 | SWE order 1 | 1.394 | 0.35 m | 0.56 | 1.21 | 0.75 |
| A5 | 0.2 scheme | 1.802 | 0.51 m | 0.48 | 1.45 | 0.40 |
| A5 | η filter only | 1.260 | 0.29 m | 0.61 | 1.07 | 1.18 |
| A5 | wave 0.05 | 1.280 | 0.31 m | 0.60 | 1.09 | 0.83 |
| A5 | wave 0.2 | 1.540 | 0.46 m | 0.53 | 1.28 | 0.60 |
| A5 | flow 0.5 | 1.504 | 0.46 m | 0.55 | 1.25 | 0.73 |
| Cartesian | SWE order 1 | 1.244 | 0.27 m | 0.62 | 0.99 | 0.88 |
| Cartesian | 0.2 scheme | 1.620 | 0.47 m | 0.54 | 1.23 | 0.50 |
| Cartesian | η filter only | 1.049 | 0.19 m | 0.68 | 0.86 | 0.99 |
| Cartesian | wave 0.05 | 1.144 | 0.21 m | 0.67 | 0.90 | 0.81 |
| Cartesian | wave 0.2 | 1.236 | 0.26 m | 0.64 | 0.96 | 0.73 |
| Cartesian | flow 0.5 | 1.168 | 0.22 m | 0.66 | 0.92 | 0.80 |

- **Every local-inertial scheme is more dissipative than the reference,
  and the momentum damping only adds to it.** The η filter alone is
  closest on both meshes (+5% storage on Cartesian, +26% on A5); each
  damping setting stores more. The 0.2 scheme stores 62–80% more.
- **On A5 the local-inertial excess is larger.** It is not at the
  boundaries (more than 20 cells from the domain edge: +668 of +691 k m³)
  and is spread over every depth class, the largest part (about 200 k m³)
  on cells the reference has dry. SWE order 1 shows the same pattern, more
  strongly (+322 k m³ on those cells). Both use a step bed with face sills
  max(z_i, z_j); the second-order reference reconstructs a sloping face bed,
  across which thin films leave shallow hollows that a step bed holds. So
  on Carlisle at ~20 m much of the remaining gap is the bed representation,
  not the dissipation.

**Resolution.** A Cartesian refinement (10, 20 and 40 m grids, all averaged
from the 5 m LiDAR), each run with the reference and the local-inertial
schemes, and a cross-mesh map to the 10 m reference (nearest cell centre):

| run | 40 m | 20 m | 10 m | CSI vs SWE 10 m (20 m run) |
|---|---|---|---|---|
| SWE, order 2 | 3.154 | 2.331 | 2.107 | 0.53 |
| η filter only | 3.268 | 2.447 | 2.278 | 0.56 |
| η filter + wave 0.05 | 3.572 | 2.682 | 2.420 | — |
| 0.2 scheme | 4.553 | 3.804 | 3.833 | — |

(6 h storage, Mm³; arrival of the reference 8 971 / 5 736 / 5 031 s.)

- The reference converges quickly (differences 0.82, then 0.22 Mm³;
  extrapolated about 2.03 Mm³). The η filter converges alongside it,
  4–8% above (extrapolated about 2.24: the local-inertial limit), and at
  20 m its flood extent matches the fine reference as well as the 20 m
  reference does (CSI 0.56 against 0.53). Momentum damping adds a further
  6–10%. The 0.2 scheme does not converge at all (3.80 → 3.83): its
  Q-centring resistance grows as Δt shrinks with the cells.
- A5 res-18 (cells of a 22 m square) sits between the 20 and 40 m Cartesian
  references: storage 2.690 and CSI 0.47 against the 10 m reference place it
  at about 28–29 m. The smooth dam break shows no such gap, so it comes from
  channels narrower than a cell (the Petteril and Caldew are one or two cells
  wide); a straight-channel series across orientations and resolutions
  tests that directly.
- Swapping beds between the meshes (each cell takes the bed of the nearest
  cell of the other mesh) breaks narrow channels whichever way it goes, so
  bed and cell geometry cannot be separated that way.

**Consequences.**
- On a given mesh and bed, the Δt-consistent scheme with the η filter alone
  is the closest local-inertial configuration on both polygon types and at
  every resolution tested; the momentum damping moves it away from the
  reference (its case rests on velocity noise alone, §15.9), and the 0.2
  scheme is far from it and not resolution-convergent.
- The remaining local-inertial excess on A5 is mostly a bed-representation
  effect shared with first-order SWE, and the A5/Cartesian difference on
  Carlisle is as large as a factor-two change in resolution. Polygon
  consistency of the hydraulics is best judged on synthetic problems, and of
  the topography on resolution sequences and channel tests.

### 15.11 Defaults in 0.3.0 (2026-10)

From 0.3.0, `--scheme recommended` for **standard flow** is the Δt-consistent
scheme, and open boundaries (standard and SGS) are normal-depth free
outflows. Against the 0.2.x default (`--scheme v0.2`), four options change:

| Option | 0.3.0 (`recommended`) | 0.2.x (`--scheme v0.2`) |
|---|---|---|
| Wet/dry limiter | positivity (one scale per donor cell), also at open boundaries (`--limiter positivity --boundary-limiter positivity`) | hybrid; legacy at boundaries |
| Momentum smoothing | none (`--q-centre-theta 1`) | consistent Q-centring, θ = 0.9 |
| Damping | η filter, C_s = 0.05 (`--eta-filter 0.05`, §15.8) | none |
| Open boundaries (also SGS) | normal-depth free outflow (`--open-boundary normal`, below) | the pre-0.3.0 extrapolation |

The diamond-wetdry face gradient, vector friction and the local CFL are
unchanged; so is SGS apart from its open boundaries (no reference exists yet
with a sub-grid bed).
`--momentum-damping` (§15.9) stays available and off. `--scheme v0.2`
reproduces 0.2.x results exactly (Carlisle, 4 h, `--cell-order file`, every
frame, A5 and Cartesian).

**Open boundaries.** A ZeroGradient (open) face extrapolates a ghost water
level from the boundary cell and an interior neighbour. Before 0.3.0 the
neighbour chosen was the most *outward* one — a comparison with the wrong
sign — which is one running along the boundary, so the slope ran along the
edge: zero at a uniform outlet (no head; water held back in the edge cells)
and, on a narrow domain, a spurious outflow through the side edges. With
the inward neighbour, extrapolating the *water* slope still has a flat pond
at the edge as a fixed point (on the 4 km rained plane the Cartesian outlet
ponded to 0.42 m). `--open-boundary normal` therefore extrapolates the
*bed* slope from the inward neighbour at the boundary cell's own depth
(zero gradient of depth: a free outflow at normal depth, as the SWE
transmissive boundary). On the rained 0.1% plane the last 250 m before the
outlet go from 0.22 m (A5) and 0.13 m (Cartesian) to 0.099 and 0.095 m;
the kinematic depth there is 0.097 m (`test/test_open_boundary.jl` T-BC11,
T-BC12). `--open-boundary legacy` keeps the old behaviour (`--scheme v0.2`,
`legacy`).

**Why.** Against the full shallow-water reference (§15.10, corrected solver;
the 0.3.0 default including the open boundaries):
- *Same mesh, Carlisle 6 h:* the 0.3.0 scheme stores 0.99× the reference on
  Cartesian (2.297 against 2.331 Mm³) and 1.25× on A5 (3.363 against
  2.690); the 0.2 scheme 1.62× and 1.80×, with deep channels moving at
  0.40–0.50 of the reference's speed. Adding momentum damping only
  increases storage.
- *Resolution (Cartesian 40/20/10 m):* 3.265 / 2.297 / 2.035 Mm³ against the
  reference's 3.154 / 2.331 / 2.107 — within 4% at every resolution and
  converging with it. The 0.2 scheme does not converge (3.80 → 3.83 Mm³
  from 20 to 10 m), because its Q-centring resistance grows as Δt shrinks
  with the cells.
- *SGS:* the open-boundary fix changes Carlisle 6 h storage by −0.8%
  (4.334 → 4.299 Mm³).
- *Time step:* halving the Courant number changes 6 h storage by 0.2% (A5)
  and 4% (Cartesian) with the η filter (§15.8), and by 10% and 16% with the
  0.2 scheme.
- *Planar tests* (0.1% slope, 20 mm/h rain or a 10 m³/s point source, 6 h,
  A5 res-18 and an equal-area Cartesian grid; `test/wp1_runs/v030/`): in
  the interior both schemes are the same and match the kinematic steady
  depth within 1–4 mm. The point-source spread agrees with the reference
  (CSI 0.89 on A5, 0.96 on Cartesian at 0.02 m; rms 9 and 4 mm). Next to
  the open outlet edge the 0.2 scheme ponds (0.34 and 0.26 m in the last
  250 m; the 0.3.0 default 0.099 and 0.095 m; kinematic 0.097 m).

**Trade-off.** Without Q-centring the standard solver keeps more grid-scale
velocity noise in deep channels: dv_rms 0.65 m/s on A5 Carlisle against
0.24 m/s for the 0.2 scheme (§15.9). The water surface is as smooth (the η
filter, §15.8). Where velocities matter more than conveyance,
`--momentum-damping 0.05` removes most of the noise (0.30 m/s) for 6–10%
more storage.

### 15.12 Channel representation: A5 against Cartesian (0.3.2, 2026-11)

*Question.* Does the A5 mesh represent a narrow channel at an arbitrary
orientation better than a Cartesian grid of the same cell area, which
staircases a diagonal channel?

*Series* (`test/channel_bench/`, `docs/TEST_CASES.md` §7):
- **Geometry.** A straight rectangular channel, 20 m wide with 2 m banks, cut into a 2 km
  valley whose floodplain rises across the valley at 0.5%. The bed slope is 0.1% along the channel, and the
  bed steepens to 2% at s = 0.3 L (a hydraulic control upstream of the
  outflow edges). Q = 20 m³/s enters at s = −0.45 L.
- **Mesh.** The cell-mean bed is sampled
  analytically by Halton points, as `sample_dem_mean!` does for a raster.
- **Cases.** 7 orientations (0–90°) × 2 cross-channel offsets (0 and ½ cell) ×
  A5 res 17–20 (cells 44, 22, 11, 5.6 m) and equal-area Cartesian grids.
- **Solvers.** Full SWE (`--flow-model swe`) and the 0.3.0 local-inertial
  default. Each run lasts 3 h; the res-20 runs were on the GPU.
- **Metric.** Δη is the mean
  water-surface error over the channel cells of the reach s ∈ [−0.3 L,
  0.05 L], against the 2-D normal flow (h_n = (Qn/(W√S))^0.6 = 0.97 m);
  lost conveyance shows as Δη > 0. The runs are steady, or close to it:
  outflow/inflow is 0.89–1.01 at 44 and 22 m and within 0.5% of 1 from
  11 m, and the stage changes by less than 4 cm between the last two frames,
  except up to 0.14 m in some 22 m SWE runs.

*Results* (mean Δη over the two offsets; figures `docs/figures/channel_vs_resolution.png`
and `channel_vs_angle.png`). Columns 0°–90° are the channel orientation (0° = channel along the x-axis); "spread" is max − min over orientations.

| solver | cell | mesh | 0° | 15° | 30° | 45° | 60° | 75° | 90° | rms | spread |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SWE | 11 m | A5 | 0.71 | 0.71 | 0.57 | 0.49 | 0.25 | 0.07 | 0.13 | 0.49 | 0.64 |
| SWE | 11 m | Cartesian | 0.18 | 0.30 | 0.43 | 0.53 | 0.43 | 0.28 | 0.17 | 0.35 | 0.36 |
| SWE | 5.6 m | A5 | 0.29 | 0.17 | 0.17 | 0.15 | 0.00 | −0.05 | −0.01 | 0.15 | 0.34 |
| SWE | 5.6 m | Cartesian | 0.03 | 0.05 | 0.11 | 0.12 | 0.12 | 0.08 | 0.03 | 0.09 | 0.09 |
| LI | 11 m | A5 | 0.50 | 0.51 | 0.45 | 0.33 | 0.10 | 0.06 | −0.06 | 0.35 | 0.57 |
| LI | 11 m | Cartesian | 0.23 | 0.12 | 0.29 | 0.48 | 0.30 | 0.27 | 0.22 | 0.29 | 0.36 |
| LI | 5.6 m | A5 | −0.08 | −0.04 | −0.02 | 0.11 | 0.00 | −0.01 | −0.06 | 0.06 | 0.19 |
| LI | 5.6 m | Cartesian | 0.09 | −0.05 | 0.00 | 0.16 | 0.00 | 0.10 | 0.09 | 0.09 | 0.21 |

- *Resolution.* Both meshes converge in the same way: Δη ≈ 1.1 m at 44 m
  cells (the channel is not resolved) and 0.8–1.0 m at 22 m, falling to
  0.3–0.5 m at 11 m and ≤ 0.15 m at 5.6 m.
- *Orientation.*
  - Cartesian behaves as expected: worst at 45°, best along the axes.
  - A5 is **not** orientation-free. It is worst at 0–15° and best at 60–90°,
    with a spread over orientations as large as Cartesian's or larger.
  - At 45° the two are equal (SWE, 11 m: 0.49 against 0.53 m).
  - With SWE, A5 has the larger rms error at every resolution from 22 m.
- *Why A5 has preferred directions.* Over a domain this size, A5 is a
  periodic tiling: every face normal lies in one of six directions (2°,
  41°, 103°, 146°, 155°, 173° from east, mod 180°; the same at res 17 and
  19). The largest gap between them is 62° (41°–103°), against 90° on a
  Cartesian grid. Face alignment alone does not explain the pattern,
  though: A5 does worst along 0–15°, where it has faces at 2° and 173°.
  The mechanism is open.
- *Local-inertial against SWE.*
  - On these in-bank cases the 0.3.0 default is
    as close to the analytic stage as SWE, or closer. At 5.6 m on A5 its rms
    is 0.06 m against SWE's 0.15 m, which keeps 0.29 m at 0°.
  - At coarse resolution, local-inertial spills onto the floodplain where the
    channel is unresolved, much more on Cartesian (up to 62 m of wet
    floodplain per metre of reach at 44 m) than on A5 (14 m). SWE barely
    spills (≤ 6 m).
  - `--momentum-damping 0.05` changes the median Δη by less than 0.05 m,
    except Cartesian at 11 m (0.29 → 0.14 m).

*Conclusion.* For a straight channel 0.5–4 cells wide, A5 gives no
systematic advantage over a Cartesian grid of the same cell area. It
trades Cartesian's 45° weakness for its own preferred directions; any
benefit of A5 for diagonal flow would have to show elsewhere than in
channel conveyance at a given resolution. This series uses one location
(51° N, 0° E); A5's local face directions vary between dodecahedron faces,
so the orientation pattern will differ elsewhere.

### 15.13 SGS conveyance (0.7.0, 2026-11)

*Question.* Until 0.6.0 SGS stored 2.5× the SWE reference at the Carlisle peak
(26.7 Mm³), with its water surface metres higher. Measurement or conveyance?
(Plan: `projectinstructions/FloodA5_SGS_Conveyance_Fix_Plan_2026-11.md`.)

*Diagnosis* (case `sgs_strip`, `docs/TEST_CASES.md` §2c,
`notebooks/04_sgs_strip.ipynb`). Storage was measured correctly; conveyance
was not. Edge curves were tabulated on each cell's own levels and stopped
growing at the cell's highest DEM sample, and the kernel averaged the two
cells' A and R: overbank a 20 m strip's faces carried 39% (channel) and 2.5%
(floodplain) of their true flow area. In bank the legacy donor cap, V/(2N)
per face per step, held channel cells brim-full.

*Changes* (§9): one table per face on the face's own levels, continued above
the profile; the 2-D conveyance form K = Σ h^(5/3); positivity limiting; a
time step from the hydraulic depth and from the faces' flow areas (cells
holding a sliver of a channel); flow-path geometry for partly wet cells;
connectivity (activation levels); tables from a per-cell working grid; DEM
profiles for boundary sides and a valley-floor outflow slope; Courant 0.5.

*Results.* Strip (analytic 2-D level −0.267 / +0.372 m): Cartesian 20 m
−0.279 / +0.362 m, A5 −0.277 / +0.340 m (2-D for the A5 strip's own width
+0.290; local-inertial on the same mesh +0.390); 0.6.0: −0.001 / +4.75 m and
still filling. Without the flow-path geometry A5 was 0.054 m high in bank;
without the face-area time-step limit, cells holding a sliver of the channel
oscillated and held it up by up to 0.5 m. A wall narrower than a cell holds
(case `planar_wall`; 0.6.0 passed 58–64% of the inflow).

*Carlisle at the peak* (`carlisle_peak`, `carlisle_peak_fine`). SGS is unstable
at Courant 0.7 in deep, fast channel flow (cell-scale level oscillations);
0.5 and 0.35 agree (5.94 / 5.96 Mm³). Against 5 m Cartesian references, with
the 5 m results aggregated onto the 22 m cells (water level where at least
half the cell is wet), at 55.5 h on A5:

| 22 m model | vs 5 m SWE: storage | level, median (rms) | vs 5 m LI: storage | level, median (rms) |
|---|---|---|---|---|
| SWE | 1.24 | +0.36 (0.54) m | 1.47 | +0.78 (0.81) m |
| local-inertial | 0.94 | −0.16 (0.25) m | 1.12 | +0.25 (0.35) m |
| SGS | 0.70 | −0.50 (0.54) m | 0.83 | −0.14 (0.24) m |

The 22 m SWE run is not a reference (24% more storage than 5 m SWE). SGS
reproduces the fine-grid local-inertial solution best, as a sub-grid scheme
for these equations should; but the local-inertial equations convey more than
SWE at 5 m (0.84 of its storage, levels 0.3–0.4 m lower) — the momentum losses
they omit — and 22 m local-inertial lands close to 5 m SWE only because its
smoothed bed under-conveys. Open: how to judge the 22 m models against the
5 m runs, and whether SGS needs those momentum losses.

### 15.14 DEM registration and the SGS critical-flow limit (0.8.0, 2026-11)

*DEM registration.* GDAL rasters are pixel-is-area: the geotransform gives
the corner of the top-left pixel, and a pixel's value belongs at its centre.
`_affine_to_pixel` returns corner-based pixel coordinates, which `_bilinear`
used as centre-based ones, so every sample until 0.7.1 was taken half a pixel
east and south of its point (on a 5 m DEM 2.5 m each way; on a slope of
1:100 that is 3.5 cm of bed). Now `_bilinear` interpolates between pixel
centres, holds the edge value within half a pixel of the raster's edge, and
falls back to the containing pixel next to nodata (`test/test_dem_registration.jl`).
Meshes record `dem_registration` (absent: sampled before 0.8.0; a warning on
`--meshload`) and SGS tables are layout 4. Raster-aligned Cartesian meshes
(`generate_cartesian_mesh_raster`, the Carlisle 20 m and 5 m grids) take the
raster's pixels and are unchanged: their runs are bit-identical. The A5 runs
move by round-off-sized amounts on real terrain: Carlisle 6 h SWE storage
2.690 → 2.641 Mm³ (−1.8%), at the peak +1.7%; local-inertial −0.1%; the
thin `planar_wall` ridge seen by the cell-mean beds 5.361 → 5.387 m.

*Critical-flow limit* (§9.2). The local-inertial equations have no advection,
so nothing stops a face passing more than critical flow where the water falls
away from it. SGS faces are now capped at √g·(2/3)^(3/2)·C(E), with
C = Σ dx·(E − a_p)^(3/2) the face's critical integral at the upstream level
(`--sgs-critical on`, default). On `planar_embankment` the notch (4 m wide,
sill 5.0 m) passed ≈ 9.8 m³/s at 1.07 m of head in 0.7.x, 30% above
critical; now 8.24 m³/s at 1.38 m between 5 and 6 h (weir 11.07): the pond
stands 0.33 m higher (6.41 m) and 26% less water is downstream at 6 h
(115 000 m³, was 155 000). On `sgs_strip` the levels in the measuring reach
are unchanged to the millimetre, but the cap binds on the 2.1% outlet reach
(supercritical; the cap neglects the velocity head): SGS holds 1.2–2.0% more
water over the whole strip, all of it below x = 3000 m (Cartesian in bank
13 361 against 13 103 m³; with `--sgs-critical off` 13 103 exactly, so the
DEM fix changes nothing there). Including the velocity head in E would
remove this; whether to do so belongs with the hydraulic formulations.

*Carlisle against the 5 m SWE benchmark* (decided 2026-11-23; case
`carlisle_peak_fine` now compares every 22 m run of `carlisle_peak` with both
5 m runs, aggregated onto the 22 m cells). At 55.5 h, A5 res 18 unless noted:

| 22 m model | vs 5 m SWE: storage | level median (rms) | CSI | vs 5 m LI: storage | level median (rms) | CSI |
|---|---|---|---|---|---|---|
| SWE | 1.26 | +0.39 (0.56) m | 0.81 | 1.49 | +0.82 (0.84) m | 0.78 |
| SWE, Cartesian 20 m | 1.12 | +0.15 (0.29) m | 0.85 | 1.33 | +0.55 (0.58) m | 0.82 |
| local-inertial | 0.94 | −0.18 (0.25) m | 0.85 | 1.12 | +0.24 (0.34) m | 0.84 |
| local-inertial, Cartesian 20 m | 0.95 | −0.22 (0.29) m | 0.85 | 1.12 | +0.20 (0.32) m | 0.85 |
| SGS, critical limit off (0.7.x hydraulics) | 0.71 | −0.49 (0.53) m | 0.77 | 0.84 | −0.14 (0.23) m | 0.81 |
| **SGS (0.8.0)** | **0.74** | **−0.47 (0.52) m** | **0.81** | **0.87** | **−0.12 (0.16) m** | **0.85** |

(5 m local-inertial stores 0.84 of 5 m SWE: 7.07 against 8.38 Mm³.) The limit
improves SGS against both 5 m runs — storage +3–4%, CSI +0.04, rms against
5 m local-inertial 0.23 → 0.16 m (against 5 m SWE within ±0.02 m either
way) — and stays on. SGS is the 22 m model
closest to the 5 m local-inertial solution in level and extent at all three
times (storage: 0.87 against 1.12–1.16 for cell-mean local-inertial; at 55.5 h
marginally the further of the two). Against 5 m SWE it fails the agreed bar
(`carlisle_peak_fine`, an expected failure): the gap is the hydraulics, not
the geometry — at 5 m the local-inertial equations themselves store 16% less
than SWE, and the cell-mean 22 m run matches SWE's storage only because its
smoothed beds store ~15% more than the 5 m local-inertial run. A sub-grid
geometry with local-inertial hydraulics cannot close that gap; a formulation
that keeps the momentum losses the local-inertial equations omit can (plan:
pluggable hydraulic formulations on the SGS geometry,
`projectinstructions/PROJECT_STATE.md`). Open: where at 5 m SWE holds the
extra 1.3 Mm³ (bends, constrictions, the bridges' reaches), which says which
momentum terms matter.

### 15.15 The SGS checkerboard: sources and the adaptive time step (0.9.0, 2026-11)

*Symptom.* Since 0.7.0 the SGS water surface on Carlisle broke into
cell-scale noise on most mesh families before or at the Eden's peak
(notebook 06 §8): at 55.5 h 36–42% of the fully wet cells were local level
extrema on Cartesian, rotated Cartesian and rHEALPix (0.1% for
local-inertial), maximum depths reached 11 m against 6.8 m, and the noise
started at times from 1.5 h (rHEALPix) to 50 h (hexagons); A5, ISEA3H and the
jittered Voronoi mesh stayed smooth. Every onset was within 20–40 m of an
inflow cell. Ruled out on the Cartesian mesh: the Courant number (0.35 gives
the same onset at a lower amplitude), the critical-flow limit, the flow-path
geometry, lumped conveyance, gradient correction, the η filter (standard
flow is smooth with and without it); component friction delays the onset on
Cartesian only. Q-centring at θ = 0.9 suppresses it but stores 50% more
water (rHEALPix, 4 h).

*Mechanism.* Sources were added to their cells before the face fluxes were
computed (`--source-order pre`). The Eden inflow (≈ 900 m³/s at the peak)
adds Q·Δt ≈ 1150 m³ a step to a 493 m² channel cell, 2.4 m of level, which
its four faces take back within the step: the level the faces see is the
stored level plus Q·Δt/a, so it moves with Δt by Q/a ≈ 1.85 m per second
of Δt. The adaptive Δt in turn follows the levels and face areas of these
cells (the SGS face rate, `_sgs_face_rate`). A step-by-step trace of the
inflow cell and its neighbours on Cartesian at 46.5 h shows an exact steady
state for hundreds of steps, then one change of Δt (1.28 → 0.83 s) that
starts an oscillation of the inflow cell against its neighbours with a
period of ≈ 7 steps (9 s) — the cell's gravity mode, ω = √(g·ΣA_f/L_f / a)
= 0.65 s⁻¹ — during which Δt oscillates with it (1.28, 1.24, 1.18, 1.11,
1.06 s …) and the amplitude grows from ±0.1 to ±0.9 m in 60 s. It spreads
downstream from there. Tests on rHEALPix (4 h, extrema at 4 h): a constant
Δt of 0.5 s or 0.7 s (below and above the adaptive steps) is smooth (0.003,
against 0.112), the adaptive step at Courant 0.25 is not (0.026 and
growing); spreading the inflow over the cell and its neighbours only delays
it (0.092). The variation of Δt drives it, not its size.

*Fix.* `--source-order post` (`_apply_sources_post!`; the default for SGS)
adds rainfall and point sources after the step's flux update, as part of the
same volume change, and brings the depth of the cells they touch up to date
(the device backends likewise). The faces then see the level after the whole
previous step, which at equilibrium does not depend on Δt; the inflow cell no
longer swings by Q·Δt/a inside each step. Test S10 (`test_sgs_general.jl`):
at steady state a 40 m³/s point source's cell holds 0.5637 / 0.5636 /
0.5636 m at Δt = 0.5 / 1 / 2 s, against 0.5551 / 0.5465 / 0.5294 m with
`pre` (Q·ΔΔt/a exactly). Standard flow keeps `pre` (its results are
unchanged; it was smooth on every family) and SWE has only `pre`. A limit on
the growth of Δt from step to step (tried as `--dt-growth`) removed most of
the noise too (0.016 on rHEALPix at 4 h), but each fall of Δt still kicks
the mode; it was not needed once the sources moved and was dropped.

*Carlisle, all eleven families at ~22 m* (`mesh_families_carlisle` jobs,
55.5 h; against 5 m SWE aggregated onto the cells, `FA.compare_to_fine`):

| mesh | extrema 0.8.x → now | r_rms (mm) | max depth (m) | storage vs 5 m SWE | level median (m) | CSI |
|---|---|---|---|---|---|---|
| Cartesian | 0.408 → 0.007 | 414 → 52 | 11.15 → 6.78 | 1.06 → 0.74 | +0.23 → −0.45 | 0.90 → 0.88 |
| rotated Cartesian | 0.424 → 0.016 | 527 → 38 | 10.81 → 6.73 | 0.96 → 0.75 | −0.07 → −0.49 | 0.88 → 0.82 |
| rHEALPix | 0.361 → 0.004 | 305 → 19 | 10.83 → 6.68 | 0.76 → 0.76 | −0.47 → −0.49 | 0.85 → 0.85 |
| hexagons | 0.083 → 0.005 | 127 → 47 | 7.24 → 6.81 | 0.77 → 0.79 | −0.46 → −0.43 | 0.87 → 0.88 |
| H3 | 0.183 → 0.004 | 285 → 20 | 7.70 → 6.67 | 0.71 → 0.72 | −0.49 → −0.50 | 0.82 → 0.82 |
| IGEO7 | 0.167 → 0.005 | 658 → 33 | 9.79 → 6.63 | 0.82 → 0.72 | −0.26 → −0.55 | 0.83 → 0.80 |
| HEALPix | 0.227 → 0.012 | 57 → 39 | 6.76 → 6.76 | 0.76 → 0.77 | −0.42 → −0.41 | 0.83 → 0.83 |
| Voronoi | 0.143 → 0.005 | 79 → 50 | 7.27 → 6.71 | 0.75 → 0.75 | −0.50 → −0.49 | 0.84 → 0.84 |
| jittered Voronoi | 0.008 → 0.008 | 65 → 48 | 6.73 → 6.73 | 0.76 → 0.76 | −0.49 → −0.49 | 0.84 → 0.84 |
| ISEA3H | 0.006 → 0.006 | 42 → 41 | 6.67 → 6.67 | 0.70 → 0.70 | −0.54 → −0.54 | 0.80 → 0.80 |
| A5 res 18 (`carlisle_peak`) | 0.009 → 0.005 | 42 → 25 | 6.68 → 6.68 | 0.735 → 0.735 | −0.47 → −0.47 | 0.81 → 0.81 |

(Local-inertial on the 22.2 m Cartesian mesh: extrema 0.001, r_rms 60 mm,
storage 0.94, median −0.08 m, CSI 0.85.) Every family is now smooth through
the peak, and the families that were smooth before are unchanged to the
third decimal. The noise had raised the levels where it was strong
(Cartesian, rotated Cartesian, IGEO7): with it gone, SGS stores 0.70–0.79 of
5 m SWE on every family with levels 0.4–0.55 m low — the hydraulic gap of
§15.14, now the same on every mesh. Rerun of notebook 06 §8 (`RERUN_SGS =
true`) gives the per-mesh picture.

### 15.16 What separates SGS from the 5 m runs: geometry, the grid, the equations (0.9.2, 2026-11)

*Question.* After 0.9.0 SGS stores 0.70–0.79 of the 5 m SWE benchmark at the
Carlisle peak on every mesh family, with levels 0.4–0.55 m low. How much of
that is the volume representation (a cell's hypsometric table against the
pixels it stands for), how much the 5 m grid, and how much the equations?

*Exact tables* (`--sgs-sampling pixels`, `build_sgs_tables!(…; sampling =
:pixels)`). On a Cartesian mesh whose cells are whole k × k blocks of the
DEM's pixels (`cases/scripts/block_dem.jl`: block means from the raster's own
origin), the tables are built from the pixels: a cell's volume at level η is
exactly Σ_p (A/k²)·max(η − z_p, 0) (bins at every pixel elevation), and a
face's bed per pixel pair is max(z_in, z_out) — what a model on the pixels
has. An SGS model on the blocks then stores, at any level, what the pixels
store (`test/test_sgs_pixels.jl`).

*Volume.* A pool at rest on rough ground (case `sgs_pool`, TEST_CASES §2e):
SGS with pixel tables reaches the exact level to 0.1 mm, the default
(bilinear) tables 2.3 mm high (0.4% of the volume: a bilinear surface is
smoother than the pixels), cell-mean beds 11.6 mm. At Carlisle, on a 25 m
grid of exact 5 × 5 blocks of the 5 m LiDAR, a flat level in each block at the
5 m runs' mean wet level holds within 0.1% of their volume. The volume
representation is not the gap.

*Carlisle, aligned 25 m grid* (55.5 h; per block, the level at which its pixel
table holds the 5 m volume):

| 25 m model | vs 5 m standard: storage | level median (rms) | vs 5 m SWE: storage | level median (rms) |
|---|---|---|---|---|
| SGS, pixel tables | 0.892 | −0.07 (0.41) m | 0.752 | −0.46 (0.67) m |
| SGS, bilinear tables | 0.899 | −0.11 (0.28) m | 0.758 | −0.49 (0.60) m |
| SGS, pixel tables, connectivity off | 0.930 | −0.07 (0.30) m | 0.784 | −0.47 (0.49) m |
| standard flow, block means | 1.091 | +0.18 (0.43) m | 0.920 | −0.22 (0.40) m |

5 m SWE holds 1.31 Mm³ more than 5 m standard flow: a backwater over the
Eden's upstream floodplain (+0.4–0.5 m, almost uniform) stepping down at the
constricted, fast reach through the city (> 3 m/s). SGS's 0.76 Mm³ less than
5 m standard flow is two-thirds in slow areas south of the Eden: in the 5 m
runs the Petteril overtops its embanked reach and floods the city and the
defended compartments behind it, which SGS keeps in bank. Turning
connectivity off fills some of those compartments, but by letting water
under walls whose crests the 5 m runs never reach from that side — the leak
connectivity exists to stop (`planar_wall`), not a fix.

*The grid.* A straight channel at an angle to the grid (case
`oblique_channel`, TEST_CASES §2f): aligned, every model holds the analytic
normal level; at 26.565° and 45° standard flow on 5 m pixels is +0.12–0.14 m
(+0.06 m on 2.5 m pixels), SGS at 25 m +0.06–0.08 m (−0.01 m from 2.5 m
pixels), and 5 m SWE +0.22–0.85 m for a 15 m channel and +0.10–0.34 m for a
45 m one, barely less on 2.5 m pixels: the staircase of the rasterised banks
costs SWE head at every step. Carlisle on a 5 m grid turned 26.565°
(`cases/scripts/rotate_dem.jl`, an oblique Mercator on the DEM's datum,
bilinear resampling; nodata cells dropped, `drop_nodata_cells.jl`):

| 55.5 h, Mm³ | 5 m standard | 5 m SWE | SWE ÷ standard | 25 m SGS (pixel tables) |
|---|---|---|---|---|
| original grid | 7.066 | 8.380 | 1.186 | 6.302 |
| grid turned 26.565° | 6.600 | 8.253 | 1.250 | 6.425 |

The 5 m standard run moves 6.6% with the grid's alignment — the Petteril's
overtopping is a threshold, the kind of grid-alignment effect Nguyen et al.
(2025, doi:10.1029/2024WR038919) measured with LISFLOOD-FP (extent ±7%,
rotation worse than translation, more where water just overtops) — and SGS
2%. SWE's backwater over the upstream floodplain is on both grids.

*The equations.* An abrupt contraction and expansion in a channel aligned
with the grid (case `channel_constriction`, TEST_CASES §2g), against 1-D
answers started from each run's own level below the opening: standard flow
and SGS follow the local-inertial answer (no convective term: no head lost at
the opening) within 2 cm; 5 m SWE stands 0.24 m above the Borda–Carnot answer
and 0.03 m below the step-backwater coefficients for abrupt transitions
(0.6 / 0.8). The head SWE loses at a constriction is physical.

*Conclusions.* (1) SGS's geometry is right: exact volumes, and it matches 5 m
standard flow where the physics is the same (`channel_constriction`,
`oblique_channel`). (2) Against 5 m standard flow — the like-for-like
reference for SGS as it is, local-inertial on its tables — SGS on A5 is
closer than cell-mean local-inertial in level (rms 0.16 against 0.34 m) and
extent at every time, and further in storage only at 55.5 h, by 0.01; the
rest of its gap there is a threshold (the Petteril) that moves with the 5 m
grid's alignment (`carlisle_peak_fine`, second bar, an expected failure).
(3) Against 5 m SWE the remaining gap is the equations: the head lost at the
Eden's constricted reach, which local-inertial hydraulics cannot lose. Part
of the 5 m SWE run's excess may also be its staircase banks wherever the
channels run obliquely; the backwater survives turning the grid, but every
orientation has oblique reaches. (4) So the hydraulic formulation for SGS
should carry convective momentum, or at least form losses where a face's flow
area changes sharply — judged on `channel_constriction` (and its oblique
variants) before Carlisle.
