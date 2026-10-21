# surfacewater/flow2d.jl
# ----------------------
# Pure physics kernels for 2D diffusive-wave surface water modelling.
#
# This file contains all functions and constants that depend only on scalar
# hydraulic arguments — no FlowState, no mesh, no I/O.  The intent is:
#
#   1. These kernels are included by FloodModel.jl after the type definitions.
#   2. They can also be included independently by external code (e.g. a 3DGeo
#      plugin or a GPU kernel wrapper) without pulling in the full application.
#   3. When the FlowState type boundary is formalised, this file becomes the
#      body of a `module Flow2D` with a clean export list.
#
# GPU note: all @inline functions here take only Float64 scalars and have no
# heap allocations — they map directly to CUDA device functions.  The step
# functions (step_standard!, step_sgs!) in FloodModel.jl that call these in
# a loop become the CUDA kernel launch wrappers.
#
# Process additions: for infiltration, evaporation, or other surface-water
# processes, add a companion file (e.g. surfacewater/infiltration.jl) that
# follows the same pattern: pure scalar kernels, no FlowState dependency,
# included by FloodModel.jl.
#
# ENV["FLOODMODEL_INCLUDE_ONLY"] note: this file contains only definitions
# and constants — no top-level executable code.  It is safe to include at
# any point without triggering side effects.

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

const _G = 9.81   # m s⁻²

# Minimum h_flow (m) below which an edge is treated as dry and flux set to zero.
# Matches the CAESAR-Lisflood default hflow_threshold = 0.001 m.  Skipping
# near-dry edges avoids carrying stale q_prev momentum on thin films and gives
# a meaningful computational saving on large meshes (most boundary edges are dry).
const HFLOW_THRESHOLD = 0.001

# Froude number limit for subcritical flow.  Matches CAESAR-Lisflood froude_limit = 0.8.
# Unit discharge is capped at q_max = h_flow × √(g × h_flow) × FROUDE_LIMIT.
# This suppresses the supercritical oscillation mode that drives checkerboarding
# in inertial models on irregular meshes.
const FROUDE_LIMIT = 0.8

# Q-centred (spatial momentum smoothing) parameter — θ in the expression:
#   q_eff = θ × q_prev_e + (1−θ)/2 × (q_prev_collinear_i + q_prev_collinear_j)
# θ = 1.0 disables smoothing (pure Bates). θ = 0.9 is the LISFLOOD-FP default.
# Lower values provide stronger damping of the checkerboard mode but reduce
# the inertial accuracy of the momentum term. Applied in step_standard! and
# step_sgs! Phase A via _q_centred(). Decoupled from FROUDE_LIMIT so both
# fixes can be tuned independently.
#
# Promoted from a module-level constant to a _q_centred argument
# (flow-direction-fixes, Step 2) so it can be set at runtime via
# --q-centre-theta and varied for A/B testing against the WLSQ gradient
# correction (FloodA5_NonOrthogonal_Correction_Plan.md §6). The value lives
# on FlowState.q_centre_theta; this default is used only by any legacy
# caller that does not pass theta explicitly.
const Q_CENTRE_THETA_DEFAULT = 0.9

"""
    _q_centred(flux, e, col_i, col_j, theta=Q_CENTRE_THETA_DEFAULT) → Float64

Return the Q-centred (spatially smoothed) unit discharge for edge `e`.

  q_eff = θ × flux[e] + (1−θ)/2 × (flux[col_i] + flux[col_j])

where `col_i` and `col_j` are the indices of the most-collinear edges on
the ci and cj sides respectively (stored in EdgeList.collinear_i/j).
0 = absent (boundary cell); handled by normalising over available neighbours.

A negative index (WP2c, 2026-10) means the partner edge is oriented against
`e` — its stored discharge has the opposite sign for the same physical
flow — so its flux enters with its sign flipped. `_orient_collinear!` sets
the signs; without it (`q_centre_mode = :legacy`) every index is
positive, which is correct on row-major Cartesian grids but on A5 reverses
about two partners in three (edges run from the lower to the higher cell
index, and A5 cell order is not geometric).

The checkerboard mode has opposite signs on alternating edges, so the average
over collinear neighbours naturally cancels it. For a coherent flow field
(all fluxes same sign and magnitude) the smoothing has negligible effect.

`theta` is normally `state.q_centre_theta` (set via --q-centre-theta, default
0.9). Passing `theta = 1.0` disables smoothing entirely (q_eff = flux[e]).

Thread safety (WP0, 2026-10): `flux` must be a READ-ONLY snapshot of the
previous step's values while the threaded edge loop runs. Passing the live
`edges.flux` array (which the same loop overwrites) is a data race: the
result then depends on edge order and thread scheduling. `step_standard!`
and `step_sgs!` pass a start-of-step copy.
"""
@inline function _q_centred(flux  :: Vector{Float64},
                              e    :: Int,
                              ci   :: Int,   # collinear_i[e]
                              cj   :: Int,   # collinear_j[e]
                              theta:: Float64 = Q_CENTRE_THETA_DEFAULT
                              )::Float64
    q0   = flux[e]
    n_nb = (ci != 0 ? 1 : 0) + (cj != 0 ? 1 : 0)
    n_nb == 0 && return q0
    q_nb = (ci != 0 ? flipsign(flux[abs(ci)], ci) : 0.0) +
           (cj != 0 ? flipsign(flux[abs(cj)], cj) : 0.0)
    return theta * q0 + (1.0 - theta) * q_nb / n_nb
end

# ---------------------------------------------------------------------------
# Friction model (WP1, 2026-10): component vs vector
# ---------------------------------------------------------------------------
# The Bates (2010) / LISFLOOD-FP "ACC" friction term uses only the face-normal
# discharge |q_n| in the denominator ("component" or decoupled friction). In
# the friction-dominated limit this makes the face conductance depend on face
# orientation: on a Cartesian grid diagonal flow is 2^¼ too conductive and the
# flow is deflected up to 12.8° from the gradient; on A5 the deflection is
# systematic and does not cancel north/south. Refinement cannot remove it
# (it is a property of the scheme's continuum limit). See
# FloodA5_Review_and_ForwardPlan_2026-10.md §3.1–§3.2 and Cozzolino et al.
# (2021, Adv. Water Resour. 152:103919).
#
# "Vector" friction uses the full discharge magnitude |q| = √(q_n² + q_t²),
# where q_t is the face-tangential discharge from the reconstructed cell
# vectors (Phase F), the mesh-general form of LISFLOOD-FP's "qvect" option.
# The same magnitude is used by the Froude limiter, which then scales q_n by
# min(1, q_max/|q|) instead of clamping the normal component alone.
#
# Selected per run via --friction-model component|vector (FlowState.
# friction_model). The kernels take `q_t` and `vector_friction` as keyword
# arguments; with the defaults (q_t = 0, vector_friction = false) every kernel
# is bit-identical to its pre-WP1 behaviour.

"""
    _friction_q(q_n, q_t, vector_friction) → Float64

Discharge magnitude used in the friction denominator: `|q_n|` (component,
Bates 2010) or `√(q_n² + q_t²)` (vector). Works for unit discharge (m²/s) and
volumetric discharge (m³/s) alike.
"""
@inline _friction_q(q_n::Float64, q_t::Float64, vector_friction::Bool) =
    vector_friction ? sqrt(q_n * q_n + q_t * q_t) : abs(q_n)

"""
    _froude_cap(q_new, q_t, q_max, vector_friction) → Float64

Froude limiter (Fix A). Component: clamp the normal discharge to ±q_max.
Vector: scale the normal discharge by min(1, q_max/|q|), |q| = √(q_new² + q_t²),
so the cap acts on the discharge vector's magnitude rather than on one
face-normal component (which is anisotropic in the same way as component
friction).
"""
@inline function _froude_cap(q_new::Float64, q_t::Float64, q_max::Float64,
                             vector_friction::Bool)::Float64
    vector_friction || return clamp(q_new, -q_max, q_max)
    q_mag = sqrt(q_new * q_new + q_t * q_t)
    return q_mag > q_max ? q_new * (q_max / q_mag) : q_new
end

"""
    _bates_flux(q_prev, wse_i, wse_j, z_sill, width, L, cos_theta, n_mann, dt) → Float64

Compute the volumetric flux Q (m³/s) from cell i to cell j using the
inertial shallow-water formulation of Bates, Horritt & Fewtrell (2010),
equation (9):

    q^t = [ q^{t-dt} - g · h_flow · dt · (dWSE / L) ]
          / [ 1 + g · h_flow · dt · n² · |q^{t-dt}| / h_flow^(10/3) ]

    Q^t = q^t · width

where:
  q^{t-dt}   unit discharge from previous timestep (m²/s), signed i→j
  h_flow     flow depth at the edge = max(WSE_i, WSE_j) - z_sill  (m)
  dWSE       WSE_i - WSE_j  (m), positive when flow is i→j
  L          centre-to-centre great-circle distance (m)
  cos_theta  cosine of angle between face normal and centre-to-centre vector
  width      shared edge length (m)
  n_mann     Manning's roughness (s m^{-1/3})

Returns volumetric flux Q (m³/s), positive = flow from i to j.
Returns 0.0 if h_flow ≤ HFLOW_THRESHOLD (no water above sill on either side).

Reference
---------
Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). A simple inertial
formulation of the shallow water equations for efficient two-dimensional
flood inundation modelling. Journal of Hydrology 387(1–2), 33–45.
https://doi.org/10.1016/j.jhydrol.2010.03.027

Sign convention (EdgeList)
---------------------------
Arguments are passed as (wse_ci, wse_cj) where ci = EdgeList.cell_i,
cj = EdgeList.cell_j, with cell_i < cell_j (canonical lower-index first).
  Q > 0  →  flow from cell_j to cell_i  (j is higher, water flows toward i)
  Q < 0  →  flow from cell_i to cell_j  (i is higher, water flows toward j)

Callers apply:  dV[ci] += Q*dt  and  dV[cj] -= Q*dt

Non-orthogonality correction
-----------------------------
L_eff = L × cos θ projects the centre-to-centre distance onto the face
normal.  cos θ is stored in EdgeList.cos_theta (1.0 = orthogonal).

⚠️  UNCORRECTED FOR DIRECTIONAL BIAS on non-orthogonal meshes.
    This L_eff = L × cos θ projection corrects gradient *magnitude* but
    not gradient *direction* — it implicitly assumes the WSE gradient
    points along the centre-to-centre vector, when on a non-orthogonal
    A5 mesh the true gradient points along the face normal. This produces
    the systematic flow-direction errors documented in the FOSS4G 2026
    paper (point-spread test: NW–SE axis bias; planar slope test: 30–60°
    deviation from downslope). Use `_bates_flux_corrected` (WLSQ gradient
    reconstruction + skewness correction) for production simulations —
    see FloodA5_NonOrthogonal_Correction_Plan.md. This function is
    retained for benchmarking and A/B comparison only
    (`--gradient-correction off`).
"""
@inline function _bates_flux(q_prev     :: Float64,
                              wse_i      :: Float64,
                              wse_j      :: Float64,
                              z_sill     :: Float64,
                              width      :: Float64,
                              L          :: Float64,
                              cos_theta  :: Float64,
                              n_mann     :: Float64,
                              dt         :: Float64;
                              q_t             :: Float64 = 0.0,
                              vector_friction :: Bool    = false)::Float64
    h_flow = max(wse_i, wse_j) - z_sill
    h_flow <= HFLOW_THRESHOLD && return 0.0
    h_flow  = max(h_flow, 1e-6)   # floor to avoid h^(10/3) underflow → NaN denominator

    dWSE  = wse_i - wse_j
    ct    = max(cos_theta, 0.1)   # hard floor: cos(84°) ≈ 0.10
    L_eff = max(L * ct, 1.0)      # guard degenerate zero-distance pairs

    # Bates et al. (2010) eq. 9 — unit discharge q (m²/s)
    numerator   = q_prev - _G * h_flow * dt * dWSE / L_eff
    denominator = 1.0 + _G * h_flow * dt * n_mann^2 *
                  _friction_q(q_prev, q_t, vector_friction) / h_flow^(10.0/3.0)
    q_new = numerator / denominator

    return q_new * width
end


"""
    _bates_flux_limited(q_prev, wse_i, wse_j, z_sill, width, L, cos_theta,
                        n_mann, dt, depth_donor) → (Q, q_stored)

Variant of `_bates_flux` used exclusively by `step_standard!`.  Applies three
stability fixes that match the CAESAR-Lisflood qroute() implementation:

**Fix A — Froude limiter** (CAESAR froude_limit = 0.8)
  Caps |q| at `h_flow × √(g × h_flow) × FROUDE_LIMIT` after Bates eq. 9.
  Prevents supercritical discharge and suppresses the checkerboard instability
  on irregular/pentagonal meshes where each edge has an independent q_prev.

**Fix B — Volume limiter** (CAESAR: depth/4 threshold → depth/5 cap per edge
  on A5, where the mesh has 5 neighbours per interior cell)
  Caps |Q × dt| at `depth_donor × width / (vol_limit_divisor × dt)`, so no
  more than `1/vol_limit_divisor` of the donor cell's water can leave via
  one edge per step. `vol_limit_divisor` defaults to `5.0` (A5's own
  neighbour count, i.e. today's exact original behaviour) but the caller
  should pass the mesh's own neighbour ceiling (`state.max_neighbours` —
  e.g. `4.0` for a Cartesian grid) so the ~20%-per-edge reasoning reflects
  this mesh's actual connectivity rather than an A5-specific assumption
  baked into the default. On any mesh, `width` is the natural spatial scale
  (shared edge length).

**Fix C — Consistent q_stored**
  Returns the *post-limiting* unit discharge as `q_stored` for writing back to
  `edges.flux[e]`.  This ensures the momentum state carried into the next step
  reflects what was actually transferred.  The previous approach stored the
  unlimited Bates q while only capping the volume in the scatter phase, causing
  divergence between stored momentum and actual hydraulic state — the root
  driver of step-to-step overshoot and checkerboarding.

Returns `(Q, q_stored)` where Q is the (limited) volumetric flux (m³/s) and
q_stored is the unit discharge (m²/s) to persist as q_prev next step.

⚠️  UNCORRECTED FOR DIRECTIONAL BIAS — see `_bates_flux`'s docstring. Use
    `_bates_flux_limited_corrected` for production simulations; this
    function is retained for `--gradient-correction off` benchmarking.
"""
@inline function _bates_flux_limited(q_prev       :: Float64,
                                      wse_i        :: Float64,
                                      wse_j        :: Float64,
                                      z_sill       :: Float64,
                                      width        :: Float64,
                                      L            :: Float64,
                                      cos_theta    :: Float64,
                                      n_mann       :: Float64,
                                      dt           :: Float64,
                                      depth_donor  :: Float64;
                                      vol_limit_divisor :: Float64 = 5.0,
                                      q_t               :: Float64 = 0.0,
                                      vector_friction   :: Bool    = false
                                      )::Tuple{Float64,Float64}
    h_flow = max(wse_i, wse_j) - z_sill
    if h_flow <= HFLOW_THRESHOLD
        # Edge is dry — zero out stale momentum so adjacent cells don't
        # carry phantom velocity in _compute_velocity!. Returning 0.0 for Q
        # without clearing edges.flux[e] was the source of spurious high
        # velocity on recently-dried cells. The caller stores q_stored in
        # edges.flux[e] after this call, so returning 0.0 here is sufficient.
        return (0.0, 0.0)
    end
    h_flow = max(h_flow, 1e-6)

    dWSE  = wse_i - wse_j
    ct    = max(cos_theta, 0.1)
    L_eff = max(L * ct, 1.0)

    # Bates et al. (2010) eq. 9 (component or vector friction — see _friction_q)
    numerator   = q_prev - _G * h_flow * dt * dWSE / L_eff
    denominator = 1.0 + _G * h_flow * dt * n_mann^2 *
                  _friction_q(q_prev, q_t, vector_friction) / h_flow^(10.0/3.0)
    q_new = numerator / denominator

    # Fix A: Froude limiter — cap at subcritical limit (Fr ≤ FROUDE_LIMIT)
    q_max = h_flow * sqrt(_G * h_flow) * FROUDE_LIMIT
    q_new = _froude_cap(q_new, q_t, q_max, vector_friction)

    # Fix B: Volume limiter — no more than 1/vol_limit_divisor of donor depth
    # per edge per step. vol_limit_divisor defaults to 5.0 (A5's neighbour
    # count, matching the original CAESAR-derived behaviour exactly) but the
    # caller should pass the mesh's own max_neighbours (e.g. 4.0 for a
    # Cartesian grid) so the cap reflects "no more than 1/n_sides of this
    # cell's donor depth", not an A5-specific assumption. Kept as a plain
    # scalar argument (not a FlowState/mesh reference) so this function
    # stays data-structure-independent — see the module-level docstring.
    if depth_donor > 0.0
        q_vol_max = (depth_donor * width) / (vol_limit_divisor * dt)
        q_new = clamp(q_new, -q_vol_max, q_vol_max)
    end

    # Fix C: return the post-limiting q so the caller stores a consistent q_prev
    Q = q_new * width
    return (Q, q_new)
end


# ---------------------------------------------------------------------------
# WLSQ non-orthogonal gradient correction — corrected flux kernels
# ---------------------------------------------------------------------------
# These supersede _bates_flux / _bates_flux_limited above for directional-
# bias-free flow on the A5 pentagonal mesh. See
# FloodA5_NonOrthogonal_Correction_Plan.md §2, §5.3, §10.6 for the full
# derivation (including the 2026-06-24 revision after two bugs were found
# via T-NOC5b — see that section for the complete account).
#
# Key difference from the legacy kernels: instead of scaling the raw WSE
# difference by L_eff = L × cos_theta (which corrects gradient *magnitude*
# but not *direction* — the root cause of the flow-direction bias documented
# in the FOSS4G 2026 paper), these kernels accept a pre-computed, fully
# corrected driving head `dWSE_n`, built by the caller from:
#   - the WLSQ-reconstructed cell-centre gradients (_build_wlsq_weights!,
#     _compute_wse_gradients!)
#   - the edge's ORIENTED cosine c and correction direction vector V̂
#     (EdgeList.cos_theta, EdgeList.skew_x/skew_y — see _edge_geometry's
#     docstring in A5Grid.jl for what these actually contain; the field
#     names skew_x/skew_y are historical and now hold a directional vector,
#     not a positional skewness offset)
#
#   ∇WSE_f = 0.5 × (grad_wse[:,ci] + grad_wse[:,cj])      (face-interpolated)
#   c      = edges.cos_theta[e]                            (oriented, ∈[0,1])
#   dWSE_n = c·(wse_i - wse_j) - L·(∇WSE_f · (skew_x[e], skew_y[e]))
#
# This reduces exactly to the legacy form (wse_i - wse_j) at orthogonality
# (c=1, V̂=0) — verified both analytically and numerically across multiple
# synthetic skewed geometries (§10.6). dWSE_n already incorporates the full
# non-orthogonal correction, so the raw L is used directly in the corrected
# kernels below (no separate cos_theta scaling at the kernel level — that
# scaling is now folded into how the caller constructs dWSE_n itself).

"""
    _bates_flux_corrected(q_prev, h_flow, dWSE_n, width, L, n_mann, dt) → Float64

WLSQ-corrected Bates (2010) inertial flux kernel — standard (non-limited)
variant, analogous to `_bates_flux` but using a skewness-corrected driving
head instead of a raw WSE difference scaled by `cos_theta`.

# Arguments
  q_prev   unit discharge from previous timestep (m²/s), signed i→j
  h_flow   flow depth at the edge = max(WSE_i, WSE_j) - z_sill (m), computed
           by the caller from the *raw* (uncorrected) WSE values — only the
           gradient term is corrected, not the flow-depth/friction term
  dWSE_n   skewness-corrected face-normal driving head (m), see module note
           above. Replaces `dWSE / L_eff` in the legacy kernel.
  width    shared edge length (m)
  L        centre-to-centre haversine distance (m) — used as-is; no
           cos_theta scaling (the WLSQ gradient already captures the
           non-orthogonality direction correctly, so this magnitude
           projection is unnecessary and would double-correct)
  n_mann   Manning's roughness (s·m⁻¹ᐟ³)
  dt       timestep (s)

Returns volumetric flux Q (m³/s). Sign convention matches `_bates_flux`:
Q < 0 when WSE_i > WSE_j (flow i→j). Returns 0.0 if h_flow ≤ HFLOW_THRESHOLD.

This function is selected by `step_standard!`/`step_sgs!` Phase A when
`state.gradient_correction == true`; see `_bates_flux` for the legacy
(uncorrected) path used when `false` (current default — see
FloodA5_NonOrthogonal_Correction_Plan.md §10.6.2 for the rationale).
"""
@inline function _bates_flux_corrected(q_prev :: Float64,
                                        h_flow :: Float64,
                                        dWSE_n :: Float64,
                                        width  :: Float64,
                                        L      :: Float64,
                                        n_mann :: Float64,
                                        dt     :: Float64;
                                        q_t             :: Float64 = 0.0,
                                        vector_friction :: Bool    = false)::Float64
    h_flow <= HFLOW_THRESHOLD && return 0.0
    h_flow_safe = max(h_flow, 1e-6)   # floor to avoid h^(10/3) underflow
    L_safe      = max(L, 1.0)         # guard degenerate zero-distance pairs

    numerator   = q_prev - _G * h_flow_safe * dt * dWSE_n / L_safe
    denominator = 1.0 + _G * h_flow_safe * dt * n_mann^2 *
                  _friction_q(q_prev, q_t, vector_friction) / h_flow_safe^(10.0/3.0)
    q_new = numerator / denominator

    return q_new * width
end


"""
    _bates_flux_limited_corrected(q_prev, h_flow, dWSE_n, width, L, n_mann,
                                  dt, depth_donor) → (Q, q_stored)

WLSQ-corrected counterpart of `_bates_flux_limited`. Identical stability
fixes (Froude limiter, volume limiter, consistent `q_stored` — see
`_bates_flux_limited`'s docstring for full rationale on each), applied on
top of the skewness-corrected driving head `dWSE_n` instead of a raw,
cos_theta-scaled WSE difference. Used exclusively by `step_standard!` when
`state.gradient_correction == true`.

See `_bates_flux_corrected` for the argument-by-argument correspondence
with the legacy kernel (`dWSE_n` replaces `wse_i`, `wse_j`, `cos_theta`).
"""
@inline function _bates_flux_limited_corrected(q_prev      :: Float64,
                                                h_flow      :: Float64,
                                                dWSE_n      :: Float64,
                                                width       :: Float64,
                                                L           :: Float64,
                                                n_mann      :: Float64,
                                                dt          :: Float64,
                                                depth_donor :: Float64;
                                                vol_limit_divisor :: Float64 = 5.0,
                                                q_t               :: Float64 = 0.0,
                                                vector_friction   :: Bool    = false
                                                )::Tuple{Float64,Float64}
    if h_flow <= HFLOW_THRESHOLD
        # See _bates_flux_limited's matching comment: zero out stale
        # momentum on dry edges rather than leaving edges.flux[e] unchanged.
        return (0.0, 0.0)
    end
    h_flow_safe = max(h_flow, 1e-6)
    L_safe      = max(L, 1.0)

    numerator   = q_prev - _G * h_flow_safe * dt * dWSE_n / L_safe
    denominator = 1.0 + _G * h_flow_safe * dt * n_mann^2 *
                  _friction_q(q_prev, q_t, vector_friction) / h_flow_safe^(10.0/3.0)
    q_new = numerator / denominator

    # Fix A: Froude limiter — cap at subcritical limit (Fr ≤ FROUDE_LIMIT).
    # h_flow_safe (not the raw h_flow) is used here, consistent with the
    # legacy kernel — the Froude cap is a function of flow depth and is
    # unaffected by the gradient-direction correction.
    q_max = h_flow_safe * sqrt(_G * h_flow_safe) * FROUDE_LIMIT
    q_new = _froude_cap(q_new, q_t, q_max, vector_friction)

    # Fix B: Volume limiter — same vol_limit_divisor convention as
    # _bates_flux_limited (see that function's comment); unchanged in
    # spirit from the legacy kernel, just no longer hardcoded to A5's 5.
    if depth_donor > 0.0
        q_vol_max = (depth_donor * width) / (vol_limit_divisor * dt)
        q_new = clamp(q_new, -q_vol_max, q_vol_max)
    end

    # Fix C: return the post-limiting q so the caller stores a consistent q_prev.
    Q = q_new * width
    return (Q, q_new)
end


# ---------------------------------------------------------------------------
# Manning R-A formulation (SGS)
# ---------------------------------------------------------------------------

"""
    _adjacency_slot(adj_matrix, ci, cj, max_nb) → Int

Return the adjacency slot (1..max_nb) of cell `cj` in cell `ci`'s neighbour list,
as stored in `adj_matrix` (shape max_nb × n_cells, 0 = unused slot).

Returns 1 as a fallback if the pair is not found — this should not occur for
valid EdgeList edges. `initialise_flow_model` now verifies at start-up (SGS
runs only) that every edge's slot is found in both directions, so this
fallback is unreachable on a validated mesh.

`max_nb` defaults to `size(adj_matrix, 1)` (the mesh's own neighbour ceiling).
WP0 fix (2026-10): the default was previously the literal `5`, and `step_sgs!`
called this without the argument, so on a mesh with more than 5 neighbours
slot 6+ was never searched and slot 1 was returned silently — i.e. the wrong
SGS edge curve was used with no error.
"""
@inline function _adjacency_slot(adj_matrix::Matrix{Int},
                                  ci::Int, cj::Int,
                                  max_nb::Int=size(adj_matrix, 1))::Int
    for s in 1:max_nb
        adj_matrix[s, ci] == cj && return s
    end
    return 1   # fallback; unreachable on a mesh validated at init (see above)
end


"""
    _manning_flux_ra(Q_prev, wse_i, wse_j, z_sill, A, R,
                     L, cos_theta, n_mann, dt) → Float64

Manning R-A inertial flux kernel for the SGS R-A solver.

Implements the LISFLOOD-FP SGC formulation:

    Q = (Q_prev - g·A·dt·Sf) / (1 + g·dt·n²·|Q_prev| / (R^(4/3)·A))

where:
  Q_prev   volumetric discharge at previous step (m³/s)
  A        cross-sectional flow area (m²) at the higher-WSE side of the edge
  R        hydraulic radius A/P (m)
  Sf       = dWSE / L_eff  (positive when wse_i > wse_j → Q is negative: flow i→j)

Sign convention: Q < 0 → flow from cell_i to cell_j (i higher);
                 Q > 0 → flow from cell_j to cell_i (j higher).
Callers apply: dV[ci] += Q*dt, dV[cj] -= Q*dt.

This matches the `_bates_flux` sign convention so the Phase B volume scatter
loop is unchanged between Bates and R-A paths.

The R-A form is self-stabilising: as |Q| grows the denominator grows as
n²|Q|/(R^(4/3)·A), providing physically correct Manning friction scaling for
confined channel flow without requiring an explicit Froude or volume limiter.

Returns 0.0 for dry edges: A ≤ 1e-6 m² or h_flow ≤ HFLOW_THRESHOLD.

⚠️  UNCORRECTED FOR DIRECTIONAL BIAS — uses the same L_eff = L × cos θ
    magnitude-only projection as `_bates_flux` (see its docstring for the
    full explanation). Use `_manning_flux_ra_corrected` for production SGS
    simulations; this function is retained for `--gradient-correction off`
    benchmarking.
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
                                   dt        :: Float64;
                                   Q_t             :: Float64 = 0.0,
                                   vector_friction :: Bool    = false)::Float64
    h_flow = max(wse_i, wse_j) - z_sill
    (h_flow <= HFLOW_THRESHOLD || A <= 1e-6) && return 0.0

    R = max(R, 1e-4)

    dWSE  = wse_i - wse_j
    ct    = max(cos_theta, 0.1)
    L_eff = max(L * ct, 1.0)

    # Q_t: face-tangential discharge (m³/s) for vector friction — see _friction_q
    numerator   = Q_prev - _G * A * dt * dWSE / L_eff
    denominator = 1.0 + _G * dt * n_mann^2 *
                  _friction_q(Q_prev, Q_t, vector_friction) / (R^(4.0/3.0) * A)
    return numerator / denominator
end


"""
    _manning_flux_ra_corrected(Q_prev, h_flow, A, R, dWSE_n, L, n_mann, dt) → Float64

WLSQ-corrected counterpart of `_manning_flux_ra`, for the SGS R-A solver.
Uses a pre-computed skewness-corrected driving head `dWSE_n` (built by the
caller from the WLSQ-reconstructed cell-centre WSE gradients and the edge's
skewness vector — see the module note above `_bates_flux_corrected` for the
exact construction) in place of a raw `dWSE` scaled by `cos_theta`.

# Arguments
  Q_prev   volumetric discharge at previous timestep (m³/s), signed i→j
  h_flow   flow depth at the edge = max(WSE_i, WSE_j) - z_sill (m), from
           the *raw* (uncorrected) WSE values — only the gradient term is
           corrected, matching `_bates_flux_corrected`'s convention
  A        cross-sectional flow area (m²) at the higher-WSE side of the edge
  R        hydraulic radius A/P (m)
  dWSE_n   skewness-corrected face-normal driving head (m)
  L        centre-to-centre haversine distance (m) — used as-is; no
           cos_theta scaling (see `_bates_flux_corrected`'s docstring for
           why this is no longer needed)
  n_mann   Manning's roughness (s·m⁻¹ᐟ³)
  dt       timestep (s)

Sign convention matches `_manning_flux_ra`: Q < 0 → flow i→j (i higher).
Returns 0.0 for dry edges: A ≤ 1e-6 m² or h_flow ≤ HFLOW_THRESHOLD.

Selected by `step_sgs!` Phase A when `state.gradient_correction == true`.
"""
@inline function _manning_flux_ra_corrected(Q_prev :: Float64,
                                             h_flow :: Float64,
                                             A      :: Float64,
                                             R      :: Float64,
                                             dWSE_n :: Float64,
                                             L      :: Float64,
                                             n_mann :: Float64,
                                             dt     :: Float64;
                                             Q_t             :: Float64 = 0.0,
                                             vector_friction :: Bool    = false)::Float64
    (h_flow <= HFLOW_THRESHOLD || A <= 1e-6) && return 0.0

    R_safe = max(R, 1e-4)
    L_safe = max(L, 1.0)

    numerator   = Q_prev - _G * A * dt * dWSE_n / L_safe
    denominator = 1.0 + _G * dt * n_mann^2 *
                  _friction_q(Q_prev, Q_t, vector_friction) / (R_safe^(4.0/3.0) * A)
    return numerator / denominator
end


# ---------------------------------------------------------------------------
# CFL timestep
# ---------------------------------------------------------------------------

"""
    _cfl_dt(state, method) → Float64

Compute the maximum stable timestep from the CFL (Courant–Friedrichs–Lewy)
condition for the inertial shallow-water equations:

    dt ≤ courant × dx / √(g × h)

where `dx` is the cell length scale (√area), `h` is the local water depth,
and `courant = 0.7` matches the CAESAR-Lisflood default.  This is the
wave-speed stability criterion appropriate for the inertial (momentum-
retaining) formulation of Bates et al. (2010) — replacing the earlier
diffusive-wave form which was too conservative at shallow depths and not
conservative enough at larger depths.

The per-cell `dx = √(cell_area)` accounts for the varying sizes of A5
pentagonal cells.  The 99th-percentile wet-cell depth is used rather than
the absolute maximum so that a single extreme cell (e.g. an overfull channel
cell extrapolating above z_max) cannot force a tiny dt for the whole domain.

**`state.cfl_mode === :local`** (WP2, 2026-10; `--cfl-mode local`). The
legacy form pairs the *smallest* cell anywhere in the mesh with the 99th-
percentile depth of *any* wet cell, which is neither conservative (the deepest
1% are ignored) nor local. The local form is a per-cell max-reduction:

    dt = courant × min_i ℓ_i / √(g × h_i)     over wet cells (h_i > 1e-4 m)

with ℓ_i = `state.cell_length[i]` = 4A_i/P_i (the hydraulic diameter: equal to
dx on a square cell, so the Courant number keeps its LISFLOOD-FP meaning;
the review's 2A/P is half of this and would halve dt at the same courant).
The SGS terrain-range depth cap still applies, so overfull channel cells do
not set dt. Falls back to √A if `cell_length` is empty (hand-built FlowState).

Note: this function reads FlowState fields and is therefore not a pure scalar
kernel.  It lives in flow2d.jl because it is tightly coupled to the physics
constants and CFL formulation, but it cannot be used independently of FlowState.
"""
function _cfl_dt(state::FlowState, method::FlowMethod; courant::Float64=0.7)::Float64
    state.cfl_mode === :local && return _cfl_dt_local(state, method, courant)
    dx_min = Inf
    for i in eachindex(state.cell_ids)
        dx = sqrt(state.cell_area[i])
        dx < dx_min && (dx_min = dx)
    end

    # Collect wet-cell depths for the CFL calculation.
    # For SGS, water_depth now holds the true hypsometric depth (uncapped) for display.
    # We apply the terrain-range cap here so that overfull channel cells do not force
    # the whole domain to tiny dt — the cap is (z_max - z_min) per cell.
    # The 99th-percentile is used instead of the absolute maximum for the same reason.
    # WP0 fix (2026-10): was `method == SGSFlow`, which compares an instance
    # (e.g. SGSFlow(0.03)) with the type and is therefore always false — the
    # SGS terrain-range depth cap below was never applied, so overfull
    # channel cells set dt for the whole domain. `isa` is the intended test.
    is_sgs = method isa SGSFlow
    wet_depths = Float64[]
    for i in eachindex(state.cell_ids)
        h = state.water_depth[i]
        h <= 1e-4 && continue
        if is_sgs && !isempty(state.sgs_tables)
            tbl = state.sgs_tables[i]
            if !isnan(tbl.z_min)
                h = min(h, max(0.0, tbl.z_max - tbl.z_min))
            end
        end
        h > 1e-4 && push!(wet_depths, h)
    end

    isempty(wet_depths) && return 60.0   # all cells dry — return safe default

    sort!(wet_depths)
    pct99_idx = max(1, round(Int, 0.99 * length(wet_depths)))
    h_cfl = wet_depths[pct99_idx]

    h_cfl < 1e-6 && return 60.0

    return courant * dx_min / sqrt(_G * h_cfl)
end


"""
    _cfl_dt_local(state, method, courant) → Float64

Per-cell CFL max-reduction for `cfl_mode === :local` — see `_cfl_dt`.
"""
function _cfl_dt_local(state::FlowState, method::FlowMethod, courant::Float64)::Float64
    is_sgs   = method isa SGSFlow && !isempty(state.sgs_tables)
    has_len  = length(state.cell_length) == length(state.cell_ids)
    rate_max = 0.0                            # max √(g h_i)/ℓ_i  (s⁻¹)
    for i in eachindex(state.cell_ids)
        h = state.water_depth[i]
        h > 1e-4 || continue                  # also skips NaN
        if is_sgs
            tbl = state.sgs_tables[i]
            isnan(tbl.z_min) || (h = min(h, max(0.0, tbl.z_max - tbl.z_min)))
            h > 1e-4 || continue
        end
        ℓ = has_len ? state.cell_length[i] : sqrt(state.cell_area[i])
        ℓ > 0.0 || continue
        rate = sqrt(_G * h) / ℓ
        rate > rate_max && (rate_max = rate)
    end
    rate_max > 0.0 || return 60.0             # all cells dry — safe default
    return courant / rate_max
end


# ---------------------------------------------------------------------------
# Positivity-preserving donor scaling (WP2)
# ---------------------------------------------------------------------------

"""
    _positivity_scale!(edge_vol, q_store, volume, cell_i, cell_j)

Positivity-preserving donor scaling (WP2, 2026-10; `--limiter positivity`),
replacing the per-edge volume limiter (Fix B, h_donor·w/(N·dt)) and the
per-edge donor cap (V/(2N)). Both of those depend on the face count and edge
width and cut individual faces, so they act anisotropically
(FloodA5_Review_and_ForwardPlan_2026-10.md §4.2).

Instead, after all face fluxes are computed, every outflow of a cell is scaled
by the same factor

    s_i = min(1, V_i / Σ_outflows(i) · dt)

so that no cell can lose more than it holds. This is mass-conservative (each
edge volume is applied to both cells after scaling), isotropic (all of a
cell's outflows are scaled together) and mesh-general — the standard wet–dry
positivity treatment in finite-volume schemes.

`edge_vol[e]` is the signed edge volume for this step (> 0: cell_j donates to
cell_i, the Phase B convention). `q_store` is the momentum array written in
Phase A (`edges.flux` for standard flow, `edges.flux_Q` for SGS R-A); scaled
edges are written back so the stored momentum matches the transfer (Fix C).
`volume` is the start-of-step volume. Serial (Phase B is serial).
"""
function _positivity_scale!(edge_vol :: Vector{Float64}, q_store :: Vector{Float64},
                            volume   :: Vector{Float64},
                            cell_i   :: Vector{Int}, cell_j :: Vector{Int})
    out = zeros(Float64, length(volume))
    @inbounds for e in eachindex(edge_vol)
        ev = edge_vol[e]
        if ev > 0.0
            out[cell_j[e]] += ev
        elseif ev < 0.0
            out[cell_i[e]] -= ev
        end
    end
    @inbounds for e in eachindex(edge_vol)
        ev = edge_vol[e]
        iszero(ev) && continue
        d = ev > 0.0 ? cell_j[e] : cell_i[e]
        o = out[d]
        V = max(volume[d], 0.0)
        if o > V
            s = V / o                          # o > V ≥ 0, so o > 0
            edge_vol[e] = ev * s
            q_store[e] *= s
        end
    end
    return edge_vol
end


# ---------------------------------------------------------------------------
# Open boundary / ghost-edge flux kernels
# ---------------------------------------------------------------------------

"""
    _bates_ghost_flux(q_prev, wse_ci, wse_ghost, z_sill, width, L,
                      n_mann, dt, depth_ci) → (Q_out, q_stored)

Compute outflow flux across a virtual ghost-cell boundary edge using the
Bates et al. (2010) inertial formulation with stability limiters.

`wse_ghost` is computed by the caller via `_ghost_wse()` based on the BC type:
  ZeroGradient: wse_ghost = wse_ci  (dWSE = 0; momentum carry-out only)
  Critical:     wse_ghost = sill + (2/3)(wse_ci - sill)

The function enforces one-way flow: Q_out is always ≥ 0. Ghost cells cannot
inject water into the domain — any negative Q (which would represent inflow
from outside) is clamped to zero and q_stored is reset.  This prevents
numerical artefacts from the zero-gradient approximation at inflow edges
(where the domain WSE may briefly dip below the ghost WSE due to routing).

Returns `(Q_out, q_stored)` where:
  Q_out    — volume flux leaving the domain through this edge (m³/s, ≥ 0)
  q_stored — unit discharge (m²/s) to persist as flux_prev next step (Fix C)

Sign convention: Q_out > 0 → water leaves cell ci.
Caller applies: state.volume[ci] -= Q_out * dt
"""
@inline function _bates_ghost_flux(q_prev    :: Float64,
                                    wse_ci    :: Float64,
                                    wse_ghost :: Float64,
                                    z_sill    :: Float64,
                                    width     :: Float64,
                                    L         :: Float64,
                                    n_mann    :: Float64,
                                    dt        :: Float64,
                                    depth_ci  :: Float64)::Tuple{Float64,Float64}
    # Closed BC sentinel
    wse_ghost == -Inf && return (0.0, 0.0)

    # One-way enforcement: only allow outflow.
    # If the ghost WSE ≥ interior WSE, water would flow inward — block it.
    # This also handles the ZeroGradient case (wse_ghost == wse_ci) where
    # the momentum term q_prev carries any remaining outflow; if q_prev is
    # also zero or negative the kernel correctly returns near-zero.
    wse_ghost >= wse_ci + 1e-10 && return (0.0, 0.0)

    # Call _bates_flux_limited with (wse_ghost, wse_ci) as the (i, j) pair.
    # wse_i = wse_ghost < wse_ci = wse_j, so dWSE = wse_ghost - wse_ci < 0
    # → flow from j to i (outward) → Q_raw > 0 in standard convention.
    # We return Q_raw directly as Q_out (positive = water leaves domain).
    Q_out, q_stored = _bates_flux_limited(q_prev, wse_ghost, wse_ci, z_sill,
                                           width, L, 1.0, n_mann, dt, depth_ci)

    # Defensive clamp: should be positive by construction above, but guard
    # against floating-point edge cases.
    Q_out < 0.0 && return (0.0, 0.0)

    return (Q_out, q_stored)
end


"""
    _manning_ghost_flux(Q_prev, wse_ci, wse_ghost, z_sill, A, R,
                        width, L, n_mann, dt) → (Q_out, Q_stored)

Ghost-edge outflow flux using the Manning R-A inertial kernel (SGS path).

`A` and `R` are evaluated at `wse_ci` (using the boundary cell's SGS table
at the ghost edge slot) since the ghost cell has no terrain data.

Returns `(Q_out, Q_stored)` in m³/s; Q_out ≥ 0 (one-way enforcement).
"""
@inline function _manning_ghost_flux(Q_prev    :: Float64,
                                      wse_ci    :: Float64,
                                      wse_ghost :: Float64,
                                      z_sill    :: Float64,
                                      A         :: Float64,
                                      R         :: Float64,
                                      L         :: Float64,
                                      n_mann    :: Float64,
                                      dt        :: Float64)::Tuple{Float64,Float64}
    wse_ghost == -Inf            && return (0.0, 0.0)
    wse_ghost >= wse_ci + 1e-10  && return (0.0, 0.0)

    # (wse_ghost, wse_ci) as (i, j): wse_ghost < wse_ci → flow outward → Q_out > 0
    Q_out = _manning_flux_ra(Q_prev, wse_ghost, wse_ci, z_sill, A, R,
                              L, 1.0, n_mann, dt)
    Q_out < 0.0 && return (0.0, 0.0)
    return (Q_out, Q_out)
end
