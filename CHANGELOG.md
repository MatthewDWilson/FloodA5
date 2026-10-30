# Changelog

All notable changes to FloodA5 are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning follows
[Semantic Versioning](https://semver.org) — see `VERSIONING.md` for the
project's specific policy and the reasoning behind starting at 0.1.0.

## Historical development (pre-0.1.0)

FloodA5 was developed over many sessions before formal versioning began.
This is a brief prose summary of that trajectory, not a version-by-version
log — no version numbers are retrofitted to this period (see
`VERSIONING.md`). Full detail is retained in the project's internal
development history.

Development proceeded roughly through these phases: core A5 mesh
generation and the Bates (2010) inertial flow solver on a flat/synthetic
domain; sub-grid sampling (SGS) using hypsometric volume curves derived
from LiDAR, added to resolve sub-cell channel features the standard
solver's cell-mean elevation cannot see; a stability investigation
(prompted by checkerboard oscillation observed on the Carlisle domain)
that identified and fixed an inconsistent momentum-storage pattern and
added Froude/volume flux limiters, informed by a line-by-line comparison
against LISFLOOD-FP and CAESAR-Lisflood reference implementations; an
alternative hydraulic-radius (R-A) flux formulation for SGS edges; dynamic
fluvial inflows and open/closed boundary conditions with LISFLOOD-FP
`.bci`/`.bdy` format support; and — the most extensive single line of
work — a multi-session investigation into a systematic directional flow
bias on the A5 pentagonal mesh, producing three independent, opt-in
correction candidates (a WLSQ gradient correction, a face-local "diamond"
reconstruction, and a cell-vector momentum representation), none of which
fully resolves the bias on its own.

## [Unreleased]

## [0.3.0] — 2026-10-30

⚠ Default results change for `--flow-model standard` (MINOR). SGS and SWE
defaults are unchanged.

### Changed
- ⚠ `--scheme recommended` (the default) for standard flow is now the
  Δt-consistent scheme: positivity limiter (also at open boundaries), no
  Q-centring (θ = 1), η filter C_s = 0.05, with the diamond-wetdry face
  gradient, vector friction and local CFL as before. It is the
  local-inertial configuration closest to the full shallow-water reference
  on A5 and Cartesian meshes and converges with resolution; the 0.2.x
  default stored 60–80% more on the Carlisle test and did not converge
  (HYDRAULICS.md §15.11).
- `--boundary-limiter`, `--eta-filter`, `--momentum-damping` and
  `--momentum-damping-scale` take their defaults from `--scheme`
  (`scheme_options` gains `eta_filter`, `momentum_damping`,
  `momentum_damping_scale`; `scheme_boundary_limiter`). Momentum damping is
  off in every scheme.

### Added
- `--scheme v0.2`: the 0.2.x default, bit-identical to 0.2.5
  (`scheme_options(:v02, …)`).

### Known trade-off
- More velocity noise in deep channels than the 0.2.x default (dv_rms
  0.65 vs 0.24 m/s on A5 Carlisle); `--momentum-damping 0.05` reduces it.

## [0.2.6] — 2026-10-30

Default results unchanged (PATCH): only the opt-in `--flow-model swe` changes.

### Fixed
- ⚠ `--flow-model swe` (the R1 reference) on thin films over a slope
  (rain-on-grid). The second-order hydrostatic reconstruction used
  cell-centre beds at the faces, removing the bed step from every downstream
  face depth; it now reconstructs the face bed (z_L = η_L − h_L, with a
  limited depth gradient) with the centred bed source of Audusse et al.
  (2004, §4), and the thin-film depth bound is no longer needed. The
  semi-implicit Manning factor took |q| from the predicted state, which
  moved its fixed point and slowed a 2 cm sheet by 36–58%; it now uses |q|
  from the stage start. Carlisle 6 h storage with the reference falls from
  3.44 to 2.69 Mm³ (A5) and 2.63 to 2.33 Mm³ (Cartesian).

### Added
- `test/test_swe_reference.jl` S10: thin sheet flow at normal depth.

### Findings (see `HYDRAULICS.md` §15.10, rewritten)
- Against the corrected reference, every local-inertial scheme is more
  dissipative; the η filter alone is closest (Carlisle storage +5%
  Cartesian, +26% A5) and the momentum damping only adds to it. The 0.2.5
  conclusion that momentum damping 0.05 matched the reference does not hold.
- Resolution (Cartesian 40/20/10 m): the reference converges (3.15/2.33/2.11
  Mm³); the η-filter scheme converges alongside it, 4–8% above; the 0.2.x
  default does not converge (3.80 → 3.83 from 20 to 10 m).

## [0.2.5] — 2026-10-29

Default results unchanged (PATCH): the default-scheme Carlisle run is
bit-identical to 0.2.4 on A5 and Cartesian (`--cell-order file`, 4 h, every
frame).

### Added
- `--flow-model swe`: a full shallow-water reference solver (R1 of the
  deep-channel damping work) on any polygon mesh. HLLC fluxes with
  hydrostatic reconstruction (well-balanced, positivity-preserving); MUSCL
  with WLSQ gradients and the Barth–Jespersen limiter plus SSP-RK2
  (`--swe-order 2`, default; `1` for first order); point-implicit Manning
  friction; optional eddy viscosity `--swe-viscosity C_E`. It uses the same
  volume, sources, rainfall and boundary types as the other solvers; open
  boundaries are outflow-only, as in the standard solver. Struct-of-arrays
  per-face/per-cell kernels without atomics, as a GPU backend needs
  (`surfacewater/swe_reference.jl`, `swe_types.jl`). FlowState field 42
  `swe`, with a 41-argument constructor shim; GhostEdge gains the outward
  face normal and midpoint offset (8-argument constructor kept).
- `test/test_swe_reference.jl` (S1–S9): lake at rest, mass and positivity,
  Stoker and Ritter dam breaks, open-boundary mass balance, uniform Manning
  flow, eddy viscosity, HDF5 velocity. Determinism fingerprints
  `r1_swe_order2`, `r1_swe_order1_visc`.
- `test/bench_swe_dambreak.jl` (circular dam break, A5 vs Cartesian vs a
  fine grid) and `test/diag_vs_reference.jl` (storage, rms depth, extent CSI
  and deep-water speed of runs against a reference run on the same mesh).

### Findings (see `HYDRAULICS.md` §15.10)
- Carlisle, against the reference on the same mesh: the recommended scheme
  stores 42–44% more and moves deep water at less than half the speed. The
  Δt-consistent scheme with the η filter and weak momentum damping (0.05,
  wave) matches storage within 1–2% on A5 and Cartesian.
- The reference converges in Courant number and is insensitive to round-off
  and to a physical eddy viscosity, but not in resolution: 40/20/10 m
  Cartesian grids store 3.56/2.63/2.18 Mm³ at 6 h. A5 res-18 behaves like a
  28–37 m Cartesian grid on Carlisle's narrow channels; on a smooth dam break
  both meshes are equally accurate.

## [0.2.4] — 2026-10-28

Default results unchanged (PATCH).

### Added
- `--momentum-damping C_m` and `--momentum-damping-scale wave|flow`
  (experimental, standard flow with the diamond gradient and vector
  friction; D6′ of the deep-channel damping work). A divergence-free damping
  of the circulating part of the face discharges: the vertex circulation
  drives a streamfunction correction, so it moves no water, leaves the
  gravity-wave dynamics alone and vanishes for uniform flow. Eddy viscosity
  ν = C_m·(√(gh) or |u|)·ℓ, a rate in physical units
  (`surfacewater/momentum_damping.jl`). Default 0 (off). FlowState fields
  40–41, with a 39-argument constructor shim. `test/test_momentum_damping.jl`
  (M1–M5); determinism fingerprint `d6_momentum_damping`; configurations in
  `test/diag_step_spectrum.jl`.

### Findings (see `HYDRAULICS.md` §15.9)
- On the Δt-consistent baseline with the η filter, momentum damping removes
  the remaining velocity noise (A5 Carlisle dv_rms 0.65 → 0.18 m/s at wave
  0.2; recommended 0.24) and keeps the runs Δt-converged. Cartesian twin
  divergence falls to 10⁻¹⁰ m.
- It is an eddy viscosity: lateral momentum exchange slows the channels and
  raises 6 h storage by 2–25% depending on C_m, at ν one to two orders above
  physical values. On A5 it adds a small planar asymmetry (up to 0.03).
  Choosing C_m needs a full shallow-water reference.

## [0.2.3] — 2026-10-28

Default results unchanged (PATCH).

### Added
- `--eta-filter C_s` (experimental, standard flow; D4 of the deep-channel
  damping work). A conservative free-surface high-pass filter: after the
  volume update, water is exchanged across wet faces in proportion to the
  jump in each cell's residual from a local planar fit of the WSE, at a rate
  C_s·√(g h)/ℓ that does not depend on Δt (`surfacewater/eta_filter.jl`). It
  leaves a lake at rest and planar surfaces untouched, conserves volume,
  keeps volumes non-negative, and has no wet/dry switch. Default 0 (off).
  `test/test_eta_filter.jl` (E1–E6); a determinism fingerprint
  (`d4_eta_filter`); filter configurations in `test/diag_step_spectrum.jl`.

### Findings (see `HYDRAULICS.md` §15.8)
- Why the square grid needs damping: the forward–backward step is neutral
  for every resolved mode. At Courant 0.7 the Cartesian checkerboard sits at
  98% of its stability limit (C* = 1/√2). A5 has no checkerboard mode
  (C* = 0.876). Damping the discharge, as Q-centring does, lowers the limit
  (κ < 4 − 2β); damping the water surface raises it (κ < (4 − 2β)/(1 − β)).
- In steady uniform flow 3 m deep at Δt = 8 s, the recommended scheme loses
  11% of the Manning conveyance (8° deflection) to the volume limiter. The
  positivity limiter is exact, with or without the filter.
- With `--limiter positivity --q-centre-theta 1 --boundary-limiter positivity
  --eta-filter 0.05` on Carlisle:
  - A5 stays Δt-converged (3.381–3.386 Mm³) and reaches the recommended
    scheme's surface-extrema level.
  - Cartesian goes from chaotic (twin divergence 56 cm) to reproducible at
    round-off (2.3 µm).
  - The remaining gap is velocity noise in shallow sheet flow, which a
    surface filter cannot reach.

## [0.2.2] — 2026-10-25

Default results unchanged (PATCH).

### Added
- `--limiter-stats` (diagnostic, D1b of the deep-channel damping work): records
  per edge and step which limiter or switch was active (Froude cap, volume
  limiter, donor cap, positivity scaling, diamond wet/dry fallback). Prints an
  hourly log line and an end-of-run table by flow-depth band, including the
  rate at which edges toggle between them. Results are unchanged; with stats
  on, Phase A visits every edge (`test/test_limiter_stats.jl`, L1–L2).
- `--boundary-limiter legacy|positivity` (experimental): `positivity` removes
  the volume limiter and the per-step donor cap from open-boundary outflow,
  keeping only "no more than the cell holds". Default `legacy`
  (`test/test_limiter_stats.jl`, L3).

### Findings (see `HYDRAULICS.md` §15.7)
- The volume limiter binds in deep water on A5, 30–90 times more often at
  Courant 0.7 than at 0.35. SGS's donor cap binds on half of all wet
  edge-steps. Both scale with 1/Δt.
- With Q-centring off and both limits replaced by positivity scaling, in the
  interior and at the boundary, A5 and SGS converge in Δt: storage at 6 h is
  3.382 / 3.382 / 3.386 Mm³ at Courant 0.7 / 0.5 / 0.35 (A5), against
  +10% to −20% trends otherwise. The boundary limits were also a major noise
  source. The square grid without damping is chaotic.

## [0.2.1] — 2026-10-23

Default results unchanged (PATCH, under the updated `VERSIONING.md`).

### Added
- Damping diagnostics (D1 of the deep-channel damping work):
  `test/diag_step_spectrum.jl` (eigenvalues of the finite-difference Jacobian
  of the real step — gravity operator, friction, limiters and Q-centring
  together — at a lake at rest or a frozen sloping flow);
  `test/diag_operator_spectrum.jl` (the free-surface operator alone:
  two-point, diamond, blends and the symmetric part of the diamond operator,
  with growth factors and a linear-field consistency residual);
  `test/diag_divergence.jl` (frame-by-frame divergence of two runs: twin runs
  with a 1e-10 perturbation, cell orders, or Courant numbers) and
  `test/diag_d2_summary.sh`.

### Changed
- `VERSIONING.md`: PATCH now covers any change that leaves default results
  unchanged (opt-in options, diagnostics, documentation, bit-identical
  performance work); MINOR is reserved for changed defaults/results or output
  format and gets a GitHub release; every merged PR bumps the version in its
  final commit and is tagged.

### Findings (see `HYDRAULICS.md` §15.6)
- On a square grid, Q-centring θ = 0.9 at Courant 0.7 is linearly unstable in
  deep water (|μ| = 1.6 per step): its damping of the checkerboard mode lowers
  the forward–backward stability limit to κ ≤ 4 − 2β, below κ ≈ 3.92.
- Twin runs: the round-off sensitivity comes from discrete switches — the
  legacy CFL's 99th-percentile depth (Cartesian twin divergence 45 cm instead
  of 3 µm) and the diamond-wetdry stencil switch — plus the diamond operator's
  slow linear growth. With the two-point gradient and the local CFL, A5 is as
  reproducible as Cartesian (6×10⁻⁷ m after 6 h).
- The Courant sensitivity of storage (4–20% between Courant 0.7 and 0.35) is
  systematic and appears in every configuration, SGS included.

## [0.2.0] — 2026-10-23

The October 2026 work packages (WP0–WP4 of the internal forward plan, plus
the Cartesian formulation and pya5 1.x support). A MINOR release under
`VERSIONING.md`: new CLI options and **changed default behaviour**.

### Upgrading from 0.1.0 — read this first
- ⚠ **Results change.** The default numerical scheme is now
  `--scheme recommended` (diamond-wetdry face gradient, vector friction,
  hybrid limiter, local CFL and consistent Q-centring for standard flow;
  vector friction, local CFL and no Q-centring for SGS). `--scheme legacy`
  reproduces the 0.1.0 numerics.
- ⚠ **Two bugs affected 0.1.0 results.** On A5, Q-centring added most of its
  partner faces with the wrong sign and acted as an artificial, dt-dependent
  resistance (Carlisle stored about twice the water of a Cartesian reference);
  and SGS edge tables were attached to the wrong neighbour for most cells.
  SGS meshes built with 0.1.0 are refused until rebuilt or repaired
  (`mesh/repair_sgs_slot_order.py`).
- **pya5 1.x** changed the A5 cell-ID encoding: meshes record their ID scheme
  and old ones are refused; migrate with `julia setup.jl --upgrade` then
  `julia setup.jl --migrate-meshes`.
- **Cell order.** Simulations keep cells in a spatial (Morton) order
  (`--cell-order spatial`, default); output is still in the mesh file's order.
  This changes results at round-off level — but the standard solver
  amplifies round-off in deep, noisy channels (HYDRAULICS.md §15.5).
- **Known open issue:** grid-scale water-surface noise in deep channels with
  the standard solver, and the sensitivity of Carlisle results to the Courant
  number (HYDRAULICS.md §15.3–15.5). SGS is not affected in the same way.

### Highlights
- Directional bias explained and corrected (diamond face gradient + vector
  friction): planar-slope N/S asymmetry on A5 from 0.87 to ~0.05.
- Mesh-general solver: Cartesian (including raster-aligned grids for
  LISFLOOD-FP comparisons), hexagonal, rotated-Cartesian, Voronoi and imported
  polygon meshes; SGS on any number of neighbours.
- CPU performance (WP4): Carlisle res-18, inflows + rain, the same options,
  about 4–5× faster simulation than before WP4 (standard 57.9 → ~12 s,
  SGS 169 → ~40 s for 2 h at 8 threads); SGS mesh load 75 → 4.5 s;
  allocation per step from MB to under 2 KB; only the wet part of the domain
  is computed.
- Test infrastructure: GitHub Actions CI (pya5 0.7.3 and 1.0.0), determinism
  checks, frozen reference steps with bit-equality tests.

The per-work-package entries below give the details and evidence.

### WP4 stage 3: active-edge set (`wp4c-active`)

Results unchanged: bit-identical to the full pass (`test/test_wp4_stage2.jl`
R4: a partly wet domain that wets, drains and dries, against the frozen
reference, mass balance included; Carlisle runs identical to `main` frame by
frame; determinism fingerprints unchanged).

#### Changed
- Phase A visits only edges that can change: those of cells holding water,
  last step's set (an edge that has just dried is visited once more, so the
  kernel itself zeroes its flux and edge volume) and edges that can carry flux
  while dry (none on normally built meshes). A full pass runs on the first
  step, when elevations/sills/areas change, and while more than half the
  cells hold water (`ACTIVE_MAX_FRACTION`). The volume gather still covers
  every edge, so no transfer can be dropped. See `_active_edges!`.
- The diamond vertex pass covers only the vertices of the active edges; the
  WSE gradient (`_compute_wse_gradients!`) is threaded.

#### Performance (Carlisle res-18, `--threads auto`, against stage 2)
- River inflows only, 6 h (ends with 6.8k of 30k cells wet): standard
  simulation 24.4 → 17.7 s, SGS 37.7 → 29.2 s.
- Nearly dry domain: standard step 2.1 → 0.6 ms (8 threads).
- Fully wet domains (e.g. uniform rain) are unchanged — they take the full pass.

### WP4 stage 2: kernels, cell order, non-allocating threads (`wp4b-stage2`)

#### Changed
- ⚠ **Cells are kept in spatial (Morton) order during the simulation**
  (`--cell-order spatial`, the new default; `--cell-order file` keeps the
  mesh file's order). Neighbouring cells sit close together in memory, which
  makes a Carlisle res-18 step ~1.4–1.7× faster. HDF5 output is still written
  in the file's cell order. Results change at round-off level only — but see
  HYDRAULICS.md §15.5: with the default standard-flow options, round-off
  differences grow to centimetres in noisy deep-channel cells (SGS: 10⁻⁵ m).
  Not applied with legacy Q-centring (`--scheme legacy`), whose result depends
  on edge orientation and hence on cell order.
- Phases B and C (volume scatter and update) are one threaded per-cell gather
  over a cell→edge CSR in edge order; the per-edge donor cap moved into Phase
  A; the positivity limiter's totals are gathered too. Bit-identical.
- The edge flux loops are kernels (`surfacewater/kernels.jl`: `StdOpts`/
  `SGSOpts` options as type parameters, arrays passed explicitly) — the form
  the GPU backend (WP5) will launch. Bit-identical.
- Threaded loops in the time step use Polyester `@batch per=core` instead of
  `Threads.@threads`: no per-step allocation (0.7–2 KB/step, was 15–20 KB).
  `per=core` matters: Polyester's default spins on every logical core and,
  under `--threads auto`, made steps ~4× slower. Polyester is now a direct
  dependency (it was already in the Manifest).

#### Added
- `surfacewater/reference_steps.jl`: frozen copies of the stage-1
  `step_standard!`/`step_sgs!` (`step_*_reference!`); `test/test_wp4_stage2.jl`
  checks the restructured steps against them bit for bit, and the cell order.
- `--cell-order spatial|file`; `spatial_cell_order`, `permute_mesh`. With a
  spatial order, SGS HDF5 output also records `/mesh/edge_cells` (each edge's
  two cells, file-order indices) because `flux_Q` follows the internal edges.

#### Performance (Carlisle res-18, 8 threads, against stage 1)
- Step (wet state): standard 4.0 → 2.4–2.9 ms, SGS 17.8 → 9.9–10.7 ms.
- Inflows + rain, 2 h: standard simulation 14.5 → 11.5–13.5 s, SGS 57 → 39–41 s.
- WP0 point source (`--threads auto`, legacy scheme, file order): unchanged
  (standard 3.2 → 2.9 s, SGS 3.9 → 4.3 s); now dominated by start-up and
  first-call compilation.

### WP4 stage 1: CPU performance (`wp4-cpu-perf`)

Results are unchanged: determinism fingerprints identical (1 and 4 threads),
Carlisle res-18 standard and SGS outputs identical frame by frame
(max |Δ| = 0), `test/test_wp4_equivalence.jl` checks each fast path against
the code it replaces.

#### Changed
- `FlowState.sgs_tables`, `ghost_edges`, `ghost_cell_bc` and `diamond_table`
  are concretely typed (were `Any`; `mesh/DiamondFlux.jl` is now included
  before `FlowState`). Fixes per-edge dynamic dispatch, e.g. SGS CFL 12 → 0.1 ms.
- Per-step work arrays live in a reused `StepWorkspace` (new last `FlowState`
  field, backward-compatible constructor); SGS edge slots precomputed; one
  bracket search per side for the SGS area and radius
  (`A5Grid.edge_area_radius`); the unused SGS wetted-area pass removed.
- Diamond face gradient: vertex values computed once per vertex per step from
  a flat (CSR) vertex table (`diamond_vertex_eta!`/`_wet!`) instead of two
  Dict lookups per edge per use.
- Dry edges (max WSE − sill ≤ HFLOW_THRESHOLD, where both kernels return
  zero) skip the momentum and gradient work.
- Velocity is computed only for written or visualised steps
  (`step_standard!/step_sgs!(…; velocity = false)` from `run_simulation!`).
- Ghost edges updated in place; `_positivity_scale!` takes a work buffer.
- GeoParquet loading passes columns as binary arrays instead of JSON text:
  Carlisle res-18 SGS mesh load 75 → 4.5 s (identical mesh).
- CUDA is loaded on first GPU-sized point-in-polygon call, not at start-up
  (`FLOODA5_NO_CUDA=1` skips it): `include("FloodModel.jl")` 11.4 → 6.4 s.
- The per-10-step debug volume line is now printed every 500 steps.

#### Performance (Carlisle res-18, 8 threads, against `main` the same day)
- Inflows + rain, 2 h, `--scheme recommended`: standard simulation
  57.9 → 14.5 s (4.0×), SGS 169 → 57 s (3.0×); SGS total wall 4 m 7 s → 1 m 4 s.
- WP0 point-source benchmark (`--scheme legacy --gradient-correction off`,
  10 h): SGS total wall 1 m 48 s → 10.5 s (10.3×; simulation 8.0×);
  standard simulation 4.8 → 3.2 s (now dominated by first-call compilation
  and HDF5 output).
- Allocation per step 1.6–3.1 MB (standard) and 13.9 MB (SGS) → 12–24 KB at 8
  threads (task overhead of the threaded loops; 2–4 KB single-threaded).

### WP2c: Q-centring on A5 (`wp2c-qcentre`)

Evidence: `docs/HYDRAULICS.md` §15.4.

#### Changed
- ⚠ **`--scheme recommended` now uses consistent Q-centring for standard flow
  and no Q-centring (θ = 1) for SGS.** `--q-centre-theta` joins the options
  taken from the scheme. `--scheme legacy`, or `--q-centre-mode legacy
  --q-centre-theta 0.9`, restores the previous Q-centring. On Carlisle res-18
  this lowers the water stored at 12 h from 12.2 to 8.9 Mm³ (standard; 20 m
  Cartesian reference 6.1) and from 11.3 to 8.0 Mm³ (SGS). Local WSE extrema
  on Carlisle rise from 2.7% to 6.6% (standard; rms residual 61 → 41 mm) and
  fall from 1.2% to 0.4% (SGS).

#### Found
- ⚠ **Q-centring partners are sign-reversed on A5.** `_q_centred` adds the
  stored fluxes of each edge's two "most collinear" faces, but edges are
  oriented from the lower to the higher cell index and A5 cell order is not
  geometric, so 67–68% of partners have the opposite orientation (0% on
  Cartesian grids). On A5 the scheme therefore acts mostly as a resistance of
  order (1−θ)/dt: in the frozen-coefficient test steady sheet-flow conveyance
  is 16% (dt 1 s) and 26% (dt 5 s) of the Manning value; on Carlisle res-18
  the domain stores about twice the water of a 20 m Cartesian grid on the
  same DEM and forcing, and more as dt is reduced. This affects every A5 run
  with θ < 1 (the default 0.9), standard and SGS.

#### Added
- `--q-centre-mode legacy|oriented|consistent` (`q_centre_mode`, new last
  `FlowState` field with a backward-compatible constructor):
  `oriented` corrects the partner signs (`_orient_collinear!`; signed indices
  in `collinear_i/j`) — on pentagons the partners are not parallel and this
  deflects the flow by up to 57°; `consistent` also maps each partner onto the
  edge normal using the cell vectors (`_q_centred_consistent`): exact for
  uniform flow (frozen test: conveyance within 0.1%, deflection < 0.4°),
  identical to Q-centring on Cartesian grids, standard and SGS kernels.
- `--courant C` (default 0.7).
- `test/test_qcentre.jl` (in CI).

### WP2b: the recommended scheme becomes the default (`wp2b-default-scheme`)

#### Changed
- ⚠ **New CLI defaults** (`--scheme recommended`, `scheme_options`): unset
  options now take the recommended values for the flow model —
  standard flow: `--gradient-correction on --face-flux-method diamond-wetdry
  --friction-model vector --limiter hybrid --cfl-mode local`; SGS:
  `--friction-model vector --cfl-mode local`. **`--scheme legacy` restores the
  previous behaviour**; explicit flags always win. `run_flood_model` follows the
  same rule (its option keywords default to `nothing` = "from the scheme");
  `initialise_flow_model`'s own defaults stay legacy. The local CFL takes
  ~20–65% more steps. Checks behind the choice (HYDRAULICS.md §15.3): open
  boundaries (Carlisle inflows + rain; planar slope under rain), SGS on the
  repaired Carlisle tables and under rain. Known trade-off: the diamond face
  gradient adds interior velocity noise in shallow sheet flow (0.05 vs
  0.015 m/s on the rained planar slope) while making the water surface
  smoother and removing the directional bias.
- `test/test_default_scheme.jl` (in CI).

### WP3 (part 2): SGS on general meshes, and the SGS slot-order fix (`wp3b-sgs-general`)

#### Fixed
- ⚠ **SGS edge sills and R–A edge curves belonged to the wrong neighbour for
  most cells** (`build_sgs_tables!`). The tables stored slots in
  `grid_disk_neighbours` order (including cells outside the mesh), while the
  solver reads slot *k* as the *k*-th entry of the cell's stored neighbour
  list. Rebuilding the build-time order with pya5 0.7.3 showed the two differ
  for 90–96% of cells (Carlisle res-16: 1,757 of 1,958; planar res-18: 15,879
  of 16,567). Effect: wrong sills acted as phantom dams — on the synthetic
  notch DEM the upstream pond rose to 2.4 m, nearly 1 m above the 1.505 m
  embankment crest, before the fix, and fills to the notch and spills after
  it. Affects every SGS result to date.
  - Slots now follow the neighbour list, for any valence; tables carry
    `static_vars["sgs_layout"] = 2`. `initialise_flow_model` refuses older SGS
    tables (`FLOODA5_ALLOW_OLD_SGS_LAYOUT=1` overrides).
  - `mesh/repair_sgs_slot_order.py` repairs existing meshes without the DEM
    (needs pya5 0.x and the original 0.x mesh, e.g. `*_a5v0.parquet`); the
    repaired Carlisle res-16 tables are bit-identical to a rebuild with the
    fixed builder.
- `test/synthetic_dem/test_sgs_synthetic.jl` now closes the domain boundaries,
  as its assertions (and the T4 mass balance, which ignores outflow) assume;
  it predates the open-boundary default. T0–T4 pass with the fixed tables.

#### Added
- SGS on meshes with more than five neighbours (`n_slots` from the mesh).
  `test/test_sgs_general.jl` (in CI): SGS tables on a hexagonal mesh — sill
  symmetry across every edge, curves within edge lengths, lake at rest, mass
  balance, refusal of the old layout.

### WP2: mesh-general momentum, limiters and CFL (`wp2-momentum-limiters`)

Work package WP2 from `FloodA5_Review_and_ForwardPlan_2026-10.md` (§4.2, §9).
All options are opt-in; the defaults are bit-identical to before.

#### Added
- `--momentum-model blend`: q_prev = θ·q_face + (1−θ)·(face-normal projection
  of the averaged adjacent cell vectors), θ = `--q-centre-theta`. Replaces
  Q-centring's "most collinear face" average, which on pentagons mixes faces
  about 144° apart and by itself deflects the steady flow (WP1 finding).
  θ = 1 reproduces edge momentum without Q-centring and θ = 0 reproduces
  `cell` exactly (tested). Standard solver only.
- `--limiter positivity`: one positivity-preserving scale factor per donor
  cell, s = min(1, V / Σ outflow), applied to all of its outflows and written
  back to the stored momentum (`_positivity_scale!`), in place of the per-edge
  volume limiter (Fix B) and the per-edge donor cap. Mass-conservative and
  isotropic. Standard and SGS solvers. **Not stable on its own** on the
  Carlisle DEM (velocity noise up to 5 m/s): the per-edge volume limiter is
  doing real stabilising work.
- `--limiter hybrid`: keeps the per-edge volume limiter and uses the
  positivity scaling only in place of the per-edge donor cap. As smooth as
  `legacy` on Carlisle res-18 and removes the face-count-dependent cap.
- `--face-flux-method diamond-wetdry`: as `diamond`, but the two-point face
  gradient is used where an edge's vertex stencil includes a dry cell (whose
  bed elevation would otherwise be read as a water level). Bit-identical to
  `diamond` when every cell is wet. `diamond_stencil_wet` in `DiamondFlux.jl`.
- `--cfl-mode local`: dt = 0.7 · min over wet cells of ℓ_i/√(g h_i), with
  ℓ_i = 4A_i/P_i (`FlowState.cell_length`; equal to dx on a square cell),
  instead of √(min area) paired with the 99th-percentile depth. The SGS
  terrain-range depth cap still applies.
- `FlowState` fields 34–36 (`limiter`, `cfl_mode`, `cell_length`), with a
  backward-compatible 33-argument constructor.
- `test/test_wp2_momentum_limiters.jl` (in CI); a WP2 case in
  `test/check_determinism.jl`.

#### Fixed
- ⚠ **Phantom ghost edges on Cartesian meshes** (`_build_ghost_edges`). Shared
  polygon sides were matched with a key ordered by longitude only, so a
  north–south side (lon1 == lon2) could keep a different vertex order in the
  two polygons and fail to match, becoming an extra "missing" side. On the
  dx = 25 m square mesh this created 423 ghost edges instead of 244 (open
  boundary width 10.6 km instead of 6.1 km): open-boundary Cartesian runs had
  too much outflow width. Sides are now matched in both vertex orders.
  A5 meshes are unaffected (Carlisle res-18 and planar res-16 ghost edges
  identical before and after). Regression test G1.

### WP1: vector friction (`wp1-vector-friction`)

Work package WP1 from `FloodA5_Review_and_ForwardPlan_2026-10.md` (§3, §9).
Opt-in; the default (`component`) is bit-identical to the previous behaviour.

#### Added
- `--friction-model component|vector` (`FlowState.friction_model`, the 33rd
  field; backward-compatible 31- and 32-argument constructors default it to
  `:component`). `vector` uses the full discharge magnitude
  |q| = √(q_n² + q_t²) in the friction term of every Bates kernel and the SGS
  R–A kernels, with q_t the face-tangential component of the averaged
  adjacent cell vectors. The Froude limiter then scales q_n by
  min(1, q_max/|q|). Removes the orientation-dependent face conductance of
  Bates/LISFLOOD-FP component friction (2^¼ diagonal excess and up to 12.8°
  flow deflection on a Cartesian grid, independent of resolution).
- Phase F (cell-vector reconstruction) factored into `_reconstruct_qvec!`
  (threaded, deterministic) and run every step under vector friction, for
  both `step_standard!` and `step_sgs!`.
- `test/test_wp1_vector_friction.jl` (in CI): kernel tests, lake at rest, and
  a frozen-coefficient friction-limit isotropy test through `step_standard!`
  that reproduces the component-friction theory exactly (12.78°, 2^¼) and
  shows none under vector friction.
- `test/frozen_isotropy_lib.jl`, `test/wp1_frozen_isotropy_mesh.jl` (the same
  test on any GeoParquet mesh), `test/wp1_metrics.jl` (§3.4 metrics:
  harmonic anisotropy, sector radii, planar N/S asymmetry) and
  `test/run_wp1_matrix.sh` (the WP1 run matrix).
- `test/check_determinism.jl`: two vector-friction cases.

#### Known limitations
- Ghost (open-boundary) edges still use component friction.

### WP3 (part 1): mesh families and raster-aligned Cartesian grids (`wp3-mesh-families`)

Work package WP3 from `FloodA5_Review_and_ForwardPlan_2026-10.md` (§6.1, §9),
the parts that do not touch the solver: SGS generalisation (per-edge-side
curves) is left until WP1/WP2 are merged.

#### Added
- `mesh/mesh_families.py`: FloodA5 GeoParquet meshes of other families —
  regular hexagons, rotated Cartesian, Voronoi (jittered, optionally Lloyd-
  relaxed to centroidal; generator or centroid centres) — and an importer for
  any polygon GeoParquet/GeoJSON/GeoPackage (reprojected to EPSG:4326).
  Edge-based adjacency; ID scheme `polygon`.
- `generate_cartesian_mesh_raster` (`mesh/CartesianGrid.jl`) and
  `--meshgen-cartesian FILE.tif|.asc`: a Cartesian grid aligned cell-for-cell
  with a raster in its own (projected) CRS, with the raster as the bed
  elevation — for LISFLOOD-FP comparisons. Note: FloodA5 measures cell areas
  on its spherical metric, 0.38% below the British National Grid pixel area
  near Carlisle (uniform).
- `test/test_mesh_families.jl` (in CI): Tier 0 on every family — closure
  Σ w·n̂ = 0, valence, orthogonality where guaranteed, lake at rest under the
  legacy, WLSQ, diamond and cell-momentum options, exact mass balance, import
  round trip, raster-aligned grid.

#### Fixed
- `a5_bridge.mesh_to_geoparquet` padded polygons with fewer vertices than the
  longest with (0, 0) vertices (a spike to the origin) and overwrote an
  unclosed ring's last vertex. Only meshes with mixed vertex counts were
  affected (not A5 or Cartesian meshes).
- Domain-boundary detection (`_build_ghost_edges`): a cell is on the boundary
  if its polygon has more sides than neighbours. It was "fewer neighbours
  than the mesh maximum", which on a mixed-valence mesh gave a max-valence
  edge cell no ghost edge. Identical on A5 and Cartesian meshes (ghost
  edges unchanged on Carlisle res-18 and planar res-16).
- Degenerate-face guards in `A5Grid._edge_geometry`: 1 mm instead of 1 m.
  A 1 m face is normal on a Voronoi mesh (and on A5 at resolution ≥ 22); such
  faces got a zero normal. A5 edge geometry is bit-identical on the Carlisle
  and planar meshes.

### WP-M1: pya5 1.x support (`maint-pya5-1x`)

Work package WP-M1 from `FloodA5_Review_and_ForwardPlan_2026-10.md` (§9).
pya5 0.7.3 behaviour is unchanged; pya5 1.x is now supported.

#### Added
- One code path per pya5 major version in `mesh/a5_bridge.py` (`A5_MAJOR`,
  `_uncompact_to`, `_grid_disk_ring1`) and `mesh/A5Grid.jl`
  (`pya5_version()`, `a5_id_scheme()`, version-aware `grid_disk_neighbours`).
  Under 1.x the AOI seed set comes from `polygon_to_cells` (overlapping
  containment): Carlisle res-18 mesh generation 1 m 02 s vs 13 m 13 s.
- Cell-ID scheme recorded in every GeoParquet mesh (schema metadata key
  `flooda5`: `a5_id_scheme`, `pya5_version`) and in the HDF5 `/mesh`
  `id_scheme` attribute; new `A5Mesh.id_scheme` field (7th; the 6-argument
  constructor defaults it to the installed pya5's scheme; Cartesian meshes
  use `"cartesian"`).
- ⚠ `load_mesh_geoparquet` refuses a mesh whose A5 ID scheme differs from the
  installed pya5's (`check_id_scheme`; override with
  `FLOODA5_ALLOW_ID_SCHEME_MISMATCH=1`). Existing meshes are `a5-v0`, so under
  pya5 1.x they must be migrated first.
- `mesh/migrate_mesh_ids.py`: rewrites `cell_id`/`neighbours` with
  `a5.migrate`, verifying every centre. Verified on Carlisle res-18: a
  simulation on the migrated mesh is bit-identical to the original, and a
  mesh generated afresh under 1.x has the same cells and adjacency graph
  (centroids within 1e-13°, areas within 6e-10 relative).
- `A5Grid.migrate_cell_id`; `test/test_pya5_compat.jl` (in CI) with the fixture
  `test/fixtures/a5_id_schemes.csv`; CI matrix over pya5 0.7.3 and 1.0.0
  (pya5 removed from `ci/python-requirements.txt`); CI also runs on `maint-*`
  branches.

### WP0: correctness and test infrastructure (`wp0-correctness`)

Correctness and test-infrastructure work package WP0 from
`FloodA5_Review_and_ForwardPlan_2026-10.md` (§4.4, §9). **Some fixes change
numerical results** (marked ⚠); re-baseline any stored comparison outputs.

#### Fixed
- ⚠ **Data race in Q-centred momentum** (`step_standard!`, `step_sgs!`): the
  threaded edge loop read neighbouring edges' momentum from the same array it
  was overwriting. Q-centring now reads a start-of-step snapshot, giving a
  proper Jacobi update that is deterministic across thread counts. Affects
  every run with `--momentum-model edge` and `--q-centre-theta` < 1 (the
  defaults), standard and SGS.
- ⚠ **SGS CFL depth cap never applied** (`_cfl_dt`): `method == SGSFlow`
  compared an instance with a type (always false); now `method isa SGSFlow`.
  SGS runs with overfull channel cells will take larger, correctly capped
  timesteps.
- **`_adjacency_slot`** searched only 5 slots by default and returned slot 1
  silently when a neighbour was not found; it now searches the mesh's own
  neighbour ceiling, `step_sgs!` passes it explicitly, and
  `initialise_flow_model` verifies SGS slot consistency at start-up (and
  errors clearly if a mesh has more neighbours than the SGS tables support).
- **Diamond degenerate vertex** (`eta_at_vertex`) returned 0.0 m instead of
  NaN, defeating the documented per-edge fallback to the legacy gradient.
- **O(n²) allocation in SGS initialisation** (`A5Grid.sgs_table`): the
  `get(dict, key, fill(NaN, n))` default was evaluated eagerly for every cell.
- ⚠ (diagnostic output only) **Cell velocity** (`_compute_velocity!`) is now the
  depth-averaged velocity in m/s from the WLSQ cell-vector reconstruction.
  The previous scatter used degree-space directions with no cos(lat) factor
  and had units of s⁻¹. HDF5 `velocity`/`vel_u`/`vel_v` magnitudes change.
- Startup warning and docs wrongly said `--gradient-correction` has no effect
  under `--flow-model sgs`; it does (legacy WLSQ + V̂, unvalidated). New
  warning for `--face-flux-method diamond` under SGS.

#### Added
- `FLOODA5_NO_VIS=1` environment variable: stubs the visualisation modules so
  FloodModel.jl loads without GLMakie/Oxygen (headless, CI, AI-assisted runs).
- `test/test_wp0_regressions.jl`: regression tests for the fixes above, plus a
  lake-at-rest (well-balancing) test across all flux/momentum options.
- `test/check_determinism.jl`: thread-count determinism fingerprint.
- Continuous integration: `.github/workflows/ci.yml`, minimal `ci/Project.toml`
  environment, `ci/run_tests.sh`, `ci/python-requirements.txt`.

### Cartesian formulation (`cartesian-formulation`)

Full account: `HYDRAULICS.md` §14 and the project's internal development
history (available on request).

#### Added
- **Cartesian mesh backend** (`mesh/CartesianGrid.jl`): a rectangular,
  4-connected grid generator, pure Julia (no Python bridge for geometry
  generation). New CLI flags `--meshgen-cartesian FILE`, `--dx METRES`,
  `--cartesian-lat-ref DEG`. Covers the AOI's bounding box (not clipped to
  its polygon shape). Returns the existing `A5Mesh`/`A5Cell` types, so
  everything downstream of mesh acquisition is unchanged and
  backend-agnostic.
- Diagnostic scripts: `test/audit_cartesian_geometry.jl`,
  `test/test_octagonal_anisotropy.jl`, `test/check_boundary_proximity.jl`.

#### Changed
- **`N_SIDES` generalisation.** A new `FlowState.max_neighbours` field,
  derived from the loaded mesh's own adjacency, replaces hardcoded
  neighbour-count assumptions in matrix sizing (`adj_matrix`,
  `wlsq_weights`, `mom_weights`, `cell_edge_index`), the donor-edge safety
  divisor (`_donor_edge_divisor(state)`), the Fix B volume limiter (new
  `vol_limit_divisor` keyword on `_bates_flux_limited`/
  `_bates_flux_limited_corrected`, default `5.0` = unchanged A5 behaviour),
  and the ghost-edge boundary threshold. **No behavioural change on A5** —
  confirmed by the full existing regression suite. `N_SIDES`/
  `DONOR_EDGE_DIVISOR` remain as A5-only legacy constants for existing test
  fixtures. A backward-compatible 31-argument `FlowState` constructor keeps
  every existing call site working unmodified.
- **Cell-ID normalisation hardened.** Six previously-duplicated, unsafe
  inline `A5Grid._to_hex(parse(UInt64, id, base=16))` calls (which threw on
  any non-hex ID) now route through `_norm_cell_id` / `_to_hex_safe`, which
  fall back to the raw string — so a non-A5 backend's own ID scheme no
  longer crashes mesh loading.

#### Known issues / open
- A real, small, **resolution-independent** cardinal-vs-diagonal
  directional bias was found on the (provably exactly orthogonal)
  Cartesian backend; it got *worse* under grid refinement, ruling out
  ordinary truncation error. Leading hypothesis is the per-edge momentum
  representation; the confirming test (`--momentum-model cell`) has not
  been run. See `HYDRAULICS.md` §14.
- SGS (`--flow-model sgs`) is **not** generalised: `build_sgs_tables!` still
  hardcodes 5 neighbour slots, so a Cartesian mesh is unsupported with SGS.
  There is currently no guard preventing the combination.
- Projected-CRS support does not exist: the Cartesian generator works in
  decimal degrees via a single reference latitude, so exact raster
  alignment against a LISFLOOD-FP grid in its own projected CRS is not yet
  possible.

## [0.1.0] — 2026-08-26

First formally versioned release. Merges the `directional-bias-reformulation`
branch to `main` and brings the documentation set up to date with the
current state of the model.

### Added
- Diamond face-flux reconstruction (`--face-flux-method diamond`) — a
  face-local, polygon-native alternative to the cell-averaged WLSQ gradient
  correction, proven exact for linear WSE fields.
- Cell-vector momentum representation (`--momentum-model cell`) — a
  Perot-style reconstruction replacing five independent per-edge momentum
  scalars with one coherent 2D vector per cell.
- `--gradient-correction-alpha` — scaling factor for the tangential
  correction term in the legacy (WLSQ) gradient-correction method.
- A startup warning when `--momentum-model cell`, `--gradient-correction
  on`, `--face-flux-method`, or `--gradient-correction-alpha` are combined
  with `--flow-model sgs` (the default flow model), since none of them
  currently have any effect there.
- `--version`/`-v` CLI flag; version now printed in the startup banner and
  `--help` header.
- `VERSION`, `VERSIONING.md`, and this `CHANGELOG.md`.

### Fixed
- `gradient_correction`'s keyword default in `initialise_flow_model` and
  `run_flood_model` was `true`, contradicting the CLI's own default of
  `false`. CLI runs were unaffected (the CLI always passes an explicit
  value), but several test files calling `initialise_flow_model` directly
  without specifying this argument were silently using the corrected
  kernel rather than the intended uncorrected baseline. Both defaults
  corrected to `false`.
- `--help` was missing documentation for the (fully working)
  `--momentum-model` and `--gradient-correction-alpha` flags.
- The resolution-to-cell-area guidance shown in `--help` was incorrect by
  roughly 16× at resolution 14, and disagreed with the (correct, internally
  consistent) table in `docs/A5_QUIRKS.md`. Corrected to match.

### Documentation
- `docs/HYDRAULICS.md` rewritten in full — the previous version predated
  the stability-limiter fixes, the SGS R-A flux path, and all of the
  directional-bias correction work.
- `docs/METHODS.md` §5.4, `README.md`: corrected statements that the
  non-orthogonal gradient correction was enabled by default and had
  resolved the directional bias — neither is accurate; the correction is
  opt-in and the bias is reduced, not eliminated, even when enabled.
- `docs/USER_GUIDE.md`: added a complete CLI reference section for the
  five directional-bias-correction flags, previously undocumented.
- `docs/DATA_FORMATS.md`: corrected the documented mesh adjacency schema
  (a single variable-length `neighbours` column, not fixed `adj_0`..`adj_4`
  columns), added missing `vel_u`/`vel_v` HDF5 datasets, corrected
  `flux_Q`'s presence condition, fixed an off-by-one in the SGS elevation
  bins array length, and documented that `.bci` `HFIX`/`HVAR` entries are
  parsed but not yet implemented.
- `docs/ARCHITECTURE.md`: `FlowState`/`EdgeList` struct references updated
  to match the current code (both were missing roughly twenty fields added
  since the doc was last touched); added the `mesh/DiamondFlux.jl` module;
  corrected the simulation step-function phase descriptions; updated the
  test-file inventory; clarified the mesh-generation thread-safety guidance
  (see below).
- `docs/A5_QUIRKS.md`: corrected references to three functions
  (`_edge_cos_theta`, `_build_adjacency_grid_disk`) that no longer exist
  under those names, having been consolidated into `_edge_geometry` and
  `_build_adjacency_shared_vertices`/`_build_adjacency_matrix!`.
- `docs/TEST_CASES.md`: corrected the Square Domain test description to
  match the actual working test script (resolution and workflow had
  drifted from what's implemented) and removed references to two files
  that don't exist in the repository; added full documentation of the
  `test/planar_embankment/` test case, previously undocumented despite
  being the primary acceptance test for the directional-bias correction
  work.
- Mesh-generation thread-safety guidance (`--threads 1`) reworded across
  all docs to be evidence-based: two specific historical crash causes are
  confirmed fixed, but the same general hazard (a GDAL/PyCall coordinate
  transform called from inside a threaded loop) remains present in two
  other code paths (standard-flow DEM sampling and SGS hypsometric-curve
  building), so `--threads 1` remains the recommended default with an
  explicit fallback note rather than being removed.
- Added a "Development Roadmap" section and an "Acknowledgements" section
  to `README.md`; removed direct references to internal-only planning and
  handover documents from the public documentation set, replaced with
  general pointers to the project's internal development history
  (retained, not published, available on request).

[0.1.0]: https://github.com/ (tag pending)
