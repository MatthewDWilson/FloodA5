# FloodA5 — Methods

This document describes the numerical methods, grid geometry, and physical
assumptions underlying FloodA5. It is intended for readers who want to understand
what the model does and why, rather than how to run it (see [USER_GUIDE.md](USER_GUIDE.md)
for that).

---

## 1. The A5 Pentagonal Grid

### 1.1 What is a DGGS?

A Discrete Global Grid System (DGGS) tessellates the surface of the Earth into
a hierarchy of cells. Unlike a projected raster, a DGGS operates natively in
geographic coordinates; cells are defined on the sphere and their geometry is
computed geodetically. FloodA5 uses the **A5 DGGS**, in which every cell is a
pentagon.

### 1.2 Key properties

**Equal area.** All cells at a given resolution level have the same plan area,
regardless of latitude. This is a fundamental property of the A5 design and means
that areal quantities (rainfall volume, storage capacity) can be computed without
any map-projection correction.

**Uniform five-connectivity.** Every interior cell has exactly five edge-sharing
neighbours. This eliminates the directional bias inherent in rectangular grids, where
north–south and east–west directions are privileged. On a square grid a circular flood
front becomes a diamond; on the A5 grid it remains circular.

**Hierarchical multi-resolution.** Each cell at resolution level L contains exactly
five children at level L+1. This strict parent–child relationship supports future
adaptive mesh refinement without re-meshing.

**Dual sublattice.** At any given resolution level, cells tile in two interleaved
sublattices (identified by the leading nibble of the 64-bit cell ID). Both sublattices
must be included in any mesh. FloodA5's mesh generation handles this automatically.

### 1.3 Cell IDs

Every A5 cell has a unique 64-bit integer identifier, serialised as a 16-character
zero-padded hexadecimal string (e.g. `6345f2518e800000`). The ID encodes both the
resolution level and the cell's hierarchical position. All cell IDs in FloodA5 are
stored and compared in this 16-character form.

### 1.4 Resolution levels

The relationship between resolution level and cell size is approximately:

| Level | Approx. cell area | Cell spacing |
|-------|-------------------|--------------|
| 10 | ~50 km² | ~7 km |
| 12 | ~2 km² | ~1.4 km |
| 14 | ~13 ha | ~355 m |
| 16 | ~8,000 m² | ~89 m |
| 18 | ~500 m² | ~22 m |
| 20 | ~31 m² | ~5.6 m |

Cell areas are exactly equal for all cells at a given level (the equal-area property).
The values above are representative; use `_polygon_area_m2` for the exact geodetic
area of a specific cell computed from its boundary polygon.

---

## 2. Mesh Generation

Given a GeoJSON area-of-interest (AOI) polygon and a target resolution level,
the mesh is built in three steps:

1. **Cell enumeration** — the A5 library (`pya5`) identifies all cells at the target
   resolution whose centres fall inside the AOI, covering both sublattices.
2. **Boundary computation** — the pentagon boundary vertices (five vertices + closing
   repeat) for each cell are computed geodetically.
3. **Adjacency construction** — for each pair of cells, the shared edge is identified
   by matching boundary vertices. Only edges where both cells are in the mesh are
   retained. The result is an undirected `EdgeList` where each edge is stored once,
   with the lower cell index as `cell_i`.

The Python bridge (`mesh/a5_bridge.py`) handles all `pya5` calls and writes the
mesh to GeoParquet. The Julia side calls the bridge as a subprocess, avoiding
thread-safety issues with PyCall on non-main threads.

Edge geometry is computed from the actual cell boundary polygons:

- **Edge width** (`w`) — haversine length of the shared edge arc.
- **Centre-to-centre distance** (`L`) — haversine distance between the two cell
  centres.
- **Non-orthogonality angle** (`cos θ`) — dot product of the centre-to-centre
  unit vector **d̂** with the edge face normal **n̂**, computed on a local
  equirectangular projection. On A5 pentagons θ ranges from ~16° to ~38° (mean ~23°).
- **WLSQ correction vector** (`V̂ = n̂ − (d̂·n̂)·d̂`) — the tangential component of
  the face normal after removing its projection onto **d̂**. This encodes both the
  direction and magnitude of the non-orthogonality and is used by the gradient
  correction described in §5.3. By construction `|V̂| = sin θ`, so the correction is
  zero for orthogonal edges and increases smoothly with skew angle.

---

## 3. DEM Ingestion

FloodA5 samples a user-supplied GeoTIFF DEM onto the mesh at initialisation time.
Two sampling methods are available:

**Mean sampling** (`--dem-method mean`, default) — distributes 256–512 quasi-random
Halton points within each cell polygon and takes the arithmetic mean of the DEM
values at those points. Halton sequences provide better coverage of the cell area than
pseudo-random points at the same count, reducing sensitivity to the placement of
individual sample points. This method gives a better estimate of the true mean cell
elevation and is recommended for all production runs.

**Centroid sampling** (`--dem-method centroid`) — a single bilinear sample at the
cell centre. Fast, but sensitive to the location of the centre point within the
terrain. Use only for quick exploratory runs.

In both cases, DEM values are sampled using bilinear interpolation after reprojecting
the query coordinates from WGS 84 into the DEM's native CRS. Out-of-bounds cells
receive `NaN` elevation (a warning is logged); they are hydraulically inert.

Sampled elevations are stored as a static column in the GeoParquet mesh file so that
subsequent runs with `--meshload` do not require re-sampling.

---

## 4. Sub-Grid Sampling (SGS)

### 4.1 Motivation

At resolution level 14, each pentagon covers approximately 13 ha. A typical 1 m
LiDAR DEM contains ~130,000 pixels in that area. Using a single representative
elevation per cell discards virtually all of this information: channels, roads,
embankments, and building footprints — the features that control urban flood routing
— are invisible to the model.

SGS addresses this by replacing the single cell elevation with a **hypsometric curve**:
a pre-computed lookup table relating water surface elevation (WSE) to stored volume
and wetted plan area. The model can then resolve partial wetting of a cell and route
water through sub-cell topographic lows even when those features are narrower than
the cell diameter.

### 4.2 Pre-processing

For each cell, SGS pre-processing:

1. Distributes 512 quasi-random Halton points within the cell polygon and samples
   the DEM at each point using bilinear interpolation.
2. Sorts the resulting elevation distribution and bins it into 100 quantile-spaced
   levels from `z_min` to `z_max`.
3. At each elevation level, computes the cumulative stored volume (treating the DEM
   surface as the bed and the current elevation as the water surface) and the wetted
   plan area.
4. Stores the resulting `(elevation, volume, area)` curves as Arrow list columns in
   the GeoParquet mesh file.

Along each shared cell boundary, the same sampling process builds cross-sectional
hypsometric curves: the sorted elevation profile of the edge arc gives cross-sectional
flow area `A(wse)` and wetted perimeter `P(wse)` at each of the 100 knot elevations.
The minimum elevation along each edge arc is also stored as the **edge sill** — the
threshold WSE below which no flow can cross that edge.

Pre-processing is automatic when `--flow-model sgs` is used with `--meshgen`, and is
stored in the mesh parquet for reuse. It is the most computationally expensive
initialisation step (~5 s at resolution 14 for a 2,900-cell domain).

### 4.3 Volume-to-WSE lookup

At each timestep, the stored volume in each cell is converted to WSE via inverse
interpolation of the volume curve (`wse_from_volume`). When a cell's volume exceeds
the top of its terrain range — i.e., all terrain is submerged — the WSE is
extrapolated linearly above `z_max`:

```
WSE = z_max + (V - V_full) / cell_area
```

This treats the cell as a vertical-walled container above its maximum terrain
elevation, preserving a positive driving head for inter-cell flux in fully-submerged
cells.

---

## 5. Flow Routing

### 5.1 Bates et al. (2010) inertial formulation

Both the standard and SGS solvers route water using the inertial shallow-water
formulation of Bates, Horritt & Fewtrell (2010), which approximates the full
Saint-Venant equations by retaining the local acceleration term while neglecting
advection. The unit discharge `q` (m²/s) across an edge is updated each timestep by:

```
q^t = [ q^{t−Δt} − g · h_flow · Δt · (ΔWSE / L_eff) ]
      / [ 1 + g · h_flow · Δt · n² · |q^{t−Δt}| / h_flow^(10/3) ]
```

where:
- `h_flow = max(WSE_i, WSE_j) − z_sill` — flow depth above the edge sill on the
  higher side
- `ΔWSE = WSE_i − WSE_j` — water surface slope (positive when cell i is higher)
- `L_eff = L × cos θ` — effective centre-to-centre distance, corrected for
  non-orthogonality
- `n` — Manning's roughness coefficient
- `z_sill` — edge sill elevation (standard: `max(elev_i, elev_j)`; SGS:
  minimum DEM elevation along the shared boundary)

The volumetric flux is `Q = q × w` (m³/s), where `w` is the edge width.

**Sign convention:** `Q < 0` means flow from cell i to cell j (i is higher);
`Q > 0` means flow from cell j to cell i. Volume changes are applied as
`ΔV[i] += Q·Δt` and `ΔV[j] -= Q·Δt`.

### 5.2 SGS R-A flux kernel

For the SGS solver, the wide-channel approximation in the denominator
(`h_flow^(10/3)`) is replaced by the Manning hydraulic radius formulation used
in LISFLOOD-FP's sub-grid channel (SGC) solver:

```
Q^t = [ Q^{t−Δt} − g · A · Δt · (ΔWSE / L_eff) ]
      / [ 1 + g · Δt · n² · |Q^{t−Δt}| / (R^(4/3) · A) ]
```

where `A` (m²) is the cross-sectional flow area and `R = A/P` (m) is the hydraulic
radius, both derived from the pre-computed edge hypsometric curves at the current WSE.

As the channel fills, `A` grows and `R` grows, making the friction term
`n²|Q| / (R^(4/3)·A)` self-stabilising. Large fluxes are naturally damped without
requiring an explicit volume limiter, and the formulation gives physically correct
friction scaling for concentrated channel flow at fine resolutions (level 18,
`dx ≈ 22 m`) where the wide-channel approximation breaks down.

`A` and `R` are averaged from both sides of each edge (`A_edge = 0.5·(A_i + A_j)`,
`R_edge = 0.5·(R_i + R_j)`) so that the flux kernel is symmetric.

### 5.3 Stability fixes

The following stability measures are applied to both solvers, informed by
comparison with the LISFLOOD-FP ACC and CAESAR-Lisflood implementations:

**Wave-speed CFL timestep.** The adaptive timestep is computed as:

```
Δt = 0.7 × dx_min / √(g × h_max)
```

where `dx_min = √(min cell area)` and `h_max` is the 99th-percentile water depth
across wet cells (to prevent a single overfull cell from collapsing the timestep
for the whole domain). The Courant number 0.7 matches the LISFLOOD-FP default.

**Dry-edge threshold.** Edges with `h_flow ≤ 0.001 m` are treated as dry and
flux is set to zero. This prevents stale momentum from carrying through thin films.

**Froude limiter.** Unit discharge is capped at `q_max = h_flow × √(g·h_flow) × 0.8`
(Froude number ≤ 0.8) for the standard solver. This suppresses supercritical
oscillation on the pentagonal mesh, where five independently-evolving momentum
states per cell provide no directional damping. The SGS R-A kernel is
self-stabilising and does not require this cap.

**Q-centred momentum smoothing.** The inertial term uses a spatially smoothed
unit discharge:

```
q_eff = θ × q_prev + (1 − θ) × mean(q_collinear)
```

where `q_collinear` is the mean of the fluxes on the most collinear edges of the
two adjacent cells (one on each side, averaged over however many are available —
boundary cells may have only one). θ = 0.9 matches the LISFLOOD-FP default and
provides additional damping of the checkerboard instability mode. The checkerboard
oscillation has opposite signs on alternating edges, so averaging with the collinear
neighbours cancels it; for a coherent flow field the smoothing has negligible effect.

**Consistent momentum state.** The stored `q_prev` value for the next step is always
the post-limiting unit discharge, not the raw Bates value. This ensures the momentum
state is consistent with the volume actually transferred, resolving the primary driver
of checkerboarding identified in comparison with LISFLOOD-FP.

**Volume limiter.** As a last-resort mass-conservation guard, volume transfer via
any single edge is capped at `V_donor / 10` per step, ensuring at most 50% of a
cell's volume can leave via all five edges in a single timestep.

---

### 5.4 Non-orthogonal gradient correction

A5 pentagon edges are generally not perpendicular to the centre-to-centre vector
**d**. Without correction, the face-normal pressure gradient — the term driving
flux in eq. (1) — is approximated along **d** rather than along the true face
normal **n̂**. On A5 cells, where the angle between **d** and **n̂** (the
non-orthogonality angle θ) ranges from ~16° to ~38° (mean ~23°), this produces
a systematic directional bias: point-source flood fronts are elongated rather
than circular, and flow on a planar slope deviates significantly from the
analytically expected downslope direction.

FloodA5 corrects this using **weighted least-squares (WLSQ) cell-centre gradient
reconstruction**, following the approach of Jasak (1996) and Moukalled et al.
(2016) for non-orthogonal finite-volume meshes.

#### Gradient reconstruction

At each timestep, before the flux loop, a 2D WSE gradient vector
∇WSE_i = (∂WSE/∂x, ∂WSE/∂y) is reconstructed for every cell _i_ by solving:

```
minimise Σ_k w_k [(WSE_jk − WSE_i) − ∇WSE_i · (x_jk − x_i)]²
```

over all mesh neighbours j₁, …, j₅. The inverse-distance weights `w_k = 1/|x_k − x_i|²`
down-weight distant neighbours. The solution (a 2×2 weighted normal equation system)
depends only on cell geometry and is pre-computed once at model initialisation as a
`(2 × 5)` projection matrix per cell (`wlsq_weights`). The per-timestep gradient
computation is then a single matrix–vector multiply per cell — O(5 × n_cells), negligible
compared to the flux loop.

#### Corrected driving head

At each edge `e = (i, j)`, the standard `ΔWSE / L` gradient term is replaced by the
non-orthogonal corrected value:

```
dWSE_n = c·(WSE_i − WSE_j) − L·(∇WSE_f · V̂)
```

where:
- `c = d̂·n̂ = cos θ` (always ≥ 0 after orienting n̂ toward j)
- `∇WSE_f = ½(∇WSE_i + ∇WSE_j)` (face-centre gradient by linear interpolation)
- `V̂ = n̂ − c·d̂` (the tangential component of the face normal — the skewness
  direction; pre-computed per edge as `skew_x`, `skew_y` in the `EdgeList`)
- `L` is the centre-to-centre distance

The first term `c·(WSE_i − WSE_j)` is the projection of the raw WSE difference
onto the face normal; the second term `L·(∇WSE_f · V̂)` corrects for the lateral
offset between the face midpoint and the point where **d** crosses the face. Together
they reconstruct the face-normal gradient without the directional bias of the
uncorrected scheme. For orthogonal edges, `V̂ = 0` and `c = 1`, so the formula
reduces exactly to the standard `WSE_i − WSE_j`, confirming backward compatibility.

This replaces `dWSE / L_eff` in eq. (1) and eq. (2), and `L` is used as-is (no
`cos θ` scaling). The same `_bates_flux_corrected` and `_manning_flux_ra_corrected`
kernels are used for standard and SGS solvers respectively.

#### Validated improvement

On A5 meshes at resolution 16 (mean θ ≈ 23°, max θ ≈ 38°):

| Test | Metric | Before correction | After correction |
|---|---|---|---|
| Point-spread (flat domain, point source) | Polsby–Popper circularity | 0.938 | **0.980** |
| Point-spread | Directional front-radius CV | 0.137 | **0.055** |
| Planar slope (0.1% E–W slope) | N/S volume asymmetry at 2 h | 0.254 | **0.066** |

The non-orthogonality range on A5 meshes (max 38°) is well within the OpenFOAM
"corrected scheme safe" envelope (≤70°), so no additional limiting of the correction
term is required.

#### Computational cost

The correction adds one O(5 × n_cells) gradient reconstruction pass per timestep — a single multiply-accumulate loop over the neighbour structure, cheaper than the edge flux loop (O(2.4 × n_cells) edges but with heavier per-edge arithmetic). In practice the overhead is modest and solver-dependent:

| Solver | Mesh | Overhead vs uncorrected |
|---|---|---|
| Standard | res-16, 1,958 cells | +15% simulation wall time |
| Standard | res-18, 29,902 cells | +30% simulation wall time |
| SGS | res-16, 1,958 cells | +11% simulation wall time |
| SGS | res-18, 29,902 cells | +4% simulation wall time |

The larger fractional overhead for the standard solver at res-18 reflects the fact that the corrected run takes slightly more steps (the correction shifts the effective CFL slightly), not a per-step cost difference. For the SGS solver the hypsometric table lookups dominate per-step cost, making the gradient pass negligible at both resolutions. All figures from Carlisle domain simulations (10 h, 50 mm/hr), validated 2026-07-02.

#### CLI control

Gradient correction is enabled by default. It can be disabled for benchmarking
or comparison via `--gradient-correction off`. The `--q-centre-theta` flag (default
0.9) adjusts the Q-centred momentum smoothing weight (set to 1.0 to disable
momentum smoothing entirely while retaining gradient correction).

---

### 5.5 Manning's roughness

Manning's `n` per edge is the arithmetic mean of the two adjacent cell values:
`n_edge = 0.5 × (n_i + n_j)`. This matches the LISFLOOD-FP convention and
standard shallow-water modelling practice.

A global value is set via `--manning-n` (default: 0.03 s/m^(1/3)). A per-cell
friction raster (GeoTIFF) can be supplied via `--friction`, which overrides the
global value at cells where the raster is finite.

---

## 6. Boundary Conditions

### 6.1 Domain boundary detection

Cells with fewer than five edge-sharing neighbours in the mesh are domain-edge cells.
Their missing edges face the exterior of the domain. FloodA5 represents these as
**ghost edges** — virtual edges with pre-computed geometry (width, centre-to-centre
distance, sill elevation) — and applies a BC flux kernel across each ghost edge at
each timestep.

Ghost edge widths are computed from the actual cell polygon vertices, not the mean
edge width, to avoid systematic under- or over-estimation at irregular boundaries.

### 6.2 Outflow boundary condition types

| Type | Description | Default? |
|---|---|---|
| `ZeroGradient` | Ghost-cell WSE = boundary cell WSE; transmissive/non-reflective | ✓ |
| `Closed` | No flux; equivalent to a solid wall | — |
| `Critical` | Ghost-cell WSE = sill + (2/3)·(WSE − sill); free outfall | — |

The zero-gradient condition sets `ΔWSE = 0` across the ghost edge, so the momentum
term carries existing flow off the edge without reflection and without any additional
pressure gradient driving flow out. This is the standard transmissive outflow
condition used in most flood models.

`--closed-boundaries` sets all boundary cells to `Closed`. This is useful for
mass-balance benchmarking and catchment studies where no water should leave the domain.

A GeoJSON file (`--bc-file`) can assign different BC types to specific boundary
segments, for example to close a levee while leaving a river channel open.

### 6.3 Inflow boundary conditions

**Uniform rainfall** (`--rainfall`) — applies a constant volumetric flux to every
cell each timestep: `ΔV = rate × Δt × cell_area`.

**Localised rainfall** (`--rainpoint`) — same as uniform rainfall but applied only
to the single nearest cell. Useful for point-source testing.

**Constant injection** (`--injection-point`) — a fixed volumetric flow rate (m³/s)
injected into the nearest cell. Represents a pipe outfall, pump, or other
constant-rate source.

**Time-varying hydrograph** (`--inflow-point`, `--inflow-bci`) — a volumetric flow
rate that varies in time according to a user-supplied hydrograph. The rate is linearly
interpolated between knot points, with flat extrapolation beyond the first and last
knot. Hydrographs are supplied as two-column CSV files or LISFLOOD-FP `.bdy` files.

Multiple sources of different types can be combined in a single run. Sources targeting
the same cell are summed.

### 6.4 LISFLOOD-FP compatibility

FloodA5 reads LISFLOOD-FP `.bci` and `.bdy` boundary condition files for direct
portability with existing LISFLOOD-FP workflows. Supported `.bci` entry types:

| Code | Meaning |
|---|---|
| `P QVAR` | Point source with time-varying discharge from a `.bdy` series |
| `P QFIX` | Point source with constant discharge |
| `P FREE` | Open outflow at that boundary cell |

Cardinal-direction edge entries (`N`, `E`, `S`, `W`) are not supported because
FloodA5 domains are arbitrary polygons with no axis-aligned boundaries. When
encountered, a helpful message is logged directing the user to use `--bc-file`
(GeoJSON) instead.

### 6.5 Mass balance accounting

Outflow through ghost edges is accumulated in `state.vol_removed`. The mass balance
at any time is:

```
input_vol − (domain_vol + vol_removed) ≈ 0
```

where `input_vol` is the cumulative volume from all sources and `domain_vol` is the
current total stored volume. This quantity is logged every 50 steps as `mb_err`.

---

## 7. Adaptive Timestep

The simulation loop maintains an adaptive timestep bounded above by `--dt-max`.
Each step:

1. Compute `Δt_cfl = 0.7 × dx_min / √(g × h_99)` where `h_99` is the 99th-percentile
   wet-cell water depth.
2. Set `Δt = min(Δt_cfl, dt_max, t_end − t)` to avoid overshooting the end time.
3. If all cells are dry, fall back to `Δt = dt_max` (the CFL denominator would be zero).

The final partial step is sized exactly to reach the requested simulation end time.

---

## 8. Coordinate System

All coordinates in FloodA5 are in **EPSG:4326** (WGS 84 geographic, longitude/latitude
in decimal degrees). There is no CRS reprojection at the grid level. Distance and
area computations use geodetic formulas:

- Haversine formula for great-circle distances (edge widths, centre-to-centre lengths)
- Shoelace formula on a local equirectangular projection centred on the polygon
  centroid for cell areas (< 0.1% error at resolution 14)
- Local equirectangular projection centred on the edge midpoint for non-orthogonality
  correction (< 0.05% error at resolution 14)

DEM ingestion reprojects sample points from WGS 84 into the DEM's native CRS using
ArchGDAL before sampling; the DEM raster itself is never warped.

When LISFLOOD-FP `.bci` files use projected coordinates (British National Grid,
etc.), the `--bc-epsg` flag converts them to longitude/latitude internally before
cell matching.

---

## 9. References

- Bates, P.D., Horritt, M.S., Fewtrell, T.J. (2010). A simple inertial formulation
  of the shallow water equations for efficient two-dimensional flood inundation
  modelling. *Journal of Hydrology* 387(1–2), 33–45.
  https://doi.org/10.1016/j.jhydrol.2010.03.027
- Jasak, H. (1996). Error analysis and estimation for the finite volume method
  with applications to fluid flows. PhD thesis, Imperial College London. *(WLSQ
  gradient reconstruction for non-orthogonal meshes.)*
- Moukalled, F., Mangani, L., Darwish, M. (2016). *The Finite Volume Method in
  Computational Fluid Dynamics.* Springer. Chapter 8: Gradient computation on
  unstructured meshes.
- Neal, J.C. et al. (2012). How much physical complexity is needed to model flood
  inundation? *Hydrological Processes* 26(15), 2264–2282.
- Weller, H. (2014). Non-orthogonal version of the arbitrary polygonal C-grid and
  a new diamond grid. *Geoscientific Model Development* 7, 779–797.
- Thacker, W.C. (1981). Some exact solutions to the nonlinear shallow-water wave
  equations. *Journal of Fluid Mechanics* 107, 499–508.
- A5 DGGS: [a5geo.org](https://a5geo.org)
