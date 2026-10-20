# FloodA5 — Test Cases

This document describes the test cases included in the repository. Each case
serves a specific purpose in the validation and benchmarking of FloodA5.
New cases will be added as the model is applied to additional domains.

---

## 1. Square Domain (Directional-Bias Benchmarks)

**Location:** `test/square/`
**Purpose:** Confirm that flood spreading from a point source on a flat domain
is radially symmetric — i.e., that the pentagonal grid does not introduce
directional bias — and provide the reference domain for the point-spread
directional-bias correction benchmark described in `docs/HYDRAULICS.md` §7–§8.

### Domain

A 1.35° × 1.35° square polygon centred at 0°N, 0°E (`test/square/square_domain.geojson`).
The domain straddles the null meridian and equator for a clean, unbiased test geometry.
No DEM is used — all cells have zero elevation.

Also present: `test/square/square.bci` and `test/square/lisflood_square.par`, a
LISFLOOD-FP-format configuration of the same domain for an external
rectangular-grid comparison run (not automated as part of this repository's
test suite — see "Comparison with rectangular grids" below).

### Point-spread circularity benchmark (`test/test_point_spread.jl`)

This is an **HDF5-output analysis tool** — it does not run the model itself.
Generate two comparison runs first, then point the script at both:

```bash
# 1. Generate the mesh once, reused for both runs (resolution 14)
julia --threads 1 --project=. FloodModel.jl \
    --meshgen test/square/square_domain.geojson --meshres 14 \
    --meshout test/square/square_mesh_res14.parquet \
    --flow-model standard --mesh-only

# 2. Baseline run (uncorrected, legacy kernel — the current default)
julia --threads auto --project=. FloodModel.jl \
    --meshload test/square/square_mesh_res14.parquet \
    --flow-model standard \
    --injection-point 0.0,0.0,50.0 \
    --gradient-correction off \
    --sim-duration 3600 --dt-max 30 \
    --output test/square/square_baseline.h5 --output-interval 45

# 3. Corrected run (WLSQ + skewness correction)
julia --threads auto --project=. FloodModel.jl \
    --meshload test/square/square_mesh_res14.parquet \
    --flow-model standard \
    --injection-point 0.0,0.0,50.0 \
    --gradient-correction on \
    --sim-duration 3600 --dt-max 30 \
    --output test/square/square_corrected.h5 --output-interval 45

# 4. Compare
julia --project=. test/test_point_spread.jl \
    --baseline  test/square/square_baseline.h5 \
    --corrected test/square/square_corrected.h5
```

The script computes, from the wet-cell centres at a matched output frame: a
convex hull, Polsby–Popper circularity (`4π·area/perimeter²`, 1.0 = perfect
circle), and a directional front-radius coefficient of variation (lower =
more circular spreading). See `--help` in the script for the `--frame` and
`--list-frames` options.

**Recorded result** (resolution 16, 331 cells, matched frame): circularity
0.938 → 0.980, front-radius CV 0.137 → 0.055, comparing
`--gradient-correction off` to `on`. See `docs/HYDRAULICS.md` §7–§8 for the
full picture, including that this is one of three tested corrections and none
is the current default.

### Expected results, uncorrected default configuration

- Water propagates outward from the injection cell in all five principal
  directions approximately equally, though with a measurable, non-zero
  directional bias (this is the phenomenon the benchmark above quantifies).
- Mass balance error (`mb_err` in logs) should be < 0.01% throughout.
- No negative volumes or NaN values.

### Comparison with rectangular grids

This comparison now exists *inside* the repository — FloodA5 has its own
Cartesian (4-connected, rectangular) mesh backend (`mesh/CartesianGrid.jl`,
`cartesian-formulation` branch), so the same solver can be run on both grid
types and compared directly; see §5 below for the procedure and results.
(`test/square/square.bci` and `test/square/lisflood_square.par` remain
available as a LISFLOOD-FP-format starting point for the separate, still
pending comparison against LISFLOOD-FP *itself*.)

**This section previously predicted** that a four-connected rectangular grid
would produce a *diamond-shaped* flood front — the classic hop-counting
argument that diagonal paths are ~1.41× longer than cardinal ones, so the
front should reach further along the four cardinal axes than the four
diagonals. **That prediction was tested and not borne out**: the measured
anisotropy on FloodA5's own Cartesian backend is small (~3% in front radius)
and goes the *opposite* way — the diagonal directions reach slightly
*further*, not less far. The hop-counting argument models a routing/
cellular-automaton scheme; a resolved finite-volume PDE solver exchanges
flux continuously across every edge of every wet cell simultaneously, and
diagonal cells fill via many overlapping short paths, not one zig-zag
route. See §5 for the numbers, and for the further, still-open finding that
this residual anisotropy does not shrink with grid refinement.

---

## 2. Planar Slope — Directional-Bias and SGS Test Domains

**Location:** `test/planar_embankment/`
**Purpose:** Two related but distinct uses of the same underlying terrain
generator, for two different tests.

### 2a. SGS sub-cell routing validation (`test_planar_embankment.jl`)

Tests SGS sub-cell routing against the standard solver on a planar 0.1%
west–east slope (4 km domain, 4 m total relief) with a 2 m rectangular
embankment ridge across the domain, containing a 3 m notch — narrower than a
resolution-18 A5 cell (~22 m) — at the domain centre.

```bash
# 1. Generate the DEM (once) — default embankment height 5 m
python test/planar_embankment/generate_planar_embankment_dem.py

# 2. Run tests (meshes auto-generated if absent — 2-5 min first time)
julia --threads auto test/planar_embankment/test_planar_embankment.jl

# 3. Optional: Makie visualisation (standard flow only)
julia --threads auto test/planar_embankment/test_planar_embankment.jl --vis
```

| Test | Description |
|------|-------------|
| T0 | Both meshes have cells on both sides of the embankment |
| T1 | Standard flow: no downstream water while upstream WSE < embankment crest |
| T2 | SGS: downstream water appears once WSE ≥ notch sill (well before crest) |
| T3 | SGS routes more water downstream than standard at the same injection history |
| T4 | Mass balance < 0.01% for both solvers (closed domain) |
| T5 | HDF5 output files valid: correct frame count, non-zero depths |

Full detail: `test/planar_embankment/README.md`.

### 2b. North/south volume-symmetry benchmark (`test/test_planar_symmetry.jl`)

The same slope generator with the embankment height set to zero
(`--emb-height 0`) gives a terrain-feature-free planar slope — a domain with
no physical north/south distinction, so any measured north/south volume
asymmetry is purely a mesh/solver artefact. This is the primary acceptance
test used throughout the directional-bias correction work (`docs/HYDRAULICS.md`
§7–§8): it is far less noisy than a per-cell velocity-direction metric,
because it's an integrated quantity (total volume on each side) rather than
an instantaneous per-cell one.

```bash
# 1. Generate a terrain-feature-free DEM
python test/planar_embankment/generate_planar_embankment_dem.py --emb-height 0.0

# 2. Generate the mesh, then run two comparison simulations (as in §1 above,
#    substituting --injection-point for the domain's actual source location
#    and --closed-boundaries so the metric isn't confounded by outflow)

# 3. Compare — this script is HDF5-output analysis only, like test_point_spread.jl
julia --project=. test/test_planar_symmetry.jl \
    --baseline  planar_baseline.h5 \
    --corrected planar_corrected.h5 \
    --source-lat <lat> --source-lon <lon> \
    --sweep 10
```

`--sweep N` reports the asymmetry metric at every Nth output frame across
the whole run, which is the form the result is usually reported in (a single
frame can be misleadingly noisy or, conversely, misleadingly clean — see
`docs/HYDRAULICS.md` §7 for examples of both).

**Recorded result** (resolution 16, closed boundaries): baseline `|asymmetry|`
rises to 0.25 by the end of a 2-hour run; with gradient correction, it stays
below 0.07 throughout. Best combination tested to date
(`--face-flux-method diamond --momentum-model cell`) reaches ~0.61 on a
harder, sustained-ponding resolution-18 configuration — see
`docs/HYDRAULICS.md` §8 for the full comparison table and the important
caveat that none of this fully eliminates the bias.

---

## 3. Carlisle 2005 Flood Event

**Location:** `test/carlisle/`  
**Purpose:** Validate FloodA5 against a well-documented historical flood event
using real LiDAR topography and measured river discharge hydrographs. The Carlisle
2005 event is a standard benchmark for 2D flood models.

### Domain

A 4.75 km × 3.06 km section of the River Eden floodplain, Carlisle, UK
(`test/carlisle/Carlisle_domain.geojson`). The domain covers the confluence of
the Eden, Caldew, and Petteril rivers. Approximate extents:
54.88°N–54.91°N, 2.96°W–2.89°W.

**DEM:** `test/carlisle/Carlisle_LiDAR_5m_mean.tif` — 5 m resolution mean LiDAR,
available in the repository (approximately 4.6 MB).

### Boundary conditions

Three upstream inflow points (Rivers Eden, Caldew, and Petteril) and one open
outflow boundary segment are specified in LISFLOOD-FP format:

- `test/carlisle/carlisle_wgs84.bci` — boundary condition entry file
  (WGS 84 decimal degree coordinates)
- `test/carlisle/carlisle_wgs84.bdy` — companion hydrograph file

The `.bdy` file contains 241 time steps at roughly 0.5-hour intervals covering the
flood event, with peak flows at approximately hour 36.

The `carlisle.bci` and `carlisle.bdy` variants use British National Grid (EPSG:27700)
projected coordinates and require `--bc-epsg 27700`.

### Simulation setup

**Resolution 18, standard solver:**

```bash
# Generate mesh (run once)
julia --threads 1 --project=. FloodModel.jl \
    --meshgen test/carlisle/Carlisle_domain.geojson --meshres 18 \
    --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif \
    --meshout test/carlisle/carlisle_mesh18_standard.parquet \
    --flow-model standard --mesh-only

# Simulate 5-day event
julia --threads auto --project=. FloodModel.jl \
    --meshload test/carlisle/carlisle_mesh18_standard.parquet \
    --flow-model standard \
    --inflow-bci test/carlisle/carlisle_wgs84.bci \
    --manning-n 0.03 \
    --sim-duration 432000 --dt-max 10 \
    --output test/carlisle/carlisle_standard_res18.h5 \
    --output-interval 3600
```

**Resolution 18, SGS solver:**

```bash
# Generate SGS mesh (run once — slower than standard due to SGS pre-processing)
julia --threads 1 --project=. FloodModel.jl \
    --meshgen test/carlisle/Carlisle_domain.geojson --meshres 18 \
    --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif \
    --meshout test/carlisle/carlisle_mesh18_sgs.parquet \
    --flow-model sgs --mesh-only

# Simulate 5-day event
julia --threads auto --project=. FloodModel.jl \
    --meshload test/carlisle/carlisle_mesh18_sgs.parquet \
    --flow-model sgs \
    --inflow-bci test/carlisle/carlisle_wgs84.bci \
    --manning-n 0.03 \
    --sim-duration 432000 --dt-max 10 \
    --output test/carlisle/carlisle_sgs_res18.h5 \
    --output-interval 3600
```

Using projected BCI coordinates with `--bc-epsg`:

```bash
    --inflow-bci test/carlisle/carlisle.bci --bc-epsg 27700
```

### Resolution 20

The FOSS4G 2026 paper also includes a resolution 20 (~31 m² cells) standard
solver run for comparison. Mesh generation at this resolution takes several
minutes and produces ~100,000+ cells:

```bash
julia --threads 1 --project=. FloodModel.jl \
    --meshgen test/carlisle/Carlisle_domain.geojson --meshres 20 \
    --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif \
    --meshout test/carlisle/carlisle_mesh20_standard.parquet \
    --flow-model standard --mesh-only
```

### Expected results

- Inundation of low-lying areas between the three rivers during and after peak flow.
- Water leaves the domain through the open (western) boundary; `vol_removed` grows
  monotonically after inundation onset.
- Mass balance error (`mb_err`) < 0.1% throughout.
- SGS run shows earlier and more extensive inundation compared to the standard solver
  due to sub-cell channel routing, particularly in the narrow inter-levee channels.

### Validation data

Observed flood outlines for the Carlisle 2005 event are available from the UK
Environment Agency and academic publications. Quantitative comparison (F2 score,
hit rate, false alarm ratio) against the observed outline is planned for a future
release.

### Performance notes

| Configuration | Approx. run time (reference hardware) |
|---|---|
| Res 18 standard, 5 days | ~10–15 min |
| Res 18 SGS, 5 days | ~20–30 min |
| Res 20 standard, 5 days | ~60–90 min |

Times are approximate and depend heavily on CPU core count and clock speed. Use
`julia --threads auto` to exploit all available cores. GPU acceleration does not
currently affect simulation speed (only SGS pre-processing).

---

## 4. Synthetic DEM SGS Validation

**Location:** `test/synthetic_dem/`  
**Purpose:** Verify that the SGS solver correctly detects and routes water through
a sub-cell topographic feature (a notch in an embankment) that is invisible to the
standard solver. This is the core scientific validation of the SGS approach.

### Domain

A synthetic 4 km × 2 km domain with:
- A parabolic bowl upstream basin
- A Gaussian embankment ridge (crest ~1.5 m, σ = 300 m)
- A notch in the embankment centred at mid-domain (width 300 m, sill 0.705 m)
- 88 A5 cells at resolution 14 (~12.6 ha each)

The notch is deliberately narrower than a resolution-14 cell diameter (~355 m),
making it invisible to the standard solver (which uses mean cell elevation) but
detectable by the SGS solver (which pre-computes the minimum elevation along each
shared cell boundary).

### Running the validation

```bash
julia --threads 1 --project=. test/synthetic_dem/test_sgs_synthetic.jl
```

### Tests

| Test | Description | Pass criterion |
|---|---|---|
| T0 | Mesh sanity | Both meshes load with cells on both sides of the embankment |
| T1 | No downstream flow before notch sill | 30 min at 50 mm/hr leaves downstream cells dry |
| T2 | Downstream flow after sill exceeded | Sustained injection drives flow through notch |
| T3 | SGS routes more water than standard | SGS downstream volume > standard at same injection history |
| T4 | Mass balance | < 0.01% error for both solvers over 600 steps |

### Key result

At 5,700 injection steps (identical injection history), the SGS solver routed
approximately 900,000 m³ downstream versus 160,000 m³ for the standard solver —
roughly 5× more. The standard solver's mean cell elevation across the embankment
face is above the notch sill, requiring upstream WSE to reach ~1.1 m before flow.
The SGS solver detects the sub-cell channel and routes water at WSE ~0.71 m (just
above the 0.705 m notch sill).

With the R-A (Manning hydraulic radius) SGS flux kernel, the ratio increases further
to approximately 14× at the same step count, reflecting the improved friction
representation for concentrated channel flow.

The synthetic DEM can be regenerated:

```bash
python test/synthetic_dem/generate_synthetic_dem.py
```

Parameters are stored in `test/synthetic_dem/synthetic_dem_params.json`.

---

## 5. Cartesian Backend Validation (`cartesian-formulation` branch)

**Purpose:** validate the new Cartesian mesh backend (`mesh/CartesianGrid.jl`)
and — using it as a provably-orthogonal reference grid — investigate whether
the solver's flux/momentum code carries a directional bias that is
independent of grid geometry. Full account, including every bug found along
the way: the project's internal development history (available on request).

### 5a. Mesh correctness (`test/audit_cartesian_geometry.jl`)

Checks adjacency (per-cell neighbour count, symmetry, and geometric
correctness of every adjacency pair) directly from `mesh.adjacency`, and
edge-geometry orthogonality (`cos_theta`, `skew_x`, `skew_y`) via the real
production `_build_edge_list`. No simulation.

```powershell
julia --project=. test\audit_cartesian_geometry.jl <mesh.parquet>
```

**Result (run, confirmed):** PASS at both 16 cells and 12,800 cells —
`cos_theta` exactly `1.0`, `skew` exactly `0.0` on every edge, zero
adjacency errors, valence distribution (2/3/4 at corners/edges/interior)
exactly as predicted for a rectangular grid.

### 5b. Flat-terrain point-spread / cardinal-vs-diagonal isotropy

Uses `test/square/square_domain.geojson` (a genuinely square 3×3 km AOI
centred on 0°N, 0°E, chosen specifically to remove any domain-shape
confound), no DEM (flat), centred point injection, closed boundaries.

```powershell
# Generate mesh (no DEM)
julia --project=. FloodModel.jl `
    --meshgen-cartesian test\square\square_domain.geojson --dx 25 `
    --meshout test\square\square_cartesian_dx25_flat.parquet --mesh-only

# Run (7200 s at a 60 s output interval gives frame 121 = t 7200 s)
julia --project=. --threads auto FloodModel.jl `
    --meshload test\square\square_cartesian_dx25_flat.parquet `
    --flow-model standard --injection-point 0.0,0.0,10.0 --closed-boundaries `
    --sim-duration 7200 --dt-max 2 `
    --output test\square\square_cartesian_dx25_flat_run.h5 --output-interval 60

# Analyse (sectors centred ON the 8 compass bearings, unlike test_point_spread.jl)
julia --project=. test\test_octagonal_anisotropy.jl `
    --input test\square\square_cartesian_dx25_flat_run.h5 `
    --frame 121 --source-lon 0.0 --source-lat 0.0
```

**Result (run, confirmed), default `--momentum-model edge`:**

| dx | n_wet (t=7200 s) | Cardinal/diagonal front-radius ratio | Deviation from 1.0 |
|---|---|---|---|
| 25 m | 2,557 | 0.9728 | 0.0272 |
| 12.5 m | 9,513 | 0.9662 | 0.0338 |

A ratio below 1.0 means *diagonal* directions reach further. **The deviation
grew by ~24% when the grid spacing was halved** — ordinary discretisation
truncation error would be expected to *shrink* by roughly 4×. This rules out
truncation error and matches the diagnostic signature A5's own bias showed
in the project's internal Phase A diagnostics. Leading hypothesis (untested
at time of writing): the per-edge momentum representation — see
`HYDRAULICS.md` §14 for the mechanism and the exact next test
(`--momentum-model cell` on the same mesh).

### 5c. Sloped-DEM planar symmetry (re-using §2b's fixture)

Same `test/planar_embankment/` fixture as §2b, with a Cartesian mesh
(`--meshgen-cartesian ... --dx 25 --dem ...`) in place of the A5 mesh.
On this exactly-orthogonal grid `--gradient-correction on` and `off` are
mathematically identical (the correction term is zero), so
`test_planar_symmetry.jl`'s baseline-vs-corrected comparison is not a real
before/after here — pass the same output file as both, and read the
reported `|asymmetry|` itself. Use `--sweep` rather than a single `--frame`
to see its trajectory over the run.

**Result (run, confirmed):** `|asymmetry|` decreased smoothly from 0.2356
(t=6,002 s) to 0.0575 (t=72,000 s); baseline and "corrected" were bit-for-bit
identical at every sampled frame (36,000 steps — a clean confirmation that
the correction term really is exactly inert on an orthogonal grid).

### 5d. Boundary-proximity check (`test/check_boundary_proximity.jl`)

Reports, per frame, how close the wetted region's edge is to each of the
four domain boundaries (N/S/E/W separately), to rule closed-wall effects
in or out as a contaminant of any asymmetry measurement.

```powershell
julia --project=. test\check_boundary_proximity.jl --input <run.h5> --sweep 20
```

**Result (run, confirmed) for §5c's run:** the northern and southern
boundary distances shrank from ~975–1000 m to ~500–525 m over 20 simulated
hours but never approached zero — the front never got close to either wall,
ruling out boundary contamination for that test. (The west distance stays
pinned near 150 m for the whole run — the source sits near that edge and
water flows east — which is expected, not a finding.)

---

## 6. Mesh Families (WP3, 2026-10)

The "arbitrary polygon" claim needs more than A5 and Cartesian meshes
(`FloodA5_Review_and_ForwardPlan_2026-10.md` §6.1). `mesh/mesh_families.py`
writes FloodA5 GeoParquet meshes that load with `--meshload`:

```bash
python mesh/mesh_families.py hex     --aoi test/square/square_domain.geojson --dx 25 --out hex25.parquet
python mesh/mesh_families.py rotcart --aoi test/square/square_domain.geojson --dx 25 --angle 30 --out rot25.parquet
python mesh/mesh_families.py voronoi --aoi test/square/square_domain.geojson --dx 25 --jitter 0.3 --lloyd 5 --out cvt25.parquet
python mesh/mesh_families.py voronoi --aoi test/square/square_domain.geojson --dx 25 --jitter 0.45 --centre centroid --out vor25.parquet
python mesh/mesh_families.py import  --input any_polygons.gpkg --out imported.parquet
```

| Family | Valence | Faces orthogonal to centre line? |
|---|---|---|
| `hex` | 6 | yes |
| `rotcart` | 4 | yes (a non-axis-aligned orthogonal grid) |
| `voronoi` (generator centres) | mixed, mostly 5–7 | yes (Voronoi property) |
| `voronoi --centre centroid` | mixed | no — a general non-orthogonal polygon mesh |
| `import` | any | depends on the input |

Generic meshes use the ID scheme `polygon` (row-index IDs). Adjacency is built
from shared polygon edges, so meshes with hanging nodes are not supported.
SGS is not yet available on these meshes (its tables are per A5 slot).

A Cartesian grid aligned with a raster in its own projected CRS (for
LISFLOOD-FP comparisons) comes from `--meshgen-cartesian dem.tif` (or `.asc`),
which also takes the bed elevation from the raster.

Tier 0 on every family is `test/test_mesh_families.jl` (in CI).

---

## Adding New Test Cases

New test cases should be placed in `test/<domain_name>/` and documented here.
Aim to include:

1. A GeoJSON AOI file and any associated DEM or boundary condition files.
2. A clear description of the domain, its purpose, and what behaviour is expected.
3. Command-line examples to reproduce the simulation from scratch.
4. Pass/fail criteria or comparison data where available.
5. A brief note on expected run time and hardware requirements.

