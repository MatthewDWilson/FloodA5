# FloodA5 — Test Cases

This document describes the test cases included in the repository. Each case
serves a specific purpose in the validation and benchmarking of FloodA5.
New cases will be added as the model is applied to additional domains.

## The test-case registry (`cases/`)

The canonical cases are registered in `cases/registry/*.jl` and run with one
command. Each case generates its own inputs (DEM, mesh, hydrograph), runs
`FloodModel.jl` with explicit flags, computes its metrics with the analysis
module (`analysis/`), and checks them: physical checks (mass conservation,
the behaviour the case exists to show), then regression against the expected
values in `cases/expected/<case>.toml`, which record the version and commit
that produced them.

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl list              # cases, tiers, runs done
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl show CASE         # description and every command, to copy
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run CASE -P 2 -T 4
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl check CASE        # analyse existing runs, check
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl baseline CASE     # record new expected values (after a deliberate change)
```

Outputs go to `test/case_runs/<case>/` (git-ignored). A run is reused when its
output records the same command line and completed; `--fresh` runs it again
(e.g. after a code change). `run smoke` / `run standard` run a whole tier.

| Case | Tier | § | What it shows | Cost (16 cores) |
|---|---|---|---|---|
| `point_spread_smoke` | smoke (CI) | 1 | the runner end to end, 900 Cartesian cells | ~2 min |
| `point_spread` | standard | 1, 5b | flat point source: A4/A6, sector ratio, circularity on A5 res 16/18 and Cartesian 25/12.5 m | ~8 min |
| `planar_slope` | standard | 2b | N/S volume symmetry on a 0.1% plane, A5 res 16/18 and Cartesian 25 m | ~10 min |
| `planar_embankment` | standard | 2a | SGS passes water through a sub-cell notch near its sill; standard only at the crest | ~6 min |
| `channel` | standard | 7 | straight channel, 22/11 m cells, 0°/45°, A5/Cartesian, SWE/LI (16 runs) | ~15 min |
| `channel_series` | full | 7 | the full 224-run channel series | hours |
| `carlisle` | full | 3 | first 6 h with a 10 mm/h *test* rain (the §15.11 record), local-inertial vs SWE on A5 res 18 and Cartesian 20 m, plus SGS; no-rain runs for contrast | ~45 min |
| `carlisle_peak` | full | 3 | the 6 h before the Eden peak (49.5–55.5 h), river inflow only: local-inertial vs SWE on A5 and Cartesian, plus SGS | ~1.5 h |

The older run scripts that these cases replace are in `test/archive/` (see
its README): they rely on the pre-0.3.0 defaults and no longer reproduce the
numbers recorded with them.

---

## 1. Square Domain (Directional-Bias Benchmarks)

**Location:** `test/square/`
**Purpose:** Confirm that flood spreading from a point source on a flat domain
is radially symmetric — i.e., that the pentagonal grid does not introduce
directional bias — and provide the reference domain for the point-spread
directional-bias correction benchmark described in `docs/HYDRAULICS.md` §7–§8.

### Domain

A 0.01349° × 0.01349° square (about 1.5 × 1.5 km) centred at 0°N, 0°E
(`test/square/square_domain.geojson`; its "3x3 km" label is wrong). The domain
straddles the null meridian and equator for a clean, unbiased test geometry.
No DEM is used — all cells have zero elevation.

**Registry cases:** `point_spread` (A5 res 16/18, Cartesian 25/12.5 m; 10 m³/s
for 1 h, closed, current default scheme) and `point_spread_smoke` (CI). The
benchmark below documents the earlier correction study, under older schemes.

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

**Terrain generator:** `test/planar_embankment/generate_planar_embankment_dem.py`
(a west→east plane with an optional N–S embankment and a notch through it;
writes `planar_embankment_dem.tif`, `_params.json` and `_aoi.geojson` into
`--out`). Both cases below are registry cases (`cases/registry/`): each
generates its own DEM and meshes in `test/case_runs/<case>/`, so they cannot
disturb each other or the legacy fixture in `test/planar_embankment/`
(an `--emb-height 0` DEM still used by older run scripts).

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl show planar_embankment   # the commands, to copy
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_embankment -P 2 -T 2
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_slope -P 3 -T 2
```

### 2a. SGS sub-cell routing through a notch (`planar_embankment`)

A 2 × 0.5 km plane sloping 0.5% west→east (2 m DEM; 10 m at the west edge),
with a 40 m-wide, 2 m-high embankment across it at x = 1000 m (crest 7.0 m),
cut by a 4 m notch to the natural ground (sill 5.0 m) on the centre line — a
fifth of a res-18 A5 cell (~22 m). 10 m³/s enters at x = 200 m (an
`--inflow-point` hydrograph with a 10-minute ramp), closed boundaries, 6 h,
output every 2 min, on res-18 SGS and standard meshes (≈ 2 min per run
including start-up; the meshes ≈ 1.5 min each). The upstream storage to the
crest (a wedge against the embankment, 1.8·10⁵ m³) is reached within the run.

Intent: standard flow (cell-mean elevations) holds the water until the
upstream water surface reaches the crest it sees; SGS sees the notch and passes
water at the sill level. Metrics per output time: the upstream WSE (the highest
WSE of the *wet* cells along the embankment's upstream foot) and the volume
stored east of the embankment; from them the time downstream water first
exceeds 1 m³ and the upstream WSE then. Results (0.4.1):

| | first downstream water | upstream WSE then | WSE at 6 h | downstream at 6 h |
|---|---|---|---|---|
| SGS | 1 082 s | 5.23 m (sill 5.0) | 6.22 m | 39 000 m³ |
| standard | 16 202 s | 6.77 m (crest 7.0) | 6.99 m | 12 600 m³ |

Standard's crest is the lowest cell-mean route across the embankment, 6.72 m:
the notch lowers the means of the cells it crosses. Checks: mass balance and
stored volume; standard passes no water until the upstream WSE reaches that
crest, which lies within 0.5 m of the real one; SGS passes water while the
upstream WSE is ≥ 0.5 m below the crest and within 0.5 m of the sill; SGS has
more water downstream at every output time before standard overtops (126 of
them). Through the 4 m notch SGS conveys ≈ 3 m³/s at 1.1–1.2 m of head (5–6 h),
so its upstream level keeps rising under 10 m³/s.

(Replaces `test_planar_embankment.jl`, which shared the `--emb-height 0` DEM,
took the upstream WSE over dry cells and used a source too small to reach the
notch.)

### 2b. North/south volume-symmetry benchmark (`planar_slope`)

The generator with `--emb-height 0 --emb-width 0` gives a pure 0.1% plane
(4 × 2 km, 4 m relief) — a domain with no physical north/south distinction,
so any north/south asymmetry of a plume is a mesh/solver artefact. This is
the primary acceptance test of the directional-bias work (`docs/HYDRAULICS.md`
§7–§8, §15): the N/S volume split is an integrated quantity, far less noisy
than per-cell velocity directions.

The case injects 0.1 m³/s at x = 400 m on the centre line (the generator's
injection point), closed boundaries, 72 000 s with Δt ≤ 2 s, standard flow and
the default scheme, on A5 res 16, 18 and a 25 m Cartesian grid (20–90 s of
simulation per run). Metrics (`FA.planar_metrics` about the source cell) at 5,
10, 15 and 20 h: `asym` = (V_N − V_S)/(V_N + V_S) excluding a ±37.5 m band,
the lateral skew `A2`, the plume deflection, and flow-direction statistics at
20 h. Results (0.4.1, 20 h):

| | asym | deflection | mean flow direction |
|---|---|---|---|
| A5 res 16 | +0.089 | +1.9° | −3.4° |
| A5 res 18 | +0.088 | +0.8° | +2.4° |
| Cartesian 25 m | −1·10⁻⁵ | −5·10⁻⁵° | 0.0° |

Checks: mass balance and stored volume; Cartesian symmetric to 10⁻³ (its grid
is mirror-symmetric about the source row); A5 |asym| < 0.2 and deflection < 5°
at every time (before WP1 the res-18 value was 0.87 with 14°); flow within 5°
of downslope. The band is 37.5 m, not the historical 50 m: on the 25 m grid
±50 m falls on cell centres, and rounding decides their side (asym −0.11 for a
symmetric plume).

**Historical results** (other source positions, older schemes): resolution 16,
2 h, baseline `|asymmetry|` 0.25, 0.07 with gradient correction; res 18,
72 000 s, baseline 0.866, `--face-flux-method diamond --momentum-model cell`
~0.61 before WP0 (0.28 since); WP1's `diamond-wetdry` + vector friction +
hybrid limiter + local CFL 0.05 (HYDRAULICS §15.1); the 0.3.0 default 0.19
(§15.8). `test/test_planar_symmetry.jl` and `test/wp1_metrics.jl planar`
compute the same metrics for arbitrary output files.

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

### Registry cases `carlisle` and `carlisle_peak`: against the SWE reference

`cases/registry/carlisle.jl` (full tier) has two cases on the same meshes. Both
start with the river dry.

**`carlisle`: the first 6 h, with a test rain.** It runs the comparison of HYDRAULICS
§15.10–15.11: the first 6 h with 10 mm/h rain, on A5 res 18 and on a 20 m
Cartesian grid aligned with the LiDAR averaged to 20 m
(`cases/scripts/average_dem.jl`: GDAL `-tr 20 20 -r average -tap`, then
`--meshgen-cartesian dem20.tif`), each with the local-inertial default and
`--flow-model swe`, plus SGS on A5. It builds its own meshes in
`test/case_runs/carlisle/` (~7 min); the runs take ~50 min on CPU, one
4-thread run at a time (SWE 21 and 17 min).

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run carlisle -P 1 -T 4
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl show carlisle    # the commands, to copy
```

Results at 6 h (0.4.1; recorded values in brackets):

| run | storage (Mm³) | vs SWE: storage | rms Δh | CSI | bias | arrival (s) |
|---|---|---|---|---|---|---|
| A5 SWE | 2.690 (2.690) | | | | | 6 685 (6 688) |
| A5 local-inertial | 3.369 | 1.252 (1.25) | 0.29 m | 0.62 | 1.06 | 9 749 |
| Cartesian SWE | 2.331 (2.331) | | | | | 5 727 (5 736) |
| Cartesian local-inertial | 2.302 | 0.988 (0.99) | 0.13 m | 0.69 | 0.85 | 5 989 |
| A5 SGS | 4.300 (4.299) | | | | | 7 717 |

(§15.11's 3.363 and 2.297 Mm³ are the log's last step line, ~100 s before
6 h.) The checks: SWE storage and arrival within 1% and 2% of the record,
local-inertial/SWE ratios within ±0.05, the A5 ratio above the Cartesian one
by > 0.15, CSI ≥ 0.55 and rms Δh ≤ 0.35 m, SGS within 2%, mass balance.

**The rain is a test forcing.** The 10 mm/h uniform rain comes from the damping
series (HYDRAULICS §15.6–15.11); it is not the 2005 event's rain. It puts 60 mm
(0.88 Mm³) on every cell in 6 h, while the inflows are still low (Eden 100–230
m³/s). It is only 15% of the water added, but it makes most of the 6 h flooded
area. Two extra runs, `a5_li_norain` and `cart_li_norain`, repeat the
local-inertial runs without it. The area deeper than 0.1 m at 6 h is 3.06 km²
with rain and 1.25 km² without on A5, and 2.32 against 0.95 km² on Cartesian
(check: without rain < 60% of with rain). The 6 h maps of this case show mostly
rain-on-grid, not river flooding.

**`carlisle_peak`: the 6 h before the Eden's peak, river inflow only.** It
simulates the event from t = 0 to the Eden's peak (1 273 m³/s at 55.5 h; the
model has no start-time option) and compares the solvers at 49.5, 52.5 and
55.5 h. It has no rain, outputs every 1800 s, and reuses the `carlisle` meshes
(run `cases/run.jl prep carlisle` first). Its runs are the same five:
local-inertial and SWE on each mesh, and SGS on A5. The two SWE runs over 55.5 h
take several hours on CPU, at 6 threads each.

Results in the window (0.4.4; 49.5 / 52.5 / 55.5 h):

| run | storage at 55.5 h (Mm³) | ÷ SWE | rms Δh | CSI | bias | area > 0.1 m at 55.5 h (km²) |
|---|---|---|---|---|---|---|
| A5 SWE | 10.51 | | | | | 5.16 |
| A5 local-inertial | 8.01 | 0.74 / 0.76 / 0.76 | 0.53–0.55 m | 0.86–0.87 | 0.87–0.89 | 4.57 |
| Cartesian SWE | 9.36 | | | | | 4.91 |
| Cartesian local-inertial | 7.92 | 0.86 / 0.86 / 0.85 | 0.29–0.34 m | 0.92–0.93 | 0.93–0.94 | 4.57 |
| A5 SGS | 26.67 | | | | | 8.12 |

**Near the peak the local-inertial default stores less than the reference.** In
the first 6 h with the test rain it stores more. Near the peak:

- **It conveys more.** Its water is faster where the reference is ≥ 3 m deep:
  1.34× on A5 and 1.15× on Cartesian at 55.5 h. It lets more water out and stands
  0.3–0.6 m lower over the floodplain.
- **It is insensitive to the mesh.** The two local-inertial runs store 8.01 and
  7.92 Mm³; the two references differ by 12%.
- **SGS holds the Eden back by metres.** Its water surface is a median 3.0 m above
  the local-inertial one, and only 0.07 Mm³ is below the cell-mean bed. This
  points to a conveyance problem at high stage, still to be investigated.

Output frames drift late of their nominal times (by up to ~110 s at 55 h, because
each write is scheduled from the previous one), so frames are matched within
300 s. `notebooks/03_carlisle.ipynb` analyses both cases.

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

The `planar_slope` case (§2b) includes this run: with the 0.3.0 default its
asymmetry is ~10⁻⁵ (about the source cell). Historically: the same
`test/planar_embankment/` fixture as §2b, with a Cartesian mesh
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

## 7. Straight-Channel Representation Series (0.3.2, 2026-11)

A straight 20 m rectangular channel at 7 orientations and 2 cross-channel
offsets on A5 meshes (res 17–20, cells 44–5.6 m) and on Cartesian grids of
the same cell area, run with the full SWE solver and the local-inertial
default. It measures how much conveyance each mesh loses as the channel
narrows to one cell or less, and how that depends on orientation. Design,
results and conclusions are in `HYDRAULICS.md` §15.12.

### Domain

A 2 km square at 51° N, 0° E with an analytic bed:
- the channel, with 2 m banks;
- a floodplain rising across the valley at 0.5%;
- a 0.1% slope along the channel, steepening to 2% near the outflow as a
  hydraulic control.

Every edge is open. The inflow is 20 m³/s, constant. The in-bank normal depth is 0.97 m.

### Running the series

Two registry cases (`cases/registry/channel.jl`, outputs in
`test/case_runs/<case>/`) run it with the generator's job names and arguments:

```bash
# `channel` (standard tier): res 18-19, 0° and 45°, offset 0, both meshes and solvers (16 runs, ~15 min)
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run channel -P 2 -T 4
# `channel_series` (full tier): the whole series below (224 runs, hours; no limd runs)
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run channel_series -P 2 -T 4
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl show channel   # the commands, to copy
```

Metrics: FA's channel metrics per run and, under `series`, the rms and
spread of Δη over the orientations (the table in HYDRAULICS §15.12). Checks:
Δη falls with resolution, the runs are steady, local-inertial is within a
margin of SWE, SWE stays in bank (details in each case's `show`). By hand, as
recorded:

```bash
# 1. Meshes (112 GeoParquet files, about 0.5 GB), inflow file, manifest and job list
FLOODA5_NO_VIS=1 julia --project=. --threads 1 test/channel_bench/make_channel_meshes.jl test/wp1_runs/chan

# 2. Runs: 3 h each, 224 jobs ("chan_<mesh>_r<res>_a<angle>_o<offset>_<swe|li>")
bash test/run_jobs.sh test/wp1_runs/chan/jobs.txt 3 4 test/wp1_runs/chan
#    (res 20 is faster on a GPU: run those lines with --backend gpu appended)

# 3. Metrics (metrics.csv) and figures
FLOODA5_NO_VIS=1 julia --project=. test/channel_bench/analyse_channel.jl test/wp1_runs/chan test/wp1_runs/chan
julia --project=. test/channel_bench/plot_channel.jl test/wp1_runs/chan
```

`--res`, `--angles` and `--offsets` on the mesh script select a subset.
Wall time per run on 16 threads of a Xeon E-2278G:
- res 19 (33 000 cells): about 35 s local-inertial;
- res 20 (130 000 cells) on an RTX 4000 Ada: 2–3 min local-inertial, 7–14 min SWE on A5, 3–5 min SWE on Cartesian.

The whole series takes several hours on that machine.

### Expected results

Δη, the mean water-surface error in the channel against normal flow, behaves as follows:
- It falls from about 1.1 m at 44 m cells to ≤ 0.15 m at 5.6 m on both meshes.
- Cartesian is worst at 45°. A5 is worst at 0–15° and best at 60–90°.
- Neither mesh is better overall (HYDRAULICS §15.12).

The series itself is a benchmark: use it to compare schemes or meshes
against these figures. The registry cases' checks encode the robust parts
of them (above).

---

## Adding New Test Cases

Add a case to the registry: a file `cases/registry/<name>.jl` that calls
`register!(Case(...))` (see `cases/registry/point_spread.jl`, and the
`Case`, `Prep`, `Job`, `Check` docstrings in `cases/FloodA5Cases.jl`).

1. **Inputs**: tracked source files (AOI, DEM, boundary conditions) go in
   `test/<domain>/`; everything derived (resampled DEMs, meshes, hydrographs)
   is a `Prep` step writing into `{case}`, never shared with another case.
2. **Runs**: one `Job` per run, with every flag that matters explicit,
   including `--flow-model` (the default changed from sgs to standard in
   0.5.0) and any scheme choice, so
   that a change of defaults cannot silently change the case.
3. **Metrics**: computed with the analysis module (`FA`); add new metrics to
   `analysis/` as pure functions, not to the case file.
4. **Checks**: the behaviour the case exists to show, as physical checks
   with stated reasons; mass conservation (`mass_balance_rel`, from the
   run's summary).
5. **Expected values**: `cases/run.jl baseline <name>` once the checks pass
   on runs you trust; commit the TOML.
6. **Documentation**: a section here (purpose, domain, what to expect,
   cost), and the `title`/`description`/`doc` fields of the case.

Helpers in a registry file share one module: prefix them with the case's
abbreviation. A `smoke` case runs in CI (`test/test_cases.jl`), so it must
take seconds and need no Python.

