# FloodA5 — Test Cases

This document describes the test cases included in the repository. Each case
serves a specific purpose in the validation and benchmarking of FloodA5.
New cases will be added as the model is applied to additional domains.

---

## 1. Square Domain (Point Spread Test)

**Location:** `test/square/`  
**Purpose:** Confirm that flood spreading from a point source on a flat domain
is radially symmetric — i.e., that the pentagonal grid does not introduce
directional bias.

### Domain

A 1.35° × 1.35° square polygon centred at 0°N, 0°E (`test/square/square_domain.geojson`).
The domain straddles the null meridian and equator for a clean, unbiased test geometry.
No DEM is used — all cells have zero elevation. The domain is geometrically simple
enough that the expected behaviour can be checked visually without an analytical
solution.

### Simulation setup

Resolution 18 (approximately 500 m² per cell), standard flow solver:

```bash
# Generate mesh
julia --threads 1 --project=. FloodModel.jl \
    --meshgen test/square/square_domain.geojson --meshres 18 \
    --meshout test/square/square_mesh18_standard.parquet \
    --flow-model standard --mesh-only

# Run with a central point source (50 mm/hr applied to one cell)
julia --threads auto --project=. FloodModel.jl \
    --meshload test/square/square_mesh18_standard.parquet \
    --flow-model standard \
    --rainpoint 0.0,0.0,50.0 \
    --sim-duration 72000 --dt-max 10 \
    --output test/square/square_rainpoint_2hr.h5 \
    --output-interval 300
```

### Expected results

- Water propagates outward from the central cell in all five principal directions
  approximately equally.
- The flood front radius at any given time should be approximately equal in all
  directions (circular front, not diamond or cross-shaped).
- Mass balance error (`mb_err` in logs) should be < 0.01% throughout.
- No negative volumes or NaN values.

### Reference output

A pre-computed GeoPackage output for comparison is provided at
`test/square/square_mesh18_standard_output.gpkg`. This was produced with the
standard solver at resolution 18 and can be used to check that results are
consistent after code changes.

### Comparison with rectangular grids

A rectangular grid model run on the same domain at equivalent resolution will
produce a diamond-shaped flood front rather than a circular one — the directional
bias inherent in four-connected rectangular grids. This comparison is a key
motivating result for the pentagonal grid approach.

### Post-processing

The provided R script (`test/square/plot_output.R`) reads the HDF5 output and
plots the flood front at selected time steps.

---

## 2. Carlisle 2005 Flood Event

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

## 3. Synthetic DEM SGS Validation

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

## Adding New Test Cases

New test cases should be placed in `test/<domain_name>/` and documented here.
Aim to include:

1. A GeoJSON AOI file and any associated DEM or boundary condition files.
2. A clear description of the domain, its purpose, and what behaviour is expected.
3. Command-line examples to reproduce the simulation from scratch.
4. Pass/fail criteria or comparison data where available.
5. A brief note on expected run time and hardware requirements.
