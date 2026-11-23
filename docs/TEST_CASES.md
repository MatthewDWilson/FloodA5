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
| `planar_wall` | standard | 2d | a wall with no gap: no water may pass below the crest; SGS passes a wall narrower than a cell at once (8 runs; the two SGS thin-wall checks are expected failures until SGS is fixed) | ~20 min |
| `sgs_strip` | standard | 2c | SGS against analytic normal flow, storage and arrival on a 4 km strip with a sub-cell channel (16 runs; several SGS checks fail by design: expected failures until SGS is fixed) | ~40 min |
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
(a west→east plane with an optional N–S embankment and a notch through it — `--notch-width 0`: none;
writes `planar_embankment_dem.tif`, `_params.json` and `_aoi.geojson` into
`--out`). The planar cases below are registry cases (`cases/registry/`): each
generates its own DEM and meshes in `test/case_runs/<case>/`, so they cannot
disturb each other or the legacy fixture in `test/planar_embankment/`
(an `--emb-height 0` DEM still used by older run scripts).

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl show planar_embankment   # the commands, to copy
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_embankment -P 2 -T 2
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_slope -P 3 -T 2
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_wall -P 2
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
exceeds 1 m³ and the upstream WSE then. Results (0.8.0; earlier SGS in
italics):

| | first downstream water | upstream WSE then | WSE at 6 h | downstream at 6 h |
|---|---|---|---|---|
| SGS | 962 s | 5.26 m (sill 5.0) | 6.41 m | 114 900 m³ |
| *SGS, 0.7.x (no critical-flow limit)* | *842 s* | *5.25 m* | *6.07 m* | *154 700 m³* |
| *SGS, 0.6.0* | *1 082 s* | *5.23 m* | *6.22 m* | *39 000 m³* |
| standard | 16 442 s | 6.78 m (crest 7.0) | 7.00 m | 9 600 m³ |

Standard's crest is the lowest cell-mean route across the embankment, 6.74 m:
the notch lowers the means of the cells it crosses. Checks: mass balance and
stored volume; standard passes no water until the upstream WSE reaches that
crest, which lies within 0.5 m of the real one; SGS passes water while the
upstream WSE is ≥ 0.5 m below the crest and within 0.5 m of the sill; SGS has
more water downstream at every output time before standard overtops (129 of
them). Until 0.6.0 SGS conveyed ≈ 3 m³/s through the notch at 1.1–1.2 m of
head and its upstream level kept rising; in 0.7.x it passed ≈ 9.8 m³/s at
≈ 1.07 m (5–6 h), more than critical flow through the notch (a 4 m
rectangle: ≈ 7.5 m³/s at that head): the local-inertial face flux had no
critical-flow limit. Since 0.8.0 SGS faces are limited to critical flow
(`--sgs-critical`, HYDRAULICS §9.2, §15.14): 8.24 m³/s at 1.38 m of head
between 5 and 6 h, and a further check holds it to the broad-crested-weir
discharge 1.705·W·H^1.5 of the 4 m notch (+2%; 11.07 m³/s at that head;
`q_down_5to6h`, `head_5to6h`).

**Where the water crosses** (metrics only, added with `planar_wall`, §2d; the
downstream volume alone does not separate the notch from a crossing
elsewhere). Cells by polygon: the first row of cells entirely east of the
embankment (24), cells wholly on its top, and cells straddling its east face;
"at the notch" = centre within 50 m of the notch's centre line, "away" =
polygon more than 100 m from it. SGS (0.6.0; the same holds in 0.7.0): at the first downstream water
(1 082 s) four row cells hold water, y = 213–316 m (the notch is at
248–252 m), none away from the notch; the wet stretch of the row then widens
continuously from the notch (y = 180–403 m at 2 h, 143–417 m at 6 h; the
first row cell more than 100 m away at 4 202 s, 3 of 13 at 6 h) with levels
4.8–5.0 m, below the ground at the upstream foot (5.1 m): water spreading
along the downstream toe. Of the 12 cells wholly on the embankment's top away
from the notch none ever holds water. **Nothing crosses the 40 m embankment
away from the notch**; the same in the standard run (first water in 5 row
cells at y = 190–293 m). Downstream volume by band (`v_down_notch_band_*`,
`v_down_off_notch_*`: cell centres within / beyond 50 m of the notch line) is
dominated by that spreading: SGS 7 100 / 31 800 m³ at 6 h.

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

### 2c. SGS on a strip with known answers (`sgs_strip`)

Written to answer whether SGS's large storage at Carlisle (26.7 Mm³ at the
Eden peak against 10.5 SWE, 8.0 LI) is how storage is measured or lost
conveyance. Notebook: `notebooks/04_sgs_strip.ipynb`.

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run sgs_strip -P 2     # ~40 min the first time
```

**Terrain** (`cases/scripts/sgs_strip_dem.jl`, 0.5 m DEM with the analytic
terrain at the pixel centres (since 0.8.0; at the corners before), EPSG:4326 in the
Cartesian generator's metric frame about 0°E 51°N; geometry and analytic
hydraulics in `cases/scripts/sgs_strip_geometry.jl`): a strip 4000 m long
and 60 m wide on a 0.1% slope, a 5 m × 1 m rectangular channel on its centre
line, Manning n = 0.03, no walls in the DEM (`--wall 0`); the mesh edges are
closed (`--bc-file`: a Closed line round the strip, a ZeroGradient line
across the downstream end). Below x = 3800 m the bed steepens to 2.1%, a
supercritical control that isolates the strip from each model's outflow
boundary (in a trial without it, the 2.5 m SWE run's transmissive outlet held
the water 0.86 m above normal depth at the outlet and was still filling at
6 h). Inflow 3 m³/s (in bank) or 20 m³/s (overbank) in the channel 31 m from
the upstream end, 10-minute ramp. Meshes: Cartesian 20 m (3 rows; the channel
a quarter of the middle row) with SGS tables or cell-mean beds; A5 res 18
likewise (every cell the AOI touches: a ragged strip 95 m wide on average); a
Cartesian 2.5 m reference. Runs (16): SGS (default, and `--sgs-conveyance
lumped`), LI and SWE on 20 m; SGS and LI on A5; LI and SWE on 2.5 m; 8 h
(2.5 m: 4 h).

**Analytic answers** (level δ above the floodplain at uniform flow): 2-D
(Σ strips h^(5/3), the fine-grid limit) −0.267 / +0.372 m with 3.57 /
27.34 m²/m of storage; 20 m cell-mean rows +0.045 / +0.409 m (7.71 /
29.52 m²/m); 20 m rows with their true section and one hydraulic radius
R = A/P per face (`--sgs-conveyance lumped`) −0.254 / +0.409 m. The SGS
default sums h^(5/3) over each face's profile and should give the 2-D level. Above the floodplain the volume above the cell-mean
bed *is* the water volume, so a cell-mean grid's storage differs from the
true volume only where cells are partly wet (in bank: 2.2× here).
`cases/scripts/sgs_strip_tables.jl` (a preparation step) compares the SGS
tables with the terrain directly.

Results (0.7.0, unchanged to the digits shown in 0.8.0; measuring reach
x = 1000–2000 m, last frame; storage in m²/m; the 0.6.0 SGS values are in
`cases/expected/v0.6.0/`):

| | 3 m³/s δ (m) | storage | 20 m³/s δ (m) | storage | dS/dt ÷ Q at 20 m³/s |
|---|---|---|---|---|---|
| analytic 2-D | −0.267 | 3.57 | 0.372 | 27.34 | 0 |
| analytic, lumped 20 m rows | −0.254 | 3.64 | 0.409 | 29.52 | 0 |
| LI / SWE 2.5 m | −0.263 / −0.264 | 3.57 / 3.56 | 0.373 / 0.373 | 27.38 | 0 |
| LI / SWE 20 m | 0.045 | 7.70 | 0.408 | 29.52 | 0 |
| **SGS 20 m** | **−0.279** | **3.50** | **0.362** | **26.71** | 0 |
| SGS 20 m, `--sgs-conveyance lumped` | −0.271 | 3.55 | 0.399 | 28.91 | 0 |
| LI A5 res 18 | 0.095 | 14.09 | 0.390 | 42.05 | 0 |
| **SGS A5 res 18** | **−0.277** | **3.50** | **0.340** | **37.35** | 0 |
| *SGS 20 m, 0.6.0* | *−0.001* | *5.05* | *4.75* | *106.7* | *0.75* |
| *SGS A5 res 18, 0.6.0* | *0.173 (filling)* | *17.43* | *0.783* | *49.33* | *0.91* |

(The 2-D analytic level for the A5 strip's own mean width, 95 m, is +0.290 m
overbank; its ragged edges hold slack water in both A5 models.)

**What the case found in 0.5.1** (and what 0.7.0 changed; HYDRAULICS §9):
(1) storage is measured correctly — the tables hold the terrain's volume to
0.3% (now: 0.3% Cartesian, 3% at worst on A5), and overbank SGS's water
volume equals its volume above the cell-mean bed; (2) overbank, an SGS face's
flow area stopped growing when the water rose above the cell's highest DEM
sample (edge curves on the cell's own levels; the top value above them; the
kernel averaged the two cells' A and R): at δ = 0.4 m faces carried 39%
(channel) and 2.5% (floodplain) of their true flow area, and the 20 m tables
at most 4.9 m³/s at the bed slope — now one table per face, on the face's own
levels and continued above them (face areas within 1–2%, conveyance within
3%; the tables carry 417 m³/s); (3) in bank, the legacy donor cap (V/(2N) per
edge and step) made a channel cell hold 2N·Q·Δt before it passed Q — now
positivity limiting. Getting the A5 strip right also needed a time-step limit
from the faces for cells holding only a sliver of the channel (they
oscillated and held the channel up by up to 0.5 m) and the flow-path geometry
of partly wet cells (the channel zig-zags through A5 cells: 0.054 m high
without it, 0.010 m low with it).

Checks: the 2.5 m references against the analytic levels (0.03 m, 3%) and
in-bank arrival times (10%); the 20 m LI/SWE runs against the cell-mean
analytic; A5 LI steady and within 0.15 m overbank; the SGS tables against the
terrain (volume 1%; face areas and conveyance within 10% / 5%; 20 m³/s at the
bed slope); SGS 20 m at the 2-D level within 0.03 m and 5% in and out of
bank; the lumped runs at the lumped level; SGS arrival within 15% of the
2.5 m reference; A5 SGS within 0.05 m in bank, and overbank no further from
the 2-D answer for its width than A5 LI; the measurement identity; mass. All
pass since 0.7.0 (eight were expected failures before).

### 2d. A wall with no gap: is a barrier inside a cell a barrier? (`planar_wall`)

SGS stores water with one hypsometric curve per cell (one level pool) and
conveys it through each face with the DEM profile sampled **on the face
line**; nothing between the faces is used. A wall that crosses the inside of
a cell therefore lowers the cell's stored volume but cannot hold water back:
water enters over the upstream faces, is stored at the cell's lowest ground
(on a slope, the downstream side of the wall) and leaves through the
downstream faces at natural ground (SGS conveyance fix plan §2.3; Phase 0,
item 1). `planar_embankment` cannot show it, because water may use the notch.

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_wall -P 2     # ~20 min the first time
```

**Terrain**: the §2a plane (2 × 0.5 km, 0.5% west→east, 2 m DEM) with a 3 m
wall across it and **no notch** (`--notch-width 0`), one DEM and four meshes
per wall in `test/case_runs/planar_wall/{wide,thin}/`:

| wall | extent | crest | against the mesh |
|---|---|---|---|
| `wide` | 40 m, x = 980–1020 m | 8.00 m | wider than a res-18 A5 cell (~22 m) and two Cartesian 20 m cells |
| `thin` | 4 m (2 pixels), x = 988–992 m | 8.05 m | narrower than any cell; inside the Cartesian column x = 981.1–1001.1 m |

**Design**: 10 m³/s (10-minute ramp) at x = 200 m on the centre line, closed
boundaries, 6 h: 213 000 m³. As a level pool against the wall's upstream face
that is a wedge, V(η) = 5·10⁴ (η − z_w)² m³, so the level ends 2.06 m above
the ground at the face: 7.16 m (wide) and 7.12 m (thin), 0.84 and 0.93 m
below the crests. In a correct model no water passes. Runs: SGS and standard
flow, A5 res 18 and Cartesian 20 m, each wall (8; ~1 min each on 2 threads).

**Metrics** per run (cells sorted by polygon: entirely west of the wall,
entirely east, straddling): the volume east of the wall through time
(`v_down_*h`, `v_down_end`, `frac_down_end`), the time it first exceeds 1 m³,
the highest water level of the wet cells along the wall's upstream foot then
and at the end (`wse_up_*`), the highest level at the wall then including the
cells the wall crosses (`wse_wall_at_first_down`), and where the water
crossed: the cells of the first row east of the wall that hold water at that
time (`cross_n_cells`, `cross_y_min/max`), against the foot cells holding
water then (`touch_*`: the stretch of wall the water has reached), and the
fraction of the row holding water through time (`row_wet_frac_*`).

Results of the scheme up to 0.6.0 (`cases/run.jl run planar_wall -P 2`;
0.7.0 below):

| run | downstream at 6 h | first downstream water | level at the wall then (crest) | upstream level at 6 h |
|---|---|---|---|---|
| SGS, A5, wide | 0 | never | — (8.00) | 6.57 m |
| standard, A5, wide | 0 | never | — (8.00) | 7.06 m |
| SGS, Cartesian, wide | 0 | never | — (8.00) | 6.64 m |
| standard, Cartesian, wide | 0 | never | — (8.00) | 7.12 m |
| **SGS, A5, thin** | **124 660 m³ (58.5%)** | **722 s** | 5.18 m (8.05) | 5.27 m |
| standard, A5, thin | 180 756 m³ (84.9%) | 2 642 s | 5.44 m (8.05) | 5.61 m |
| **SGS, Cartesian, thin** | **136 972 m³ (64.3%)** | **962 s** | 5.10 m (8.05) | 5.15 m |
| standard, Cartesian, thin | 175 467 m³ (82.4%) | 3 842 s | 5.68 m (8.05) | 5.72 m |

- **SGS passes the thin wall as soon as the water reaches it**, 2.9 m below
  the crest and ≈ 0.1 m above the ground at its foot (5.06 m), on both meshes.
  It crosses wherever the water touches the wall, not at particular cells: at
  the first output with water downstream, 10 row cells hold water on A5
  (y = 324–531 m; 13 foot cells hold water, y = 213–507 m — the plume has
  drifted north on A5) and 7 on the Cartesian grid (y = 180–320 m; 9 foot
  cells, y = 160–340 m, symmetric about the source row); within an hour 93%
  (A5) and 100% (Cartesian) of the row holds water. The upstream level never
  rises above 5.3 m: the wall is not there.
- **The wide wall holds on both meshes and both models**: every route across
  it passes a face that lies on the wall. (SGS reaches only 6.6 m at the wall
  against 7.1 m in standard flow, because its conveyance defects, §2c, keep
  water near the source: 5–8 m deep at x < 200 m. The wall is therefore
  tested by SGS only up to that level.)
- **Standard flow overtops the thin wall early**, at the level of the ridge
  that the wall makes in the cell-mean beds (the lowest cell-mean route:
  5.36 m on A5, 5.63 m Cartesian, against a crest of 8.05 m). That is the
  known limit of cell-mean beds, not the SGS defect; the check for these two
  runs says so (no water downstream until the upstream level reaches that
  ridge).

**0.7.0: SGS holds both walls.** With connectivity (HYDRAULICS §9.3) a
point of a face conveys only from the level at which it is connected, below
the water, to the main basins of both cells, so the faces behind the wall are
closed up to its crest. No water passes either wall in any SGS run; the pond
rises to 7.07 m (A5) and 7.14 m (Cartesian) against the 4 m wall and 7.08 /
7.15 m against the 40 m wall (level pool 7.12 / 7.16 m, crests 8.05 / 8.00 m).
The wall may be moved by up to one cell (the faces that close are those of
the cells it crosses). 0.8.0 (DEM at pixel centres): the same, 7.06 / 7.14 m
and 7.08 / 7.15 m; the thin wall's cell-mean ridge is 5.39 m on A5 (was 5.36)
and 5.64 m Cartesian (5.63), and standard flow on A5 crosses it at 2 882 s
with the water at 5.48 m.

Checks: mass balance and stored volume; the design (level-pool level ≥ 0.5 m
below the crest; the upstream level stays below the crest and the level-pool
level in every run); no water downstream (< 1 m³ at 6 h) for the four wide
runs and the two SGS thin runs (expected failures until 0.7.0); the standard
thin runs as above.

**A wall on a face is a different test.** In a trial with the thin wall at
x = 998–1002 m the Cartesian face at x = 1001.1 m lay on it, and the
Cartesian SGS run held the water until the upstream level reached 6.31 m
(16 442 s; 15 000 m³ downstream at 6 h), then passed it along the whole wall
at once. 6.31 m is the lowest elevation of that face's profile as sampled
up to 0.7.1: the DEM was sampled with pixel values placed at pixel corners
rather than centres (`_affine_to_pixel` returns corner-based coordinates and
`_bilinear` treated them as centre-based), so the profile was taken half a
pixel (1 m) east of the face, on the ramp that interpolation makes of the
wall's east side (8.0 → 4.99 m over one pixel). Sampled at the face itself
it is 7.82 m. The half-pixel registration applied to every DEM sample (SGS
tables, face profiles and cell-mean elevations); fixed in 0.8.0, with every
case rerun.

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

Results at 6 h (0.8.0; recorded values in brackets):

| run | storage (Mm³) | vs SWE: storage | rms Δh | CSI | bias | arrival (s) |
|---|---|---|---|---|---|---|
| A5 SWE | 2.641 (2.641) | | | | | 6 676 (6 676) |
| A5 local-inertial | 3.297 | 1.248 (1.25) | 0.28 m | 0.61 | 1.06 | 9 647 |
| Cartesian SWE | 2.331 (2.331) | | | | | 5 727 (5 736) |
| Cartesian local-inertial | 2.302 | 0.988 (0.99) | 0.13 m | 0.69 | 0.85 | 5 989 |
| A5 SGS | 1.775 (1.775) | | | | | 3 570 |

The A5 values moved in 0.8.0 with the DEM sampled at pixel centres (0.4.1–0.7.1:
A5 SWE 2.690 Mm³ and 6 688 s, local-inertial 3.369, SGS 1.755); the Cartesian
mesh takes the raster's pixels and its runs are bit-identical.

(§15.11's 3.363 and 2.297 Mm³ are the log's last step line, ~100 s before
6 h.) The checks: SWE storage and arrival within 1% and 2% of the record,
local-inertial/SWE ratios within ±0.05, the A5 ratio above the Cartesian one
by > 0.15, CSI ≥ 0.55 and rms Δh ≤ 0.35 m, SGS within 2% of its record, mass
balance. (SGS until 0.6.0: 4.300 Mm³, arrival 7 717 s — conveyance held back.)

**The rain is a test forcing.** The 10 mm/h uniform rain comes from the damping
series (HYDRAULICS §15.6–15.11); it is not the 2005 event's rain. It puts 60 mm
(0.88 Mm³) on every cell in 6 h, while the inflows are still low (Eden 100–230
m³/s). It is only 15% of the water added, but it makes most of the 6 h flooded
area. Two extra runs, `a5_li_norain` and `cart_li_norain`, repeat the
local-inertial runs without it. The area deeper than 0.1 m at 6 h is 3.02 km²
with rain and 1.19 km² without on A5, and 2.32 against 0.95 km² on Cartesian
(check: without rain < 60% of with rain). The 6 h maps of this case show mostly
rain-on-grid, not river flooding.

**`carlisle_peak`: the 6 h before the Eden's peak, river inflow only.** It
simulates the event from t = 0 to the Eden's peak (1 273 m³/s at 55.5 h; the
model has no start-time option) and compares the solvers at 49.5, 52.5 and
55.5 h. It has no rain, outputs every 1800 s, and reuses the `carlisle` meshes
(run `cases/run.jl prep carlisle` first). Its runs are the same five —
local-inertial and SWE on each mesh, and SGS on A5 — and SGS without the
critical-flow limit (`a5_sgs_nocrit`, `--sgs-critical off`: the 0.7.x
hydraulics on the same tables). The two SWE runs over 55.5 h
take several hours on CPU, at 6 threads each.

Results in the window (0.8.0; 49.5 / 52.5 / 55.5 h):

| run | storage at 55.5 h (Mm³) | ÷ SWE | rms Δh | CSI | bias | area > 0.1 m at 55.5 h (km²) |
|---|---|---|---|---|---|---|
| A5 SWE | 10.68 | | | | | 5.21 |
| A5 local-inertial | 8.00 | 0.73 / 0.74 / 0.75 | 0.57–0.59 m | 0.86–0.88 | 0.87–0.89 | 4.63 |
| Cartesian SWE | 9.36 | | | | | 4.91 |
| Cartesian local-inertial | 7.92 | 0.86 / 0.86 / 0.85 | 0.29–0.34 m | 0.92–0.93 | 0.93–0.94 | 4.57 |
| A5 SGS | 6.25 | 0.55 / 0.57 / 0.58 | 0.99–1.02 m* | 0.68–0.71* | 0.70–0.73* | 3.55* |
| A5 SGS, no critical-flow limit | 5.98 | 0.53 / 0.55 / 0.56 | 1.00–1.03 m* | 0.64–0.67* | 0.66–0.69* | 3.36* |
| *A5 SGS (0.7.x)* | *5.94* | *0.52 / 0.55 / 0.57* | *0.97–0.99 m* | *0.61–0.67* | *0.63–0.69* | |
| *A5 SGS (0.6.0)* | *26.67* | | | | | |

(0.4.4–0.7.1, before the DEM registration fix: A5 SWE 10.51 Mm³, A5
local-inertial 8.01, 0.74–0.76; the Cartesian runs are bit-identical.)

\* SGS against SWE: mean depth V/A for CSI and bias (wet ≥ 0.1 m), water level
where both are wet for rms Δh (`_carl_sgs_compare`): an SGS depth is measured
from the cell's lowest ground, so its depth field is not comparable with a
cell-mean model's.

**Near the peak the local-inertial default stores less than the reference.** In
the first 6 h with the test rain it stores more. Near the peak:

- **It conveys more.** Its water is faster where the reference is ≥ 3 m deep:
  1.34× on A5 and 1.15× on Cartesian at 55.5 h. It lets more water out and stands
  0.3–0.6 m lower over the floodplain.
- **It is insensitive to the mesh.** The two local-inertial runs store 8.00 and
  7.92 Mm³; the two references differ by 14%.
- **SGS conveys more than local-inertial.** Until 0.6.0 it held the Eden
  back by metres (26.7 Mm³; the conveyance defects of §2c). It stores 0.58 of
  the 22 m SWE run (0.56 without the critical-flow limit of 0.8.0) where
  local-inertial stores 0.75. The 22 m SWE run is not a sound reference
  (below): since 0.8.0 the SGS bar is checked against the 5 m SWE benchmark in
  `carlisle_peak_fine` (until 0.7.1 it was an expected failure here).

Output frames drift late of their nominal times (by up to ~110 s at 55 h, because
each write is scheduled from the previous one), so frames are matched within
300 s. `notebooks/03_carlisle.ipynb` analyses both cases.

**`carlisle_peak_fine`: the 5 m benchmark (GPU).** The `carlisle_peak` event on
the benchmark suite's 5 m Cartesian mesh (582 012 cells; `julia --project=.
benchmarks/run.jl` builds it), local-inertial (GPU Float64, 1 h 10 min) and
SWE (GPU Float32, 3 h 51 min). Not part of `carlisle_peak`, which must run
on a CPU. **The 5 m SWE run is the benchmark** for the 22 m models (decided
2026-11-23; it replaces the 22 m SWE runs). The 5 m mesh takes the raster's
pixels, so the DEM registration fix of 0.8.0 does not change it. The case's
analysis compares every 22 m run of `carlisle_peak` (run that first) with
both 5 m runs at 49.5, 52.5 and 55.5 h (`_carl_fine_metrics`, entries
`vs5:<job>`).

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run carlisle_peak_fine -P 2
```

Compared on the 22 m cells (5 m results aggregated onto each coarse cell;
water level where at least half the cell is wet; CSI of the extent at 0.1 m;
SGS wetted area from its tables), at 55.5 h (0.8.0; HYDRAULICS §15.14):

| 22 m model | vs 5 m SWE: storage | level median (rms) | CSI | vs 5 m LI: storage | level median (rms) | CSI |
|---|---|---|---|---|---|---|
| SWE, A5 | 1.26 | +0.39 (0.56) m | 0.81 | 1.49 | +0.82 (0.84) m | 0.78 |
| SWE, Cartesian | 1.12 | +0.15 (0.29) m | 0.85 | 1.33 | +0.55 (0.58) m | 0.82 |
| local-inertial, A5 | 0.94 | −0.18 (0.25) m | 0.85 | 1.12 | +0.24 (0.34) m | 0.84 |
| local-inertial, Cartesian | 0.95 | −0.22 (0.29) m | 0.85 | 1.12 | +0.20 (0.32) m | 0.85 |
| SGS, A5, no critical-flow limit | 0.71 | −0.49 (0.53) m | 0.77 | 0.84 | −0.14 (0.23) m | 0.81 |
| **SGS, A5** | **0.74** | **−0.47 (0.52) m** | **0.81** | **0.87** | **−0.12 (0.16) m** | **0.85** |

The 22 m SWE runs store 12–26% more than 5 m SWE, so they are not
references. SGS reproduces the fine local-inertial solution best (level and
extent at all three times); but the local-inertial equations themselves store
0.84 of SWE's volume at 5 m, and 22 m local-inertial lands close to 5 m SWE
through compensating errors (its smoothed beds store ~1.15× the 5 m
local-inertial run). Checks: mass; every 22 m run compared; and the agreed
SGS bar — no further from 5 m SWE than local-inertial on the same mesh in
storage, rms level and CSI at all three times — an **expected failure**: the
gap is the hydraulics (local-inertial), not the SGS geometry. Closing it
needs a hydraulic formulation with momentum on the SGS geometry.

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
SGS works on them (tables per neighbour slot, `test/test_sgs_general.jl`);
from the command line, `--meshload mesh.parquet --dem dem.tif --flow-model sgs
--meshout mesh_sgs.parquet --mesh-only` builds the tables (0.8.1). The
discrete global grids (H3, IGEO7, ISEA3H, HEALPix, rHEALPix) are in
`docs/MESHES.md`; §8 runs the benchmarks on all the families.

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

## 8. The Benchmarks on Every Mesh Family (0.8.1)

Four registry cases (`cases/registry/mesh_families.jl`) run the point-spread,
planar-slope, straight-channel and Carlisle benchmarks on eleven meshes at
matched mean cell area, √area ≈ 22.2 m (A5 res 18):

| Mesh | Generator | Cells (√ mean area) |
|---|---|---|
| `cart` | `--meshgen-cartesian --dx 22.2` | 22.2 m |
| `rotcart` | `mesh_families.py rotcart --angle 30` | 22.2 m |
| `hex` | `mesh_families.py hex` (side 13.8 m) | 22.2 m |
| `voronoi` | `mesh_families.py voronoi --jitter 0.3 --lloyd 5` (centroidal) | 22.1 m |
| `voronoi_jit` | `mesh_families.py voronoi --jitter 0.3 --lloyd 0` | 22.1 m |
| `a5` | `--meshgen-dggs --dggs a5 --dx 22.2` (res 18) | 22.2 m |
| `h3` | `--meshgen-dggs --dggs h3 --dx 22.2` (level 12) | 16.6 m (15.1 at the equator) |
| `igeo7` | `… --dggs igeo7` (level 13) | 22.9 m |
| `isea3h` | `… --dggs isea3h` (level 23) | 23.3 m |
| `healpix` | `… --dggs healpix` (level 18) | 24.9 m |
| `rhealpix` | `… --dggs rhealpix` (level 12) | 17.3 m |

The global grids' levels are discrete, so H3 and rHEALPix are about 25% finer
and HEALPix 12% coarser than the target; every job records its `cell_size_m`.
Each case generates its own meshes in `test/case_runs/<case>/`.

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run mesh_families_point_spread -P 3 -T 4   # ~5 min
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run mesh_families_planar_slope -P 3 -T 4   # ~10 min
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run mesh_families_channel -P 3 -T 4        # ~45 min
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run mesh_families_carlisle -P 3 -T 4       # ~10 h (full tier)
```

The numbers below were recorded with the branch of 0.8.1 (the expected values
in `cases/expected/mesh_families_*.toml`). A5 and Cartesian reproduce the
`point_spread` (A5 res 18) and `channel` expected values bit for bit.

### `mesh_families_point_spread`: isotropy on a flat bed

10 m³/s for 1 h at the centre of the closed, flat 1.5 km square at (0°, 0°).
A4 and A6 are the four- and six-fold harmonics of the depth about the source
(0 for a circle); circularity is the hull's.

| Mesh | A4 | A6 | Front radius CV | Circularity |
|---|---|---|---|---|
| cart | 0.034 | 0 | 0.027 | 0.987 |
| rotcart | −0.017 | 0 | 0.012 | 0.987 |
| hex | 0 | −0.0070 | 0 | 0.996 |
| voronoi | 0.017 | −0.0003 | 0.013 | 0.993 |
| voronoi_jit | 0.012 | −0.0006 | 0.008 | 0.992 |
| a5 | −0.0015 | 0.00003 | 0.028 | 0.990 |
| h3 | 0.0001 | 0.0052 | 0.005 | 0.993 |
| igeo7 | 0.0000 | 0.0030 | 0.013 | 0.991 |
| isea3h | 0.0002 | 0.0072 | 0.010 | 0.992 |
| healpix | −0.037 | 0.0006 | 0.061 | 0.983 |
| rhealpix | 0.028 | 0.0002 | 0.016 | 0.990 |

Every mesh spreads nearly circularly (> 0.98). The bias follows the face
directions: the four-direction meshes have a four-fold term (Cartesian,
rHEALPix, and HEALPix, whose equatorial cells are near-squares turned 45°:
A4 < 0), the hexagonal ones a weaker six-fold term, and A5, with its spread
face directions, neither. HEALPix has the least round front (CV 0.061).

### `mesh_families_planar_slope`: directional bias on a slope

`planar_slope` (0.1 m³/s on a 0.1% west→east plane at 51° N, 20 h) on every
mesh. asym is the north/south volume asymmetry, deflection the plume's
direction against east.

| Mesh | asym 5 h | asym 20 h | deflection 20 h (°) | mean flow direction (°) | cells within ±30° |
|---|---|---|---|---|---|
| cart | −0.00001 | −0.00001 | −0.0001 | 0.0 | 0.97 |
| rotcart | −0.081 | −0.055 | −0.65 | −1.9 | 0.94 |
| hex | −0.00001 | −0.00001 | −0.0000 | 0.0 | 0.89 |
| voronoi | −0.024 | −0.013 | 0.04 | −0.5 | 0.94 |
| voronoi_jit | 0.010 | 0.005 | 0.07 | 1.1 | 0.98 |
| a5 | −0.083 | 0.088 | 0.81 | 2.4 | 0.96 |
| h3 | −0.061 | −0.044 | −0.42 | −1.0 | 0.93 |
| igeo7 | 0.050 | 0.046 | 0.46 | −0.1 | 0.95 |
| isea3h | −0.128 | −0.087 | −0.64 | −1.0 | 0.97 |
| healpix | −0.104 | −0.270 | −3.34 | −0.7 | 0.77 |
| rhealpix | −0.133 | −0.128 | −0.78 | −0.5 | 0.99 |

The two meshes that are mirror-symmetric about the source row (Cartesian,
hexagons) are symmetric to 1e-5, as they must be (a check). The Voronoi
meshes, irregular but orthogonal, are the least biased of the rest
(|asym| ≤ 0.024). The global grids sit between 0.04 and 0.13, as A5 does —
except HEALPix, whose skewed polar-region rhombi (33° non-orthogonal at
51° N, `docs/MESHES.md` §3) deflect the plume by 3.3° and give 0.27. Mean
flow directions are within 2.4° of east on every mesh (a check). The
fraction of moving cells within ±30° of east, which `planar_slope` checks at
≥ 90%, is recorded but not checked here: at 22 m the plume moves only 22–160
cells, and the radial flow next to the source decides it — the hexagons'
four cells at ±40° (two symmetric pairs) give 0.89; HEALPix's five, lopsided,
go with its bias.

### `mesh_families_channel`: orientation dependence in a one-cell channel

`channel` at 22 m (res 18, the 20 m channel about one cell wide) at 0°, 15°,
30° and 45°, no offset, SWE and local-inertial; meshes from
`make_channel_meshes.jl --families`. Δη is the in-channel stage error against
normal flow (m); spread = max − min of Δη over the four orientations, the
mesh's orientation dependence.

| Mesh | LI Δη at 0 / 15 / 30 / 45° | LI rms | LI spread | SWE rms | SWE spread |
|---|---|---|---|---|---|
| cart | 0.74 / 0.82 / 0.99 / 1.29 | 0.98 | 0.55 | 0.91 | 0.34 |
| rotcart | 0.93 / 0.70 / −0.09 / 0.95 | 0.76 | 1.04 | 0.80 | 0.91 |
| hex | 0.42 / 0.70 / 0.92 / 0.87 | 0.75 | 0.50 | 0.82 | 0.45 |
| voronoi | 0.72 / 0.70 / 0.81 / 1.14 | 0.86 | 0.44 | 0.89 | 0.28 |
| voronoi_jit | 0.69 / 0.55 / 0.93 / 1.07 | 0.84 | 0.53 | 0.92 | 0.30 |
| a5 | 1.06 / 1.07 / 0.97 / 1.00 | 1.03 | 0.10 | 1.06 | 0.09 |
| h3 (16.6 m) | 0.46 / 0.57 / 0.61 / 0.57 | 0.56 | 0.15 | 0.65 | 0.18 |
| igeo7 | 0.58 / 0.36 / 0.83 / 1.05 | 0.75 | 0.70 | 0.81 | 0.63 |
| isea3h | 0.68 / 0.74 / 0.63 / 0.97 | 0.77 | 0.34 | 0.83 | 0.36 |
| healpix (24.9 m) | 1.16 / 1.16 / 1.09 / 1.12 | 1.14 | 0.07 | 1.07 | 0.06 |
| rhealpix (17.3 m) | 0.28 / 0.71 / 0.93 / 1.07 | 0.80 | 0.79 | 0.78 | 0.76 |

At one cell per channel width every mesh loses conveyance (Δη > 0, the water
too high), and each loses least where its faces line up with the channel:
the rotated grid at 30°, Cartesian, rHEALPix and the hexagons at 0°. Hence
the regular meshes' large spreads (0.5–1.0 m). A5 and HEALPix are nearly
independent of orientation (0.07–0.10 m), but at about the Cartesian mean; H3
does best overall, helped by its 25% smaller cells. SWE and local-inertial
agree on all of this (both columns), so it is the meshes, not the
local-inertial scheme. Water leaves the channel in few runs: local-inertial
on Cartesian at 45° (32 m of wet floodplain per metre of reach) and on HEALPix
at 15–45° (10–12 m; SWE 3.5 m at 30–45°).

### `mesh_families_carlisle`: a real event on every mesh

`carlisle_peak` (the Eden peak, 55.5 h, inflows only) on the families that
case does not run, local-inertial and SGS (tables from the 5 m DEM); the
planar meshes keep the cells meeting the domain polygon (`--clip`). Every
run, with `carlisle_peak`'s A5 and 20 m Cartesian runs, is compared on its
own cells with `carlisle_peak_fine`'s 5 m runs (HYDRAULICS §15.14). Needs
those two cases run first; 18 runs, about 10 h on 3 × 4 threads. At 55.5 h,
against the 5 m SWE benchmark (storage ratio, rms water-level difference,
CSI at 0.1 m); "rms vs 5 m LI" is SGS against the 5 m local-inertial run,
which SGS approximates; "extrema" is the fraction of fully wet cells that are
local water-level extrema (`FA.noise_metrics`; a checkerboard).

| Mesh (√ mean area) | LI storage | LI rms (m) | LI CSI | SGS storage | SGS rms (m) | SGS CSI | SGS rms vs 5 m LI | SGS extrema |
|---|---|---|---|---|---|---|---|---|
| cart (20.0 m) | 0.946 | 0.289 | 0.853 | – | – | – | – | – |
| rotcart (22.2 m) | 0.953 | 0.291 | 0.845 | 0.959 | 0.417 | 0.882 | 0.526 | **0.42** |
| hex (22.2 m) | 0.920 | 0.265 | 0.841 | 0.774 | 0.492 | 0.872 | 0.159 | 0.08 |
| voronoi (22.2 m) | 0.929 | 0.275 | 0.840 | 0.749 | 0.544 | 0.842 | 0.209 | 0.14 |
| voronoi_jit (22.2 m) | 0.970 | 0.268 | 0.831 | 0.756 | 0.530 | 0.843 | 0.204 | 0.008 |
| a5 (22.2 m) | 0.943 | 0.254 | 0.848 | 0.735 | 0.515 | 0.813 | 0.159 | 0.009 |
| h3 (16.4 m) | 0.868 | 0.283 | 0.877 | 0.714 | 0.564 | 0.815 | 0.325 | 0.18 |
| igeo7 (22.9 m) | 0.943 | 0.269 | 0.832 | 0.823 | 0.598 | 0.829 | 0.564 | 0.17 |
| isea3h (23.2 m) | 0.915 | 0.277 | 0.833 | 0.703 | 0.547 | 0.799 | 0.234 | 0.006 |
| healpix (24.8 m) | 0.995 | 0.300 | 0.823 | 0.761 | 0.446 | 0.825 | 0.217 | 0.23 |
| rhealpix (17.3 m) | 0.939 | 0.261 | 0.874 | 0.763 | 0.549 | 0.850 | 0.267 | **0.36** |

Local-inertial: on a real event the mesh family matters little at matched
cell size — storage 0.92–0.97 of the benchmark's, level rms 0.25–0.30 m, CSI
0.82–0.85 on every ~22 m mesh. The finer H3 and rHEALPix meshes have the
best extents (CSI 0.88); A5 the lowest level rms, by a few centimetres. All
of them are smooth (extrema ≤ 0.008, storage rising throughout).

SGS stores 0.70–0.82 of the benchmark's volume, as on A5 (§15.14: it is
local-inertial on the 5 m geometry, and close to the 5 m local-inertial
run). **But SGS checkerboards on most of the other meshes** — a finding of
this case, recorded as an expected failure (`xfail`) until it is fixed:
smooth on A5, ISEA3H and the jittered Voronoi mesh; a checkerboard over
8–23% of the wet cells on hexagons, centroidal Voronoi, H3, IGEO7 and
HEALPix; and over 36–42% on rHEALPix and the rotated grid, where the stored
volume also oscillates on the rising limb (by 36% of the peak storage on
rHEALPix, frame to frame) and Δt and the maximum depth jump (to 11 m, against
7 m for local-inertial). The rotated grid's SGS storage (0.96) is that
noise, not accuracy. On the smooth meshes SGS stays within 0.16–0.23 m rms
of the 5 m local-inertial run; on the noisy ones up to 0.56 m. The pattern
follows neither regularity nor orthogonality (ISEA3H's hexagons are smooth,
the planar ones are not); it is open (PROJECT_STATE), for the SGS work.

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

