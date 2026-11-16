# SGS Phase 0 audit: hypsometry, sub-basins, time step (2026-11-16)

Plan: `projectinstructions/FloodA5_SGS_Conveyance_Fix_Plan_2026-11.md`, Phase 0
items 2–4. Measurement only; no model change. Code: `analysis/sgs_audit.jl`
(pure functions, tested in `test/analysis/test_sgs_audit.jl`) and the CLI
`test/diag_sgs_tables_audit.jl`. Meshes are the 0.6.0 tables (layout 2).

## Summary

- **Sampling (H1) is the main hypsometry error and is large enough to matter
  on Carlisle.** In the lowest 0.5 m the table volume is off by 1–2 % in the
  median cell, 4–11 % at the 95th percentile, over 50 % at worst; the level for
  a 1–20 cm layer is off by 2 mm (median), 2 cm (p95), 0.3–0.7 m (worst).
  4096 candidates instead of 512 cut these about four-fold.
- **Extremes (H2) matter.** Carlisle: the table z_min is more than 5 cm above
  the terrain minimum in 11 % of cells (10 cm: 5 %; worst 2.7 m). 13 % of
  cells have a face point more than 5 cm below their own z_min (10 cm: 6 %;
  worst 3.15 m), so C3 is real. A5 strip: 2.1 % of cells miss the channel
  altogether (z_min 0.25–0.99 m too high).
- **Uniform levels (H3) do not matter much.** With sampling error removed,
  100 uniform levels give under 1 % volume error (p95) at every level tested
  and under 1 mm level error (p95) on all terrains, although the spacing
  exceeds 1 cm in 39 % of Carlisle cells. The worst cells reach 26 mm.
- **S1 ("storage within 0.3 %") holds for the Cartesian strip and the median
  A5 cell, not for A5 cells that clip the channel**: p95 3–7 %, worst 20–33 %
  of the (small) volume they hold; level error up to 0.19 m for a 1 cm layer.
- **Sub-basins.** Carlisle (5 m DEM, a lower bound): 7.1 % of cells hold a
  second basin with over 5 % of the area and deeper than 0.25 m (7.9 % of cells
  wet at the peak); 1.8 % at 20 % and 0.5 m (2.3 % wet). For a
  single-storage rule: most cells have one basin, 5–8 % need a choice, 1–2 %
  hold two comparable pools and are the D2 candidates.
- **The seed must not be the largest-area basin.** On the A5 strip the channel
  is the smaller-area basin beside a 2 cm floodplain shelf in 9.8 % of cells;
  an area seed would hold the channel shut to bank level there. On Carlisle
  the lowest point is outside the largest-area basin in 13.6 % of cells.
- **Time step (C1): smaller exposure than the plan expected at Carlisle.**
  With the true depth the largest Courant number at the peak is 0.95 against
  0.75 with the cap; the true depth needs 1.13 s against the 1.43 s taken.
  No cell exceeds 1.0. Cells with large relief (cap up to 11.3 m) set the step.

## Method and assumptions

- Reference: the DEM's bilinear surface (`A5Grid._bilinear`, the surface the
  model samples, with its CRS transform) on a regular grid inside each
  polygon, spacing h = min(pixel/8, √A/64): 0.35 m on Carlisle, 0.25 m on the
  embankment, 0.0625 m on the strip. Pixel/2 was not enough on the strip (the
  channel walls are ramps one pixel wide; 2–4 % error). Each point stands for
  A/n, as in the table build.
- Reference accuracy, h against h/2 on 200 cells, volume in the lowest 0.5 m:
  Carlisle p95 0.4–1.2 % (worst 5.6 %), embankment p95 < 0.1 %, strip p95
  < 0.5 %. On the strip the reference agrees with the analytic terrain
  (0.05 m grid) to 0.5 % (p95).
- Table read as the model does: `FA.sgs_table_volume` / `sgs_table_wse` (the
  arithmetic of `wse_from_volume`). Relative volume errors use
  max(V, 1 cm × cell area). Levels are measured from the table's own z_min.
- Check of the script: tables rebuilt from the model's 512 Halton candidates
  equal the stored ones exactly (0.0 difference) on all four meshes.
- Cells the DEM covers only in part are left out: 947 of 29 902 (Carlisle),
  327 of 2 189 (embankment), none on the strip.
- **Basin definition (an assumption).** As the level rises, ground below it
  forms 4-connected regions on the working grid. When two meet, the one with
  the smaller area ends and is recorded with its area share and its depth
  (merge level − its lowest point). The region that survives every merge is
  the main basin. Counted: basins that merge directly with the main one.
  "Two-sided" also requires the main basin to be deeper than the threshold
  below the merge level (two pools, not a pool beside a shelf). A region
  shallower than 1 mm never outlasts one with depth (plateaus are not basins).
- Carlisle's 5 m DEM gives 4.4 pixels across a cell. The working grid
  (57 × 105 points) only interpolates them, so walls, kerbs and ditches under
  about 10 m are absent and the basin counts are a lower bound.

## A. Hypsometry

Stored table against the reference, all audited cells. Volume: |error| in %
of max(V, 1 cm·A); level: |η_table(V) − η_ref(V)| in mm. Median / p95 / max.

| | Carlisle A5 r18 | Embankment A5 r18 | Strip A5 r18 | Strip Cartesian 20 m |
|---|---|---|---|---|
| V at z_min + 0.05 m | 0.97 / 4.3 / 74 | 1.33 / 1.5 / 7.2 | 0.22 / 2.7 / 21 | 0.05 / 0.43 / 0.47 |
| V at z_min + 0.10 m | 1.73 / 7.1 / 51 | 0.39 / 0.89 / 4.4 | 0.10 / 4.6 / 22 | 0.02 / 0.25 / 1.9 |
| V at z_min + 0.25 m | 1.47 / 11.3 / 52 | 0.10 / 0.21 / 8.8 | 0.04 / 6.6 / 25 | 0.01 / 0.06 / 1.8 |
| V at z_min + 0.50 m | 0.70 / 8.7 / 57 | 0.04 / 0.09 / 13 | 0.02 / 7.3 / 24 | 0.00 / 0.14 / 0.40 |
| V at z_min + 1.00 m | 0.28 / 4.5 / 72 | 0.02 / 0.04 / 41 | 0.01 / 7.2 / 33 | 0.00 / 0.15 / 0.17 |
| V at z_max | 0.56 / 1.8 / 23 | 0.45 / 0.72 / 48 | 0.60 / 6.9 / 22 | 0.20 / 0.20 / 0.22 |
| level, mean depth 1 cm | 1.9 / 23.6 / 730 | 0.32 / 0.45 / 254 | 0.09 / 15.7 / 192 | 0.02 / 0.22 / 2.6 |
| level, mean depth 5 cm | 1.7 / 21.0 / 449 | 0.19 / 0.40 / 49 | 0.08 / 15.4 / 76 | 0.02 / 0.04 / 1.0 |
| level, mean depth 20 cm | 1.7 / 16.8 / 302 | 0.19 / 0.40 / 49 | 0.08 / 6.4 / 44 | 0.02 / 0.27 / 1.0 |

The errors have almost no bias (median signed error under 0.2 % and 0.2 mm,
except −0.7 % at z_min + 0.05 m on the embankment). Cells wet at the Carlisle peak are the same as all cells (p95
4.3–9.5 % volume, 17–22 mm level). Against the analytic strip terrain the
table gives the same numbers as against the reference (A5: p95 2.7–7.3 %).

Extremes, sample counts, relief and faces:

| | Carlisle | Embankment | Strip A5 | Strip Cartesian |
|---|---|---|---|---|
| Halton points inside (of 512): median (min–max) | 353 (306–357) | 305 (303–355) | 353 (303–355) | 512 |
| table z_min − terrain z_min: median / p95 / max (m) | 0.008 / 0.102 / 2.71 | 0.001 / 0.003 / 0.37 | 0.000 / 0.009 / 0.99 | 0 / 0 / 0.003 |
| z_min too high by > 1 / 5 / 10 / 25 cm (% of cells) | 45.1 / 11.4 / 5.2 / 1.6 | 0.2 / 0.1 / 0.1 / 0.1 | 4.8 / 2.4 / 2.3 / 2.1 | 0 |
| table z_max − terrain z_max: median / p5 / min (m) | −0.010 / −0.115 / −1.90 | 0.000 / −0.002 / −0.35 | 0.000 / −0.002 / −0.009 | 0 |
| z_max too low by > 1 / 5 / 10 / 25 cm (% of cells) | 51.3 / 14.1 / 6.0 / 1.7 | 0.2 / 0.1 / 0.1 / 0.1 | 0 | 0 |
| cells with a face point below their z_min by > 1 / 5 / 10 / 25 cm (%) | 43.3 / 13.1 / 6.2 / 2.0 | 0.2 / 0.1 / 0.1 / 0.1 | 5.0 / 2.6 / 2.4 / 2.1 | 0 |
| face sides below their cell's z_min by > 1 / 5 / 10 cm (%) | 16.0 / 5.1 / 2.4 | 0.1 / 0 / 0 | 2.1 / 1.2 / 1.2 | 0 |
| largest face deficit (m) | 3.15 | 0.50 | 0.99 | 0.003 |
| relief z_max − z_min: median / p95 / max (m) | 0.74 / 4.59 / 21.5 | 0.12 / 0.12 / 2.2 | 0.02 / 1.02 / 1.49 | 0.02 / 1.02 / 1.42 |
| level spacing: median / p95 / max (mm) | 7.5 / 46 / 217 | 1.2 / 1.2 / 22 | 0.24 / 10 / 15 | 0.20 / 10 / 14 |
| spacing above 1 / 2 / 5 cm (% of cells) | 38.9 / 17.6 / 3.7 | 2.2 / 1.2 / 0 | 40.5 / 0 / 0 | 33.3 / 0 / 0 |

Level interpolation alone (a 100-level table built from the dense samples,
against the dense curve). Median / p95 / max.

| | Carlisle | Embankment | Strip A5 | Strip Cartesian |
|---|---|---|---|---|
| V at z_min + 0.05 m (%) | 0.03 / 0.88 / 78 | 0.01 / 0.02 / 7.0 | 0 / 0.01 / 0.11 | 0 / 0.01 / 0.12 |
| V at z_min + 0.25 m (%) | 0.02 / 0.28 / 12 | 0 / 0 / 0.02 | 0 / 0.01 / 0.08 | 0 / 0 / 0.08 |
| V at z_min + 1.00 m (%) | 0 / 0.03 / 2.1 | 0 / 0 / 0.04 | 0 / 1.3 / 5.2 | 0 / 0.05 / 0.05 |
| level, mean depth 1 cm (mm) | 0.04 / 0.84 / 26 | 0 / 0.01 / 1.5 | 0 / 0.06 / 2.3 | 0 / 0.01 / 0.1 |
| level, mean depth 20 cm (mm) | 0 / 0.13 / 5.1 | 0 / 0 / 0.1 | 0 / 0.01 / 1.4 | 0 / 0 / 0.1 |

Sampling convergence: the same cells rebuilt with 4096 candidates (about
2450–2850 inside), p95 of |error|, 512 → 4096.

| | Carlisle (491 cells) | Embankment (437) | Strip A5 (500) | Strip Cartesian (500) |
|---|---|---|---|---|
| V at z_min + 0.10 m (%) | 7.4 → 2.0 | 0.89 → 0.24 | 4.6 → 1.1 | 0.25 → 0.22 |
| V at z_min + 0.25 m (%) | 9.7 → 3.1 | 0.21 → 0.06 | 6.5 → 1.9 | 0.06 → 0.11 |
| V at z_min + 0.50 m (%) | 8.2 → 2.0 | 0.09 → 0.03 | 7.0 → 2.3 | 0.14 → 0.13 |
| level, mean depth 1 cm (mm) | 20.1 → 5.7 | 0.45 → 0.13 | 12.3 → 3.6 | 0.22 → 0.21 |
| level, mean depth 20 cm (mm) | 17.3 → 4.4 | 0.40 → 0.11 | 6.9 → 1.2 | 0.27 → 1.1 |
| z_min − terrain z_min, p95 (m) | 0.095 → 0.015 | 0.003 → 0.000 | 0.006 → 0.002 | 0 |
| z_min − terrain z_min, max (m) | 0.64 → 0.09 | 0.003 → 0.000 | 0.99 → 0.50 | 0 |

4096 points still miss a channel sliver by 0.50 m on the A5 strip: more
points reduce H2 but do not remove it. A pixel-based working grid shared by
the cell and its faces (plan §3.4) would.

Side finding: `build_sgs_tables!` converts a geographic DEM's pixel to metres
without cos φ, so faces are sampled every 0.40 m on the strip (0.5 m pixels)
and 1.59 m on the embankment (2 m pixels), not at half a pixel.

## B. Sub-basins

Share of cells (%) with at least one basin, beside the main one, above both
thresholds. One-sided as the task defines it / two-sided.

Carlisle, all 28 955 audited cells:

| area > | depth > 0.10 m | 0.25 m | 0.50 m |
|---|---|---|---|
| 1 % | 22.4 / 20.0 | 9.2 / 6.9 | 4.5 / 2.8 |
| 5 % | 16.0 / 14.1 | 7.1 / 5.2 | 3.9 / 2.3 |
| 20 % | 5.1 / 4.6 | 2.9 / 2.2 | 1.8 / 1.1 |

Carlisle, 16 760 audited cells wet at the peak (last frame of `a5_sgs.h5`,
t = 199 800 s, depth > 10⁻⁴ m):

| area > | depth > 0.10 m | 0.25 m | 0.50 m |
|---|---|---|---|
| 1 % | 22.6 / 20.2 | 9.8 / 7.5 | 5.2 / 3.2 |
| 5 % | 17.3 / 15.3 | 7.9 / 5.8 | 4.6 / 2.7 |
| 20 % | 6.0 / 5.5 | 3.5 / 2.6 | 2.3 / 1.3 |

At (1 %, 0.10 m), 2.2 % of Carlisle cells have two or more such basins and
0.3 % three or more.

Where the cell's lowest point lies (matters for the seed rule):

| | Carlisle all | Carlisle wet | Strip A5 | Strip Cartesian | Embankment |
|---|---|---|---|---|---|
| lowest point outside the largest-area basin (its basin > 1 % area, > 2 cm deep), % of cells | 13.6 | 15.0 | 9.8 | 0 | 0.05 |
| … and that basin passes (5 %, 0.25 m) | 3.9 | 4.3 | 5.2 | 0 | 0 |
| … and passes (20 %, 0.50 m) | 1.2 | 1.5 | 0 | 0 | 0 |

A5 strip: 9.8 % of cells show a basin over 1 % of the area and over 0.5 m
deep, but only 0.6 % two-sided at 0.10 m and none at 0.5 m. These are cells
where the channel covers less area than the floodplain shelf (2 cm of
relief). There the lowest point is always in the channel.

Embankment (x = 980–1020 m, crest 7.0 m, notch y = 248–252 m). The widest
cell is 24.3 m in x against 40 m, so **no cell spans the embankment** and
none is split by it into two pools:

| class (x-range against the embankment) | cells | on the notch | relief median (m) | second basin (1 %, 0.10 m) |
|---|---|---|---|---|
| spans both faces | 0 | – | – | – |
| straddles the west face | 19 | 2 | 1.89 | 1 |
| straddles the east face | 21 | 2 | 2.14 | 0 |
| wholly on the crest | 19 | 1 | 0.00 | 0 |

Notch cells: #400 (x 962.2–986.5), #554 (971.6–995.9), #2035 (984.2–1003.7),
#2151 (995.9–1020.2), #612 (1017.9–1037.4). The one cell with a second basin
is #554 (y 216.6–251.1), which touches the notch at a corner. Whole mesh:
0.05 % of cells. The present case therefore does not test a barrier inside a
cell; the narrow-wall variant of Phase 0 item 1 is needed for that.

Implications for one storage per cell:

1. A lowest-point seed is right for a channel beside a shelf (all the A5
   strip cases) and wrong only where the lowest point is in a minor hollow:
   at most 13.6 % of Carlisle cells, 3.9 % if the hollow has to pass (5 %,
   0.25 m). An area seed is wrong on 9.8 % of A5 strip cells and would dam the
   channel. Neither rule alone is safe; the build-time check of faces whose
   sill is raised between two through-flowing regions (plan §3.5) is needed.
2. Two real pools (two-sided) at (5 %, 0.25 m): 5.2 % of Carlisle cells, 5.8 %
   of wet ones. At (20 %, 0.50 m): 1.1 % and 1.3 %. D1 has to choose in the
   former; the latter are where a wrong choice costs most and are the
   candidates for D2. On a finer DEM these shares will rise.

## C. Time-step exposure (C1)

Last frame of the `carlisle_peak` runs (t = 199 800 s), wet cells
(depth > 10⁻⁴ m), Courant = Δt √(g h) / ℓ with ℓ = 4A/P (median 21.58 m) and
the run's mean step. SGS: 131 053 steps, Δt = 1.5246 s. Standard: 97 439
steps, Δt = 2.0505 s. Both ran with `--dt-max 10`, Courant 0.7, local CFL.

| run, depth measure | wet cells | depth median / p95 / max (m) | Courant median / p95 / p99 / max | > 0.7 (%) | > 1.0 (%) | step required at 0.7 (s) |
|---|---|---|---|---|---|---|
| SGS, capped min(η − z_min, z_max − z_min) (the model's) | 16 942 | 0.50 / 4.02 / 11.30 | 0.16 / 0.44 / 0.52 / 0.75 | 0.03 | 0 | 1.429 |
| SGS, η − z_min | 16 945 | 3.31 / 9.11 / 18.21 | 0.40 / 0.67 / 0.78 / 0.95 | 2.87 | 0 | 1.127 |
| SGS, V / A_wet | 16 945 | 2.81 / 8.24 / 18.14 | 0.37 / 0.64 / 0.74 / 0.95 | 1.94 | 0 | 1.128 |
| standard, water depth | 9 977 | 1.40 / 5.58 / 7.13 | 0.35 / 0.70 / 0.76 / 0.80 | 5.35 | 0 | 1.799 |

- The capped rule reproduces the run: 1.429 s required, 1.43 s in the log at
  the last steps. At that step (not the mean) the true-depth maximum is 0.89.
- 84.0 % of wet cells have water above z_max (plan S5: 85 %).
- The true depth needs a step 21 % shorter than the one taken (1.13 s
  against 1.43 s), not "several times". The step is set by a cell with
  11.3 m of relief, close to the deepest water (18.2 m).
- The exposure is still real: the median wet cell's wave speed is
  understated 2.5-fold (0.16 against 0.40), and on low-relief terrain (the
  strip: 2 cm of relief under metres of water) nothing bounds the step. C1
  should be fixed with the conveyance as planned; at Carlisle it will cost
  about a quarter more steps.

## Against the plan's expectations

- H1 "small at res 18 on a 0.5–2 m DEM": true for the embankment and the
  Cartesian strip, not for A5 strip cells that clip the 5 m channel, and not
  on Carlisle's 5 m DEM.
- H3: smaller than H1 and H2 everywhere; changing the level spacing alone
  would gain little.
- C1 at Carlisle: a factor 1.27, not several.
- S1: see the summary.

## Not verified

- Sub-basins on a finer Carlisle DEM (none available here).
- The Carlisle reference is converged to about 1 % (p95) in volume; table
  errors below that are not resolved.
- Courant numbers use the mean step and one frame; no other frames checked.
- The strip's time-step exposure and face-based depth A_f/T_f (plan §3.3)
  were not computed.

## Commands

From the repository root; outputs (JSON and log) in `test/wp1_runs/sgs_audit/`
(git-ignored). About 3 minutes each on one thread.

```bash
export FLOODA5_NO_VIS=1
O=test/wp1_runs/sgs_audit
SD=D:/FloodA5-wpc7/test/case_runs/sgs_strip
PE=D:/FloodA5/test/case_runs/planar_embankment
julia --project=. test/diag_sgs_tables_audit.jl D:/FloodA5/test/case_runs/carlisle/carlisle_a5_res18_sgs.parquet \
    D:/FloodA5/test/carlisle/Carlisle_LiDAR_5m_mean.tif $O/carlisle.json \
    --run D:/FloodA5/test/case_runs/carlisle_peak/runs/a5_sgs.h5 \
    --compare-run D:/FloodA5/test/case_runs/carlisle_peak/runs/a5_li.h5 > $O/carlisle.log
julia --project=. test/diag_sgs_tables_audit.jl $PE/a5_res18_sgs.parquet $PE/planar_embankment_dem.tif \
    $O/embankment.json --embankment-params $PE/planar_embankment_params.json > $O/embankment.log
julia --project=. test/diag_sgs_tables_audit.jl $SD/a5r18_sgs.parquet $SD/sgs_strip_dem.tif \
    $O/strip_a5r18.json --strip-params $SD/sgs_strip_params.json > $O/strip_a5r18.log
julia --project=. test/diag_sgs_tables_audit.jl $SD/cart20_sgs.parquet $SD/sgs_strip_dem.tif \
    $O/strip_cart20.json --strip-params $SD/sgs_strip_params.json > $O/strip_cart20.log
julia --project=. test/test_analysis.jl      # includes test/analysis/test_sgs_audit.jl
```

Defaults used: `--subset 500 --halton-hi 4096 --seed 1`. The run commands of
the two Carlisle outputs are in the files (`FO.run_info(p)["argv"]`).
