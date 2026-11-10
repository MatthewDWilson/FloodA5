# Planar slope / embankment terrain

`generate_planar_embankment_dem.py` builds a west→east planar slope with an
optional N–S embankment and a notch through it, and writes into `--out`:
`planar_embankment_dem.tif` (EPSG:4326), `planar_embankment_params.json` (key
elevations, injection point, extent) and `planar_embankment_aoi.geojson`.

```
West (high)                         East (low)
      │           ▓▓▓▓▓            │
      │           ▓▓ ▓▓  ← notch   │   embankment: --emb-height, --emb-width,
 source ─→        ▓▓ ▓▓            │   at --emb-pos (fraction of the width)
(x = 10%)         ▓▓▓▓▓            │
                 (x = 50%)
```

Generator semantics worth knowing:
- The embankment top is flat at `z(x_emb) + emb_height` (the crest); the notch
  floor (the sill) is `z(x_emb)`, flat along the embankment.
- Embankment and notch are the pixels whose centres lie within ±width/2, so
  widths are effectively whole pixels: with the default 5 m pixel a 3 m notch
  is empty and a 5 m notch is 10 m wide. Choose `--pixel-size` so that
  width/2 does not fall just short of a pixel centre.
- `--emb-height 0` leaves a flat strip of `--emb-width` at `z(x_emb)`; add
  `--emb-width 0` for a pure plane.
- Defaults: 4 × 2 km, 0.1% slope, 5 m pixels, embankment 5 m high.

## Test cases

The planar tests are registry cases (`cases/registry/`, `docs/TEST_CASES.md`
§2). Each generates its own DEM and meshes in `test/case_runs/<case>/`:

```bash
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl show planar_embankment   # commands, to copy
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_embankment -P 2 -T 2
FLOODA5_NO_VIS=1 julia --project=. cases/run.jl run planar_slope -P 3 -T 2
```

- `planar_embankment` (§2a): a 4 m notch through a 2 m embankment on a 0.5%
  slope; SGS passes water near the sill, standard flow only at the crest it
  sees. It replaces `test_planar_embankment.jl`.
- `planar_slope` (§2b): the N/S volume-symmetry benchmark on a pure 0.1%
  plane, A5 res 16, 18 and Cartesian 25 m.

## Files in this directory

- `planar_embankment_params.json`, `planar_embankment_aoi.geojson` (and the
  git-ignored `planar_embankment_dem.tif`): the legacy fixture, generated with
  `--emb-height 0` (a 20 m flat strip at x = 2 km), still used by older run
  scripts (`test/run_wp1_matrix.sh`, `test/run_FOSS4G2026_simulations.ps1`, …).
  The cases do not use it.
- `planar_embankment_dem_OLDwithWALL.tif`: an earlier DEM with the embankment.
- `planar_slope_tests/`: results and animations of the earlier symmetry experiments.
