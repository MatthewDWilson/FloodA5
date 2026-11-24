# Globe viewers: CesiumJS and Giro3D

FloodA5 has the same web globe viewer on two rendering engines, to choose
between them on evidence (0.8.4): `--vis cesium` (CesiumJS 1.147) and
`--vis giro3d` (Giro3D 2.0.7 on three.js r180). Usage: docs/USER_GUIDE.md
§6.2; protocol: docs/DATA_FORMATS.md §6.

## What is shared and what is not

Everything except drawing is shared (`visualisation/viewer-common/`): the
server and protocol, the controls and their behaviour (`js/app.js`,
`js/ui.js`), the mesh decoding, triangulation, cell lookup and colour maps,
and the surface model — one mesh per surface, each cell on its own vertices,
per-frame values in float textures, the vertex shader moving each vertex along
the ellipsoid normal to `datum + depth (+ offset)` (`js/surface.js`, GLSL used
by both engines). An engine is a thin adapter (`cesium/js/engine.js`,
`giro3d/js/engine.js`, about 200 lines each): viewer, basemaps, camera,
picking, basemap heights, and its version of `water.js`.

## Measurements

Carlisle, `test/case_runs/carlisle_peak_fine/runs/cart5_swe.h5` (582,012
cells, 112 frames), headless Edge on an RTX 4000 SFF Ada (Direct3D 11),
2026-11-25, `visualisation/viewer-common/tools/smoke.mjs`:

| | CesiumJS 1.147 | Giro3D 2.0.7 |
|---|---|---|
| Page start → first frame (empty browser cache, OSM) | 10.1–10.2 s | 6.6–7.9 s |
| Mesh decode + surfaces | 4.3–7.5 s | 4.0–5.5 s |
| Change of frame, data already in the browser | 6.2 ms | 5.8 ms |
| JS heap after loading | 446 MB | 379 MB |
| Vertical offset estimate, Google tiles (A5 run, 302 cells) | 52.2 m (p10 52.0, p50 52.6) | 52.4 m (p10 51.5, p50 53.0) |

Live (Carlisle A5 SGS, 29,902 cells, `-t 4`): both follow the run one frame
per second (frames 10 → 34–54 over 20 s, each rendered), status "running".

## Differences

| | CesiumJS | Giro3D |
|---|---|---|
| Delivery | One script from cesium.com (`FLOODA5_CESIUM_URL` for a local copy: offline) | 15 ES-module packages from esm.sh, pinned in an import map; offline would need a local mirror or a bundle |
| Google Photorealistic 3D Tiles | Native (`createGooglePhotorealistic3DTileset`), Google Maps logo shown | 3d-tiles-renderer auth plugin (Cesium ion or Google key); the viewer shows the words "Google Maps" and the data providers (the logo is a to-do) |
| Terrain | Cesium World Terrain (ion token), quantized mesh | Terrarium (AWS open data, no key) — in England from EA LiDAR, so in the DEM's datum (offset ≈ 0 at Carlisle) |
| Basemap heights for the offset | `sampleHeightMostDetailed` loads the detail it needs | Rays onto the loaded tiles: the viewer looks down over the domain and waits for the tiles first |
| Custom drawing | `CustomShader` on a glTF model built in the browser | Plain three.js mesh + `ShaderMaterial`: simpler, full three.js available |
| Licence / governance | Apache-2.0; Cesium (Bentley); ion services commercial | MIT; Oslandia, OSGeo incubation; three.js ecosystem |
| Maturity of the globe | Long-standing | Globe mode since 2024; two major versions in 2026 (3.0 beta) |

## Problems found while building them

Giro3D / three.js / 3d-tiles-renderer:
- **`Ellipsoid.toGeodetic` in Giro3D 2.0.7 is wrong by ~130 m in height** (and
  ~0.0008° in latitude) at mid-latitudes: Bowring's formula with the first
  eccentricity where the second belongs. `toCartesian` is right. The viewer
  uses its own `ecefToGeodetic` (`viewer-common/js/geometry.js`, tested to
  < 1 mm). Worth reporting upstream.
- Giro3D turns off three.js's shader error checks: a shader that does not
  compile draws nothing, silently. (Ours did: three.js declares no fragment
  output for an explicit GLSL3 `ShaderMaterial`.) `water.js` now checks the
  program after the first draw and reports.
- The globe controls raycast the whole scene every frame; a brute-force pass
  over the surfaces' millions of vertices starved tile loading (15 s of a
  15 s profile). The surfaces are not raycastable now.
- The Google auth plugin's "recommended settings" set `errorTarget` 20
  (Giro3D's default is 8): coarser tiles, fewer requests. Kept.

CesiumJS:
- The pre-0.8.3 viewer passed the removed `imageryProvider` option, so CesiumJS
  loaded Bing imagery unseen (against Cesium ion's terms with Google tiles).
- `CustomShader` compiles the vertex and fragment texts separately: shared
  GLSL helpers have to be in both.

Both (server side): the first request after the page was busy for seconds
building a large mesh is sometimes aborted by the browser (on a reused
keep-alive connection to the HTTP.jl server); the viewers retry network
errors with back-off. The root cause is not yet known.

## Choosing

Both viewers meet the requirement — water at its height over photorealistic
3D tiles, live and replay, 500k+ cells — with the same controls and similar
speed. CesiumJS is the more complete globe out of the box (tile handling,
height sampling, attribution, offline copy). Giro3D starts faster, uses less
memory, gives plain three.js for custom drawing, and has no vendor in the
middle, at the price of more integration work (three bugs above) and a less
settled API. The shared code keeps either choice cheap: dropping one engine
deletes one directory.
