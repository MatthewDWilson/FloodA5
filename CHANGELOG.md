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

## [0.8.15] — 2026-11-26

Output files on the web in the globe viewers; shared palettes in the
repository (branch `wpv-remote`; PATCH: viewer only, results unchanged).

- Files on a web server: `?file=URL` (absolute, or relative to the page:
  a static host with the viewer and the outputs) or *Open…* → *or from a
  link*. Read in the browser like a local file (0.8.10), through HTTP range
  requests where the server serves them: `js/local/rangefile.js`
  (`RangeCache`: 64 KiB blocks fetched on demand by synchronous range
  requests in the worker, adjacent blocks in one request, an LRU cache of at
  most 256 MiB) behind a read-only file in h5wasm's file system; otherwise
  the whole file is downloaded (with progress in the log). The 400 MB
  Carlisle A5 run opened from another origin with 17 MB fetched (1.7 s);
  the maximum depth over all 112 frames brought it to 33 MB. Errors name
  the cause (HTTP status; CORS/offline). The page's address then carries
  `?file=…` (a shareable link to the run); a failed link stays in the
  *Open* panel to be corrected. `/api/info` adds `remote`.
- The top bar's *Open…* shows the *Open* panel (file or link) instead of
  the file chooser directly.
- `visualisation/viewer-common/palettes.json` is in the repository (no
  longer git-ignored; the template `palettes.example.json` became it, with
  its two palettes); without a server the viewer reads it directly, so a
  static host offers the *Shared* palettes too.
- Fixes from testing 0.8.14:
  - CesiumJS (Firefox): "e._va is undefined" when opening a file — from
    CesiumJS's dynamic environment map, which every glTF model computes by
    default (image-based lighting) and whose clean-up can run twice when
    models are replaced; the surfaces are unlit: it is now off for them
    (`environmentMapOptions: { enabled: false }`), which also saves the work.
  - Terrain basemaps: water hidden at a few hundred metres — the depth bias
    (0.8.6) came in only from 800 m, for 3D tiles; terrain tiles are coarse
    at middle distances too: there it comes in from 30 m, full from 600 m
    (`surface.js` `DEPTH_BIAS_RANGE`, uniforms biasNear/biasFar, set per
    basemap).
  - Vertical offset: an estimate beyond ±120 m (from tiles too coarse at the
    view, e.g. a basemap switched while zoomed out) is refused, not saved;
    such a value saved before is ignored.
  - Giro3D: a sky — a gradient near the ground (Giro3D's physical SkyDome is
    almost white at the low angles of an oblique view, and takes the
    camera's up, which these controls keep on the Earth's axis), the
    Atmosphere entity's halo from space only (inside its shell it whitens
    the view; its veil washes out light basemaps and is off).
- Tests: `viewer-common/test/rangefile.test.js` (reads across blocks and at
  the end, one request for adjacent blocks, the request size limit, cache
  hits, the cache limit and LRU order, short answers and failures, file
  names from URLs, which links may be opened — http(s) and relative only).
  Browser smoke: a cross-origin range server, a server without ranges (full
  download), a 404, a static host with a relative link and the shared
  palettes.

## [0.8.14] — 2026-11-26

Smooth water and bed surfaces in the globe viewers (branch `wpv-smooth`;
PATCH: viewer only, results unchanged).

- *Smooth surface* (Display) and *Smooth* (Bed surface), off by default:
  each cell's ring corners are drawn at the mean level of the wet cells
  sharing them (the bed: the mean elevation), each centre at the cell's own
  level, so the surface follows the slope instead of stepping; colours
  blend (not with classed palettes), and the surface is shaded by its slope
  (exaggerated 20× for water, 3× for the bed, for the shading only; lit
  from the north-west). For viewing — the flat cells remain the exact view.
- How: `geometry.js` `cornerIndex` finds the shared corners by sorting
  (no hash map: affordable for millions of cells), once, when smoothing is
  first on; `cornerValues` / `bedCornerValues` give each corner's means
  (one RGBA32F texel per corner: level, bed, depth, colour value) for every
  frame; the shaders (`surface.js` `fa5Corner`, `fa5LevelC`,
  `fa5ColorValueC`, `fa5Shade`) read them for ring vertices. Giro3D adds a
  `corner` attribute; CesiumJS rebuilds the surface once with a
  `_FA5CORNER` glTF attribute (`buildGlb(surface, corner)`).
- Cost, measured: the 30 000-cell A5 Carlisle mesh — corners 57 ms, nothing
  noticeable per frame; the 582 000-cell 5 m mesh — corners 0.8–0.9 s
  (CesiumJS: 1.9 s with the rebuild), a frame 90–125 ms instead of 30–75 ms
  (23 ms of it the corner means on the CPU).
- About: a *Smooth surface* note (`io/variables.toml`).
- Tests: `viewer-common/test/geometry.test.js` (corners shared and welded
  within ~1 cm, centres −1, corner → cells; means over wet cells, missing
  colour values left out, dry corners zero; the bed's means; the glTF
  attribute). Browser smoke on both engines (flat vs smooth, water and bed,
  frames stepping).

## [0.8.13] — 2026-11-26

Statistics over the run written by the model (branch `wpv-stats`; PATCH:
results unchanged — frames bit-identical with and without them).

- The output gets `/statistics` (DATA_FORMATS §2): maximum depth and its
  time, maximum speed and D·V, the ARR class, time of first flooding and
  duration, gathered while the model runs and rewritten with every output
  frame. `--stats-interval S` sets how often (default: at the output frames,
  at no extra cost; `0` = every step, catching peaks between frames — cell
  speeds are then computed, and a GPU state copied to the host, that often);
  `--flood-depth M` the threshold (default 0.1 m). Hazard rating and degree
  maxima depend on the viewer's parameters and are not stored.
- One implementation: `analysis/hazard.jl`'s `run_statistics` is now a loop
  over `StatsAccumulator` (`accumulate_frame!`, `accumulated_statistics`;
  exported by FA), which the model uses during the run, on the same Float32
  values the viewers read — at the default cadence the model's statistics
  equal those derived from the frames bit for bit (tested), and
  `run_statistics` gives the same values as before (the local-file fixture
  is unchanged).
- Globe viewers: a replay (server, and a local file in the browser) serves
  those statistics from the file — on the 400 MB Carlisle SGS run the first
  *Maximum depth* took 2.5 s, now 0.4 s, with identical values; the hazard
  maxima are still derived. The *Run* tab says how the statistics were
  gathered.
- `FA.write_run_statistics!(path)` adds `/statistics` to an earlier file
  from its frames (18 s for that run; +0.8 MB).
- `FO.write_statistics!`, `FO.read_statistics`; `run_info` gains `statistics`.
- Tests: `test/test_run_statistics.jl` (new, in CI): default cadence equals
  the frames' statistics bit for bit; `--stats-interval 0` and
  `--flood-depth` change no frame (CPU and the ka-cpu device backend) and
  give maxima ≥ the frames'; the accumulator part-way and in one pass;
  `write_run_statistics!`; files without statistics. The local-file fixture
  has a `/statistics` group with values of its own (the answers show which
  source was used).

## [0.8.12] — 2026-11-26

Globe viewer fixes from testing 0.8.8–0.8.9 (branch `wpv-feedback`; PATCH:
viewer only, results unchanged).

- The right-click menu could be covered: the browser's own menu is now
  suppressed over the whole map (not only the engine's canvas, which other
  layers can cover), the right click is taken anywhere on the map, the
  hover tooltip and palette menus close while the menu is open, and the
  menu is above every other panel.
- *Parameters* (was *Hazard parameters*) shows only what the layer shown —
  and the plot's variable, while the plot is open — depends on: land use
  for the hazard rating, the person and Chen et al. calibration for the
  hazard degree, the flooded threshold for the time of first flooding and
  the duration; hidden for the other layers. The heading names the layers.
- New statistic *Time of maximum depth* (`time_max_water_depth`: the first
  output frame with the cell's largest depth; `run_statistics`, the
  server's `/api/stats` and the browser's port alike; `io/variables.toml`).
- A shorter side panel: the bed's controls show only while the bed is
  shown (its checkbox is in the heading); basemap and vertical offset in one
  section, the offset's explanation in a tooltip (ⓘ); tighter spacing; the
  panel ends higher, clear of the logos and credits.
- Tests: `test/analysis/test_hazard.jl` (time of the maximum),
  `test_vis_server.jl` (depth-only statistics), the local-file fixture
  regenerated (`test/vis_local_fixture.jl --write`), `local.test.js`.

## [0.8.11] — 2026-11-26

Palettes in the globe viewers: scientific colour maps, user-defined palettes
and a group's shared palettes (branch `wpv-palettes`; PATCH: viewer only,
results unchanged). Step 5 of the viewer plan
(`projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md` §4).

- Scientific colour maps as fixed 256-entry tables in
  `visualisation/viewer-common/colormaps.json` (35 KB, served with the
  viewer): cmocean *deep, dense, turbid, haline, speed, amp, tempo, thermal,
  ice, rain* and Crameri's Scientific colour maps 8.0.0 *batlow, oslo,
  devon, lajolla, davos, imola, lapaz, hawaii, roma, vik*. Generated by
  `viewer-common/tools/make_palettes.jl` from the data files of
  ColorSchemes.jl 3.31.0 (in the Manifest through Makie; read, not loaded —
  no dependency added). Both families are MIT-licensed (cmocean © 2015
  Kristen M. Thyng; Scientific colour maps © 2023 Fabio Crameri — checked in
  the package's `data/cmoceanlicense.txt` and `data/scientificcolormaps.txt`);
  the JSON carries each family's source, version, reference and licence
  text.
- The palette menus (water and bed) are grouped: built-in, cmocean,
  Crameri, shared, yours. The built-in palettes keep their names, so
  choices stored by earlier versions (`flooda5.palette.<var>`,
  `flooda5.palette.bed`) still apply.
- User-defined palettes (*Edit palettes…* in a palette menu): colour stops
  (colour and position, 2–32), add/remove, reverse, even spacing, start
  from any palette, live preview; colours interpolated in OKLab; saved in
  the browser, usable for the water and the bed; export to / import from a
  JSON file in the browser (nothing uploaded), checked strictly (types,
  2–32 stops, positions 0–1 sorted, `#rrggbb`, names 1–40 characters,
  unique, ≤ 64 palettes, ≤ 256 KiB; names shown only as text). Editing a
  palette in use recolours the map (each saved table gets a new key, so
  both engines upload it).
- Optional `visualisation/viewer-common/palettes.json` (git-ignored;
  template `palettes.example.json`): a group's shared palettes, validated
  by the server with the same rules (`VisualisationServer.validate_palettes`,
  `shared_palettes`) and offered read-only by `GET /api/palettes`; an
  invalid file is ignored with a warning. Never served as a file (nor are
  `config.json` and `keys.json`).
- Code: `js/palettes.js` (pure: OKLab, stops → table, file checks;
  `test/palettes.test.js`), `js/palette-editor.js` (loading, menu, editor;
  `window.__flooda5api.palettes` for scripting and the smoke test),
  `css/palettes.css`; `js/colormaps.js` registers run-time tables
  (`registerColormap`, `colormapKey`). Tests: `test/test_vis_server.jl`
  (palettes.json absent, valid, invalid, oversized; config/keys/palettes
  files not served). Docs: USER_GUIDE §6.2 *Palettes*, DATA_FORMATS §6
  *Palette files*.
## [0.8.10] — 2026-11-26

Output files from the user's own computer in the globe viewers, read in the
browser (branch `wpv-local`; PATCH: no solver change, results unchanged).
Step 4 of the viewer plan (`projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md`).

- *Open file…* (top bar), or a file dropped on the page, opens a FloodA5
  output (format 2) in this page. h5wasm 0.10.3 (HDF5 in WebAssembly) in a
  module worker reads it through Emscripten's WORKERFS — in place, a frame
  or a cell's series at a time, never uploaded — and
  `viewer-common/js/local/` answers the viewer's `/api` requests as the
  replay server does: `localrun.js` ports `web_replay.jl`, `encode_mesh`,
  `depth_datum`, `run_ranges` and `FO.run_record`, `hazardcalc.js` ports
  `analysis/hazard.jl` (hazard variables, `run_statistics`). Everything a
  replay offers works: variables, hazards, statistics, colour ranges, the
  time-series plot, the *Run* tab. On the Carlisle A5 SGS run (400 MB,
  29 902 cells, 112 frames) every answer equals the server's for the same
  file (binary to the last bit, JSON exactly; the mean cell area within
  rounding); the worker opens it in 0.2–1.3 s, a statistic takes about
  0.5 s. A file that cannot be opened (e.g. format 1, no cell polygons)
  leaves the one open before it in place; the file replaces the run in that
  page only.
- No server needed: the pages work from a static web host (`/api/config`
  missing → `viewer-common/defaults.json`; *About* from
  `viewer-common/about.json`, generated from `io/variables.toml` by
  `viewer-common/tools/make_about.jl`), starting with a prompt to open a
  file. `web_replay.jl` without a file starts a server with no run
  (`web_replay_empty`; config `no_data`) that does the same.
- Configuration: `h5wasm_url` (`FLOODA5_H5WASM_URL`; offline: a copy in
  the git-ignored `viewer-common/vendor/`). The library URLs' defaults
  (CesiumJS, Chart.js, h5wasm) moved to `viewer-common/defaults.json`, read
  by the server and by the viewer without one.
- Plan: the `--dir` fallback (a server-side file list) is not needed —
  WORKERFS reads large files without loading them — and was not built.
- Tests: `test/test_vis_local.jl` (new, in CI): two generated runs (SGS with
  NaNs, non-finite options, open and closed rings; depth-only standard)
  served by the replay server and recorded with their contents in
  `viewer-common/test/fixtures/local_run.json` (`test/vis_local_fixture.jl
  --write`); the committed fixture must still be what the server answers,
  and the JavaScript port must answer a fresh one alike (node);
  `about.json` = `variables.toml`; `defaults.json` = the server's
  defaults; the no-run server. JavaScript: `test/local.test.js` (the
  fixture's 168 answers, frame picks incl. Julia's round-half-even ties,
  JSON with Infinity/NaN, Float32 shortest decimals, mesh encoding, hazard
  edge cases). `viewer-common/tools/smoke.mjs --file RUN.h5 [--file …]`
  opens local files in the page in turn.

## [0.8.9] — 2026-11-25

Hazard to people, statistics over the run and classed palettes in the globe
viewers (branch `wpv-hazards`; PATCH: no solver change, results unchanged).
Step 3 of the viewer plan (`projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md`).

- `analysis/hazard.jl` (`FloodA5Analysis`; pure, also included by the
  viewers' server so both use one implementation): `depth_velocity` (D·V);
  `hazard_rating`, the Defra/EA HR = d (v + 0.5) + DF of FD2321/TR1 §3.5
  with the debris factor of its Table 3.1 by land use (`pasture`,
  `woodland`, `urban`), FD2320's precautionary 0.5 / 1 for any land use
  (`fd2320`, the default) or none; `arr_hazard_class`, the Australian
  combined hazard curves H1–H6 (AIDR Guideline 7-3 Table 2);
  `chen_critical_velocity` and `hazard_degree`, HD = min(1, v / U_c) with
  the incipient toppling velocity of Chen, Xia, Falconer & Guo (2019)
  eq. (10) (person mass and height, θ, α, β; checked against their Fig. 6
  and eq. (11) at θ = 90°); `derive_hazard` (elementwise, so a frame and a
  cell's series derive alike); `run_statistics` — per-cell maxima of depth,
  speed and each hazard, time of first flooding and duration of flooding
  (depth ≥ `flood_depth`, default 0.1 m), at frame resolution; `HazardParams`
  / `hazard_params`. On a run: `FA.hazard_field`, `FA.statistics_maps`.
  Sources: FD2321/TR1 Table 3.1 and 3.2, FD2321/TR2 Table 3.1, and the
  EA/HR Wallingford supplementary note of May 2008 (FD2320 thresholds
  0.75 / 1.25 / 2.0; FD2321 used 2.5 for the last) were read from the
  published PDFs; the ARR limits from AIDR Guideline 7-3 (2017). Reading of
  FD2321's "d > 0.75 m and/or v > 2": the third debris row applies above
  0.25 m when v > 2 m/s; shallower water keeps the first row.
- Viewer server: hazard variables served as frame variables
  (`/api/frames/{id}/hazard_rating?landuse=urban…`), in `/api/ranges` and
  `/api/cell`; statistics at `/api/stats/{name}` (computed once per parameter
  set off the server's thread; live runs: over the frames kept, recomputed
  as frames arrive); `/api/info` lists `derived` and `stats`
  (DATA_FORMATS §6). `/api/series` (0.8.8) gives hazard variables too,
  derived from the cell's depth and speed series with the same parameters.
- Viewers: the variable menu in groups (each frame / hazard / over the run);
  *Hazard parameters* (land use, person preset or own mass and height,
  Chen et al. calibration, flooded threshold), kept in the browser; a
  statistic disables the time slider and playback and draws the water at
  each cell's maximum depth. Classed palettes (`viewer-common/js/classes.js`)
  with a class legend: depth bands, D·V (ARR limits), HR (FD2320 classes;
  default for HR), ARR H1–H6 (always); editable breaks. The time-series
  plot (0.8.8) offers the hazard variables (not the statistics) and
  follows the hazard parameters.
- `io/variables.toml`: the new variables and statistics (definitions,
  classes, references, limitations — Chen et al. model toppling only, and
  sit above some urban low-depth/high-velocity data; a land-use raster is a
  later development), shown in *About* and in the static plots' captions;
  `FloodA5Plots` maps accept `:depth_velocity`, `:hazard_rating`,
  `:hazard_degree`, `:hazard_class_arr` (statistics: pass
  `FA.statistics_maps(r)[name]` as the values).

## [0.8.8] — 2026-11-25

Globe viewers: time series at cells (branch `wpv-series`; PATCH: no solver
change, results unchanged). Step 2 of the viewer plan
(`projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md`).

- Right-click on the map (a click: a right *drag* still turns or zooms the
  view) → *Add cell N to plot*. The first cell opens a plot panel docked at
  the bottom (⇆ moves it to the right; remembered); each further cell adds a
  line, up to eight, each in its own colour with a numbered marker on the
  map. The plot shows a variable over the whole run (its own selector:
  depth, level — where wet —, speed, wetted fraction), a dashed line at the
  frame shown; a click on the plot shows the nearest frame; ✕ in the legend
  removes a line (or the right-click menu on the cell); *CSV* downloads the
  lines (time in seconds, one column per cell). Live runs: the lines grow
  as frames arrive. Plots use Chart.js 4.5.1 (an ES module from jsDelivr,
  loaded when the first plot opens; `FLOODA5_CHARTJS_URL` for a local copy).
- Server: `GET /api/series/{i}?var=V` — one cell's values over the frames
  kept (`VisualisationServer.series_data`). A replay reads the cell's row
  from the file in one go (`set_replay!(…; series = (var, cell) -> …)`,
  `FO.read_series`); live runs read the frames in memory. `/api/config`
  gives `chartjs_url`. Engines gain `project(lon, lat, h)` (screen
  position of a point; the markers).

## [0.8.7] — 2026-11-25

Globe viewers: Run and About tabs, bed colours, *brown flood water*; run
and variable captions under the static plots (branch `wpv-run-about`;
PATCH: no solver change, default results unchanged). Step 1 of the viewer
plan (`projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md`, decisions
of 2026-11-25).

- The side panel has tabs: *Layers* (the controls), *Run* and *About*.
  *Run* shows what the run records — FloodA5 version and commit (and
  whether the tree was modified), flow model, backend, the command line
  (with a copy button), the mesh (family, DGGS level, cells, cell areas,
  file), frames, inputs with sizes and checksums, all options and the
  end-of-run summary — for a replay from the file, for a live run from the
  record its output file gets. *About* describes each variable (definition,
  unit, how the viewer derives or draws it) and the colour range, vertical
  offset, wet threshold and depth bias.
- One source for those descriptions: `io/variables.toml`
  (`FloodA5Output.variable_docs`); the run record for display:
  `FloodA5Output.run_record` (`run_info` without the long texts, plus mesh
  metadata, cell count and area range, frame times; JSON-safe). New
  endpoints `/api/run` and `/api/about`, `VisualisationServer.set_run!`
  (DATA_FORMATS §6).
- Static plots (`FloodA5Plots`): a caption under each figure — one line per
  run (label, flow model, mesh, FloodA5 version and commit) and the
  variable's definition, from the same sources (`run_caption`,
  `variable_caption`, `figure_caption`). Off for one figure with
  `caption = false` (export, publication), everywhere with
  `FP.CAPTIONS[] = false`.
- Bed surface: its own section with a palette menu (default *terrain*), a
  range (auto: 1st–99th percentiles of the bed, or fixed) and an opacity;
  below 1 the bed is drawn translucent (CesiumJS: the shader is rebuilt
  when that changes).
- *Sediment* is now **brown flood water** (no sediment variable yet): a
  narrow range around the colour of open brown floodwater in aerial
  photographs, ≈ rgb(112, 100, 92) at the middle. A palette chosen as
  *sediment* falls back to the default.
- `VERSION` had been left at 0.8.5 in 0.8.6; it is in step again.

## [0.8.6] — 2026-11-25

Globe viewers: fixes and small features from the first round of testing
(branch `wpv-fixes`; PATCH: no solver change). Larger requests are planned
in `projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md` (supplementary repository).

- Water no longer disappears when zoomed out over 3D tiles (both viewers):
  distant tiles are coarse and stand metres above the ground, hiding shallow
  water in large triangles. The water is drawn nearer the camera by up to
  1.5 % of its distance, only beyond 800 m (close views keep exact
  occlusion by buildings and trees).
- Vertical offset per basemap: `--offset` applies only to basemaps with
  ellipsoidal heights (Google tiles, Cesium World Terrain), no longer to
  Terrarium terrain or the flat maps (where it lifted the flood into the
  air). Without `--offset` the viewer estimates the offset from the basemap
  when it is first shown and saves it; a flat basemap puts the model's lowest
  ground on the map. The panel says where the value came from; its help text
  is no longer specific to Great Britain.
- Colour ranges fixed for the whole run by default: the 2nd–98th percentiles
  over wet cells and up to 24 frames (new endpoint `/api/ranges`,
  `VisualisationServer.run_ranges`; DATA_FORMATS §6); *frame* and *fixed*
  remain.
- Palettes: click the colour bar to choose one, remembered per variable;
  new *sediment* (light to dark brown), *cividis*, *greys*.
- *Cell borders*: white hairlines along the cells (off by default). The
  surfaces are now a fan from each cell's centre (one more vertex per cell).
- Timeline: a replay (or a finished run) opens at the first frame and plays
  only on *Play*; `web_replay.jl --autoplay` (or `?autoplay=1`) plays on
  opening. The slider follows the frame shown and no longer jumps back while
  it is dragged.
- CesiumJS: switching basemaps (e.g. OpenStreetMap → Google tiles) could
  fail in the render loop ("e._va is undefined"): basemap layers are now kept
  and shown or hidden, not destroyed and rebuilt.
- The opacity slider fits the panel.
- Tests: `/api/ranges` (live, replay, caching, sampling) in
  `test/test_vis_server.jl`; `--autoplay` and the ranges of a recorded run in
  `test/test_vis_live.jl`; the centre fan and the `_FA5EDGE` attribute in the
  viewer's JavaScript tests.

## [0.8.5] — 2026-11-25

Access keys from a key file (branch `maint-keys`; PATCH: no solver change).

- `io/FloodA5Keys.jl` (`flooda5_key`, `key_source`): a key such as
  `FLOODA5_CESIUM_ION_TOKEN` or `FLOODA5_GOOGLE_TILES_KEY` is taken from the
  environment variable of that name, then from JSON key files: the file
  named by `FLOODA5_KEYS`, `keys.json` in the repository root, and
  `~/.flooda5/keys.json` (one per user, shared by every clone and worktree).
  Template `keys.example.json`; `keys.json` is git-ignored everywhere and
  never served by the viewers' server. Same on every OS, and unlike a
  `setx` variable it does not depend on the environment a long-running
  VS Code was started with.
- The viewers' server uses it and logs where the Cesium ion token came from
  (never its value), or that there is none.
- `test/test_keys.jl` (CI): lookup order, placeholders, a malformed file,
  the template, the server's `/api/config`, `keys.json` not served and
  git-ignored.

## [0.8.4] — 2026-11-25

The globe viewer on a second engine, Giro3D (branch `wpv-giro3d`, on 0.8.3;
PATCH: no solver change, default results unchanged).

- **`--vis giro3d`**: the same viewer on Giro3D 2.0.7 / three.js r180, loaded
  as pinned ES modules from esm.sh through an import map (no build step).
  Google Photorealistic 3D Tiles through Cesium ion (or a Google key),
  Terrarium terrain + OSM (open AWS elevation tiles, no key), OSM, none.
  Comparison with the CesiumJS viewer, measurements and the problems found:
  docs/GLOBE_VIEWERS.md.
- The viewers share `visualisation/viewer-common/`: the app (controls, live
  and replay, playback, colour range, statistics, hover, offset estimate),
  the controls' markup, protocol, geometry, colour maps, and the surface
  model with the GLSL both engines use. `visualisation/cesium/` and
  `visualisation/giro3d/` are thin engine adapters. The JS tests, the smoke
  tool and `config.example.json` moved there; `config.json` is now
  `visualisation/viewer-common/config.json` (git-ignored).
- The server serves both viewers (`/viz/cesium/index.html`,
  `/viz/giro3d/index.html`) and only the viewer directories; `/` opens the
  one chosen (`start(; viewer)`, `/api/info` `viewer`).
- `visualisation/cesium_replay.jl` → `visualisation/web_replay.jl`
  (`--viewer cesium|giro3d`); `--replay RUN.h5 --vis giro3d`.
- The vertical offset is remembered per basemap (Google tiles are
  ellipsoidal, Terrarium is in the DEM's datum: 52.2 m vs ≈ 0 m at Carlisle).
- `geometry.ecefToGeodetic` (Giro3D 2.0.7's `Ellipsoid.toGeodetic` is ~130 m
  off in height); the viewers retry API requests on network errors;
  `smoke.mjs --profile N` and `--netlog`.

## [0.8.3] — 2026-11-24

The web globe viewer, rebuilt (branch `wpv-cesium`; PATCH: no solver change,
default results unchanged).

- **Globe viewer** (`--vis cesium`; visualisation/cesium/, CesiumJS 1.147,
  was 1.114): the water is one surface for the whole mesh, each cell flat at
  its water level (`datum + depth`: the cell elevation, or `sgs_z_min` for
  SGS), translucent and depth-tested — over Google Photorealistic 3D Tiles,
  buildings and trees stand out of the flood. A frame is one float-texture
  upload read by a custom shader (no per-cell entities, no geometry rebuild):
  the frame-update flash is gone, and a 582,012-cell mesh loads and changes
  frame in tens of milliseconds in the browser. Colour by depth, water
  level, speed or wetted fraction; per-frame or fixed range; depth ×
  exaggeration mode; the mesh's bed as a surface for the flat basemaps;
  hover read-out; time-based playback; LIVE follows a running simulation.
  Basemaps: Google Photorealistic 3D, Cesium World Terrain + OSM, OSM,
  none; no Bing imagery (it was loaded unseen before, against Cesium ion's
  terms when combined with Google tiles or OSM).
- **Vertical offset** from model heights (the DEM's datum) to ellipsoidal
  heights, with an *Estimate* from the basemap's heights: Carlisle (Google
  tiles) 52.2 m, 10th–50th percentile 52.0–52.6 m over 302 cells; remembered
  per domain in the browser.
- **Replay of output files in the globe viewer**:
  `visualisation/cesium_replay.jl RUN.h5` or `FloodModel.jl --replay RUN.h5
  --vis cesium` (frames read from the file on demand).
- **Server** (`visualisation/VisualisationServer.jl`) rewritten on plain
  HTTP.jl; **Oxygen.jl dropped** from the project. Protocol 2
  (docs/DATA_FORMATS.md §6): binary mesh, Float32 frames on demand, stable
  frame ids; a client that connects or reloads after the run learns that it
  is complete (before, it stayed "live"); the page's port follows
  `--vis-port` (it was fixed to 8080); requests must be addressed to
  localhost (DNS-rebinding guard; the ion token is served to the page);
  static paths are checked.
- ⚠ **Live frames thinned by wall-clock time** (default one per second,
  `FLOODA5_VIS_INTERVAL`) within a memory budget (1 GB,
  `FLOODA5_VIS_BUDGET_MB`; when full, every second frame is dropped). Before,
  every time step was stored in memory for the whole run (tens of GB on long
  runs) and, on the GPU backend, the device state was copied to the host
  every step. The end state is always sent. Diagnostics only: no result
  changes.
- With `--vis cesium` the time loop yields once per step: in Julia 1.12 the
  main task runs on interactive thread 1, the pool the server uses, and
  requests otherwise waited for the run (measured: 22 s per request on
  Carlisle; after: 0.09 s; solver speed unchanged, 44.4 steps/s).
- **A bare `--vis` now opens the GLMakie viewer** (it was the CesiumJS one).
- `FLOODA5_NO_VIS=1` stubs only GLMakie: the globe viewer works headless.
  `FLOODA5_VIS_LINGER=S` stops serving S s after a run (scripted runs).
- Cesium ion token from `FLOODA5_CESIUM_ION_TOKEN` (never in the
  repository); `visualisation/cesium/config.json` is now git-ignored.
- Tests: `test/test_vis_server.jl` (protocol, thinning, budget, late
  clients, shutdown, Host guard, paths), `test/test_vis_live.jl` (a live run
  followed over the WebSocket against its HDF5 output; replay; bare `--vis`),
  and the viewer's JavaScript (`cd visualisation/cesium && node --test test/*.test.js`,
  run by `ci/run_tests.sh`); a shared binary fixture keeps the Julia
  encoder and the JavaScript decoder in step. `ci/Project.toml` gains HTTP
  (1.x). Browser check (local): `visualisation/cesium/tools/smoke.mjs`.

## [0.8.2] — 2026-11-24

DGGS mesh families, Phase C (WP3 part 3, branch `wp3e-mesh-notebook`; PATCH:
default results unchanged, no solver change).

- Notebook `notebooks/06_mesh_families.ipynb`: the eleven meshes of the
  mesh-family cases — what they are and look like, their geometry as the
  solver sees it (`FA.mesh_quality`), isotropy, directional bias, orientation
  in a one-cell channel, Carlisle against the 5 m references, cost,
  advantages and disadvantages — and §8, the SGS checkerboard per mesh:
  provenance of the runs, per-mesh verdicts, the extrema through the event,
  residual maps, storage on the rising limb. `RERUN_SGS = true` runs the
  Carlisle SGS simulations again, as the check after SGS work.
- `mesh_families_carlisle` also runs Cartesian at 22.2 m, local-inertial
  and SGS (`carlisle_peak`'s is 20 m and local-inertial only); expected
  values re-recorded (the other runs' values unchanged).
- Finding: **SGS checkerboards on the Cartesian grid too** (41% of the wet
  cells local level extrema at the Eden peak, from about 46 h; the stored
  volume then runs away to 1.06 of the 5 m SWE benchmark). The open SGS
  finding of 0.8.1 is therefore not about non-A5 meshes: SGS is smooth only
  on A5, ISEA3H and jittered Voronoi (TEST_CASES §8; notebook 06 §8).
- `FA.wse_residual`: per cell, the water level minus its neighbours' mean and
  the local-extremum flag, which `FA.noise_metrics` aggregates (for maps).
- `FP.cellmap`/`cellmap!`: any per-cell values, or the cells alone (with
  outlines), on a subset of cells; `FP.cells_within(run, lon, lat, m)`
  picks a window; `FP.fieldmap(r; var = values)` takes a vector.
- Case runner: a run made in another checkout of the repository (a `git
  worktree`) is the same job — paths under the recording checkout (its `cwd`)
  are mapped onto this one before comparing commands, except those under the
  case root. Before, every run looked changed from another worktree and was
  run again.
- `notebooks/nbtools.jl execute`: `FLOODA5_NB_KERNEL` picks the kernel (e.g.
  one installed for a worktree).
- `FC` public names for the notebook: `MESH_FAMILIES`, `MESH_FAMILIES_DX`,
  `MESH_FAMILIES_CARLISLE`, `CHECKERBOARD_EXTREMA`, `CARLISLE_RISE`,
  `CARLISLE_PEAK_WINDOW`.
- Tests: `test/analysis/test_noise.jl` (`wse_residual` = `noise_metrics` per
  cell), `test/test_plots.jl` (vector values, `cells_within`),
  `test/test_cases.jl` C3 (worktree paths), `test/test_plots_render.jl`
  (`cellmap`).

## [0.8.1] — 2026-11-24

DGGS mesh families, Phase B (WP3 part 3, branch `wp3d-dggs-phaseb`; PATCH:
default results unchanged, no solver change).

- `--meshgen-dggs AOI --dggs FAMILY (--meshres N | --dx M)`: a mesh of H3,
  IGEO7, ISEA3H, HEALPix or rHEALPix in any run (`mesh/DGGSGrid.jl` over
  `mesh/dggs_families.py`); `--dx` picks the level whose mean cell area is
  nearest M². `--dggs a5` is the A5 generator (`--dx` → the nearest A5
  resolution). Both flags are mesh-generation flags for `--rerun`
  (`FO._DROP_FLAGS`).
- The mesh file's `flooda5` metadata (DGGS family, library, level, refinement
  ratio, equal-area flag, target size; the planar family) is kept on the mesh
  (`A5Mesh.meta`), through re-saves (DEM sampling, SGS tables — before, a
  re-save dropped it), and written to the output's `/mesh` attributes;
  `FO.mesh_info(path)` reads them.
- A mesh file that records no cell-ID scheme is still taken for an old A5 or
  Cartesian mesh, now with a warning. `build_sgs_tables!` refuses a non-A5
  mesh without neighbour lists (it used to call pya5's `grid_disk` with its
  IDs).
- `--meshload MESH --dem DEM --flow-model sgs` builds the SGS tables of a
  loaded mesh that has none (it used to refuse it): the command-line route
  to SGS on the planar families and imported meshes (with `--meshout …
  --mesh-only` to save them). Without `--dem` such a mesh is still refused.
- `PolyMesh`/`PolyCell`: neutral aliases of `A5Mesh`/`A5Cell`.
- `mesh/mesh_families.py --clip`: keep the cells that meet the AOI polygon
  (the A5/DGGS rule) instead of covering its bounding box.
  `test/channel_bench/make_channel_meshes.jl --families`: the channel series
  on other mesh families.
- Cases (`cases/registry/mesh_families.jl`, docs/TEST_CASES.md §8): point
  spread, planar slope, straight channel (0–45°) and Carlisle to the Eden peak
  (against the 5 m references) on eleven meshes at ~22 m — Cartesian, rotated
  Cartesian, hexagons, centroidal and jittered Voronoi, A5, H3, IGEO7,
  ISEA3H, HEALPix, rHEALPix.
- Finding (open): SGS checkerboards at Carlisle on most non-A5 meshes —
  rHEALPix and the rotated grid worst, with the storage oscillating;
  smooth on A5, ISEA3H and jittered Voronoi; local-inertial is smooth on
  every mesh (`mesh_families_carlisle`, an `xfail` check; TEST_CASES §8).
- Tests: `test/test_dggs_meshes.jl` C (command line, metadata, `--rerun`),
  `test/test_sgs_general.jl` S8 (SGS on H3, IGEO7, HEALPix), S9 (SGS tables
  for a loaded mesh from the command line).
- The compiled `a5_bridge` `.pyc` files are no longer tracked.

## [0.8.0] — 2026-11-23

DEM registration and the SGS critical-flow limit (MINOR, result-changing;
branch `wps3-dem-critical`; HYDRAULICS §15.14). **Breaking for SGS meshes**:
tables are layout 4; SGS meshes must be rebuilt. Meshes sampled from a DEM
before 0.8.0 load with a warning and should be regenerated. Every case whose
results moved is re-recorded (deliberately); the 0.7.1 values are kept in
`cases/expected/v0.7.1/` (and tag `v0.7.1`).

### Fixed
- ⚠ **DEM sampled at pixel centres.** GDAL rasters are pixel-is-area (the
  geotransform is the top-left pixel's corner; a value belongs at the pixel's
  centre). `_bilinear` used `_affine_to_pixel`'s corner-based coordinates as
  centre-based ones, so every DEM sample — cell means and centroids, SGS
  tables and face profiles, in every engine — was taken half a pixel east and
  south of its point. It now interpolates between pixel centres, holds the
  edge value within half a pixel of the raster's edge, and next to nodata
  takes the pixel containing the point. Meshes record `dem_registration`
  (`docs/DATA_FORMATS.md`); `--meshload` warns when it is missing.
  Raster-aligned Cartesian meshes (`--meshgen-cartesian dem.tif`) read the
  pixels directly and are unchanged. Test: `test/test_dem_registration.jl`.
- `cases/scripts/sgs_strip_dem.jl` writes the analytic terrain at pixel
  centres (it wrote it at the corners, matching the old sampling).

### Added
- ⚠ **SGS critical-flow limit** (`--sgs-critical on|off`, default on). A face
  passes at most √g·(2/3)^(3/2)·C(E) with C = Σ dx·(E − a_p)^(3/2), its
  critical integral at the upstream cell's level (each strip of the face at
  its own critical depth; for a rectangle the broad-crested-weir discharge
  1.705·W·H^1.5), times the flow-path factor; SGS open-boundary outflow is
  capped the same way. The local-inertial equations have no advection, so
  nothing limited the discharge where water falls away from a face: 0.7.x
  passed ≈ 9.8 m³/s through the 4 m notch of `planar_embankment` at 1.07 m of
  head (critical ≈ 7.5). C is tabulated per face (`sgs_face_C`, SGS layout 4,
  `A5Grid.face_crit`). The limit is an option of the hydraulics, not of the
  SGS geometry; `off` gives the 0.7.x flux. CPU, device and the frozen
  `reference_steps.jl` loop (a deliberate change). Tests: `test/test_sgs_general.jl`
  S7, `test/test_sgs_edge_ra_tables.jl` F9.
- Case `carlisle_peak`: job `a5_sgs_nocrit` (SGS with `--sgs-critical off`).
- Case `planar_embankment`: check that SGS passes no more than critical flow
  through the notch (metrics `q_down_5to6h`, `head_5to6h`).

### Changed
- **The 5 m SWE run is the Carlisle benchmark** (decided 2026-11-23).
  `carlisle_peak_fine` compares every 22 m run of `carlisle_peak` with both
  5 m runs (`vs5:<job>`) and holds the agreed SGS bar; `carlisle_peak` no
  longer checks SGS against the 22 m SWE run.
- `carlisle` 6 h: the A5 records are re-recorded (SWE 2.641 Mm³, arrival
  6 676 s; SGS 1.775 Mm³).

### Results
- Cartesian runs on raster-aligned meshes (Carlisle 20 m and 5 m) are
  bit-identical. A5 Carlisle: SWE 6 h storage −1.8%, at the peak +1.7%;
  local-inertial −0.1%. Synthetic cases move by the half-pixel shift
  (`planar_wall` thin ridge 5.36 → 5.39 m on A5); `sgs_strip` levels are
  unchanged to the millimetre, but SGS holds 1–2% more water over the whole
  strip: the critical-flow cap binds on its steep (2.1%, supercritical)
  outlet reach, since it neglects the velocity head.
- `planar_embankment`: SGS passes 8.24 m³/s at 1.38 m of head through the
  notch (weir 11.07); the pond stands 0.33 m higher and 26% less water is
  downstream at 6 h.
- Carlisle at the peak, against 5 m SWE: SGS stores 0.74 (0.71 without the
  limit; 22 m local-inertial 0.94, 22 m SWE 1.26); against 5 m local-inertial
  SGS is 0.87 with the lowest rms level error of any 22 m model (0.16 m;
  local-inertial 0.34 m). **The SGS bar fails** and is an expected failure:
  at 5 m the local-inertial equations themselves store 0.84 of SWE, so the
  gap is the hydraulics, not the SGS geometry.

### Known issues (found, not changed)
- The SGS gap to 5 m SWE needs a hydraulic formulation with momentum on the
  SGS geometry (PROJECT_STATE, "Next").
- A local-inertial run can overrun its end time by up to 0.1 s: the 0.1 s
  floor on Δt (`FloodModel.jl`, run loop) is applied after the clamp to the
  remaining time (`planar_wall` `a5_sgs_wide` ends at 21 600.086 s).
- Output files and expected values written by these reruns record version
  0.7.1 (dirty tree): `FLOODA5_VERSION` was bumped after them.

## [0.7.1] — 2026-11-23

DGGS mesh families (WP3 part 3, branch `wp3c-dggs-meshes`; PATCH: default
results unchanged, no solver change).

- `mesh/dggs_families.py`: FloodA5 GeoParquet meshes of five discrete
  global grids, loaded with `--meshload`: H3 (`h3`), and through DGGAL IGEO7
  (ISEA7H, Z7 index), ISEA3H, HEALPix and rHEALPix. Level by `--res` or
  `--target-dx`, `--list-levels`. Cells intersecting the AOI; corner vertices
  snapped to one copy per vertex; shared-side adjacency checked against the
  library's neighbour lists.
- `docs/MESHES.md`: the mesh contract (what any mesh must provide), the
  families, levels near flood-model cell sizes, and each family's geometry at
  Carlisle. None of the DGGS is orthogonal at 55° N in FloodA5's metric
  (H3 4.7°, rHEALPix 7°, IGEO7 10–16°, ISEA3H 13.5°, HEALPix 33°; A5 20–37°).
- `FA.mesh_quality` (`analysis/mesh_quality.jl`): valence, area uniformity,
  non-orthogonality, skew, centre offset, face-orientation moments, vertex
  valence of any mesh.
- Tests: `test/test_dggs_meshes.jl` (in CI: generator, geometry, lake at rest
  under four options, mass conservation, HEALPix against healpix-geo),
  `test/analysis/test_mesh_quality.jl`.
- Python packages: `h3`, `dggal==0.0.6`, `ecrt==0.0.6`, `healpix-geo`
  (`setup.jl`, `ci/python-requirements.txt`). dggal 0.0.7's Windows wheel
  crashes on setup (two copies of its eC runtime in one process).
- Replay: only A5 meshes get an A5 level in the diagnostics panel.
- Plan and findings: `docs/notes/MeshFamilies_DGGS_Plan_2026-11.md`.

## [0.7.0] — 2026-11-18

SGS conveyance (MINOR; plan: `projectinstructions/FloodA5_SGS_Conveyance_Fix_Plan_2026-11.md`).
**Breaking for SGS meshes**: tables are layout 3; SGS meshes must be rebuilt
(`--meshgen … --flow-model sgs --mesh-only`). Standard flow and SWE are
unchanged. The results of the old scheme are kept in `cases/expected/v0.6.0/`
(and tag `v0.6.0`).

### Changed
- ⚠ **SGS face tables (layout 3).** Each face has one table, built from the
  DEM profile along the face on the face's own levels and held identically by
  both cells: flow area A, wetted width P and the conveyance integral
  K = Σ dx·h^(5/3) at 64 levels (denser near the sill) that reach above the
  highest point of the profile, continued analytically above them with the
  full face width (`A5Grid.face_table_from_profile`, `face_apk`; mesh columns
  `sgs_face_A/P/K`, `sgs_face_zr`, with `sgs_edge_sills` their lowest level).
  They replace the per-cell edge curves (`sgs_edge_area_curve`,
  `sgs_edge_perim_curve`), which were tabulated on each cell's own levels and
  stopped growing at the cell's highest DEM sample, and the ½(A_i + A_j),
  ½(R_i + R_j) average in the kernel (case `sgs_strip`: faces carried 39% and
  2.5% of their flow area 0.4 m above the floodplain). `edge_area_radius`,
  `flow_area_from_wse`, `wetted_perim_from_wse`, `hydraulic_radius_from_wse`
  are removed. Layout ≤ 2 meshes are refused (`FLOODA5_ALLOW_OLD_SGS_LAYOUT`
  is gone); `--rerun` of an SGS output written by ≤ 0.6.0 stops with that
  message. The SGS loop of the frozen `surfacewater/reference_steps.jl` reads
  the face tables too (a deliberate change of that file).
- ⚠ **SGS conveyance form**: `--sgs-conveyance 2d` (default) uses
  R = (K/A)^(3/2), i.e. A·R^(2/3) = Σ h^(5/3) over the face profile — the
  uniform-flow limit of the 2-D equations. `lumped` keeps one hydraulic radius
  A/P per face (for the record: overbank it gives the level of a cell-mean
  grid, +0.04 m on the strip).
- ⚠ **SGS limiter**: `--scheme recommended` gives SGS `--limiter positivity`
  and `--boundary-limiter positivity` (were `legacy`). The legacy donor cap,
  V/(2N) per face per step, held in-bank channels brim-full.
- ⚠ **SGS time step**: the depth of an SGS cell for the time-step limit is
  the hydraulic depth V/A_wet (was: capped at the cell's relief
  z_max − z_min, which ignored all water above a cell's highest ground), and
  a second limit follows from the faces: rate ½√(g·Σ(A_f/L_f)/A_wet)
  (`_sgs_face_rate`), which for a fully wet square cell is the usual
  √(g h)/Δx and is larger for a cell that holds only a sliver of a channel
  its faces convey. Without it such cells oscillated and held the channel up
  by up to 0.5 m on the A5 strip. The wetted area is floored at 2% of the
  cell (`SGS_CFL_AREA_FLOOR`). The A5 strip in bank takes 24 300 steps (0.6.0: 5 650).

- ⚠ **SGS tables from a working grid.** `build_sgs_tables!` covers each cell
  with a regular grid in a local metric frame (`--sgs-samples`, default 1600
  points; was 512 Halton candidates whatever the cell or DEM) and takes the
  storage curve from four samples per grid cell. The cell's lowest and
  highest ground are taken over the grid and the profiles of its sides; faces
  are sampled at midpoints, half a pixel apart in metres (a geographic DEM's
  pixel was converted without cos φ). Volume tables on the strip: within 0.3%
  of the terrain on Cartesian (was 0.3%), 3% at worst on A5 (was 13%).
- ⚠ **SGS connectivity** (`--sgs-connectivity on|off` at mesh generation,
  default on). A point of a face conveys only from the level at which it is
  connected, below the water, to the main basins of both cells
  (`_cell_connectivity`: the first region of ground below a rising level to
  reach two sides of the cell). A barrier inside a cell closes the faces
  behind it up to its crest, or to the sill of a gap through it. Until now a
  wall narrower than a cell passed water as soon as the water reached it
  (case `planar_wall`: 58–64% of the inflow, 2.9 m below the crest); now
  nothing passes. `sgs_face_raise` in the mesh records the change per face.
- ⚠ **SGS flow-path geometry** (`--sgs-path on|off`, default on). For cells
  less than half wet the slope across a face is taken over the distance
  between the two cells' wetted-area centroids and the face's flow area is
  projected on that line (`_sgs_path`, `wet_centroid`; mesh columns
  `sgs_wet_cx/cy`); fully wet cells are unchanged. A channel crossing A5
  cells in a zig-zag was 0.054 m high in bank on the strip; now 0.010 m low.
- ⚠ **SGS open boundaries.** Polygon sides on the mesh boundary have face
  tables like any face (was: a rectangle of the side's width above the
  cell's lowest ground; still the fallback where the DEM does not reach).
  The bed slope of the normal-depth outflow is the slope of the local valley
  floor (`_sgs_boundary_slopes`), at least 1e-4 (was: the difference of two
  cells' mean elevations, 21% at the Eden's outlet from the Carlisle domain).
- ⚠ **SGS Courant number**: 0.5 by default with `--flow-model sgs` (0.7
  otherwise). At 0.7 the SGS water surface in deep, fast channel flow breaks
  into cell-scale oscillations (Carlisle at the peak); 0.5 and 0.35 agree.
- `mesh/repair_sgs_slot_order.py` is removed (it repaired tables to layout 2,
  which is now refused). `visualisation/visualise_mesh.jl` reads the face
  tables.
- Case runner: a run's identity no longer depends on the case of a Windows
  path or on its slash direction.

### Results
- `sgs_strip` (8 h; analytic 2-D level δ −0.267 / +0.372 m). Cartesian 20 m:
  in bank −0.279 m, storage 3.50 m²/m (3.57); overbank +0.362 m, 26.7
  (27.3). A5 res 18: in bank −0.277 m, 3.50; overbank +0.340 m (2-D for its
  own width +0.290, local-inertial on the same mesh +0.390). All steady.
  0.6.0: −0.008 m (brim-full), +4.75 m and still filling, out of bank on A5.
  All eight expected failures of the case now pass and are ordinary checks.
- `planar_wall`: nothing passes the 4 m wall (was 58–64% of the inflow); the
  pond reaches 7.07–7.14 m (level pool 7.12, crest 8.05). Its two expected
  failures are ordinary checks again.
- `planar_embankment`: water passes at the notch (first at 842 s, upstream
  level 5.25 m; sill 5.0, crest 7.0) and nowhere else.
- `carlisle` (6 h, rain): SGS stores 1.755 Mm³ (0.6.0: 4.30; local-inertial
  3.37, SWE 2.69). The regression pin is re-recorded.
- `carlisle_peak` (55.5 h): SGS stores 5.94 Mm³ (0.6.0: 26.7), 0.57 of the
  SWE reference where local-inertial stores 0.76; its water surface is a
  median 0.89 m below SWE (local-inertial 0.53 m). SGS is smooth and
  time-step converged but conveys more than local-inertial. **The check
  "SGS no further from SWE than local-inertial" fails** and is kept as an
  expected failure. Against 5 m references (case `carlisle_peak_fine`,
  plan §9.3a, HYDRAULICS §15.13) the 22 m SWE run stores 24% more than 5 m
  SWE, so it is not a reference; SGS is the 22 m model closest to the 5 m
  local-inertial solution (level −0.14 m median, rms 0.24 m), but
  local-inertial conveys more than SWE at 5 m, and against 5 m SWE SGS
  stores 0.70 (22 m local-inertial 0.94, by compensating errors).
- Checks changed where they referred to the old scheme: the strip's default
  run is compared with the 2-D analytic (was the lumped one), tolerances
  0.03 m / 5% overbank (were 0.10 m / 10%); its `--limiter positivity` jobs
  are replaced by `--sgs-conveyance lumped` jobs; its A5 overbank check
  compares SGS and local-inertial with the 2-D answer for the strip's own
  width. Expected values of `sgs_strip`, `planar_wall`, `planar_embankment`,
  `carlisle` and `carlisle_peak` are re-recorded.
- Cost: the A5 strip in bank takes 24 300 steps (5 650 in 0.6.0); Carlisle
  peak 181 000 (131 000).

### Known issues (found, not changed)
- DEM sampling is half a pixel off for pixel-is-area rasters, in every
  engine (`_affine_to_pixel` / `_bilinear`); see plan §9.4.
- Standard flow's open-boundary slope is taken from two cells' mean
  elevations, as SGS's was.
- The η filter is not implemented for SGS (hence its Courant 0.5).
- SGS faces have no critical-flow limit: `planar_embankment` passes ≈ 9.8 m³/s
  through the 4 m notch at ≈ 1.07 m of head, above critical flow (≈ 7.5 m³/s
  for a 4 m rectangle).

### Added
- `cases/expected/v0.6.0/`: the SGS metrics, table diagnostics and strip
  series of the scheme up to 0.6.0, for before/after comparison.
- **Case `carlisle_peak_fine`**: the `carlisle_peak` event on the 5 m
  Cartesian mesh (582 012 cells), local-inertial and SWE on the GPU — the
  fine reference for the 22 m runs (docs/TEST_CASES.md §3).
- `notebooks/04_sgs_strip.ipynb` rewritten as diagnosis, fix and result
  against the 0.6.0 record.
- `notebooks/01_planar.ipynb`: the `planar_wall` case (new section) and the
  embankment against the 0.6.0 record; `notebooks/03_carlisle.ipynb`: SGS
  compared like for like (storage, mean depth, level where both are wet,
  wetted area from its tables) and an optional section against the 5 m
  references.
- **Analysis: fine-to-coarse comparison** (`analysis/aggregate.jl`):
  `FA.assign_to_cells`, `aggregate_fine` / `aggregate_onto`, `compare_to_fine`
  (a fine run aggregated onto a coarse mesh's cells: volume, wetted fraction,
  level of the wet fine cells), and SGS tables read from the mesh GeoParquet
  without Python (`load_sgs_tables`, `sgs_tables`, `wetted_fraction`,
  `load_mesh_columns`). Tests in `test/analysis/test_aggregate.jl`.
- ⚠ `carlisle` / `carlisle_peak`: the SGS flooded area (`wet_area_0p1*`) is
  the area wetter than 0.1 m from its tables (was whole cells by SGS depth,
  measured from the lowest ground: 4.96 km² at 6 h, now 1.79). Re-recorded.

SGS conveyance fix, Phase 0 (measurement; no model change).

### Added
- **Registry case `planar_wall`** (`cases/registry/planar_wall.jl`,
  docs/TEST_CASES.md §2d): a 3 m wall with no notch across the planar slope —
  40 m wide, and 4 m wide (narrower than a cell, placed inside the cells of
  both meshes) — with an inflow that provably leaves the level-pool level
  0.84–0.93 m below the crest (wedge storage). SGS and standard flow on A5
  res 18 and Cartesian 20 m (8 runs). Checks state the correct behaviour (no
  water downstream of the wall). Result (0.6.0): the 40 m wall holds in every
  run; **SGS passes the 4 m wall as soon as the water reaches it** (A5: 722 s,
  124 660 m³ = 58.5% of the inflow downstream at 6 h; Cartesian: 962 s,
  136 972 m³ = 64.3%), 2.9 m below the crest, wherever the water touches the
  wall. Those two checks are expected failures (`xfail`: SGS conveys through
  faces from the DEM profile on the face line only; a barrier inside a cell
  is not seen) until SGS connectivity (fix plan, Phase 3). Standard flow on
  the thin wall is checked against the cell-mean ridge it sees (5.36 /
  5.63 m), which it overtops early. Expected values:
  `cases/expected/planar_wall.toml` (the six passing runs and the design).
- `planar_embankment`: metrics (no new checks) of where water crosses the
  embankment — the first-row cells downstream that hold water at the first
  crossing and through time, by distance from the notch, cells on the
  embankment's top, and downstream volume by band. Result: nothing crosses
  the 40 m embankment away from the notch in either run. The 61 existing
  expected values are unchanged; the new keys were added with `baseline`.
- `test/test_cases.jl` C5: `planar_wall` is registered, well formed, its
  design holds and its expected-failure marks are in place (no simulation).

### Fixed
- `generate_planar_embankment_dem.py --notch-width 0` now guarantees no
  notch (the `|y − y_c| ≤ w/2` test would open a pixel row whose centre lies
  on the centre line, i.e. with an odd number of rows). Output for existing
  arguments is unchanged.

### Known issues (found, not changed)
- DEM sampling is registered half a pixel off: `_affine_to_pixel`
  (`mesh/A5Grid.jl`) returns corner-based pixel coordinates and `_bilinear`
  treats them as centre-based, so every DEM sample (cell means, SGS tables,
  face profiles) is taken half a pixel east and south of its position
  (docs/TEST_CASES.md §2d, "A wall on a face").


## [0.6.0] — 2026-11-15

Float32 on the device backends and a benchmark suite (MINOR). Float64
results do not change: end states of every solver and scheme (CPU and
`ka-cpu`, A5 and Cartesian) and the determinism fingerprints are bit-identical
to 0.5.1.

### Added
- **`--precision float64|float32`** (default `float64`; `run_flood_model(…;
  precision = :float32)`), for `--backend gpu` and `ka-cpu` with
  `--flow-model standard` and `swe`. Refused with `--backend cpu` and with
  `--flow-model sgs` (SGS follows its review). Depths, fluxes and velocities
  are Float32; cell volumes and the mass totals stay Float64, and water
  levels are formed relative to a local datum from bed differences taken in
  Float64, so the result does not depend on the elevation datum
  (`surfacewater/precision.jl`: `@fp`, `AbsLevels` / `RelLevels`). Output
  files stay Float64; the precision is recorded in the run record and
  `--rerun` reproduces it. Design: review plan §7.3; ARCHITECTURE.md §6a.
- `test/test_precision.jl` (CI, on `ka-cpu`): Float32 against Float64 for
  SWE and standard flow on Cartesian and A5 meshes, mass conservation, lake
  at rest, datum independence, relative levels in Float64, the refusals and
  the run record.
- **Benchmark suite** `benchmarks/` (`julia --project=. benchmarks/run.jl`,
  `benchmarks/README.md`): step throughput on 256²–2048² grids, CPU thread
  scaling, and whole Carlisle runs (A5 res 18, Cartesian 20 m and 5 m;
  standard and SWE) on CPU, GPU Float64 and GPU Float32, each device run
  with its agreement with the reference (storage, extent CSI, depth
  differences) and its conservation error. Results are kept per machine in
  `benchmarks/results/<host>-<os>.json` (committed); first file:
  `c002it-windows.json`.
- **`notebooks/05_performance.ipynb`**: runs the sections the machine it is
  opened on has no results for, and compares all machines with results.
- `test/benchmark_backend.jl`: `--configs cpu,gpu64,gpu32[,kacpu64,kacpu32]`,
  `--json FILE`, accuracy and mass columns, and `--min-cellsteps N` (the
  suite times at least 10⁸ cell-steps per grid: 20 steps of a 512² grid take
  a GPU 0.1 s, and gave rates wrong by up to a factor of three).
- Run summary: `wall_first5_steps_s` (set-up and the first five steps; on a
  GPU this is kernel compilation, about 40–50 s per process).

### Changed
- Kernel bodies (`kernels.jl`, `std_device.jl`, `swe_reference.jl`,
  `eta_filter.jl`, the flux functions in `flow2d.jl`) are generic in the
  float type. `KADevice` carries the precision and the level mode.
- `test/benchmark_sim.jl` and `test/run_benchmarks.*` are superseded by the
  suite's `threads` section (they need a configuration file that is no
  longer in the repository).

### Measured (Xeon E-2278G, 16 threads; RTX 4000 SFF Ada; `benchmarks/run.jl --force`)
- Step throughput at 2048² (million cell-steps per second), CPU / GPU
  Float64 / GPU Float32: standard 8.3 / 113 / 214; SWE 3.4 / 36 / 83.
  Float32 differs from the CPU by 1·10⁻⁶ m there (26 steps) and by at most
  8·10⁻⁶ m on any grid (256², 1 528 steps). Below
  about 3·10⁵ cells Float32 gains little (standard, 512²: 64 / 60).
- Carlisle 5 m (582 012 cells, 1 h), SWE: GPU Float64 893 s, Float32 388 s;
  standard: CPU 328 s, GPU Float64 103 s, Float32 94 s (about 50 s of each
  GPU run is compilation). At 30 000 cells the GPU is slower than the CPU
  for standard flow in either precision.
- Float32 against Float64 at the end of the Carlisle runs: storage within
  1.4·10⁻⁴, extent CSI ≥ 0.997, conservation error ≤ 1·10⁻⁸. On A5 with the
  standard solver, the differences are those of GPU Float64 against the CPU
  (round-off amplification, HYDRAULICS.md §15.5), not of the precision.

## [0.5.1] — 2026-11-12

SGS diagnosis (PATCH): the model does not change.

### Added
- Registry checks can be marked as **expected failures**: `Check(…; xfail =
  "reason")` marks a check of the correct behaviour that the model is known to
  fail. It is reported as XFAIL (or XPASS once it passes) and does not fail
  the case. `sgs_strip` uses it for its eight SGS checks (`test/test_cases.jl`
  C4).
- **Case `sgs_strip`** (`cases/registry/sgs_strip.jl`, docs/TEST_CASES.md
  §2c) and **`notebooks/04_sgs_strip.ipynb`**: the SGS model on a synthetic
  4 km strip (0.1% slope, a 5 × 1 m channel a quarter of a 20 m cell wide) at
  3 m³/s (in bank) and 20 m³/s (overbank), against analytic normal levels,
  storage per metre and kinematic arrival times, and a 2.5 m reference.
  SGS, local-inertial and SWE on Cartesian 20 m, SGS and local-inertial on
  A5 res 18. The SGS tables are also checked against the terrain directly
  (`cases/scripts/sgs_strip_tables.jl`). DEM generator
  `cases/scripts/sgs_strip_dem.jl`; terrain and analytic hydraulics in
  `cases/scripts/sgs_strip_geometry.jl`; public names `FC.sgs_strip_*`.
  The model does not change.
- Findings (no fix yet; the case's SGS checks state the correct behaviour
  and fail): SGS storage is measured correctly (the tables' volumes match
  the terrain to 0.3%; overbank, the volume above the cell-mean bed equals
  the water volume), but SGS **conveyance** is wrong:
  - overbank, an edge's flow area stops growing once the water is above the
    cell's highest DEM sample (`build_sgs_tables!` tabulates edge curves on
    the cell's own z_min–z_max levels; `edge_area_radius` returns the top
    value above them), and face points above a cell's z_max never wet in its
    table. The 20 m strip's tables carry at most 4.9 m³/s at the bed slope,
    so at 20 m³/s SGS piles water up towards the inflow (δ 4.8 m above the
    floodplain against 0.37 m; still filling at 8 h). At the `carlisle_peak`
    SGS peak 85% of wet cells, holding 26.1 of 26.7 Mm³, are in this regime;
  - in bank, the legacy donor cap (SGS's default limiter) lets a cell pass
    at most V/(2N) per step, so a channel cell carries Q only once it holds
    2N·Q·Δt: at 3 m³/s the channel fills to the brim (storage 1.4× the true
    value). `--limiter positivity` gives the analytic level within 1 cm.

## [0.5.0] — 2026-11-12

Default engine change (MINOR): results change for runs that did not pass
`--flow-model`.

### Changed
- ⚠ **`--flow-model` defaults to `standard`** (local-inertial on cell-mean
  elevations with the `--scheme recommended` corrections); it was `sgs`. SGS
  becomes opt-in until its behaviour at high stage is understood: on
  Carlisle near the Eden peak (`carlisle_peak`, 0.4.4) it stores 2.5× the
  SWE reference, with the river about 3 m high (PROJECT_STATE). A synthetic
  test is in preparation. Pass `--flow-model sgs` to keep the old default.
  Every test and registry case already passes `--flow-model` explicitly.
- `--rerun` pins the recorded flow model and scheme when the original
  command left them to the default. A run made before 0.5.0 without
  `--flow-model` therefore reruns as SGS, as recorded, not as today's
  default (`test/test_output_format.jl` O7).
- **`--flow-model swe` is a full option**, not a benchmark only. It keeps
  inertia and advection, so it is the engine for dam breaks, rapidly varying
  or supercritical flow and flow around obstacles, and remains the reference
  for the local-inertial engines. It typically takes 5–10× the time steps of
  `standard`. Help text, README and USER_GUIDE updated.
- The `--flow-model` error message lists all three engines.

## [0.4.4] — 2026-11-11

Notebooks (PATCH). Neither the model nor the output format changes.

### Added
- **Jupyter notebooks** (`notebooks/`, Julia kernel "Julia FloodA5 1.12";
  setup in `notebooks/README.md`). Each describes every simulation, prints
  the exact commands from the case registry and documents the Python steps
  (DEM generator; A5 meshes through pya5). It runs the cases through the
  registry (reusing completed runs), checks them, and analyses and draws
  them with `FA`/`FP`:
  - `01_planar.ipynb`: `planar_slope`, `planar_embankment`
  - `02_channel.ipynb`: `channel`; the full series from its recorded
    metrics
  - `03_carlisle.ipynb`: `carlisle_peak` (the 6 h before the Eden peak)
    and `carlisle` (the recorded first-6 h window), each with where the
    local-inertial storage difference sits by depth class
- Committed without outputs. `notebooks/nbtools.jl` provides `strip` (writes
  JSON exactly as Jupyter does), `check` and `execute` (headless).
  `test/test_notebooks.jl` (in CI) checks the kernel, that there are no
  outputs, and that every code cell parses.
- Public registry names for notebooks and scripts (`FC.python_command`,
  `FC.planar_params`, `FC.log_mass_balance`, …), checked in
  `test/test_cases.jl`.
- **Case `carlisle_peak`**: the event from t = 0 to the River Eden's peak
  (1273 m³/s at 55.5 h), river inflow only. It compares the 6 h before the
  peak, local-inertial vs SWE on A5 res 18 and Cartesian 20 m, plus SGS.
  - At high flow the local-inertial default stores **less** than SWE: 0.76×
    on A5 and 0.85× on Cartesian at 55.5 h, with CSI 0.87 / 0.92. That is the
    opposite of the early window's 1.25× / 0.99×.
  - SGS stores 2.5× SWE, with the Eden about 3 m high. This is flagged for
    investigation and has no check.

### Changed
- Case `carlisle`: its 10 mm/h uniform rain is labelled as a test forcing (60
  mm in 6 h, not the event's rain). Two no-rain local-inertial runs show
  that the rain-on-grid sheet makes most of the 6 h flooded area: 3.06 vs
  1.25 km² on A5, 2.32 vs 0.95 km² on Cartesian. There is a new
  `wet_area_0p1_km2` metric and a check.

## [0.4.3] — 2026-11-11

Visualisation (PATCH). Simulation results and the output format do not change.

### Added
- **Static figures** (`visualisation/FloodA5Plots.jl`, module `FloodA5Plots`, `FP`;
  CairoMakie, from output files alone — no model, no Python): `fieldmap`/`fieldmaps`
  (maps at a frame or a runs × times grid, elevation basemap, sources from
  the run record), `difference_map`, `extent_map` (hits/misses/false alarms,
  CSI), `timeseries`, `storage_plot`, `mass_balance_plot`; `…!` forms compose
  panels for notebooks. A Carlisle res-18 map takes ~0.9 s. Makie-free parts
  (`visualisation/plots/core.jl`: colour ranges, frame matching, extent
  classes, projections) are tested in CI (`test/test_plots.jl`); every plot
  type is rendered by the local `test/test_plots_render.jl`.
  `test/channel_bench/plot_channel.jl` uses its house style (byte-identical
  PNGs).
- **Replay of an output file in the GLMakie viewer**: `visualisation/replay.jl
  RUN.h5` (or `FloodModel.jl --replay RUN.h5`) opens the `--vis makie` window
  on a recorded run, with a frame slider, a playback-speed slider and
  Pause/Play; `--fps`, `--var`, `--from`, `--paused`, `--png`, `--mesh`
  (format-1 files). The time-series strip shows the whole run: domain volume
  and wet cells from the file, input volume rebuilt from the run record and
  checked against `/summary`, removed volume per frame for closed boundaries
  only (open: added − domain plus the recorded end values). Data layer
  `visualisation/replay_data.jl` (`FloodA5Replay`, Makie-free, tested by
  `test/test_replay.jl` in CI); `MakieVisualiser.replay`, `replay_goto!`;
  `MakieVisualiser.start` gains `show_window`. The live path is unchanged
  (`push_frame!` now shares `_set_fields!`/`_update_ring_bars!` with replay).

### Changed
- `visualisation/FloodViewer.jl` reads the mesh in pure Julia (`FA`) — from
  the output file for format 2, so the parquet argument is optional — instead
  of through Python; its emoji button labels (which Makie 0.24 cannot render)
  are plain text; `main()` runs only from the command line.
- `visualisation/view_h5_output.jl --animate` uses the mesh embedded in a
  format-2 file; `--mesh` is needed only for format-1 files.

## [0.4.2] — 2026-11-11

Test-case registry (PATCH). Simulation results do not change. The output's
`/summary` gains two attributes, and its volumes are now filled in.

### Added
- **Test-case registry** (`cases/`, module `FloodA5Cases`; `docs/TEST_CASES.md`):
  - `julia --project=. cases/run.jl list | show | prep | run | check |
    baseline CASE`.
  - Each case generates its own inputs (DEM, meshes, hydrographs) in
    `test/case_runs/<case>/`, runs `FloodModel.jl` with explicit flags, and
    computes its metrics with `analysis/`.
  - It then checks physical expectations and compares with the expected
    values in `cases/expected/<case>.toml`, recorded with the version and
    commit that produced them.
  - A run is reused when its output records the same command line.
  - Cases:
    - `point_spread_smoke` (CI)
    - `point_spread`: A5 res 16/18, Cartesian 25/12.5 m
    - `planar_slope`: N/S symmetry on a 0.1% plane
    - `planar_embankment`: SGS passes water through a sub-cell notch near
      its sill, 1.8 m below the crest; standard flow only at its cell-mean
      crest
    - `channel`: 16 runs of the channel series
    - `channel_series`: the full 224 runs; its jobs and meshes reproduce
      the recorded series bit for bit
    - `carlisle`: 6 h, local-inertial / SWE storage 1.25 on A5 and 0.99 on
      Cartesian, as recorded in HYDRAULICS §15.11
  - `test/test_cases.jl` (in CI) runs the smoke tier end to end.
- `/summary` attributes `vol_applied_m3` (what the time stepping added) and
  `conservation_error_m3` (= applied − stored − removed).
  `mass_balance_error_m3` (= exact input − stored − removed, as the run log
  reports it) also contains the sources' first-order time integration: rates
  are evaluated at the start of each step. That is 4.7% on a 600 s triangular
  hydrograph at default Δt, and zero for constant rates.

### Fixed
- ⚠ The end-of-run mass balance was never computed. `/summary`'s
  `vol_added_m3`, `vol_removed_m3`, `vol_domain_m3` and
  `mass_balance_error_m3` were 0 in every 0.4.0/0.4.1 output, and the timing
  report's mass-balance table was never printed.
- `FA.storage` and `FA.compare` reported depth·area as storage for SGS runs.
  SGS depth is measured from the cell's lowest point, so the store was
  overstated (+24% on the embankment case); they now use the stored volume.
- Mesh generation used fixed temporary file names in `mesh/`, so two mesh
  generations in one checkout could corrupt each other. The names are now
  unique per process.

### Changed
- `test/planar_embankment/test_planar_embankment.jl` is removed, replaced by
  the `planar_embankment` case. It was broken by design: its DEM had been
  regenerated without the embankment for the planar-slope tests; it measured
  the upstream water level over dry cells too; and its source was too small
  for its depth check.
- Twelve superseded run scripts moved to `test/archive/` (with a README).
  They rely on the pre-0.3.0 defaults and no longer reproduce their recorded
  numbers; the README maps each to its registry case.

## [0.4.1] — 2026-11-11

Analysis tooling only (PATCH). Neither the model nor the output format changes.

### Added
- `analysis/FloodA5Analysis.jl` (module `FloodA5Analysis`, `FA`): one module
  for loading runs and computing metrics, in pure Julia. It does not include
  `FloodModel.jl` and needs no Python; it loads in about 3 s.
  - `Run` / `MeshGeom`: frame times, mesh geometry and run record; fields by
    frame (`depth`, `wse` SGS-aware, `speed`, …); frames at matched times;
    storage; local coordinates.
  - Format-1 (≤ 0.3.x) runs: geometry from the mesh's GeoParquet
    (`load_run(p; parquet = …)`), read through ArchGDAL's Parquet driver and
    matching the model's areas, ring convention and neighbour order. Runs on
    a re-labelled mesh copy are matched by row (`match = :row`).
  - Metrics: run vs reference (storage ratio, rms Δh, CSI, bias, deep-water
    speed) and twin-run divergence; point-spread and planar-slope shape and
    direction; checkerboard and velocity noise; straight-channel series.
- `test/test_analysis.jl` (in CI) with `test/analysis/test_*.jl`: the
  metrics checked on synthetic fields with known answers; the GeoParquet
  mesh checked against a format-2 file's mesh.

### Changed
- The analysis scripts are now thin command-line wrappers over `FA`, and
  none of them includes `FloodModel.jl` any more: `wp1_metrics.jl`,
  `diag_vs_reference.jl`, `diag_divergence.jl`, `wp2_checkerboard.jl`,
  `channel_bench/analyse_channel.jl`, `test_point_spread.jl`,
  `test_octagonal_anisotropy.jl`, `test_planar_slope.jl` and
  `test_planar_symmetry.jl`. Same flags, same printed output and CSVs on
  existing runs. The channel series' `metrics.csv` is byte-identical over
  all 336 runs, and the analysis takes 3m41 instead of 5m44. The mesh
  argument of `diag_vs_reference.jl` and `wp2_checkerboard.jl` is optional
  for format-2 runs.

## [0.4.0] — 2026-11-10

Output format change (MINOR). Simulation results are unchanged: every
time step is computed as in 0.3.3; only the HDF5 file layout and its
contents change.

### Changed
- ⚠ **HDF5 output format 2** (`docs/DATA_FORMATS.md` §2). Time series are
  one dataset per variable, `/results/<var>` (n_frames × n_cells in C
  order), instead of one `/frames/NNNNNN` group per output time. Against
  format 1 on existing runs (`test/bench_output_format.jl`):
  - time series of cells read 2.4–28× faster
  - whole variables read 3.0–5.8× faster
  - file size within +2.5% on 25-frame runs; 4.6× smaller with 241 frames
  - writing 4–11% slower, tens of ms per frame

  Scripts that read `/frames` directly no longer work on new files; read
  through `io/FloodA5Output.jl`, which reads both formats. Re-run earlier
  runs to get format-2 files.

### Added
- The output records how to reproduce its run:
  - `/run`: version, git commit (and the patch if the tree was dirty), the
    command line, every option after the `--scheme` defaults, backend,
    threads, cell order, and the environment's `Project.toml` and
    `Manifest.toml`
  - the input files (hydrographs, `.bci`/`.bdy`, `--bc-file`, friction
    raster; files over 64 MB by checksum only)
  - the mesh: its GeoParquet file (`/mesh/source`), cell polygons, neighbour
    lists and cell areas
  - `/summary` at the end: completed, steps, mass balance, wall time
- `--rerun FILE.h5 [flags…]`: runs the recorded simulation again from the
  embedded mesh and inputs. Flags given with it replace the recorded ones.
  Same code and backend: bit-identical frames.
- `--embed-mesh full|geometry` (default `full`): `geometry` stores the
  polygons and the mesh file's checksum but not the file, for sweeps on a
  large SGS mesh (Carlisle res-18 SGS: 179 MB).
- `io/FloodA5Output.jl` (module `FloodA5Output`, `FO`): writer, reader
  (`read_frame`, `read_var`, `read_series`, `frame_times`, `cell_polygons`,
  `run_info`, …) and reproduction (`extract_mesh`, `rerun_argv`), without
  depending on `FloodModel.jl`. All scripts that read output now use it.
- A mesh without elevation data is recorded with the z = 0 the model runs it
  at (format 1 wrote NaN).
- `run_simulation!` called directly (not through `run_flood_model`) creates
  its output file, with the mesh but no run record, replacing any file at
  that path. Format 1 needed the file to exist already and appended to it.
- `test/test_output_format.jl` (O1–O6, in CI): run record and geometry,
  `--rerun` bit-identical with the original mesh and inputs moved away (also
  for a mesh generated in the run, and with `--embed-mesh geometry`), the
  format-1 reader, and a direct `run_simulation!` call.
- `SHA` (standard library) as a dependency.

## [0.3.3] — 2026-11-10

Documentation and housekeeping only (PATCH); no code changes.

### Changed
- README "Development Roadmap": what was done in 0.2.1–0.3.2, and the next
  phase (consolidation: self-describing HDF5 output, an analysis module, a
  test-case registry, HDF5 replay in the GLMakie viewer, notebooks).
- `CLAUDE.md` status updated; `projectinstructions` (supplementary repo)
  updated to its 2026-11-10 `PROJECT_STATE.md` revision.

## [0.3.2] — 2026-11-02

Default results unchanged (PATCH): a benchmark series and documentation only.

### Added
- `test/channel_bench/`: the straight-channel representation series — a
  20 m channel at 7 orientations and 2 offsets on A5 (res 17–20) and
  equal-area Cartesian grids, with SWE and the local-inertial default
  (mesh and job generation, metrics, figures). `docs/TEST_CASES.md` §7.

### Findings (see `HYDRAULICS.md` §15.12, figures in `docs/figures/`)
- Both meshes lose conveyance in the same way as the channel narrows
  below a cell (stage error about 1.1 m at 44 m cells, ≤ 0.15 m at 5.6 m).
- A5 is not orientation-free. Cartesian is worst at 45°, A5 at 0–15° (at
  51° N, 0° E). The spread over orientations is as large on A5 or larger,
  and neither mesh is better overall. Locally, A5 faces lie in six fixed
  directions.
- On these in-bank cases the local-inertial default is as close to normal
  depth as SWE, or closer.

## [0.3.1] — 2026-11-02

Default results unchanged (PATCH): the CPU path is bit-identical to 0.3.0
(Carlisle 6 h, A5 and Cartesian standard flow and A5 SGS, `--cell-order
file`, every frame).

### Added
- `--backend cpu|gpu|ka-cpu` (WP5): every flow model (standard, SGS, SWE)
  runs its time steps on CPU threads or an NVIDIA GPU from one code base.
  Each pass is a kernel body launched by `foreach_index!` — Polyester on
  the CPU, KernelAbstractions on a device; the state stays on the device
  and is copied back for output and log lines. CUDA.jl loads only for
  `--backend gpu`; `ka-cpu` runs the device path on the CPU, for testing.
  The device path needs edge momentum and the local CFL (the defaults).
  Open-boundary outflow (Phase D) runs in dependency levels that reproduce
  the serial CPU order. `docs/ARCHITECTURE.md` §6a.
- GPU throughput against 16 CPU threads (all cells wet, 10⁶ cells): standard
  flow 93 vs 8.4, SWE 29 vs 2.7 million cell-steps per second; the GPU is
  faster from about 10⁵ cells (`test/benchmark_backend.jl`; per-phase
  timing in `test/profile_backend.jl`).
- `test/test_backend.jl` (K1–K5, in CI): the device path on the CPU is
  bit-identical to CPU threads for all three solvers; on a GPU
  (`FLOODA5_TEST_GPU=1`) single steps agree to 10⁻¹⁵. Over a 6 h Carlisle
  run the GPU matches the CPU to 3·10⁻¹¹ m on Cartesian; on A5 the
  solvers amplify round-off at wet/dry fronts, and the GPU run departs from
  the CPU run exactly as much as a CPU run in another cell order does
  (storage within 3·10⁻⁴).
- KernelAbstractions and Adapt as dependencies.

### Fixed
- Start-up time with the diamond face gradient (the default): building its
  vertex table was quadratic in the number of vertices (51 s for 65 000
  Cartesian cells, hours for 10⁶); now linear (4 s). Same table.

## [0.3.0] — 2026-10-30

⚠ Default results change for `--flow-model standard` (MINOR), and for
`--flow-model sgs` through its open boundaries. SWE is unchanged.

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

### Fixed
- ⚠ Open (ZeroGradient) boundaries extrapolated the ghost water level from
  the most *outward* neighbour (a comparison with the wrong sign) — one
  along the boundary — so outlets held water back and narrow domains leaked
  through their sides. `--open-boundary normal` (the default for standard
  and SGS flow under `--scheme recommended`) extrapolates the bed slope
  from the inward neighbour at the boundary cell's depth: a free outflow at
  normal depth. `--open-boundary legacy` (and `--scheme v0.2`/`legacy`)
  keeps the old behaviour. FlowState field 43 `open_boundary`, with a
  42-argument constructor shim. Tests T-BC11/12.

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
