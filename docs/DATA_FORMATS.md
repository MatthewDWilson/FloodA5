# FloodA5 — Data Formats

This document specifies all file formats used by FloodA5: GeoParquet mesh schema,
HDF5 simulation output, hydrograph input formats, and the CesiumJS binary wire
protocol.

---

## 1. GeoParquet Mesh File

Produced by `--meshout` — via `mesh/a5_bridge.py` for A5 meshes (`--meshgen`),
or via the pure-Julia generator `mesh/CartesianGrid.jl` for Cartesian meshes
(`--meshgen-cartesian`, `cartesian-formulation` branch; the final GeoParquet
*write* still goes through the same existing Python-side serialiser in both
cases — "no Python bridge" for the Cartesian backend means no Python is needed
to *generate* the geometry, not that file I/O avoids Python). Loaded by
`--meshload`, which doesn't distinguish between backends — everything
downstream of load is identical. Extension: `.parquet`.

### Geometry column

| Property | Value |
|---|---|
| Column name | `geometry` |
| Encoding | WKB via geopandas |
| CRS | EPSG:4326 (WGS 84, longitude/latitude) |
| Shape | A5: pentagon polygon — 5 vertices + closing repeat (6 coordinate pairs). Cartesian: square polygon — 4 vertices + closing repeat (5 coordinate pairs), counter-clockwise. |

### Required scalar columns

| Column | dtype | Description |
|---|---|---|
| `cell_id` | string | 16-character zero-padded hex ID. A5: the real A5 cell ID, e.g. `08a2a1072b59ffff`. Cartesian: same *shape* but a different meaning — `(row << 32) \| col`, hex-encoded, so e.g. row 40 col 6 is `0000002800000006`; trivially invertible, not a real A5 cell. The loader normalises both identically (see `A5_QUIRKS.md` §2); a non-hex ID from some future backend would also pass through unmodified rather than failing, since the normalisation falls back to the raw string. |
| `center_lon` | float64 | Cell centre longitude (degrees) |
| `center_lat` | float64 | Cell centre latitude (degrees) |
| `resolution` | int32 | A5 resolution level. `0` for Cartesian meshes (the concept doesn't apply — see the `cell_size_m` static column below for the actual cell size; `.resolution` is only ever read for display/logging, never branched on) |

### Static variable columns (added by FloodA5)

| Column | dtype | Added by | Description |
|---|---|---|---|
| `elevation` | float64 | `--dem` | Mean cell bed elevation (m). `NaN` if outside DEM extent. |
| `sgs_cell_area` | float64 | SGS build | Geodetic plan area of cell polygon (m²) |
| `sgs_z_min` | float64 | SGS build | Minimum DEM sample within cell (m) |
| `sgs_z_max` | float64 | SGS build | Maximum DEM sample within cell (m) |
| `sgs_n_bins` | float64 | SGS build | Number of hypsometric bins (same for all cells) |
| `cell_size_m` | float64 | Cartesian generator | Nominal cell edge length (m) used to generate the mesh — present only on Cartesian meshes. Constant across all cells in a given mesh, stored per-cell for consistency with how every other per-cell scalar is stored. |

### SGS array columns

Stored as Apache Arrow list columns (variable-length arrays). Each row contains a
1-D array. Length notation: `n_bins` = number of elevation bins (default 100);
`max_nb` = 5 (maximum neighbours per cell). **SGS is currently A5-only** —
`build_sgs_tables!` still hardcodes `max_nb = 5` in three places and has not
yet been generalised to a mesh-derived neighbour ceiling the way the
standard-flow path has (see `FlowState.max_neighbours` in `ARCHITECTURE.md`).
A Cartesian mesh should therefore be treated as **unsupported** with
`--flow-model sgs` — not tested, and the hardcoded 5-slot matrices would not
match a 4-neighbour mesh. (There is currently no explicit guard blocking the
combination; it isn't prevented, just not expected to work.) See the
Development Roadmap in `README.md`.

| Column | Length | Description |
|---|---|---|
| `sgs_elev_bins` | n_bins | Elevation bin knot points (m), quantile-spaced from z_min to z_max |
| `sgs_vol_curve` | n_bins | Cumulative stored volume at each bin edge (m³) |
| `sgs_area_curve` | n_bins | Wetted plan area at each bin elevation (m²) |
| `sgs_edge_sills` | max_nb | Minimum DEM elevation along the shared edge to each neighbour (m), indexed by adjacency slot |
| `sgs_edge_area_curve` | n_bins × max_nb | Cross-sectional flow area at each elevation knot, per adjacency slot (m²). Stored as a flat array of length n_bins × 5; reshape to (n_bins, 5) on load. |
| `sgs_edge_perim_curve` | n_bins × max_nb | Wetted perimeter at each elevation knot, per adjacency slot (m). Same layout as `sgs_edge_area_curve`. |

### Cell-ID scheme metadata (WP-M1)

The Parquet schema metadata carries a `flooda5` key with a JSON object, e.g.
`{"a5_id_scheme": "a5-v1", "pya5_version": "1.0.0"}` (plus `migrated_from` /
`migrated_by` after `mesh/migrate_mesh_ids.py`). `a5_id_scheme` is `a5-v0`
(pya5 0.x cell IDs), `a5-v1` (pya5 1.x) or `cartesian`. Files without the key
were written before WP-M1 by pya5 0.7.x. The same cell has different IDs in the
two A5 schemes, so FloodA5 refuses to load a mesh of the other scheme — see
`A5_QUIRKS.md` §11.

### Adjacency column

| Column | dtype | Description |
|---|---|---|
| `neighbours` | list of string | Neighbour cell IDs (edge-sharing adjacency, up to 5 on A5 or 4 on Cartesian — i.e. up to the mesh's own `max_neighbours`). Variable-length — boundary cells with fewer neighbours than the mesh's maximum simply have a shorter list, rather than a fixed 5-slot layout with empty placeholders. Cell ID normalisation (16-char zero-padded hex) is applied on load; see `A5_QUIRKS.md` §2. |

### Reading in Python

```python
import geopandas as gpd
import pyarrow.parquet as pq

# Geometry and scalar columns
gdf = gpd.read_file("mesh.parquet")

# All columns including array columns
tbl = pq.read_table("mesh.parquet")
sgs_vol = tbl["sgs_vol_curve"].to_pylist()   # list of lists (one per cell)
```

### Reading in Julia

```julia
mesh = A5Grid.load_mesh_geoparquet("mesh.parquet")
# mesh.static_vars["elevation"]          → Vector{Float64}
# mesh.array_vars["sgs_vol_curve"]       → Matrix{Float64} (n_bins × n_cells)
```

---

## 2. HDF5 Simulation Output

Written by `--output`. Extension: `.h5` or `.hdf5`. Since 0.4.0 the file is
**format 2**: it records everything needed to reproduce the run
(`--rerun FILE.h5`) and to replay or analyse it without the original mesh
file. The writer and the reader live in one module, `io/FloodA5Output.jl`;
read files through it (`FO.read_frame`, `FO.read_var`, …, below) rather than
by dataset path, so scripts keep working across format versions. Files from
FloodA5 ≤ 0.3.x (format 1, `/frames/NNNNNN` groups, no run record or
geometry) are still readable through the same functions.

```
/                 attrs: format = "FloodA5 output", format_version = 2,
                         flooda5_version, created_utc
/run              the run record
/mesh             static mesh data, including the mesh file itself
/results          the time series: one dataset per variable
/summary          attrs written when the run ends (absent if it did not finish)
```

### `/run` — run record

| Item | Kind | Description |
|---|---|---|
| `argv` | dataset, JSON array | The command line, as given (empty for runs not started from the CLI) |
| `options` | dataset, JSON object | Every option after the `--scheme` defaults are applied (flow model, scheme, limiter, Courant number, η filter, open boundary, sources, …) |
| `git_commit`, `git_dirty` | attrs | Code version; `git_dirty` = uncommitted changes to tracked files |
| `git_diff` | dataset, text | The uncommitted patch (only when `git_dirty`) |
| `project`, `manifest` | datasets, text | `Project.toml` and `Manifest.toml` of the active environment (exact package versions) |
| `julia_version`, `hostname`, `cwd`, `backend`, `threads`, `cell_order`, `flow_model`, `started_utc` | attrs | Environment and run settings |
| `inputs/NN/` | groups | One per input file the run read (`--friction`, `--inflow-point` hydrographs, `--inflow-bci` and its `.bdy`, `--bc-file`). Attrs `role`, `flag`, `argv_value` (as written in `argv`; empty for a `.bdy` found next to the `.bci`), `original_path`, `sha256`, `size`, `embedded`; dataset `bytes` (the file) when `embedded`. Files over 64 MB are referenced by path and checksum only. |

The DEM is not stored: the mesh (below) already holds what the simulation
used from it (cell elevations, SGS tables).

### `/mesh` — static data

| Dataset | Shape (C order) | dtype | Description |
|---|---|---|---|
| `cell_ids` | (n_cells,) | string | Cell IDs. Group attr `id_scheme` (`a5-v0`, `a5-v1` or `cartesian`) |
| `center_lons`, `center_lats` | (n_cells,) | float64 | Cell centre (degrees) |
| `elevations` | (n_cells,) | float64 | Bed elevation (m); `NaN` if not sampled |
| _static variables_ | (n_cells,) | float64 | The mesh's other per-cell scalars, e.g. `friction`, `sgs_z_min` |
| `cell_area` | (n_cells,) | float64 | Plan area (m²) |
| `polygon_offsets` | (n_cells+1,) | int64 | Cell *i*'s vertices are `polygon_lon/lat[offsets[i]:offsets[i+1]]` (0-based, CSR) |
| `polygon_lon`, `polygon_lat` | (n_vertices,) | float64 | Polygon vertices as stored in the mesh file (degrees) |
| `neighbours` | (n_cells, max_nb) | int32 | Neighbour cells in neighbour-list order, 0-based, `-1` = none |
| `edge_cells` | (n_edges, 2) | int32 | The two cells of each edge, 0-based (only with `flux_Q`) |
| `source/` | group | | The mesh file. Attrs `embedded` (`full` or `geometry`), `kind` (`file`: the loaded or `--meshout` file; `resaved`: a mesh generated in the run and not saved, saved for the output), `original_path`, `sha256`, `size`; datasets `parquet` (the GeoParquet file's bytes; `full` only) and `prj` |

All cell arrays are in the **mesh file's cell order** (the simulation may
reorder cells internally, `--cell-order spatial`; output is written back in
file order). `--embed-mesh geometry` leaves out `source/parquet` — for
sweeps on one large SGS mesh, where the copy would dominate the file — and
`--rerun` then uses the original mesh file after checking its checksum.

### `/results` — time series

| Dataset | Shape (C order) | dtype | Description |
|---|---|---|---|
| `t` | (n_frames,) | float64 | Simulation time (s) |
| `water_depth` | (n_frames, n_cells) | float64 | Depth above local bed (m) |
| `volume` | (n_frames, n_cells) | float64 | Stored water volume (m³) — primary state variable |
| `saturation` | (n_frames, n_cells) | float64 | Wetted fraction 0–1 (SGS: meaningful; standard: binary) |
| `velocity` | (n_frames, n_cells) | float64 | Velocity magnitude (m/s): WLSQ cell-vector reconstruction of face discharges divided by depth (since WP0, 2026-10; earlier files were not in m/s) |
| `vel_u`, `vel_v` | (n_frames, n_cells) | float64 | Eastward / northward velocity (m/s) |
| `flux_Q` | (n_frames, n_edges) | float64 | Flux per edge (m³/s), SGS runs with the R-A kernel only (`HYDRAULICS.md` §9.3); frames before its first non-zero value read as 0 |

In Julia the shapes read transposed: `(n_cells, n_frames)`. Datasets grow
by one frame per output time; chunks are one frame × up to 4096 cells,
deflate level 4. A whole frame is one contiguous read, and a cell's time
series one small read per frame. Against format 1 on existing runs
(`test/bench_output_format.jl`): time series 2–30× and whole-variable reads
3–6× faster; files within 3% of the size, or several times smaller with many
frames; writing 4–11% slower (tens of ms per frame).

### `/summary`

Attrs `completed`, `n_frames`, `n_steps`, `wall_simulation_s`,
`wall_since_start_s`, `vol_added_m3`, `vol_removed_m3`, `vol_domain_m3`,
`mass_balance_error_m3`, `vol_applied_m3`, `conservation_error_m3`,
`finished_utc`.

- `vol_added_m3`: the exact time integral of the inputs (rainfall,
  hydrographs, point sources); `vol_applied_m3`: what the time stepping
  actually added (Σ rate·Δt, rates evaluated at the start of each step).
  They differ by the sources' first-order time integration (a few % on a
  steep 10-minute hydrograph ramp; negligible for constant rates).
- `conservation_error_m3` = applied − in domain − removed: conservation
  alone (round-off). `mass_balance_error_m3` = added − in domain − removed,
  as in the run log, includes the time-integration difference.
- Before 0.4.2 the volume attrs and the mass-balance error were written as 0.

### Reading in Julia

```julia
include("io/FloodA5Output.jl"); const FO = FloodA5Output   # (FloodModel.jl defines FO)

nf    = FO.n_frames("sim.h5");  t = FO.frame_times("sim.h5")
h_end = FO.read_frame("sim.h5", "water_depth")             # last frame; or (…, k)
h     = FO.read_var("sim.h5", "water_depth")               # n_cells × n_frames
ts    = FO.read_series("sim.h5", "water_depth", [12, 40])  # n_frames × 2
k     = FO.frame_at("sim.h5", 3600.0)                      # frame nearest t = 1 h
cells = FO.cell_table("sim.h5")                            # ids, lon, lat, elevation
polys = FO.cell_polygons("sim.h5")                         # 2 × n_vertices per cell
info  = FO.run_info("sim.h5")                              # argv, options, git, inputs, summary
```

Every function also takes an open `HDF5.File`. Indices are 1-based.

### Reproducing a run

```
julia --project=. FloodModel.jl --rerun sim.h5                   # → sim_rerun.h5
julia --project=. FloodModel.jl --rerun sim.h5 --output b.h5 --backend gpu
```

The recorded command line is run again with the embedded mesh
(`--meshload`) and the embedded input files; mesh-generation flags are
dropped, and any flag given with `--rerun` replaces the recorded one. With
the same code (commit) and backend the frames are bit-identical
(`test/test_output_format.jl` O2–O4); a warning is printed when the code
differs from the recorded commit.

### Reading in Python

```python
import h5py, json
import numpy as np

with h5py.File("sim.h5") as f:
    t      = f["results/t"][:]
    depth  = f["results/water_depth"]            # (n_frames, n_cells)
    final  = depth[-1, :]                        # one frame
    gauge  = depth[:, 1234]                      # one cell's time series
    lons   = f["mesh/center_lons"][:]
    lats   = f["mesh/center_lats"][:]
    offs   = f["mesh/polygon_offsets"][:]
    plon, plat = f["mesh/polygon_lon"][:], f["mesh/polygon_lat"][:]
    poly0  = np.column_stack([plon[offs[0]:offs[1]], plat[offs[0]:offs[1]]])
    argv   = json.loads(f["run/argv"][()])
    opts   = json.loads(f["run/options"][()])
```

With xarray:

```python
import xarray as xr
with h5py.File("sim.h5") as f:
    ds = xr.Dataset(
        {"water_depth": (["time", "cell"], f["results/water_depth"][:])},
        coords={"time": f["results/t"][:],
                "lon": ("cell", f["mesh/center_lons"][:]),
                "lat": ("cell", f["mesh/center_lats"][:])})
```

### Format 1 (FloodA5 ≤ 0.3.x)

`/mesh` with `cell_ids`, `elevations`, `center_lons`, `center_lats` (and
static variables); one group per output time, `/frames/000001/`, …, each
holding `t` and the per-cell variables above as `(n_cells,)` datasets
(`flux_Q` only in frames where it was non-zero). No run record, geometry or
summary. `io/FloodA5Output.jl` reads it; re-run such runs with 0.4.0 to get
a format-2 file.

---

## 3. Hydrograph Input Formats

### 3.1 Two-column CSV

The simplest format for `--inflow-point`. Two columns: simulation time in seconds
and discharge in m³/s. An optional header row is auto-detected (if the first row
is non-numeric, it is skipped).

```
t_s,Q_m3s
0,0.0
3600,5.2
7200,18.1
10800,42.7
```

Plain (no header):

```
0    0.0
3600 5.2
```

Delimiter: comma or whitespace, auto-detected.

### 3.2 LISFLOOD-FP `.bdy` file

One or more named time series in a single file. Each series begins with a name line,
followed by a count, a time unit, and then paired rows of time and discharge. Time
units are `seconds`, `hours`, or `days`; FloodA5 converts all to seconds on read.
Comments (lines beginning `#` or `!`) and blank lines are ignored.

```
UPSTREAM_GAUGE
241
hours
0.0    0.0
1.0    5.2
2.0    18.1
...

TRIBUTARY
48
hours
0.0    0.0
1.0    0.8
...
```

When multiple series are present in a `.bdy` file, use the LABEL argument to
`--inflow-point` to select the desired series by name. `--inflow-bci` with `QVAR`
entries references series by name automatically.

### 3.3 LISFLOOD-FP `.bci` file

Specifies boundary condition entries for the domain. FloodA5 supports `P`-type
(point source) entries and `FREE` entries. Column layout:

| Col | Description |
|---|---|
| 1 | Entry type: `P` (point source), `N/E/S/W` (unsupported — see below), `F` (unsupported) |
| 2 | Longitude (or easting if `--bc-epsg` is set) |
| 3 | Latitude (or northing) |
| 4 | BC code: `QVAR`, `QFIX`, `FREE`, or `HFIX`/`HVAR` (see below) |
| 5 | For `QVAR`: series name in the companion `.bdy` file. For `QFIX`: discharge in m³/s. |

**`HFIX`/`HVAR` (fixed/time-varying water surface elevation — tide-style
boundary) entries are recognised by the parser but not yet implemented.**
They are logged with a warning and skipped rather than silently ignored or
misapplied. This is planned future work (see the roadmap in the top-level
`README.md`).

Example:

```
P   -2.8957   54.9090   QVAR   main
P   -2.9372   54.8838   QVAR   cummersdale
P   -2.9146   54.8843   QFIX   1.5
```

FloodA5 looks for a `.bdy` file with the same name stem as the `.bci` file in the
same directory. Use `--inflow-bdy` to specify a different path.

`N`, `E`, `S`, `W` (cardinal edge boundaries) are not supported, because FloodA5
domains can be arbitrary polygons with no axis-aligned boundaries. A helpful warning
is logged if these are encountered, directing the user to use `--bc-file` (GeoJSON)
instead.

---

## 4. GeoJSON Boundary Condition File (`--bc-file`)

A GeoJSON `FeatureCollection` where each feature is a `LineString` or `Polygon`
that marks a boundary segment. Boundary cells whose centres fall within 1.5× the
cell diameter of each feature geometry are assigned the specified BC type.

| Property | Value | Description |
|---|---|---|
| `bc_type` | `Closed`, `ZeroGradient`, or `Critical` | BC type to apply to matched cells |
| `label` | any string (optional) | Used in log output for traceability |

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "LineString",
        "coordinates": [[172.55, -43.60], [172.60, -43.60]]
      },
      "properties": {
        "bc_type": "Closed",
        "label": "Northern levee"
      }
    },
    {
      "type": "Feature",
      "geometry": {
        "type": "LineString",
        "coordinates": [[172.70, -43.55], [172.75, -43.55]]
      },
      "properties": {
        "bc_type": "ZeroGradient",
        "label": "Eastern outlet"
      }
    }
  ]
}
```

Cells not matched by any feature use the default BC type (`ZeroGradient` unless
`--closed-boundaries` is set).

---

## 5. AOI GeoJSON Input

Used with `--meshgen`. Must be a valid GeoJSON `Feature` or `FeatureCollection`
with a `Polygon` geometry in EPSG:4326.

```json
{
  "type": "Feature",
  "properties": {},
  "geometry": {
    "type": "Polygon",
    "coordinates": [[
      [172.55, -43.60],
      [172.75, -43.60],
      [172.75, -43.45],
      [172.55, -43.45],
      [172.55, -43.60]
    ]]
  }
}
```

The polygon must be closed (first and last coordinate pair identical). Coordinates
are `[longitude, latitude]` in decimal degrees (GeoJSON convention). Multi-polygon
AOIs are not currently supported.

---

## 6. CesiumJS Binary Wire Protocol

The visualisation server serves the CesiumJS viewer via HTTP and WebSocket. The
key design principle is **one variable fetched at a time** — only the actively
displayed variable crosses the wire, scaling to 1M+ cells.

### HTTP endpoints

| Endpoint | Response | Description |
|---|---|---|
| `GET /viz/{file}` | `text/html` / JSON | Static files: `index.html`, `config.json` |
| `GET /mesh` | JSON | GeoJSON FeatureCollection + `cell_order` array |
| `GET /frames/count` | JSON | `{count: N, vars: [...]}` |
| `GET /frames/{idx}` | JSON | `{t: 123.0, vars: [...]}` — metadata only |
| `GET /frames/{idx}/{varname}` | `application/octet-stream` | Raw Float32 LE array, n_cells × 4 bytes |
| `GET /status` | JSON | Server diagnostics |

### Binary frame format

`GET /frames/{idx}/{varname}` returns a raw binary body:

- Encoding: little-endian `float32`
- Length: exactly `n_cells × 4` bytes
- Ordering: matches the `cell_order` array from `/mesh` — index `i` in the Float32
  array corresponds to the cell at `cell_order[i]`

**JavaScript (receive):**
```javascript
const resp = await fetch(`/frames/${idx}/depth`);
const buf  = await resp.arrayBuffer();
const vals = new Float32Array(buf);   // vals[i] = depth of cell_order[i]
```

### WebSocket messages

The WebSocket endpoint (`WS /live`) pushes notifications from the server.

| `type` field | Payload | Sent when |
|---|---|---|
| `"mesh"` | `{type, data: GeoJSON, cell_order: [...]}` | On connect |
| `"framecount"` | `{type, count: N, vars: [...]}` | On connect (after mesh) |
| `"newframe"` | `{type, idx: N, t: 123.0, vars: [...]}` | Each new frame |
| `"simcomplete"` | `{type, frames: N}` | Simulation end |

`newframe` carries only metadata — the client fetches binary data on demand based
on the currently displayed variable.

---

## 7. Cesium Viewer Configuration

```json
{
  "cesium_ion_token": "YOUR_TOKEN_FROM_ion.cesium.com",
  "home_lon": 0.0,
  "home_lat": 0.0,
  "home_alt": 25000,
  "default_basemap": "google3d"
}
```

| Key | Default | Description |
|---|---|---|
| `cesium_ion_token` | `""` | Cesium Ion token for Google 3D Photorealistic Tiles. Falls back to OSM if absent. |
| `home_lon` | 0.0 | Camera home longitude |
| `home_lat` | 0.0 | Camera home latitude |
| `home_alt` | 25000 | Camera home altitude (metres) |
| `default_basemap` | `"google3d"` | `"google3d"` / `"osm"` / `"none"` |

Copy `visualisation/cesium/config.example.json` to `visualisation/cesium/config.json`
and add `config.json` to `.gitignore`. The viewer validates the token at startup and
falls back gracefully if it is missing.
