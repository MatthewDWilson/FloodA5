# FloodA5 — Data Formats

This document specifies all file formats used by FloodA5: GeoParquet mesh schema,
HDF5 simulation output, hydrograph input formats, and the CesiumJS binary wire
protocol.

---

## 1. GeoParquet Mesh File

Produced by `--meshout` (via `mesh/a5_bridge.py`). Loaded by `--meshload`.
Extension: `.parquet`.

### Geometry column

| Property | Value |
|---|---|
| Column name | `geometry` |
| Encoding | WKB via geopandas |
| CRS | EPSG:4326 (WGS 84, longitude/latitude) |
| Shape | Pentagon polygon — 5 vertices + closing repeat (6 coordinate pairs) |

### Required scalar columns

| Column | dtype | Description |
|---|---|---|
| `cell_id` | string | 16-character zero-padded hex A5 cell ID, e.g. `08a2a1072b59ffff` |
| `center_lon` | float64 | Cell centre longitude (degrees) |
| `center_lat` | float64 | Cell centre latitude (degrees) |
| `resolution` | int32 | A5 resolution level |

### Static variable columns (added by FloodA5)

| Column | dtype | Added by | Description |
|---|---|---|---|
| `elevation` | float64 | `--dem` | Mean cell bed elevation (m). `NaN` if outside DEM extent. |
| `sgs_cell_area` | float64 | SGS build | Geodetic plan area of cell polygon (m²) |
| `sgs_z_min` | float64 | SGS build | Minimum DEM sample within cell (m) |
| `sgs_z_max` | float64 | SGS build | Maximum DEM sample within cell (m) |
| `sgs_n_bins` | float64 | SGS build | Number of hypsometric bins (same for all cells) |

### SGS array columns

Stored as Apache Arrow list columns (variable-length arrays). Each row contains a
1-D array. Length notation: `n_bins` = number of elevation bins (default 100);
`max_nb` = 5 (maximum neighbours per cell).

| Column | Length | Description |
|---|---|---|
| `sgs_elev_bins` | n_bins | Elevation bin knot points (m), quantile-spaced from z_min to z_max |
| `sgs_vol_curve` | n_bins | Cumulative stored volume at each bin edge (m³) |
| `sgs_area_curve` | n_bins | Wetted plan area at each bin elevation (m²) |
| `sgs_edge_sills` | max_nb | Minimum DEM elevation along the shared edge to each neighbour (m), indexed by adjacency slot |
| `sgs_edge_area_curve` | n_bins × max_nb | Cross-sectional flow area at each elevation knot, per adjacency slot (m²). Stored as a flat array of length n_bins × 5; reshape to (n_bins, 5) on load. |
| `sgs_edge_perim_curve` | n_bins × max_nb | Wetted perimeter at each elevation knot, per adjacency slot (m). Same layout as `sgs_edge_area_curve`. |

### Adjacency column

| Column | dtype | Description |
|---|---|---|
| `neighbours` | list of string | Neighbour cell IDs (edge-sharing adjacency, up to 5). Variable-length — boundary cells with fewer than 5 neighbours simply have a shorter list, rather than a fixed 5-slot layout with empty placeholders. Cell ID normalisation (16-char zero-padded hex) is applied on load; see `A5_QUIRKS.md` §2. |

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

Written by `--output`. Extension: `.h5` or `.hdf5`.

### `/mesh` group — static data (written once at simulation start)

| Dataset | Shape | dtype | Description |
|---|---|---|---|
| `/mesh/cell_ids` | (n_cells,) | string | 16-char zero-padded hex cell IDs |
| `/mesh/elevations` | (n_cells,) | float64 | Bed elevation (m). `NaN` if not sampled. |
| `/mesh/center_lons` | (n_cells,) | float64 | Cell centre longitude (degrees) |
| `/mesh/center_lats` | (n_cells,) | float64 | Cell centre latitude (degrees) |

### `/frames` group — time series

Each snapshot is a numbered subgroup:

```
/frames/000001/   t=60.0 s
/frames/000002/   t=120.0 s
...
```

| Dataset | Shape | dtype | Description |
|---|---|---|---|
| `t` | scalar | float64 | Simulation time (seconds) |
| `water_depth` | (n_cells,) | float64 | Depth above local bed (m) |
| `volume` | (n_cells,) | float64 | Stored water volume (m³) — primary state variable |
| `saturation` | (n_cells,) | float64 | Wetted fraction 0–1 (SGS: meaningful; standard: binary) |
| `velocity` | (n_cells,) | float64 | Scalar velocity magnitude (m/s) |
| `vel_u` | (n_cells,) | float64 | Eastward velocity component (m/s) |
| `vel_v` | (n_cells,) | float64 | Northward velocity component (m/s) |
| `flux_Q` | (n_edges,) | float64 | Volumetric flux per edge (m³/s). Only present in the frame group when at least one edge is non-zero — i.e. only for SGS runs using the R-A flux kernel (see `HYDRAULICS.md` §9.3). Absent (not zero-filled) otherwise. |

Datasets are chunked (`min(n_cells, 4096)`) and gzip-compressed (level 4). Frame
groups are zero-padded to six digits.

### Reading in Python

```python
import h5py
import numpy as np

with h5py.File("sim.h5") as f:
    # Static mesh
    cell_ids = f["mesh/cell_ids"][:].astype(str)
    lons     = f["mesh/center_lons"][:]
    lats     = f["mesh/center_lats"][:]

    # All frames
    frames = sorted(f["frames"].keys())
    times  = np.array([f[f"frames/{fr}/t"][()] for fr in frames])
    depths = np.stack([f[f"frames/{fr}/water_depth"][:] for fr in frames])

# depths.shape = (n_frames, n_cells)
print(depths[-1].max())   # peak depth in final frame
```

### Reading with xarray (recommended for analysis)

```python
import xarray as xr
import h5py
import numpy as np

with h5py.File("sim.h5") as f:
    frames = sorted(f["frames"].keys())
    times  = [f[f"frames/{fr}/t"][()] for fr in frames]
    lons   = f["mesh/center_lons"][:]
    lats   = f["mesh/center_lats"][:]
    depths = np.stack([f[f"frames/{fr}/water_depth"][:] for fr in frames])

ds = xr.Dataset(
    {"water_depth": (["time", "cell"], depths)},
    coords={"time": times, "lon": ("cell", lons), "lat": ("cell", lats)}
)
```

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
