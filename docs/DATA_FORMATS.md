# FloodA5 — Data Formats Reference

_Specification of all file formats used by FloodA5: GeoParquet mesh schema, HDF5 output layout, and the CesiumJS binary wire protocol. Intended as project context for visualisation, post-processing, and I/O development conversations._

---

## 1. GeoParquet Mesh File

Produced by `save_mesh_geoparquet` (via `a5_bridge.py`). Read by `load_mesh_geoparquet` into an `A5Mesh` struct.

### Geometry column
- Column name: `geometry`
- Type: GeoParquet polygon geometry (WKB-encoded via geopandas)
- CRS: EPSG:4326 (WGS84 geographic, lon/lat)
- Each row: a single pentagon polygon (5 vertices + closing repeat = 6 coordinate pairs)

### Required scalar columns

| Column | dtype | Description |
|--------|-------|-------------|
| `cell_id` | string (hex) | 16-character zero-padded hex A5 cell ID (e.g. `"08a2a1072b59ffff"`) |
| `center_lon` | float64 | Cell centre longitude (degrees, EPSG:4326) |
| `center_lat` | float64 | Cell centre latitude (degrees, EPSG:4326) |
| `resolution` | int32 | A5 resolution level (e.g. 14) |

### Optional scalar columns (added by FloodA5)

| Column | dtype | Added by | Description |
|--------|-------|----------|-------------|
| `elevation` | float64 | `sample_dem_*` | Mean cell bed elevation (m above datum). NaN if out-of-DEM. |
| `sgs_cell_area` | float64 | `build_sgs_tables!` | Geodetic plan area of cell polygon (m²) |
| `sgs_z_min` | float64 | `build_sgs_tables!` | Minimum DEM sample within cell (m) |
| `sgs_z_max` | float64 | `build_sgs_tables!` | Maximum DEM sample within cell (m) |
| `sgs_n_bins` | float64 | `build_sgs_tables!` | Number of hypsometric bins (scalar, same for all cells) |

### Array columns (SGS hypsometric tables)

Stored as Apache Arrow list columns (variable-length arrays). Each row contains a 1-D array of length `n_bins` (or `n_bins+1` for bin edges).

| Column | Array length | Description |
|--------|-------------|-------------|
| `sgs_elev_bins` | n_bins + 1 | Elevation bin edges (m) — quantile-spaced from z_min to z_max |
| `sgs_vol_curve` | n_bins | Cumulative stored volume at each bin edge (m³) |
| `sgs_area_curve` | n_bins | Wetted plan area at each bin elevation (m²) |
| `sgs_edge_sills` | max_nb (= 5) | Minimum DEM elevation along shared edge to each neighbour (m) — indexed by adjacency slot |

### Adjacency columns

| Column | dtype | Description |
|--------|-------|-------------|
| `adj_0` … `adj_4` | string (hex) | Neighbour cell IDs in each of the up to 5 adjacency slots. Empty string if slot unused (boundary cells have <5 neighbours). |

### Reading in Python
```python
import geopandas as gpd
import pyarrow.parquet as pq

gdf = gpd.read_file("mesh.parquet")          # geometry + scalar columns
tbl = pq.read_table("mesh.parquet")          # all columns including arrays
sgs_vol = tbl["sgs_vol_curve"].to_pylist()  # list of lists
```

### Reading in Julia
```julia
mesh = A5Grid.load_mesh_geoparquet("mesh.parquet")
# mesh.static_vars["elevation"]          → Vector{Float64}
# mesh.array_vars["sgs_vol_curve"]       → Matrix{Float64} (n_bins × n_cells)
# mesh.adjacency                         → Dict{String, Vector{String}}
```

---

## 2. HDF5 Simulation Output

Produced by `_write_frame!` and `_write_mesh_metadata!`. Extension `.h5` or `.hdf5`.

### `/mesh` group — static data written once at simulation start

| Dataset | shape | dtype | Description |
|---------|-------|-------|-------------|
| `/mesh/cell_ids` | (n_cells,) | string | 16-char zero-padded hex cell IDs |
| `/mesh/elevations` | (n_cells,) | float64 | Bed elevation (m) |
| `/mesh/center_lons` | (n_cells,) | float64 | Cell centre longitude (degrees) |
| `/mesh/center_lats` | (n_cells,) | float64 | Cell centre latitude (degrees) |

Additional static variables from `mesh.static_vars` may be written if present.

### `/frames` group — time series

Each timestep snapshot is stored as a numbered subgroup:

```
/frames/000001/   t=60.0 s
/frames/000002/   t=120.0 s
...
```

| Dataset | shape | dtype | Description |
|---------|-------|-------|-------------|
| `t` | scalar | float64 | Simulation time (seconds) |
| `water_depth` | (n_cells,) | float64 | Depth above local bed (m) |
| `volume` | (n_cells,) | float64 | Stored water volume (m³) — primary state |
| `saturation` | (n_cells,) | float64 | Wetted fraction 0–1 (SGS: meaningful; standard: binary) |
| `velocity` | (n_cells,) | float64 | Scalar velocity magnitude (m/s) — currently always zero |

Frame groups are zero-padded to 6 digits. Chunk size: `min(n_cells, 4096)`. Compression: gzip level 4.

### Reading in Python
```python
import h5py, numpy as np

with h5py.File("sim.h5") as f:
    cell_ids = f["mesh/cell_ids"][:].astype(str)
    lons     = f["mesh/center_lons"][:]
    lats     = f["mesh/center_lats"][:]
    
    frames = sorted(f["frames"].keys())
    depths = np.stack([f[f"frames/{fr}/water_depth"][:] for fr in frames])
    times  = np.array([f[f"frames/{fr}/t"][()] for fr in frames])

# depths shape: (n_frames, n_cells)
```

### Reading in Julia
```julia
using HDF5

h5open("sim.h5", "r") do f
    cell_ids = read(f["mesh/cell_ids"])
    frames   = sort(keys(f["frames"]))
    depth_0  = read(f["frames/$(frames[1])/water_depth"])
end
```

### Reading with xarray (recommended for analysis)
```python
import xarray as xr, h5py, numpy as np

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

## 3. CesiumJS Binary Wire Protocol

The VisualisationServer serves the CesiumJS viewer via HTTP and WebSocket. The key design principle is **one variable fetched at a time** — only the actively displayed variable crosses the wire, scaling to 1M+ cells.

### HTTP endpoints

| Endpoint | Response | Description |
|----------|----------|-------------|
| `GET /viz/{file}` | `text/html` or `application/json` | Static files: `index.html`, `config.json` |
| `GET /mesh` | `application/json` | GeoJSON FeatureCollection + `cell_order` array |
| `GET /frames/count` | `application/json` | `{count: N, vars: ["depth", "saturation", ...]}` |
| `GET /frames/{idx}` | `application/json` | `{t: 123.0, vars: ["depth", ...]}` — metadata only |
| `GET /frames/{idx}/{varname}` | `application/octet-stream` | Raw Float32 LE array, n_cells × 4 bytes |
| `GET /status` | `application/json` | Server diagnostics |

### Binary frame format

`GET /frames/{idx}/{varname}` returns a raw binary body:
- Encoding: little-endian `float32` (`Float32` in Julia, `Float32Array` in JS)
- Length: exactly `n_cells × 4` bytes
- Ordering: matches the `cell_order` array from `/mesh` — index `i` in the Float32 array corresponds to the cell at `cell_order[i]`

**Julia side (push):**
```julia
data = Float32.(state.water_depth)   # n_cells Float32 values
body = Vector{UInt8}(reinterpret(UInt8, data))   # raw bytes
```

**JavaScript side (receive):**
```javascript
const resp = await fetch(`/frames/${idx}/depth`);
const buf  = await resp.arrayBuffer();
const vals = new Float32Array(buf);   // vals[i] = depth of cell_order[i]
```

### WebSocket messages (JSON)

The WebSocket endpoint (`WS /live`) pushes notifications from the server. The client does not send messages.

| `type` field | Payload | When sent |
|-------------|---------|-----------|
| `"mesh"` | `{type, data: GeoJSON, cell_order: [...]}` | On WebSocket connect |
| `"framecount"` | `{type, count: N, vars: [...]}` | On connect (after mesh) |
| `"newframe"` | `{type, idx: N, t: 123.0, vars: [...]}` | Each time `push_frame!` is called |
| `"simcomplete"` | `{type, frames: N}` | When `notify_complete!` is called |

The `newframe` message contains only metadata — the client decides whether to fetch the binary data based on its current variable selection and whether a render is already in progress.

### Variable definitions (client-side)

```javascript
const VAR_DEFS = {
    depth:      { colormap: "turbo",  fixedMax: null },
    saturation: { colormap: "blues",  fixedMax: 1.0  },
    volume:     { colormap: "turbo",  fixedMax: null },
    velocity:   { colormap: "plasma", fixedMax: null },
};
```

Unknown variables sent by the server are handled by a `varDef(name)` fallback that returns default colormap settings — new variables appear in the dropdown automatically without client changes.

---

## 4. GeoJSON Mesh Format (alternative to GeoParquet)

When `--meshout` is given a `.geojson` extension, the mesh is saved as a GeoJSON FeatureCollection. Less efficient than GeoParquet for large meshes (no columnar compression) but human-readable and compatible with any GIS tool.

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": { "type": "Polygon", "coordinates": [[[lon, lat], ...]] },
      "properties": {
        "cell_id": "08a2a1072b59ffff",
        "center_lon": 172.636,
        "center_lat": -43.531,
        "resolution": 14,
        "elevation": 3.4
      }
    }
  ]
}
```

GeoJSON meshes do **not** include SGS array columns (hypsometric tables) — these require Parquet's columnar array support. SGS runs must use `.parquet`.

---

## 5. Config File (`viz/config.json`)

Gitignored. Copy from `viz/config.example.json`:

```json
{
  "cesium_ion_token": "YOUR_TOKEN_FROM_ion.cesium.com",
  "home_lon": 172.636,
  "home_lat": -43.531,
  "home_alt": 25000,
  "default_basemap": "google3d"
}
```

| Key | Default | Description |
|-----|---------|-------------|
| `cesium_ion_token` | `""` | Cesium Ion token for Google 3D Photorealistic Tiles. Viewer falls back gracefully if missing. |
| `home_lon` | 172.636 | Camera home longitude (Christchurch) |
| `home_lat` | -43.531 | Camera home latitude |
| `home_alt` | 25000 | Camera home altitude (metres) |
| `default_basemap` | `"google3d"` | `"google3d"` / `"osm"` / `"none"` |

---

## 6. AOI GeoJSON Input

Used with `--meshgen`. Must be a valid GeoJSON `Feature` or `FeatureCollection` with a `Polygon` geometry in EPSG:4326.

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

The polygon should be closed (first and last coordinate identical). Coordinates are `[longitude, latitude]` (GeoJSON convention). Multi-polygon AOIs are not currently supported — use the outer ring only.
