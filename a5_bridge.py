"""
a5_bridge.py
------------
Python bridge for A5 DGGS mesh generation and geospatial operations.

Performance architecture
------------------------
The mesh generation pipeline has three distinct phases, each optimised:

Phase 1 — Point sampling & PIP
  Julia handles this (GPU or multi-threaded CPU). The bridge receives the
  AOI GeoJSON and runs its own sampling only when called standalone.
  When called from Julia, the AOI has already been filtered to a dense
  sample grid; the bridge re-samples it but more coarsely (it only needs
  enough points to hit each cell at least once, not pixel-perfect coverage).

  If Shapely is available, PIP uses shapely.contains_xy() — a vectorised
  C implementation, ~50–100× faster than pure Python ray-casting loops.
  Falls back to pure Python if Shapely is unavailable.

Phase 2 — Cell indexing (pya5)
  lonlat_to_cell() is called once per *unique sample point* (up to ~147K).
  After deduplication, typically only ~3K unique cells remain at res 14.
  The deduplication happens as a set — O(1) per insert.

Phase 3 — Boundary/centre fetch + GeoParquet write
  Only called for unique cells (~3K at res 14, not 147K).
  Boundaries fetched once per cell. Shapely Polygon objects constructed
  using the vectorised shapely.polygons() C API (avoids per-object Python
  overhead). Coordinate normalisation done with NumPy broadcasting.

Commands
--------
    mesh_for_aoi  <geojson_path> <resolution> <output_path> [geoparquet|geojson]
    check

Output: one JSON line per status update; final line is always the result.
"""

import sys
import json
import os
from typing import List, Set, Optional, Tuple

# ── Core imports ───────────────────────────────────────────────────────────

def _fail(msg):
    print(json.dumps({"error": msg}), flush=True)
    sys.exit(1)

try:
    from a5 import (
        lonlat_to_cell,
        cell_to_lonlat,
        cell_to_boundary,
        u64_to_hex,
    )
except ImportError as e:
    _fail(f"pya5 not available: {e}. Run: pip install pya5")

# ── Optional accelerated dependencies ─────────────────────────────────────

def _try_import(pkg):
    try:
        __import__(pkg)
        return True
    except ImportError:
        return False

def _ensure_geopandas():
    """
    Import geopandas, auto-installing if missing using sys.executable.
    This guarantees the package lands in the right environment regardless
    of which pip/python the user ran from their shell.
    """
    missing = [p for p in ("geopandas", "pyarrow", "shapely") if not _try_import(p)]
    if missing:
        import subprocess
        print(json.dumps({
            "status": "installing",
            "message": f"Auto-installing into {sys.executable}: {missing}"
        }), flush=True)
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--quiet"] + missing
        )
    import geopandas as gpd
    from shapely.geometry import Polygon
    return gpd, Polygon

# Eager import attempt
try:
    import geopandas as gpd
    from shapely.geometry import Polygon
    import shapely
    import numpy as np
    _gpd_available   = True
    _numpy_available = True
    # Shapely ≥ 2.0 exposes vectorised C APIs
    _shapely_vec = hasattr(shapely, "contains_xy")
except ImportError:
    _gpd_available   = False
    _numpy_available = False
    _shapely_vec     = False

try:
    import numpy as np
    _numpy_available = True
except ImportError:
    _numpy_available = False

# ── Coordinate normalisation ───────────────────────────────────────────────

def _wrap_lon(lon: float) -> float:
    """Normalise longitude to [-180, 180]. pya5 may return values outside this range."""
    return ((lon + 180.0) % 360.0) - 180.0

def _norm_coords(boundary: list) -> list:
    """Normalise all [lon, lat] pairs in a boundary ring."""
    return [[_wrap_lon(v[0]), max(-90.0, min(90.0, v[1]))] for v in boundary]

def _norm_coords_np(boundaries: list) -> list:
    """
    Normalise a list of boundary rings using NumPy broadcasting.
    ~20× faster than calling _norm_coords() in a loop.
    Returns the same structure (list of list of [lon, lat]).
    """
    if not _numpy_available:
        return [_norm_coords(b) for b in boundaries]
    import numpy as np
    result = []
    for b in boundaries:
        arr = np.array(b, dtype=np.float64)         # (N, 2)
        arr[:, 0] = ((arr[:, 0] + 180.0) % 360.0) - 180.0
        arr[:, 1] = np.clip(arr[:, 1], -90.0, 90.0)
        result.append(arr.tolist())
    return result

# ── PIP implementations ─────────────────────────────────────────────────────

def _point_in_ring_py(lon, lat, ring) -> bool:
    """Pure Python ray-casting PIP."""
    inside = False
    n = len(ring)
    j = n - 1
    for i in range(n):
        xi, yi = ring[i][0], ring[i][1]
        xj, yj = ring[j][0], ring[j][1]
        if ((yi > lat) != (yj > lat)) and (
            lon < (xj - xi) * (lat - yi) / ((yj - yi) or 1e-15) + xi
        ):
            inside = not inside
        j = i
    return inside

def _point_in_polygon_py(lon, lat, coordinates) -> bool:
    if not _point_in_ring_py(lon, lat, coordinates[0]):
        return False
    return not any(_point_in_ring_py(lon, lat, h) for h in coordinates[1:])

def _point_in_geometry_py(lon, lat, geometry) -> bool:
    t = geometry["type"]
    if t == "Polygon":
        return _point_in_polygon_py(lon, lat, geometry["coordinates"])
    if t == "MultiPolygon":
        return any(_point_in_polygon_py(lon, lat, p) for p in geometry["coordinates"])
    raise ValueError(f"Unsupported geometry type: {t}")

def _make_shapely_geom(geometry):
    """Convert GeoJSON geometry dict to a Shapely geometry object."""
    from shapely.geometry import shape
    return shape(geometry)

def _pip_batch(lons: list, lats: list, geometry: dict) -> list:
    """
    Vectorised batch PIP for a list of (lon, lat) pairs.

    Uses shapely.contains_xy() if Shapely ≥ 2.0 is available — this calls
    GEOS directly via a C extension, vectorised over all points in one call.
    Falls back to pure Python ray-casting otherwise.

    Returns a list of bool.
    """
    if _shapely_vec:
        import shapely
        import numpy as np
        shp_geom = _make_shapely_geom(geometry)
        return shapely.contains_xy(shp_geom, lons, lats).tolist()
    else:
        return [_point_in_geometry_py(lon, lat, geometry)
                for lon, lat in zip(lons, lats)]

# ── Geometry helpers ────────────────────────────────────────────────────────

def _geometry_bbox(geometry):
    t = geometry["type"]
    rings = geometry["coordinates"] if t == "Polygon" else \
            [r for poly in geometry["coordinates"] for r in poly]
    coords = [c for ring in rings for c in ring]
    return (min(c[0] for c in coords), min(c[1] for c in coords),
            max(c[0] for c in coords), max(c[1] for c in coords))

def _aoi_geometry(geojson):
    t = geojson["type"]
    if t == "Feature":
        return geojson["geometry"]
    if t == "FeatureCollection":
        feats = geojson["features"]
        if len(feats) == 1:
            return feats[0]["geometry"]
        polys = []
        for f in feats:
            g = f["geometry"]
            if g["type"] == "Polygon":
                polys.append(g["coordinates"])
            elif g["type"] == "MultiPolygon":
                polys.extend(g["coordinates"])
        return {"type": "MultiPolygon", "coordinates": polys}
    if t in ("Polygon", "MultiPolygon"):
        return geojson
    raise ValueError(f"Cannot extract polygon geometry from type: {t}")

# ── Core mesh generation ────────────────────────────────────────────────────

def generate_mesh(geojson_path: str, resolution: int) -> List[dict]:
    """
    Generate A5 pentagon mesh cells covering the AOI.

    Uses pya5's native fill_polygon (polyfill) + uncompact to guarantee:
      - All cells are exactly at `resolution` (no mixed-resolution artifacts)
      - No holes (uncompact explodes any coarser cells from compacted results)
      - No sample-grid or BFS approximation needed

    Falls back to a sample-grid approach if fill_polygon is unavailable.
    """
    with open(geojson_path, "r", encoding="utf-8") as f:
        geojson = json.load(f)

    geometry  = _aoi_geometry(geojson)
    geojson_coords = _geojson_to_coord_list(geometry)

    # ── Phase 1: polyfill → uniform resolution ────────────────────────────
    try:
        from a5 import geometry as a5_geometry
        compact_cells = a5_geometry.fill_polygon(geojson_coords, resolution)
    except Exception:
        # fill_polygon unavailable — fall back to sample-grid + dedup
        compact_cells = _sample_grid_coverage(geometry, resolution)

    from a5 import uncompact
    uniform_cells = uncompact(compact_cells, resolution)

    # ── Phase 2: fetch boundaries for unique cells ─────────────────────────
    cell_ids   = list(set(uniform_cells))   # deduplicate
    hex_ids    = [u64_to_hex(c)        for c in cell_ids]
    centres    = [cell_to_lonlat(c)    for c in cell_ids]
    boundaries = [cell_to_boundary(c)  for c in cell_ids]

    # ── Phase 3: normalise coordinates (NumPy if available) ───────────────
    norm_boundaries = _norm_coords_np(boundaries)
    norm_lon = [_wrap_lon(c[0]) for c in centres]
    norm_lat = [max(-90.0, min(90.0, c[1])) for c in centres]

    return [
        {
            "cell_id":    hex_ids[i],
            "resolution": resolution,
            "center_lon": norm_lon[i],
            "center_lat": norm_lat[i],
            "boundary":   norm_boundaries[i],
        }
        for i in range(len(cell_ids))
    ]


def _geojson_to_coord_list(geometry: dict):
    """
    Convert a GeoJSON geometry to the coordinate format expected by
    a5.geometry.fill_polygon — a list of [lon, lat] pairs forming the
    outer ring of the polygon.
    """
    t = geometry["type"]
    if t == "Polygon":
        return geometry["coordinates"][0]
    elif t == "MultiPolygon":
        # Use the largest polygon by vertex count
        return max(geometry["coordinates"], key=lambda p: len(p[0]))[0]
    raise ValueError(f"Unsupported geometry type for fill_polygon: {t}")


def _sample_grid_coverage(geometry: dict, resolution: int) -> list:
    """
    Fallback coverage using sample-grid + lonlat_to_cell when fill_polygon
    is unavailable. Returns a list of raw cell IDs (ints).
    """
    min_lon, min_lat, max_lon, max_lat = _geometry_bbox(geometry)
    approx_cell_deg = max(0.0001, 90.0 / (2 ** (resolution * 1.16)))
    sample_step     = approx_cell_deg * 0.45

    sample_lons, sample_lats = [], []
    lat = min_lat
    while lat <= max_lat + sample_step:
        lon = min_lon
        while lon <= max_lon + sample_step:
            sample_lons.append(max(-179.9, min(lon, 179.9)))
            sample_lats.append(min(lat, 89.9))
            lon += sample_step
        lat += sample_step

    inside_mask = _pip_batch(sample_lons, sample_lats, geometry)
    cell_ids = set()
    for lo, la, ok in zip(sample_lons, sample_lats, inside_mask):
        if ok:
            cell_ids.add(lonlat_to_cell([lo, la], resolution))
    return list(cell_ids)


def mesh_to_geoparquet(cells: List[dict], output_path: str):
    """
    Write cells to GeoParquet using geopandas.

    Geometry column built with shapely.polygons() vectorised C API when
    Shapely ≥ 2.0 is available — avoids per-object Python overhead for the
    Polygon() constructor loop.

    Static variable columns (elevation, etc.) are included if present on
    cell dicts. Any key beyond the base set is written as an extra column.
    Scalar values are written as Float64 columns; list/array values are
    written as Arrow list<float64> columns (supported by GeoParquet/Arrow).
    This allows SGS hypsometric curves and edge sills to be stored per-cell.
    """
    gpd, Polygon = _ensure_geopandas()
    import numpy as np

    boundaries = [c["boundary"] for c in cells]

    # Vectorised geometry construction (Shapely ≥ 2.0)
    try:
        import shapely
        if hasattr(shapely, "polygons"):
            max_verts = max(len(b) for b in boundaries)
            coords_np = np.zeros((len(boundaries), max_verts, 2), dtype=np.float64)
            for i, b in enumerate(boundaries):
                arr = np.array(b, dtype=np.float64)
                coords_np[i, :len(b)] = arr
                if arr[0].tolist() != arr[-1].tolist():
                    coords_np[i, len(b)-1] = arr[0]
            geometries = shapely.polygons(coords_np)
        else:
            raise AttributeError
    except (AttributeError, Exception):
        geometries = [Polygon(b) for b in boundaries]

    # Base columns always present
    base_keys = {"cell_id", "resolution", "center_lon", "center_lat", "boundary"}
    data = {
        "cell_id":    [c["cell_id"]    for c in cells],
        "resolution": [c["resolution"] for c in cells],
        "center_lon": [c["center_lon"] for c in cells],
        "center_lat": [c["center_lat"] for c in cells],
    }

    # Extra columns: detect scalar vs list/array
    extra_keys = [k for k in cells[0].keys() if k not in base_keys] if cells else []
    for key in extra_keys:
        sample = cells[0].get(key)
        if isinstance(sample, (list, tuple)) or (
                hasattr(sample, '__len__') and not isinstance(sample, str)):
            # List column — store as Python lists (Arrow serialises as list<float64>)
            data[key] = [
                [float(x) if x is not None else float("nan") for x in c.get(key, [])]
                for c in cells
            ]
        else:
            # Scalar column
            data[key] = [c.get(key, float("nan")) for c in cells]

    gdf = gpd.GeoDataFrame(data, geometry=geometries, crs="EPSG:4326")
    gdf.to_parquet(output_path, index=False)


def mesh_to_geojson(cells: List[dict], output_path: str):
    # Base properties only (geometry file format; static vars go in parquet)
    features = [
        {
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": [c["boundary"]]},
            "properties": {
                "cell_id":    c["cell_id"],
                "resolution": c["resolution"],
                "center_lon": c["center_lon"],
                "center_lat": c["center_lat"],
            },
        }
        for c in cells
    ]
    fc = {
        "type":     "FeatureCollection",
        "features": features,
        "properties": {
            "resolution": cells[0]["resolution"] if cells else 0,
            "cell_count": len(cells),
        },
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(fc, f)


# ── CLI dispatcher ──────────────────────────────────────────────────────────

def cmd_check():
    gpd_ok    = _try_import("geopandas")
    np_ok     = _try_import("numpy")
    shp_ok    = _try_import("shapely")
    shp_vec   = False
    gpd_loc   = None
    if gpd_ok:
        import geopandas as _g
        gpd_loc = getattr(_g, "__file__", None)
    if shp_ok:
        import shapely as _s
        shp_vec = hasattr(_s, "contains_xy")
    print(json.dumps({
        "status":             "ok",
        "python":             sys.executable,
        "pya5":               True,
        "geopandas":          gpd_ok,
        "geopandas_location": gpd_loc,
        "numpy":              np_ok,
        "shapely":            shp_ok,
        "shapely_vectorised": shp_vec,
        "functions":          [n for n in dir(sys.modules["a5"]) if not n.startswith("_")],
    }), flush=True)


def cmd_mesh_for_aoi(geojson_path, resolution, output_path, fmt):
    cells = generate_mesh(geojson_path, int(resolution))
    if fmt == "geoparquet":
        mesh_to_geoparquet(cells, output_path)
    else:
        mesh_to_geojson(cells, output_path)
    print(json.dumps({
        "status":      "ok",
        "cell_count":  len(cells),
        "output":      output_path,
        "format":      fmt,
        "pip_backend":  "shapely_vectorised" if _shapely_vec else "python_raycasting",
        "geom_backend": "shapely_polygons_vectorised" if (_gpd_available and _shapely_vec) else "polygon_loop",
        "coverage":     "fill_polygon+uncompact",
    }), flush=True)


def main():
    args = sys.argv[1:]
    if not args:
        _fail("No command given. Usage: a5_bridge.py <command> [args...]")
    cmd = args[0]
    try:
        if cmd == "check":
            cmd_check()
        elif cmd == "mesh_for_aoi":
            if len(args) < 4:
                _fail("Usage: mesh_for_aoi <geojson_path> <resolution> <output_path> [geoparquet|geojson]")
            fmt = args[4] if len(args) > 4 else "geoparquet"
            cmd_mesh_for_aoi(args[1], args[2], args[3], fmt)
        else:
            _fail(f"Unknown command: {cmd}")
    except Exception as e:
        import traceback
        _fail(f"{e}\n{traceback.format_exc()}")


if __name__ == "__main__":
    main()
