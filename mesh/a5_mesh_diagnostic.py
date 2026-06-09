"""
a5_mesh_diagnostic.py
---------------------
Standalone Python diagnostic for A5 mesh generation and neighbour computation.

Usage:
    python a5_mesh_diagnostic.py <aoi.geojson> <resolution> [output.parquet]

Key design decisions:
  Cell generation: sample-grid (or fill_polygon if available) + grid_disk
    expansion pass. AOI inclusion test uses shapely.intersects (ST_Intersects
    equivalent) so boundary-straddling cells are included.
  Adjacency: shared-boundary-vertex detection, NOT grid_disk. Two cells are
    edge-sharing neighbours iff they share exactly 2 boundary vertices.
    This avoids pya5 grid_disk compaction artefacts which cause some cells to
    report fewer than 5 neighbours.
"""

import sys
import json
from collections import defaultdict, Counter
from typing import List, Dict, Set, Tuple

try:
    from a5 import lonlat_to_cell, cell_to_lonlat, cell_to_boundary, u64_to_hex, get_resolution
    print("[OK] pya5 base functions available")
except ImportError as e:
    print(f"[ERROR] pya5: {e}"); sys.exit(1)

try:
    from a5 import grid_disk
    HAS_GRID_DISK = True
    print("[OK] pya5.grid_disk available")
except ImportError:
    HAS_GRID_DISK = False
    print("[WARN] pya5.grid_disk not available")

try:
    from a5 import uncompact
    HAS_UNCOMPACT = True
except ImportError:
    HAS_UNCOMPACT = False

try:
    from a5 import geometry as a5_geometry
    _ = a5_geometry.fill_polygon
    HAS_FILL_POLYGON = True
    print("[OK] pya5.geometry.fill_polygon available")
except Exception:
    HAS_FILL_POLYGON = False
    print("[WARN] pya5.geometry.fill_polygon not available - will use sample grid")

try:
    import geopandas as gpd
    from shapely.geometry import Polygon
    import shapely, pandas as pd
    print("[OK] geopandas / shapely available")
except ImportError as e:
    print(f"[ERROR] {e}"); sys.exit(1)


# ── Geometry helpers ──────────────────────────────────────────────────────────

def to_hex(cell_id) -> str:
    return f"{int(cell_id) & 0xFFFFFFFFFFFFFFFF:016x}"

def wrap_lon(lon):
    return ((float(lon) + 180) % 360) - 180

def load_aoi(path):
    with open(path, "r", encoding="utf-8") as f:
        gj = json.load(f)
    t = gj.get("type", "")
    if t == "FeatureCollection": return gj["features"][0]["geometry"]
    if t == "Feature":           return gj["geometry"]
    if t in ("Polygon", "MultiPolygon"): return gj
    raise ValueError(f"Unsupported GeoJSON type: {t}")

def bbox(geometry):
    t = geometry["type"]
    rings = geometry["coordinates"] if t == "Polygon" else \
            [r for poly in geometry["coordinates"] for r in poly]
    coords = [c for ring in rings for c in ring]
    return (min(c[0] for c in coords), min(c[1] for c in coords),
            max(c[0] for c in coords), max(c[1] for c in coords))

def make_shapely(geometry):
    try:
        return shapely.from_geojson(json.dumps(geometry))
    except Exception:
        from shapely.geometry import shape as _shape
        return _shape(geometry)

def cell_polygon(cell_id) -> Polygon:
    bnd = cell_to_boundary(cell_id)
    return Polygon([(wrap_lon(v[0]), float(v[1])) for v in bnd])

def lonlat_to_cell_safe(lon, lat, resolution):
    try:    return int(lonlat_to_cell(lon, lat, resolution)) & 0xFFFFFFFFFFFFFFFF
    except TypeError: return int(lonlat_to_cell([lon, lat], resolution)) & 0xFFFFFFFFFFFFFFFF

def cell_res(cell_id: int) -> int:
    try:    return int(get_resolution(cell_id))
    except: return -1


# ── Cell generation ───────────────────────────────────────────────────────────

def generate_cells(geometry, resolution) -> List[int]:
    """
    Generate all cells at `resolution` that intersect the AOI.

    Two-pass approach:
    Pass 1 (broad): sample-grid (or fill_polygon) using a buffered AOI to seed
      initial cells. The buffer = ~1 cell width ensures we don't miss any cell
      whose centre is just outside the tight AOI boundary.
    Pass 2 (expansion): iteratively add grid_disk(c,1) neighbours that:
      - are at the correct resolution, and
      - their polygon intersects the ORIGINAL AOI (shapely.intersects)
    This combination catches every cell with any overlap of the AOI.
    """
    aoi_shp  = make_shapely(geometry)
    min_lon, min_lat, max_lon, max_lat = bbox(geometry)

    # Approximate cell width in degrees for buffer distance
    import math
    cell_deg = 90.0 / (2 ** (resolution * 1.16))
    buf_deg  = cell_deg * 1.5   # ~1.5 cell widths

    def intersects_aoi(cid: int) -> bool:
        """True if cell polygon intersects the (unbuffered) AOI."""
        try:
            return bool(aoi_shp.intersects(cell_polygon(cid)))
        except Exception:
            return False

    # ── Pass 1: seed cells from buffered AOI sample ───────────────────────
    if HAS_FILL_POLYGON and HAS_UNCOMPACT:
        t = geometry["type"]
        coords = geometry["coordinates"][0] if t == "Polygon" else \
                 geometry["coordinates"][0][0]
        try:
            compact = a5_geometry.fill_polygon(coords, resolution)
            cells = set(int(c) & 0xFFFFFFFFFFFFFFFF
                        for c in uncompact(compact, resolution))
            print(f"  fill_polygon + uncompact: {len(cells)} seed cells")
        except Exception as e:
            print(f"  fill_polygon failed ({e}), using sample grid")
            cells = _sample_grid(min_lon - buf_deg, min_lat - buf_deg,
                                 max_lon + buf_deg, max_lat + buf_deg,
                                 geometry, resolution, buf_deg)
    else:
        cells = _sample_grid(min_lon - buf_deg, min_lat - buf_deg,
                             max_lon + buf_deg, max_lat + buf_deg,
                             geometry, resolution, buf_deg)

    if not cells:
        return []
    print(f"  Seed cells: {len(cells)}")

    # ── Pass 2: expand via grid_disk + ST_Intersects ──────────────────────
    if not HAS_GRID_DISK:
        print("  [WARN] grid_disk not available — skipping expansion")
        return list(cells)

    rounds = 0
    while True:
        candidates: Set[int] = set()
        for c in cells:
            try:
                for nb in grid_disk(c, 1):
                    nb_id = int(nb) & 0xFFFFFFFFFFFFFFFF
                    if nb_id != c and nb_id not in cells:
                        candidates.add(nb_id)
            except Exception:
                pass

        new_cells = set()
        for cid in candidates:
            if cell_res(cid) != resolution:
                continue  # wrong resolution (compaction artefact)
            if intersects_aoi(cid):
                new_cells.add(cid)

        if not new_cells:
            break
        cells |= new_cells
        rounds += 1
        print(f"  Expansion round {rounds}: +{len(new_cells)} cells (total {len(cells)})")

    # Final filter: keep only cells that actually intersect the AOI
    before = len(cells)
    cells = {c for c in cells if intersects_aoi(c)}
    if len(cells) < before:
        print(f"  Final AOI filter removed {before - len(cells)} non-overlapping cells")

    print(f"  Final cell count: {len(cells)}")
    return list(cells)


def _sample_grid(min_lon, min_lat, max_lon, max_lat, geometry, resolution, step_override=None):
    """Sample a grid of points and collect containing cells."""
    step = step_override if step_override else \
           max(0.0001, 90.0 / (2 ** (resolution * 1.16))) * 0.45

    aoi_shp = make_shapely(geometry)

    lons, lats = [], []
    lat = min_lat
    while lat <= max_lat + step:
        lon = min_lon
        while lon <= max_lon + step:
            lons.append(max(-179.9, min(lon, 179.9)))
            lats.append(min(lat, 89.9))
            lon += step
        lat += step

    print(f"  Sample grid: {len(lons)} points (step={step:.5f}°)")

    try:
        inside = shapely.contains_xy(aoi_shp, lons, lats).tolist()
    except Exception:
        # Buffered version: include all points in bbox
        inside = [True] * len(lons)

    n_in = sum(inside)
    print(f"  Points inside: {n_in}")
    if n_in == 0:
        # Try all points in the buffered bbox
        print("  [WARN] No points inside AOI — using all bbox points as fallback")
        inside = [True] * len(lons)

    cells = set()
    for lo, la, ok in zip(lons, lats, inside):
        if not ok: continue
        try:
            cells.add(lonlat_to_cell_safe(lo, la, resolution))
        except Exception:
            pass
    return cells


# ── Adjacency via shared vertices ─────────────────────────────────────────────

def build_adjacency_shared_vertices(cells: List[int]) -> Dict[int, List[int]]:
    """
    Build exact edge-sharing adjacency using shared-vertex detection.

    Two cells are edge-sharing neighbours iff their boundaries share
    exactly 2 vertices (the two endpoints of the shared edge).

    This is geometrically correct and avoids pya5 grid_disk compaction
    artefacts, which cause some cells to return fewer than 5 neighbours.

    Algorithm:
      1. Fetch all cell boundaries.
      2. Build a vertex -> [cell_ids] index (rounded to 8 decimal places
         to handle floating-point differences in pya5 output).
      3. For each cell, find all cells sharing any vertex.
      4. Keep those sharing exactly 2 vertices = true edge-sharing neighbours.
    """
    PREC = 7   # decimal places for vertex rounding (~1cm precision)

    # Fetch boundaries
    print(f"  Fetching boundaries for {len(cells)} cells...")
    boundaries: Dict[int, List[Tuple[float, float]]] = {}
    for c in cells:
        try:
            bnd = cell_to_boundary(c)
            # Store as (lon, lat) tuples, normalised, excluding closing vertex
            verts = [(round(wrap_lon(v[0]), PREC), round(float(v[1]), PREC))
                     for v in bnd[:-1]]   # drop closing duplicate
            boundaries[c] = verts
        except Exception as e:
            boundaries[c] = []

    # Build vertex → cell index
    vertex_to_cells: Dict[Tuple, List[int]] = defaultdict(list)
    for c, verts in boundaries.items():
        for v in verts:
            vertex_to_cells[v].append(c)

    # Find neighbours: cells sharing exactly 2 vertices
    adj: Dict[int, List[int]] = {}
    for c in cells:
        verts = set(boundaries[c])
        # Candidate neighbours: cells that share at least one vertex
        candidates: Dict[int, int] = defaultdict(int)
        for v in verts:
            for nb in vertex_to_cells[v]:
                if nb != c:
                    candidates[nb] += 1
        # Keep those sharing exactly 2 vertices (one shared edge)
        adj[c] = [nb for nb, count in candidates.items() if count == 2]

    return adj


# ── Diagnostics ───────────────────────────────────────────────────────────────

def print_diagnostics(cells, adj, grid_disk_raw=None):
    print()
    print("=" * 60)
    print("DIAGNOSTIC SUMMARY")
    print("=" * 60)
    if not cells:
        print("[!] Zero cells."); return

    sublattices = Counter(to_hex(c)[0] for c in cells)
    print(f"\nTotal cells: {len(cells)}")
    print("Cells per sublattice (leading hex nibble):")
    for nibble, cnt in sorted(sublattices.items()):
        print(f"  '{nibble}...': {cnt} ({100*cnt/len(cells):.1f}%)")

    adj_counts = Counter(len(v) for v in adj.values())
    print(f"\nAdjacency (shared-vertex method):")
    for k, v in sorted(adj_counts.items()):
        print(f"  {k} neighbours: {v} cells")

    zero = [c for c, v in adj.items() if not v]
    if zero:
        print(f"\n[!] {len(zero)} cells have ZERO neighbours:")
        for c in zero:
            print(f"    {to_hex(c)}")
    else:
        print(f"\n[OK] All cells have >= 1 neighbour.")

    total_edges = sum(len(v) for v in adj.values()) // 2
    print(f"\nTotal undirected edges: {total_edges}")
    expected = len(cells) * 5 // 2
    print(f"Expected if all interior (~5n/2): {expected}")

    if grid_disk_raw:
        raw_counts = Counter(len(v) for v in grid_disk_raw.values())
        print(f"\ngrid_disk(c,1) raw counts (for comparison):")
        for k, v in sorted(raw_counts.items()):
            print(f"  {k} raw: {v} cells")
        cells_set = set(cells)
        cross = sum(1 for c in cells for nb in grid_disk_raw.get(c,[]) if nb in cells_set)
        total_raw = sum(len(v) for v in grid_disk_raw.values())
        print(f"  {cross}/{total_raw} ({100*cross/max(total_raw,1):.0f}%) raw neighbours in mesh")

    print("=" * 60)


# ── Save ──────────────────────────────────────────────────────────────────────

def save_output(cells, adj, out_path, grid_disk_raw=None):
    if not cells:
        print("[SKIP] No cells."); return
    rows = []
    for c in cells:
        try:
            hx  = to_hex(c)
            ct  = cell_to_lonlat(c)
            lon = wrap_lon(ct[0])
            lat = max(-90.0, min(90.0, float(ct[1])))
            bnd = cell_to_boundary(c)
            coords = [(wrap_lon(v[0]), float(v[1])) for v in bnd]
            geom = Polygon(coords)
            nbrs = [to_hex(nb) for nb in adj.get(c, [])]
            row = {
                "cell_id"      : hx,
                "sublattice"   : hx[0],
                "center_lon"   : lon,
                "center_lat"   : lat,
                "neighbours"   : nbrs,
                "n_neighbours" : len(nbrs),
                "geometry"     : geom,
            }
            if grid_disk_raw:
                raw = [to_hex(nb) for nb in grid_disk_raw.get(c, [])]
                row["neighbours_raw"]   = raw
                row["n_neighbours_raw"] = len(raw)
            rows.append(row)
        except Exception as e:
            print(f"  [WARN] {to_hex(c)}: {e}")

    if not rows:
        print("[ERROR] No rows."); return

    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    gdf.to_parquet(out_path, index=False)
    print(f"[OK] {len(rows)} cells -> {out_path}")

    vis = out_path.replace(".parquet", "_visual.geojson")
    drop = ["neighbours"]
    if "neighbours_raw" in gdf.columns: drop.append("neighbours_raw")
    gdf.drop(columns=drop).to_file(vis, driver="GeoJSON")
    print(f"[OK] Visual -> {vis}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 3:
        print("Usage: python a5_mesh_diagnostic.py <aoi.geojson> <resolution> [out.parquet]")
        sys.exit(1)
    aoi_path   = sys.argv[1]
    resolution = int(sys.argv[2])
    out_path   = sys.argv[3] if len(sys.argv) > 3 else \
                 aoi_path.replace(".geojson", f"_mesh{resolution}_diag.parquet")

    print(f"\nA5 Mesh Diagnostic  |  {aoi_path}  |  res={resolution}\n")
    geometry = load_aoi(aoi_path)

    print("Step 2: Generating cells (two-pass: buffered seed + ST_Intersects expansion)...")
    cells = generate_cells(geometry, resolution)

    # Also get grid_disk raw for comparison
    grid_disk_raw = None
    if HAS_GRID_DISK and cells:
        print(f"\nStep 3: grid_disk raw neighbours (for comparison only)...")
        grid_disk_raw = {}
        cell_set = set(cells)
        for c in cells:
            try:
                grid_disk_raw[c] = [int(nb) & 0xFFFFFFFFFFFFFFFF
                                     for nb in grid_disk(c, 1)
                                     if (int(nb) & 0xFFFFFFFFFFFFFFFF) != c]
            except Exception:
                grid_disk_raw[c] = []

    print(f"\nStep 4: Building adjacency via shared-vertex detection...")
    adj = build_adjacency_shared_vertices(cells)

    print_diagnostics(cells, adj, grid_disk_raw)

    print("\nStep 5: Saving...")
    save_output(cells, adj, out_path, grid_disk_raw)

if __name__ == "__main__":
    main()
