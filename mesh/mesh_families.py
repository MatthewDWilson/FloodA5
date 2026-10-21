"""
mesh_families.py
----------------
Generic polygon meshes for FloodA5 (WP3, 2026-10 — FloodA5_Review_and_ForwardPlan
_2026-10.md §6.1, §9 WP3). Writes FloodA5 GeoParquet meshes (the same format
as A5 meshes: cell_id, resolution, center_lon, center_lat, geometry, neighbours,
cell_area; schema metadata a5_id_scheme = "polygon") that load with
`--meshload`, so every solver option can be tested on mesh families other than
A5 and Cartesian:

    hex       regular hexagons (6 neighbours; pointy-top), side from --dx
    rotcart   Cartesian squares rotated by --angle degrees
    voronoi   Voronoi cells of jittered lattice points (--jitter, fraction of
              --dx); --lloyd N relaxes them towards a centroidal Voronoi
              tessellation. --centre generator (default) uses the Voronoi
              generators as cell centres, which makes every face exactly
              orthogonal to its centre-to-centre line; --centre centroid uses
              the polygon centroids instead (a non-orthogonal mesh)
    import    any polygon GeoParquet/GeoJSON (--input), reprojected to
              EPSG:4326 if it has a CRS; centres are polygon centroids

hex/rotcart/voronoi cover the AOI's bounding box (or --bbox) in a local
equirectangular plane about its centre: hex and rotcart keep whole cells whose
centre lies inside the box; voronoi cells are clipped to it. Cell IDs are the
16-hex-digit row index. Adjacency is built from shared polygon edges (vertex
pairs matched after rounding), so meshes with hanging nodes (T-junctions) are
not supported — such cells would get no neighbour across the split edge.

Usage (from the repo root, with the Python that FloodA5 uses):
    python mesh/mesh_families.py hex     --aoi AOI.geojson --dx 50 --out mesh.parquet
    python mesh/mesh_families.py voronoi --bbox LON0,LAT0,LON1,LAT1 --dx 50 \
                                         --jitter 0.3 --lloyd 5 --seed 1 --out m.parquet
    python mesh/mesh_families.py rotcart --aoi AOI.geojson --dx 50 --angle 30 --out m.parquet
    python mesh/mesh_families.py import  --input polygons.gpkg --out mesh.parquet
Prints one JSON line with a summary.
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

_R = 6_371_000.0          # matches A5Grid's equirectangular helpers
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ── Local plane ─────────────────────────────────────────────────────────────

class LocalPlane:
    """Equirectangular projection about (lon0, lat0), metres."""

    def __init__(self, lon0, lat0):
        self.lon0, self.lat0 = lon0, lat0
        self.kx = math.radians(1.0) * _R * math.cos(math.radians(lat0))
        self.ky = math.radians(1.0) * _R

    def to_xy(self, lon, lat):
        return (np.asarray(lon) - self.lon0) * self.kx, (np.asarray(lat) - self.lat0) * self.ky

    def to_lonlat(self, x, y):
        return self.lon0 + np.asarray(x) / self.kx, self.lat0 + np.asarray(y) / self.ky


def _bbox_from_args(args):
    if args.bbox:
        lon0, lat0, lon1, lat1 = (float(v) for v in args.bbox.split(","))
        return lon0, lat0, lon1, lat1
    with open(args.aoi, "r", encoding="utf-8") as f:
        gj = json.load(f)
    import shapely
    geoms = ([shapely.geometry.shape(ft["geometry"]) for ft in gj["features"]]
             if gj.get("type") == "FeatureCollection" else
             [shapely.geometry.shape(gj.get("geometry", gj))])
    return shapely.union_all(geoms).bounds


# ── Families (polygons in local metres) ─────────────────────────────────────

def hex_cells(w, h, dx):
    """Pointy-top regular hexagons of side `dx` covering [-w/2,w/2]×[-h/2,h/2]."""
    s = dx
    cw, ch = math.sqrt(3.0) * s, 1.5 * s            # column / row spacing
    corners = [(s * math.cos(math.radians(30 + 60 * k)),
                s * math.sin(math.radians(30 + 60 * k))) for k in range(6)]
    cells = []
    nr = int(math.ceil(h / ch)) + 2; nc = int(math.ceil(w / cw)) + 2
    for r in range(-nr, nr + 1):
        for c in range(-nc, nc + 1):
            cx = c * cw + (cw / 2.0 if r % 2 else 0.0); cy = r * ch
            if abs(cx) <= w / 2 and abs(cy) <= h / 2:
                cells.append(((cx, cy), [(cx + a, cy + b) for a, b in corners]))
    return cells


def rotcart_cells(w, h, dx, angle_deg):
    """Squares of side `dx` on a lattice rotated by `angle_deg` about the origin."""
    t = math.radians(angle_deg); ct, st = math.cos(t), math.sin(t)
    rot = lambda u, v: (u * ct - v * st, u * st + v * ct)
    n = int(math.ceil(math.hypot(w, h) / dx / 2.0)) + 1
    cells = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            cx, cy = rot(i * dx, j * dx)
            if abs(cx) <= w / 2 and abs(cy) <= h / 2:
                corners = [rot(i * dx + a, j * dx + b) for a, b in
                           ((-dx / 2, -dx / 2), (dx / 2, -dx / 2), (dx / 2, dx / 2), (-dx / 2, dx / 2))]
                cells.append(((cx, cy), corners))
    return cells


def voronoi_cells(w, h, dx, jitter, lloyd, seed, centre):
    """Voronoi tessellation of jittered square-lattice points, clipped to the box."""
    import shapely
    rng = np.random.default_rng(seed)
    nx, ny = max(2, int(round(w / dx))), max(2, int(round(h / dx)))
    xs = (np.arange(nx) + 0.5) * (w / nx) - w / 2
    ys = (np.arange(ny) + 0.5) * (h / ny) - h / 2
    pts = np.array([(x, y) for y in ys for x in xs])
    pts += rng.uniform(-jitter * dx, jitter * dx, pts.shape)
    pts[:, 0] = np.clip(pts[:, 0], -w / 2 * 0.999, w / 2 * 0.999)
    pts[:, 1] = np.clip(pts[:, 1], -h / 2 * 0.999, h / 2 * 0.999)
    box = shapely.box(-w / 2, -h / 2, w / 2, h / 2)
    for it in range(lloyd + 1):
        mp = shapely.multipoints(pts)
        vor = shapely.voronoi_polygons(mp, extend_to=box)
        polys = [shapely.intersection(p, box) for p in shapely.get_parts(vor)]
        # voronoi_polygons does not keep input order: match each cell to its generator
        order = shapely.STRtree(polys).query(shapely.points(pts), predicate="within")
        cell_of = dict(zip(order[0], order[1]))
        polys = [polys[cell_of[k]] for k in range(len(pts))]
        if it < lloyd:
            pts = np.array([[p.centroid.x, p.centroid.y] for p in polys])
    cells = []
    for g, p in zip(pts, polys):
        ring = list(p.exterior.coords)[:-1]
        c = (g[0], g[1]) if centre == "generator" else (p.centroid.x, p.centroid.y)
        cells.append((c, ring))
    return cells


# ── Common output ───────────────────────────────────────────────────────────

def _dedupe(ring):
    """Drop consecutive duplicate vertices (and a closing repeat)."""
    out = []
    for v in ring:
        if not out or v != out[-1]:
            out.append(v)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def _orient_ccw(ring):
    a = sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]))
    return ring if a > 0 else ring[::-1]


def edge_adjacency(rings, decimals):
    """Neighbours = cells sharing a polygon edge (rounded vertex pair)."""
    key = lambda v: (round(float(v[0]), decimals), round(float(v[1]), decimals))
    edge_cells = defaultdict(list)
    for i, ring in enumerate(rings):
        for a, b in zip(ring, ring[1:] + ring[:1]):
            ka, kb = key(a), key(b)
            if ka != kb:
                edge_cells[(min(ka, kb), max(ka, kb))].append(i)
    nbrs = [set() for _ in rings]
    bad = 0
    for cells in edge_cells.values():
        if len(cells) == 2:
            i, j = cells; nbrs[i].add(j); nbrs[j].add(i)
        elif len(cells) > 2:
            bad += 1
    if bad:
        raise ValueError(f"{bad} edges are shared by more than two cells")
    return [sorted(s) for s in nbrs]


def write_mesh(cells_xy, plane, out_path, family, extra_md):
    """cells_xy: list of ((cx, cy), [(x, y), ...]) in the local plane."""
    import a5_bridge as bridge
    # Snap every vertex to one canonical copy (1 µm grid key), so a vertex
    # shared by several cells is bitwise identical in each of them. GEOS can
    # produce copies that differ by ~1e-16, and FloodA5 matches shared
    # vertices after rounding (A5Grid._shared_edge), where such copies can
    # straddle a rounding boundary and drop the face.
    canon = {}
    def _snap(v):
        k = (round(v[0], 6), round(v[1], 6))
        return canon.setdefault(k, (float(v[0]), float(v[1])))
    rings_xy = [_orient_ccw(_dedupe([_snap(v) for v in ring])) for _, ring in cells_xy]
    # Convert to lon/lat once; neighbours share the same float vertices.
    rings_ll = []
    for ring in rings_xy:
        lon, lat = plane.to_lonlat([v[0] for v in ring], [v[1] for v in ring])
        rings_ll.append([(float(a), float(b)) for a, b in zip(lon, lat)])
    nbrs = edge_adjacency(rings_ll, 10)
    ids = [f"{k:016x}" for k in range(len(cells_xy))]
    cells = []
    for k, ((cx, cy), ring_xy) in enumerate(cells_xy):
        clon, clat = plane.to_lonlat(cx, cy)
        ring = rings_ll[k]
        xs, ys = zip(*rings_xy[k])
        area = 0.5 * abs(sum(xs[m] * ys[(m + 1) % len(xs)] - xs[(m + 1) % len(xs)] * ys[m]
                             for m in range(len(xs))))
        cells.append({
            "cell_id": ids[k], "resolution": 0,
            "center_lon": float(clon), "center_lat": float(clat),
            "boundary": [list(v) for v in ring] + [list(ring[0])],
            "neighbours": [ids[j] for j in nbrs[k]],
            "cell_area": area,
        })
    bridge.mesh_to_geoparquet(cells, out_path, id_scheme="polygon")
    bridge.write_id_scheme_metadata(out_path, "polygon", mesh_family=family, **extra_md)
    valence = np.bincount([len(n) for n in nbrs])
    return {"status": "ok", "output": out_path, "family": family, "cells": len(cells),
            "valence_histogram": valence.tolist(),
            "max_neighbours": int(max(len(n) for n in nbrs)),
            "mean_area_m2": float(np.mean([c["cell_area"] for c in cells]))}


def import_polygons(path):
    """Read any polygon file; return (cells_xy, plane) with centroid centres."""
    import geopandas as gpd
    gdf = gpd.read_file(path) if not path.endswith(".parquet") else gpd.read_parquet(path)
    if gdf.crs is not None and gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(4326)
    gdf = gdf.explode(index_parts=False)
    gdf = gdf[gdf.geometry.geom_type == "Polygon"]
    minx, miny, maxx, maxy = gdf.total_bounds
    plane = LocalPlane(0.5 * (minx + maxx), 0.5 * (miny + maxy))
    import shapely
    cells = []
    for g in gdf.geometry:
        lon, lat = zip(*list(g.exterior.coords)[:-1])
        x, y = plane.to_xy(lon, lat)
        ring = list(zip(x.tolist(), y.tolist()))
        c = shapely.Polygon(ring).centroid
        cells.append(((c.x, c.y), ring))
    return cells, plane


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("family", choices=["hex", "rotcart", "voronoi", "import"])
    ap.add_argument("--aoi"); ap.add_argument("--bbox"); ap.add_argument("--input")
    ap.add_argument("--dx", type=float, default=50.0)
    ap.add_argument("--angle", type=float, default=30.0)
    ap.add_argument("--jitter", type=float, default=0.3)
    ap.add_argument("--lloyd", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--centre", choices=["generator", "centroid"], default="generator")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.family == "import":
            if not args.input:
                raise ValueError("import needs --input FILE")
            cells, plane = import_polygons(args.input)
            md = {"source": os.path.basename(args.input)}
        else:
            if not (args.aoi or args.bbox):
                raise ValueError("give --aoi FILE or --bbox LON0,LAT0,LON1,LAT1")
            lon0, lat0, lon1, lat1 = _bbox_from_args(args)
            plane = LocalPlane(0.5 * (lon0 + lon1), 0.5 * (lat0 + lat1))
            (x0, x1), (y0, y1) = plane.to_xy([lon0, lon1], [lat0, lat1])
            w, h = x1 - x0, y1 - y0
            md = {"dx_m": args.dx}
            if args.family == "hex":
                cells = hex_cells(w, h, args.dx)
            elif args.family == "rotcart":
                cells = rotcart_cells(w, h, args.dx, args.angle); md["angle_deg"] = args.angle
            else:
                cells = voronoi_cells(w, h, args.dx, args.jitter, args.lloyd, args.seed, args.centre)
                md.update(jitter=args.jitter, lloyd=args.lloyd, seed=args.seed, centre=args.centre)
        print(json.dumps(write_mesh(cells, plane, args.out, args.family, md)), flush=True)
    except Exception as e:
        import traceback
        print(json.dumps({"error": f"{e}\n{traceback.format_exc()}"}), flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
