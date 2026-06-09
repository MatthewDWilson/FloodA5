#!/usr/bin/env python3
"""
a5_cell_geometry.py
-------------------
For a given reference point (default: Christchurch, NZ), compute the
geodetic plan area and characteristic spacing (sqrt(area)) of A5 cells
at each resolution level from 5 to 19.

Usage:
    python a5_cell_geometry.py                    # uses default point
    python a5_cell_geometry.py 172.636 -43.531    # lon lat

Output:
    Tab-separated table printed to stdout.
    Also prints values formatted for the LaTeX paper table.

Method:
    For each resolution:
      1. Get the A5 cell ID containing the reference point.
      2. Retrieve the cell boundary vertices (lon/lat degrees).
      3. Project vertices to metres using a local equirectangular
         projection centred on the cell centroid (error < 0.1 % at
         mid-latitudes).
      4. Apply the Shoelace formula for polygon area.
      5. Characteristic spacing = sqrt(area).

Dependencies: a5  (pip install pya5)
"""

import sys
import math


def polygon_area_m2(lons, lats):
    """Geodetic plan area of a polygon via Shoelace on local equirectangular."""
    R = 6_371_000.0  # mean Earth radius, m
    cx = sum(lons) / len(lons)
    cy = sum(lats) / len(lats)
    cos_cy = math.cos(math.radians(cy))
    xs = [(lo - cx) * math.radians(1) * R * cos_cy for lo in lons]
    ys = [(la - cy) * math.radians(1) * R           for la in lats]
    n  = len(xs)
    area = sum(xs[i] * ys[(i+1) % n] - xs[(i+1) % n] * ys[i]
               for i in range(n))
    return abs(area) / 2.0


def fmt_area(m2):
    """Human-readable area with appropriate unit."""
    if m2 >= 1e6:
        return f"{m2/1e6:>7.2f} km²"
    elif m2 >= 1e4:
        return f"{m2/1e4:>7.2f} ha"
    else:
        return f"{m2:>8.1f} m²"


def main():
    import a5

    lon = float(sys.argv[1]) if len(sys.argv) > 1 else 172.636
    lat = float(sys.argv[2]) if len(sys.argv) > 2 else -43.531

    print(f"\nReference point: lon={lon}, lat={lat}")
    print(f"{'Res':>4}  {'Area':>14}  {'Spacing (√A)':>14}  {'Cell ID':>18}")
    print("-" * 60)

    rows = []
    for res in range(5, 21):
        cell  = a5.lonlat_to_cell((lon, lat), res)
        bnd   = a5.cell_to_boundary(cell)
        lons  = [v[0] for v in bnd]
        lats  = [v[1] for v in bnd]
        area  = polygon_area_m2(lons, lats)
        spac  = math.sqrt(area)
        cid   = a5.u64_to_hex(cell)
        rows.append((res, area, spac, cid))
        print(f"{res:>4}  {fmt_area(area):>14}  {spac:>11.1f} m  {cid:>18}")

    print()
    print("LaTeX table rows (for \\tabular):")
    print("-" * 60)
    for res, area, spac, _ in rows:
        if area >= 1e6:
            area_str = f"\\SI{{{area/1e6:.2f}}}{{\\km\\squared}}"
        elif area >= 1e4:
            area_str = f"\\SI{{{area/1e4:.1f}}}{{\\ha}}"
        else:
            area_str = f"\\SI{{{area:.0f}}}{{\\m\\squared}}"
        if spac >= 1000:
            spac_str = f"\\SI{{{spac/1000:.2f}}}{{\\km}}"
        else:
            spac_str = f"\\SI{{{spac:.0f}}}{{\\m}}"
        print(f"{res} & {area_str} & {spac_str} \\\\")


if __name__ == "__main__":
    main()
