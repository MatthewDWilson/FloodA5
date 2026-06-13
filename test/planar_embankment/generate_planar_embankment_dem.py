#!/usr/bin/env python3
"""
generate_planar_embankment_dem.py — Planar slope + thin embankment DEM
=======================================================================

Generates a reproducible GeoTIFF DEM for the planar embankment SGS test.

Terrain description
-------------------
  1. Planar slope: linear gradient from west (high) to east (low).
     Default slope: 0.1% (1 m rise per 1000 m). Total relief across
     the 4 km domain: 4.0 m (west edge at z=4.0 m, east edge at z=0.0 m).

  2. Thin rectangular embankment: a sharp-edged ridge running N-S across
     the full domain at x = embankment_pos * Lx (default 50% of the way
     down the slope). Default height 2.0 m above the local slope, width 20 m.

  3. A single gap (hole) centred at y = 0.5*Ly, default 3 m wide —
     narrower than a res-18 A5 cell (~22 m diameter). Standard flow will
     not detect it (it uses mean cell elevation); SGS will route water
     through it because it pre-computes the minimum DEM elevation along
     each shared cell boundary.

Key design intent
-----------------
  - The planar slope is perfectly predictable: z(x) = z_west - slope * x.
    This makes analytical verification straightforward.
  - The embankment is rectangular (not Gaussian), giving sharp, well-defined
    cell-mean elevations that make it easy to reason about standard flow.
  - The 3 m notch is well below the res-18 cell scale (~22 m), ensuring
    it is genuinely sub-cell for the standard solver.
  - With closed boundaries and a constant upstream injection, the system
    fills predictably and the embankment blocks standard flow until
    WSE reaches the embankment crest, while SGS routes through the notch.

Key elevations (metres above datum):
  West edge (high side) z:  slope_pct/100 * domain_m_ew  (default 4.0 m)
  East edge (low side)  z:  0.0 m
  Embankment crest z:       z_at_emb + embankment_h (default varies with pos)
  Notch sill z:             z_at_emb  (hole goes through to the local slope level)
  Injection point z:        z at (0.1*Lx, 0.5*Ly) ≈ 3.6 m

Usage:
    cd <project_root>
    python test/planar_embankment/generate_planar_embankment_dem.py [options]

Options:
    --out DIR           Output directory (default: test/planar_embankment)
    --pixel-size M      DEM pixel size in metres (default: 5)
    --domain-km W H     Domain width and height in km (default: 4 2)
    --centre LON LAT    Domain centre lon lat (default: -0.017 51.0)
    --slope-pct S       East-west slope in percent (default: 0.1)
    --emb-height M      Embankment height above local slope (default: 2.0)
    --emb-width M       Embankment width in metres (default: 20.0)
    --emb-pos X         Embankment x-position as fraction of domain (default: 0.5)
    --notch-width M     Notch width in metres (default: 5.0)
    --help              Print this message

Output files (in --out directory):
    planar_embankment_dem.tif        GeoTIFF DEM, EPSG:4326, Float32
    planar_embankment_params.json    Key elevations and geometry for Julia assertions
    planar_embankment_aoi.geojson    AOI polygon for --meshgen
"""

import argparse, json, math, os
import numpy as np


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Planar slope + thin embankment DEM generator")
    p.add_argument("--out",          default="test/planar_embankment")
    p.add_argument("--pixel-size",   type=float, default=5.0,      metavar="M")
    p.add_argument("--domain-km",    type=float, nargs=2, default=[4.0, 2.0],
                   metavar=("W", "H"))
    p.add_argument("--centre",       type=float, nargs=2, default=[-0.017, 51.0],
                   metavar=("LON", "LAT"))
    p.add_argument("--slope-pct",    type=float, default=0.1,      metavar="S",
                   help="East-west slope in percent (default: 0.1)")
    p.add_argument("--emb-height",   type=float, default=2.0,      metavar="M",
                   help="Embankment height above local slope in metres (default: 2.0)")
    p.add_argument("--emb-width",    type=float, default=20.0,     metavar="M",
                   help="Embankment width in metres (default: 20.0)")
    p.add_argument("--emb-pos",      type=float, default=0.5,      metavar="X",
                   help="Embankment x-position as fraction of domain (default: 0.5)")
    p.add_argument("--notch-width",  type=float, default=5.0,      metavar="M",
                   help="Notch width in metres (default: 5.0)")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# DEM construction
# ---------------------------------------------------------------------------

def build_dem(domain_m_ew, domain_m_ns, pixel_m,
              slope_pct, emb_height, emb_width, emb_pos, notch_width):
    """
    Return (elev_array, params_dict).

    elev_array: 2-D float32, shape (ny, nx), north-up (row 0 = northernmost).

    Coordinate convention:
      x increases west→east  (col direction)
      y increases south→north (row direction; row 0 is northernmost in stored array)

    Elevation:
      z(x) = z_west - (slope_pct/100) * x
      z_west = slope_pct / 100 * domain_m_ew   (so east edge = 0.0 m)
    """
    nx = int(round(domain_m_ew / pixel_m))
    ny = int(round(domain_m_ns / pixel_m))
    Lx = domain_m_ew
    Ly = domain_m_ns

    # Pixel-centre coordinates in metres from SW corner
    xs = (np.arange(nx) + 0.5) * pixel_m    # E-W (col direction)
    ys = (np.arange(ny) + 0.5) * pixel_m    # S-N (row direction)
    xg, yg = np.meshgrid(xs, ys)            # shape (ny, nx)

    # 1. Planar slope: z decreases linearly west→east.
    #    z_west is set so the east edge is at z=0 (natural datum).
    slope   = slope_pct / 100.0             # m/m
    z_west  = slope * Lx                    # elevation at west edge
    z       = z_west - slope * xg           # (ny, nx)

    # 2. Rectangular embankment at x_emb ± emb_width/2.
    x_emb    = emb_pos * Lx
    half_w   = emb_width / 2.0
    in_emb   = (xg >= x_emb - half_w) & (xg <= x_emb + half_w)
    z_at_emb = z_west - slope * x_emb      # slope elevation at embankment centreline
    z[in_emb] = z_at_emb + emb_height      # flat top of embankment

    # Embankment crest elevation (uniform across full N-S extent before notch)
    emb_crest_z = z_at_emb + emb_height

    # 3. Notch: remove the embankment in a narrow N-S band at y = 0.5*Ly.
    #    The notch goes through to the natural slope level (not below it).
    #    Width is deliberately well below the res-18 cell scale (~22 m).
    y_centre  = 0.5 * Ly
    in_notch  = np.abs(yg - y_centre) <= 0.5 * notch_width
    # Within the notch, restore the natural slope elevation (no embankment)
    z[in_notch & in_emb] = z_at_emb        # = natural slope at this x

    # Notch sill is the slope elevation at the embankment centreline.
    notch_sill_z = z_at_emb

    z = np.clip(z, 0.0, None).astype(np.float32)

    # 4. Injection point: upstream near x=0.1*Lx, centred N-S at y=0.5*Ly.
    #    Positioned well upstream of the embankment (10% of domain width).
    inj_x_m  = 0.10 * Lx
    inj_y_m  = 0.50 * Ly
    inj_xi   = int(inj_x_m / pixel_m)
    inj_yi   = int(inj_y_m / pixel_m)
    inj_xi   = min(inj_xi, nx - 1)
    inj_yi   = min(inj_yi, ny - 1)

    params = {
        "domain_m_ew":      domain_m_ew,
        "domain_m_ns":      domain_m_ns,
        "pixel_m":          pixel_m,
        "nx":               nx,
        "ny":               ny,
        "slope_pct":        slope_pct,
        "z_west":           float(z_west),
        "z_east":           0.0,
        "emb_pos_frac":     emb_pos,
        "emb_x_m":          float(x_emb),
        "emb_width_m":      emb_width,
        "emb_height_m":     emb_height,
        "emb_crest_z":      float(emb_crest_z),
        "notch_width_m":    notch_width,
        "notch_y_centre_m": float(y_centre),
        "notch_sill_z":     float(notch_sill_z),
        "injection_x_m":    float(inj_x_m),
        "injection_y_m":    float(inj_y_m),
        "injection_z":      float(z[inj_yi, inj_xi]),
        "z_min_global":     float(z.min()),
        "z_max_global":     float(z.max()),
    }
    return z, params


# ---------------------------------------------------------------------------
# GeoTIFF writer — rasterio preferred, stdlib fallback (same as existing test)
# ---------------------------------------------------------------------------

def write_geotiff(path, z_arr, lon0, lat0, lon1, lat1):
    """Write a Float32 GeoTIFF with EPSG:4326 CRS and north-up orientation."""
    ny, nx = z_arr.shape
    try:
        import rasterio
        from rasterio.transform import from_bounds
        transform = from_bounds(lon0, lat0, lon1, lat1, nx, ny)
        with rasterio.open(path, "w", driver="GTiff",
                           height=ny, width=nx, count=1, dtype="float32",
                           crs="EPSG:4326", transform=transform,
                           nodata=float("nan")) as dst:
            dst.write(np.flipud(z_arr), 1)
        print(f"  [rasterio] {path}")
        return
    except ImportError:
        pass

    import struct
    pixel_lon = (lon1 - lon0) / nx
    pixel_lat = (lat1 - lat0) / ny
    data = np.flipud(z_arr).astype("<f4").tobytes()

    SHORT, LONG, DOUBLE = 3, 4, 12

    def pack_entry(tag, typ, count, val_or_off):
        if typ == SHORT:
            return struct.pack("<HHI", tag, typ, count) + struct.pack("<HH",
                               val_or_off & 0xFFFF, 0)
        return struct.pack("<HHII", tag, typ, count, val_or_off)

    DATA_OFFSET = 4096
    extra_off   = DATA_OFFSET + len(data)
    geo_scale   = struct.pack("<3d", pixel_lon, pixel_lat, 0.0)
    geo_tie     = struct.pack("<6d", 0.0, 0.0, 0.0, lon0, lat1, 0.0)
    geo_keys    = struct.pack("<16H", 1, 1, 0, 3,
                              1024, 0, 1, 2, 1025, 0, 1, 1, 2048, 0, 1, 4326)
    off_scale, off_tie = extra_off, extra_off + len(geo_scale)
    off_keys    = off_tie + len(geo_tie)

    tags = sorted([
        pack_entry(256,   LONG,  1, nx),
        pack_entry(257,   LONG,  1, ny),
        pack_entry(258,   SHORT, 1, 32),
        pack_entry(259,   SHORT, 1, 1),
        pack_entry(262,   SHORT, 1, 1),
        pack_entry(273,   LONG,  1, DATA_OFFSET),
        pack_entry(278,   LONG,  1, ny),
        pack_entry(279,   LONG,  1, len(data)),
        pack_entry(284,   SHORT, 1, 1),
        pack_entry(339,   SHORT, 1, 3),
        pack_entry(33550, DOUBLE, 3, off_scale),
        pack_entry(33922, DOUBLE, 6, off_tie),
        pack_entry(34735, SHORT, len(geo_keys) // 2, off_keys),
    ], key=lambda b: struct.unpack("<H", b[:2])[0])

    ifd_data = struct.pack("<H", len(tags)) + b"".join(tags) + struct.pack("<I", 0)
    pad = bytes(max(0, DATA_OFFSET - 8 - len(ifd_data)))
    with open(path, "wb") as f:
        f.write(b"II")
        f.write(struct.pack("<H", 42))
        f.write(struct.pack("<I", 8))
        f.write(ifd_data)
        f.write(pad)
        f.write(data)
        f.write(geo_scale)
        f.write(geo_tie)
        f.write(geo_keys)
    print(f"  [stdlib TIFF] {path}")


# ---------------------------------------------------------------------------
# AOI GeoJSON
# ---------------------------------------------------------------------------

def write_aoi(path, lon0, lat0, lon1, lat1):
    geojson = {
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[lon0, lat0], [lon1, lat0], [lon1, lat1],
                              [lon0, lat1], [lon0, lat0]]]
        }
    }
    with open(path, "w") as f:
        json.dump(geojson, f, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args    = parse_args(argv)
    os.makedirs(args.out, exist_ok=True)

    domain_m_ew   = args.domain_km[0] * 1000.0
    domain_m_ns   = args.domain_km[1] * 1000.0
    pixel_m       = args.pixel_size
    centre_lon, centre_lat = args.centre

    m_per_deg_lon = 111320.0 * math.cos(math.radians(centre_lat))
    m_per_deg_lat = 111320.0
    half_lon = (domain_m_ew / 2.0) / m_per_deg_lon
    half_lat = (domain_m_ns / 2.0) / m_per_deg_lat
    lon0, lon1 = centre_lon - half_lon, centre_lon + half_lon
    lat0, lat1 = centre_lat - half_lat, centre_lat + half_lat

    print("Planar embankment DEM generator")
    print(f"  Domain       : {domain_m_ew:.0f} m × {domain_m_ns:.0f} m")
    print(f"  Extent       : lon [{lon0:.5f}, {lon1:.5f}]  lat [{lat0:.5f}, {lat1:.5f}]")
    print(f"  Pixel size   : {pixel_m} m")
    print(f"  Slope        : {args.slope_pct}%  ({args.slope_pct/100:.4f} m/m)")
    print(f"  Embankment   : h={args.emb_height} m, w={args.emb_width} m "
          f"at x={args.emb_pos:.0%}")
    print(f"  Notch        : {args.notch_width} m wide at domain centre")

    z, params = build_dem(domain_m_ew, domain_m_ns, pixel_m,
                          args.slope_pct, args.emb_height,
                          args.emb_width, args.emb_pos, args.notch_width)

    # Compute injection lon/lat from pixel coordinates
    inj_lon = lon0 + params["injection_x_m"] / m_per_deg_lon
    inj_lat = lat0 + params["injection_y_m"] / m_per_deg_lat
    params.update({
        "lon0": lon0, "lon1": lon1, "lat0": lat0, "lat1": lat1,
        "centre_lon": centre_lon, "centre_lat": centre_lat,
        "injection_lon": inj_lon, "injection_lat": inj_lat,
        # Downstream partition boundary: just east of embankment
        # Used by Julia tests to split cells into upstream/downstream.
        "partition_x_m": params["emb_x_m"] + params["emb_width_m"] / 2.0 + pixel_m,
    })
    params["partition_lon"] = lon0 + params["partition_x_m"] / m_per_deg_lon

    print(f"  DEM size     : {params['nx']} × {params['ny']} pixels  "
          f"({params['nx']*params['ny']//1000}k)")
    print(f"  z range      : {params['z_min_global']:.3f} – {params['z_max_global']:.3f} m")
    print(f"  z west/east  : {params['z_west']:.3f} / {params['z_east']:.3f} m")
    print(f"  Emb crest z  : {params['emb_crest_z']:.3f} m")
    print(f"  Notch sill z : {params['notch_sill_z']:.3f} m  "
          f"(= slope at x={params['emb_x_m']:.0f} m)")
    print(f"  Injection    : z={params['injection_z']:.3f} m  "
          f"(lon={inj_lon:.5f}, lat={inj_lat:.5f})")

    dem_path    = os.path.join(args.out, "planar_embankment_dem.tif")
    params_path = os.path.join(args.out, "planar_embankment_params.json")
    aoi_path    = os.path.join(args.out, "planar_embankment_aoi.geojson")

    write_geotiff(dem_path, z, lon0, lat0, lon1, lat1)
    with open(params_path, "w") as f:
        json.dump(params, f, indent=2)
    print(f"  Params JSON  : {params_path}")
    write_aoi(aoi_path, lon0, lat0, lon1, lat1)
    print(f"  AOI GeoJSON  : {aoi_path}")
    print("Done.")


if __name__ == "__main__":
    main()
