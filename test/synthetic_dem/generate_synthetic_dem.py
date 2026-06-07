#!/usr/bin/env python3
"""
generate_synthetic_dem.py — Synthetic DEM for FloodA5 SGS validation
======================================================================

Generates a reproducible GeoTIFF DEM with known geometry:

  1. A parabolic bowl (lowest on the upstream/west side, rising to east).
  2. A transverse embankment crossing the full domain N-S at x ≈ 0.45*Lx,
     rising 1.5 m above the local bowl floor.
  3. A notch (gap + lowered floor) in the embankment centred at y = 0.5*Ly,
     narrow (default 300 m), with sill at params["notch_sill_z"] (below embankment crest).
     The notch is carved only in the embankment body (x >= x_emb - sigma) to avoid
     flattening the upstream basin to the notch-sill elevation.
  4. Injection point at x = 0.2*Lx, y = 0.5*Ly (upstream bowl, clear of notch).

Key elevations (metres above datum):
  Injection point z:   params["injection_z"]        (≈ 0.20 m)
  Notch sill:          params["notch_sill_z"]        (default 0.70 m)
  Embankment crest:    params["embankment_crest_z"]  (default 1.50 m)

These values are written to synthetic_dem_params.json so the Julia test
can read them without re-parsing the GeoTIFF.

Usage:
    cd <project_root>
    python test/synthetic_dem/generate_synthetic_dem.py [options]

Options:
    --out DIR           Output directory (default: test/synthetic_dem)
    --pixel-size M      DEM pixel size in metres (default: 5)
    --domain-km W H     Domain width and height in km (default: 4 2)
    --centre LON LAT    Domain centre lon lat (default: -0.017 51.0)
    --embankment-m Z    Embankment crest height above local bowl (default: 1.5)
    --notch-sill-m Z    Notch sill height above local bowl (default: 0.7)
    --notch-width-m W   Notch width in metres (default: 300)
    --embankment-pos X  Embankment x-position as fraction of domain (default: 0.45)
    --help              Print this message

Output files (in --out directory):
    synthetic_dem.tif           GeoTIFF DEM, EPSG:4326, Float32
    synthetic_dem_params.json   Key elevations and geometry for Julia assertions
    synthetic_aoi.geojson       AOI polygon for --meshgen
"""

import argparse, json, math, os, sys
import numpy as np


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    p.add_argument("--out",            default="test/synthetic_dem")
    p.add_argument("--pixel-size",     type=float, default=5.0,    metavar="M")
    p.add_argument("--domain-km",      type=float, nargs=2, default=[4.0, 2.0],
                   metavar=("W","H"))
    p.add_argument("--centre",         type=float, nargs=2, default=[-0.017, 51.0],
                   metavar=("LON","LAT"))
    p.add_argument("--embankment-m",   type=float, default=1.5,   metavar="Z")
    p.add_argument("--notch-sill-m",   type=float, default=0.7,   metavar="Z")
    p.add_argument("--notch-width-m",  type=float, default=300.0, metavar="W")
    p.add_argument("--embankment-sigma-m", type=float, default=300.0, metavar="S",
                   help="Embankment Gaussian sigma in metres (default: 300). "
                        "Keep 3*sigma < injection_x_m to avoid overwriting injection point.")
    p.add_argument("--embankment-pos", type=float, default=0.45,  metavar="X")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# DEM construction
# ---------------------------------------------------------------------------

def build_dem(domain_m_ew, domain_m_ns, pixel_m,
              embankment_m, notch_sill_m, notch_width_m, embankment_pos,
              embankment_sigma_m=300.0):
    """
    Return (elev_array, params_dict).
    elev_array: 2-D float32, shape (ny, nx), north-up (row 0 = northernmost).
    """
    nx = int(round(domain_m_ew / pixel_m))
    ny = int(round(domain_m_ns / pixel_m))
    Lx, Ly = domain_m_ew, domain_m_ns

    # Pixel-centre coordinates in metres from SW corner
    xs = (np.arange(nx) + 0.5) * pixel_m   # E-W (col direction)
    ys = (np.arange(ny) + 0.5) * pixel_m   # S-N (row direction, south = low index)
    xg, yg = np.meshgrid(xs, ys)           # shape (ny, nx)

    # 1. Parabolic bowl: minimum at x=Lx/2, rises symmetrically to east and west.
    #    Max height at edges = 0.5 m.  Slight N-S tilt: +0.05 m north-to-south.
    bowl_ew = 0.5 * (2.0 * xg / Lx - 1.0) ** 2
    bowl_ns = 0.05 * (yg / Ly)
    z = bowl_ew + bowl_ns

    # 2. Gaussian embankment ridge at x = embankment_pos * Lx, sigma=75m.
    x_emb = embankment_pos * Lx
    sigma  = embankment_sigma_m
    ridge  = embankment_m * np.exp(-0.5 * ((xg - x_emb) / sigma) ** 2)
    z     += ridge

    # 3. Notch: suppress the ridge in a N-S band centred at y = 0.5*Ly.
    #    Carve only within the embankment body itself (x >= x_emb - sigma).
    #    The old code used |x - x_emb| < 3*sigma, which also overwrote the
    #    upstream bowl for 900m west of the embankment, flattening SGS cells
    #    to z_min=z_max=notch_sill.  That prevented WSE from ever rising above
    #    the sill (T2 deadlock).  The fix restricts carving to the eastern half
    #    of the Gaussian ridge where the embankment actually rises.
    y_centre   = 0.5 * Ly
    in_notch   = np.abs(yg - y_centre) < 0.5 * notch_width_m
    in_emb     = (xg >= x_emb - sigma) & (xg <= x_emb + 3.0 * sigma)
    bowl_at_emb = 0.5 * (2.0 * x_emb / Lx - 1.0) ** 2  # bowl height at emb centreline
    notch_elev  = bowl_at_emb + notch_sill_m
    z[in_notch & in_emb] = notch_elev

    z = np.clip(z, 0.0, None).astype(np.float32)

    # 4. Injection point: upstream bowl at x=0.2*Lx, y=0.5*Ly (clear of notch/emb).
    inj_x_m = 0.2 * Lx
    inj_y_m = 0.5 * Ly
    inj_xi  = int(inj_x_m / pixel_m)
    inj_yi  = int(inj_y_m / pixel_m)

    params = {
        "domain_m_ew":          domain_m_ew,
        "domain_m_ns":          domain_m_ns,
        "pixel_m":              pixel_m,
        "nx":                   nx,
        "ny":                   ny,
        "embankment_pos_frac":  embankment_pos,
        "embankment_x_m":       float(x_emb),
        "embankment_sigma_m":   float(sigma),
        "embankment_m":         embankment_m,
        "embankment_crest_z":   float(bowl_at_emb + embankment_m),
        "notch_sill_z":         float(notch_elev),
        "notch_width_m":        notch_width_m,
        "notch_y_centre_m":     float(y_centre),
        "injection_x_m":        float(inj_x_m),
        "injection_y_m":        float(inj_y_m),
        "injection_z":          float(z[inj_yi, inj_xi]),
        "z_min_global":         float(z.min()),
        "z_max_global":         float(z.max()),
    }
    return z, params


# ---------------------------------------------------------------------------
# GeoTIFF writer — uses rasterio if available, else a minimal stdlib TIFF
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
            # rasterio expects row 0 = north: flip our south-up array
            dst.write(np.flipud(z_arr), 1)
        print(f"  [rasterio] {path}")
        return
    except ImportError:
        pass

    # Minimal stdlib GeoTIFF (GDAL-readable, no compression)
    import struct

    pixel_lon = (lon1 - lon0) / nx
    pixel_lat = (lat1 - lat0) / ny   # positive step size

    # Data is stored north-up (row 0 = northernmost = lat1)
    data = np.flipud(z_arr).astype("<f4").tobytes()

    # TIFF tag constants
    SHORT, LONG, DOUBLE = 3, 4, 12

    def pack_entry(tag, typ, count, val_or_off):
        if typ == SHORT:
            return struct.pack("<HHI", tag, typ, count) + struct.pack("<HH", val_or_off & 0xFFFF, 0)
        else:
            return struct.pack("<HHII", tag, typ, count, val_or_off)

    DATA_OFFSET = 4096  # image data after this fixed offset

    # Extra data appended after DATA_OFFSET + len(data)
    extra_off   = DATA_OFFSET + len(data)

    geo_scale  = struct.pack("<3d", pixel_lon, pixel_lat, 0.0)  # ModelPixelScaleTag
    geo_tie    = struct.pack("<6d", 0.0, 0.0, 0.0, lon0, lat1, 0.0)  # ModelTiepointTag (UL corner)
    # GeoKeyDirectory: 3 keys (GTModelType=2, GTRasterType=1, GeographicType=4326)
    geo_keys   = struct.pack("<16H",
                     1, 1, 0, 3,
                     1024, 0, 1, 2,
                     1025, 0, 1, 1,
                     2048, 0, 1, 4326)

    off_scale  = extra_off
    off_tie    = off_scale + len(geo_scale)
    off_keys   = off_tie   + len(geo_tie)

    tags = sorted([
        pack_entry(256,   LONG,  1, nx),                    # ImageWidth
        pack_entry(257,   LONG,  1, ny),                    # ImageLength
        pack_entry(258,   SHORT, 1, 32),                    # BitsPerSample
        pack_entry(259,   SHORT, 1, 1),                     # Compression=none
        pack_entry(262,   SHORT, 1, 1),                     # PhotometricInterp=BlackIsZero
        pack_entry(273,   LONG,  1, DATA_OFFSET),           # StripOffsets
        pack_entry(278,   LONG,  1, ny),                    # RowsPerStrip
        pack_entry(279,   LONG,  1, len(data)),             # StripByteCounts
        pack_entry(284,   SHORT, 1, 1),                     # PlanarConfig=chunky
        pack_entry(339,   SHORT, 1, 3),                     # SampleFormat=IEEE float
        pack_entry(33550, DOUBLE, 3, off_scale),            # ModelPixelScaleTag
        pack_entry(33922, DOUBLE, 6, off_tie),              # ModelTiepointTag
        pack_entry(34735, SHORT, len(geo_keys)//2, off_keys), # GeoKeyDirectoryTag
    ], key=lambda b: struct.unpack("<H", b[:2])[0])

    n_tags   = len(tags)
    ifd_data = struct.pack("<H", n_tags) + b"".join(tags) + struct.pack("<I", 0)

    pad = bytes(max(0, DATA_OFFSET - 8 - len(ifd_data)))

    with open(path, "wb") as f:
        f.write(b"II")                    # little-endian
        f.write(struct.pack("<H", 42))    # magic
        f.write(struct.pack("<I", 8))     # IFD at byte 8
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
            "coordinates": [[[lon0,lat0],[lon1,lat0],[lon1,lat1],[lon0,lat1],[lon0,lat0]]]
        }
    }
    with open(path, "w") as f:
        json.dump(geojson, f, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out, exist_ok=True)

    domain_m_ew = args.domain_km[0] * 1000.0
    domain_m_ns = args.domain_km[1] * 1000.0
    pixel_m     = args.pixel_size
    centre_lon, centre_lat = args.centre

    m_per_deg_lon = 111320.0 * math.cos(math.radians(centre_lat))
    m_per_deg_lat = 111320.0
    half_lon = (domain_m_ew / 2.0) / m_per_deg_lon
    half_lat = (domain_m_ns / 2.0) / m_per_deg_lat
    lon0, lon1 = centre_lon - half_lon, centre_lon + half_lon
    lat0, lat1 = centre_lat - half_lat, centre_lat + half_lat

    print(f"Synthetic DEM generator")
    print(f"  Domain       : {domain_m_ew:.0f} m × {domain_m_ns:.0f} m")
    print(f"  Extent       : lon [{lon0:.5f}, {lon1:.5f}]  lat [{lat0:.5f}, {lat1:.5f}]")
    print(f"  Pixel size   : {pixel_m} m")
    print(f"  Embankment   : {args.embankment_m} m crest @ x={args.embankment_pos:.0%}, sigma={args.embankment_sigma_m:.0f}m")
    print(f"  Notch        : sill={args.notch_sill_m} m, width={args.notch_width_m} m")

    z, params = build_dem(domain_m_ew, domain_m_ns, pixel_m,
                          args.embankment_m, args.notch_sill_m,
                          args.notch_width_m, args.embankment_pos,
                          args.embankment_sigma_m)

    # Compute injection lon/lat
    inj_lon = lon0 + params["injection_x_m"] / m_per_deg_lon
    inj_lat = lat0 + params["injection_y_m"] / m_per_deg_lat
    params.update({
        "lon0": lon0, "lon1": lon1, "lat0": lat0, "lat1": lat1,
        "centre_lon": centre_lon, "centre_lat": centre_lat,
        "injection_lon": inj_lon, "injection_lat": inj_lat,
    })

    print(f"  DEM size     : {params['nx']} × {params['ny']} pixels  "
          f"({params['nx']*params['ny']//1000}k)")
    print(f"  z range      : {params['z_min_global']:.3f} – {params['z_max_global']:.3f} m")
    print(f"  Injection z  : {params['injection_z']:.3f} m  "
          f"(lon={inj_lon:.5f}, lat={inj_lat:.5f})")
    print(f"  Notch sill   : {params['notch_sill_z']:.3f} m")
    print(f"  Emb crest    : {params['embankment_crest_z']:.3f} m")

    dem_path    = os.path.join(args.out, "synthetic_dem.tif")
    params_path = os.path.join(args.out, "synthetic_dem_params.json")
    aoi_path    = os.path.join(args.out, "synthetic_aoi.geojson")

    write_geotiff(dem_path, z, lon0, lat0, lon1, lat1)
    with open(params_path, "w") as f:
        json.dump(params, f, indent=2)
    print(f"  Params JSON  : {params_path}")
    write_aoi(aoi_path, lon0, lat0, lon1, lat1)
    print(f"  AOI GeoJSON  : {aoi_path}")
    print("Done.")


if __name__ == "__main__":
    main()
