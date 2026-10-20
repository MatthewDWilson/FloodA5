"""
migrate_mesh_ids.py
-------------------
Rewrite the cell IDs of a FloodA5 GeoParquet mesh from the pya5 0.x ID scheme
("a5-v0") to the pya5 1.x scheme ("a5-v1") with `a5.migrate`, so existing
meshes (and their sampled elevation, SGS tables, …) do not need regenerating.
WP-M1, 2026-10 — see FloodA5_Review_and_ForwardPlan_2026-10.md §9.

Only `cell_id` and `neighbours` change. Every other column, the geometry,
the row order and the GeoParquet metadata are kept, so a simulation on the
migrated mesh is bit-identical to one on the original. Each migrated cell's
centre is checked against the stored centre (pya5 1.x `cell_to_lonlat`), and
the schema metadata records `{"a5_id_scheme": "a5-v1", "migrated_from":
"a5-v0", ...}`. A `.prj` sidecar is copied if present.

Requires pya5 >= 1.0. Usage:
    python mesh/migrate_mesh_ids.py IN.parquet OUT.parquet [--tol-deg 1e-9]

Prints one JSON line with a summary; exits non-zero on any failure.
"""

import json
import os
import shutil
import sys


def _fail(msg):
    print(json.dumps({"error": msg}), flush=True)
    sys.exit(1)


def migrate_mesh(in_path: str, out_path: str, tol_deg: float = 1e-9) -> dict:
    import a5
    import pyarrow as pa
    import pyarrow.parquet as pq
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import a5_bridge as bridge

    if bridge.A5_MAJOR < 1:
        raise RuntimeError(f"pya5 {bridge.PYA5_VERSION} is installed; migration needs pya5 >= 1.0")
    if os.path.abspath(in_path) == os.path.abspath(out_path):
        raise RuntimeError("OUT must differ from IN (the original is kept)")

    scheme = bridge.read_id_scheme(in_path) or "a5-v0"
    if scheme != "a5-v0":
        raise RuntimeError(f"{in_path} has ID scheme '{scheme}'; only 'a5-v0' meshes can be migrated")

    tbl = pq.read_table(in_path)
    old_hex = tbl.column("cell_id").to_pylist()
    mapping = {h: bridge._cell_int_to_hex(a5.migrate(int(h, 16))) for h in old_hex}
    new_hex = [mapping[h] for h in old_hex]
    if len(set(new_hex)) != len(new_hex):
        raise RuntimeError("migrated IDs are not unique")

    # Centre check: the migrated ID must name the same cell.
    lons = tbl.column("center_lon").to_pylist()
    lats = tbl.column("center_lat").to_pylist()
    max_dev = 0.0
    for h, lon, lat in zip(new_hex, lons, lats):
        c_lon, c_lat = a5.cell_to_lonlat(int(h, 16))
        dlon = abs(((c_lon - lon) + 180.0) % 360.0 - 180.0)
        max_dev = max(max_dev, dlon, abs(c_lat - lat))
    if max_dev > tol_deg:
        raise RuntimeError(f"migrated cell centres differ from the stored centres by up to "
                           f"{max_dev:.3e} deg (> {tol_deg:g}); is the input really 'a5-v0'?")

    cols = {"cell_id": pa.array(new_hex, type=tbl.schema.field("cell_id").type)}
    if "neighbours" in tbl.column_names:
        nb_new = []
        for nbrs in tbl.column("neighbours").to_pylist():
            nb_new.append(None if nbrs is None else
                          [mapping.get(n) or bridge._cell_int_to_hex(a5.migrate(int(n, 16)))
                           for n in nbrs])
        cols["neighbours"] = pa.array(nb_new, type=tbl.schema.field("neighbours").type)
    for name, arr in cols.items():
        tbl = tbl.set_column(tbl.column_names.index(name), tbl.schema.field(name), arr)

    pq.write_table(tbl, out_path)
    bridge.write_id_scheme_metadata(out_path, "a5-v1", migrated_from="a5-v0",
                                    migrated_by="mesh/migrate_mesh_ids.py")
    prj_in = os.path.splitext(in_path)[0] + ".prj"
    if os.path.isfile(prj_in):
        shutil.copyfile(prj_in, os.path.splitext(out_path)[0] + ".prj")
    return {"status": "ok", "input": in_path, "output": out_path, "cells": len(new_hex),
            "max_centre_dev_deg": max_dev, "pya5_version": bridge.PYA5_VERSION,
            "a5_id_scheme": "a5-v1"}


def main():
    args = sys.argv[1:]
    tol = 1e-9
    if "--tol-deg" in args:
        i = args.index("--tol-deg"); tol = float(args[i + 1]); del args[i:i + 2]
    if len(args) != 2:
        _fail("Usage: migrate_mesh_ids.py IN.parquet OUT.parquet [--tol-deg 1e-9]")
    try:
        print(json.dumps(migrate_mesh(args[0], args[1], tol)), flush=True)
    except Exception as e:
        _fail(str(e))


if __name__ == "__main__":
    main()
